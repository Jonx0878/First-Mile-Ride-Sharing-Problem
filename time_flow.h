#pragma once
#pragma once
#include <algorithm>
#include <gurobi_c++.h>
#include <set>

#include "cvrpsep/cnstrmgr.h"
#include "col_gen.h"
#include "load.h"
#include "rc.h"
#include "utils.h"


class TF {
private:

public:
    // Instance Properties
    std::string name;
    int capacity;
    int no_veh;
    int no_cust;
    int no_cust_d;
    int time_ub;
    int no_dir_veh_vars = 0;
    int no_total_veh_vars = 0;
    int no_total_vars = 0;
    double max_cap;
    GRBEnv env = GRBEnv();

    // Vehicle Data
    std::vector<Coords> veh_coords;
    std::vector<int> veh_occ;
    std::vector<int> veh_arr_time;
    std::vector<std::vector<int>> veh_dist; // Distance to each customer
    std::set<int> vehicles;

    // Customer Data
    std::vector<Coords> cust_coords;
    std::vector<int> demand;
    std::vector<int> cust_arr_time;
    std::vector<std::vector<int>> cust_dist; // Upper triangular matrix of distances - customer zero is the destination
    std::set<int> customers;
    std::set<int> dest_and_cust;

    // Model Properties
    GRBModel* model = nullptr;

    std::vector<GRBVar> flow;
    std::vector<GRBVar> vars;
    std::vector<std::vector<GRBVar>> xi;
    int no_base_constraints = 0;

    //Valid inequalities
    std::vector<std::set<int>> rounded_cap;
    std::set<int> active_rc_ineqs;
    int total_rcc = 0;


    // Edges
    std::set<int> var_indices;
    std::set<int> var_indices_subset;
    std::set<int> flow_indices;

    // CVRPSEP
    CnstrMgrPointer MyCutsCMP;
    CnstrMgrPointer MyOldCutsCMP;

    // Timings
    double solve_time = 0;
    double model_construction_time = 0;
    double graph_construction_time = 0;
    double col_gen_time = 0;
    double RC_separation_time = 0;
    double add_vars_time = 0;


    // Initializes instance
    TF(const std::string& filename) {
        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        capacity = inst_data.capacity[0];
        no_veh = inst_data.no_veh;
        no_cust = inst_data.no_cust;
        no_cust_d = no_cust + 1;
        time_ub = inst_data.time_ub;
        max_cap = capacity - *std::ranges::min_element(inst_data.veh_occ);

        veh_coords = inst_data.veh_coords;
        veh_occ = inst_data.veh_occ;
        veh_arr_time = inst_data.veh_arr_time;
        veh_dist = inst_data.veh_dist;

        cust_coords = inst_data.cust_coords;
        demand = inst_data.demand;
        cust_arr_time = inst_data.cust_arr_time;
        cust_dist = inst_data.cust_dist;

        // Populate vertices, customers, and deste_and_cust
        for (int k = 0; k < no_veh; k++) {
            vehicles.insert(k);
        }
        for (int i = 0; i <= no_cust; i++) {
            dest_and_cust.insert(i);
            if (i > 0) customers.insert(i);
        }

        rounded_cap.resize(no_cust * (no_cust - 1) / 2);

        // Setup No customers
        no_dir_veh_vars = no_veh * (no_cust_d);
        no_total_veh_vars = no_dir_veh_vars + no_veh * (no_cust_d * (no_cust_d - 1) / 2);
        no_total_vars = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2);
        for (int idx = 0; idx < no_total_vars; idx++) {
            var_indices.insert(var_indices.end(), idx);
        }

        for (const int& i : dest_and_cust) {
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                const int idx = i * no_cust_d + j;
                flow_indices.insert(flow_indices.end(), idx);
            }
        }

        // Setup Constraint Managers
        CMGR_CreateCMgr(&MyCutsCMP, 100);
        CMGR_CreateCMgr(&MyOldCutsCMP, 100);
    }

    void preprocess() {
        auto start = std::chrono::high_resolution_clock::now();
        // Preprocess YPS Vars
        std::cout << "Potential # Edge Variables: " << var_indices.size() << std::endl;
        std::cout << "Preprocesing..." << std::endl;
        std::set<int> indices_removed;
        std::set<int> flow_removed;

        // Remove vehicle edges with too large demand or time constraint
        for (const int& k : vehicles) {
            int remaning_cap = capacity - veh_occ[k];
            int k_times_cust_d = k * no_cust_d;
            for (const int& j : customers) {
                if (veh_dist[k][j] + cust_dist[j][0] > std::min(veh_arr_time[k], cust_arr_time[j]) || demand[j] > remaning_cap) {
                    int idx = k_times_cust_d + j;
                    indices_removed.insert(idx);

                    for (const int& i : dest_and_cust) {
                        if (i < j) {
                            idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                            indices_removed.insert(idx);
                        }
                        else if (i > j) {
                            idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
                            indices_removed.insert(idx);
                        }
                    }
                }
            }
        }

        // Remove customer edges
        for (int i : customers) {
            for (int j = i + 1; j < no_cust_d; j++) {
                int dem = demand[i] + demand[j];
                int min_time = std::min(cust_arr_time[i], cust_arr_time[j]);

                bool feasible = false;
                int min_veh_dist = time_ub;
                for (const int& k : vehicles) {
                    if (dem <= capacity - veh_occ[k] && veh_dist[k][i] + cust_dist[i][j] + cust_dist[j][0] <= std::min(veh_arr_time[k], min_time)) {
                        min_veh_dist = std::min(min_veh_dist, veh_dist[k][i]);
                        feasible = true;
                    }
                    else if (dem <= capacity - veh_occ[k] && veh_dist[k][j] + cust_dist[j][i] + cust_dist[i][0] <= std::min(veh_arr_time[k], min_time)) {
                        min_veh_dist = std::min(min_veh_dist, veh_dist[k][j]);
                        feasible = true;
                    }
                }
                double max_time = fmax(cust_arr_time[i], cust_arr_time[j]);
                if (!feasible) {
                    for (int k : vehicles) {
                        int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                        indices_removed.insert(idx);
                    }
                }
                else {
                    for (int k : vehicles) {
                        int min_dist = std::min(veh_dist[k][i] + cust_dist[i][j] + cust_dist[j][0], veh_dist[k][j] + cust_dist[j][i] + cust_dist[i][0]);
                        if (dem > capacity - veh_occ[k] || min_dist > fmax(veh_arr_time[k], max_time)) {
                            int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                            indices_removed.insert(idx);
                        }
                    }
                }
            }
        }

        // Remove Symmetric edge variables
        for (const int& i : dest_and_cust) {
            for (int j = i + 1; j < no_cust_d; j++) {
                bool exists = false;
                for (const int& k : vehicles) {
                    int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (std::find(indices_removed.begin(), indices_removed.end(), idx) == indices_removed.end()) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) {
                    int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    indices_removed.insert(indices_removed.end(), idx);
                    if (i > 0) {
                        idx = i * no_cust_d + j;
                        flow_removed.insert(idx);
                        idx = j * no_cust_d + i;
                        flow_removed.insert(idx);
                    }
                }
            }
        }

        std::set<int> temp_set;
        std::set_difference(var_indices.begin(), var_indices.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices = temp_set;
        temp_set.clear();
        std::set_difference(var_indices_subset.begin(), var_indices_subset.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices_subset = temp_set;
        temp_set.clear();
        std::set_difference(flow_indices.begin(), flow_indices.end(), flow_removed.begin(), flow_removed.end(), std::inserter(temp_set, temp_set.begin()));
        flow_indices = temp_set;
        std::cout << "Variables Removed: " << indices_removed.size() + flow_removed.size() << " # Variables left: " << var_indices.size() + flow_indices.size() << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Creates the TCF formulation
    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);

        // Add Variables
        if (var_indices_subset.size() == 0) var_indices_subset = var_indices;

        double obj = 0;
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        std::string name;

        // flow variables
        flow.resize(no_cust_d * no_cust_d);

        for (const int& idx : flow_indices) {
            const int i = floor(idx / float(no_cust_d));
            const int j = idx - i * no_cust_d;
            name = "flow[" + itos(i) + "][" + itos(j) + "]";
            int idx = i * no_cust_d + j;
            double ub = (j > 0) ? max_cap - demand[i] : 2 * max_cap - 2;
            flow[idx] = model->addVar(0, time_ub, 0, GRB_CONTINUOUS, name);
        }

        vars.resize(no_total_vars);
        xi.resize(no_cust_d);
        for (const int& i : dest_and_cust) {
            xi[i].resize(no_veh);
        }
        for (const int& k : vehicles) {
            name = "xi[" + itos(0) + "][" + itos(k) + "]";
            xi[0][k] = model->addVar(0.0, 1.0, 0.0, vtype, name);
        }
        for (const int& idx : var_indices) {
            if (idx >= no_dir_veh_vars) break;
            // xi
            int k = floor(idx / no_cust_d);
            int i = idx - k * no_cust_d;
            obj = -cust_dist[0][i] * demand[i];
            name = "xi[" + itos(i) + "][" + itos(k) + "]";
            xi[i][k] = model->addVar(0.0, 1.0, obj, vtype, name);
        }

        for (const int& idx : var_indices_subset) {
            if (idx < no_dir_veh_vars) {
                // gamma_veh
                int k = floor(idx / no_cust_d);
                int j = idx - k * no_cust_d;
                obj = veh_dist[k][j];
                name = "gamma[" + itos(k) + "][" + itos(j) + "]";
                vars[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
            }
            else if (idx < no_total_veh_vars) {
                // chi
                int k = floor((idx - no_dir_veh_vars) / (no_cust_d * (no_cust_d - 1) / 2));
                int sym_idx = idx - no_dir_veh_vars - k * (no_cust_d * (no_cust_d - 1) / 2);
                int i = no_cust_d - 2 - int(floor(sqrt(-8 * sym_idx + 4 * no_cust_d * (no_cust_d - 1) - 7) / 2.0 - 0.5));
                int j = sym_idx + i * (i + 3) / 2 - no_cust_d * i + 1;
                obj = 0;
                name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, name);
            }
            else {
                // gamma_cust
                int sym_idx = idx - no_total_veh_vars;
                int i = no_cust_d - 2 - int(floor(sqrt(-8 * sym_idx + 4 * no_cust_d * (no_cust_d - 1) - 7) / 2.0 - 0.5));
                int j = sym_idx + i * (i + 3) / 2 - no_cust_d * i + 1;
                obj = cust_dist[i][j];
                name = "psi[" + itos(i) + "][" + itos(j) + "]";
                vars[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
            }
        }

        // Add Constraints
        GRBLinExpr lexpr, lexpr2, lexpr3, lexpr4, lexpr5;
        model->update();

        // 2.2
// Edge from vehicle if chosen
        for (const int& k : vehicles) {
            for (const int& i : dest_and_cust) {
                int idx = k * no_cust_d + i;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
            }
            model->addConstr(lexpr == xi[0][k]);
            lexpr.clear();
        }

        // 2.3
        // Capacity Constraints
        for (int k : vehicles) {
            for (const int& i : customers) {
                int idx = k * no_cust_d + i;
                if (var_indices.find(idx) != var_indices.end()) lexpr += demand[i] * xi[i][k];
            }
            int idx = k * no_cust_d + 0;
            if (var_indices_subset.find(idx) != var_indices_subset.end()) model->addConstr(lexpr <= (capacity - veh_occ[k]) * (xi[0][k] - vars[idx]));
            else model->addConstr(lexpr <= (capacity - veh_occ[k]) * (xi[0][k]));
            lexpr.clear();
        }

        // 2.4
        //Vehicle must reach destination if it already has passengers
        for (int k : vehicles) {
            model->addConstr(veh_occ[k] <= veh_occ[k] * xi[0][k]);
            lexpr.clear();
        }

        // 2.5
        // Subtour Elimination
        for (const int& k : vehicles) {
            for (const int& i : dest_and_cust) {
                int idx = k * no_cust_d + i;
                if (i > 0 && var_indices.find(idx) != var_indices.end()) lexpr2 += xi[i][k];
                for (int j = i + 1; j < no_cust_d; j++) {
                    int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                }
            }
        }
        model->addConstr(lexpr == lexpr2);
        lexpr.clear();
        lexpr2.clear();

        // 2.6
        // Customer serviced at most once
        for (int i : customers) {
            for (const int& k : vehicles) {
                int idx = k * no_cust_d + i;
                if (var_indices.find(idx) != var_indices.end()) lexpr += xi[i][k];
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // 2.7
        // 2 Edges per customers
        for (const int& i : customers) {
            for (const int& j : dest_and_cust) {
                if (i < j) {
                    int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                }
                else if (i > j) {
                    int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                }
            }
            for (const int& k : vehicles) {
                int veh_idx = k * no_cust_d + i;
                if (var_indices_subset.find(veh_idx) != var_indices_subset.end()) {
                    lexpr += vars[veh_idx];
                    lexpr2 += 2 * xi[i][k];
                }
            }
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
        }

        // 5.10
        for (const int& i : customers) {
            // Flow in node
            for (const int& k : vehicles) {
                int idx = k * no_cust_d + i;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) {
                    lexpr += (cust_arr_time[i] - veh_dist[k][i]) * vars[idx];
                    if (cust_arr_time[k] < veh_arr_time[k]) {
                        lexpr5 += 2 * cust_arr_time[i] * xi[i][k];
                    }
                }
            }
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                int idx = j * no_cust_d + i;
                if (flow_indices.find(idx) != flow_indices.end()) lexpr2 += flow[idx];
                idx = i * no_cust_d + j;
                if (flow_indices.find(idx) != flow_indices.end()) lexpr3 += flow[idx];
            }

            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                if (i < j) {
                    int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr4 += cust_dist[i][j] * vars[idx];
                }
                else if (i > j) {
                    int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr4 += cust_dist[j][i] * vars[idx];
                }
            }

            model->addConstr(lexpr + lexpr2 <= lexpr3 + lexpr4);
            lexpr.clear();
            lexpr2.clear();
            lexpr3.clear();
            lexpr4.clear();
            lexpr5.clear();

            // flow edge constraints
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                int idx;
                if (i < j) {
                    idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                }
                else if (i > j) {
                    idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr4 += cust_dist[j][i] * vars[idx];
                }
                if (var_indices_subset.find(idx) == var_indices_subset.end()) continue;

                model->addConstr(cust_dist[j][0] * vars[idx] <= flow[i * no_cust_d + j]);
                model->addConstr(cust_dist[i][0] * vars[idx] <= flow[j * no_cust_d + i]);

                model->addConstr(flow[i * no_cust_d + j] <= cust_arr_time[i] * vars[idx]);
                model->addConstr(flow[j * no_cust_d + i] <= cust_arr_time[j] * vars[idx]);

            }

            // Flow from depot
            int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - 0) * ((no_cust_d - 0) - 1) / 2 + i - 0 - 1;

            for (const int& k : vehicles) {
                int veh_idx = k * no_cust_d + i;
                lexpr3 -= veh_arr_time[k];
                if (var_indices_subset.find(veh_idx) != var_indices_subset.end()) {
                    lexpr += veh_arr_time[k] * xi[i][k];
                    lexpr3 += veh_arr_time[k] * xi[i][k];
                    lexpr4 += vars[veh_idx];
                }
                if (var_indices_subset.find(idx) != var_indices_subset.end()) {
                    lexpr2 += veh_arr_time[k] * vars[idx];
                    lexpr3 += veh_arr_time[k] * vars[idx];
                }

            }
            model->addConstr(flow[i] <= lexpr);
            model->addConstr(flow[i] <= lexpr2);
            model->addConstr(flow[i] >= lexpr3);
            lexpr.clear();
            lexpr2.clear();
            lexpr3.clear();

            model->addConstr(flow[i * no_cust_d] <= cust_arr_time[i] * (vars[idx] + lexpr4));
            lexpr4.clear();

        }

        for (const int& i : customers) {
            int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - 0) * ((no_cust_d - 0) - 1) / 2 + i - 0 - 1;
            if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
        }
        for (const int& k : vehicles) {
            int idx = k * no_cust_d;
            if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += xi[0][k];
        }
        model->addConstr(lexpr == lexpr2);
        lexpr.clear();
        lexpr2.clear();

        //// 5.11
        //for (const int& i : customers) {
        //    for (int j = i + 1; j < no_cust_d; j++) {
        //        int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
        //        if (var_indices_subset.find(idx) == var_indices_subset.end()) continue;
        //        for (const int& k : vehicles) {
        //            int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;;
        //            if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += (capacity - veh_occ[k]) * vars[idx];
        //        }
        //        lexpr2 += flow[i * no_cust_d + j];
        //        lexpr2 += flow[j * no_cust_d + i];
        //        model->addConstr(lexpr2 == lexpr);
        //        lexpr.clear();
        //        lexpr2.clear();
        //    }
        //}

        //// 5.12-5.13
        //for (const int& i : customers) {
        //    for (int j = i + 1; j < no_cust_d; j++) {
        //        int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
        //        if (var_indices_subset.find(idx) == var_indices_subset.end()) continue;
        //        lexpr += vars[idx];
        //        for (const int& k : vehicles) {
        //            int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
        //            if (var_indices_subset.find(idx) != var_indices_subset.end()) {
        //                lexpr2 += (capacity - veh_occ[k] - demand[i]) * vars[idx];
        //                lexpr3 += (capacity - veh_occ[k] - demand[j]) * vars[idx];
        //            }
        //            idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
        //            if (var_indices_subset.find(idx) != var_indices_subset.end()) {
        //                lexpr2 += (capacity - veh_occ[k] - demand[i]) * vars[idx];
        //                lexpr3 += (capacity - veh_occ[k] - demand[j]) * vars[idx];
        //            }
        //        }
        //        // 5.12
        //        model->addConstr(demand[j] * lexpr <= flow[i * no_cust_d + j]);
        //        model->addConstr(flow[i * no_cust_d + j] <= lexpr2);
        //        lexpr2.clear();
        //        // 5. 13
        //        model->addConstr(demand[i] * lexpr <= flow[j * no_cust_d + i]);
        //        model->addConstr(flow[j * no_cust_d + i] <= lexpr3);
        //        lexpr.clear();
        //        lexpr3.clear();
        //    }
        //}

        //// 5.14
        //for (const int& i : customers) {
        //    for (const int& k : vehicles) {
        //        int idx = k * no_cust_d + i;
        //        if (var_indices_subset.find(idx) != var_indices_subset.end()) {
        //            lexpr += (capacity - veh_occ[k] - demand[i]) * vars[idx];
        //        }
        //        idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - 0) * ((no_cust_d - 0) - 1) / 2 + i - 0 - 1;
        //        if (var_indices_subset.find(idx) != var_indices_subset.end()) {
        //            lexpr += (capacity - veh_occ[k] - demand[i]) * vars[idx];
        //        }
        //    }
        //    model->addConstr(flow[i * no_cust_d] <= lexpr);
        //    lexpr.clear();
        //}

        std::cout << "Model has been built" << std::endl;
        model->update();
        no_base_constraints = model->get(GRB_IntAttr_NumConstrs);
        auto stop = std::chrono::high_resolution_clock::now();
        model_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void optimise_model() {
        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        auto stop = std::chrono::high_resolution_clock::now();
        solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }
};