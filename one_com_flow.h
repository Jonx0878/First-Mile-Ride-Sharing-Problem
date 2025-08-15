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


class OCF {
private:

public:
    // Instance Properties
    std::string name;
    int capacity;
    int no_veh;
    int no_cust;
    int no_cust_d;
    int time_ub;
    int no_veh_vars = 0;
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
    OCF(const std::string& filename) {
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

        // Setup Constraint Managers
        CMGR_CreateCMgr(&MyCutsCMP, 100);
        CMGR_CreateCMgr(&MyOldCutsCMP, 100);

        no_veh_vars = no_veh * (no_cust + 1);
        no_total_vars = no_veh * (no_cust + 1) + no_veh * ((no_cust + 1) * (no_cust + 1) - no_cust - 1);
        for (int idx = 1; idx <= no_total_vars; idx++) {
            if (idx > no_veh_vars) {
                int k = floor((idx - 1 - no_veh * no_cust_d) / (pow(no_cust_d, 2) - no_cust - 1));
                int i = floor(
                    (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1)) / no_cust_d
                    + 1);
                int j = (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1) - (i - 1) * no_cust_d);
                if (i != j) var_indices.insert(idx);
            }
            else var_indices.insert(idx);
        }

        for (const int& i : customers) {
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                const int idx = i * no_cust_d + j;
                flow_indices.insert(flow_indices.end(), idx);
            }
        }
    }

    void set_starting_edges() {
        int no_edges = 15;

        // For YPS
        no_veh_vars = no_veh * (no_cust + 1);
        no_total_vars = no_veh * (no_cust + 1) + no_veh * ((no_cust + 1) * (no_cust + 1) - no_cust - 1);
        std::set<std::pair<int, int>> veh_loc;


        for (int k : vehicles) {
            int k_times_cust_d = k * no_cust_d;
            if (veh_occ[k] > 0) {
                var_indices_subset.insert(k_times_cust_d + 1);
            }
            std::vector<DistanceLocationPair> distances;
            for (int j : customers) {
                distances.push_back(DistanceLocationPair(veh_dist[k][j], j));
            }
            // Sort distances in ascending order
            std::sort(distances.begin(), distances.end(), [](const DistanceLocationPair& a, const DistanceLocationPair& b) {
                return a.distance < b.distance;
                });
            for (int l = 0; l < fmin(no_edges, distances.size()); l++) {
                int j = distances[l].location;
                if (j > 0) veh_loc.insert({ k, j });
                var_indices_subset.insert(k_times_cust_d + j + 1);
            }
        }

        std::vector<std::set<std::tuple<int, int, int>>> veh_loc2(capacity);
        for (auto& pair : veh_loc) {
            int k = pair.first;
            int i = pair.second;
            std::vector<DistanceLocationPair> distances;
            for (int j : customers) {
                if (i != j) distances.push_back(DistanceLocationPair(cust_dist[i][j], j));
            }
            std::sort(distances.begin(), distances.end(), [](const DistanceLocationPair& a, const DistanceLocationPair& b) {
                return a.distance < b.distance;
                });
            for (int l = 0; l < fmin(no_edges, distances.size()); l++) {
                int j = distances[l].location;
                int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                var_indices_subset.insert(idx);
                if (j > 0) {
                    veh_loc2[0].insert({ i, j, k });
                    idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + i + 1;
                    var_indices_subset.insert(idx);
                }
            }
        }

        for (int n = 0; n < capacity - 1; n++) {
            for (auto& pair : veh_loc2[n]) {
                int i = std::get<0>(pair);
                int m = std::get<1>(pair);
                int k = std::get<2>(pair);
                std::vector<DistanceLocationPair> distances;
                for (int j : customers) {
                    if (i != j && j != m) distances.push_back(DistanceLocationPair(cust_dist[i][j], j));
                }
                std::sort(distances.begin(), distances.end(), [](const DistanceLocationPair& a, const DistanceLocationPair& b) {
                    return a.distance < b.distance;
                    });
                for (int l = 0; l < fmin(no_edges, distances.size()); l++) {
                    int j = distances[l].location;
                    if (j > 0) veh_loc2[n + 1].insert({ i, j, k });
                    int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                    var_indices_subset.insert(idx);
                    idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + i + 1;
                    var_indices_subset.insert(idx);
                    if (n == capacity - 2 && j > 0) {
                        idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + 1;
                        var_indices_subset.insert(idx);
                    }
                }
            }
        }
    }

    void preprocess() {
        auto start = std::chrono::high_resolution_clock::now();
        // Preprocess YPS Vars
        std::cout << var_indices.size() << std::endl;
        std::cout << "Preprocesing..." << std::endl;
        std::set<int> indices_removed;
        std::set<int> flow_removed;
        std::vector<std::set<int>> flow_ind_veh(no_cust_d * no_cust_d);

        // Remove vehicle edges with too large demand or time constraint
        for (const int& k : vehicles) {
            int remaning_cap = capacity - veh_occ[k];
            int k_times_cust_d = k * no_cust_d;
            for (const int& j : customers) {
                //double t_max = fmax(veh_arr_time[k], cust_arr_time[j]);
                if (veh_dist[k][j] + cust_dist[j][0] > std::min(veh_arr_time[k], cust_arr_time[j]) || demand[j] > remaning_cap) {
                    int idx = k_times_cust_d + j + 1;
                    indices_removed.insert(idx);

                    idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + 1;
                    indices_removed.insert(idx);

                    idx = j * no_cust_d;
                    flow_ind_veh[idx].insert(k);

                    for (const int& i : customers) {
                        if (i != j) {
                            idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                            indices_removed.insert(idx);

                            idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + i + 1;
                            indices_removed.insert(idx);

                            idx = j * no_cust_d + i;
                            flow_ind_veh[idx].insert(k);
                            idx = i * no_cust_d + j;
                            flow_ind_veh[idx].insert(k);
                        }
                    }
                }
            }
        }

        // Remove customer edges
        for (const int& i : customers) {
            for (const int& j : customers) {
                if (i != j) {
                    int dist = cust_dist[i][j] + cust_dist[j][0];
                    int dem = demand[i] + demand[j];
                    int min_time = std::min(cust_arr_time[i], cust_arr_time[j]);

                    bool feasible = false;
                    int min_veh_dist = time_ub;
                    for (const int& k : vehicles) {
                        if (dem <= capacity - veh_occ[k] && veh_dist[k][i] + dist <= std::min(veh_arr_time[k], min_time)) {
                            min_veh_dist = std::min(min_veh_dist, veh_dist[k][i]);
                            feasible = true;
                        }
                    }
                    //int dist = cust_dist[i][j] + cust_dist[j][0];
                    double max_time = fmax(cust_arr_time[i], cust_arr_time[j]);
                    if (!feasible) {
                        for (int k : vehicles) {
                            int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                            indices_removed.insert(idx);
                            idx = i * no_cust_d + j;
                            flow_ind_veh[idx].insert(k);
                        }
                    }
                    else {
                        for (int k : vehicles) {
                            int new_dist = dist + veh_dist[k][i];
                            if (new_dist > fmax(veh_arr_time[k], max_time)) {
                                int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                                indices_removed.insert(idx);

                                idx = i * no_cust_d + j;
                                flow_ind_veh[idx].insert(k);
                            }
                        }
                    }
                }
            }
        }

        // Remove Flow variables
        for (const int& idx : flow_indices) {
            if (flow_ind_veh[idx].size() == no_veh) {
                flow_removed.insert(flow_removed.end(), idx);
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
        std::cout << flow_removed.size() << std::endl;
        std::cout << "Variables Removed: " << indices_removed.size() + flow_removed.size() << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Creates the GP formulation
    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);

        // Add Variables
        no_veh_vars = no_veh * (no_cust + 1);
        no_total_vars = no_veh * (no_cust + 1) + no_veh * ((no_cust + 1) * (no_cust + 1) - no_cust - 1);
        if (var_indices.size() == 0) {
            for (int idx = 1; idx <= no_total_vars; idx++) {
                if (idx > no_veh_vars) {
                    int k = floor((idx - 1 - no_veh * no_cust_d) / (pow(no_cust_d, 2) - no_cust - 1));
                    int i = floor(
                        (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1)) / no_cust_d
                        + 1);
                    int j = (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1) - (i - 1) * no_cust_d);
                    if (i != j) var_indices.insert(idx);
                }
                else var_indices.insert(idx);
            }
        }
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
            flow[idx] = model->addVar(0, max_cap - demand[i], 0, GRB_CONTINUOUS, name);
        }

        // rest of variables
        vars.resize(no_total_vars + 1);
        for (int idx : var_indices_subset) {
            if (idx <= no_veh_vars) {
                int k = floor((idx - 1) / (no_cust + 1));
                int j = idx - k * (no_cust + 1) - 1;
                obj = veh_dist[k][j];
                name = "gamma[" + itos(k) + "][" + itos(j) + "]";
                vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, name);
            }
            else {
                int k = floor((idx - 1 - no_veh * no_cust_d) / (pow(no_cust_d, 2) - no_cust - 1));
                int i = floor(
                    (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1)) / no_cust_d
                    + 1);
                int j = (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1) - (i - 1) * no_cust_d);
                if (i != j) {
                    obj = cust_dist[i][j] - demand[i] * cust_dist[i][0];
                    name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                    vars[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
                }
            }
        }

        // Add Constraints
        GRBLinExpr lexpr, lexpr2, lexpr3, lexpr4;

        // 1.2
        for (int k : vehicles) {
            for (int j : dest_and_cust) {
                int idx = k * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // 1.3
        for (int k : vehicles) {
            for (int i : customers) {
                int idx1 = k * no_cust_d + i + 1;
                if (var_indices_subset.find(idx1) != var_indices_subset.end()) lexpr += vars[idx1];
                int idx2 = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + 0 + 1;
                if (var_indices_subset.find(idx2) != var_indices_subset.end()) lexpr2 += vars[idx2];
            }
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }

        // 1.4
        for (int k : vehicles) {
            for (int j : customers) {
                int idx = k * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += vars[idx];
            }
            for (int i : customers) {
                for (int j : dest_and_cust) {
                    if (i != j) {
                        int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx] * demand[i];
                    }
                }
            }
            model->addConstr(lexpr <= (capacity - veh_occ[k]) * lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }

        //1.5
        for (int k : vehicles) {
            for (int j : dest_and_cust) {
                int idx = k * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += veh_occ[k] * vars[idx];
            }
            model->addConstr(veh_occ[k] <= lexpr);
            lexpr.clear();
        }

        // 1.6
        for (int i : customers) {
            for (int j : dest_and_cust) {
                if (i != j) {
                    for (int k : vehicles) {
                        int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                    }
                }
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // 1.7
        for (int k : vehicles) {
            for (int j : customers) {
                int idx = k * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                for (int i : customers) {
                    if (i != j) {
                        idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                    }
                }
                for (int i : dest_and_cust) {
                    if (i != j) {
                        idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + i + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += vars[idx];
                    }
                }
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }

        // 1.8
        for (int k : vehicles) {
            // Vehicle edges
            for (int j : dest_and_cust) {
                int idx = k * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += (veh_dist[k][j] - veh_arr_time[k]) * vars[idx];
            }
            // All edges
            for (int l : customers) {
                for (int j : dest_and_cust) {
                    if (l != j) {
                        int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (l - 1) * no_cust_d + j + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += cust_dist[l][j] * vars[idx];
                    }
                }
            }
            for (int i : customers) {
                if (cust_arr_time[i] < veh_arr_time[k]) {
                    double coeff = veh_arr_time[k] - cust_arr_time[i];
                    // Customer Edges
                    for (int j : dest_and_cust) {
                        int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr3 += coeff * vars[idx];
                    }
                }
                model->addConstr(lexpr + lexpr2 + lexpr3 <= 0);
                lexpr3.clear();
            }
            lexpr.clear();
            lexpr2.clear();
        }

        model->update();
        for (const int& i : customers) {
            for (const int& k : vehicles) {
                int idx = k * no_cust_d + i + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += (capacity -  veh_occ[k]) * vars[idx];
            }

            for (const int& j : customers) {
                if (i == j) continue;
                int idx = j * no_cust_d + i;
                if (flow_indices.find(idx) != flow_indices.end()) lexpr2 += flow[idx];
            }

            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                int idx = i * no_cust_d + j;
                if (flow_indices.find(idx) != flow_indices.end()) lexpr3 += flow[idx];
            }

            for (const int& k : vehicles) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr4 += demand[i] * vars[idx];
                }
            }

            model->addConstr(lexpr + lexpr2 == lexpr3 + lexpr4);
            lexpr.clear();
            lexpr2.clear();
            lexpr3.clear();
            lexpr4.clear();
        }

        for (const int& idx : flow_indices) {
            const int i = floor(idx / float(no_cust_d));
            const int j = idx - i * no_cust_d;
            for (const int& k : vehicles) {
                int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) {
                    lexpr += (capacity - veh_occ[k] - demand[i]) * vars[idx];
                    if (j != 0) lexpr2 += demand[j] * vars[idx];
                }
            }
            model->addConstr(flow[i * no_cust_d + j] <= lexpr);
            if (j != 0) model->addConstr(lexpr2 <= flow[i * no_cust_d + j]);
            lexpr.clear();
            lexpr2.clear();
        }

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