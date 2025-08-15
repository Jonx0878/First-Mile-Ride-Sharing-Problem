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


class LB {
private:
    void calculate_reduced_cost_gamma(
        double& rc,
        std::vector<GRBConstr>& col_constrs,
        std::vector<double>& col_coeffs,
        const int& k,
        const int& j,
        GRBConstr* constrs,
        double* dual_values
    ) {
        int idx;
        double coeff;
        // 2.2
        idx = k;
        rc -= dual_values[idx];
        col_constrs.push_back(constrs[idx]);
        col_coeffs.push_back(1);

        // 2.3
        if (j == 0) {
            idx = no_veh + k;
            coeff = capacity - veh_occ[k];
            rc -= coeff * dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(coeff);
        }

        // 2.7
        if (j > 0) {
            idx = 4 * no_veh + (1 + k) * no_cust + j - 1;
            rc -= dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(1);
        }

        // 2.8
        for (const int& i : customers) {
            idx = 4 * no_veh + (1 + no_veh + k) * no_cust + i - 1;
            coeff = veh_dist[k][j] - veh_arr_time[k];
            rc -= coeff * dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(coeff);
        }
        // 2.14 - 2.15
        if (j > 0) {
            idx = 4 * no_veh + (1 + 2 * k + 2 * no_veh) * no_cust + (4 * no_veh + 1) * (no_cust_d * (no_cust_d - 1)) / 2 + 2 * j - 2;
            rc -= dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(1);
            idx++;
            rc -= dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(1);
        }
    }

    void calculate_reduced_cost_chi(
        double& rc,
        std::vector<GRBConstr>& col_constrs,
        std::vector<double>& col_coeffs,
        const int& k,
        const int& i,
        const int& j,
        GRBConstr* constrs,
        double* dual_values
    ) {
        int idx;
        double coeff;
        int sym_idx = (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;

        // 2.5 
        idx = 3 * no_veh + k;
        rc -= dual_values[idx];
        col_constrs.push_back(constrs[idx]);
        col_coeffs.push_back(1);


        // 2.7
        for (const int& l : customers) {
            if (i == l || j == l) {
                idx = 4 * no_veh + (1 + k) * no_cust + l - 1;
                rc -= dual_values[idx];
                col_constrs.push_back(constrs[idx]);
                col_coeffs.push_back(1);
            }
        }

        // 2.8
        for (const int& l : customers) {
            idx = 4 * no_veh + (1 + k + no_veh) * no_cust + l - 1;
            coeff = cust_dist[i][j];
            rc -= coeff * dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(coeff);
        }

        // 2.9
        idx = 4 * no_veh + (1 + 2 * no_veh) * no_cust + sym_idx;
        rc -= -dual_values[idx];
        col_constrs.push_back(constrs[idx]);
        col_coeffs.push_back(-1);

        // 2.10 - 2.13
        idx = 4 * no_veh + (1 + 2 * no_veh) * no_cust + (4 * k + 1) * (no_cust_d * (no_cust_d - 1)) / 2 + 4 * sym_idx - 1;
        for (int l = 0; l < 4; l++) {
            idx++;
            coeff = (l == 0) ? -1 : 1;
            rc -= coeff * dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(coeff);
        }
    }

    void calculate_reduced_cost_psi(
        double& rc,
        std::vector<GRBConstr>& col_constrs,
        std::vector<double>& col_coeffs,
        const int& i,
        const int& j,
        GRBConstr* constrs,
        double* dual_values,
        std::unordered_map<int, int>& rc_indices
    ) {
        int idx;
        double coeff;
        int sym_idx = (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;

        // 2.9
        idx = 4 * no_veh + (1 + 2 * no_veh) * no_cust + sym_idx;
        rc -= dual_values[idx];
        col_constrs.push_back(constrs[idx]);
        col_coeffs.push_back(1);


        // 2.10 + 2.13
        for (const int& k : vehicles){
            idx = 4 * no_veh + (1 + 2 * no_veh) * no_cust + (4 * k + 1) * (no_cust_d * (no_cust_d - 1)) / 2 + 4 * sym_idx;
            rc -= dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(1);
            idx += 3;
            rc -= -dual_values[idx];
            col_constrs.push_back(constrs[idx]);
            col_coeffs.push_back(-1);
        }

        // RCC
        if (i > 0) {
            idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (j - 1) - 1;
            for (int rc_idx : rounded_cap[idx]) {
                if (rc_indices.count(rc_idx) == 1) {
                    int constr_idx = rc_indices[rc_idx];
                    rc -= dual_values[constr_idx];
                    col_constrs.push_back(constrs[constr_idx]);
                    col_coeffs.push_back(1);
                }
            }
        }
    }

    void convert_to_cvrp_solution(
        std::vector<int>& cvrp_to_fmrsp,
        int& cvrp_no_cust,
        std::set<int>& cvrp_customers,
        std::vector<int>& cvrp_demand,
        std::set<std::tuple<int, int, double>>& cvrp_solution
    ) {
        std::vector<int> fmrsp_to_cvrp;
        fmrsp_to_cvrp.resize(customers.size() + 1);
        cvrp_to_fmrsp.emplace_back(0);
        cvrp_demand.emplace_back(0);
        std::vector<int> mapping_cvrp_fmrsp;
        std::vector<double> customer_sum(no_cust_d);
        for (int i : dest_and_cust) {
            customer_sum[i] = 0;
        }

        cvrp_no_cust = 0;

        std::vector<std::vector<double>> cust_edges(no_cust_d);
        for (int i : customers) {
            cust_edges[i].resize(no_cust_d);
            for (int j : customers) {
                cust_edges[i][j] = 0;
            }
        }
        std::vector<double> cust_edges_missing(no_cust_d);
        for (int i : customers) {
            cust_edges_missing[i] = 2;
        }

        for (const int& idx : var_indices_subset) {
            if (idx < no_dir_veh_vars) {
                int k = floor(idx / no_cust_d);
                int j = idx - k * no_cust_d;
                if (j > 0) customer_sum[j] += vars[idx].get(GRB_DoubleAttr_X);
            }
            else if (idx >= no_total_veh_vars) {
                int sym_idx = idx - no_total_veh_vars;
                int i = no_cust_d - 2 - int(floor(sqrt(-8 * sym_idx + 4 * no_cust_d * (no_cust_d - 1) - 7) / 2.0 - 0.5));
                int j = sym_idx + i * (i + 3) / 2 - no_cust_d * i + 1;
                double val = vars[idx].get(GRB_DoubleAttr_X);
                if (i > 0) customer_sum[i] += val;
                customer_sum[j] += val;
                if (i > 0) cust_edges[i][j] += val;
            }
        }

        for (const int& i : customers) {
            if (customer_sum[i] >= 0.001) {
                cvrp_no_cust++;
                cvrp_to_fmrsp.emplace_back(i);
                cvrp_demand.emplace_back(demand[i]);
                cvrp_customers.insert(cvrp_no_cust);
                fmrsp_to_cvrp[i] = cvrp_no_cust;
            }
            else fmrsp_to_cvrp[i] = 0;
        }

        for (const int& i : customers) {
            for (const int& j : customers) {
                if (i < j) cust_edges_missing[i] -= cust_edges[i][j];
                else if (i > j) cust_edges_missing[i] -= cust_edges[j][i];
            }
        }

        for (const int& i : customers) {
            for (int j = i + 1; j < no_cust_d; j++) {
                if (cust_edges[i][j] >= 0.001) {
                    cvrp_solution.insert({ fmrsp_to_cvrp[i], fmrsp_to_cvrp[j], cust_edges[i][j] });
                }
            }
        }

        for (int i : cvrp_customers) {
            if (0.001 <= cust_edges_missing[cvrp_to_fmrsp[i]]) {
                cvrp_solution.insert({ 0, i, cust_edges_missing[cvrp_to_fmrsp[i]] });
            }
        }
        //std::cout << cvrp_solution.size() << std::endl;
        //for (auto& edge : cvrp_solution) {
        //	std::cout << std::fixed << std::setprecision(3) << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
        //}
    }

    void add_rc_inequalities(
        int& no_new_cuts,
        const std::vector<RC_Inequality>& new_cuts
    ) {
        for (RC_Inequality cut : new_cuts) {
            GRBLinExpr lexpr;
            std::vector<int> indices;
            for (const int& i : cut.set) {
                for (const int& j : cut.set) {
                    if (i < j) {
                        int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                        int sym_idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (j - 1) - 1;
                        indices.emplace_back(sym_idx);
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                    }
                }
            }
            if (cut.check_violation) {
                if (lexpr.getValue() >= cut.rhs + 0.001) {
                    goto addcut;
                }
            }
            else {
            addcut:
                model->addConstr(lexpr <= cut.rhs, "RC_" + itos(total_rcc));
                for (const int& idx : indices) {
                    rounded_cap[idx].insert(total_rcc);
                }
                no_new_cuts++;
                total_rcc++;
            }
        }
    }

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

    // Column Generation
    int outer_iterations = 0;

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
    LB(const std::string& filename) {
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

        // Setup Constraint Managers
        CMGR_CreateCMgr(&MyCutsCMP, 100);
        CMGR_CreateCMgr(&MyOldCutsCMP, 100);
    }

    void set_starting_edges() {
        int no_edges = 10;

        std::set<std::pair<int, int>> veh_loc;

        for (const int& k : vehicles) {
            if (veh_occ[k] > 0) {
                var_indices_subset.insert(k * no_cust_d);
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
                int idx = k * no_cust_d + j;
                if (var_indices.find(idx) != var_indices.end()) var_indices_subset.insert(idx);
            }
        }

        for (const int& i : dest_and_cust) {
            std::vector<DistanceLocationPair> distances;
            for (const int& j : customers) {
                if (i == j) continue;
                distances.push_back(DistanceLocationPair(cust_dist[i][j], j));
            }
            // Sort distances in ascending order
            std::sort(distances.begin(), distances.end(), [](const DistanceLocationPair& a, const DistanceLocationPair& b) {
                return a.distance < b.distance;
                });
            for (int l = 0; l < fmin(no_edges, distances.size()); l++) {
                int j = distances[l].location;
                int idx;
                for (const int& k : vehicles) {
                    if (i < j) {
                        idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    }
                    else {
                        idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
                    }
                    if (var_indices.find(idx) != var_indices.end()) var_indices_subset.insert(idx);
                }
            }
        }

        // Add psi variables
        for (const int& i : dest_and_cust) {
            for (int j = i + 1; j < no_cust_d; j++) {
                bool exists = false;
                for (const int& k : vehicles) {
                    int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices.find(idx) != var_indices.end()) {
                        exists = true;
                        break;
                    }
                }
                if (exists) {
                    int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    var_indices_subset.insert(var_indices_subset.end(), idx);
                }
            }
        }

    }

    void preprocess() {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Potential # Edge Variables: " << var_indices.size() << std::endl;
        std::cout << "Preprocesing..." << std::endl;
        std::set<int> indices_removed;

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
            for (int j = i + 1; j < no_cust_d; j++){
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
                }
            }
        }

        std::set<int> temp_set;
        std::set_difference(var_indices.begin(), var_indices.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices = temp_set;
        temp_set.clear();
        std::set_difference(var_indices_subset.begin(), var_indices_subset.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices_subset = temp_set;
        std::cout << "Variables Removed: " << indices_removed.size()  << " # Variables left: " << var_indices.size() << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Creates the JLS formulation
    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);

        // Add Variables
        double obj;
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        std::string name;

        if (var_indices_subset.size() == 0) var_indices_subset = var_indices;
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
                vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, name);
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
        GRBLinExpr lexpr, lexpr2, lexpr3;



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
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }

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
        for (const int& k : vehicles) {
            for (const int& i : customers) {
                int veh_idx = k * no_cust_d + i;
                if (std::find(var_indices.begin(), var_indices.end(), veh_idx) == var_indices.end()) model->addConstr(lexpr == 0); // Removing constr may be difficult for Column Generation
                else {
                    if (var_indices_subset.find(veh_idx) != var_indices_subset.end()) lexpr += vars[veh_idx];
                    for (const int& j : dest_and_cust) {
                        if (i < j) {
                            int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                            if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                        }
                        else if (i > j) {
                            int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + i - j - 1;
                            if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                        }
                    }
                    model->addConstr(lexpr == 2 * xi[i][k]);
                    lexpr.clear();
                }
            }
        }

        // 2.8
        // Time Constraints
        for (const int& k : vehicles) {
            // Vehicle edges
            for (const int& j : dest_and_cust) {
                int idx = k * no_cust_d + j;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += (veh_dist[k][j] - veh_arr_time[k]) * vars[idx];
            }
            // All edges
            for (const int& l : dest_and_cust) {
                for (int j = l + 1; j < no_cust_d; j++) {
                    int idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - l) * ((no_cust_d - l) - 1) / 2 + j - l - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += cust_dist[l][j] * vars[idx];
                }
            }
            for (const int& i : customers) {
                if (cust_arr_time[i] < veh_arr_time[k]) {
                    double coeff = veh_arr_time[k] - cust_arr_time[i];
                    int idx = k * no_cust_d + i;
                    if (var_indices.find(idx) != var_indices.end()) lexpr3 += coeff * xi[i][k];
                }
                model->addConstr(lexpr + lexpr2 + lexpr3 <= 0);
                lexpr3.clear();
            }
            lexpr.clear();
            lexpr2.clear();
        }

        // 2.9
        // x constraints
        for (const int& i : dest_and_cust) {
            for (int j = i + 1; j < no_cust_d; j++) {
                int idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                for (const int& k : vehicles) {
                    idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += vars[idx];
                }
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }

        // 2.10 - 2.13
        for (const int& k : vehicles) {
            for (const int& i : dest_and_cust) {
                int idx = k * no_cust_d + i;
                if (var_indices.find(idx) == var_indices.end()) {
                    for (int j = i + 1; j < no_cust_d; j++) {
                        model->addConstr(-2 <= lexpr);
                        if (i == 0) model->addConstr(lexpr <= xi[0][k]);
                        else model->addConstr(lexpr <= 0);
                        model->addConstr(lexpr <= 0);
                        model->addConstr(lexpr <= 0);
                    }
                    continue;
                }
                for (int j = i + 1; j < no_cust_d; j++) {
                    // chi var
                    idx = no_dir_veh_vars + (k + 1) * (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices.find(idx) == var_indices.end()) {
                        model->addConstr(-2 <= lexpr);
                        if (i == 0) model->addConstr(lexpr <= xi[0][k]);
                        else model->addConstr(lexpr <= 0);
                        model->addConstr(lexpr <= 0);
                        model->addConstr(lexpr <= 0);
                        continue;
                    }
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                    //gamma_var
                    idx = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr2 += vars[idx];

                    model->addConstr(xi[i][k] + xi[j][k] + lexpr2 - 2 <= lexpr);
                    model->addConstr(lexpr <= xi[i][k]);
                    model->addConstr(lexpr <= xi[j][k]);
                    model->addConstr(lexpr <= lexpr2);
                    lexpr.clear();
                    lexpr2.clear();
                    lexpr3.clear();
                }
            }
        }

        // 2.14 - 2.15
        for (const int& k : vehicles) {
            for (const int& j : customers) {
                int idx = k * no_cust_d + j;
                if (var_indices.find(idx) == var_indices.end()) {
                    model->addConstr(lexpr == 0);
                    model->addConstr(lexpr == 0);
                    continue;
                }
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                model->addConstr(lexpr <= xi[j][k]);
                //model->addConstr(lexpr <= xi[0][k]);
                lexpr.clear();
            }
        }
        model->update();
        no_base_constraints = model->get(GRB_IntAttr_NumConstrs);

        // 2.16
        // RCC for |S| == 3 - CG Currently requires RCC - Change later
        for (int i : customers) {
            for (int j = i + 1; j <= no_cust; j++) {
                int idx1 = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + j - i - 1;
                for (int l = j + 1; l <= no_cust; l++) {
                    int idx2 = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + l - j - 1;
                    int idx3 = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + l - i - 1;
                    if (capacity == 4) {
                        if (var_indices_subset.find(idx1) != var_indices_subset.end()) lexpr += vars[idx1];
                        if (var_indices_subset.find(idx2) != var_indices_subset.end()) lexpr += vars[idx2];
                        if (var_indices_subset.find(idx3) != var_indices_subset.end()) lexpr += vars[idx3];
                        model->addConstr(lexpr <= 3 - ceil((demand[i] + demand[j] + demand[l]) / float(max_cap)), "RC_" + itos(total_rcc));
                        int sym_idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (j - 1) - 1;
                        rounded_cap[sym_idx].insert(total_rcc);
                        sym_idx = (2 * no_cust - 3 - (j - 1)) * (j - 1) / 2.0 + (l - 1) - 1;
                        rounded_cap[sym_idx].insert(total_rcc);
                        sym_idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (l - 1) - 1;
                        rounded_cap[sym_idx].insert(total_rcc);
                        lexpr.clear();
                        total_rcc++;
                    }
                    else {
                        for (int m = l + 1; m <= no_cust; m++) {
                            int idx4 = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - i) * ((no_cust_d - i) - 1) / 2 + m - i - 1;
                            int idx5 = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - j) * ((no_cust_d - j) - 1) / 2 + m - j - 1;
                            int idx6 = no_total_veh_vars + (no_cust_d * (no_cust_d - 1) / 2) - (no_cust_d - l) * ((no_cust_d - l) - 1) / 2 + m - l - 1;
                            if (var_indices_subset.find(idx1) != var_indices_subset.end()) lexpr += vars[idx1];
                            if (var_indices_subset.find(idx2) != var_indices_subset.end()) lexpr += vars[idx2];
                            if (var_indices_subset.find(idx3) != var_indices_subset.end()) lexpr += vars[idx3];
                            if (var_indices_subset.find(idx4) != var_indices_subset.end()) lexpr += vars[idx4];
                            if (var_indices_subset.find(idx5) != var_indices_subset.end()) lexpr += vars[idx5];
                            if (var_indices_subset.find(idx6) != var_indices_subset.end()) lexpr += vars[idx6];
                            int sym_idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (j - 1) - 1;
                            rounded_cap[sym_idx].insert(total_rcc);
                            sym_idx = (2 * no_cust - 3 - (j - 1)) * (j - 1) / 2.0 + (l - 1) - 1;
                            rounded_cap[sym_idx].insert(total_rcc);
                            sym_idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (l - 1) - 1;
                            rounded_cap[sym_idx].insert(total_rcc);
                            sym_idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (m - 1) - 1;
                            rounded_cap[sym_idx].insert(total_rcc);
                            sym_idx = (2 * no_cust - 3 - (j - 1)) * (j - 1) / 2.0 + (m - 1) - 1;
                            rounded_cap[sym_idx].insert(total_rcc);
                            sym_idx = (2 * no_cust - 3 - (l - 1)) * (l - 1) / 2.0 + (m - 1) - 1;
                            rounded_cap[sym_idx].insert(total_rcc);
                            lexpr.clear();
                            total_rcc++;
                        }
                    }
                }
            }
        }
        auto stop = std::chrono::high_resolution_clock::now();
        model_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void optimise_model() {
        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        auto stop = std::chrono::high_resolution_clock::now();
        solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    /*------------------------------------------------------------------------------------------*/
    /*									Column Generation										*/
    /*------------------------------------------------------------------------------------------*/

    void add_variables(
        bool& new_vars,
        bool add_all = false
    ) {
        auto start = std::chrono::high_resolution_clock::now();
        int no_new_vars = 0;
        std::unordered_map<int, int> rc_indices;

        outer_iterations++;
        std::set<int> var_indices_unused;
        std::set_difference(var_indices.begin(), var_indices.end(), var_indices_subset.begin(), var_indices_subset.end(), std::inserter(var_indices_unused, var_indices_unused.begin()));
        GRBConstr* constrs = model->getConstrs();
        int no_constr = model->get(GRB_IntAttr_NumConstrs);
        double* dual_values = model->get(GRB_DoubleAttr_Pi, constrs, no_constr);

        std::string* constr_names = model->get(GRB_StringAttr_ConstrName, constrs, no_constr);

        process_valid_inequalties(rc_indices, constr_names, no_base_constraints, no_constr);
        std::vector<int> new_vars_veh(no_veh + 1);
        for (const int& k : vehicles) {
            new_vars_veh[k] = 0;
        }
        new_vars_veh[no_veh] = 0;

        for (int idx : var_indices_unused) {
            double obj;
            double rc;
            int k = - 1, i, j;
            std::string name;
            std::vector<GRBConstr> col_constrs;
            std::vector<double> col_coeffs;
            if (idx < no_dir_veh_vars) {
                k = floor(idx / no_cust_d);
                //if (new_vars_veh[k] >= no_cust_d) continue;
                j = idx - k * no_cust_d;
                rc = veh_dist[k][j];
                obj = rc;
                name = "gamma[" + itos(k) + "][" + itos(j) + "]";
                calculate_reduced_cost_gamma(rc, col_constrs, col_coeffs, k, j, constrs, dual_values);
            }
            else if (idx < no_total_veh_vars) {
                k = floor((idx - no_dir_veh_vars) / (no_cust_d * (no_cust_d - 1) / 2));
                //if (new_vars_veh[k] >= no_cust_d) continue;
                int sym_idx = idx - no_dir_veh_vars - k * (no_cust_d * (no_cust_d - 1) / 2);
                i = no_cust_d - 2 - int(floor(sqrt(-8 * sym_idx + 4 * no_cust_d * (no_cust_d - 1) - 7) / 2.0 - 0.5));
                j = sym_idx + i * (i + 3) / 2 - no_cust_d * i + 1;
                rc = 0;
                obj = rc;
                name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                calculate_reduced_cost_chi(rc, col_constrs, col_coeffs, k, i, j, constrs, dual_values);
            }
            else {
                //if (new_vars_veh[no_veh] >= no_cust_d) continue;
                int sym_idx = idx - no_total_veh_vars;
                i = no_cust_d - 2 - int(floor(sqrt(-8 * sym_idx + 4 * no_cust_d * (no_cust_d - 1) - 7) / 2.0 - 0.5));
                j = sym_idx + i * (i + 3) / 2 - no_cust_d * i + 1;
                rc = cust_dist[i][j];
                obj = rc;
                name = "psi[" + itos(i) + "][" + itos(j) + "]";
                calculate_reduced_cost_psi(rc, col_constrs, col_coeffs, i, j, constrs, dual_values, rc_indices);
            }

            //std::cout << obj << " " << rc << std::endl;
            if (rc <= -0.001 || add_all) {
                if (k >= 0) new_vars_veh[k]++;
                else new_vars_veh[no_veh]++;
                GRBColumn col;
                col.addTerms(&col_coeffs[0], &col_constrs[0], int(col_coeffs.size()));
                vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, col, name);
                var_indices_subset.insert(idx);
                no_new_vars++;
            }
            //if (no_new_vars >= no_cust * no_veh) break;
        }

        new_vars = (no_new_vars > 0);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
        add_vars_time += duration.count();
    }

    void add_all_variables() {
        bool new_vars = false;
        add_variables(new_vars, true);
    }

    void remove_constraints() { remove_constraints_from_model(model, no_base_constraints); };

    void make_integral() {
        optimise_model();
        model->getEnv().set(GRB_IntParam_OutputFlag, 1);

        for (const int& k : vehicles) {
            xi[0][k].set(GRB_CharAttr_VType, GRB_BINARY);
        }
        for (const int& idx : var_indices) {
            if (idx >= no_dir_veh_vars) break;
            // xi
            int k = floor(idx / no_cust_d);
            int i = idx - k * no_cust_d;
            xi[i][k].set(GRB_CharAttr_VType, GRB_BINARY);
        }

        for (const int& idx : var_indices_subset) {
            if (idx > no_total_veh_vars){
                vars[idx].set(GRB_CharAttr_VType, GRB_BINARY);
            }
        }
    }

    void separate_rc() {
        int iteration = 0;
        while (true) {
            iteration++;
            std::vector<int> cvrp_to_fmrsp;
            int cvrp_no_cust;
            std::vector<int> cvrp_demand;
            std::set<int> cvrp_customers;
            std::set<std::tuple<int, int, double>> cvrp_solution;
            if (iteration % 10 == 0) {
                remove_constraints();
                optimise_model();
            }


            convert_to_cvrp_solution(
                cvrp_to_fmrsp,
                cvrp_no_cust,
                cvrp_customers,
                cvrp_demand,
                cvrp_solution
            );
            double max_violation = 0;
            std::vector<RC_Inequality> new_cuts;
            double threshold = 0.001;
            if (cvrp_no_cust >= 2) {

                auto start = std::chrono::high_resolution_clock::now();

                separate_rcc_heuristically(
                    new_cuts,
                    max_violation,
                    cvrp_solution,
                    cvrp_no_cust + 1,
                    cvrp_demand,
                    capacity,
                    0.001,
                    0.001,
                    MyOldCutsCMP,
                    MyCutsCMP
                );

                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
                RC_separation_time += duration.count();
            }
            if (max_violation < threshold) break;

            if (cvrp_no_cust < no_cust) {
                convert_cuts_to_fmrsp(cvrp_to_fmrsp, new_cuts);
            }


            int no_new_cuts = 0;
            add_rc_inequalities(no_new_cuts, new_cuts);
            if (no_new_cuts == 0) break;
            optimise_model();
            if (iteration % 10 == 0) std::cout << model->get(GRB_DoubleAttr_ObjVal) << std::endl;

            for (int i = 0; i < MyCutsCMP->Size; i++)
            {
                CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP, i, 0);
            }
            MyCutsCMP->Size = 0;

        }
    }

};