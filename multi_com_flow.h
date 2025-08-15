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


class MCF {
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

    std::vector<GRBVar> f_flow;
    std::vector<GRBVar> g_flow;
    std::vector<GRBVar> vars;
    int no_base_constraints = 0;

    //Valid inequalities
    std::vector<std::set<int>> rounded_cap;
    std::set<int> active_rc_ineqs;
    int total_rcc = 0;


    // Edges
    std::set<int> var_indices;
    std::set<int> var_indices_subset;
    std::set<int> f_flow_indices;
    std::set<int> g_flow_indices;

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
    MCF(const std::string& filename) {
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

        for (const int& l : customers) {
            for (const int& i : dest_and_cust) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                    f_flow_indices.insert(f_flow_indices.end(), idx);
                    g_flow_indices.insert(g_flow_indices.end(), idx);
                }
            }
        }

        rounded_cap.resize(no_cust * (no_cust - 1) / 2);

        // Setup Constraint Managers
        CMGR_CreateCMgr(&MyCutsCMP, 100);
        CMGR_CreateCMgr(&MyOldCutsCMP, 100);
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
        std::set<int> f_flow_indices_removed;
        std::set<int> g_flow_indices_removed;

        // Remove flow_variables fixed at 0
        for (const int& i : customers) {
            for (const int& j : dest_and_cust) {
                int idx = i * pow(no_cust_d, 2) + j * no_cust_d + 0;
                f_flow_indices_removed.insert(idx);
                idx = i * pow(no_cust_d, 2) + i * no_cust_d + j;
                f_flow_indices_removed.insert(idx);
                idx = i * pow(no_cust_d, 2) + j * no_cust_d + i;
                g_flow_indices_removed.insert(idx);
                idx = i * pow(no_cust_d, 2) + 0 * no_cust_d + j;
                g_flow_indices_removed.insert(idx);
            }
        }


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

                    for (int i : customers) {
                        if (i != j) {
                            idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                            indices_removed.insert(idx);

                            idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (j - 1) * no_cust_d + i + 1;
                            indices_removed.insert(idx);
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
                        for (const int& k : vehicles) {
                            int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                            indices_removed.insert(idx);
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            int new_dist = dist + veh_dist[k][i];
                            if (new_dist > fmax(veh_arr_time[k], max_time)) {
                                int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                                indices_removed.insert(idx);
                            }
                        }
                    }
                }
            }
        }


        // Remove direct flow edges
        for (const int& i : dest_and_cust) {
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                bool exists = false;
                for (const int& k : vehicles) {
                    int idx;
                    if (i == 0) idx = k * no_cust_d + j + 1;
                    else idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                    if (std::find(indices_removed.begin(), indices_removed.end(), idx) == indices_removed.end()) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) {
                    for (const int& l : customers) {
                        int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        f_flow_indices_removed.insert(idx);
                        g_flow_indices_removed.insert(idx);
                    }
                }
            }
        }

        // Remove indirect flow edges
        for (const int& l : customers) {
            for (const int& i : dest_and_cust) {
                bool exists = false;
                for (const int& k : vehicles) {
                    int idx;
                    if (i == 0) idx = k * no_cust_d + l + 1;
                    else idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + l + 1;
                    if (std::find(indices_removed.begin(), indices_removed.end(), idx) == indices_removed.end()) {
                        exists = true;
                        break;
                    }
                }
                if (!exists) {
                    for (const int& j : dest_and_cust) {
                        int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        f_flow_indices_removed.insert(idx);
                        idx = l * pow(no_cust_d, 2) + j * no_cust_d + i;
                        f_flow_indices_removed.insert(idx);
                        idx = i * pow(no_cust_d, 2) + j * no_cust_d + l;
                        g_flow_indices_removed.insert(idx);
                        idx = i* pow(no_cust_d, 2) + l * no_cust_d + j;
                        g_flow_indices_removed.insert(idx);
                    }
                }
            }

            for (const int& i : dest_and_cust) {
                if (i == l) continue;
                for (const int& j : customers) {
                    if (j == l || j == i) continue;
                    if (demand[l] + demand[i] + demand[j] > max_cap) {
                        int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        f_flow_indices_removed.insert(idx);
                        g_flow_indices_removed.insert(idx);
                    }
                }
            }
        }
        int count = f_flow_indices.size() + g_flow_indices.size();
        std::set<int> temp_set;
        std::set_difference(var_indices.begin(), var_indices.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices = temp_set;
        temp_set.clear();
        std::set_difference(var_indices_subset.begin(), var_indices_subset.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices_subset = temp_set;
        temp_set.clear();
        std::set_difference(f_flow_indices.begin(), f_flow_indices.end(), f_flow_indices_removed.begin(), f_flow_indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        f_flow_indices = temp_set;
        temp_set.clear();
        std::set_difference(g_flow_indices.begin(), g_flow_indices.end(), g_flow_indices_removed.begin(), g_flow_indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        g_flow_indices = temp_set;
        count = count - f_flow_indices.size() - g_flow_indices.size();
        std::cout << count << std::endl;
        std::cout << "Variables Removed: " << indices_removed.size() << std::endl;
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
        if (var_indices_subset.size() == 0) var_indices_subset = var_indices;

        double obj = 0;
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        std::string name;

        // flow variables
        f_flow.resize(pow(no_cust_d, 3));
        g_flow.resize(pow(no_cust_d, 3));

        for (const int& idx : f_flow_indices) {
            int l = floor(idx / (pow(no_cust_d, 2)));
            int i = floor((idx - l * pow(no_cust_d, 2)) / no_cust_d);
            int j = idx - l * pow(no_cust_d, 2) - i * no_cust_d;
            name = "f[" + itos(l) + "][" + itos(i) + "][" + itos(j) + "]";
            f_flow[idx] = model->addVar(0, 1, 0, GRB_CONTINUOUS, name);
        }
        for (const int& idx : g_flow_indices) {
            int l = floor(idx / (pow(no_cust_d, 2)));
            int i = floor((idx - l * pow(no_cust_d, 2)) / no_cust_d);
            int j = idx - l * pow(no_cust_d, 2) - i * no_cust_d;
            name = "g[" + itos(l) + "][" + itos(i) + "][" + itos(j) + "]";
            g_flow[idx] = model->addVar(0, 1, 0, GRB_CONTINUOUS, name);
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

        std::vector<std::vector<GRBLinExpr>> f_out(no_cust_d);
        std::vector<std::vector<GRBLinExpr>> f_in(no_cust_d);
        std::vector<std::vector<GRBLinExpr>> g_out(no_cust_d);
        std::vector<std::vector<GRBLinExpr>> g_in(no_cust_d);

        for (const int& l : customers) {
            f_out[l].resize(no_cust_d);
            f_in[l].resize(no_cust_d);
            g_out[l].resize(no_cust_d);
            g_in[l].resize(no_cust_d);
            for (const int& i : dest_and_cust) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    int out_index = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                    int in_index = l * pow(no_cust_d, 2) + j * no_cust_d + i;
                    if (f_flow_indices.find(out_index) != f_flow_indices.end()) f_out[l][i] += f_flow[out_index];
                    if (f_flow_indices.find(in_index) != f_flow_indices.end())f_in[l][i] += f_flow[in_index];
                    if (g_flow_indices.find(out_index) != g_flow_indices.end())g_out[l][i] += g_flow[out_index];
                    if (g_flow_indices.find(in_index) != g_flow_indices.end())g_in[l][i] += g_flow[in_index];
                }
            }
        }

        for (const int& i : customers) {
            for (int j : dest_and_cust) {
                for (const int& k : vehicles) {
                    int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                }
            }
            model->addConstr(f_out[i][0] == f_in[i][i]);
            model->addConstr(f_in[i][i] == g_out[i][i]);
            model->addConstr(g_out[i][i] == g_in[i][0]);
            model->addConstr(g_in[i][0] == lexpr);
            lexpr.clear();
            for (const int& j : customers) {
                if (i == j) continue;
                model->addConstr(f_out[i][j] == f_in[i][j]);
                model->addConstr(f_in[i][j] == g_out[j][i]);
                model->addConstr(g_out[j][i] == g_in[j][i]);
            }
        }

        //model->update();
        for (const int& i : dest_and_cust) {
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                bool exists = false;
                for (const int& k : vehicles) {
                    int idx;
                    if (i > 0) idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                    else idx = k * no_cust_d + j + 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) {
                        lexpr2 += vars[idx];
                        exists = true;
                    }
                }
                if (!exists) continue;
                for (const int& l : customers) {
                    int flow_idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                    if (f_flow_indices.find(flow_idx) == f_flow_indices.end() && g_flow_indices.find(flow_idx) == g_flow_indices.end()) continue;
                    if (f_flow_indices.find(flow_idx) != f_flow_indices.end()) lexpr += f_flow[flow_idx];
                    if (g_flow_indices.find(flow_idx) != g_flow_indices.end()) lexpr += g_flow[flow_idx];

                    model->addConstr(lexpr <= lexpr2);
                    lexpr.clear();
                }
                lexpr2.clear();
            }
        }

        for (const int& i : dest_and_cust) {
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                bool exists = false;
                for (const int& k : vehicles) {
                    int idx;
                    if (i > 0) idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                    else idx = k * no_cust_d + j + 1;
                    if (var_indices_subset.find(idx) != var_indices_subset.end()) {
                        lexpr2 += (capacity - veh_occ[k] - demand[i] - demand[j]) * vars[idx];
                        exists = true;
                    }
                }
                if (!exists) continue;
                for (const int& l : customers) {
                    if (l == i || l == j) continue;
                    int flow_idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                    if (f_flow_indices.find(flow_idx) != f_flow_indices.end()) lexpr += demand[l] * f_flow[flow_idx];
                    if (g_flow_indices.find(flow_idx) != g_flow_indices.end()) lexpr += demand[l] * g_flow[flow_idx];
                }
                model->addConstr(lexpr <= lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
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