#pragma once
#include <algorithm>
#include <gurobi_c++.h>
#include <set>

#include "cvrpsep/cnstrmgr.h"
#include "col_gen.h"
#include "load.h"
#include "rc.h"
#include "utils.h"


class FMRSP2 {
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
    std::vector<std::vector<std::vector<GRBVar>>> x;
    std::vector<std::vector<GRBVar>> y_veh;
    std::vector<std::vector<GRBVar>> y_cust;
    std::vector<std::vector<GRBVar>> xi;
    std::vector<GRBVar> q;


    std::vector<GRBVar> phi;
    std::vector<GRBVar> vars;
    int no_base_constraints = 0;

    //Valid inequalities
    std::vector<std::set<int>> rounded_cap;
    std::set<int> active_rc_ineqs;
    int total_rcc = 0;

    // Columns
    std::vector<std::vector<GRBColumn>> y_veh_col;
    std::vector<std::vector<GRBColumn>> y_cust_col;

    // Edges
    std::set<int> var_indices;
    std::set<int> var_indices_subset;
    std::set<std::pair<int, int>> veh_edges;
    std::set<std::pair<int, int>> veh_edges_subset;
    std::set<std::pair<int, int>> cust_edges;
    std::set<std::pair<int, int>> cust_edges_subset;
   
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
    FMRSP2(const std::string& filename) {
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

        // Setup Edges - Only Symmetrical
        for (int k : vehicles) {
            for (int i : dest_and_cust) {
                veh_edges.insert({ k, i });
            }
        }

        for (int i : dest_and_cust) {
            for (int j = i + 1; j <= no_cust; j++) {
                cust_edges.insert({ i, j });
            }
        }
        veh_edges_subset = veh_edges;
        cust_edges_subset = cust_edges;

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
                    if (j > 0) veh_loc2[n+1].insert({ i, j, k });
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

    void preprocess_JLS() {
        //DOES NOT CURRENTLY WORK
        // Preprocess JLS Vars - Needs to remove xi_i,k's also - Can also be further improved by looking at individual vehicles - DOES NOT CURRENTLY WORK
        std::cout << "Preprocesing..." << std::endl;
        std::set<std::pair<int, int>> veh_edges_removed;
        std::set<std::pair<int, int>> cust_edges_removed;

        // Remove vehicle edges with too large demand or time constraint
        for (const auto& edge : veh_edges) {
            int k = edge.first;
            int j = edge.second;
            double t_max = std::max(veh_arr_time[k], cust_arr_time[j]);
            if (veh_dist[k][j] + cust_dist[j][0] > t_max || demand[j] > capacity - veh_occ[k]) {
                veh_edges_removed.insert(edge);
            }
        }
        
        // Remove customer edges
        int max_veh_time = 0;
        for (const int& k : vehicles) {
            max_veh_time = std::max(max_veh_time, veh_arr_time[k]);
        }
        int max_cap = capacity - *std::ranges::min_element(veh_occ);;
        for (const auto& edge : cust_edges) {
            int i = edge.first;
            int j = edge.second;
            int min_veh_dist = time_ub;
            for (const int& k : vehicles) {
                min_veh_dist = std::min(min_veh_dist, veh_dist[k][i]);
            }
            int dist = min_veh_dist + cust_dist[i][j] + cust_dist[j][0];
            double max_time = std::max({ cust_arr_time[i], cust_arr_time[j], max_veh_time });
            if (dist > max_time ||demand[i] + demand[j] > max_cap) {
                cust_edges_removed.insert(edge);
            }
        }

        std::set<std::pair<int, int>> temp_set;
        std::set_difference(veh_edges.begin(), veh_edges.end(), veh_edges_removed.begin(), veh_edges_removed.end(), std::inserter(temp_set, temp_set.begin()));
        veh_edges = temp_set;
        temp_set.clear();
        std::set_difference(cust_edges.begin(), cust_edges.end(), cust_edges_removed.begin(), cust_edges_removed.end(), std::inserter(temp_set, temp_set.begin()));
        cust_edges = temp_set;
        temp_set.clear();
        std::set_difference(veh_edges_subset.begin(), veh_edges_subset.end(), veh_edges_removed.begin(), veh_edges_removed.end(), std::inserter(temp_set, temp_set.begin()));
        veh_edges_subset = temp_set;
        temp_set.clear();
        std::set_difference(cust_edges_subset.begin(), cust_edges_subset.end(), cust_edges_removed.begin(), cust_edges_removed.end(), std::inserter(temp_set, temp_set.begin()));
        cust_edges_subset = temp_set;
        temp_set.clear();
        std::cout << "Edges Removed: " << veh_edges_removed.size() + cust_edges_removed.size() << std::endl;
    }

    void preprocess() {
        auto start = std::chrono::high_resolution_clock::now();

        // Preprocess YPS Vars
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

        std::cout << var_indices.size() << std::endl;
        std::cout << "Preprocesing..." << std::endl;
        std::set<int> indices_removed;

        // Remove vehicle edges with too large demand or time constraint
        for (int k : vehicles) {
            int remaning_cap = capacity - veh_occ[k];
            int k_times_cust_d = k * no_cust_d;
            for (int j : customers) {
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
        for (int i : customers) {
            for (int j : customers) {
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
                        }
                    }
                    else {
                        for (int k : vehicles) {
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
        std::set<int> temp_set;
        std::set_difference(var_indices.begin(), var_indices.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices = temp_set;
        temp_set.clear();
        std::set_difference(var_indices_subset.begin(), var_indices_subset.end(), indices_removed.begin(), indices_removed.end(), std::inserter(temp_set, temp_set.begin()));
        var_indices_subset = temp_set;
        std::cout << "Variables Removed: " << indices_removed.size() << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Creates the JLS formulation
    void create_jls_model(GRBEnv ENV, bool relax = true, bool silent = true) {
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);

        // Add Variables
        double obj;
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        std::string name;

        // y_cust variables
        y_cust.resize(no_cust + 1);
        y_cust_col.resize(no_cust + 1);
        for (int i : dest_and_cust) {
            y_cust[i].resize(no_cust + 1);
            y_cust_col[i].resize(no_cust + 1);
        }

        for (auto& edge : cust_edges_subset) {
            int i = edge.first;
            int j = edge.second;
            obj = cust_dist[i][j];
            name = "y_cust[" + itos(i) + "][" + itos(j) + "]";
            y_cust[i][j] = model->addVar(0.0, 1.0, obj, vtype, name);
        }


        // y_veh variables
        y_veh.resize(no_veh);
        y_veh_col.resize(no_veh);
        for (int k : vehicles) {
            y_veh[k].resize(no_cust + 1);
            y_veh_col[k].resize(no_cust + 1);
        }

        for (auto& edge : veh_edges_subset) {
            int k = edge.first;
            int i = edge.second;
            obj = veh_dist[k][i];
            name = "y_veh[" + itos(k) + "][" + itos(i) + "]";
            y_veh[k][i] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, name);
        }


        // x variables
        x.resize(no_cust + 1);
        for (int i : dest_and_cust) {
            x[i].resize(no_cust + 1);
            for (int j : customers) {
                if (i < j) {
                    x[i][j].resize(no_veh);
                }
            }
        }
        obj = 0;
        for (auto& edge : cust_edges_subset) {
            for (int k : vehicles) {
                int i = edge.first;
                int j = edge.second;
                name = "x[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                x[i][j][k] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, name);
            }
        }

        // xi variables
        xi.resize(no_cust + 1);
        for (int i : dest_and_cust) {
            xi[i].resize(no_veh);
        }

        for (auto& edge : veh_edges) {
            int k = edge.first;
            int i = edge.second;
            obj = -cust_dist[0][i]*demand[i];
            name = "xi[" + itos(i) + "][" + itos(k) + "]";
            xi[i][k] = model->addVar(0.0, 1.0, obj, vtype, name);
         }

        // Add Constraints
        GRBLinExpr lexpr, lexpr2;
        model->update();

        // 2 Edges per customers
        std::set<std::tuple<int, int, double>> veh_edges_not_in_constr;
        for (int k : vehicles) {
            for (int i : customers) {
                for (auto& edge : delta_edges(cust_edges_subset, {i}, dest_and_cust)) {
                    lexpr += x[edge.first][edge.second][k];
                }
                if (veh_edges_subset.find({ k, i }) != veh_edges_subset.end()) lexpr += y_veh[k][i];
                else veh_edges_not_in_constr.insert({ k, i, 1.0 });
                GRBConstr constr = model->addConstr(lexpr == 2 * xi[i][k]);
                for (auto& edge : veh_edges_not_in_constr) {
                    y_veh_col[k][i].addTerm(std::get<2>(edge), constr);
                }
                lexpr.clear();
                veh_edges_not_in_constr.clear();
            }
        }


        // Edge from vehicle if chosen
        for (int k : vehicles) {
            for (int i : dest_and_cust) {
                if (veh_edges_subset.find({ k, i }) != veh_edges_subset.end()) lexpr += y_veh[k][i];
                else veh_edges_not_in_constr.insert({ k, i, 1.0 });
            }
            GRBConstr constr = model->addConstr(lexpr == xi[0][k]);
            for (auto& edge : veh_edges_not_in_constr) {
                y_veh_col[k][std::get<1>(edge)].addTerm(std::get<2>(edge), constr);
            }
            lexpr.clear();
            veh_edges_not_in_constr.clear();
        }


        // Customer serviced at most once
        for (int i : customers) {
            for (int k : vehicles) {
                lexpr += xi[i][k];
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // Capacity Constraints
        for (int k : vehicles) {
            for (int i : customers) {
                lexpr += demand[i] * xi[i][k];
            }
            model->addConstr(lexpr <= (capacity - veh_occ[k]) * (xi[0][k] - y_veh[k][0]));
            lexpr.clear();
        }

        // Vehicle must reach destination if it already has passengers
        for (int k : vehicles) {
            model->addConstr(veh_occ[k] <= veh_occ[k] * xi[0][k]);
        }

        // Time Constraints
        for (int k : vehicles) {
            for (int l : dest_and_cust) {
                if (veh_edges_subset.find({ k, l }) != veh_edges_subset.end()) lexpr += (veh_dist[k][l] - veh_arr_time[k]) * y_veh[k][l];
                else veh_edges_not_in_constr.insert({ k, l, veh_dist[k][l] });
            }
            for (auto& edge : cust_edges_subset) {
                lexpr += cust_dist[edge.first][edge.second] * x[edge.first][edge.second][k];
            }
            for (int i : customers) {
                GRBConstr constr;
                if (veh_arr_time[k] > cust_arr_time[i]) {
                    constr = model->addConstr(lexpr + (veh_arr_time[k] - cust_arr_time[i]) * xi[i][k] <= 0);
                }
                else {
                    constr = model->addConstr(lexpr <= 0);
                }
                for (auto& edge : veh_edges_not_in_constr) {
                    y_veh_col[k][std::get<1>(edge)].addTerm(std::get<2>(edge), constr);
                }
            }
            lexpr.clear();
            veh_edges_not_in_constr.clear();
        }

        // Subtour Elimination
        for (int k : vehicles) {
            for (auto& edge : veh_edges_subset) {
                if (edge.first == k) lexpr += y_veh[k][edge.second];
            }
            for (auto& edge : cust_edges_subset) {
                lexpr += x[edge.first][edge.second][k];
            }
            for (int i : dest_and_cust) {
                lexpr2 += xi[i][k];
            }
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }

        // RCC for |S| == 3
        for (int i : customers) {
            for (int j = i + 1; j <= no_cust; j++) {
                for (int l = j + 1; l <= no_cust; l++) {
                    model->addConstr(y_cust[i][j] + y_cust[j][l] + y_cust[i][l] <= 3 - ceil((demand[i] + demand[j] + demand[l]) / float(capacity)));
                }
            }
        }


        // Calculate total demand
        int dem = 0;
        for (int i : customers) {
            dem += demand[i];
        }
        int vehicles_for_customers = int(ceil(dem / float(capacity)));
        for (auto& edge : edges_in(cust_edges_subset, customers)) {
            for (int k : vehicles) {
                lexpr += (capacity - veh_occ[k]) * x[edge.first][edge.second][k];
            }
        }
        for (int i : customers) {
            for (int k : vehicles) {
                lexpr2 += (capacity - veh_occ[k] - demand[i]) * xi[i][k];
            }
        }
        model->addConstr(lexpr <= lexpr2);
        lexpr.clear();
        lexpr2.clear();


        // x constraints
        for (auto& edge : cust_edges_subset) {
            int i = edge.first;
            int j = edge.second;
            for (int k : vehicles) {
                lexpr += x[i][j][k];
            }
            model->addConstr(y_cust[i][j] == lexpr);
            lexpr.clear();
        }


        for (auto& edge : cust_edges_subset) {
            int i = edge.first;
            int j = edge.second;
            for (int k : vehicles) {
                model->addConstr(x[i][j][k] >= xi[i][k] + xi[j][k] + y_cust[i][j] - 2);
                model->addConstr(x[i][j][k] <= xi[i][k]);
                model->addConstr(x[i][j][k] <= xi[j][k]);
                model->addConstr(x[i][j][k] <= y_cust[i][j]);
            }
        }


        for (auto& edge : veh_edges_subset) {
            int k = edge.first;
            int i = edge.second;
            model->addConstr(y_veh[k][i] <= xi[i][k]);
            model->addConstr(y_veh[k][i] <= xi[0][k]);
        }
    }

    // Creates the GP formulation
    void create_gp_model(GRBEnv ENV, bool relax = true, bool silent = true) {
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);

        // Add Variables
        double obj = 0;
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        std::string name;

        // x variables
        x.resize(no_cust + 1);
        for (int i : dest_and_cust) {
            x[i].resize(no_cust + 1);
            for (int j : dest_and_cust) {
                x[i][j].resize(no_veh);
            }
        }

        for (int i : dest_and_cust) {
            double ub = (i == 0) ? 0.0 : 1.0;
            for (int j : dest_and_cust) {
                if (i != j) {
                    for (int k : vehicles) {
                        name = "x[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                        x[i][j][k] = model->addVar(0.0, ub, obj, vtype, name);
                    }
                }
            }
        }

        // y_veh variables
        y_veh.resize(no_veh);
        for (int k : vehicles) {
            y_veh[k].resize(no_cust + 1);
        }

        for (int k : vehicles) {
            for (int i : dest_and_cust) {
                name = "y_veh[" + itos(k) + "][" + itos(i) + "]";
                y_veh[k][i] = model->addVar(0.0, 1.0, obj, vtype, name);
            }
        }

        // q variables
        q.resize(no_cust + 1);
        double max_demand = fmin(no_cust, capacity - *std::ranges::min_element(veh_occ));

        for (int i : dest_and_cust) {
            name = "q[" + itos(i) + "]";
            q[i] = model->addVar(0.0, max_demand, obj, GRB_CONTINUOUS, name);
        }


        model->update();
        // Add Constraints
        GRBLinExpr lexpr, lexpr2, lexpr3;

        // 4c
        for (int k : vehicles) {
            for (int j : dest_and_cust) {
                lexpr += y_veh[k][j];
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // 4d
        for (int k : vehicles) {
            for (int j : dest_and_cust) {
                lexpr += y_veh[k][j];
            }
            lexpr2 += y_veh[k][0];
            for (int i : customers) {
                lexpr2 += x[i][0][k];
            }
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }

        // 4e
        for (int k : vehicles) {
            for (int j : customers) {
                lexpr += y_veh[k][j];
                for (int i : customers) {
                    if (i != j) lexpr += x[i][j][k];
                }
                for (int i : dest_and_cust) {
                    if (i != j) lexpr2 += x[j][i][k];
                }
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }

        // 4f1
        for (int i : customers) {
            for (int j : customers) {
                if (i != j) {
                    for (int k : vehicles) {
                        lexpr += x[i][j][k];
                    }
                    model->addConstr(q[j] >= q[i] - (no_cust + 1) * (1 - lexpr) + demand[j]);
                    lexpr.clear();
                }
            }
        }

        // 4f2
        for (int k : vehicles) {
            for (int j : customers) {
                model->addConstr(q[j] >= veh_occ[k] - (no_cust + 1) * (1 - y_veh[k][j]) + demand[j]);
            }
        }


        // 4g
        for (int i : customers) {
            for (int j : dest_and_cust) {
                if (i != j) {
                    for (int k : vehicles) {
                        lexpr += x[i][j][k];
                    }
                }
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // 4h
        for (int k : vehicles) {
            for (int i : customers) {
                for (int j : dest_and_cust) {
                    if (i != j) lexpr += x[i][j][k] * demand[i];
                    }
            }
            model->addConstr(lexpr <= (capacity - veh_occ[k]));
            lexpr.clear();
            lexpr2.clear();
        }

        // 4i
        for (int k : vehicles) {
            for (int j : dest_and_cust) {
                lexpr += veh_occ[k] * y_veh[k][j];
            }
            model->addConstr(veh_occ[k] <= lexpr);
            lexpr.clear();
        }

        // 4j
        for (int k : vehicles) {
            for (int i : customers) {
                for (int j : dest_and_cust) {
                    lexpr += veh_dist[k][j] * y_veh[k][j];
                }
                for (int l : customers) {
                    for (int j : dest_and_cust){
                        if (l != j) {
                            int time = (l < j) ? cust_dist[l][j] : cust_dist[j][l];
                            lexpr += time * x[l][j][k];
                        }
                    }
                }
                for (int j : dest_and_cust) {
                    if (i != j) lexpr2 += x[i][j][k];
                }
                model->addConstr(lexpr <= fmin(veh_arr_time[k], cust_arr_time[i]) + time_ub * (1 - lexpr2));
                lexpr.clear();
                lexpr2.clear();
            }
        }



        // Set Objective

        for (int k : vehicles) {
            for (int i : dest_and_cust) {
                lexpr += veh_dist[k][i] * y_veh[k][i];
                if (i > 0) {
                    for (int j : dest_and_cust) {
                        if (i != j) {
                            int time = (i < j) ? cust_dist[i][j] : cust_dist[j][i];
                            lexpr += time * x[i][j][k];
                        }
                    }
                }
            }
        }

        for (int i : customers) {
            for (int j : dest_and_cust) {
                if (i != j) {
                    for (int k : vehicles) {
                        lexpr2 += x[i][j][k];
                    }
                }
            }

            lexpr -= demand[i] * cust_dist[0][i] * lexpr2;
            lexpr2.clear();
        }
        model->setObjective(lexpr, GRB_MINIMIZE);
    }
   
    // Creates the GP formulation
    void create_model(GRBEnv ENV, bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);

        // Add Variables
        no_veh_vars = no_veh * (no_cust + 1);
        no_total_vars = no_veh * (no_cust + 1) + no_veh * ((no_cust + 1)* (no_cust + 1) - no_cust - 1);
        if (var_indices.size() == 0){
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

        // phi variables
        phi.resize(no_cust + 1);

        for (int i : customers) {
            name = "phi[" + itos(i) + "]";
            phi[i] = model->addVar(demand[i], std::max(double(demand[i]), max_cap), obj, GRB_CONTINUOUS, name);
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
                int k = floor((idx - 1 - no_veh * no_cust_d) / (pow(no_cust_d,2) - no_cust - 1));
                int i = floor(
                    (idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1)) / no_cust_d
                    + 1);
                int j = (idx - 1 - no_veh * no_cust_d -k * (pow(no_cust_d, 2) - no_cust - 1) - (i-1)*no_cust_d);
                if (i != j) {
                    obj = cust_dist[i][j] - demand[i]*cust_dist[i][0];
                    name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
                    vars[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
                }
            }
        }

        // Add Constraints
        GRBLinExpr lexpr, lexpr2, lexpr3;

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
            for (int i: customers) {
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

        // 1.9
        for (int k : vehicles) {
            for (int j : customers) {
                int idx = k * no_cust_d + j + 1;
                if (var_indices_subset.find(idx) != var_indices_subset.end()) model->addConstr(phi[j] >= (veh_occ[k] + demand[j]) * vars[idx]);
                else model->addConstr(phi[j] >= 0);
            }
        }

        // 1.10
        for (int i : customers) {
            for (int j : customers) {
                if (i != j) {
                    for (int k : vehicles) {
                        int idx = no_veh * no_cust_d + k * (pow(no_cust_d, 2) - no_cust - 1) + (i - 1) * no_cust_d + j + 1;
                        if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
                    }
                    //model->addConstr(phi[i] - (max_cap + demand[j]) * (1 - lexpr) + demand[j] <= phi[j]);
                    model->addConstr(phi[i] + demand[j] - (1 - lexpr) * std::max(double(demand[i]), max_cap) <= phi[j]);
                    lexpr.clear();
                }
            }
        }

        model->update();
        no_base_constraints = model->get(GRB_IntAttr_NumConstrs);
        std::cout << "Model Constructed." << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        model_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void optimise_model() {
        auto start = std::chrono::high_resolution_clock::now(); 
        model->optimize();
        auto stop = std::chrono::high_resolution_clock::now();
        solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void add_variables(bool& new_vars) {
        auto start = std::chrono::high_resolution_clock::now();

        add_variables_cg(new_vars, max_cap, veh_dist, cust_dist, veh_occ, veh_arr_time, cust_arr_time,
            demand, customers, no_veh, no_cust, no_cust_d, capacity, no_veh_vars, no_total_vars, no_base_constraints, model, vars, var_indices, var_indices_subset, rounded_cap, false);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
        add_vars_time += duration.count();
    }
    void add_all_variables(bool& new_vars) {
        add_variables_cg(new_vars, max_cap, veh_dist, cust_dist, veh_occ, veh_arr_time, cust_arr_time,
            demand, customers, no_veh, no_cust, no_cust_d, capacity, no_veh_vars, no_total_vars, no_base_constraints, model, vars, var_indices, var_indices_subset, rounded_cap, true);
    }

    void remove_constraints() { remove_constraints_from_model(model, no_base_constraints);  };

    void make_integral() {
        optimise_model();
        model->getEnv().set(GRB_IntParam_OutputFlag, 1);
        for (int idx : var_indices_subset) {
            if (idx >= no_veh_vars) vars[idx].set(GRB_CharAttr_VType, GRB_BINARY);
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


            //convert_fmrsp_to_cvrp_solution(
            //    cvrp_to_fmrsp,
            //    cvrp_no_cust,
            //    cvrp_demand,
            //    cvrp_customers,
            //    cvrp_solution,
            //    demand,
            //    customers,
            //    vehicles,
            //    y_cust,
            //    y_veh,
            //    xi
            //);

            //convert_fmrsp_to_cvrp_solution(
            //    cvrp_to_fmrsp,
            //    cvrp_no_cust,
            //    cvrp_demand,
            //    cvrp_customers,
            //    cvrp_solution,
            //    demand,
            //    customers,
            //    dest_and_cust,
            //    vehicles,
            //    x,
            //    y_veh
            //);

            convert_YPS_to_cvrp_solution(
                cvrp_to_fmrsp,
                cvrp_no_cust,
                cvrp_customers,
                cvrp_demand,
                cvrp_solution,
                demand,
                customers,
                dest_and_cust,
                vehicles,
                vars,
                no_veh_vars,
                no_total_vars,
                no_veh,
                no_cust,
                no_cust_d,
                var_indices_subset
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
                //separate_rcc_exactly(new_cuts, max_violation, cvrp_solution, cvrp_customers, capacity, cvrp_demand, 1, 0.001, 0.001, RC_separation_time, env);
            }
            if (max_violation < threshold) break;

            if (cvrp_no_cust < no_cust) {
                convert_cuts_to_fmrsp(cvrp_to_fmrsp, new_cuts);
            }

            //add_rc_inequalities(model, y_cust, customers, new_cuts);
            int no_new_cuts = 0;
            add_rc_inequalities(no_new_cuts, rounded_cap, model, vars, customers, vehicles, new_cuts, no_veh, no_cust, no_cust_d, var_indices_subset, total_rcc);
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