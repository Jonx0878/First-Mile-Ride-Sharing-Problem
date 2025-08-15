#pragma once
#include <any>
#include <gurobi_c++.h>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "load.h"
#include "rc.h"
#include "utils.h"


class RB {
private:

    const std::string vert_to_str(const std::vector<int>& vertex) {
        std::stringstream str;
        for (const int& j : vertex) {
            str << std::setfill('0') << std::setw(cust_magn) << j;
        }
        return str.str();
    }
    const std::vector<int> str_to_vert(const std::string& str) {
        const int size = str.size() / cust_magn;
        std::vector<int> vertex(size);  // Pre-allocate space to avoid resizing
        for (int i = 0; i < size; ++i) {
            int value = 0;
            // Manually parse the integer without calling std::stoi or std::from_chars
            for (int j = 0; j < cust_magn; ++j) {
                value = value * 10 + (str[i * cust_magn + j] - '0');
            }
            vertex[i] = value;  // Directly set the value in the pre-allocated vector
        }
        return vertex;
    }

    void routes_pr_vehicle(const int& k) {
        const int min_prof_dist = (veh_occ[k] == 0) ? 0 : veh_dist[k][0];
        std::set<int> routes_to_remove;
        std::vector<int> best_child_cost;
        std::vector<int> stop_index;
        stop_index.emplace_back(routes[k].size());
        std::vector<int> depot = { 0 };
        routes[k].emplace_back(depot);
        route_demand[k].emplace_back(0);
        route_dist[k].emplace_back(veh_dist[k][0]);
        route_stops[k].emplace_back(0);
        route_min_time[k].emplace_back(veh_arr_time[k]);
        route_cost[k].emplace_back(veh_dist[k][0]);
        best_child_cost.emplace_back(veh_dist[k][0]);
        const int size = capacity[k] - veh_occ[k];
        // Determine possible routes of 1 stop
        stop_index.emplace_back(routes[k].size());
        if (size == 0) return;
        for (const int& i : customers) {
            std::vector<int> route;
            int dist = veh_dist[k][i] + cust_dist[i][0];
            int min_time = std::min(veh_arr_time[k], cust_arr_time[i]);
            if (demand[i] > size || dist > min_time) continue;
            route.emplace_back(i);
            routes[k].emplace_back(route);
            route_demand[k].emplace_back(demand[i]);
            route_dist[k].emplace_back(dist - cust_dist[i][0]);
            route_stops[k].emplace_back(1);
            route_min_time[k].emplace_back(min_time);
            route_cost[k].emplace_back(dist - demand[i] * cust_dist[i][0]);
            best_child_cost.emplace_back(dist - demand[i] * cust_dist[i][0]);
            if (dist - demand[i] * cust_dist[i][0] > min_prof_dist) routes_to_remove.insert(routes[k].size() - 1);

        }

        stop_index.emplace_back(routes[k].size());
        // Determine rest of possible routes
        for (int l = 2; l <= size; l++) {
            std::unordered_map<std::string, int> ordered_best_costs;
            std::unordered_map<std::string, std::unordered_set<int>> index_to_routes;
            for (int m = stop_index[l - 1]; m < stop_index[l]; m++) {
                for (const int& i : customers) {
                    if (std::find(routes[k][m].begin(), routes[k][m].end(), i) != routes[k][m].end()) continue;
                    int dem = route_demand[k][m] + demand[i];
                    int dist = route_dist[k][m] + cust_dist[routes[k][m][routes[k][m].size() - 1]][i] + cust_dist[i][0];
                    int min_time = std::min(route_min_time[k][m], cust_arr_time[i]);
                    if (dem > size || dist > min_time) continue;
                    std::vector<int> route = routes[k][m];
                    route.emplace_back(i);
                    routes[k].emplace_back(route);
                    route_demand[k].emplace_back(dem);
                    route_dist[k].emplace_back(dist - cust_dist[i][0]);
                    route_stops[k].emplace_back(l);
                    route_min_time[k].emplace_back(min_time);
                    int cost = dist - demand[i] * cust_dist[i][0];
                    for (const int& j : routes[k][m]) cost -= demand[j] * cust_dist[j][0];
                    route_cost[k].emplace_back(cost);
                    std::vector<int> cust_sorted = route;
                    std::sort(cust_sorted.begin(), cust_sorted.end());
                    std::string index = vert_to_str(cust_sorted);

                    bool exists = false;
                    if (cost > best_child_cost[m]) {
                        exists = true;
                        routes_to_remove.insert(routes[k].size() - 1);
                        if (ordered_best_costs.find(index) == ordered_best_costs.end()) ordered_best_costs[index] = best_child_cost[m];
                        else ordered_best_costs[index] = std::min(ordered_best_costs[index], best_child_cost[m]);
                    }
                    else if (cost > min_prof_dist) routes_to_remove.insert(routes[k].size() - 1);

                    // If routes with same costumers in different order
                    if (!exists && ordered_best_costs.find(index) == ordered_best_costs.end()) ordered_best_costs[index] = cost;
                    // Remove if better route with same set of customers
                    else if (ordered_best_costs[index] < cost) routes_to_remove.insert(routes[k].size() - 1);
                    // Remove other previously created worse routes
                    else if (ordered_best_costs[index] > cost) {
                        for (const int& r : index_to_routes[index]) {
                            routes_to_remove.insert(r);
                        }
                        // Update Cost
                        ordered_best_costs[index] = cost;
                    }
                    index_to_routes[index].insert(routes[k].size() - 1);
                    // Update best_child
                    best_child_cost.emplace_back(std::min(cost, best_child_cost[m]));
                }
            }
            stop_index.emplace_back(routes[k].size());
        }

        {
            std::vector<std::vector<int>> temp_routes;
            std::vector<int> temp_route_stops;
            std::vector<int> temp_route_costs;
            for (int idx = 0; idx < routes[k].size(); idx++) {
                if (routes_to_remove.find(idx) != routes_to_remove.end()) continue;
                temp_routes.emplace_back(routes[k][idx]);
                temp_route_stops.emplace_back(route_stops[k][idx]);
                temp_route_costs.emplace_back(route_cost[k][idx]);
            }
            routes[k] = temp_routes;
            route_stops[k] = temp_route_stops;
            route_cost[k] = temp_route_costs;
        }
    }

    void worker_remaining(std::atomic<int>& task_counter) {
        while (true) {
            int task = task_counter.fetch_add(1);   // Get the next task and increment the counter
            if (task >= no_veh) break;              // Stop if all tasks are processed
            routes_pr_vehicle(task);                // Process the task
        }
    }

    void retrieve_solution(std::set<std::tuple<int, int, double>>& solution) {
        // Create initial customer edge array
        std::vector<std::vector<double>> sol_values(dest_and_cust.size());
        for (const int& i : dest_and_cust) {
            sol_values[i].resize(dest_and_cust.size());
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                sol_values[i][j] = 0;
            }
        }
        // Add solution values
        for (const int& idx : var_indices) {
            const double val = vars[idx].get(GRB_DoubleAttr_X);
            if (val < 1e-3) continue;
            const int k = idx_to_route[idx].first;
            const int route_num = idx_to_route[idx].second;
            const int route_length = routes[k][route_num].size();
            const int first_cust = routes[k][route_num][0];
            sol_values[first_cust][0] += val;
            if (first_cust == 0) continue;
            sol_values[routes[k][route_num][route_length - 1]][0] += val;
            for (int i = 0; i < route_length - 1; i++) {
                sol_values[routes[k][route_num][i]][routes[k][route_num][i + 1]] += val;
            }
        }

        // Add positive solution values to solution
        for (const int& i : customers) {
            for (const int& j : dest_and_cust) {
                if (i == j) continue;
                if (sol_values[i][j] >= 1e-3) {
                    if (j == 0) solution.insert({ j, i, sol_values[i][j] });
                    else solution.insert({ i, j, sol_values[i][j] });
                }
            }
        }
    }

    void add_rc_inequalities(std::vector<RC_Inequality>& violated_inequalties, bool& cuts_added) {
        cuts_added = false;
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        for (RC_Inequality& cut : violated_inequalties) {
            GRBLinExpr lexpr;
            for (const int& i : cut.set) {
                for (const int& idx : route_cust_idx[i]) {
                    const int k = idx_to_route[idx].first;
                    const int route_num = idx_to_route[idx].second;
                    if (routes[k][route_num].size() == 1) continue;
                    auto it = std::find(routes[k][route_num].begin(), routes[k][route_num].end() - 1, i);
                    if (it == routes[k][route_num].end() - 1) continue;
                    if (cut.set.find(*(it + 1)) != cut.set.end()) lexpr += vars[idx];
                }
            }

            if (lexpr.getValue() >= cut.rhs + 1e-3) {
                cuts_added = true;
                model->addConstr(lexpr <= cut.rhs);
            }
        }
    }

public:
    // Instance Properties
    std::string name;
    int no_veh;
    int no_cust;
    int no_cust_d;
    int time_ub;
    int no_veh_vars = 0;
    int no_total_vars = 0;
    int max_cap;
    GRBEnv env = GRBEnv();
    int cust_magn = 1;
    int num_threads = 12;
    int time_limit = 3600;

    // Solution Values
    double initial_lp;
    double strengthened_lp;
    double mip;

    // Vehicle Data
    std::vector<int> capacity;
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
    std::unordered_set<int> var_indices;
    std::unordered_set<int> var_indices_subset;

    // Graph Data
    std::vector<std::vector<std::vector<int>>> routes;
    std::vector<std::pair<int, int>> idx_to_route;
    std::vector<std::vector<int>> route_cost;
    std::vector<std::vector<int>> route_demand;
    std::vector<std::vector<int>> route_dist;
    std::vector<std::vector<int>> route_stops;
    std::vector<std::vector<int>> route_min_time;
    std::vector<int> first_route_index;
    std::vector<std::set<int>> route_cust_idx;
    
    // Column Generation
    int outer_iterations = 0;

    // Timings
    double total_lp_solve_time = 0;
    double lp_solve_time = 0;
    double lp_rcc_solve_time = 0;
    double mip_solve_time = 0;
    double model_construction_time = 0;
    double graph_construction_time = 0;
    double RC_separation_time = 0;

    // Determine vertices
    void determine_routes() {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Constructing Routes..." << std::endl;
        std::set<int> routes_to_remove;
        std::vector<int> best_child_cost;
        routes.resize(no_veh);
        route_cost.resize(no_veh);
        route_demand.resize(no_veh);
        route_dist.resize(no_veh);
        route_stops.resize(no_veh);
        route_min_time.resize(no_veh);

        std::atomic<int> task_counter(0); // Atomic counter to assign tasks
        std::vector<std::thread> threads;
        for (int i = 0; i < std::min(num_threads, no_veh); i++) {
            threads.emplace_back([this, &task_counter]() {
                this->worker_remaining(task_counter);  // Call member function on the instance
                });
        }

        // Join threads
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }

        int idx = 0;
        for (const int& k : vehicles) {
            first_route_index.emplace_back(idx);
            for (int i = 0; i < routes[k].size(); i++) {
                idx_to_route.emplace_back(std::make_pair(k, i));
                var_indices.insert(idx);
                idx++;
                //std::cout << idx - 1  << ": (";
                //for (int j : routes[k][i]) {
                //    std::cout << " " << j << ",";
                //}
                //std::cout << ") " << route_cost[k][i] << std::endl;
            }
        }
        first_route_index.emplace_back(idx);
        //for (int j = 0; j < routes.size();j++) {
        //    std::cout << j << ": (";
        //    for (int i : routes[j]) {
        //        std::cout << " " << i << ",";
        //    }
        //    std::cout << ")" << std::endl;
        //} 
        std::cout << "# Routes: " << idx << std::endl;

        route_demand.clear();
        route_dist.clear();
        route_stops.clear();
        route_min_time.clear();
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Initializes instance
    RB(const std::string& filename, const int& threads, const int& time) {
        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        no_veh = inst_data.no_veh;
        no_cust = inst_data.no_cust;
        no_cust_d = no_cust + 1;
        time_ub = inst_data.time_ub;
        num_threads = threads;
        time_limit = time;

        capacity = inst_data.capacity;
        veh_coords = inst_data.veh_coords;
        veh_occ = inst_data.veh_occ;
        veh_arr_time = inst_data.veh_arr_time;
        veh_dist = inst_data.veh_dist;

        cust_coords = inst_data.cust_coords;
        demand = inst_data.demand;
        cust_arr_time = inst_data.cust_arr_time;
        cust_dist = inst_data.cust_dist;
        cust_magn = ceil(log10(no_cust_d));

        // Populate vertices, customers, and deste_and_cust
        for (int k = 0; k < no_veh; k++) {
            vehicles.insert(k);
        }
        for (int i = 0; i <= no_cust; i++) {
            dest_and_cust.insert(i);
            if (i > 0) customers.insert(i);
        }

        rounded_cap.resize(no_cust * (no_cust - 1) / 2.0);
        max_cap = 0;
        for (const int& k : vehicles) {
            max_cap = std::max(capacity[k] - veh_occ[k], max_cap);
        }

        determine_routes();
    }

    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Get subset indices
        std::cout << "Building Model..." << std::endl;
        if (var_indices_subset.size() == 0) var_indices_subset = var_indices;
        // Determine routes of customer i
        route_cust_idx.resize(no_cust_d);
        for (const int& k : vehicles) {
            for (int idx = 0; idx < routes[k].size(); idx++) {
                const int route_idx = first_route_index[k] + idx;
                for (const int& i : routes[k][idx]) {
                    route_cust_idx[i].insert(route_idx);
                }
            }
        }
        // Create Model
        model = new GRBModel(env);
        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);
        model->getEnv().set(GRB_IntParam_Threads, num_threads);

        // Create Variables
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        vars.resize(var_indices.size());

        {
            const int size = int(var_indices.size());
            double* ub = new double[size];
            char* vtypes = new char[size];
            double* objs = new double[size];
            std::string* names = new std::string[size];

            for (int idx = 0; idx < var_indices.size(); idx++) {
                ub[idx] = 1.0;
                objs[idx] = route_cost[idx_to_route[idx].first][idx_to_route[idx].second];
                vtypes[idx] = vtype;
                names [idx] ="lamda[" + itos(idx) + "]";
            }
            GRBVar* var_list = model->addVars(NULL, ub, objs, vtypes, names, size);

            vars.assign(var_list, var_list + size);
        }
        // Add constraints
        GRBLinExpr lexpr, lexpr2, lexpr3;
        // At most service customer once
        for (const int& i : customers) {
            for (const int& idx : route_cust_idx[i]) {
                lexpr += vars[idx];
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }
        // One vehicle can only take one path
        for (const int& k : vehicles) {
            for (int idx = first_route_index[k]; idx < first_route_index[k + 1]; idx++) {
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
            }
            if (veh_occ[k] == 0) model->addConstr(lexpr <= 1);
            else model->addConstr(lexpr == 1);
            //model->addConstr(veh_occ[k] <= veh_occ[k] * lexpr);
            lexpr.clear();
        }
        std::cout << "Model Constructed." << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        model_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void optimise_model() {
        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        auto stop = std::chrono::high_resolution_clock::now();
        total_lp_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }


    void remove_variables() { 
        int k = 0;
        for (int idx : var_indices_subset) {
            if (idx >= first_route_index[k + 1]) k++;
            if (veh_occ[k] > 0 && std::find(first_route_index.begin(), first_route_index.end(), idx) != first_route_index.end()) // Do not remove route directly to dest if vehicle must move
            if (vars[idx].get(GRB_DoubleAttr_X) <= 0.0001 && vars[idx].get(GRB_DoubleAttr_RC) > 0.01) {
                model->remove(vars[idx]);
                var_indices_subset.erase(idx);
            }
        }
    }

    void solve_root_relaxation(bool separate_rci = true, bool silent = true) {
        if (model_construction_time == 0) create_model(true, silent);

        auto start = std::chrono::high_resolution_clock::now();
        optimise_model();
        initial_lp = model->get(GRB_DoubleAttr_ObjVal);
        auto stop = std::chrono::high_resolution_clock::now();
        lp_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();

        if (separate_rci) {
            std::cout << initial_lp << std::endl;
            std::set<std::tuple<int, int, double>> solution;
            std::vector<double> no_vehicles_used;
            std::vector<RC_Inequality> violated_inequalties;
            double max_viol = 0;
            bool cuts_added = false;
            while (true) {
                retrieve_solution(solution);
                separate_rcc_exactly(violated_inequalties, max_viol, solution, customers, max_cap, demand, 1, RC_separation_time, env);
                add_rc_inequalities(violated_inequalties, cuts_added);
                if (!cuts_added) break;

                solution.clear();
                no_vehicles_used.clear();
                violated_inequalties.clear();
                optimise_model();
                std::cout << model->get(GRB_DoubleAttr_ObjVal) << std::endl;

            }
        }
        strengthened_lp = model->get(GRB_DoubleAttr_ObjVal);
        stop = std::chrono::high_resolution_clock::now();
        lp_rcc_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void solve_to_optimality(bool separate_rci = true, bool silent = true) {
        if (total_lp_solve_time == 0) {
            solve_root_relaxation(separate_rci, silent);
            mip_solve_time += lp_rcc_solve_time;
        }
        model->getEnv().set(GRB_IntParam_OutputFlag, 1);
        model->getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
        model->getEnv().set(GRB_DoubleParam_MIPGap, 0);
        model->getEnv().set(GRB_DoubleParam_MIPGapAbs, 0.99);

        for (const int& idx : var_indices) {
            vars[idx].set(GRB_CharAttr_VType, GRB_BINARY);
        }
        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        mip = model->get(GRB_DoubleAttr_ObjVal);
        auto stop = std::chrono::high_resolution_clock::now();
        mip_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Prints all routes
    void print_all_routes() {
        for (const int k : vehicles) {
            std::cout << "Vehicle " << k << ":" << std::endl;
            for (int r = 0; r < routes[k].size(); r++) {
                for (const int& i : routes[k][r]) {
                    std::cout << i << "->";
                }
                std::cout << "d" << std::endl;
            }
        }
    }
};

