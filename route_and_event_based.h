#pragma once
#include <any>
#include <gurobi_c++.h>
#include <set>
#include <unordered_map>
#include <thread>

#include "load.h"
#include "rc.h"
#include "utils.h"



class REB {
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

    void remaining_vertices(const int& k, std::vector<std::vector<int>>& min_costs, std::vector<std::vector<int>>& vert_costs,
        std::vector<std::vector<int>>& min_dist_prev, std::vector<std::vector<int>>& no_stops_indices,
        std::vector<std::vector<std::vector<int>>>& min_costs_n_step) {
        auto start_time = std::chrono::high_resolution_clock::now();
        const int size = capacity[k] - veh_occ[k];
        no_stops_indices[k].emplace_back(0);
        std::vector<int> vertex = { 0 };
        veh_vertices[k].emplace_back(vertex);
        no_stops[k].emplace_back(0);
        demand_verts[k].emplace_back(0);
        min_times[k].emplace_back(time_ub);
        min_dist_prev[k].emplace_back(0);
        no_stops_indices[k].emplace_back(veh_vertices[k].size());
        if (size == 0) return;
        for (const int& i : customers) {
            const int min_time = std::min(cust_arr_time[i], veh_arr_time[k]);
            if (demand[i] > size || veh_dist[k][i] + cust_dist[i][0] > min_time) continue;
            vertex = { i };
            veh_vertices[k].emplace_back(vertex);
            no_stops[k].emplace_back(1);
            demand_verts[k].emplace_back(demand[i]);
            min_times[k].emplace_back(min_time);
            min_dist_prev[k].emplace_back(veh_dist[k][i]);
        }
        no_stops_indices[k].emplace_back(veh_vertices[k].size());
        for (int l = 2; l <= size; l++) {
            std::unordered_map<std::string, int> duplicate_detector;
            for (const int& i : customers) {
                for (int n = no_stops_indices[k][l - 1]; n < no_stops_indices[k][l]; n++) {
                    if (std::find(veh_vertices[k][n].begin(), veh_vertices[k][n].end(), i) != veh_vertices[k][n].end()) continue;
                    const int demand_sum = demand_verts[k][n] + demand[i];
                    const int min_time = std::min(min_times[k][n], cust_arr_time[i]);
                    const int min_dist_prev_vert = min_dist_prev[k][n] + cust_dist[veh_vertices[k][n][0]][i];
                    if (demand_sum > size || min_dist_prev_vert + cust_dist[i][0] > min_time) continue;
                    // Add vertex
                    vertex = { i };
                    int first_cust = veh_vertices[k][n][0];
                    bool inserted = false;
                    for (int m = 1; m < l - 1; m++) {
                        if (!inserted && first_cust < veh_vertices[k][n][m]) {
                            vertex.emplace_back(first_cust);
                            inserted = true;
                        }
                        vertex.emplace_back(veh_vertices[k][n][m]);
                    }
                    if (!inserted) vertex.emplace_back(first_cust);
                    const std::string str = vert_to_str(vertex);
                    if (duplicate_detector.find(str) != duplicate_detector.end()) {
                        int idx = duplicate_detector[str];
                        min_dist_prev[k][idx] = std::min(min_dist_prev[k][idx], min_dist_prev_vert);
                    }
                    else {
                        duplicate_detector[str] = veh_vertices[k].size();
                        veh_vertices[k].emplace_back(vertex);
                        no_stops[k].emplace_back(l);
                        demand_verts[k].emplace_back(demand_sum);
                        min_times[k].emplace_back(min_time);
                        min_dist_prev[k].emplace_back(min_dist_prev_vert);
                    }
                }
            }
            no_stops_indices[k].emplace_back(veh_vertices[k].size());
        }
        no_stops_indices[k].emplace_back(veh_vertices[k].size());

        // Determine minimum possible cost for each number of possible steps and max and min next distance
        min_costs[k].resize(veh_vertices[k].size());
        vert_costs[k].resize(veh_vertices[k].size());
        min_costs_n_step[k].resize(veh_vertices[k].size());
        std::vector<int> max_next_distance(veh_vertices[k].size());
        std::vector<int> min_next_distance(veh_vertices[k].size());

        for (int idx = veh_vertices[k].size() - 1; idx >= 0; idx--) {
            // Calculate the minimum possible cost of having this event as the last one
            const int idx_cust = veh_vertices[k][idx][0];
            int min_cost = min_dist_prev[k][idx] + cust_dist[idx_cust][0];
            for (const int& i : veh_vertices[k][idx]) {
                min_cost -= demand[i] * cust_dist[i][0];
            }
            vert_costs[k][idx] = min_cost;
            max_next_distance[idx] = -time_ub;
            min_next_distance[idx] = time_ub;

            // Set the min cost of 0 step equal vert_costs and rest equal to an upper bound
            min_costs_n_step[k][idx].resize(size - demand_verts[k][idx] + 1);
            min_costs_n_step[k][idx][0] = min_cost;
            if (demand_verts[k][idx] == size) continue;

            for (int n = 1; n <= size - demand_verts[k][idx]; n++) {
                min_costs_n_step[k][idx][n] = time_ub;
            }
            // Update future min costs by using the ones from events reachable in 1 step
            const int start = no_stops_indices[k][no_stops[k][idx] + 1];
            const int end = no_stops_indices[k][no_stops[k][idx] + 2];
            for (int i = start; i < end; i++) {
                int same_cust_counter = 0;
                for (const int& l : veh_vertices[k][idx]) {
                    if (std::find(veh_vertices[k][i].begin() + 1, veh_vertices[k][i].end(), l) != veh_vertices[k][i].end()) same_cust_counter++;
                    else break;
                }
                if (same_cust_counter < no_stops[k][idx]) continue;
                const int next_dist = cust_dist[idx_cust][veh_vertices[k][i][0]];
                const int min_cost_diff = min_dist_prev[k][idx] + next_dist - min_dist_prev[k][i];
                for (int n = 1; n <= size - demand_verts[k][i] + 1; n++) {
                    min_costs_n_step[k][idx][n] = std::min(min_costs_n_step[k][idx][n], min_costs_n_step[k][i][n - 1] + min_cost_diff);
                }
                max_next_distance[idx] = std::max(max_next_distance[idx], min_dist_prev[k][idx] + next_dist);
                min_next_distance[idx] = std::min(min_next_distance[idx], min_dist_prev[k][idx] + next_dist);
            }
        }

        // Calculate min_cost and events with visiting same customers
        std::unordered_map<int, std::string> vert_to_num;
        std::unordered_map<std::string, std::vector<int>> num_to_verts;
        for (int idx = 1; idx < veh_vertices[k].size(); idx++) {
            min_costs[k][idx] = *std::ranges::min_element(min_costs_n_step[k][idx]);
            std::vector<int> vertex;
            int first_cust = veh_vertices[k][idx][0];
            bool inserted = false;
            for (int m = 1; m < veh_vertices[k][idx].size(); m++) {
                if (!inserted && first_cust < veh_vertices[k][idx][m]) {
                    vertex.emplace_back(first_cust);
                    inserted = true;
                }
                vertex.emplace_back(veh_vertices[k][idx][m]);
            }
            if (!inserted) vertex.emplace_back(first_cust);
            const std::string str = vert_to_str(vertex);
            vert_to_num[idx] = str;
            num_to_verts[str].emplace_back(idx);
        }

        // Remove vertices that cannot be optimal choices due to a better order
        const int min_profitable_cost = (veh_occ[k] == 0) ? 0 : veh_dist[k][0];
        std::vector<bool> keep_vertex(veh_vertices[k].size());
        keep_vertex[0] = true;
        for (int idx = 1; idx < veh_vertices[k].size(); idx++) {
            keep_vertex[idx] = true;

            if (min_costs[k][idx] > min_profitable_cost) {
                keep_vertex[idx] = false;
                continue;
            }
            std::string vert_as_num = vert_to_num[idx];
            if (num_to_verts[vert_as_num].size() == 1) {
                continue;
            }
            int max_other_next_dist = -time_ub;
            int min_other_current_dist = time_ub;
            for (const int i : num_to_verts[vert_as_num]) {
                if (i == idx) continue;
                max_other_next_dist = std::max(max_other_next_dist, max_next_distance[i]);
                if (demand_verts[k][idx] == size) {
                    min_other_current_dist = std::min(min_other_current_dist,
                        min_dist_prev[k][i] + cust_dist[veh_vertices[k][i][0]][0]);
                }
            }
        }

        const int no_verts = std::count(keep_vertex.begin(), keep_vertex.end(), true);
        std::vector<std::vector<int>> temp_vertices(no_verts);
        std::vector<int> temp_no_stops(no_verts);
        std::vector<int> temp_demand_verts(no_verts);
        std::vector<int> temp_min_times(no_verts);
        std::vector<int> temp_min_dist_prev(no_verts);
        std::vector<int> temp_min_costs(no_verts);
        std::vector<int> temp_vert_costs(no_verts);
        int no_stop = 0;
        int i = 0;
        for (int idx = 0; idx < veh_vertices[k].size(); idx++) {
            if (no_stops[k][idx] > no_stop) {
                no_stop++;
                no_stops_indices[k][no_stop] = veh_vertices[k].size();
            }
            if (keep_vertex[idx]) {
                temp_vertices[i] = veh_vertices[k][idx];
                temp_no_stops[i] = no_stops[k][idx];
                temp_demand_verts[i] = demand_verts[k][idx];
                temp_min_times[i] = min_times[k][idx];
                temp_min_dist_prev[i] = min_dist_prev[k][idx];
                temp_min_costs[i] = min_costs[k][idx];
                temp_vert_costs[i] = vert_costs[k][idx];
                i++;
            }
        }
        veh_vertices[k] = temp_vertices;
        no_stops[k] = temp_no_stops;
        demand_verts[k] = temp_demand_verts;
        min_times[k] = temp_min_times;
        min_dist_prev[k] = temp_min_dist_prev;
        min_costs[k] = temp_min_costs;
        vert_costs[k] = temp_vert_costs;

        int idx_stop = 0;
        // Change the index for when no of stops changes
        for (int i = 1; i <= capacity[k] - veh_occ[k] + 1; i++) {
            bool changed = false;
            for (int idx = idx_stop; idx < veh_vertices[k].size(); idx++) {
                if (no_stops[k][idx] == i) {
                    idx_stop = idx;
                    no_stops_indices[k][i] = idx_stop;
                    changed = true;
                    break;
                }
            }
            if (!changed) no_stops_indices[k][i] = veh_vertices[k].size();
        }
    }

    void worker_remaining(std::atomic<int>& task_counter, std::vector<std::vector<int>>& min_costs, std::vector<std::vector<int>>& vert_costs,
        std::vector<std::vector<int>>& min_dist_prev, std::vector<std::vector<int>>& no_stops_indices,
        std::vector<std::vector<std::vector<int>>>& min_costs_n_step) {
        while (true) {
            int task = task_counter.fetch_add(1); // Get the next task and increment the counter
            if (task >= no_veh) break;        // Stop if all tasks are processed
            remaining_vertices(task, min_costs, vert_costs, min_dist_prev, no_stops_indices, min_costs_n_step);               // Process the task
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
            int k = 0;
            while (idx >= first_route_index[k + 1]) k++;
            const int route_length = route_stops[idx];
            const int first_cust = veh_vertices[k][routes[idx][0]][0];
            sol_values[first_cust][0] += val;
            if (first_cust == 0) continue;
            const int last_cust = veh_vertices[k][routes[idx][route_length - 1]][0];
            sol_values[last_cust][0] += val;
            for (int i = 0; i < route_length - 1; i++) {
                sol_values
                    [veh_vertices[k][routes[idx][i]][0]]
                    [veh_vertices[k][routes[idx][i + 1]][0]]
                    += val;
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
                    if (route_stops[idx] == 1) continue;
                    int k = 0;
                    while (idx >= first_route_index[k + 1]) k++;
                    int cust_i;
                    for (cust_i = 0; cust_i < route_stops[idx]; cust_i++) {
                        if (veh_vertices[k][routes[idx][cust_i]][0] == i) break;
                    }
                    if (cust_i == route_stops[idx] - 1) continue;
                    int j = veh_vertices[k][routes[idx][cust_i + 1]][0];
                    if (cut.set.find(j) != cut.set.end()) lexpr += vars[idx];
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

    // Edges
    std::unordered_set<int> var_indices;
    std::unordered_set<int> var_indices_subset;

    // Graph Data
    // Graph Data
    std::vector<std::vector<int>> vertices;
    std::vector<std::vector<std::vector<int>>> veh_vertices;
    std::vector<std::unordered_map<int, int>> vertex_mapping;
    std::vector<std::pair<int, int>> edges;
    std::vector<int> first_var_index;
    std::vector<std::vector<int>> current_costs;
    std::vector<std::vector<int>> no_stops;
    std::vector<std::vector<int>> demand_verts;
    std::vector<std::vector<int>> min_times;

    std::vector<std::vector<int>> routes;
    std::vector<int> route_cost;
    std::vector<int> route_demand;
    std::vector<int> route_dist;
    std::vector<int> route_stops;
    std::vector<int> route_min_time;
    std::vector<int> first_route_index;
    std::vector<std::unordered_set<int>> route_cust_idx;


    // Timings
    double total_lp_solve_time = 0;
    double lp_solve_time = 0;
    double lp_rcc_solve_time = 0;
    double mip_solve_time = 0;
    double model_construction_time = 0;
    double graph_construction_time = 0;
    double RC_separation_time = 0;

    // Determine vertices
    void determine_graph() {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<int>> min_dist_prev(no_veh);
        std::vector<std::vector<int>> no_stops_indices(no_veh + 2);
        veh_vertices.resize(no_veh);
        no_stops.resize(no_veh);
        demand_verts.resize(no_veh);
        min_times.resize(no_veh);
        int cap_magn = ceil(log10(max_cap));
        int max_dist_to_dest = *std::ranges::max_element(cust_dist[0]);
        std::vector<std::vector<std::vector<int>>> min_costs_n_step(no_veh);
        std::vector<std::vector<int>> min_costs(no_veh);
        std::vector<std::vector<int>> vert_costs(no_veh);


        std::cout << "Constructing Vertices..." << std::endl;

        std::atomic<int> task_counter(0); // Atomic counter to assign tasks
        std::vector<std::thread> threads;
        for (int i = 0; i < std::min(num_threads, no_veh); i++) {
            threads.emplace_back([this, &task_counter, &min_dist_prev, &no_stops_indices, &min_costs, &vert_costs, &min_costs_n_step]() {
                this->worker_remaining(task_counter, min_costs, vert_costs, min_dist_prev, no_stops_indices, min_costs_n_step);  // Call member function on the instance
                });
        }

        // Join threads
        for (auto& t : threads) {
            if (t.joinable()) {
                t.join();
            }
        }

        //int idx = 0;
        //for (const int& k : vehicles) {
        //    for (int j = 0; j < veh_vertices[k].size(); j++) {
        //        auto end = veh_vertices[k][j].end();
        //            std::cout << idx << " " << k << " " << j << ": (";
        //            for (int i : veh_vertices[k][j]) {
        //                std::cout << " " << i << ",";
        //            }
        //            std::cout << ")"  << std::endl;
        //            idx++;
        //    }
        //}

        //std::cout << "# Vertices: " << vertices.size() << std::endl;

        std::cout << "Constructing Routes..." << std::endl;

        
        std::set<int> routes_to_remove;
        std::vector<int> best_child_cost;
        // Determine possible routes pr vehicle
        for (const int& k : vehicles) {
            int size = capacity[k] - veh_occ[k];
            std::vector<std::set<int>> possible_next_events(veh_vertices[k].size());
            std::unordered_map<int, int> index_mapping;
            for (int idx = 1; idx < no_stops_indices[k][std::min(int(size), 2)]; idx++) {
                index_mapping[veh_vertices[k][idx][0]] = idx;
            }
            if (size >= 2) {
                for (int idx = no_stops_indices[k][2]; idx < no_stops_indices[k][size + 1]; idx++) {
                    int index = 0;
                    int mult = 1;
                    for (const int& i : veh_vertices[k][idx]) {
                        index += mult * i;
                        mult *= pow(10, cust_magn);
                    }
                    index_mapping[index] = idx;
                    for (int i = 1; i < no_stops[k][idx]; i++) {
                        int index = veh_vertices[k][idx][i];
                        int mult = pow(10, cust_magn);
                        for (int j = 1; j < no_stops[k][idx]; j++) {
                            if (i == j) continue;
                            index += mult * veh_vertices[k][idx][j];
                            mult *= pow(10, cust_magn);
                        }
                        if (index_mapping.find(index) != index_mapping.end()) {
                            possible_next_events[index_mapping[index]].insert(idx);
                        }
                    }
                }
            }



            const int min_profitable_cost = (veh_occ[k] == 0) ? 0 : veh_dist[k][0];
            std::vector<int> stop_index;
            stop_index.emplace_back(routes.size());
            first_route_index.emplace_back(routes.size());
            // Add route directly to depot
            std::vector<int> depot = { 0 };
            routes.emplace_back(depot);
            route_demand.emplace_back(0);
            route_dist.emplace_back(veh_dist[k][0]);
            route_stops.emplace_back(0);
            route_min_time.emplace_back(veh_arr_time[k]);
            route_cost.emplace_back(veh_dist[k][0]);
            best_child_cost.emplace_back(veh_dist[k][0]);
            if (size == 0) continue;
            // Determine possible routes of 1 stop
            stop_index.emplace_back(routes.size());
            for (int idx = no_stops_indices[k][1]; idx < no_stops_indices[k][2]; idx++) {
                int i = veh_vertices[k][idx][0];
                std::vector<int> route;
                int dist = veh_dist[k][i] + cust_dist[i][0];
                int min_time = min_times[k][idx];
                if (demand[i] > size || dist > min_time) continue;
                route.emplace_back(idx);
                routes.emplace_back(route);
                route_demand.emplace_back(demand[i]);
                route_dist.emplace_back(dist - cust_dist[i][0]);
                route_stops.emplace_back(1);
                route_min_time.emplace_back(min_time);
                route_cost.emplace_back(dist - demand[i] * cust_dist[i][0]);
                best_child_cost.emplace_back(dist - demand[i] * cust_dist[i][0]);
                if (route_cost[routes.size() - 1] > min_profitable_cost) {
                    routes_to_remove.insert(routes.size() - 1);
                }
            }
            stop_index.emplace_back(routes.size());

            // Determine rest of possible routes
            for (int l = 2; l <= size; l++) {
                std::unordered_map<std::string, int> ordered_best_costs;
                std::unordered_map<std::string, std::unordered_set<int>> index_to_routes;
                for (int m = stop_index[l - 1]; m < stop_index[l]; m++) {
                   for (const int& idx : possible_next_events[routes[m][route_stops[m] - 1]]) {
                        int i = veh_vertices[k][idx][0];
                        int dem = route_demand[m] + demand[i];

                        int prev_cust = veh_vertices[k][routes[m][route_stops[m] - 1]][0];
                        int dist = route_dist[m] + cust_dist[prev_cust][i] + cust_dist[i][0];
                        int min_time = std::min(route_min_time[m], cust_arr_time[i]);
                        if (dem > size || dist > min_time) continue;

                        std::vector<int> route = routes[m];
                        route.emplace_back(idx);
                        routes.emplace_back(route);
                        route_demand.emplace_back(dem);
                        route_dist.emplace_back(dist - cust_dist[i][0]);
                        route_stops.emplace_back(l);
                        route_min_time.emplace_back(min_time);
                        int cost = dist - demand[i] * cust_dist[i][0];
                        for (const int& j : routes[m]) {
                            cost -= demand[veh_vertices[k][j][0]] * cust_dist[veh_vertices[k][j][0]][0];
                        }
                        route_cost.emplace_back(cost);

                        std::vector<int> cust_sorted;
                        for (const int& i : veh_vertices[k][route[route_stops[m]]]) {
                            cust_sorted.emplace_back(i);
                        }
                        std::sort(cust_sorted.begin(), cust_sorted.end());

                        std::string index = vert_to_str(cust_sorted);
                        int mult = 1;
                        bool exists = false;
                        if (cost > best_child_cost[m]) {
                            exists = true;
                            routes_to_remove.insert(routes.size() - 1);
                            if (ordered_best_costs.find(index) == ordered_best_costs.end()) ordered_best_costs[index] = best_child_cost[m];
                            else ordered_best_costs[index] = std::min(ordered_best_costs[index], best_child_cost[m]);
                        }
                        else if (cost > min_profitable_cost) routes_to_remove.insert(routes.size() - 1);

                        // If routes with same costumers in different order
                        if (!exists && ordered_best_costs.find(index) == ordered_best_costs.end()) ordered_best_costs[index] = cost;
                        // Remove if better route with same set of customers
                        else if (ordered_best_costs[index] < cost) routes_to_remove.insert(routes.size() - 1);
                        else if (ordered_best_costs[index] > cost) {
                            for (const int& r : index_to_routes[index]) {
                                routes_to_remove.insert(r);
                            }
                            // Update Cost
                            ordered_best_costs[index] = cost;
                        }
                        index_to_routes[index].insert(routes.size() - 1);
                        // Update best_child
                        best_child_cost.emplace_back(std::min(cost, best_child_cost[m]));
                    }
                }
                stop_index.emplace_back(routes.size());
            }
        }
        first_route_index.emplace_back(routes.size());

        {
            std::vector<std::vector<int>> temp_routes;
            std::vector<int> temp_route_stops;
            std::vector<int> temp_first_route_index;
            std::vector<int> temp_route_costs;
            for (const int& k : vehicles) {
                temp_first_route_index.emplace_back(temp_routes.size());
                for (int idx = first_route_index[k]; idx < first_route_index[k + 1]; idx++) {
                    if (routes_to_remove.find(idx) != routes_to_remove.end()) continue;
                    temp_routes.emplace_back(routes[idx]);
                    temp_route_stops.emplace_back(route_stops[idx]);
                    temp_route_costs.emplace_back(route_cost[idx]);
                }
            }
            temp_first_route_index.emplace_back(temp_routes.size());
            routes = temp_routes;
            route_stops = temp_route_stops;
            first_route_index = temp_first_route_index;
            route_cost = temp_route_costs;
            route_demand.clear();
            route_dist.clear();
            route_min_time.clear();
        }

        for (int idx = 0; idx < routes.size(); idx++) {
            var_indices.insert(var_indices.end(), idx);
        }
/*        for (int j = 0; j < routes.size();j++) {
            std::cout << j << ": (";
            for (int i : routes[j]) {
                std::cout << " " << i << ",";
            }
            std::cout << ") " << route_cost[j] << " " << route_stops[j] << std::endl;
        }*/ 
        std::cout << "# Routes: " << routes.size() << std::endl;

        //for (int idx = 0; idx < var_indices.size(); idx++) {
        //    int k = 0;
        //    while (idx >= first_route_index[k + 1]) k++;
        //        std::cout << k << ": ";
        //        for (const int& j : routes[idx]) {
        //            std::cout << veh_vertices[k][j][0] << "->";
        //        }
        //        std::cout << "d" << std::endl;
        //    }
        //}

        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Initializes instance
    REB(const std::string& filename, const int& threads, const int& time) {
        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        capacity = inst_data.capacity;
        no_veh = inst_data.no_veh;
        no_cust = inst_data.no_cust;
        no_cust_d = no_cust + 1;
        time_ub = inst_data.time_ub;
        num_threads = threads;
        time_limit = time;

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

        determine_graph();
    }

    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Get subset indices
        std::cout << "Building Model..." << std::endl;
        if (var_indices_subset.size() == 0) var_indices_subset = var_indices;

        // Determine routes of customer i
        route_cust_idx.resize(no_cust_d);
        for (const int& idx : var_indices_subset) {
            if (route_stops[idx] == 0) continue;
            int k = 0;
            while (idx >= first_route_index[k + 1]) k++;
            for (const int& i : veh_vertices[k][routes[idx][route_stops[idx] - 1]]) {
                route_cust_idx[i].insert(idx);
            }
        }

        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);
        model->getEnv().set(GRB_IntParam_Threads, num_threads);

        // Create Variables
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;
        vars.resize(var_indices.size());

        {
            const int size = int(var_indices.size());
            double* ub = new double[size];
            char* vtypes = new char[size];
            double* objs = new double[size];
            std::string* names = new std::string[size];

            for (int idx = 0; idx < var_indices.size(); idx++) {
                ub[idx] = 1.0;
                objs[idx] = -route_cost[idx];
                vtypes[idx] = vtype;
                names[idx] = "lamda[" + itos(idx) + "]";
            }
            GRBVar* var_list = model->addVars(NULL, ub, objs, vtypes, names, size);

            vars.assign(var_list, var_list + size);
        }

        // Add constraints
        GRBLinExpr lexpr, lexpr2, lexpr3;

        // RB.2
        for (const int& i : customers) {
            for (const int& idx : route_cust_idx[i]) {
                lexpr += vars[idx];
            }
            model->addConstr(lexpr <= 1);
            lexpr.clear();
        }

        // RB.3/4
        for (const int& k : vehicles) {
            for (int idx = first_route_index[k]; idx < first_route_index[k + 1]; idx++) {
                if (var_indices_subset.find(idx) != var_indices_subset.end()) lexpr += vars[idx];
            }
            if (veh_occ[k] == 0) model->addConstr(lexpr <= 1);
            else model->addConstr(lexpr == 1);
            lexpr.clear();
        }

        model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

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

    void solve_root_relaxation(bool separate_rci = true, bool silent = true, bool _ = true) {
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

    void solve_to_optimality(bool separate_rci = true, bool silent = true, bool _ = true) {
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
            for (int r = first_route_index[k]; r < first_route_index[k + 1]; r++) {
                for (const int& e : routes[r]) {
                    std::cout << "(";
                    for (const int& i : veh_vertices[k][e]) {
                        std::cout << i << ",";
                    }
                    std::cout << ")->";
                }
                std::cout << "d" << std::endl;

            }
        }
    }
};