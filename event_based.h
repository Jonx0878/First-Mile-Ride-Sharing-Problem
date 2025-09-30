#pragma once
#include <any>
#include <algorithm>
#include <execution>
#include <gurobi_c++.h>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <functional>

#include "load.h"
#include "rc.h"
#include "utils.h"


class EB {
private:
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
            if (demand[i] > size || veh_dist[k][i] + shortest_path_to_dest[i] > min_time) continue;
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
                    if (demand_sum > size || min_dist_prev_vert + shortest_path_to_dest[i] > min_time) continue;
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
            }
        }

        // Remove vertices that cannot be optimal choices due to a better order
        const int min_profitable_cost = (veh_occ[k] == 0) ? 0 : veh_dist[k][0];
        std::vector<bool> keep_vertex(veh_vertices[k].size());
        keep_vertex[0] = true;
        for (int idx = 1; idx < veh_vertices[k].size(); idx++) {
            min_costs[k][idx] = *std::ranges::min_element(min_costs_n_step[k][idx]);

            keep_vertex[idx] = true;

            // Event cannot be part of a profitable route
            if (min_costs[k][idx] > min_profitable_cost) {
                keep_vertex[idx] = false;
                continue;
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
            const int i = vertices[edges[idx].first][0];
            const int j = vertices[edges[idx].second][0];
            if (i == 0 && j == 0) continue;
            else if (i == 0) sol_values[j][i] += vars[idx].get(GRB_DoubleAttr_X);
            else sol_values[i][j] += vars[idx].get(GRB_DoubleAttr_X);
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
                for (const int& k : vehicles) {
                    for (const int& idx : edges_in_cust[k][i]) {
                        if (cut.set.find(vertices[edges[idx].first][0]) != cut.set.end()) {
                            lexpr += vars[idx];
                        }
                    }
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
    std::string time_constr;

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
    std::vector<int> profit;
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
    std::unordered_map<int, GRBVar> tau;
    std::unordered_map<int, GRBVar> Gamma;
    std::unordered_map<int, GRBVar> sigma;
    std::unordered_map<int, GRBVar> omega;
    std::unordered_set<int> tau_indices;
    std::unordered_set<int> sigma_indices;
    std::unordered_set<int> omega_indices;

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


    std::vector<std::vector<std::unordered_set<int>>> edges_out;
    std::vector<std::vector<std::unordered_set<int>>> edges_in;
    std::vector<std::vector<std::unordered_set<int>>> edges_out_cust;
    std::vector<std::vector<std::unordered_set<int>>> edges_in_cust;

    std::vector<std::vector<std::unordered_set<int>>> edges_i_to_j;
    std::vector<std::vector<std::vector<std::unordered_set<int>>>> edges_cust_k_i_to_j;


    // Shortest Path
    std::vector<std::vector<int>> shortest_path_from_veh;
    std::vector<int> shortest_path_to_dest;
    std::vector<std::vector<std::vector<int>>> u;


    // Timings
    double total_lp_solve_time = 0;
    double lp_solve_time = 0;
    double lp_rcc_solve_time = 0;
    double mip_solve_time = 0;
    double model_construction_time = 0;
    double graph_construction_time = 0;
    double RC_separation_time = 0;

    // Utils
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


    // Determine vertices
    void determine_graph() {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<int>> min_dist_prev(no_veh);
        veh_vertices.resize(no_veh);
        no_stops.resize(no_veh);
        demand_verts.resize(no_veh);
        min_times.resize(no_veh);
        int cap_magn = ceil(log10(max_cap));
        int max_dist_to_dest = *std::ranges::max_element(cust_dist[0]);
        std::vector<std::vector<std::vector<int>>> min_costs_n_step(no_veh);
        std::vector<std::vector<int>> min_costs(no_veh);
        std::vector<std::vector<int>> vert_costs(no_veh);
        std::vector<std::vector<int>> no_stops_indices(no_veh + 2);


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

        // Construct edges pr. vehicle
        first_var_index.emplace_back(edges.size());
        for (const int& k : vehicles) {
            std::vector<int> edges_in(veh_vertices[k].size());
            std::vector<int> edges_out(veh_vertices[k].size());
            for (int i = 0; i < veh_vertices[k].size(); i++) {
                edges_in[i] = 0;
                edges_out[i] = 0;
            }
            const int size = capacity[k] - veh_occ[k];
            int min_profitable_dist = 0;
            if (veh_occ[k] > 0) {
                min_profitable_dist = veh_dist[k][0];
                edges.emplace_back(std::make_pair(0, 0));
                edges_in[0]++;
                edges_out[0]++;
            }
            for (int i = 1; i < veh_vertices[k].size(); i++) {
                if (min_costs[k][i] > min_profitable_dist) continue;
                const int i_cust = veh_vertices[k][i][0];
                if (no_stops[k][i] == 1) {
                    edges.emplace_back(std::make_pair(0, i));
                    edges_in[i]++;
                    edges_out[0]++;
                }
                if (edges_in[i] == 0) continue;
                if (vert_costs[k][i] <= min_profitable_dist) {
                    edges.emplace_back(std::make_pair(i, 0));
                    edges_in[0]++;
                    edges_out[i]++;
                }
                if (vert_costs[k][i] == min_costs[k][i] || demand_verts[k][i] == size) continue;
                for (int j = no_stops_indices[k][no_stops[k][i] + 1]; j < no_stops_indices[k][no_stops[k][i] + 2]; j++) {
                    const int j_cust = veh_vertices[k][j][0];
                    int same_cust_counter = 0;
                    for (const int& l : veh_vertices[k][i]) {
                        if (std::find(veh_vertices[k][j].begin() + 1, veh_vertices[k][j].end(), l) != veh_vertices[k][j].end()) same_cust_counter++;
                        else break;
                    }
                    if (same_cust_counter < no_stops[k][i]) continue;
                    const int min_total_dist = min_dist_prev[k][i] + cust_dist[i_cust][j_cust] + shortest_path_to_dest[j_cust];
                    if ((min_dist_prev[k][i] + cust_dist[i_cust][j_cust] == min_dist_prev[k][j]) && (vert_costs[k][i] > min_costs[k][j]) && (min_total_dist <= min_times[k][j])) {
                        edges.emplace_back(std::make_pair(i, j));
                        edges_out[i]++;
                        edges_in[j]++;
                    }
                }
            }
            while (true) {
                int edges_removed = 0;
                for (int idx = edges.size() - 1; idx >= first_var_index[k]; idx--) {
                    int i = edges[idx].first;
                    int j = edges[idx].second;
                    if ((edges_in[i] == 0) || (edges_out[j] == 0)) {
                        edges_removed++;
                        edges_out[i]--;
                        edges_in[j]--;
                        edges.erase(edges.begin() + idx);
                    }
                }
                if (edges_removed == 0) break;
            }
            first_var_index.emplace_back(edges.size());


            std::vector<int> vertex_map;
             
            std::vector<std::vector<int>> temp_vertices;
            for (int idx = 0; idx < veh_vertices[k].size(); idx++) {
                vertex_map.emplace_back(temp_vertices.size());
                if (edges_in[idx] > 0 && edges_out[idx] > 0) {
                    temp_vertices.emplace_back(veh_vertices[k][idx]);
                }
            }

            veh_vertices[k] = temp_vertices;
            for (int i = first_var_index[k]; i < first_var_index[k + 1]; i++) {
                edges[i].first = vertex_map[edges[i].first];
                edges[i].second = vertex_map[edges[i].second];
            }

            first_var_index[k + 1] = edges.size();
            int idx_stop = 0;
            for (int i = 1; i <= capacity[k] - veh_occ[k] + 1; i++) {
                bool changed = false;
                for (int idx = idx_stop; idx < veh_vertices[k].size(); idx++) {
                    if (no_stops[k][idx] == i) {
                        idx_stop = idx;
                        no_stops_indices[k][i] = idx_stop;
                        changed = true;
                        break;
                    }
                    if (!changed) no_stops_indices[k][i] = veh_vertices[k].size();
                }
            }

        }
        vertices = veh_vertices[0];
        std::unordered_map<std::string, int> num_to_index;
        vertex_mapping.resize(no_veh);

        int idx = 0;
        for (const int& k : vehicles) {
            for (int j = 0; j < veh_vertices[k].size(); j++) {
                std::string index = vert_to_str(veh_vertices[k][j]);
                if (k == 0) {
                    num_to_index[index] = idx;
                    vertex_mapping[k][j] = idx;
                    idx++;
                }
                else if (num_to_index.find(index) != num_to_index.end()) {
                    vertex_mapping[k][j] = num_to_index[index];
                }
                else {
                    num_to_index[index] = vertices.size();
                    vertex_mapping[k][j] = vertices.size();
                    vertices.emplace_back(veh_vertices[k][j]);
                }

            }
        }
        //for (int j = 0; j < vertices.size(); j++) {
        //    std::cout << j << ": (";
        //    for (int i : vertices[j]) {
        //        std::cout << " " << i << ",";
        //    }
        //    std::cout << ") " << std::endl;
        //}
        std::cout << "# Vertices: " << vertices.size() << std::endl;

        for (const int& k : vehicles) {
            for (int idx = first_var_index[k]; idx < first_var_index[k + 1]; idx++) {
                edges[idx].first = vertex_mapping[k][edges[idx].first];
                edges[idx].second = vertex_mapping[k][edges[idx].second];
            }
        }
        std::cout << "#Edges: " << edges.size() << std::endl;

        for (int idx = 0; idx < edges.size(); idx++) {
            var_indices.insert(idx);
        }

        current_costs.clear();
        no_stops.clear();
        demand_verts.clear();
        min_times.clear();
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Initializes instance
    EB(const std::string& filename, const std::string time_constr, const int& threads, const int& time) {
        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        capacity = inst_data.capacity;
        no_veh = inst_data.no_veh;
        no_cust = inst_data.no_cust;
        no_cust_d = no_cust + 1;
        time_ub = inst_data.time_ub;
        this->time_constr = time_constr;
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

        max_cap = 0;
        for (const int& k : vehicles) {
            max_cap = std::max(capacity[k] - veh_occ[k], max_cap);
        }

        // Shortest path to destination
        shortest_path_to_dest.resize(no_cust_d);
        for (const int& i : customers) {
            shortest_path_to_dest[i] = cust_dist[i][0];
        }
        while (true) {
            bool any_change = false;
            for (const int& i : customers) {
                for (const int& j : customers) {
                    if (i == j) continue;
                    if (shortest_path_to_dest[i] > cust_dist[i][j] + shortest_path_to_dest[j]) {
                        shortest_path_to_dest[i] = cust_dist[i][j] + shortest_path_to_dest[j];
                        any_change = true;
                    }
                }
            }
            if (!any_change) break;
        }

        determine_graph();

		// Determine shortest constants for time constraints
        if (time_constr == "TTCF") {
            // Shortest path to customer
            shortest_path_from_veh.resize(no_veh);
            for (const int& k : vehicles) {
                shortest_path_from_veh[k].resize(no_cust_d);
                for (const int& i : dest_and_cust) {
                    shortest_path_from_veh[k][i] = veh_dist[k][i];
                    //std::cout << k << "->" << i << " " << veh_dist[k][i] << std::endl;
                }
            }
            while (true) {
                bool any_change = false;
                for (const int& k : vehicles) {
                    for (const int& i : customers) {
                        for (const int& j : dest_and_cust) {
                            if (i == j) continue;
                            if (shortest_path_from_veh[k][j] > cust_dist[i][j] + shortest_path_from_veh[k][i]) {
                                shortest_path_from_veh[k][j] = cust_dist[i][j] + shortest_path_from_veh[k][i];
                                any_change = true;
                            }
                        }
                    }
                }
                if (!any_change) break;
            }


            // u
            u.resize(no_veh);
            for (const int& k : vehicles) {
                u[k].resize(no_cust_d);
                for (const int& i : customers) {
                    u[k][i].resize(no_cust_d);
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        u[k][i][j] = std::min(std::min(cust_arr_time[i], cust_arr_time[j]), veh_arr_time[k])
                            - cust_dist[i][j] - shortest_path_to_dest[j];
                    }
                }
            }

            //for (const int& k : vehicles) {
            //    for (const int& i : customers) {
            //        for (const int& j : dest_and_cust) {
            //            if (i == j) continue;
            //                if (u[k][i][j] < shortest_path_from_veh[k][i] - 1e-9) {
            //                    std::cerr << "Violation: u(" << k << "," << i << "," << j << ")="
            //                        << u[k][i][j]
            //                        << " < sp(" << k << "," << i << ")=" << shortest_path_from_veh[k][i]
            //                        << "\n";
            //                }
            //                
            //            
            //        }
            //    }
            //}
        }

        // Determine profits
        profit.resize(no_cust_d);
        for (const int& i : customers) {
            profit[i] = demand[i] * cust_dist[i][0];
        }
    }

    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        // Get subset indices
        std::cout << "Building Model..." << std::endl;

        edges_out.resize(no_veh);
        edges_in.resize(no_veh);
        edges_out_cust.resize(no_veh);
        edges_in_cust.resize(no_veh);
        edges_i_to_j.resize(no_cust_d);
        edges_cust_k_i_to_j.resize(no_veh);

        for (const int& k : vehicles) {
            edges_cust_k_i_to_j[k].resize(no_cust_d);
            for (const int& i : dest_and_cust) {
                edges_cust_k_i_to_j[k][i].resize(no_cust_d);
            }
        }
        for (const int& i : dest_and_cust) edges_i_to_j[i].resize(no_cust_d);

        for (const int& k : vehicles) {
            edges_out[k].resize(vertices.size() + 1);
            edges_in[k].resize(vertices.size() + 1);
            edges_out_cust[k].resize(no_cust_d);
            edges_in_cust[k].resize(no_cust_d);
            for (int idx = first_var_index[k]; idx < first_var_index[k + 1]; idx++) {
                const int i = edges[idx].first;
                const int j = edges[idx].second;
                const int i_cust = vertices[i][0];
                const int j_cust = vertices[j][0];
                edges_out[k][i].insert(idx);
                edges_in[k][j].insert(idx);
                edges_out_cust[k][i_cust].insert(idx);
                edges_in_cust[k][j_cust].insert(idx);
                edges_i_to_j[i_cust][j_cust].insert(idx);
                edges_cust_k_i_to_j[k][i_cust][j_cust].insert(idx);
            }
        }

        // Create Model
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);
        model->getEnv().set(GRB_IntParam_Threads, num_threads);

        // Create Variables
        double obj;
        char vtype = relax ? GRB_CONTINUOUS : GRB_BINARY;;
        std::string name;
        vars.resize(var_indices.size());
        for (const auto& idx : var_indices){
            int k = 0;
            while (idx >= first_var_index[k + 1]) k++;
            int i = edges[idx].first;
            int i_cust = vertices[i][0];
            int j = edges[idx].second;
            int j_cust = vertices[j][0];
            if (i == 0) {
                obj = -veh_dist[k][j_cust];
                name = "gamma[" + itos(k) + "][" + itos(j) + "]";
            }
            else {
                obj = profit[i_cust] - cust_dist[i_cust][j_cust];
                name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
            }
            vars[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
        }

        // Tau
        for (const int& i : customers) {
            bool need_variable = false;
            for (const int& k : vehicles) {
                if (!edges_in_cust[k][i].empty() && !edges_out_cust[k][i].empty()) {
                    need_variable = true;
                    break;
                }
            }
            if (need_variable) {
                tau_indices.insert(i);
				const std::string name = "tau[" + itos(i) + "]";
                tau[i] = model->addVar(0, 1, 0, GRB_CONTINUOUS, name);
            }
        }

        if (time_constr == "TIME") {
            // Gamma
            for (const int& k : vehicles) {
                const std::string name = "Gamma[" + itos(k) + "]";
                Gamma[k] = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
            }
        }
        else if (time_constr == "TTCF") {
			// Construct Indices for sigma and omega
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
  
                    omega_indices.insert(i * no_cust_d + j);
                    sigma_indices.insert(i* no_cust_d + j);
                }
            }

            // Sigma
            for (const int& idx : sigma_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                const std::string name = "sigma[" + itos(i) + "][" + itos(j) + "]";
                if (j > 0) sigma[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
                else sigma[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
            }

            // Omega
            for (const int& idx : omega_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                const std::string name = "omega[" + itos(i) + "][" + itos(j) + "]";
                omega[idx] = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
            }
        }

        // Add constraints
        GRBLinExpr lexpr, lexpr2, lexpr3, lexpr4;
        
        // 2.2/2.3
        for (const int& k : vehicles) {
            for (const int& idx : edges_out[k][0]) lexpr += vars[idx]; // LHS
            if (veh_occ[k] == 0) model->addConstr(lexpr <= 1);
            else model->addConstr(lexpr == 1);
            lexpr.clear();
        }

        // 2.4
        for (const int& k : vehicles) {
            for (const int& idx : edges_out[k][0]) lexpr += vars[idx]; // LHS
            for (const int& idx : edges_in[k][0]) lexpr2 += vars[idx]; // RHS
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }


        // 2.5
        for (const int& i : customers) {
            // LHS
            for (const int& k : vehicles) {
                for (const int& idx : edges_out_cust[k][i]) lexpr += vars[idx];
            }
            // RHS
			if (tau_indices.contains(i)) lexpr2 += tau[i];
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();
        }

        // 2.6
        for (int j = 1; j < vertices.size(); j++) {
            for (const int& k : vehicles) {
                for (const int& idx : edges_in[k][j]) lexpr += vars[idx]; // LHS
                for (const int& idx : edges_out[k][j]) lexpr2 += vars[idx]; // RHS
                if (lexpr.size() + lexpr2.size() > 0) model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }

        if (time_constr == "TIME"){
            // TIME.1
            for (const int& k : vehicles) {
                // RHS
                for (const int& v : dest_and_cust) {
                    for (const int& idx : edges_out_cust[k][v]) {
                        int dist = 0;
                        int j = vertices[edges[idx].second][0];
                        if (v == 0) dist = veh_dist[k][j];
                        else dist = cust_dist[v][j];
                        lexpr += dist * vars[idx];
                    }
                }
                model->addConstr(Gamma[k] == lexpr);
                lexpr.clear();
            }

            // TIME.2
			int max_cust_arr = *std::ranges::max_element(cust_arr_time);
            for (const int& k : vehicles) {
                if (max_cust_arr >= veh_arr_time[k]) {
                    model->addConstr(Gamma[k] <= veh_arr_time[k]);
                }
                else {
                    for (const int& i : customers) {
                        if (veh_arr_time[k] <= cust_arr_time[i]) continue;
                        // RHS
                        for (const int& idx : edges_out_cust[k][i]) lexpr += vars[idx];
                        model->addConstr(
                            Gamma[k] <= veh_arr_time[k] + (cust_arr_time[i] - veh_arr_time[k]) * lexpr
                        );
                        lexpr.clear();
                    }
				}
            }            
        }
        else if (time_constr == "TTCF") {
            // TTCF.1
            for (const int& i : customers) {
                // LHS
                // First term
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (sigma_indices.contains(idx)) lexpr += sigma[idx];
                }
                // Second term
                for (const int& k : vehicles) {
                    for (const int& v_idx : edges_cust_k_i_to_j[k][0][i]) {
                        lexpr2 += veh_dist[k][i] * vars[v_idx];
                    }
                }
                for (const int& j : customers) {
                    if (i == j) continue;
                    for (const int& v_idx : edges_i_to_j[j][i]) {
                        lexpr2 += cust_dist[j][i] * vars[v_idx];
                    }
                }
                // RHS
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (sigma_indices.contains(idx)) lexpr3 += sigma[idx];
                }
                model->addConstr(lexpr + lexpr2 == lexpr3);
                lexpr.clear();
                lexpr2.clear();
                lexpr3.clear();
            }
            // TTCF.2
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (!sigma_indices.contains(idx)) continue;
                    // RHS
                    for (const int& k : vehicles) {
						for (const int& v_idx : edges_cust_k_i_to_j[k][i][j]) {
                            lexpr += (u[k][i][j]) * vars[v_idx];
                            lexpr2 += shortest_path_from_veh[k][i] * vars[v_idx];
                        }
                    }
                    model->addConstr(sigma[idx] <= lexpr);
                    //model->addConstr(sigma[idx] >= lexpr2);

                    lexpr.clear();
                    lexpr2.clear();
                }
            }
            // TTCF.3 (2.18)
            for (const int& i : customers) {
                // LHS
                for (const int& k : vehicles) {
                    for (const int& v_idx : edges_cust_k_i_to_j[k][0][i]) {
                        lexpr += veh_arr_time[k] * vars[v_idx];
                    }
                }
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (omega_indices.contains(idx)) lexpr += omega[idx];
                }
                // RHS
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (omega_indices.contains(idx)) lexpr2 += omega[idx];
                }
                model->addConstr(lexpr >= lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
            // TTCF.3 (2.19)
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (!omega_indices.contains(idx)) continue;
                    for (const int& k : vehicles) {
                        for (const int& v_idx : edges_cust_k_i_to_j[k][i][j]) {
                            lexpr += cust_arr_time[i] * vars[v_idx];
                        }
                    }
                    model->addConstr(omega[idx] <= lexpr);
                    lexpr.clear();
                }
            }
            // TTCF.4
            for (const int& i : customers) {
                // LHS
                const int idx = i * no_cust_d + 0;
                if (sigma_indices.contains(idx)) lexpr += sigma[idx];
                for (const int& v_idx : edges_i_to_j[i][0]) {
                    lexpr += cust_dist[i][0] * vars[v_idx];
                }
                // RHS
                if (omega_indices.contains(idx)) lexpr2 += omega[idx];
                model->addConstr(lexpr <= lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
            
        }
        model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
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
            bool cuts_added;
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

    void print_solution_routes() {
        for (const int& k : vehicles) {
            std::vector<int> route_k;
            int current_event = 0;
            for (int idx = first_var_index[k]; idx < first_var_index[k + 1]; idx++) {
                if (vars[idx].get(GRB_DoubleAttr_X) < 0.5) continue;
                if (edges[idx].first == current_event) {
                    current_event = edges[idx].second;
                    route_k.emplace_back(current_event);
                    if (current_event == 0) break;
                }
            }
            if (route_k.size() == 0) {
                std::cout << k << ": " << std::endl;
            }
            else {
                std::cout << k << ": " << 0 << "->";
                for (int i = 0; i < route_k.size() - 1; i++) {
                    std::cout << vertices[route_k[i]][0] << "->";
                }
                 std::cout << "d" << std::endl;
            }
        }
    }
};