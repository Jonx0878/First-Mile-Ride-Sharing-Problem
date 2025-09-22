#pragma once
#include <algorithm>
#include <gurobi_c++.h>
#include <set>
#include <unordered_set>

#include "load.h"
#include "rc.h"
#include "utils.h"


class FMRSP {
private:
    // Construct indices of variables
    void construct_indices() {
		// ----------------------------- Main Variables -----------------------------
        // Tau
        if (main_type == "3I" || main_type == "2I") {
            for (const int& i : customers) {
                tau_indices.insert(i);
            }
        }

        // Gamma
        if (main_type == "3I" ||  main_type == "2I") {
            for (const int& k : vehicles) {
                for (const int& i : dest_and_cust) {
                    gamma_indices.insert(k * no_cust_d + i);
                }
            }
        }

        // Chi
        if (main_type == "3I") {
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    for (const int j : dest_and_cust) {
                        if (i == j) continue;
                        chi_indices.insert(k * no_cust_d_squared + i * no_cust_d + j);
                    }
                }
            }
        }

        // Psi
        if (main_type == "2I") {
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    psi_indices.insert(i * no_cust_d + j);
                }
            }
        }

        // ----------------------------- Capacity Variables -----------------------------
        // Phi
        if (cap_type == "DL") {
            for (const int& i : customers) {
                phi_indices.insert(i);
            }
        }

        // Zeta
        if (cap_type == "COCF") {
            for (const int& i : customers) {
                for (const int j : customers) {
                    if (i == j) continue;
                    zeta_indices.insert(i * no_cust_d + j);
                }
            }
        }

        // Zeta Multi
        if (cap_type == "CMCF") {
            for (const int& l : customers) {
                for (const int& i : dest_and_cust) {
                    if (i == l) continue;
                    for (const int& j : customers) {
                        if (i == j) continue;
                        const int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        zeta_multi_indices.insert(idx);
                    }
                }
            }
        }

        // Rho Multi
        if (cap_type == "CMCF") {
            for (const int& l : customers) {
                for (const int& i : customers) {
                    for (const int& j : dest_and_cust) {
                        if (j == i || j == l) continue;
                        const int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        rho_multi_indices.insert(idx);
                    }
                }
            }
        }
        // ----------------------------- Time Variables -----------------------------
        // Gamma_Capital
        if (time_type == "TIME") {
            for (const int& k : vehicles) {
                Gamma_indices.insert(k);
            }
        }

        // Sigma
        if (time_type == "TTCF") {
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    sigma_indices.insert(i * no_cust_d + j);
                }
            }
        }

        // Omega
        if (time_type == "TTCF") {
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    omega_indices.insert(i * no_cust_d + j);
                }
            }
        }
    };

    // Eliminate Variables
    void preprocess() {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Eliminating Variables..." << std::endl;
        std::unordered_set<int> tau_indices_removed;
        std::unordered_set<int> gamma_indices_removed;
        std::unordered_set<int> chi_indices_removed;
        std::unordered_set<int> phi_indices_removed;
        std::unordered_set<int> psi_indices_removed;
        std::unordered_set<int> zeta_indices_removed;
        std::unordered_set<int> zeta_multi_indices_removed;
        std::unordered_set<int> rho_multi_indices_removed;
        std::unordered_set<int> sigma_indices_removed;
        std::unordered_set<int> omega_indices_removed;

		// Keep track of percentage of eliminated variables
        int tau_before = tau_indices.size();
		int gamma_before = gamma_indices.size();
		int chi_before = chi_indices.size();
		int psi_before = psi_indices.size();
		int phi_before = phi_indices.size();
		int zeta_before = zeta_indices.size();
		int zeta_multi_before = zeta_multi_indices.size();
		int rho_multi_before = rho_multi_indices.size();
		int sigma_before = sigma_indices.size();
		int omega_before = omega_indices.size();
        // Eliminate Variables
        while (true) {
            int total_removed = 0;
            // ----------------------------- Main Variables -----------------------------
            // Tau
            for (const int& i : tau_indices) {
                bool can_eliminate1 = true, can_eliminate2 = true;
                if (main_type == "3I") {
                    // 1.1
                    for (const int& k : vehicles) {
                        for (const int& j : dest_and_cust) {
							int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) {
                                can_eliminate1 = false;
                                break;
                            }
                        }
						if (!can_eliminate1) break;
                    }
                    // 2.1
                    for (const int& k : vehicles) {
                        if (gamma_indices.contains(k * no_cust_d + i)) {
                            can_eliminate2 = false;
                            break;
						}
                        for (const int& j : customers) {
                            int idx = k * no_cust_d_squared + j * no_cust_d + i;
                            if (chi_indices.contains(idx)) {
                                can_eliminate2 = false;
                                break;
                            }
                        }
                        if (!can_eliminate2) break;
                    }
                }
                else {
                    // 1.2
                    for (const int& j : dest_and_cust) {
                        int idx = i * no_cust_d + j;
                        if (psi_indices.contains(idx)) {
                            can_eliminate1 = false;
                            break;
                        }
                    }
                    // 2.2
                    for (const int& k : vehicles) {
                        if (gamma_indices.contains(k * no_cust_d + i)) {
                            can_eliminate2 = false;
                            break;
                        }
                    }
                    for (const int& j : customers) {
                        int idx = j * no_cust_d + i;
                        if (chi_indices.contains(idx)) {
                            can_eliminate2 = false;
                            break;
                        }
                    }
                }
                if (can_eliminate1 || can_eliminate2) {
                    tau_indices_removed.insert(i);
                    total_removed++;
				}
            }
            for (const int& i : tau_indices_removed) {
                tau_indices.erase(i);
			}
            tau_indices_removed.clear();

            // Gamma
            for (const int& idx : gamma_indices) {
                const int k = int(floor(idx / no_cust_d));
                const int j = idx - k * no_cust_d;
                bool can_eliminate = false;
                // 1
                if (j == 0 && veh_occ[k] == 0) can_eliminate = true;
                // 2
                else if (demand[j] > capacity[k] - veh_occ[k]) can_eliminate = true;
                // 3
                else if (std::min(veh_arr_time[k], cust_arr_time[j]) <
                    veh_dist[k][j] + shortest_path_to_dest[j]) can_eliminate = true;
                // 4
                else if (!tau_indices.contains(j) && j > 0) can_eliminate = true;
                
                if (can_eliminate) {
                    gamma_indices_removed.insert(idx);
                    total_removed++;
				}
            }
            for (const int& idx : gamma_indices_removed) {
                gamma_indices.erase(idx);
            }
            gamma_indices_removed.clear();

            // Chi
            for (const int& idx : chi_indices) {
                const int k = int(floor(idx / no_cust_d_squared));
                const int i = int(floor((idx - k * no_cust_d_squared) / no_cust_d));
                const int j = idx - k * no_cust_d_squared - i * no_cust_d;
                bool can_eliminate = false;
                // 1
                if (demand[i] + demand[j] > capacity[k] - veh_occ[k]) can_eliminate = true;
                // 2
                else if (std::min(veh_arr_time[k], std::min(cust_arr_time[i], cust_arr_time[j])) <
                    shortest_path_from_veh[k][i] + cust_dist[i][j] + shortest_path_to_dest[j]) can_eliminate = true;
                // 3
                else if (!tau_indices.contains(i) || (!tau_indices.contains(j) && j > 0)) can_eliminate = true;
                if (can_eliminate) {
                    chi_indices_removed.insert(idx);
                    total_removed++;
                }
            }
            for (const int& idx : chi_indices_removed) {
                chi_indices.erase(idx);
            }
            chi_indices_removed.clear();

            // Psi
            for (const int& idx : psi_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                bool can_eliminate = false;
                // 1
                if (demand[i] + demand[j] > max_cap) can_eliminate = true;
                // 3
                else if (!tau_indices.contains(i) || (!tau_indices.contains(j) && j > 0)) can_eliminate = true;
                // 2
                else {
                    int min_time = std::min(cust_arr_time[i], cust_arr_time[j]);
                    int max_time = 0;
					int shortest_path_to_i = std::numeric_limits<int>::max();
                    for (const int& k : vehicles) {
						max_time = std::max(veh_arr_time[k], max_time);
						shortest_path_to_i = std::min(shortest_path_to_i, shortest_path_from_veh[k][i]);
                    }
					min_time = std::min(min_time, max_time);
                    if (min_time < shortest_path_to_i + cust_dist[i][j] + shortest_path_to_dest[j]) can_eliminate = true;
                }
                if (can_eliminate) {
                    psi_indices_removed.insert(idx);
                    total_removed++;
                }
            }
            for (const int& idx : psi_indices_removed) {
                psi_indices.erase(idx);
            }
            psi_indices_removed.clear();

            // ----------------------------- Capacity Variables -----------------------------
            // Phi
            for (const int& i : phi_indices) {
                int can_eliminate = false;
                // 1
				if (!tau_indices.contains(i)) can_eliminate = true;
                if (can_eliminate) {
                    phi_indices_removed.insert(i);
                    total_removed++;
                }
            }
            for (const int& idx : phi_indices_removed) {
                phi_indices.erase(idx);
            }
            phi_indices_removed.clear();

            // Zeta
            for (const int& idx : zeta_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                bool can_eliminate = false;
                // 1
                if (max_cap - demand[i] - demand[j] <= 0) can_eliminate = true;
                else {
                    // 2.1
                    if (main_type == "3I") {
                        bool all_edges_eliminated = true;
                        for (const int& k : vehicles) {
                            int k_idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(k_idx)) {
                                all_edges_eliminated = false;
                                break;
							}
                        }
                        if (all_edges_eliminated) can_eliminate = true;
                    }
                    // 2.2
                    else {
                        if (!psi_indices.contains(idx)) can_eliminate = true;
                    }
                }
                if (can_eliminate) {
                    zeta_indices_removed.insert(idx);
                    total_removed++;
                }
            }
            for (const int& idx : zeta_indices_removed) {
                zeta_indices.erase(idx);
            }
            zeta_indices_removed.clear();

            // Zeta_Multi
            for (const int& idx : zeta_multi_indices) {
                const int l = int(floor(idx / no_cust_d_squared));
                const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
                const int j = idx - l * no_cust_d_squared - i * no_cust_d;
                bool can_eliminate = false;
                // i = d
                if (i == 0) {
                    // 1
					if (!tau_indices.contains(l)) can_eliminate = true;
                    // 2
                    else {
						bool all_edges_eliminated = true;
                        for (const int& k : vehicles) {
                            if (gamma_indices.contains(k * no_cust_d + j)) {
                                all_edges_eliminated = false;
                                break;
                            }
                        }
						if (all_edges_eliminated) can_eliminate = true;
                    }
                }
                else {
                    // 1
                    if (l !=j && demand[i] + demand[j] + demand[l] > max_cap) can_eliminate = true;
                    else {
                        // 2
						bool all_edges_eliminated = true;
                        for (const int& v : dest_and_cust) {
                            if (v == i || v == l) continue;
							int z_idx = l * no_cust_d_squared + v * no_cust_d + i;
                            if (zeta_multi_indices.contains(z_idx)) {
                                all_edges_eliminated = false;
                                break;
							}
                        }
						if (all_edges_eliminated) can_eliminate = true;
                        // 3
                        all_edges_eliminated = true;
                        for (const int& v : customers) {
                            if (v == l) continue;
                            int r_idx = i * no_cust_d_squared + v * no_cust_d + l;
                            if (rho_multi_indices.contains(r_idx)) {
                                all_edges_eliminated = false;
                                break;
                            }
                        }
                        if (all_edges_eliminated) can_eliminate = true;
                        // 4.1
                        if (main_type == "3I") {
                            bool all_edges_eliminated = true;
                            for (const int& k : vehicles) {
                                int k_idx = k * no_cust_d_squared + i * no_cust_d + j;
                                if (chi_indices.contains(k_idx)) {
                                    all_edges_eliminated = false;
                                    break;
                                }
                            }
                            if (all_edges_eliminated) can_eliminate = true;
                        }
                        // 4.2
                        else {
                            if (!psi_indices.contains(i * no_cust_d + j)) can_eliminate = true;
                        }
                    }
                }

                if (can_eliminate) {
                    zeta_multi_indices_removed.insert(idx);
                    total_removed++;
				}
            }
            for (const int& idx : zeta_multi_indices_removed) {
                zeta_multi_indices.erase(idx);
            }
            zeta_multi_indices_removed.clear();

            // Rho_Multi
            for (const int& idx : rho_multi_indices) {
                const int l = int(floor(idx / no_cust_d_squared));
                const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
                const int j = idx - l * no_cust_d_squared - i * no_cust_d;
                bool can_eliminate = false;
                // 1
                if (l != i && demand[i] + demand[j] + demand[l] > max_cap) can_eliminate = true;
                else if (j > 0) {
                    // 2
                    bool all_edges_eliminated = true;
                    for (const int& v : dest_and_cust) {
                        if (v == j || v == l) continue;
                        int r_idx = l * no_cust_d_squared + j * no_cust_d + v;
                        if (rho_multi_indices.contains(r_idx)) {
                            all_edges_eliminated = false;
                            break;
                        }
                    }
                    // 3
                    bool all_edges_eliminated2 = true;
                    for (const int& v : customers) {
                        if (v == l || v == j) continue;
                        int z_idx = j * no_cust_d_squared + v * no_cust_d + l;
                        if (zeta_multi_indices.contains(z_idx)) {
                            all_edges_eliminated2 = false;
                            break;
                        }
                    }
                    if (all_edges_eliminated || all_edges_eliminated2) {
                        //can_eliminate = true;
                    }
                }
                if (!can_eliminate) {
                    // 4.1
                    if (main_type == "3I") {
                        bool all_edges_eliminated = true;
                        for (const int& k : vehicles) {
                            int k_idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(k_idx)) {
                                all_edges_eliminated = false;
                                break;
                            }
                        }
                        if (all_edges_eliminated) can_eliminate = true;
                    }
                    // 4.2
                    else {
                        if (!psi_indices.contains(i * no_cust_d + j)) can_eliminate = true;
                    }
                }


                if (can_eliminate) {
                    rho_multi_indices_removed.insert(idx);
                    total_removed++;
                }
            }
            for (const int& idx : rho_multi_indices_removed) {
                rho_multi_indices.erase(idx);
            }
            rho_multi_indices_removed.clear();

            // ----------------------------- Time Variables -----------------------------
            // Sigma
            for (const int& idx : sigma_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                bool can_eliminate = false;
                // 1
                int u_ij = 0;
                int max_cust_time = std::min(cust_arr_time[i], cust_arr_time[j]);
				int shortest_path_to_i = std::numeric_limits<int>::max();
				int max_ub = std::numeric_limits<int>::min();
                for (const int& k : vehicles) {
					int u_ij_k = std::min(max_cust_time, veh_arr_time[k])
                        - cust_dist[i][j] - shortest_path_to_dest[j];
                    u_ij = std::max(u_ij, u_ij_k);
					shortest_path_to_i = std::min(shortest_path_to_i, shortest_path_from_veh[k][i]);
                    max_ub = std::max(max_ub, u_ij_k - shortest_path_from_veh[k][i]);
				}
                if (main_type == "3I" && max_ub <= 0) can_eliminate = true;
				else if (main_type == "2I" && u_ij - shortest_path_to_i <= 0) can_eliminate = true;
                // 2
                else {
                    // 2.1
                    if (main_type == "3I") {
						bool all_edges_eliminated = true;
                        for (const int& k : vehicles) {
                            int k_idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(k_idx)) {
                                all_edges_eliminated = false;
                                break;
                            }
						}
						if (all_edges_eliminated) can_eliminate = true;
                    }
                    // 2.2
                    else {
						if (!psi_indices.contains(idx)) can_eliminate = true;
                    }
                }

                if (can_eliminate) {
                    sigma_indices_removed.insert(idx);
                    total_removed++;
                }
            }
            for (const int& idx : sigma_indices_removed) {
                sigma_indices.erase(idx);
            }
            sigma_indices_removed.clear();

            // Omega
            for (const int& idx : omega_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                bool can_eliminate = false;
                // 1.1
                if (main_type == "3I") {
                    bool all_edges_eliminated = true;
                    for (const int& k : vehicles) {
                        int k_idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(k_idx)) {
                            all_edges_eliminated = false;
                            break;
                        }
                    }
                    if (all_edges_eliminated) can_eliminate = true;
                }
                // 1.2
                else {
                    if (!psi_indices.contains(idx)) can_eliminate = true;
                }
                if (can_eliminate) {
                    omega_indices_removed.insert(idx);
                    total_removed++;
                }
            }
            for (const int& idx : omega_indices_removed) {
                omega_indices.erase(idx);
            }
            omega_indices_removed.clear();
            if (total_removed == 0) break;
            break;
        }
        if (tau_before > 0) std::cout << "Tau variables eliminated: " << std::fixed << std::setprecision(2) << (tau_before - tau_indices.size()) / float(tau_before) * 100 << "%" << std::endl;
        if (gamma_before > 0) std::cout << "Gamma variables eliminated: " << std::fixed << std::setprecision(2) << (gamma_before - gamma_indices.size()) / float(gamma_before) * 100 << "%" << std::endl;
        if (chi_before > 0) std::cout << "Chi variables eliminated: " << std::fixed << std::setprecision(2) << (chi_before - chi_indices.size()) / float(chi_before) * 100 << "%" << std::endl;
        if (psi_before > 0) std::cout << "Psi variables eliminated: " << std::fixed << std::setprecision(2) << (psi_before - psi_indices.size()) / float(psi_before) * 100 << "%" << std::endl;
        if (phi_before > 0) std::cout << "Phi variables eliminated: " << std::fixed << std::setprecision(2) << (phi_before - phi_indices.size()) / float(phi_before) * 100 << "%" << std::endl;
        if (zeta_before > 0) std::cout << "Zeta variables eliminated: " << std::fixed << std::setprecision(2) << (zeta_before - zeta_indices.size()) / float(zeta_before) * 100 << "%" << std::endl;
		if (zeta_multi_before > 0) std::cout << "Zeta Multi variables eliminated: " << std::fixed << std::setprecision(2) << (zeta_multi_before - zeta_multi_indices.size()) / float(zeta_multi_before) * 100 << "%" << std::endl;
		if (rho_multi_before > 0) std::cout << "Rho Multi variables eliminated: " << std::fixed << std::setprecision(2) << (rho_multi_before - rho_multi_indices.size()) / float(rho_multi_before) * 100 << "%" << std::endl;
		if (sigma_before > 0) std::cout << "Sigma variables eliminated: " << std::fixed << std::setprecision(2) << (sigma_before - sigma_indices.size()) / float(sigma_before) * 100 << "%" << std::endl;
		if (omega_before > 0) std::cout << "Omega variables eliminated: " << std::fixed << std::setprecision(2) << (omega_before - omega_indices.size()) / float(omega_before) * 100 << "%" << std::endl;
    }

    // Construct the variables of the model
    void construct_variables(const char& vtype) {
        // ----------------------------- Main Variables -----------------------------
        // Tau
        for (const int& idx : tau_indices) {
            const std::string name = "tau[" + itos(idx) + "]";
            tau[idx] = model->addVar(0, 1, profit[idx], GRB_CONTINUOUS, name);
        }

        // Gamma
        for (const int& idx : gamma_indices) {
            const int k = int(floor(idx / no_cust_d));
            const int j = idx - k * no_cust_d;
            const std::string name = "gamma[" + itos(k) + "][" + itos(j) + "]";
            gamma[idx] = model->addVar(0.0, GRB_INFINITY, -veh_dist[k][j], vtype, name);
        }

        // Chi
        for (const int& idx : chi_indices) {
            const int k = int(floor(idx / no_cust_d_squared));
            const int i = int(floor((idx - k * no_cust_d_squared) / no_cust_d));
            const int j = idx - k * no_cust_d_squared - i * no_cust_d;
            const std::string name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
            chi[idx] = model->addVar(0.0, 1.0, -cust_dist[i][j], vtype, name);
        }

        // Psi
        for (const int& idx : psi_indices) {
            const int i = int(floor(idx / no_cust_d));
            const int j = idx - i * no_cust_d;
            const std::string name = "psi[" + itos(i) + "][" + itos(j) + "]";
            psi[idx] = model->addVar(0.0, 1.0, -cust_dist[i][j], vtype, name);
        }
        // ----------------------------- Capacity Variables -----------------------------
        // Phi
        for (const int& idx : phi_indices) {
            name = "phi[" + itos(idx) + "]";
            phi[idx] = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        }

        // Zeta
        for (const int& idx : zeta_indices) {
            const int i = int(floor(idx / no_cust_d));
            const int j = idx - i * no_cust_d;
            const std::string name = "zeta[" + itos(i) + "][" + itos(j) + "]";
            zeta[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        }
        
        // Zeta_Multi
        for (const int& idx : zeta_multi_indices) {
            const int l = int(floor(idx / no_cust_d_squared));
            const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
            const int j = idx - l * no_cust_d_squared - i * no_cust_d;
            const std::string name = "zeta_multi[" + itos(l) + "][" + itos(i) + "][" + itos(j) + "]";
            zeta_multi[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        }

        // Rho_Multi
        for (const int& idx : rho_multi_indices) {
            const int l = int(floor(idx / no_cust_d_squared));
            const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
            const int j = idx - l * no_cust_d_squared - i * no_cust_d;
            const std::string name = "rho_multi[" + itos(l) + "][" + itos(i) + "][" + itos(j) + "]";
            rho_multi[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        }
        // ----------------------------- Time Variables -----------------------------
        // Gamma_Capital
        for (const int& idx : Gamma_indices) {
            const std::string name = "Gamma[" + itos(idx) + "]";
            Gamma[idx] = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
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

    // Adds Constraint to the model
    void add_constraints() {
        GRBLinExpr lexpr, lexpr2, lexpr3, lexpr4;

        // Main Constraints
        if (main_type == "3I") {
            for (const int& k : vehicles) {
                for (const int& i : dest_and_cust) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += gamma[idx];
                }
                if (veh_occ[k] == 0) model->addConstr(lexpr <= 1); // 2.2
                else model->addConstr(lexpr == 1); // 2.3
                lexpr.clear();
            }
            // 2.4
            for (const int& k : vehicles) {
                // LHS
                for (const int& i : customers) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += gamma[idx];
                }
                // RHS
                for (const int& i : customers) {
                    const int idx = k * pow(no_cust_d, 2) + i * no_cust_d + 0;
                    if (chi_indices.contains(idx)) lexpr2 += chi[idx];
                }
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
            // 2.5
            for (const int& i : customers) {
                // LHS
                for (const int& k : vehicles) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr += chi[idx];
                    }
                }
                // RHS
                if (tau_indices.contains(i)) lexpr2 += tau[i];
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
            // 2.6
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    // LHS
                    int idx = k * no_cust_d + i;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                    for (const int& j : customers) {
                        if (i == j) continue;
                        idx = k * no_cust_d_squared + j * no_cust_d + i;
                        if (chi_indices.contains(idx)) lexpr += chi[idx];
                    }
                    // RHS
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr2 += chi[idx];
                    }
                    model->addConstr(lexpr == lexpr2);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }

        }
        else if (main_type == "2I") {
            for (const int& k : vehicles) {
                for (const int& i : dest_and_cust) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += gamma[idx];
                }
                if (veh_occ[k] == 0) model->addConstr(lexpr <= 1); // 2.25
                else model->addConstr(lexpr == 1); // 2.26
                lexpr.clear();
            }

            // 2.27
            // LHS
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += gamma[idx];
                }
            }
            // RHS
            for (const int& i : customers) {
                const int idx = i * no_cust_d + 0;
                if (psi_indices.contains(idx)) lexpr2 += psi[idx];
            }
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();

            // 2.28
            for (const int& i : customers) {
                // LHS
                for (const int& j : dest_and_cust) {
                    const int idx = i * no_cust_d + j;
                    if (psi_indices.contains(idx)) lexpr += psi[idx];
                }
                // RHS
				if (tau_indices.contains(i)) lexpr2 += tau[i];
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }

            // 2.29
            for (const int& i : customers) {
                // LHS
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += gamma[idx];
                }

                for (const int& j : customers) {
                    const int idx = j * no_cust_d + i;
                    if (psi_indices.contains(idx)) lexpr += psi[idx];
                }
                // RHS
                if (tau_indices.contains(i)) lexpr2 += tau[i];
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }

        }

        // Capacity Constraints
        if (cap_type == "DL") {
            // DL.1
            for (const int& i : customers) {
                for (const int& j : customers) {
                    if (i == j) continue;
                    // LHS
					if (phi_indices.contains(i)) lexpr += phi[i];
                    if (phi_indices.contains(j)) lexpr -= phi[j];
                    for (const int& k : vehicles) {
						const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr += max_cap * chi[idx];
                    }
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                        if (chi_indices.contains(idx)) lexpr += (max_cap - demand[i]) * chi[idx];
                    }
                    // RHS
                    lexpr2 += max_cap;
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr2 -= demand[j] * chi[idx];
                    }
					model->addConstr(lexpr <= lexpr2);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }
            // DL.2
            for (const int& i : customers) {
                if (!phi_indices.contains(i)) continue;
                // RHS
                for (const int& k : vehicles) {
					int idx = k * no_cust_d + i;
					if (gamma_indices.contains(idx)) lexpr2 += demand[i] * gamma[idx];
                    for (const int& j :customers) {
                        if (i == j) continue;
                        idx = k * no_cust_d_squared + j * no_cust_d + i;
                        if (chi_indices.contains(idx)) lexpr2 += (demand[i] + demand[j]) * chi[idx];
					}
                }
				model->addConstr(phi[i] >= lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }

            // DL.3
            for (const int& i : customers) {
                if (!phi_indices.contains(i)) continue;
                // RHS
                for (const int& k : vehicles) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr2 += (capacity[k] - veh_occ[k] - demand[j]) * chi[idx];
                    }
                }
                model->addConstr(phi[i] <= lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }
        else if (cap_type == "COCF") {
            // COCF.1
            for (const int& i : customers) {
                // LHS
                for (const int& k : vehicles) {
                    int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += (capacity[k] - veh_occ[k]) * gamma[idx];
                }

                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (zeta_indices.contains(idx)) lexpr += zeta[idx];
                    if (main_type == "2I") {
                        const int idx = j * no_cust_d + i;
                        if (psi_indices.contains(idx)) {
                            lexpr += demand[i] * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                            if (chi_indices.contains(idx)) lexpr += demand[i] * chi[idx];
                        }
                    }
                }
                // RHS
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (zeta_indices.find(idx) != zeta_indices.end()) lexpr2 += zeta[idx];
                    if (main_type == "2I") {
                        const int idx = i * no_cust_d + j;
                        if (psi_indices.contains(idx)) lexpr2 += demand[j] * psi[idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) lexpr2 += demand[j] * chi[idx];
                        }
                    }
                }
                if (tau_indices.contains(i)) lexpr2 += demand[i] * tau[i];

                model->addConstr(lexpr >= lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
            // COCF.2
            for (const int& i : customers) {
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int zeta_idx = i * no_cust_d + j;
                    if (!zeta_indices.contains(zeta_idx)) continue;
                    // RHS
                    if (main_type == "2I") {
                        const int idx = i * no_cust_d + j;
                        if (psi_indices.contains(idx)) {
                            lexpr2 += (max_cap - demand[i] - demand[j]) * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) {
                                lexpr2 += (capacity[k] - veh_occ[k] - demand[i] - demand[j]) * chi[idx];
                            }
                        }
                    }
                    model->addConstr(zeta[zeta_idx] <= lexpr2 );
                    lexpr2.clear();
                }
            }
        }
        else if (cap_type == "CMCF") {
            std::vector<std::vector<GRBLinExpr>> zeta_out(no_cust_d);
            std::vector<std::vector<GRBLinExpr>> zeta_in(no_cust_d);
            std::vector<std::vector<GRBLinExpr>> rho_out(no_cust_d);
            std::vector<std::vector<GRBLinExpr>> rho_in(no_cust_d);

            for (const int& l : customers) {
                zeta_out[l].resize(no_cust_d);
                zeta_in[l].resize(no_cust_d);
                rho_out[l].resize(no_cust_d);
                rho_in[l].resize(no_cust_d);
            }
            for (const int& idx : zeta_multi_indices) {
                const int l = int(floor(idx / no_cust_d_squared));
                const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
                const int j = idx - l * no_cust_d_squared - i * no_cust_d;
				zeta_out[l][i] += zeta_multi[idx];
				zeta_in[l][j] += zeta_multi[idx];
            }
            for (const int& idx : rho_multi_indices) {
                const int l = int(floor(idx / no_cust_d_squared));
                const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
                const int j = idx - l * no_cust_d_squared - i * no_cust_d;
                rho_out[l][i] += rho_multi[idx];
                rho_in[l][j] += rho_multi[idx];
            }

            for (const int& l : customers) {
				if (tau_indices.contains(l)) lexpr += tau[l];
                // Zeta
                model->addConstr(zeta_out[l][0] == lexpr); // CMCF.1
                model->addConstr(zeta_in[l][l] == lexpr); // CMCF.2
	            for (const int& i : customers) {
                    if (i == l) continue;
			        model->addConstr(zeta_in[l][i] == zeta_out[l][i]); // CMCF.3
	            }
                // Rho
                model->addConstr(rho_in[l][0] == lexpr); // CMCF.4
                model->addConstr(rho_out[l][l] == lexpr); // CMCF.5
                for (const int& i : customers) {
                    if (i == l) continue;
                    model->addConstr(rho_in[l][i] == rho_out[l][i]); // CMCF.6
                }
                // CMCF.7
	            for (const int& i : customers) {
		            if (i == l) continue;
		            model->addConstr(zeta_in[l][i] == rho_in[i][l]);
	            }
                lexpr.clear();
            }

            // CMCF.8-11
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    if (main_type == "2I") {
                        const int psi_idx = i * no_cust_d + j;
                        if (psi_indices.contains(psi_idx)) lexpr += psi[psi_idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) lexpr += chi[idx];
                        }
                    }
                    for (const int& l : customers) {
                        const int idx = l * no_cust_d_squared + i * no_cust_d + j;
                        if (zeta_multi_indices.contains(idx)) lexpr2 += zeta_multi[idx];
                        if (rho_multi_indices.contains(idx)) lexpr2 += rho_multi[idx];
                        if (lexpr2.size() > 0) model->addConstr(lexpr2 <= lexpr);
                        lexpr2.clear();
                    }
                    lexpr.clear();
                }
            }

            // CMCF.12
            for (const int& i : customers) {
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += gamma[idx];
                }
                for (const int& l : customers) {
                    const int idx = l * no_cust_d_squared + 0 * no_cust_d + i;
                    if (zeta_multi_indices.contains(idx)) lexpr2 += zeta_multi[idx];
                    if (lexpr2.size() > 0) model->addConstr(lexpr2 <= lexpr);
                    lexpr2.clear();
                }
                lexpr.clear();
            }
            // CMCF.13-14
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    if (main_type == "2I") {
                        const int psi_idx = i * no_cust_d + j;
                        if (psi_indices.contains(psi_idx)) lexpr += (max_cap - demand[i] - demand[j]) * psi[psi_idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) lexpr += (capacity[k] - veh_occ[k] - demand[i] - demand[j]) * chi[idx];
                        }
                    }
                    for (const int& l : customers) {
                        if (l == i || l == j) continue;
                        const int idx = l * no_cust_d_squared + i * no_cust_d + j;
                        if (zeta_multi_indices.contains(idx)) lexpr2 += demand[l] * zeta_multi[idx];
                        if (rho_multi_indices.contains(idx)) lexpr2 += demand[l] * rho_multi[idx];
                    }
                    model->addConstr(lexpr2 <= lexpr);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }
            // CMCF.15
            for (const int& i : customers) {
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += (capacity[k] - veh_occ[k] - demand[i]) * gamma[idx];
                }
                for (const int& l : customers) {
                    if (l == i) continue;
                    const int idx = l * no_cust_d_squared + 0 * no_cust_d + i;
                    if (zeta_multi_indices.contains(idx)) lexpr2 += demand[l] * zeta_multi[idx];
                }
                model->addConstr(lexpr2 <= lexpr);
                lexpr.clear();
                lexpr2.clear();
            }
        }

        // Time Constraints
        if (time_type == "TIME") {
            // TIME.1
            for (const int& k : vehicles) {
                // RHS
                for (const int& i : dest_and_cust) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += veh_dist[k][i] * gamma[idx];
                }
                for (const int& i : customers) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr += cust_dist[i][j] * chi[idx];
                    }
                }
				model->addConstr(Gamma[k] == lexpr);
                lexpr.clear();
                lexpr2.clear();
            }
            // TIME.2
            int max_cust_arr_time = 0;
			for (const int& i : customers) max_cust_arr_time = std::max(max_cust_arr_time, cust_arr_time[i]);
            for (const int& k : vehicles) {
                if (max_cust_arr_time > veh_arr_time[k]) {
                    model->addConstr(Gamma[k] <= veh_arr_time[k]);
                }
                else {
                    for (const int& i : customers) {
                        if (cust_arr_time[k] >= veh_arr_time[k]) continue;
                        // RHS
                        lexpr += veh_arr_time[k];
                        for (const int& j : dest_and_cust) {
                            if (i == j) continue;
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) {
                                lexpr += (cust_arr_time[i] - veh_arr_time[k]) * chi[idx];
                            }
                        }
                        model->addConstr(Gamma[k] <= lexpr);
                        lexpr.clear();
                    }
                }
            }
        }
        else if (time_type == "TTCF") {
            // Determine bounds
            std::vector<std::vector<std::vector<int>>> u(no_veh);
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
            // TTCF.1
            for (const int& i : customers) {
                // LHS
                // First term
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (sigma_indices.contains(idx)) lexpr += sigma[idx];
                    if (main_type == "2I") {
                        if (!psi_indices.contains(idx)) continue;
                        int k_min = std::numeric_limits<int>::max();
                        for (const int& k : vehicles) {
                            k_min = std::min(k_min, shortest_path_from_veh[k][j]);
                        }
                        lexpr += k_min * psi[idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                            if (chi_indices.contains(idx)) {
                                lexpr += shortest_path_from_veh[k][j] * chi[idx];
                            }
                        }
                    }
                }
                // Second term
                for (const int& k : vehicles) {
                    int veh_idx = k * no_cust_d + i;
                    if (gamma_indices.contains(veh_idx)) lexpr2 += veh_dist[k][i] * gamma[veh_idx];
                }
                if (main_type == "2I") {
                    for (const int& j : customers) {
                        if (i == j) continue;
                        const int idx = j * no_cust_d + i;
                        if (psi_indices.contains(idx)) lexpr2 += cust_dist[j][i] * psi[idx];
                    }
                }
                else {
                    for (const int& k : vehicles) {
						for (const int& j : customers) {
							if (i == j) continue;
							const int idx = k * no_cust_d_squared + j * no_cust_d + i;
							if (chi_indices.contains(idx)) lexpr2 += cust_dist[j][i] * chi[idx];
						}
                    }
                }
                // RHS
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (sigma_indices.contains(idx)) lexpr3 += sigma[idx];
                    if (main_type == "2I") {
                        if (!psi_indices.contains(idx)) continue;
                        int k_min = std::numeric_limits<int>::max();
                        for (const int& k : vehicles) {
                            k_min = std::min(k_min, shortest_path_from_veh[k][i]);
                        }
                        lexpr3 += k_min * psi[idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) {
                                lexpr3 += shortest_path_from_veh[k][i] * chi[idx];
                            }
                        }
                    }
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
                    if (main_type == "2I") {
                        if (!psi_indices.contains(idx)) continue;
                        int k_min_path = std::numeric_limits<int>::max();
                        int k_max_u = std::numeric_limits<int>::min();
                        for (const int& k : vehicles) {
                            k_min_path = std::min(k_min_path, shortest_path_from_veh[k][i]);
                            k_max_u = std::max(k_max_u, u[k][i][j]);
                        }
                        lexpr += (k_max_u - k_min_path) * psi[idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(idx)) {
                                lexpr += (u[k][i][j] - shortest_path_from_veh[k][i]) * chi[idx];
                            }
                        }
                    }
                    model->addConstr(sigma[idx] <= lexpr);
                    lexpr.clear();
                }
            }
            // TTCF.3 (2.18)
            for (const int& i : customers) {
                // LHS
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr += veh_arr_time[k] * gamma[idx];
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
                    if (main_type == "2I") {
                        if (psi_indices.contains(idx)) lexpr += cust_arr_time[i] * psi[idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int chi_idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.contains(chi_idx)) lexpr += cust_arr_time[i] * chi[chi_idx];
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
                if (main_type == "2I") {
                    if (psi_indices.contains(idx)) lexpr += cust_dist[i][0] * psi[idx];
                }
                else {
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + 0;
                        if (chi_indices.contains(idx))  lexpr += cust_dist[i][0] * chi[idx];
                    }
                }
                // RHS
                if (omega_indices.contains(idx)) lexpr2 += omega[idx];
                if (main_type == "2I") {
                    if (!psi_indices.contains(idx)) continue;
                    int k_min = std::numeric_limits<int>::max();
                    for (const int& k : vehicles) {
                        k_min = std::min(k_min, shortest_path_from_veh[k][i]);
                    }
                    lexpr2 -= k_min * psi[idx];
                }
                else {
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + 0;
                        if (chi_indices.contains(idx)) {
                            lexpr2 -= shortest_path_from_veh[k][i] * chi[idx];
                        }
                    }
                }
                model->addConstr(lexpr <= lexpr2);
                lexpr.clear();
				lexpr2.clear();
            }


        }

		// Standard Capacity constraints
        if (main_type == "3I") {
            // CAP
            for (const int& k : vehicles) {
                // LHS
                for (const int& i : customers) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.contains(idx)) lexpr += demand[i] * chi[idx];
                    }
                }
                // RHS
                for (const int& i : customers) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.contains(idx)) lexpr2 += gamma[idx];
                }
                model->addConstr(lexpr <= (capacity[k] - veh_occ[k]) * lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }
    }

    // Retrieves solution
    void retrieve_solution(std::set<std::tuple<int, int, double>>& x_solution, 
        std::set<std::tuple<int, int, double>>& gamma_solution, 
        std::unordered_map<int, double>& tau_solution
    ){
        if (main_type == "2I") {
            for (const int& idx : psi_indices) {
                const int i = int(floor(idx / no_cust_d));
                const int j = idx - i * no_cust_d;
                double val = psi[idx].get(GRB_DoubleAttr_X);
                if (j == 0) {
                    for (const int& k : vehicles) {
                        const int gamma_idx = k * no_cust_d + i;
                        if (gamma_indices.find(gamma_idx) != gamma_indices.end()) val += gamma[gamma_idx].get(GRB_DoubleAttr_X);
                    }
                }
                if (val >= 1e-3) {
                    if (j == 0) x_solution.insert({ j, i, val });
                    else x_solution.insert({ i, j, val });
                }
            }
        }
        else if (main_type == "3I") {
            const int no_cust_d_squared = int(pow(no_cust_d, 2));
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
					//if (i == j) continue;
                    double val = 0;
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(idx) != chi_indices.end()) val += chi[idx].get(GRB_DoubleAttr_X);
                    }
                    if (j == 0) {
                        for (const int& k : vehicles) {
                            const int gamma_idx = k * no_cust_d + i;
                            if (gamma_indices.find(gamma_idx) != gamma_indices.end()) val += gamma[gamma_idx].get(GRB_DoubleAttr_X);
                        }
                    }
                    if (val >= 1e-3) {
                        if (j == 0) x_solution.insert({ j, i, val });
                        else x_solution.insert({ i, j, val });
                    }
                }
            }
        }

		for (const auto& [idx, var] : gamma) {
			double val = var.get(GRB_DoubleAttr_X);
			if (val >= 1e-3) {
				const int k = int(floor(idx / no_cust_d));
				const int i = idx - k * no_cust_d;
				gamma_solution.insert({ k, i, val });
			}
		}

		for (const auto& [i, var] : tau) {
			double val = var.get(GRB_DoubleAttr_X);
            if (val >= 1e-3) tau_solution[i] = val;
		}
    }

    void add_rc_inequalities(std::vector<RC_Inequality>& violated_inequalties, bool& cuts_added) {
        cuts_added = false;
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        for (RC_Inequality& cut : violated_inequalties) {
            GRBLinExpr lexpr;
            for (const int& i : cut.set) {
                for (const int& j : cut.set) {
                    if (i == j) continue;
                    if (main_type == "2I") {
                        const int idx = i * no_cust_d + j;
                        if (psi_indices.find(idx) != psi_indices.end()) lexpr += psi[idx];

                    }
                    else if (main_type == "3I") {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
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
    int no_cust_d_squared;
    int time_ub;
    int no_veh_vars = 0;
    int no_total_vars = 0;
    int max_cap;
    GRBEnv env = GRBEnv();
    std::string graph_type;
    std::string main_type;
    std::string cap_type;
    std::string time_type;
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
    std::vector<int> profit;
    std::vector<int> cust_arr_time;
    std::vector<std::vector<int>> cust_dist; // Upper triangular matrix of distances - customer zero is the destination
    std::set<int> customers;
    std::set<int> dest_and_cust;

    // Shortest Path
    std::vector<std::vector<int>> shortest_path_from_veh;
    std::vector<int> shortest_path_to_dest;

    // Model Properties
    GRBModel* model = nullptr;
    std::unordered_map<int, GRBVar> tau;
    std::unordered_map<int, GRBVar> gamma;
    std::unordered_map<int, GRBVar> chi;
    std::unordered_map<int, GRBVar> phi;
    std::unordered_map<int, GRBVar> psi;
    std::unordered_map<int, GRBVar> zeta;
    std::unordered_map<int, GRBVar> zeta_multi;
    std::unordered_map<int, GRBVar> rho_multi;
    std::unordered_map<int, GRBVar> sigma;
    std::unordered_map<int, GRBVar> omega;
    std::unordered_map<int, GRBVar> Gamma;


    int no_base_constraints = 0;
    bool preprocessed = false;

    //Valid inequalities
    std::vector<std::set<int>> rounded_cap;
    std::set<int> active_rc_ineqs;
    int total_rcc = 0;

    // Indices
    std::unordered_set<int> tau_indices;
    std::unordered_set<int> gamma_indices;
    std::unordered_set<int> chi_indices;
    std::unordered_set<int> phi_indices;
    std::unordered_set<int> psi_indices;
    std::unordered_set<int> zeta_indices;
    std::unordered_set<int> zeta_multi_indices;
    std::unordered_set<int> rho_multi_indices;
    std::unordered_set<int> Gamma_indices;
    std::unordered_set<int> sigma_indices;
    std::unordered_set<int> omega_indices;

    // Separation
    std::vector<int> vehicle_sizes;
    std::vector<std::vector<int>> vehicles_per_type;

    // Timings
    double total_lp_solve_time = 0;
    double lp_solve_time = 0;
    double lp_rcc_solve_time = 0;
    double mip_solve_time = 0;
    double model_construction_time = 0;
    double graph_construction_time = 0;
    double RC_separation_time = 0;

    // Initializes instance
    FMRSP(const std::string& filename, const std::string main, const std::string cap, const std::string time, const int& threads, const int& time_lim) {
        std::cout << main << "_" << cap << "_" << time << std::endl;

        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        no_veh = inst_data.no_veh;
        no_cust = inst_data.no_cust;
        no_cust_d = no_cust + 1;
        no_cust_d_squared = int(pow(no_cust_d, 2));
        time_ub = inst_data.time_ub;
        num_threads = threads;
        time_limit = time_lim;

        capacity = inst_data.capacity;
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
        max_cap = 0;
        for (const int& k : vehicles) {
            max_cap = std::max(capacity[k] - veh_occ[k], max_cap);
        }

        // Determine profits
		profit.resize(no_cust_d);
        for (const int& i : customers) {
            profit[i] = demand[i] * cust_dist[i][0];
        }

        // Determine Shortest Paths
        shortest_path_from_veh.resize(no_veh);
        for (const int& k : vehicles) {
			shortest_path_from_veh[k].resize(no_cust_d);
            for (const int& i : dest_and_cust) {
				shortest_path_from_veh[k][i] = veh_dist[k][i];
            }
		}
        while (true) {
            bool any_change = false;
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    for (const int& j : customers) {
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

        rounded_cap.resize(no_cust * (no_cust - 1) / 2);

        // Setup Separation Data
        std::set<int> veh_sizes_set;
        for (const int& k : vehicles) {
            if (capacity[k] - veh_occ[k] > 0) veh_sizes_set.insert(capacity[k] - veh_occ[k]);
        }
        vehicle_sizes.resize(veh_sizes_set.size());
        int i = vehicle_sizes.size() - 1;
        for (const int size : veh_sizes_set) {
            vehicle_sizes[i] = size;
            i--;
        }
        vehicles_per_type.resize(veh_sizes_set.size());
        for (const int& k : vehicles) {
            if (capacity[k] - veh_occ[k] == 0) continue;
            int type = 0;
            while (capacity[k] - veh_occ[k] != vehicle_sizes[type]) type++;
            vehicles_per_type[type].emplace_back(k);
        }

        // Set model types
        main_type = main;
        cap_type = cap;
        time_type = time;

        construct_indices();
    }

    // Creates the model
    void create_model(bool relax = true, bool silent = true) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Constructing model..." << std::endl;
        // Create Model Object
        model = new GRBModel(env);

        // Change parameters
        if (silent) model->getEnv().set(GRB_IntParam_OutputFlag, 0);
        model->getEnv().set(GRB_IntParam_Threads, num_threads);

        const char vtype = (relax) ? GRB_CONTINUOUS : GRB_BINARY;
        construct_variables(vtype);
        add_constraints();
        model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

auto stop = std::chrono::high_resolution_clock::now();
model_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Optimises the model
    void optimise_model() {
        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        auto stop = std::chrono::high_resolution_clock::now();
        total_lp_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Creates and solves root node relaxation
    void solve_root_relaxation(bool separate_rci = true, bool silent = true, bool preprocess_instance = true) {
        if (preprocess_instance) preprocess();
        if (model_construction_time == 0) create_model(true, silent);
        model->set(GRB_DoubleParam_BestObjStop, 4761);
        auto start = std::chrono::high_resolution_clock::now();
        optimise_model();
        initial_lp = model->get(GRB_DoubleAttr_ObjVal);
        auto stop = std::chrono::high_resolution_clock::now();
        lp_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();

        if (separate_rci) {
            std::cout << initial_lp << std::endl;
            std::set<std::tuple<int, int, double>> x_solution;
            std::set<std::tuple<int, int, double>> gamma_solution;
            std::unordered_map<int, double> tau_solution;
            std::vector<RC_Inequality> violated_inequalties;

            double max_viol = 0;
            bool cuts_added = false, cuts_added2 = false;
            while (true) {
                retrieve_solution(x_solution, gamma_solution, tau_solution);
                separate_rcc_exactly(violated_inequalties, max_viol, x_solution, customers, max_cap, demand, 1, RC_separation_time, env);
                add_rc_inequalities(violated_inequalties, cuts_added);

                if (!cuts_added) break;

                x_solution.clear();
                gamma_solution.clear();
                tau_solution.clear();
                violated_inequalties.clear();
                optimise_model();
                std::cout << model->get(GRB_DoubleAttr_ObjVal) << std::endl;

            }
        }
        optimise_model();
        strengthened_lp = model->get(GRB_DoubleAttr_ObjVal);
        stop = std::chrono::high_resolution_clock::now();
        lp_rcc_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    void solve_to_optimality(bool separate_rci = true, bool silent = true, bool preprocess_instance = true) {
        if (total_lp_solve_time == 0) {
            solve_root_relaxation(separate_rci, silent, preprocess_instance);
            mip_solve_time += lp_rcc_solve_time;
        }
        model->getEnv().set(GRB_IntParam_OutputFlag, 1);
        model->getEnv().set(GRB_DoubleParam_TimeLimit, time_limit);
        model->getEnv().set(GRB_DoubleParam_MIPGap, 0);
        model->getEnv().set(GRB_DoubleParam_MIPGapAbs, 0.99);

        if (main_type == "2I") {
            for (const int& idx : psi_indices) {
                psi[idx].set(GRB_CharAttr_VType, GRB_BINARY);
            }
            for (const int& idx : gamma_indices) {
                gamma[idx].set(GRB_CharAttr_VType, GRB_BINARY);
            }
        }
        else if (main_type == "3I") {
            for (const int& idx : chi_indices) {
                chi[idx].set(GRB_CharAttr_VType, GRB_BINARY);
            }
        }

        auto start = std::chrono::high_resolution_clock::now();
        model->optimize();
        mip = model->get(GRB_DoubleAttr_ObjVal);
        auto stop = std::chrono::high_resolution_clock::now();
        mip_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
    }

    // Prints routes
    void print_routes() {
        std::unordered_set<int> psi_indices_in_routes;
        std::unordered_set<int> customers_visited;
        int total_profit = 0;
        for (const int k : vehicles) {
            std::cout << "Vehicle " << k << " (" << capacity[k] - veh_occ[k] << "): ";
            int current_vert = 0;
            int next_vert = -1;
            int time = 0;
            int dem = 0;
			int latest_arr_time = veh_arr_time[k];
            int profit = 0;
            for (const int& i : dest_and_cust) {
                int idx = k * no_cust_d + i;
                if (!gamma_indices.contains(idx)) continue;
                if (gamma[idx].get(GRB_DoubleAttr_X) > 0.5) {
                    next_vert = i;
					time += veh_dist[k][i];
                    break;
                }
            }
            if (next_vert == -1) {
                std::cout << std::endl;
                continue;
            }
            while (true) {
                std::cout << next_vert << " ";
                if (next_vert == 0) {
                    profit -= time;
                    total_profit += profit;
                    std::cout << " -- " << "(" << dem << ", " << time << ") -- " << latest_arr_time << " " << profit << std::endl;
                    break;
                }
				current_vert = next_vert;
				dem += demand[current_vert];
				latest_arr_time = std::min(latest_arr_time, cust_arr_time[current_vert]);
				profit += demand[current_vert] * cust_dist[current_vert][0];
		        customers_visited.insert(current_vert);
                if (main_type == "2I") {
                    for (const int& i : dest_and_cust) {
                        int idx = current_vert * no_cust_d + i;
                        if (!psi_indices.contains(idx)) continue;
                        if (psi[idx].get(GRB_DoubleAttr_X) > 0.5) {
                            next_vert = i;
							time += cust_dist[current_vert][i];
							psi_indices_in_routes.insert(idx);
                            break;
                        }
                    }
                }
            }
        }
        for (const int& idx : psi_indices) {
            if (psi[idx].get(GRB_DoubleAttr_X) > 0.1 && !psi_indices_in_routes.contains(idx)) {
				std::cout << "Subtour: " << int(floor(idx / no_cust_d)) << " -> " << idx - int(floor(idx / no_cust_d)) * no_cust_d << std::endl;
            }
        }
        for (const int& idx : tau_indices) {
            if (tau[idx].get(GRB_DoubleAttr_X) > 0.1 && !customers_visited.contains(idx)) {
                std::cout << "Customer " << idx << " not serviced." << std::endl;
            }
		}
        std::cout << total_profit << std::endl;
    }
};