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
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        // Gamma
        if (main_type == "3I" ||  main_type == "2I") {
            for (const int& k : vehicles) {
                for (const int& i : dest_and_cust) {
                    gamma_indices.insert(k * no_cust_d + i);
                }
            }
        }

        // Chi
        if (main_type == "3I" || time_type == "TIME") {
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    for (const int j : dest_and_cust) {
                        if (i == j) continue;
                        chi_indices.insert(k * no_cust_d_squared + i * no_cust_d + j);
                    }
                }
            }
        }

        // Phi
        if (sec_type == "DL") {
            for (const int& i : customers) {
                phi_indices.insert(i);
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

        // Zeta
        if (sec_type == "COCF") {
            for (const int& i : customers) {
                for (const int j : dest_and_cust) {
                    if (i == j) continue;
                    zeta_indices.insert(i * no_cust_d + j);
                }
            }
        }

        // Rho
        if (time_type == "TOCF") {
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    rho_indices.insert(i * no_cust_d + j);
                }
            }
        }

        // Rho_veh
        if (time_type == "TOCF") {
            for (const int& k : vehicles) {
                for (const int& j : customers) {
                    rho_veh_indices.insert(k * no_cust_d + j);
                }
            }
        }

        // Zeta + Rho Multi
        if (sec_type == "CMCF") {
            for (const int& l : customers) {
                for (const int& i : dest_and_cust) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        zeta_multi_indices.insert(idx);
                        rho_multi_indices.insert(idx);
                    }
                }
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

        // Gamma_Capital
        if (time_type == "TIME") {
            for (const int& k : vehicles) {
                Gamma_indices.insert(k);
            }
        }

        // tau
        for (const int& i : customers) {
            tau_indices.insert(i);
        }
    };

    void preprocess_main(std::unordered_set<int>& gamma_indices_removed, std::unordered_set<int>& chi_indices_removed,
        std::unordered_set<int>& phi_indices_removed, std::unordered_set<int>& xi_indices_removed, std::unordered_set<int>& psi_indices_removed, std::unordered_set<int>& zeta_indices_removed,
        std::unordered_set<int>& rho_indices_removed, std::unordered_set<int>& rho_veh_indices_removed, std::unordered_set<int>& zeta_multi_indices_removed, std::unordered_set<int>& rho_multi_indices_removed,
        std::unordered_set<int>& lambda_indices_removed, std::unordered_set<int>& sigma_indices_removed, std::unordered_set<int>& omega_indices_removed
    ) {
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        for (const int& k : vehicles) {
            if (veh_occ[k] == 0) {
                int idx = k * no_cust_d;
                gamma_indices_removed.insert(idx);
            }
        }
    };

    void preprocess_sec(std::unordered_set<int>& gamma_indices_removed, std::unordered_set<int>& chi_indices_removed,
        std::unordered_set<int>& phi_indices_removed, std::unordered_set<int>& xi_indices_removed, std::unordered_set<int>& psi_indices_removed, std::unordered_set<int>& zeta_indices_removed,
        std::unordered_set<int>& rho_indices_removed, std::unordered_set<int>& rho_veh_indices_removed, std::unordered_set<int>& zeta_multi_indices_removed, std::unordered_set<int>& rho_multi_indices_removed,
        std::unordered_set<int>& lambda_indices_removed, std::unordered_set<int>& sigma_indices_removed, std::unordered_set<int>& omega_indices_removed
    ) {
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        if (main_type == "3I") {
            for (const int& i : customers) {
                for (const int& k : vehicles) {
                    if (demand[i] + veh_occ[k] > capacity[k]) {
                        int idx = k * no_cust_d + i;
                        gamma_indices_removed.insert(idx);
                        for (const int& j : dest_and_cust) {
                            if (i == j) continue;
                            idx = k * no_cust_d_squared + i * no_cust_d + j;
                            chi_indices_removed.insert(idx);
                            if (j == 0) continue;
                            idx = k * no_cust_d_squared + j * no_cust_d + i;
                            chi_indices_removed.insert(idx);
                        }
                    }
                }
                for (const int& j : customers) {
                    for (const int& k : vehicles) {
                        if (demand[i] + demand[j] > capacity[k] - veh_occ[k]) {
                            int idx = k * no_cust_d_squared + j * no_cust_d + i;
                            chi_indices_removed.insert(idx);
                            idx = k * no_cust_d_squared + i * no_cust_d + j;
                            chi_indices_removed.insert(idx);
                        }
                    }
                }
            }
        }
        else {
            for (const int& i : customers) {
                for (const int& k : vehicles) {
                    if (demand[i] + veh_occ[k] > capacity[k]) {
                        int idx = k * no_cust_d + i;
                        gamma_indices_removed.insert(idx);
                    }
                }
                for (const int& j : customers) {
                    if (demand[i] + demand[j] > max_cap) {
                        int idx = j * no_cust_d + i;
                        psi_indices_removed.insert(idx);
                        idx = i * no_cust_d + j;
                        psi_indices_removed.insert(idx);
                    }
                }
            }
        }

        if (sec_type == "DL") {
            for (const int& i : customers) {
                if (max_cap - demand[i] <= 0) phi_indices_removed.insert(i);
                else {
                    bool not_servicable = true;
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d + i;
                        if (gamma_indices_removed.find(idx) == gamma_indices_removed.end()) {
                            not_servicable = false;
                            break;
                        }
                    }
                    if (not_servicable) phi_indices_removed.insert(i);
                }
            }
        }
        else if (sec_type == "COCF") {
            for (const int& i : customers) {
                // Remove Zeta Variables
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    if (max_cap - demand[i] - demand[j] <= 0) {
                        const int zeta_idx = i * no_cust_d + j;
                        zeta_indices_removed.insert(zeta_idx);
                    }
                    else {
                        bool no_edge = true;
                        if (main_type == "3I") {
                            for (const int& k : vehicles) {
                                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                                if (chi_indices_removed.find(idx) == chi_indices_removed.end()) {
                                    no_edge = false;
                                    break;
                                }
                            }
                        }
                        else {
                            const int idx = i * no_cust_d + j;
                            no_edge = psi_indices_removed.find(idx) != psi_indices_removed.end();
                        }
                        if (no_edge) {
                            const int zeta_idx = i * no_cust_d + j;
                            zeta_indices_removed.insert(zeta_idx);
                        }
                    }
                }
            }
        }
        else if (sec_type == "CMCF") {
            // Remove Flow Variables
            for (const int& i : dest_and_cust) {
                for (const int& j : dest_and_cust) {
                    // Forced to zero
                    int idx = i * no_cust_d_squared + j * no_cust_d + 0;
                    zeta_multi_indices_removed.insert(idx);
                    idx = i * no_cust_d_squared + 0 * no_cust_d + j;
                    rho_multi_indices_removed.insert(idx);
                    if (i == j) continue;
                    idx = i * no_cust_d_squared + i * no_cust_d + j;
                    zeta_multi_indices_removed.insert(idx);
                    idx = i * no_cust_d_squared + j * no_cust_d + i;
                    rho_multi_indices_removed.insert(idx);
                }
            }

            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    bool no_edge = true;
                    if (main_type == "3I") {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices_removed.find(idx) == chi_indices_removed.end()) {
                                no_edge = false;
                                break;
                            }
                        }
                    }
                    else {
                        const int idx = i * no_cust_d + j;
                        no_edge = psi_indices_removed.find(idx) != psi_indices_removed.end();
                    }
                    if (no_edge) {
                        if (j > 0) {
                            // 8.24
                            for (const int& l : dest_and_cust) {
                                if (l == i || l == j) continue;
                                int idx = j * no_cust_d_squared + i * no_cust_d + l;
                                zeta_multi_indices_removed.insert(idx);
                                idx = j * no_cust_d_squared + l * no_cust_d + i;
                                zeta_multi_indices_removed.insert(idx);
                                idx = i * no_cust_d_squared + j * no_cust_d + l;
                                rho_multi_indices_removed.insert(idx);
                                idx = i * no_cust_d_squared + l * no_cust_d + j;
                                rho_multi_indices_removed.insert(idx);
                            }
                        }
                        // 8.25
                        for (const int& l : customers) {
                            const int idx = l * no_cust_d_squared + i * no_cust_d + j;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                        }
                    }
                }
            }

            // 8.26
            for (const int& j : customers) {
                bool no_edge = true;
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices_removed.find(idx) == gamma_indices_removed.end()) {
                        no_edge = false;
                        break;
                    }
                }
                if (no_edge) {
                    for (const int& l : customers) {
                        const int idx = l * no_cust_d_squared + 0 * no_cust_d + j;
                        zeta_multi_indices_removed.insert(idx);
                        rho_multi_indices_removed.insert(idx);
                    }
                }
            }
            
            // 8.27 - 8.28
            for (const int& i : dest_and_cust) {
                for (int l = i + 1; l < no_cust_d; l++) {
                    for (int j = l + 1; j < no_cust_d; j++) {
                        if (demand[l] + demand[i] + demand[j] > max_cap) {
                            // Any combination can be removed
                            int idx = l * no_cust_d_squared + i * no_cust_d + j;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                            idx = l * no_cust_d_squared + j * no_cust_d + i;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                            idx = j * no_cust_d_squared + i * no_cust_d + l;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                            idx = j * no_cust_d_squared + l * no_cust_d + i;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                            if (i == 0) continue;
                            idx = i * no_cust_d_squared + l * no_cust_d + j;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                            idx = i * no_cust_d_squared + j * no_cust_d + l;
                            zeta_multi_indices_removed.insert(idx);
                            rho_multi_indices_removed.insert(idx);
                        }
                    }
                }

            }
            
        }

    };

    void preprocess_time(std::unordered_set<int>& gamma_indices_removed, std::unordered_set<int>& chi_indices_removed,
        std::unordered_set<int>& phi_indices_removed, std::unordered_set<int>& xi_indices_removed, std::unordered_set<int>& psi_indices_removed, std::unordered_set<int>& zeta_indices_removed,
        std::unordered_set<int>& rho_indices_removed, std::unordered_set<int>& rho_veh_indices_removed, std::unordered_set<int>& zeta_multi_indices_removed, std::unordered_set<int>& rho_multi_indices_removed,
        std::unordered_set<int>& lambda_indices_removed, std::unordered_set<int>& sigma_indices_removed, std::unordered_set<int>& omega_indices_removed
    ) {
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        if (main_type == "3I") {
            for (const int& k : vehicles) {
                for (const int& i : dest_and_cust) {
                    const int min_veh_cust = std::min(veh_arr_time[k], cust_arr_time[i]);
                    if (veh_dist[k][i] + cust_dist[i][0] > min_veh_cust) {
                        int idx = k * no_cust_d + i;
                        gamma_indices_removed.insert(idx);
                        for (const int& j : dest_and_cust) {
                            if (i == j) continue;
                            idx = k * no_cust_d_squared + i * no_cust_d + j;
                            chi_indices_removed.insert(idx);
                            if (j == 0) continue;
                            idx = k * no_cust_d_squared + j * no_cust_d + i;
                            chi_indices_removed.insert(idx);
                        }
                    }
                    if (i == 0) continue;
                    for (const int& j : customers) {
                        if (i == j) continue;
                        if (veh_dist[k][i] + cust_dist[i][j] + cust_dist[j][0] > std::min(min_veh_cust, cust_arr_time[j])) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            chi_indices_removed.insert(idx);
                        }
                    }
                }
            }
        }
        else {
            for (const int& i : dest_and_cust) {
                for (const int& k :vehicles){
                    const int min_veh_cust = std::min(veh_arr_time[k], cust_arr_time[i]);
                    if (veh_dist[k][i] + cust_dist[i][0] > min_veh_cust) {
                        const int idx = k * no_cust_d + i;
                        gamma_indices_removed.insert(idx);
                    }
                }
                if (i == 0) continue;
                for (const int& j : customers) {
                    if (i == j) continue;
                    bool not_servicable = true;
                    for (const int& k : vehicles) {
                        const int min_veh_cust = std::min(std::min(veh_arr_time[k], cust_arr_time[i]), cust_arr_time[j]);
                        if (veh_dist[k][i] + cust_dist[i][j] + cust_dist[j][0] <= min_veh_cust) {
                            not_servicable = false;
                            break;
                        }
                    }
                    if (not_servicable) {
                        const int idx = i * no_cust_d + j;
                        psi_indices_removed.insert(idx);
                    }
                }
            }
        }

        if (time_type == "TOCF") {
            for (const int& k : vehicles) {
                for (const int& i : dest_and_cust) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices_removed.find(idx) != gamma_indices_removed.end()) rho_veh_indices_removed.insert(idx);
                }
            }
            // Remove Rho
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    int min_dist = cust_dist[j][0];
                    for (const int& l : customers) {
                        if (l == i || l == j) continue;
                        min_dist = std::min(min_dist, cust_dist[j][l]);
                    }
                    if (std::min(cust_arr_time[i], cust_arr_time[j]) - cust_dist[i][j] - min_dist <= 0) {
                        const int rho_idx = i * no_cust_d + j;
                        rho_indices_removed.insert(rho_idx);
                    }
                    else {
                        bool no_edge = true;
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices_removed.find(idx) == chi_indices_removed.end()) {
                                no_edge = false;
                                break;
                            }
                        }
                        if (no_edge) {
                            const int rho_idx = i * no_cust_d + j;
                            rho_indices_removed.insert(rho_idx);
                        }
                    }
                }
            }
        }
        else if (time_type == "TTCF") {
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    bool no_edge = true;
                    if (main_type == "3I") {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices_removed.find(idx) == chi_indices_removed.end()) {
                                no_edge = false;
                                break;
                            }
                        }
                    }
                    else {
                        const int idx = i * no_cust_d + j;
                        no_edge = psi_indices_removed.find(idx) != psi_indices_removed.end();
                    }
                    if (no_edge) {
                        const int flow_idx = i * no_cust_d + j;
                        sigma_indices_removed.insert(flow_idx);
                        omega_indices_removed.insert(flow_idx);
                    }
                    else if (std::min(cust_arr_time[i], cust_arr_time[j]) - cust_dist[i][j] - cust_dist[j][0] <= 0) {
                        const int flow_idx = i * no_cust_d + j;
                        sigma_indices_removed.insert(flow_idx);
                    }
                }
            }
        }
    };

    void preprocess_cap(std::unordered_set<int>& gamma_indices_removed, std::unordered_set<int>& chi_indices_removed,
        std::unordered_set<int>& phi_indices_removed, std::unordered_set<int>& xi_indices_removed, std::unordered_set<int>& psi_indices_removed, std::unordered_set<int>& zeta_indices_removed,
        std::unordered_set<int>& rho_indices_removed, std::unordered_set<int>& rho_veh_indices_removed, std::unordered_set<int>& zeta_multi_indices_removed, std::unordered_set<int>& rho_multi_indices_removed,
        std::unordered_set<int>& lambda_indices_removed, std::unordered_set<int>& sigma_indices_removed, std::unordered_set<int>& omega_indices_removed
    ) {
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        for (const int& i : customers) {
            int vehicles_removed_counter = 0;
            for (const int& k : vehicles) {
                if (demand[i] + veh_occ[k] > capacity[k]) {
                    vehicles_removed_counter++;
                    int idx = k * no_cust_d + i;
                    gamma_indices_removed.insert(idx);

                    if (main_type == "3I") {
                        for (const int& j : dest_and_cust) {
                            if (i == j) continue;
                            idx = k * no_cust_d_squared + j * no_cust_d + i;
                            chi_indices_removed.insert(idx);
                            idx = k * no_cust_d_squared + i * no_cust_d + j;
                            chi_indices_removed.insert(idx);
                        }
                    }
                }
            }
            if (main_type == "3I" && no_veh == vehicles_removed_counter) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    int idx = j * no_cust_d + i;
                    psi_indices_removed.insert(idx);
                    idx =  i * no_cust_d + j;
                    psi_indices_removed.insert(idx);
                }
            }
            
        }
    };

    // Construct the variables of the model
    void construct_variables(const char& vtype) {
        const int no_cust_d_squared = int(pow(no_cust_d, 2));
        // Gamma
        for (const int& idx : gamma_indices) {
            const int k = int(floor(idx / no_cust_d));
            const int j = idx - k * no_cust_d;
            const std::string name = "gamma[" + itos(k) + "][" + itos(j) + "]";
            if (main_type == "2I") gamma[idx] = model->addVar(0.0, 1.0, veh_dist[k][j], vtype, name);
            else gamma[idx] = model->addVar(0.0, 1.0, -veh_dist[k][j], vtype, name);
        }

        // Chi
        for (const int& idx : chi_indices) {
            const int k = int(floor(idx / no_cust_d_squared));
            const int i = int(floor((idx - k * no_cust_d_squared) / no_cust_d));
            const int j = idx - k * no_cust_d_squared - i * no_cust_d;
            const std::string name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
            const int obj = -cust_dist[i][j]; // Change to penalty
            if (main_type == "LB") chi[idx]= model->addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
            else chi[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
        }

        // Phi
        for (const int& idx : phi_indices) {
            name = "phi[" + itos(idx) + "]";
            phi[idx] = model->addVar(0, max_cap - demand[idx], 0, GRB_CONTINUOUS, name);
        }

        // Xi
        for (const int& idx : xi_indices) {
            const int k = int(floor(idx / no_cust_d));
            const int j = idx - k * no_cust_d;
            const std::string name = "xi[" + itos(k) + "][" + itos(j) + "]";
            xi[idx] = model->addVar(0.0, 1.0, - demand[j] * cust_dist[j][0], vtype, name); // Change to Penalty
        }

        // Psi
        for (const int& idx : psi_indices) {
            const int i = int(floor(idx / no_cust_d));
            const int j = idx - i * no_cust_d;
            const std::string name = "psi[" + itos(i) + "][" + itos(j) + "]";
            double obj = -cust_dist[i][j];
            psi[idx] = model->addVar(0.0, 1.0, obj, vtype, name);
        }

        // Zeta
        for (const int& idx : zeta_indices) {
            const int i = int(floor(idx / no_cust_d));
            const int j = idx - i * no_cust_d;
            const std::string name = "zeta[" + itos(i) + "][" + itos(j) + "]";
            zeta[idx] = model->addVar(0.0, max_cap - demand[i], 0, GRB_CONTINUOUS, name);
        }

        // Rho
        for (const int& idx : rho_indices) {
            const int i = int(floor(idx / no_cust_d));
            const int j = idx - i * no_cust_d;
            const std::string name = "rho[" + itos(i) + "][" + itos(j) + "]";
            rho[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        }

        // Rho_veh
        for (const int& idx : rho_veh_indices) {
            const int k = int(floor(idx / no_cust_d));
            const int j = idx - k * no_cust_d;
            const std::string name = "rho_veh[" + itos(k) + "][" + itos(j) + "]";
            rho_veh[idx] = model->addVar(0.0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        }
        
        // Zeta_Multi
        for (const int& idx : zeta_multi_indices) {
            const int l = int(floor(idx / no_cust_d_squared));
            const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
            const int j = idx - l * no_cust_d_squared - i * no_cust_d;
            const std::string name = "zeta_multi[" + itos(l) + "][" + itos(i) + "][" + itos(j) + "]";
            zeta_multi[idx] = model->addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
        }

        // Rho_Multi
        for (const int& idx : rho_multi_indices) {
            const int l = int(floor(idx / no_cust_d_squared));
            const int i = int(floor((idx - l * no_cust_d_squared) / no_cust_d));
            const int j = idx - l * no_cust_d_squared - i * no_cust_d;
            const std::string name = "rho_multi[" + itos(l) + "][" + itos(i) + "][" + itos(j) + "]";
            rho_multi[idx] = model->addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
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

        // T
        for (const int& idx : Gamma_indices) {
            const std::string name = "T[" + itos(idx) + "]";
            Gamma[idx] = model->addVar(0, veh_arr_time[idx], 0, GRB_CONTINUOUS, name);
        }

        // Tau
        for (const int& idx : tau_indices) {
            const std::string name = "tau[" + itos(idx) + "]";
            tau[idx] = model->addVar(0, 1, demand[idx] * cust_dist[idx][0], GRB_CONTINUOUS, name);
        }
    }

    // Adds Constraint to the model
    void add_constraints() {
        GRBLinExpr lexpr, lexpr2, lexpr3, lexpr4;
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        if (main_type == "3I") {
            for (const int& k : vehicles) {
                for (const int& j : dest_and_cust) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                }
                if (veh_occ[k] == 0) model->addConstr(lexpr <= 1);
                else model->addConstr(lexpr == 1);
                lexpr.clear();
            }

            for (const int& k : vehicles) {
                for (const int& j : customers) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                }
                for (const int& i : customers) {
                    const int idx = k * pow(no_cust_d, 2) + i * no_cust_d + 0;
                    if (chi_indices.find(idx) != chi_indices.end()) lexpr2 += chi[idx];
                }
                model->addConstr(lexpr == lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }

            for (const int& i : customers) {
                if (tau_indices.find(i) == tau_indices.end()) continue;
                for (const int& k : vehicles) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
                    }
                }
                model->addConstr(lexpr == tau[i]);
                lexpr.clear();
            }

            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    if (tau_indices.find(i) == tau_indices.end()) continue;
                    int idx = k * no_cust_d + i;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];

                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        if (j > 0) {
                            idx = k * no_cust_d_squared + j * no_cust_d + i;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
                        }
                        idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr2 += chi[idx];
                    }
                    model->addConstr(lexpr == lexpr2);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }

        }
        else if (main_type == "2I") {
            // 2.16/2.17
            for (const int& k : vehicles) {
                for (const int& j : dest_and_cust) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                }
                if (veh_occ[k] == 0) model->addConstr(lexpr <= 1);
                else model->addConstr(lexpr == 1);
                lexpr.clear();
            }

            // Same no vehicle starting and ending
            for (const int& k : vehicles) {
                for (const int& j : customers) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                }
            }
            for (const int& i : customers) {
                const int idx = i * no_cust_d + 0;
                if (psi_indices.find(idx) != psi_indices.end()) lexpr2 += psi[idx];
            }
            model->addConstr(lexpr == lexpr2);
            lexpr.clear();
            lexpr2.clear();

            // Max one service - Chceck if necessary after SEC
            for (const int& i : customers) {
                if (tau_indices.find(i) == tau_indices.end()) continue;

                for (const int& j : dest_and_cust) {
                    const int idx = i * no_cust_d + j;
                    if (psi_indices.find(idx) != psi_indices.end()) lexpr += psi[idx];
                }
                model->addConstr(lexpr == tau[i]);
                lexpr.clear();
            }

            // In == Out
            for (const int& i : customers) {
                if (tau_indices.find(i) == tau_indices.end()) continue;

                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                }

                for (const int& j : customers) {
                    const int idx = j * no_cust_d + i;
                    if (psi_indices.find(idx) != psi_indices.end()) lexpr += psi[idx];
                }

                for (const int& j : dest_and_cust) {
                    const int idx = i * no_cust_d + j;
                    if (psi_indices.find(idx) != psi_indices.end()) lexpr2 += psi[idx];
                }
                model->addConstr(lexpr == tau[i]);

                lexpr.clear();
                lexpr2.clear();
            }

        }

        if (sec_type == "DL") {
            for (const int& i : customers) {
                if (phi_indices.find(i) == phi_indices.end()) continue;
                for (const int& k : vehicles) {
                    for (const int& j : dest_and_cust) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr += (capacity[k] - veh_occ[k] - demand[j]) * chi[idx];
                    }
                    const int veh_idx = k * no_cust_d + i;
                    if (gamma_indices.find(veh_idx) != gamma_indices.end()) lexpr -= demand[i] * gamma[veh_idx];
                    for (const int& j : customers) {
                        const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr -= (demand[i] + demand[j]) * chi[idx];
                    }
                }
                model->addConstr(phi[i] <= lexpr);
                lexpr.clear();
            }

            for (const int& i : customers) {
                //if (phi_indices.find(i) == phi_indices.end()) continue;
                if (phi_indices.find(i) != phi_indices.end()) lexpr = phi[i];
                for (const int& k : vehicles) {
                    for (const int& l : customers) {
                        if (l == i) continue;
                        int idx = k * no_cust_d_squared + l * no_cust_d + i;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr += (demand[i] + demand[l]) * chi[idx];
                    }
                    const int veh_idx = k * no_cust_d + i;
                    if (gamma_indices.find(veh_idx) != gamma_indices.end()) lexpr += demand[i] * gamma[veh_idx];

                }
                for (const int& j : customers) {
                    if (i == j) continue;
                    if (phi_indices.find(j) == phi_indices.end() && phi_indices.find(i) == phi_indices.end()) continue;
                    if (phi_indices.find(j) != phi_indices.end()) lexpr2 = phi[j];
                    for (const int& k : vehicles) {
                        int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr3 += (max_cap + demand[j]) * chi[idx];

                        idx = k * no_cust_d_squared + j * no_cust_d + i;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr3 += (max_cap - demand[i]) * chi[idx];
                        for (const int& l : customers) {
                            if (l == j) continue;
                            int idx = k * no_cust_d_squared + l * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr2 += (demand[j] + demand[l]) * chi[idx];
                        }
                        const int veh_idx = k * no_cust_d + j;
                        if (gamma_indices.find(veh_idx) != gamma_indices.end()) lexpr2 += demand[j] * gamma[veh_idx];
                    }
                    model->addConstr(lexpr - lexpr2 + lexpr3  <= max_cap);
                    lexpr2.clear();
                    lexpr3.clear();
                }
                lexpr.clear();

            }
        }
        if (sec_type == "COCF") {
            for (const int& i : customers) {
                for (const int& k : vehicles) {
                    int idx = k * no_cust_d + i;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += (capacity[k] - veh_occ[k]) * gamma[idx];
                }

                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (zeta_indices.find(idx) != zeta_indices.end()) lexpr += zeta[idx];
                    if (main_type == "2I") {
                        const int idx = j * no_cust_d + i;
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr += demand[i] * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                            if (chi_indices.find(idx) != chi_indices.end()) {
                                lexpr += demand[i] * chi[idx];
                            }
                        }
                    }
                }

                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (zeta_indices.find(idx) != zeta_indices.end()) lexpr2 += zeta[idx];
                    if (main_type == "2I") {
                        const int idx = i * no_cust_d + j;
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr2 += demand[j] * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) {
                                lexpr2 += demand[j] * chi[idx];
                            }
                        }
                    }
                }
                if (tau_indices.find(i) != tau_indices.end()) lexpr3 += tau[i];
                //if (main_type == "2I") {
                //    for (const int& j : dest_and_cust) {
                //        const int idx = i * no_cust_d + j;
                //        if (psi_indices.find(idx) != psi_indices.end()) lexpr3 += psi[idx];
                //    }
                //}
                //else {
                //    for (const int& k : vehicles) {
                //        for (const int& j : dest_and_cust) {
                //            if (i == j) continue;
                //            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                //            if (chi_indices.find(idx) != chi_indices.end()) lexpr3 += chi[idx];
                //        }
                //    }
                //}
                model->addConstr(lexpr == lexpr2 + demand[i] * lexpr3);
                lexpr.clear();
                lexpr2.clear();
                lexpr3.clear();
            }
            
            for (const int& i : customers) {
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int zeta_idx = i * no_cust_d + j;
                    if (zeta_indices.find(zeta_idx) == zeta_indices.end()) continue;
                    if (main_type == "2I") {
                        const int idx = i * no_cust_d + j;
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr += psi[idx];
                            lexpr2 += (max_cap - demand[i]) * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) {
                                lexpr += chi[idx];
                                lexpr2 += (capacity[k] - veh_occ[k] - demand[i]) * chi[idx];
                            }
                        }
                    }
                    model->addConstr(zeta[zeta_idx] <= lexpr2 - demand[j] * lexpr);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }

            for (const int& i : customers) {
                const int zeta_idx = i * no_cust_d;
                if (zeta_indices.find(zeta_idx) == zeta_indices.end()) continue;
                if (main_type == "2I") {
                    const int idx = i * no_cust_d + 0;
                    if (psi_indices.find(idx) != psi_indices.end()) {
                        lexpr += (max_cap - demand[i]) * psi[idx];
                    }
                }
                else {
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + 0;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr += (capacity[k] - veh_occ[k] - demand[i]) * chi[idx];
                    }
                }
                model->addConstr(zeta[zeta_idx] <= lexpr);
                lexpr.clear();
            }
        }
        if (sec_type == "CMCF") {
            std::vector<std::vector<GRBLinExpr>> zeta_out(no_cust_d);
            std::vector<std::vector<GRBLinExpr>> zeta_in(no_cust_d);
            std::vector<std::vector<GRBLinExpr>> rho_out(no_cust_d);
            std::vector<std::vector<GRBLinExpr>> rho_in(no_cust_d);

            for (const int& l : customers) {
                zeta_out[l].resize(no_cust_d);
                zeta_in[l].resize(no_cust_d);
                rho_out[l].resize(no_cust_d);
                rho_in[l].resize(no_cust_d);
                for (const int& i : dest_and_cust) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = l * pow(no_cust_d, 2) + i * no_cust_d + j;
                        if (zeta_multi_indices.find(idx) != zeta_multi_indices.end()) {
                            zeta_out[l][i] += zeta_multi[idx];
                            zeta_in[l][j] += zeta_multi[idx];
                        }
                        if (rho_multi_indices.find(idx) != rho_multi_indices.end()) {
                            rho_out[l][i] += rho_multi[idx];
                            rho_in[l][j] += rho_multi[idx];
                        }
                    }
                }
            }
            for (const int& l : customers) {
                if (main_type == "2I") {
                    for (int j : dest_and_cust) {
                        if (l == j) continue;
                        const int idx = l * no_cust_d + j;
                        if (psi_indices.find(idx) != psi_indices.end()) lexpr += psi[idx];
                    }
                }
                else {
                    for (const int& k : vehicles) {
                        for (int j : dest_and_cust) {
                            if (l == j) continue;
                            const int idx = k * no_cust_d_squared + l * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
                        }
                    }
                }
                // Zeta
                model->addConstr(zeta_out[l][0] == lexpr);
				for (const int& i : customers) {
                    if (i == l) {
                        model->addConstr(zeta_in[l][i] == zeta_out[l][i] + lexpr);
                    }
                    else {
						model->addConstr(zeta_in[l][i] == zeta_out[l][i]);
                    }
				}
                // Rho
                model->addConstr(rho_in[l][0] == lexpr);
                for (const int& i : customers) {
                    if (i == l) {
                        model->addConstr(rho_in[l][i] == rho_out[l][i] - lexpr);
                    }
                    else {
                        model->addConstr(rho_in[l][i] == rho_out[l][i]);
                    }
                }
				for (const int& i : customers) {
					if (i == l) continue;
					model->addConstr(zeta_in[l][i] == rho_in[i][l]);
				}
                lexpr.clear();
            }
            //for (const int& i : customers) {
            //    if (main_type == "2I") {
            //        for (int j : dest_and_cust) {
            //            if (i == j) continue;
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) lexpr += psi[idx];
            //        }
            //    }
            //    else {
            //        for (const int& k : vehicles) {
            //            for (int j : dest_and_cust) {
            //                if (i == j) continue;
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
            //            }
            //        }
            //    }
            //    model->addConstr(zeta_out[i][0] == lexpr);
            //    model->addConstr(zeta_in[i][i] == lexpr);
            //    model->addConstr(rho_out[i][i] == lexpr);
            //    model->addConstr(rho_in[i][0] == lexpr);
            //    lexpr.clear();
            //    if (!preprocessed) {
            //        model->addConstr(zeta_in[i][0] == 0);
            //        model->addConstr(zeta_out[i][i] == 0);
            //        model->addConstr(rho_in[i][i] == 0);
            //        model->addConstr(rho_out[i][0] == 0);
            //    }
            //    for (const int& j : customers) {
            //        if (i == j) continue;
            //        model->addConstr(zeta_out[i][j] == zeta_in[i][j]);
            //        model->addConstr(zeta_in[i][j] == rho_out[j][i]);
            //        model->addConstr(rho_out[j][i] == rho_in[j][i]);
            //    }
            //}

            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    if (main_type == "2I") {
                        const int psi_idx = i * no_cust_d + j;
                        if (psi_indices.find(psi_idx) != psi_indices.end()) lexpr += psi[psi_idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
                        }
                    }
                    for (const int& l : customers) {
                        const int idx = l * no_cust_d_squared + i * no_cust_d + j;
                        if (zeta_multi_indices.find(idx) != zeta_multi_indices.end()) lexpr2 += zeta_multi[idx];
                        if (rho_multi_indices.find(idx) != rho_multi_indices.end()) lexpr2 += rho_multi[idx];
                        if (lexpr2.size() > 0) model->addConstr(lexpr2 <= lexpr);
                        lexpr2.clear();
                    }
                    lexpr.clear();
                }
            }

            for (const int& j : dest_and_cust) {
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
                }
                for (const int& l : customers) {
                    const int idx = l * no_cust_d_squared + 0 * no_cust_d + j;
                    if (zeta_multi_indices.find(idx) != zeta_multi_indices.end()) lexpr2 += zeta_multi[idx];
                    if (rho_multi_indices.find(idx) != rho_multi_indices.end()) lexpr2 += rho_multi[idx];
                    if (lexpr2.size() > 0) model->addConstr(lexpr2 <= lexpr);
                    lexpr2.clear();
                }
                lexpr.clear();
            }

            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    if (main_type == "2I") {
                        const int psi_idx = i * no_cust_d + j;
                        if (psi_indices.find(psi_idx) != psi_indices.end()) lexpr += (max_cap - demand[i] - demand[j]) * psi[psi_idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr += (capacity[k] - veh_occ[k] - demand[i] - demand[j]) * chi[idx];
                        }
                    }
                    for (const int& l : customers) {
                        if (l == i || l == j) continue;
                        const int idx = l * no_cust_d_squared + i * no_cust_d + j;
                        if (zeta_multi_indices.find(idx) != zeta_multi_indices.end()) lexpr2 += demand[l] * zeta_multi[idx];
                        if (rho_multi_indices.find(idx) != rho_multi_indices.end()) lexpr2 += demand[l] * rho_multi[idx];
                    }
                    model->addConstr(lexpr2 <= lexpr);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }

            for (const int& j : dest_and_cust) {
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += (capacity[k] - veh_occ[k] - demand[j]) * gamma[idx];
                }
                for (const int& l : customers) {
                    if (l == j) continue;
                    const int idx = l * no_cust_d_squared + 0 * no_cust_d + j;
                    if (zeta_multi_indices.find(idx) != zeta_multi_indices.end()) lexpr2 += demand[l] * zeta_multi[idx];
                    if (rho_multi_indices.find(idx) != rho_multi_indices.end()) lexpr2 += demand[l] * rho_multi[idx];
                }
                model->addConstr(lexpr2 <= lexpr);
                lexpr.clear();
                lexpr2.clear();
            }

        }

        if (time_type == "TIME") {
            for (const int& k : vehicles) {
                for (const int& j : dest_and_cust) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += veh_dist[k][j] * gamma[idx];
                }

                for (const int& i : customers) {
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr2 += cust_dist[i][j] * chi[idx];
                    }
                }
				model->addConstr(Gamma[k] == lexpr + lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }

            for (const int& k : vehicles) {
                int no_impl_for_k = 0;
                for (const int& i : customers) {
                    if (cust_arr_time[i] < veh_arr_time[k]) {
						for (const int& j : dest_and_cust) {
							if (i == j) continue;
							const int idx = k * no_cust_d_squared + i * no_cust_d + j;
							if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
						}
                        model->addConstr(Gamma[k] <= veh_arr_time[k] + (cust_arr_time[i] - veh_arr_time[k]) * lexpr);
                        lexpr.clear();
						no_impl_for_k++;
                    }
                }
            }
            //for (const int& k : vehicles) {
            //    for (const int& j : dest_and_cust) {
            //        const int idx = k * no_cust_d + j;
            //        if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += (veh_dist[k][j] - veh_arr_time[k]) * gamma[idx];
            //    }

            //    for (const int& l : customers) {
            //        for (const int& j : dest_and_cust) {
            //            if (l == j) continue;
            //            const int idx = k * no_cust_d_squared + l * no_cust_d + j;
            //            if (chi_indices.find(idx) != chi_indices.end()) lexpr2 += cust_dist[l][j] * chi[idx];
            //        }
            //    }
            //    
            //    int no_impl_for_k = 0;
            //    for (const int& i : customers) {
            //        if (no_impl_for_k > 0 && veh_arr_time[k] <= cust_arr_time[i]) continue;

            //        for (const int& j : dest_and_cust) {
            //            if (i == j) continue;
            //            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //            if (chi_indices.find(idx) != chi_indices.end()) lexpr3 += chi[idx];
            //        }
            //        model->addConstr(lexpr + lexpr2 + (veh_arr_time[k] - cust_arr_time[i]) * lexpr3 <= 0);
            //        no_impl_for_k++;
            //        lexpr3.clear();
            //    }
            //    lexpr.clear();
            //    lexpr2.clear();
            //}
        }
        else if (time_type == "TOCF") {
            // Determine bounds
            std::vector<std::vector<std::vector<int>>> L(no_veh);
            std::vector<std::vector<std::vector<int>>> U(no_veh);
            for (const int& k : vehicles) {
                L[k].resize(no_cust_d);
                U[k].resize(no_cust_d);
                for (const int& i : customers) {
                    L[k][i].resize(no_cust_d);
                    U[k][i].resize(no_cust_d);
                    for (const int& j : dest_and_cust) {
                        if (i == j) continue;
                        int l_min = std::numeric_limits<int>::max();
                        int u_min = std::numeric_limits<int>::max();
                        for (const int& l : customers) {
                            if (l == i || l == j) continue;
                            l_min = std::min(l_min, cust_dist[j][l]);// +veh_dist[k][l]);
                            u_min = std::min(u_min, cust_dist[l][i]);// +cust_dist[l][0]);
                        }
                        if (j == 0) {
                            L[k][i][j] = cust_dist[i][j];
                        }
                        else {
                            L[k][i][j] = cust_dist[i][j] + std::min(cust_dist[j][0], l_min);
                        }
                        U[k][i][j] = std::min(std::min(cust_arr_time[i], cust_arr_time[j]), veh_arr_time[k])
                            - std::min(veh_dist[k][i], u_min);
                    }
                }
            }
            // TOCF2.1
            for (const int& i : customers) {
                // First term
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (rho_veh_indices.find(idx) != rho_veh_indices.end()) lexpr += rho_veh[idx];
                }
                //Second term
                for (const int& j : customers) {
					if (i == j) continue;
					const int idx = j * no_cust_d + i;
					if (rho_indices.find(idx) != rho_indices.end()) lexpr2 += rho[idx];
					for (const int& k : vehicles) {
						const int chi_idx = k * no_cust_d_squared + j * no_cust_d + i;
						if (chi_indices.find(chi_idx) != chi_indices.end()) lexpr2 += L[k][j][i] * chi[chi_idx];
					}
                }
				// Third term
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (rho_indices.find(idx) != rho_indices.end()) lexpr3 += rho[idx];
                    for (const int& k : vehicles) {
                        const int chi_idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(chi_idx) != chi_indices.end()) lexpr3 += L[k][i][j] * chi[chi_idx];
                    }
                }
                // Fourth term
                for (const int& k : vehicles) {
                    for (const int& j : customers) {
                        if (i == j) continue;
						const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                        if (chi_indices.find(idx) != chi_indices.end()) lexpr4 += cust_dist[j][i] * chi[idx];
                    }
                    const int veh_idx = k * no_cust_d + i;
                    if (gamma_indices.find(veh_idx) != gamma_indices.end()) {
                        lexpr4 += veh_dist[k][i] * gamma[veh_idx];
                    }
                }
                model->addConstr(lexpr + lexpr2 == lexpr3 + lexpr4);
                lexpr.clear();
                lexpr2.clear();
                lexpr3.clear();
                lexpr4.clear();
            }

			// TOCF2.2
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
					const int idx = i * no_cust_d + j;
                    if (rho_indices.find(idx) == rho_indices.end()) continue;
					for (const int& k : vehicles) {
						const int chi_idx = k * no_cust_d_squared + i * no_cust_d + j;
                        if (chi_indices.find(chi_idx) != chi_indices.end()) {
                            lexpr += (U[k][i][j] - L[k][i][j]) * chi[chi_idx];
                        }
					}
                    model->addConstr(rho[idx] <= lexpr);
                    lexpr.clear();
                }
            }
            // TOCF.4
            for (const int& k : vehicles) {
                for (const int& i : customers) {
					const int idx = k * no_cust_d + i;
					if (rho_veh_indices.find(idx) == rho_veh_indices.end()) continue;
                    if (gamma_indices.find(idx) != gamma_indices.end()) {
                        lexpr += std::min(veh_arr_time[k], cust_arr_time[k]) * gamma[idx];
                    }
					model->addConstr(rho_veh[idx] <= lexpr);
                    lexpr.clear();
                }
            }
            // TOCF.5
            //for (const int& k : vehicles) {
            //    for (const int& i : customers) {
            //        const int idx = k * no_cust_d + i;
            //        if (rho_veh_indices.find(idx) != rho_veh_indices.end()) lexpr2 += rho_veh[idx];
            //        for (const int& j : customers) {
            //            if (i == j) continue;
            //            for (const int& l : dest_and_cust) {
            //                if (j == l) continue;
            //                const int idx = k * no_cust_d_squared + j * no_cust_d + l;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    lexpr += chi[idx];
            //                }
            //            }
            //            GRBLinExpr RHS = cust_arr_time[j] * lexpr;
            //            +std::min(veh_arr_time[k], cust_arr_time[i]) * (1 - lexpr);
            //            model->addConstr(lexpr2 <= RHS);
            //            lexpr.clear();
            //            lexpr2.clear();
            //        }
            //    }
            //}
            // TOCF.5 - sum i
            int max_cust_arr_time = 0;
            for (const int& i : customers) max_cust_arr_time = std::max(max_cust_arr_time, cust_arr_time[i]);
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    const int idx = k * no_cust_d + i;
                    if (rho_veh_indices.find(idx) != rho_veh_indices.end()) lexpr2 += rho_veh[idx];
                }
                for (const int& j : customers) {
                    for (const int& l : dest_and_cust) {
						if (j == l) continue;
						const int idx = k * no_cust_d_squared + j * no_cust_d + l;
						if (chi_indices.find(idx) != chi_indices.end()) {
							lexpr += chi[idx];
						}
                    }
                    GRBLinExpr RHS = cust_arr_time[j] * lexpr
                        + std::min(veh_arr_time[k], max_cust_arr_time) * (1 - lexpr);
                    model->addConstr(lexpr2 <= RHS);
                    lexpr.clear();
                    lexpr2.clear();
                }
            }
            //// TOCF.5 - Sum k
            //int max_veh_arr_time = 0;
            //for (const int& k : vehicles) max_veh_arr_time = std::max(max_veh_arr_time, veh_arr_time[k]);
            //for (const int& i : customers) {
            //    for (const int& k : vehicles) {
            //        const int idx = k * no_cust_d + i;
            //        if (rho_veh_indices.find(idx) != rho_veh_indices.end()) lexpr2 += rho_veh[idx];
            //    }
            //    for (const int& j : customers) {
            //        if (i == j) continue;
            //        if (tau_indices.find(j) != tau_indices.end()) {
            //            lexpr += tau[j];
            //        }
            //        GRBLinExpr RHS = cust_arr_time[j] * lexpr
            //            + std::min(max_veh_arr_time, cust_arr_time[i]) * (1 - lexpr);
            //        model->addConstr(lexpr2 <= RHS);
            //        lexpr.clear();
            //        lexpr2.clear();
            //    }
            //}
            //for (const int& k : vehicles) {
            //    for (const int& i : customers) {
            //        const int idx = k * no_cust_d + i;
            //        if (rho_veh_indices.find(idx) == rho_veh_indices.end()) continue;
            //        if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += gamma[idx];
            //        
            //        const int min_time = std::min(veh_arr_time[k], cust_arr_time[i]);
            //        // TOCF.4
            //        model->addConstr(rho_veh[idx] <= min_time * lexpr);
            //        lexpr.clear();
            //        for (const int& j : customers) {
            //            if (i == j) continue;
            //            if (cust_arr_time[j] >= min_time) continue;
            //            for (const int& l : dest_and_cust) {
            //                if (l == j) continue;
            //                const int idx = k * no_cust_d_squared + j * no_cust_d + l;
            //                if (chi_indices.find(idx) != chi_indices.end()) lexpr += chi[idx];
            //            }
            //            // TOCF.5
            //            model->addConstr(rho_veh[idx] <= cust_arr_time[j] * lexpr + min_time * (1 - lexpr));
            //            lexpr.clear();
            //        }
            //    }
            //}
            //for (const int& i : customers) {
            //    // First term
            //    for (const int& k : vehicles) {
            //        const int idx = k * no_cust_d + i;
            //        if (rho_veh_indices.find(idx) != rho_veh_indices.end()) lexpr += rho_veh[idx];
            //        if (gamma_indices.find(idx) != gamma_indices.end()) lexpr -= veh_dist[k][i] * gamma[idx];
            //    }
            //    // Second Term
            //    for (const int& j : customers) {
            //        if (i == j) continue;
            //        const int idx = j * no_cust_d + i;
            //        if (rho_indices.find(idx) != rho_indices.end()) lexpr2 += rho[idx];
            //        int min_dist = cust_dist[i][0];
            //        for (const int& l : customers) {
            //            if (l == i || l == j) continue;
            //            min_dist = std::min(min_dist, cust_dist[i][l]);
            //        }
            //        if (main_type == "2I") {
            //            const int idx = j * no_cust_d + i;
            //            if (psi_indices.find(idx) != psi_indices.end()) {
            //                lexpr2 += min_dist * psi[idx];
            //            }
            //        }
            //        else {
            //            for (const int& k : vehicles) {
            //                const int idx = k * no_cust_d_squared + j * no_cust_d + i;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    lexpr2 += min_dist * chi[idx];
            //                }
            //            }
            //        }
            //    }
            //    // Third Term
            //    for (const int& j : dest_and_cust) {
            //        if (i == j) continue;
            //        const int idx = i * no_cust_d + j;
            //        if (rho_indices.find(idx) != rho_indices.end()) lexpr3 += rho[idx];
            //        int min_dist = cust_dist[j][0];
            //        for (const int& l : customers) {
            //            if (l == i || l == j) continue;
            //            min_dist = std::min(min_dist, cust_dist[j][l]);
            //        }
            //        if (main_type == "2I") {
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) {
            //                lexpr3 += min_dist * psi[idx];
            //            }
            //        }
            //        else {
            //            for (const int& k : vehicles) {
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    lexpr3 += min_dist * chi[idx];
            //                }
            //            }
            //        }
            //    }
            //    //Fourth Term
            //    if (main_type == "2I") {
            //        for (const int& j : dest_and_cust) {
            //            if (i == j) continue;
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) lexpr4 += cust_dist[i][j] * psi[idx];
            //        }
            //    }
            //    else {
            //        for (const int& k : vehicles) {
            //            for (const int& j : dest_and_cust) {
            //                if (i == j) continue;
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) lexpr4 += cust_dist[i][j] * chi[idx];
            //            }
            //        }
            //    }
            //    // TOCF2.1
            //    model->addConstr(lexpr + lexpr2 == lexpr3 + lexpr4);
            //    lexpr.clear();
            //    lexpr2.clear();
            //    lexpr3.clear();
            //    lexpr4.clear();
            //}

            //for (const int& i : customers) {
            //    for (const int& j : customers) {
            //        if (i == j) continue;
            //        int min_dist = cust_dist[j][0];
            //        for (const int& l : customers) {
            //            if (l == i || l == j) continue;
            //            min_dist = std::min(min_dist, cust_dist[j][l]);
            //        }
            //        const int rho_idx = i * no_cust_d + j;
            //        if (rho_indices.find(rho_idx) == rho_indices.end()) continue;
            //        const int min_time = std::min(cust_arr_time[i], cust_arr_time[j]);
            //        if (main_type == "2I") {
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) {
            //                //lexpr += psi[idx];
            //                lexpr2 += (min_time - cust_dist[i][j] - min_dist) * psi[idx];
            //            }
            //        }
            //        else {
            //            for (const int& k : vehicles) {
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    //lexpr += chi[idx];
            //                    lexpr2 += (min_time - cust_dist[i][j] - min_dist) * chi[idx];
            //                }
            //            }
            //        }
            //        // TOCF2.2
            //        //model->addConstr(min_dist * lexpr <= rho[rho_idx]);
            //        model->addConstr(rho[rho_idx] <= lexpr2);
            //        lexpr.clear();
            //        lexpr2.clear();
            //    }
            //}
        }
        else if (time_type == "TTCF") {
            // Determine bounds
            std::vector<std::vector<std::vector<int>>> L(no_veh);
            std::vector<std::vector<std::vector<int>>> U(no_veh);
			for (const int& k : vehicles) {
				L[k].resize(no_cust_d);
                U[k].resize(no_cust_d);
				for (const int& i : customers) {
					L[k][i].resize(no_cust_d);
                    U[k][i].resize(no_cust_d);
					for (const int& j : dest_and_cust) {
						if (i == j) continue;
						int l_min = std::numeric_limits<int>::max();
                        int u_min = std::numeric_limits<int>::max();
                        for (const int& l : customers) {
						    if (l == i || l == j) continue;
                            l_min = std::min(l_min, cust_dist[l][i]);// +veh_dist[k][l]);
                            u_min = std::min(u_min, cust_dist[j][l]);// +cust_dist[l][0]);
					    }
                        L[k][i][j] = std::min(veh_dist[k][i], l_min);
                        U[k][i][j] = std::min(std::min(cust_arr_time[i], cust_arr_time[j]), veh_arr_time[k])
                            - std::min(cust_dist[j][0], u_min) - cust_dist[i][j];
					}
				}
			}
            // Time flow - TTCF2.1
            for (const int& i : customers) {
                // First term
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (sigma_indices.find(idx) != sigma_indices.end()) lexpr += sigma[idx];
                    if (main_type == "2I") {
                        const int idx = j * no_cust_d + i;
						int k_min_l = std::numeric_limits<int>::max();
						for (const int& k : vehicles) {
							k_min_l = std::min(k_min_l, L[k][j][i]);
						}
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr += k_min_l * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + j * no_cust_d + i;
                            if (chi_indices.find(idx) != chi_indices.end()) {
                                lexpr += L[k][j][i] * chi[idx];
                            }
                        }
                    }
                }
                // Second term
                if (main_type == "2I") {
                    for (const int& k : vehicles) {
                        int veh_idx = k * no_cust_d + i;
                        if (gamma_indices.find(veh_idx) != gamma_indices.end()) {
                            lexpr2 += veh_dist[k][i] * gamma[veh_idx];
                        }
                    }
                    for (const int& j : customers) {
                        if (i == j) continue;
                        const int idx = j * no_cust_d + i;
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr2 += cust_dist[j][i] * psi[idx];
                        }
                    }
                }
                else {
                    for (const int& k : vehicles) {
						int veh_idx = k * no_cust_d + i;
						if (gamma_indices.find(veh_idx) != gamma_indices.end()) {
							lexpr2 += veh_dist[k][i] * gamma[veh_idx];
						}
						for (const int& j : customers) {
							if (i == j) continue;
							const int idx = k * no_cust_d_squared + j * no_cust_d + i;
							if (chi_indices.find(idx) != chi_indices.end()) {
								lexpr2 += cust_dist[j][i] * chi[idx];
							}
						}
                    }
                }
                // Third term
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (sigma_indices.find(idx) != sigma_indices.end()) lexpr3 += sigma[idx];
                    if (main_type == "2I") {
                        const int idx = i * no_cust_d + j;
                        int k_min_l = std::numeric_limits<int>::max();
                        for (const int& k : vehicles) {
                            k_min_l = std::min(k_min_l, L[k][i][j]);
                        }
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr3 += k_min_l * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) {
                                lexpr3 += L[k][i][j] * chi[idx];
                            }
                        }
                    }
                }
                model->addConstr(lexpr + lexpr2 == lexpr3);
                lexpr.clear();
                lexpr2.clear();
                lexpr3.clear();
            }
            // Upper Bound - TTCF2.2
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (sigma_indices.find(idx) == sigma_indices.end()) continue;
                    if (main_type == "2I") {
                        const int idx = j * no_cust_d + i;
                        int k_min_l = std::numeric_limits<int>::max();
						int k_max_u = std::numeric_limits<int>::min();
                        for (const int& k : vehicles) {
                            k_min_l = std::min(k_min_l, L[k][i][j]);
                            k_max_u = std::max(k_max_u, U[k][i][j]);
                        }
                        if (psi_indices.find(idx) != psi_indices.end()) {
                            lexpr += (k_max_u - k_min_l) * psi[idx];
                        }
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) {
                                lexpr += (U[k][i][j] - L[k][i][j]) * chi[idx];
                            }
                        }
                    }
                    model->addConstr(sigma[idx] <= lexpr);
                    lexpr.clear();
                }
            }
            
            // TTCF2.3
            for (const int& i : customers) {
                // LHS
				const int idx = i * no_cust_d + 0;
                if (sigma_indices.find(idx) != sigma_indices.end()) lexpr += sigma[idx];
                if (main_type == "2I") {
                    if (psi_indices.find(idx) != psi_indices.end()) lexpr += cust_dist[i][0] * psi[idx];
                }
                else {
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + 0;
                        if (chi_indices.find(idx) != chi_indices.end()) {
                            lexpr += cust_dist[i][0] * chi[idx];
                        }
                    }
                }
                // RHS
                if (omega_indices.find(idx) != omega_indices.end()) lexpr2 += omega[idx];
                if (main_type == "2I") {
                    int k_min_l = std::numeric_limits<int>::max();
                    for (const int& k : vehicles) {
                        k_min_l = std::min(k_min_l, L[k][i][0]);
                    }
                    if (psi_indices.find(idx) != psi_indices.end()) lexpr2 -= k_min_l * psi[idx];
                }
                else {
                    for (const int& k : vehicles) {
                        const int idx = k * no_cust_d_squared + i * no_cust_d + 0;
                        if (chi_indices.find(idx) != chi_indices.end()) {
                            lexpr2 -= L[k][i][0] * chi[idx];
                        }
                    }
                }
                model->addConstr(lexpr <= lexpr2);
                lexpr.clear();
				lexpr2.clear();
            }
            //for (const int& i : customers) {
            //    for (const int& k : vehicles) {
            //        const int idx = k * no_cust_d + i;
            //        if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += veh_dist[k][i] * gamma[idx];
            //    }

            //    for (const int& j : customers) {
            //        if (i == j) continue;
            //        const int idx = j * no_cust_d + i;
            //        if (sigma_indices.find(idx) != sigma_indices.end()) lexpr2 += sigma[idx];
            //        if (main_type == "2I") {
            //            const int idx = j * no_cust_d + i;
            //            if (psi_indices.find(idx) != psi_indices.end()) {
            //                lexpr2 += cust_dist[j][i] * psi[idx];
            //            }
            //        }
            //        else {
            //            for (const int& k : vehicles) {
            //                const int idx = k * no_cust_d_squared + j * no_cust_d + i;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    lexpr2 += cust_dist[j][i] * chi[idx];
            //                }
            //            }
            //        }
            //    }

            //    for (const int& j : dest_and_cust) {
            //        if (i == j) continue;
            //        const int idx = i * no_cust_d + j;
            //        if (sigma_indices.find(idx) != sigma_indices.end()) lexpr3 += sigma[idx];
            //        if (main_type == "2I") {
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) {
            //                lexpr3 += cust_dist[i][j] * psi[idx];
            //            }
            //        }
            //        else {
            //            for (const int& k : vehicles) {
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    lexpr3 += cust_dist[i][j] * chi[idx];
            //                }
            //            }
            //        }
            //    }

            //    if (main_type == "2I") {
            //        for (const int& j : dest_and_cust) {
            //            if (i == j) continue;
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) lexpr4 += cust_dist[i][j] * psi[idx];
            //        }
            //    }
            //    else {
            //        for (const int& k : vehicles) {
            //            for (const int& j : dest_and_cust) {
            //                if (i == j) continue;
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) lexpr4 += cust_dist[i][j] * chi[idx];
            //            }
            //        }
            //    }
            //    model->addConstr(lexpr + lexpr2 == lexpr3 - lexpr4);
            //    lexpr.clear();
            //    lexpr2.clear();
            //    lexpr3.clear();
            //    lexpr4.clear();
            //}

            //for (const int& i : customers) {
            //    for (const int& j : dest_and_cust) {
            //        if (i == j) continue;
            //        const int sigma_idx = i * no_cust_d + j;
            //        if (sigma_indices.find(sigma_idx) == sigma_indices.end()) continue;
            //        const int min_time = std::min(cust_arr_time[i], cust_arr_time[j]);
            //        if (main_type == "2I") {
            //            const int idx = i * no_cust_d + j;
            //            if (psi_indices.find(idx) != psi_indices.end()) {
            //                lexpr += psi[idx];
            //            }
            //        }
            //        else {
            //            for (const int& k : vehicles) {
            //                const int idx = k * no_cust_d_squared + i * no_cust_d + j;
            //                if (chi_indices.find(idx) != chi_indices.end()) {
            //                    lexpr += chi[idx];
            //                }
            //            }
            //        }
            //        //model->addConstr(cust_dist[i][j] * lexpr <= sigma[sigma_idx]);
            //        if (j > 0) model->addConstr(sigma[sigma_idx] <= (min_time - cust_dist[i][j] - cust_dist[j][0]) * lexpr);
            //        else model->addConstr(sigma[sigma_idx] <= (min_time - cust_dist[i][j]) * lexpr);
            //        lexpr.clear();
            //    }
            //}

            // Deadline flow - TTCF.4
            for (const int& i : customers) {
                for (const int& k : vehicles) {
                    const int idx = k * no_cust_d + i;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr += veh_arr_time[k] * gamma[idx];
                }
                for (const int& j : customers) {
                    if (i == j) continue;
                    const int idx = j * no_cust_d + i;
                    if (omega_indices.find(idx) != omega_indices.end()) lexpr += omega[idx];
                }
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (omega_indices.find(idx) != omega_indices.end()) lexpr2 += omega[idx];
                }
                model->addConstr(lexpr >= lexpr2);
                lexpr.clear();
                lexpr2.clear();
                lexpr3.clear();
            }
            // TTCF.5
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
                    if (i == j) continue;
                    const int idx = i * no_cust_d + j;
                    if (omega_indices.find(idx) == omega_indices.end()) continue;
                    if (main_type == "2I") {
                        if (psi_indices.find(idx) != psi_indices.end()) lexpr += psi[idx];
                    }
                    else {
                        for (const int& k : vehicles) {
                            const int chi_idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(chi_idx) != chi_indices.end()) lexpr += chi[chi_idx];
                        }
                    }
                    model->addConstr(omega[idx] <= cust_arr_time[i] * lexpr);
                    lexpr.clear();
                }
            }

            //for (const int& i : customers) {
            //    const int idx = i * no_cust_d + 0;
            //    if (sigma_indices.find(idx) != sigma_indices.end()) lexpr += sigma[idx];
            //    if (main_type == "2I") {
            //        const int idx = i * no_cust_d + 0;
            //        if (psi_indices.find(idx) != psi_indices.end()) {
            //            lexpr += cust_dist[i][0] * psi[idx];
            //        }
            //    }
            //    else {
            //        for (const int& k : vehicles) {
            //            const int idx = k * no_cust_d_squared + i * no_cust_d + 0;
            //            if (chi_indices.find(idx) != chi_indices.end()) {
            //                lexpr += cust_dist[i][0] * chi[idx];
            //            }
            //        }
            //    }
            //    if (omega_indices.find(idx) != omega_indices.end()) lexpr2 += omega[idx];
            //    model->addConstr(lexpr <= lexpr2);
            //    lexpr.clear();
            //    lexpr2.clear();
            //}
        }

        if (cap_type == "CAP") {
            for (const int& k : vehicles) {
                for (const int& i : customers) {
                    if (main_type == "2I") {
                        const int idx = k * no_cust_d + i;
                        if (xi_indices.find(idx) != xi_indices.end()) lexpr += demand[i] * xi[idx];
                    }
                    else {
                        for (const int& j : dest_and_cust) {
                            if (i == j) continue;
                            const int idx = k * no_cust_d_squared + i * no_cust_d + j;
                            if (chi_indices.find(idx) != chi_indices.end()) lexpr += demand[i] * chi[idx];

                        }
                    }
                }

                for (const int& j : customers) {
                    const int idx = k * no_cust_d + j;
                    if (gamma_indices.find(idx) != gamma_indices.end()) lexpr2 += gamma[idx];
                }

                model->addConstr(lexpr <= (capacity[k] - veh_occ[k]) * lexpr2);
                lexpr.clear();
                lexpr2.clear();
            }
        }
    }

    // Retrieves solution
    void retrieve_solution(std::set<std::tuple<int, int, double>>& solution, std::unordered_map<int, double>& tau_solution){
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
                    if (j == 0) solution.insert({ j, i, val });
                    else solution.insert({ i, j, val });
                }
            }
        }
        else if (main_type == "3I") {
            const int no_cust_d_squared = int(pow(no_cust_d, 2));
            for (const int& i : customers) {
                for (const int& j : dest_and_cust) {
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
                        if (j == 0) solution.insert({ j, i, val });
                        else solution.insert({ i, j, val });
                    }
                }
            }
        }
		for (const auto& [i, var] : tau) {
			double val = var.get(GRB_DoubleAttr_X);
			//if (val >= 1e-3) {
				tau_solution[i] = val;
			//}
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

    void add_rc_inequalities2(std::vector<RC_Inequality2>& violated_inequalties, bool& cuts_added) {
        cuts_added = false;
        const int no_cust_d_squared = int(pow(no_cust_d, 2));

        for (RC_Inequality2& cut : violated_inequalties) {
            GRBLinExpr lexpr, lexpr2;
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
                if (tau_indices.contains(i)) lexpr2 += (max_cap - demand[i]) * tau[i];
            }
            if (max_cap * lexpr.getValue() >= lexpr2.getValue() + 1e-3) {
                cuts_added = true;
                model->addConstr(max_cap * lexpr <= lexpr2);
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
    std::string graph_type;
    std::string main_type;
    std::string sec_type;
    std::string time_type;
    std::string cap_type;
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
    std::unordered_map<int, GRBVar> gamma;
    std::unordered_map<int, GRBVar> chi;
    std::unordered_map<int, GRBVar> phi;
    std::unordered_map<int, GRBVar> xi;
    std::unordered_map<int, GRBVar> psi;
    std::unordered_map<int, GRBVar> zeta;
    std::unordered_map<int, GRBVar> rho;
    std::unordered_map<int, GRBVar> rho_veh;
    std::unordered_map<int, GRBVar> zeta_multi;
    std::unordered_map<int, GRBVar> rho_multi;
    std::unordered_map<int, GRBVar> lambda;
    std::unordered_map<int, GRBVar> sigma;
    std::unordered_map<int, GRBVar> omega;
    std::unordered_map<int, GRBVar> Gamma;
    std::unordered_map<int, GRBVar> tau;


    int no_base_constraints = 0;
    bool preprocessed = false;

    //Valid inequalities
    std::vector<std::set<int>> rounded_cap;
    std::set<int> active_rc_ineqs;
    int total_rcc = 0;

    // Indices
    std::unordered_set<int> gamma_indices;
    std::unordered_set<int> chi_indices;
    std::unordered_set<int> phi_indices;
    std::unordered_set<int> xi_indices;
    std::unordered_set<int> psi_indices;
    std::unordered_set<int> zeta_indices;
    std::unordered_set<int> rho_indices;
    std::unordered_set<int> rho_veh_indices;
    std::unordered_set<int> zeta_multi_indices;
    std::unordered_set<int> rho_multi_indices;
    std::unordered_set<int> lambda_indices;
    std::unordered_set<int> sigma_indices;
    std::unordered_set<int> omega_indices;
    std::unordered_set<int> Gamma_indices;
    std::unordered_set<int> tau_indices;

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
    FMRSP(const std::string& filename, const std::string main, const std::string sec, const std::string time, const std::string cap, const int& threads, const int& time_lim) {
        std::cout << main << "_" << sec << "_" << time << "_" << cap << std::endl;

        // Load Data
        InstanceData inst_data = load_instance(filename);

        // Instance Properties
        name = filename;
        no_veh = inst_data.no_veh;
        no_cust = inst_data.no_cust;
        no_cust_d = no_cust + 1;
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
        sec_type = sec;
        time_type = time;
        cap_type = cap;

        construct_indices();
    }

    // Eliminate Variables
    void preprocess() {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Eliminating Variables..." << std::endl;
        std::unordered_set<int> gamma_indices_removed;
        std::unordered_set<int> chi_indices_removed;
        std::unordered_set<int> phi_indices_removed;
        std::unordered_set<int> xi_indices_removed;
        std::unordered_set<int> psi_indices_removed;
        std::unordered_set<int> zeta_indices_removed;
        std::unordered_set<int> rho_indices_removed;
        std::unordered_set<int> rho_veh_indices_removed;
        std::unordered_set<int> zeta_multi_indices_removed;
        std::unordered_set<int> rho_multi_indices_removed;
        std::unordered_set<int> lambda_indices_removed;
        std::unordered_set<int> sigma_indices_removed;
        std::unordered_set<int> omega_indices_removed;

        preprocessed = true;
        for (int i = 0; i < 2; i++) {

            preprocess_main(gamma_indices_removed, chi_indices_removed, phi_indices_removed, xi_indices_removed, psi_indices_removed, zeta_indices_removed, rho_indices_removed,
                rho_veh_indices_removed, zeta_multi_indices_removed, rho_multi_indices_removed, lambda_indices_removed, sigma_indices_removed, omega_indices_removed);

            preprocess_sec(gamma_indices_removed, chi_indices_removed, phi_indices_removed, xi_indices_removed, psi_indices_removed, zeta_indices_removed, rho_indices_removed,
                rho_veh_indices_removed, zeta_multi_indices_removed, rho_multi_indices_removed, lambda_indices_removed, sigma_indices_removed, omega_indices_removed);

            preprocess_time(gamma_indices_removed, chi_indices_removed, phi_indices_removed, xi_indices_removed, psi_indices_removed, zeta_indices_removed, rho_indices_removed,
                rho_veh_indices_removed, zeta_multi_indices_removed, rho_multi_indices_removed, lambda_indices_removed, sigma_indices_removed, omega_indices_removed);

            preprocess_cap(gamma_indices_removed, chi_indices_removed, phi_indices_removed, xi_indices_removed, psi_indices_removed, zeta_indices_removed, rho_indices_removed,
                rho_veh_indices_removed, zeta_multi_indices_removed, rho_multi_indices_removed, lambda_indices_removed, sigma_indices_removed, omega_indices_removed);
        }

        // Gamma
        int variables_removed = 0;
        if (gamma_indices.size() > 0) std::cout << "Gamma variables eliminated: " << std::fixed << std::setprecision(2) << gamma_indices_removed.size() / double((no_veh * no_cust_d)) * 100 << "%" << std::endl;
        if (gamma_indices_removed.size() > 0) {
            for (const int& idx : gamma_indices_removed) {
                gamma_indices.erase(idx);
            }
            variables_removed += gamma_indices_removed.size();
        }

        // Chi
        if (chi_indices.size() > 0) std::cout << "Chi variables eliminated: " << std::fixed << std::setprecision(2) << chi_indices_removed.size() / double(no_veh * (pow(no_cust, 2))) * 100 << "%" << std::endl;
        if (chi_indices_removed.size() > 0) {
            for (const int& idx : chi_indices_removed) {
                chi_indices.erase(idx);
            }
            variables_removed += chi_indices_removed.size();
        }

        // Phi
        if (phi_indices.size() > 0) std::cout << "Phi variables eliminated: " << std::fixed << std::setprecision(2) << phi_indices_removed.size() / double(no_cust) * 100 << "%" << std::endl;
        if (phi_indices_removed.size() > 0) {
            for (const int& idx : phi_indices_removed) {
                phi_indices.erase(idx);
            }
            variables_removed += phi_indices_removed.size();
        }

        // xi
        if (xi_indices.size() > 0) std::cout << "Xi variables eliminated: " << std::fixed << std::setprecision(2) << xi_indices_removed.size() / double(no_veh * no_cust_d) * 100 << "%" << std::endl;
        if (xi_indices_removed.size() > 0) {
            for (const int& idx : xi_indices_removed) {
                xi_indices.erase(idx);
            }
            variables_removed += xi_indices_removed.size();
        }

        //Psi
        if (psi_indices.size() > 0 ) std::cout << "Psi variables eliminated: " << std::fixed << std::setprecision(2) << psi_indices_removed.size() / double(pow(no_cust_d, 2) - no_cust_d) * 100 << "%" << std::endl;
        if (psi_indices_removed.size() > 0) {
            for (const int& idx : psi_indices_removed) {
                psi_indices.erase(idx);
            }
            variables_removed += psi_indices_removed.size();
        }

        // Zeta
        if (zeta_indices.size() > 0 && sec_type == "COCF") std::cout << "Zeta variables eliminated: " << std::fixed << std::setprecision(2) << zeta_indices_removed.size() / double(pow(no_cust, 2)) * 100 << "%" << std::endl;
        else if (zeta_indices.size() > 0 && sec_type == "TCF") std::cout << "Zeta variables eliminated: " << std::fixed << std::setprecision(2) << zeta_indices_removed.size() / double(pow(no_cust_d, 2) - no_cust_d) * 100 << "%" << std::endl;
        if (zeta_indices_removed.size() > 0) {
            for (const int& idx : zeta_indices_removed) {
                zeta_indices.erase(idx);
            }
            variables_removed += zeta_indices_removed.size();
        }

        // Rho
        if (rho_indices.size() > 0) std::cout << "Rho variables eliminated: " << std::fixed << std::setprecision(2) << rho_indices_removed.size() / double(pow(no_cust, 2)) * 100 << "%" << std::endl;
        if (rho_indices_removed.size() > 0) {
            for (const int& idx : rho_indices_removed) {
                rho_indices.erase(idx);
            }
            variables_removed += rho_indices_removed.size();
        }

        // Rho_veh
        if (rho_veh_indices.size() > 0) std::cout << "Rho_veh variables eliminated: " << std::fixed << std::setprecision(2) << rho_veh_indices_removed.size() / double((no_veh * no_cust)) * 100 << "%" << std::endl;
        if (rho_veh_indices_removed.size() > 0) {
            for (const int& idx : rho_veh_indices_removed) {
                rho_veh_indices.erase(idx);
            }
            variables_removed += rho_veh_indices_removed.size();
        }

        // Zeta_multi
        if (zeta_multi_indices.size() > 0) {
            std::cout << "Zeta_multi variables eliminated: " << std::fixed << std::setprecision(2);
            std::cout << zeta_multi_indices_removed.size() / double(no_cust * (pow(no_cust_d, 2) - no_cust_d)) * 100 << "%" << std::endl;
        }
        if (zeta_multi_indices_removed.size() > 0) {
            for (const int& idx : zeta_multi_indices_removed) {
                zeta_multi_indices.erase(idx);
            }
            variables_removed += zeta_multi_indices_removed.size();
        }

        // Rho_multi
        if (rho_multi_indices.size() > 0) {
            std::cout << "Rho_multi variables eliminated: " << std::fixed << std::setprecision(2);
            std::cout << rho_multi_indices_removed.size() / double(no_cust * (pow(no_cust_d, 2) - no_cust_d)) * 100 << "%" << std::endl;
        }
        if (rho_multi_indices_removed.size() > 0) {
            for (const int& idx : rho_multi_indices_removed) {
                rho_multi_indices.erase(idx);
            }
            variables_removed += rho_multi_indices_removed.size();
        }

        // Sigma
        if (sigma_indices.size() > 0) std::cout << "Sigma variables eliminated: " << std::fixed << std::setprecision(2) << sigma_indices_removed.size() / double(pow(no_cust, 2)) * 100 << "%" << std::endl;
        if (sigma_indices_removed.size() > 0) {
            for (const int& idx : sigma_indices_removed) {
                sigma_indices.erase(idx);
            }
            variables_removed += sigma_indices_removed.size();
        }

        // Omega
        if (omega_indices.size() > 0) std::cout << "Omega variables eliminated: " << std::fixed << std::setprecision(2) << omega_indices_removed.size() / double(pow(no_cust, 2)) * 100 << "%" << std::endl;
        if (omega_indices_removed.size() > 0) {
            for (const int& idx : omega_indices_removed) {
                omega_indices.erase(idx);
            }
            variables_removed += omega_indices_removed.size();
        }

        std::cout << "# Total variables removed: " << variables_removed << std::endl;
        auto stop = std::chrono::high_resolution_clock::now();
        graph_construction_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();
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
    void solve_root_relaxation(bool separate_rci = true, bool silent = true) {
        if (graph_construction_time == 0) preprocess();
        if (model_construction_time == 0) create_model(true, silent);

        auto start = std::chrono::high_resolution_clock::now();
        optimise_model();
        initial_lp = model->get(GRB_DoubleAttr_ObjVal);
        auto stop = std::chrono::high_resolution_clock::now();
        lp_solve_time += std::chrono::duration_cast<std::chrono::duration<double>>(stop - start).count();

        if (separate_rci) {
            std::cout << initial_lp << std::endl;
            std::set<std::tuple<int, int, double>> x_solution;
			std::unordered_map<int, double> tau_solution;
            std::vector<double> no_vehicles_used;
            std::vector<RC_Inequality> violated_inequalties;
            std::vector<RC_Inequality2> violated_inequalties2;

            double max_viol = 0;
            bool cuts_added, cuts_added2 = false;
            while (true) {
                retrieve_solution(x_solution, tau_solution);
                separate_rcc_exactly(violated_inequalties, max_viol, x_solution, customers, max_cap, demand, 1, RC_separation_time, env);
                //separate_rcc_exactly2(violated_inequalties2, max_viol, x_solution, tau_solution, customers, max_cap, demand, RC_separation_time, env);
                add_rc_inequalities(violated_inequalties, cuts_added);
                //add_rc_inequalities2(violated_inequalties2, cuts_added2);

                if (!cuts_added && !cuts_added2) break;

                x_solution.clear();
                tau_solution.clear();
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
};