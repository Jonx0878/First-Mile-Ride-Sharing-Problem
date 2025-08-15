#pragma once

#include <algorithm>
#include <set>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>
#include <gurobi_c++.h>

#include "utils.h"

void calculate_rc_gamma_kd(
	double& rc,
	std::vector<GRBConstr>& col_constrs,
	std::vector<double>& col_coeffs,
	const int& k,
	const std::vector<std::vector<int>>& veh_dist,
	const std::set<int>& customers,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const int& no_veh,
	const int& no_cust,
	GRBConstr* constrs,
	double* dual_values
)
{
	// Cost
	rc = veh_dist[k][0];

	// 1.2
	rc -= dual_values[k];
	col_constrs.push_back(constrs[k]);
	col_coeffs.push_back(1);

	// 1.5
	int idx = 3 * no_veh + k;
	rc -= -veh_occ[k] * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(-veh_occ[k]);

	// 1.8
	for (int i : customers) {
		idx = 4 * no_veh + no_cust * (1 + no_veh + k) + i - 1;
		double coeff = veh_dist[k][0] - veh_arr_time[k];
		rc -= coeff * dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(coeff);
	}
}

void calculate_rc_gamma_kj(
	double& rc,
	std::vector<GRBConstr>& col_constrs,
	std::vector<double>& col_coeffs,
	const int& k,
	const int& j,
	const std::vector<std::vector<int>>& veh_dist,
	const std::vector<int>& demand,
	const std::set<int>& customers,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const int& no_veh,
	const int& no_cust,
	const int& capacity,
	GRBConstr* constrs,
	double* dual_values
)
{
	// Cost
	rc = veh_dist[k][j]; 

	// 1.2
	rc -= dual_values[k];
	col_constrs.push_back(constrs[k]);
	col_coeffs.push_back(1);

	// 1.3
	int idx = no_veh + k;
	rc -= dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(1);

	// 1.4
	idx = 2 * no_veh + k;
	double coeff = -(capacity - veh_occ[k]);
	rc -= coeff * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(coeff);

	// 1.5
	idx = 3 * no_veh + k;
	rc -= -veh_occ[k] * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(-veh_occ[k]);

	// 1.7
	idx = 4 * no_veh + (1 + k) * no_cust + j - 1;
	rc -= dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(1);


	// 1.8
	for (int i : customers) {
		idx = 4 * no_veh + no_cust * (1 + no_veh + k) + i - 1;
		coeff = veh_dist[k][j] - veh_arr_time[k];
		rc -= coeff * dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(coeff);
	}

	// 1.9
	idx = 4 * no_veh + (1 + 2 * no_veh + k) * no_cust + j;
	coeff = (veh_occ[k] + demand[j]);
	rc -= coeff * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(coeff);
}


void calculate_rc_chi_id(
	double& rc,
	std::vector<GRBConstr>& col_constrs,
	std::vector<double>& col_coeffs,
	const int& i,
	const int& k,
	const std::vector<std::vector<int>>& cust_dist,
	const std::vector<int>& demand,
	const std::set<int>& customers,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const int& no_veh,
	const int& no_cust,
	GRBConstr* constrs,
	double* dual_values
) {
	// Cost 
	rc = cust_dist[i][0] - demand[i] * cust_dist[i][0]; // Cost

	// 1.3
	int idx = no_veh + k;
	rc -= -dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(-1);

	// 1.4
	idx = 2 * no_veh + k;
	rc -= demand[i] * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(demand[i]);

	// 1.6
	idx = 4 * no_veh + i - 1;
	rc -= dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(1);

	// 1.7
	idx = 4 * no_veh + (1 + k) * no_cust + i - 1;
	rc -= -dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(-1);

	// 1.8
	for (int l : customers) {
		idx = 4 * no_veh + no_cust * (1 + no_veh + k) + l - 1;
		double coeff = cust_dist[i][0];
		if (l == i && cust_arr_time[i] < veh_arr_time[k]) coeff += veh_arr_time[k] - cust_arr_time[i];
		rc -= coeff * dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(coeff);
	}
}

void calculate_rc_chi_ij(
	double& rc,
	std::vector<GRBConstr>& col_constrs,
	std::vector<double>& col_coeffs,
	const int& i,
	const int& j,
	const int& k,
	const double& max_cap,
	const std::vector<std::vector<int>>& cust_dist,
	const std::vector<int>& demand,
	const std::set<int>& customers,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const int& no_veh,
	const int& no_cust,
	GRBConstr* constrs,
	double* dual_values,
	const std::vector<std::set<int>>& rounded_cap,
	std::unordered_map<int, int>& rc_indices
	) {
	// Cost 
	rc = cust_dist[i][j] - demand[i] * cust_dist[i][0]; // Cost

	// 1.4
	int idx = 2 * no_veh + k;
	rc -= demand[i] * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(demand[i]);

	// 1.6
	idx = 4 * no_veh + i - 1;
	rc -= dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(1);

	// 1.7
	idx = 4 * no_veh + (1 + k) * no_cust + j - 1;
	rc -= dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(1);

	idx = 4 * no_veh + (1 + k) * no_cust + i - 1;
	rc -= -dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(-1);

	// 1.8
	for (int l : customers) {
		idx = 4 * no_veh + no_cust * (1 + no_veh + k) + l - 1;
		double coeff = cust_dist[i][j];
		if (l == i && cust_arr_time[i] < veh_arr_time[k]) coeff += veh_arr_time[k] - cust_arr_time[i];
		rc -= coeff * dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(coeff);
	}

	// 1.10
	idx = 4 * no_veh + (1 + 3 * no_veh) * no_cust + (i - 1) * (no_cust - 1) + j - 1;
	if (j > i) idx--;
	double coeff = std::max(max_cap, double(demand[i]));
	rc -= coeff * dual_values[idx];
	col_constrs.push_back(constrs[idx]);
	col_coeffs.push_back(coeff);

	// RC
	if (i < j)	idx = (2 * no_cust - 3 - (i - 1)) * (i - 1) / 2.0 + (j - 1) - 1;
	else idx = (2 * no_cust - 3 - (j - 1)) * (j - 1) / 2.0 + (i - 1) - 1;
	for (int rc_idx : rounded_cap[idx]) {
		if (rc_indices.count(rc_idx) == 1) {
			int constr_idx = rc_indices[rc_idx];
			rc -= dual_values[constr_idx];
			col_constrs.push_back(constrs[constr_idx]);
			col_coeffs.push_back(1);
		}
	}
}

void process_valid_inequalties(std::unordered_map<int, int>& rc_indices, std::string* constr_names, const int& no_base_constraints, const int& no_constrs) {
	for (int i = no_base_constraints; i < no_constrs; i++) {
		std::string name = constr_names[i];
		if (name.rfind("RC_", 0) == 0) {
			int idx = std::stoi(name.substr(3));
			rc_indices[idx] = i;
		}
	}
}

void add_variables_cg(
	bool& new_vars,
	const double& max_cap,
	const std::vector<std::vector<int>>& veh_dist,
	const std::vector<std::vector<int>>& cust_dist,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const std::vector<int>& demand,
	const std::set<int>& customers,
	const int& no_veh,
	const int& no_cust,
	const int& no_cust_d,
	const int& capacity,
	const int& no_veh_vars,
	const int& no_total_vars,
	const int& no_base_constraints,
	GRBModel* model,
	std::vector<GRBVar>& vars,
	const std::set<int>& var_indices,
	std::set<int>& var_indices_subset,
	const std::vector<std::set<int>>& rounded_cap,
	const bool& add_all
	) {
	int no_new_vars = 0;
	std::unordered_map<int, int> rc_indices;

	std::set<int> var_indices_unused;
	std::set_difference(var_indices.begin(), var_indices.end(), var_indices_subset.begin(), var_indices_subset.end(), std::inserter(var_indices_unused, var_indices_unused.begin()));
	GRBConstr* constrs = model->getConstrs();
	int no_constr = model->get(GRB_IntAttr_NumConstrs);
	double* dual_values = model->get(GRB_DoubleAttr_Pi, constrs, no_constr);
	std::string* constr_names = model->get(GRB_StringAttr_ConstrName, constrs, no_constr);

	process_valid_inequalties(rc_indices, constr_names, no_base_constraints, no_constr);

	for (int idx : var_indices_unused) {
		int i, j, k;
		double rc = 0;
		double obj = 0;
		std::string name;
		std::vector<GRBConstr> col_constrs;
		std::vector<double> col_coeffs;
		if (idx <= no_veh_vars) {
			k = int(floor((idx - 1) / (no_cust + 1.0)));
			j = idx - k * (no_cust + 1) - 1;
			obj = veh_dist[k][j];
			name = "gamma[" + itos(k) + "][" + itos(j) + "]";
			if (j == 0) {
				calculate_rc_gamma_kd(rc, col_constrs, col_coeffs, k, veh_dist, customers, veh_occ, veh_arr_time, cust_arr_time, no_veh, no_cust, constrs, dual_values);
			}
			else {
				calculate_rc_gamma_kj(rc, col_constrs, col_coeffs, k, j, veh_dist, demand, customers, veh_occ, veh_arr_time, cust_arr_time, no_veh, no_cust, capacity, constrs, dual_values);
			}
		}
		else {
			k = int(floor((idx - 1 - no_veh * no_cust_d) / (pow(no_cust_d, 2) - no_cust - 1.0)));
			i = int(floor(
				(idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1.0)) / no_cust_d
				+ 1));
			j = int((idx - 1 - no_veh * no_cust_d - k * (pow(no_cust_d, 2) - no_cust - 1) - (i - 1) * no_cust_d));
			obj = cust_dist[i][j] - demand[i] * cust_dist[i][0];
			name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
			if (j == 0) {
				calculate_rc_chi_id(rc, col_constrs, col_coeffs, i, k, cust_dist, demand, customers, veh_arr_time, cust_arr_time, no_veh, no_cust, constrs, dual_values);
			}
			else {
				calculate_rc_chi_ij(rc, col_constrs, col_coeffs, i, j, k, max_cap, cust_dist, demand, customers, veh_arr_time, cust_arr_time, no_veh, no_cust, constrs, dual_values, rounded_cap, rc_indices);
			}
		}
		if (rc <= -0.001 || add_all) {
			GRBColumn col;
			col.addTerms(&col_coeffs[0], &col_constrs[0], int(col_coeffs.size()));
			vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, col, name);
			var_indices_subset.insert(idx);
			no_new_vars++;
		}
		

		//if (no_new_vars >= no_total_vars / 10) break;
	}
	new_vars = (no_new_vars > 0);
}

/*------------------------------------------------------------------------------------------*/
/*										EVENT BASED											*/
/*------------------------------------------------------------------------------------------*/

void calculate_reduced_cost_EB(
	double& rc,
	std::vector<GRBConstr>& col_constrs,
	std::vector<double>& col_coeffs,
	const int& k,
	const int& i,
	const int& j,
	const int& no_veh,
	const int& no_cust,
	const int& no_cust_d,
	const std::vector<std::vector<int>>& veh_dist,
	const std::vector<std::vector<int>>& cust_dist,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const std::set<int>& customers,
	std::vector<std::vector<int>>& vertices,
	GRBConstr* constrs,
	double* dual_values
) {
	int idx;
	double coeff;
	// 1.2
	if (i == 0) {
		idx = k;
		rc -= dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(1);
	}

	// 1.3
	if (i == 0 && j > 0) {
		idx = no_veh + k;
		rc -= dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(1);
	}
	else if (i > 0 && j == 0) {
		idx = no_veh + k;
		rc -= -dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(-1);
	}

	// 1.5
	if (i == 0) {
		idx = 2 * no_veh + k;
		rc -= -veh_occ[k] * dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(-veh_occ[k]);
	}

	// 1.6
	if (i > 0) {
		idx = 3 * no_veh + vertices[i][0] - 1;
		rc -= dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(1);
	}

	// 1.7
	if (j > 0) {
		idx = 3 * no_veh + no_cust + k * (vertices.size() - 1) + j - 1;
		rc -= dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(1);
	}
	if (i > 0) {
		idx = 3 * no_veh + no_cust + k * (vertices.size() - 1) + i - 1;
		rc -= -dual_values[idx];
		col_constrs.push_back(constrs[idx]);
		col_coeffs.push_back(-1);
	}

	// 1.8
	if (i == 0) {
		for (int l : customers) {
			idx = 3 * no_veh + (1 + k) * no_cust + no_veh * (vertices.size() - 1) + l - 1;
			coeff = veh_dist[k][vertices[j][0]] - veh_arr_time[k];
			rc -= coeff * dual_values[idx];
			col_constrs.push_back(constrs[idx]);
			col_coeffs.push_back(coeff);
		}
	}
	else {
		for (int l : customers) {
			idx = 3 * no_veh + (1 + k) * no_cust + no_veh * (vertices.size() - 1) + l - 1;
			coeff = cust_dist[vertices[i][0]][vertices[j][0]]; 
			if (l == vertices[i][0] && cust_arr_time[l] < veh_arr_time[k]) coeff += veh_arr_time[k] - cust_arr_time[l];
			rc -= coeff * dual_values[idx];
			col_constrs.push_back(constrs[idx]);
			col_coeffs.push_back(coeff);
		}
	}
}

void add_variables_EB(
	bool& new_vars,
	const int& no_base_constraints,
	const std::vector<std::vector<int>>& veh_dist,
	const std::vector<std::vector<int>>& cust_dist,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_arr_time,
	const std::vector<int>& cust_arr_time,
	const std::vector<int>& demand,
	const std::set<int>& customers,
	const int& no_veh,
	const int& no_cust,
	const int& no_cust_d,
	std::vector<std::vector<int>>& vertices,
	std::vector<std::pair<int, int>>& edges,
	std::vector<int>& first_var_index,
	GRBModel* model,
	std::vector<GRBVar>& vars,
	const std::set<int>& var_indices,
	std::set<int>& var_indices_subset,
	const std::vector<std::set<int>>& rounded_cap,
	const bool& add_all,
	int& outer_iterations
) {
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
	std::vector<int> new_vars_veh;
	for (int k = 0; k < first_var_index.size() - 1; k++) {
		new_vars_veh.emplace_back(0);
	}

	std::vector<GRBColumn> cols(var_indices.size() + 1);
	std::vector<double> rcs(var_indices.size() + 1);
	std::vector<int> indices;

	for (int idx : var_indices_unused) {
		int k = -1;
		while (idx >= first_var_index[k + 1]) k++;
		//if (new_vars_veh[k] >= 3000) continue;
		int i = edges[idx].first;
		int j = edges[idx].second;
		int i_cust = vertices[i][0];
		int j_cust = vertices[j][0];

		double obj;
		std::string name;
		if (i == 0) {
			obj = veh_dist[k][j_cust];
			name = "gamma[" + itos(k) + "][" + itos(j) + "]";
		}
		else {
			obj = cust_dist[i_cust][j_cust] - demand[i_cust] * cust_dist[i_cust][0];
			name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
		}

		double rc = obj;
		std::vector<GRBConstr> col_constrs;
		std::vector<double> col_coeffs;

		calculate_reduced_cost_EB(rc, col_constrs, col_coeffs, k, i, j, no_veh, no_cust, no_cust_d, veh_dist, cust_dist, veh_occ, veh_arr_time, cust_arr_time, customers, vertices, constrs, dual_values);

		if (rc <= -0.001 || add_all) {
			GRBColumn col;
			col.addTerms(&col_coeffs[0], &col_constrs[0], int(col_coeffs.size()));
			cols[idx] = col;
			rcs[idx] = rc;
			indices.emplace_back(idx);
			//double obj;
			//std::string name;
			//if (i == 0) {
			//	obj = veh_dist[k][j_cust];
			//	name = "gamma[" + itos(k) + "][" + itos(j) + "]";
			//}
			//else {
			//	obj = cust_dist[i_cust][j_cust] - demand[i_cust] * cust_dist[i_cust][0];
			//	name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
			//}

			//GRBColumn col;
			//col.addTerms(&col_coeffs[0], &col_constrs[0], int(col_coeffs.size()));
			//vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, col, name);
			//var_indices_subset.insert(idx);
			//no_new_vars++;
			//new_vars_veh[k]++;
		}
		//if (no_new_vars > 10000) break;
	}

	std::sort(indices.begin(), indices.end(),
		[&](int A, int B) {
			return rcs[A] < rcs[B];
		});

	int max_new_vars = indices.size();
	// std::max(- 50 * outer_iterations * outer_iterations + 1000 * outer_iterations + 1000, 1000)
	if (!add_all) max_new_vars = std::min(max_new_vars, int(pow(outer_iterations, 0.5) * 1000) );
	for (int l = 0; l < max_new_vars; l++) {
		int idx = indices[l];
		int k = -1;
		while (idx >= first_var_index[k + 1]) k++;
		//if (new_vars_veh[k] >= 3000) continue;
		int i = edges[idx].first;
		int j = edges[idx].second;
		int i_cust = vertices[i][0];
		int j_cust = vertices[j][0];
		double obj = 0;
		std::string name;
		if (i == 0) {
			obj = veh_dist[k][j_cust];
			name = "gamma[" + itos(k) + "][" + itos(j) + "]";
		}
		else {
			obj = cust_dist[i_cust][j_cust] - demand[i_cust] * cust_dist[i_cust][0];
			name = "chi[" + itos(i) + "][" + itos(j) + "][" + itos(k) + "]";
		}
		vars[idx] = model->addVar(0.0, 1.0, obj, GRB_CONTINUOUS, cols[idx], name);
		var_indices_subset.insert(idx);
		no_new_vars++;
		new_vars_veh[k]++;
	}
	new_vars = (no_new_vars > 0);
}

void remove_constraints_from_model(GRBModel* model, const int& no_base_constraints) {
	GRBConstr* constrs = model->getConstrs();
	GRBConstr constr;
	for (int i = no_base_constraints; i < model->get(GRB_IntAttr_NumConstrs); i++) {
		constr = constrs[i];
		double slack = constr.get(GRB_DoubleAttr_Slack);
		if (slack >= 0.001) {
			model->remove(constr);
		}
	}
}

void remove_variables_EB(GRBModel* model, std::vector<GRBVar> vars, std::set<int>& var_indices_subset, std::vector<int>& first_var_index) {
	for (int idx : var_indices_subset) {
		if (std::find(first_var_index.begin(), first_var_index.end(), idx) != first_var_index.end()) break; // EB Specific
		if (vars[idx].get(GRB_DoubleAttr_X) <= 0.0001 && vars[idx].get(GRB_DoubleAttr_RC) > 0) {
			model->remove(vars[idx]);
			var_indices_subset.erase(idx);
		}
	}

}