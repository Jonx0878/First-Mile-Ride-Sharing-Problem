#pragma once
#include <gurobi_c++.h>
#include <iostream>
#include <unordered_set>
#include <set>
#include <vector>


struct RC_Inequality {
	std::unordered_set<int> set;
	double rhs;
};

struct RC_Inequality2 {
	std::unordered_set<int> set;
};

GRBModel* setup_rcc_model(
	GRBVar& alpha,
	std::vector<GRBVar>& xi,
	const std::set<std::tuple<int, int, double>>& solution,
	const std::set<int>& customers,
	const int& max_cap,
	const std::vector<int>& demand,
	const int& gcd,
	GRBEnv env
) {
	// Create Model
	GRBModel* model = new GRBModel(env);

	xi.resize(customers.size() + 1);
	std::vector<std::vector<GRBVar>> beta;
	beta.resize(customers.size() + 2);
	for (int i = 0; i <= customers.size();i++) {
		// Resize each inner vector
		beta[i].resize(customers.size() + 1);
	}

	// Change parameters
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	// Add Variables
	int total_demand = 0;
	for (const int& i : customers) total_demand += demand[i];
	const int ub = std::ceil(total_demand / float(max_cap));
	alpha = model->addVar(1, ub, 1, GRB_INTEGER, "alpha");
	for (int i : customers) {
		xi[i] = model->addVar(0, 1, -1, GRB_BINARY, "xi[" + std::to_string(i) + "]");
	}

	for (auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double obj = std::get<2>(edge);
		if (i > 0) {
			beta[i][j] = model->addVar(0, 1, obj, GRB_CONTINUOUS, "beta[" + std::to_string(i) + "," + std::to_string(j) + "]");
		}
	}

	// Add Constraints
	GRBLinExpr expr;
	for (int i : customers) {
		expr += demand[i] * xi[i];
	}
	model->addConstr(max_cap * (alpha - 1) + gcd, GRB_LESS_EQUAL, expr);

	expr.clear();

	for (int i : customers) {
		expr += xi[i];
	}
	model->addConstr(expr >= 2);
	expr.clear();

	for (auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);

		if (i > 0) {
			model->addConstr(beta[i][j], GRB_LESS_EQUAL, xi[i]);
			model->addConstr(beta[i][j], GRB_LESS_EQUAL, xi[j]);
		}
	}

	model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

	return model;
}

GRBModel* setup_fci_model(
	std::vector<GRBVar>& delta,
	const std::set<std::tuple<int, int, double>>& x_solution,
	const std::set<std::tuple<int, int, double>>& gamma_solution,
	const std::unordered_map<int, double>& tau_solution,
	const std::set<int>& customers,
	const std::set<int>& vehicles,
	const std::vector<int>& capacity,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_order,
	const std::vector<int>& demand,
	GRBEnv env
) {
	// Create Model
	GRBModel* model = new GRBModel(env);

	delta.resize(customers.size() + 1);
	std::vector<GRBVar> beta;
	std::vector<GRBVar> theta;
	std::vector<std::vector<GRBVar>> lambda;
	beta.resize(vehicles.size() + 1);
	theta.resize(vehicles.size() + 1);
	lambda.resize(customers.size() + 1);
	for (int i = 0; i <= customers.size(); i++) {
		// Resize each inner vector
		lambda[i].resize(customers.size() + 1);
	}

	// Change parameters
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	// Add Variables
	for (int i : customers) {
		delta[i] = model->addVar(0, 1, NULL, GRB_BINARY, "delta[" + std::to_string(i) + "]");
	}
	for (const int& k : vehicles) {
		beta[k] = model->addVar(0, 1, NULL, GRB_BINARY, "beta[" + std::to_string(k) + "]");
		theta[k] = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "theta[" + std::to_string(k) + "]");
	}

	for (auto& [i, j, sol] : x_solution) {
		if (i > 0) {
			lambda[i][j] = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "lambda[" + std::to_string(i) + "," + std::to_string(j) + "]");
		}
	}

	// Add Constraints
	model->update();
	GRBLinExpr lexpr, lexpr2;
	for (const int& k : vehicles) {
		for (const auto& [k2, i, sol] : gamma_solution) {
			if (veh_order[k2] != k || i == 0) continue;
			lexpr += (beta[k] + theta[k]) * (capacity[veh_order[k]] - veh_occ[veh_order[k]]) * sol;
		}
	}
	
	for (const auto& [i, sol] : tau_solution) {
		lexpr2 += demand[i] * sol * delta[i];
	}
	//std::cout << lexpr << " <= " << lexpr2 << std::endl;
	model->addConstr(lexpr, GRB_LESS_EQUAL, lexpr2); // 3.19
	lexpr.clear();
	lexpr2.clear();

	for (const int& k : vehicles) {
		model->addConstr(theta[k] <= 1 - beta[k]); // 3.22
		if (k == vehicles.size() - 1) continue;
		model->addConstr(beta[k] + theta[k] >= beta[k + 1] + theta[k + 1]); // 3.20
	}

	// Add SOS
	std::vector<double> weights(vehicles.size(), 1);
	model->addSOS(&theta[0], &weights[0], vehicles.size(), GRB_SOS_TYPE1);

	// Psi(S, S_bar)
	for (int i : customers) {
		lexpr += delta[i];
	}
	model->addConstr(lexpr >= 2);
	lexpr.clear();

	for (auto& [i, j, sol] : x_solution) {
		if (i > 0) {
			model->addConstr(lambda[i][j], GRB_GREATER_EQUAL, delta[j] - delta[i]);
			//model->addConstr(lambda[i][j], GRB_LESS_EQUAL, delta[j]);
			//model->addConstr(lambda[i][j], GRB_LESS_EQUAL, 1 - delta[i]);

		}
	}


	model->update();
	GRBLinExpr obj;
	for (const int& k : vehicles) {
		for (const auto& [k2, i, sol] : gamma_solution) {
			if (veh_order[k2] != k || i == 0) continue;
			obj += (beta[k] + theta[k]) * sol;
			obj -= delta[i] * sol;
		}
	}
	for (auto& [i, j, sol] : x_solution) {
		if (i > 0) {
			obj -= sol * lambda[i][j];
		}
	}
	model->setObjective(obj, GRB_MAXIMIZE);
	//model->optimize();
	//
	//for (auto& [i, j, sol] : gamma_solution) {
	//	std::cout << i << " " << j << " " << sol << std::endl;
	//}
	//for (int i : customers) {
	//	if (tau_solution.contains(i)) std::cout << tau_solution.at(i) << " ";
	//	std::cout << delta[i].get(GRB_StringAttr_VarName) << " " << delta[i].get(GRB_DoubleAttr_X) << std::endl;
	//}
	//for (const int& k : vehicles) {
	//	std::cout << beta[k].get(GRB_StringAttr_VarName) << " " << beta[k].get(GRB_DoubleAttr_X) << std::endl;
	//	std::cout << theta[k].get(GRB_StringAttr_VarName) << " " << theta[k].get(GRB_DoubleAttr_X) << std::endl;
	//}
	//for (auto& [i, j, sol] : x_solution) {
	//	if (i > 0) {
	//		std::cout << sol << " " << lambda[i][j].get(GRB_StringAttr_VarName) << " " << lambda[i][j].get(GRB_DoubleAttr_X) << std::endl;
	//	}
	//}
	return model;

}


GRBModel* setup_fci_model2(
	std::vector<GRBVar>& delta,
	const std::set<std::tuple<int, int, double>>& x_solution,
	const std::set<std::tuple<int, int, double>>& gamma_solution,
	const std::unordered_map<int, double>& tau_solution,
	const std::set<int>& customers,
	const std::set<int>& vehicles,
	const std::vector<int>& capacity,
	const std::vector<int>& veh_occ,
	const std::vector<int>& veh_order,
	const std::vector<int>& demand,
	GRBEnv env
) {
	// Create Model
	GRBModel* model = new GRBModel(env);

	delta.resize(customers.size() + 1);
	GRBVar alpha;
	GRBVar b;

	std::vector<GRBVar> beta;
	std::vector<GRBVar> theta;
	std::vector<std::vector<GRBVar>> lambda;
	std::vector<std::vector<GRBVar>> xi;
	beta.resize(vehicles.size() + 1);
	theta.resize(vehicles.size() + 1);
	lambda.resize(customers.size() + 1);
	xi.resize(customers.size() + 1);

	for (int i = 0; i <= customers.size(); i++) {
		// Resize each inner vector
		lambda[i].resize(customers.size() + 1);
		xi[i].resize(vehicles.size() + 1);
	}

	// Change parameters
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	// Add Variables
	alpha = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "alpha");
	b = model->addVar(0, 1, NULL, GRB_BINARY, "b");

	for (int i : customers) {
		delta[i] = model->addVar(0, 1, NULL, GRB_BINARY, "delta[" + std::to_string(i) + "]");
	}
	for (const int& k : vehicles) {
		beta[k] = model->addVar(0, 1, NULL, GRB_BINARY, "beta[" + std::to_string(k) + "]");
		theta[k] = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "theta[" + std::to_string(k) + "]");
	}

	for (auto& [i, j, sol] : x_solution) {
		if (i > 0) {
			lambda[i][j] = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "lambda[" + std::to_string(i) + "," + std::to_string(j) + "]");
		}
	}

	for (const int& i : customers) {
		for (const int& k : vehicles) {
			xi[i][veh_order[k]] = model->addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS, "xi[" + std::to_string(i) + "," + std::to_string(veh_order[k]) + "]");
		}
	}
	// Add Constraints
	model->update();
	GRBLinExpr lexpr, lexpr2;
	// 3.40
	for (const int& k : vehicles) {
		for (const auto& [k2, i, sol] : gamma_solution) {
			if (veh_order[k2] != k || i == 0) continue;
			lexpr += (delta[i] + xi[i][k]) * (capacity[veh_order[k]] - veh_occ[veh_order[k]]) * sol;
			//std::cout << i << " " << k << " " << k2 << " " << (delta[i] + xi[i][k]) * (capacity[k] - veh_occ[k]) * sol << std::endl;
		}
	}
	//std::cout << lexpr << std::endl;
	model->addConstr(lexpr, GRB_LESS_EQUAL, alpha);
	lexpr.clear();

	// 3.41 / 3.42
	double sum_1 = 0, sum_2 = 0;
	for (const auto& [i, sol] : tau_solution) {
		sum_1 += demand[i] * sol;
		lexpr += demand[i] * sol * delta[i];
	}
	for (const int& k : vehicles) {
		for (const auto& [k2, i, sol] : gamma_solution) {
			if (veh_order[k2] != k || i == 0) continue;
			sum_2 += (capacity[veh_order[k]] - veh_occ[veh_order[k]]) * sol;
			lexpr2 += delta[i] * (capacity[veh_order[k]] - veh_occ[veh_order[k]]) * sol;
		}
	}
	double M = std::max(sum_1, sum_2);
	model->addConstr(alpha <= lexpr + M * b);
	model->addConstr(alpha <= lexpr2 + M * (1 - b));
	lexpr.clear();
	lexpr2.clear();

	for (const int& k : vehicles) {
		model->addConstr(theta[k] <= 1 - beta[k]); // 3.47
		for (const int& i : customers) {
			model->addConstr(xi[i][k] <= 1 - delta[i]); // 3.45
			model->addConstr(xi[i][k] <= beta[k] + theta[k]); // 3.45
			model->addConstr(xi[i][k] >= beta[k] + theta[k] - delta[i]); // 3.45

		}
		if (k == vehicles.size() - 1) continue;
		// 3.43
		model->addConstr(beta[k] + theta[k] >= beta[k + 1] + theta[k + 1]); // 3.20
	}

	// Add SOS - 3.44
	std::vector<double> weights(vehicles.size());
	std::fill(weights.begin(), weights.end(), 1);
	model->addSOS(&theta[0], &weights[0], vehicles.size(), GRB_SOS_TYPE1);

	// Psi(S, S_bar)
	for (int i : customers) {
		lexpr += delta[i];
	}
	model->addConstr(lexpr >= 2);
	lexpr.clear();

	for (auto& [i, j, sol] : x_solution) {
		if (i > 0) {
			model->addConstr(lambda[i][j], GRB_GREATER_EQUAL, delta[j] - delta[i]);
			//model->addConstr(lambda[i][j], GRB_LESS_EQUAL, delta[j]);
			//model->addConstr(lambda[i][j], GRB_LESS_EQUAL, 1 - delta[i]);

		}
	}


	model->update();
	GRBLinExpr obj;
	for (const int& k : vehicles) {
		for (const auto& [k2, i, sol] : gamma_solution) {
			if (veh_order[k2] != k || i == 0) continue;
			obj += xi[i][k] * sol;
			//obj -= delta[i] * sol;
		}
	}
	for (auto& [i, j, sol] : x_solution) {
		if (i > 0) {
			obj -= sol * lambda[i][j];
		}
	}
	//std::cout << obj << std::endl;
	model->setObjective(obj, GRB_MAXIMIZE);

	//model->optimize();
	//
	//for (auto& [i, j, sol] : gamma_solution) {
	//	std::cout << i << " " << j << " " << sol << std::endl;
	//}
	//std::cout << alpha.get(GRB_StringAttr_VarName) << " " << alpha.get(GRB_DoubleAttr_X) << std::endl;

	//for (int i : customers) {
	//	if (tau_solution.contains(i)) std::cout << tau_solution.at(i) << " ";
	//	std::cout << delta[i].get(GRB_StringAttr_VarName) << " " << delta[i].get(GRB_DoubleAttr_X) << std::endl;
	//}
	//for (const int& k : vehicles) {
	//	std::cout << beta[k].get(GRB_StringAttr_VarName) << " " << beta[k].get(GRB_DoubleAttr_X) << std::endl;
	//	std::cout << theta[k].get(GRB_StringAttr_VarName) << " " << theta[k].get(GRB_DoubleAttr_X) << std::endl;
	//}
	//for (auto& [i, j, sol] : x_solution) {
	//	if (i > 0) {
	//		std::cout << sol << " " << lambda[i][j].get(GRB_StringAttr_VarName) << " " << lambda[i][j].get(GRB_DoubleAttr_X) << std::endl;
	//	}
	//}
	//for (const int& i : customers) {
	//	for (const int& k : vehicles) {
	//		std::cout << xi[i][k].get(GRB_StringAttr_VarName) << " " << xi[i][k].get(GRB_DoubleAttr_X) << std::endl;
	//	}
	//}

	return model;
}


void separate_rcc_exactly(
	std::vector<RC_Inequality>& new_cuts,
	double& max_violation,
	const std::set<std::tuple<int, int, double>>& solution,
	const std::set<int>& customers,
	const int& max_cap,
	const std::vector<int>& demand,
	const int& gcd,
	double& time,
	GRBEnv env
) {
	auto start = std::chrono::high_resolution_clock::now();
	max_violation = 0;
	GRBVar alpha;
	std::vector<GRBVar> xi;

	// Setup separation model
	GRBModel* model = setup_rcc_model(alpha, xi, solution, customers, max_cap, demand, gcd, env);

	//mycallback cb = mycallback(new_cuts, 25, customers, alpha, xi, max_violation);
	//model->setCallback(&cb);
	//model->set(GRB_IntParam_LazyConstraints, 1);
	model->optimize();

	max_violation = std::fmax(max_violation, model->get(GRB_DoubleAttr_ObjVal));
	for (int sol = 0; sol < model->get(GRB_IntAttr_SolCount); ++sol) {
		model->set(GRB_IntParam_SolutionNumber, sol);
		if (model->get(GRB_DoubleAttr_PoolObjVal) >= 1e-3) {
			std::unordered_set<int> S;
			for (int i : customers) {
				if (xi[i].get(GRB_DoubleAttr_Xn) >= 1 - 1e-3) {
					S.insert(i);
				}
			}
			new_cuts.emplace_back(RC_Inequality(S, std::round(S.size() - alpha.get(GRB_DoubleAttr_Xn))));
		}
	}
	std::cout << max_violation << std::endl;

	delete model;

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
}


double separate_fci_exactly(
	std::vector<RC_Inequality2>& new_cuts,
	double& max_violation,
	const std::set<std::tuple<int, int, double>>& x_solution,
	const std::set<std::tuple<int, int, double>>& gamma_solution,
	const std::unordered_map<int, double>& tau_solution,
	const std::set<int>& customers,
	const std::set<int>& vehicles,
	const std::vector<int>& capacity,
	const std::vector<int>& veh_occ,
	const std::vector<int>& demand,
	double& time,
	GRBEnv env
) {
	auto start = std::chrono::high_resolution_clock::now();
	max_violation = 0;
	GRBVar alpha;
	std::vector<GRBVar> delta;

	// Determine vehicle order
	std::vector<std::pair<int, int>> veh_order_pairs;
	for (const int& k : vehicles) {
		veh_order_pairs.emplace_back(k, capacity[k] - veh_occ[k]);
	}
	std::sort(veh_order_pairs.begin(), veh_order_pairs.end(),
		[](const std::pair<int, int>& a, const std::pair<int, int>& b) {
			return a.second > b.second; // Sort in descending order
		});
	std::vector<int> veh_order(vehicles.size());
	for (int k = 0; k < vehicles.size(); k++) {
		veh_order[k] = veh_order_pairs[k].first;
		//std::cout << k << " " << veh_order[k] << std::endl;
	}


	// Setup separation model
	GRBModel* model = setup_fci_model(delta, x_solution, gamma_solution, tau_solution, customers, vehicles, capacity, veh_occ, veh_order, demand, env);

	model->optimize();

	max_violation = std::fmax(max_violation, model->get(GRB_DoubleAttr_ObjVal));
	for (int sol = 0; sol < model->get(GRB_IntAttr_SolCount); ++sol) {
		model->set(GRB_IntParam_SolutionNumber, sol);
		if (model->get(GRB_DoubleAttr_PoolObjVal) >= 1e-3) {
			std::unordered_set<int> S;
			for (int i : customers) {
				if (delta[i].get(GRB_DoubleAttr_Xn) >= 1 - 1e-3) {
					S.insert(i);
					//std::cout << i << " ";
				}
			}
			//std::cout << std::endl;
			new_cuts.emplace_back(RC_Inequality2(S));
		}
	}
	std::cout << max_violation << std::endl;
	delete model;

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
	return max_violation;
}