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

GRBModel* setup_rcc_model2(
	GRBVar& alpha,
	std::vector<GRBVar>& xi,
	const std::set<std::tuple<int, int, double>>& x_solution,
	const std::unordered_map<int, double>& tau_solution,
	const std::set<int>& customers,
	const int& max_cap,
	const std::vector<int>& demand,
	GRBEnv env
) {
	// Create Model
	GRBModel* model = new GRBModel(env);

	xi.resize(customers.size() + 1);
	std::vector<std::vector<GRBVar>> beta;
	beta.resize(customers.size() + 2);
	for (int i = 0; i <= customers.size(); i++) {
		// Resize each inner vector
		beta[i].resize(customers.size() + 1);
	}

	// Change parameters
	model->getEnv().set(GRB_IntParam_OutputFlag, 0);

	// Add Variables
	for (int i : customers) {
		xi[i] = model->addVar(0, 1, -(max_cap - demand[i])*tau_solution.at(i), GRB_BINARY, "xi[" + std::to_string(i) + "]");
	}

	for (auto& edge : x_solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double obj = std::get<2>(edge);
		if (i > 0) {
			beta[i][j] = model->addVar(0, 1, obj * max_cap, GRB_CONTINUOUS, "beta[" + std::to_string(i) + "," + std::to_string(j) + "]");
		}
	}

	// Add Constraints
	GRBLinExpr expr;
	for (int i : customers) {
		expr += xi[i];
	}
	model->addConstr(expr >= 2);
	expr.clear();

	for (auto& edge : x_solution) {
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


void separate_rcc_exactly2(
	std::vector<RC_Inequality2>& new_cuts,
	double& max_violation,
	const std::set<std::tuple<int, int, double>>& x_solution,
	const std::unordered_map<int, double>& tau_solution,
	const std::set<int>& customers,
	const int& max_cap,
	const std::vector<int>& demand,
	double& time,
	GRBEnv env
) {
	auto start = std::chrono::high_resolution_clock::now();
	max_violation = 0;
	GRBVar alpha;
	std::vector<GRBVar> xi;

	// Setup separation model
	GRBModel* model = setup_rcc_model2(alpha, xi, x_solution, tau_solution, customers, max_cap, demand, env);

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
			new_cuts.emplace_back(RC_Inequality2(S));
		}
	}
	std::cout << max_violation << std::endl;
	delete model;

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
	time += duration.count();
}