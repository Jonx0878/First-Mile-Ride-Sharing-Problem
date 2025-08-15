#pragma once
#include <set>
#include <unordered_set>

#include "gurobi_c++.h"

struct BBV_RC_Inequality {
	const std::unordered_set<int> set;
	const int constant;
	const int veh_type;
};

void create_separation_model(
	GRBModel* model,
	GRBVar& alpha,
	std::vector<GRBVar>& xi,
	const std::set<std::tuple<int, int, double>>& solution,
	const std::vector<double>& no_vehicles_used, // Chceck whether int or fractional values
	const std::vector<int>& vehicle_sizes,
	const std::set<int>& customers,
	const int& veh_type,
	const std::vector<int>& demand,
	const int& gcd
) {

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
	int total_demand = 0;
	for (const int& i : customers) total_demand += demand[i];
	const int ub = std::ceil(total_demand / float(vehicle_sizes[vehicle_sizes.size() - 1]));
	alpha = model->addVar(1, ub, 1, GRB_INTEGER, "alpha");

	for (const int& i : customers) {
		xi[i] = model->addVar(0, 1, -1, GRB_BINARY, "xi[" + std::to_string(i) + "]");
	}

	for (auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		double obj_val = std::get<2>(edge);
		if (j > 0) {
			beta[i][j] = model->addVar(0, 1, obj_val, GRB_CONTINUOUS, "beta[" + std::to_string(i) + "," + std::to_string(j) + "]");
			//obj -= obj_val * (xi[i] + xi[j] - 2 * beta[i][j]);
		}
		//else {
		//	obj -= obj_val * xi[i];
		//}

	}

	// Add Constraints
	GRBLinExpr expr;
	for (int i : customers) {
		expr += demand[i] * xi[i];
	}
	model->addConstr(vehicle_sizes[veh_type +1] * (alpha - 1) + gcd, GRB_LESS_EQUAL, expr);
	expr.clear();

	for (int i : customers) {
		expr += xi[i];
	}
	model->addConstr(expr >= 2);
	expr.clear();

	for (auto& edge : solution) {
		int i = std::get<0>(edge);
		int j = std::get<1>(edge);
		if (j > 0) {
			model->addConstr(beta[i][j], GRB_LESS_EQUAL, xi[i]);
			model->addConstr(beta[i][j], GRB_LESS_EQUAL, xi[j]);
			model->addConstr(beta[i][j], GRB_GREATER_EQUAL, xi[i] + xi[j] - 1);

		}
	}

	// Set Objective
	model->update();
	GRBQuadExpr obj = model->getObjective();

	for (int h = 0; h <= veh_type; h++) {
		obj -= std::ceil(vehicle_sizes[h] / float(vehicle_sizes[veh_type + 1]) - 1) * no_vehicles_used[h];
	}
	//std::cout << obj << std::endl;
	model->setObjective(obj, GRB_MAXIMIZE);
}

void separate_BBV_RCI(
	const std::set<std::tuple<int, int, double>>& solution,
	const std::vector<double>& no_vehicles_used,
	const std::vector<int>& vehicle_sizes,
	const std::set<int>& customers,
	const std::vector<int>& demand,
	const int& gcd,
	GRBEnv env,
	std::vector<BBV_RC_Inequality>& violated_inequalities
) {
	for (int h = 0; h < vehicle_sizes.size() - 1; h++) {
		GRBModel* model = new GRBModel(env);
		GRBVar alpha;
		std::vector<GRBVar> xi;

		create_separation_model(model, alpha, xi, solution, no_vehicles_used, vehicle_sizes,
			customers, h, demand, gcd);

		model->optimize();

		for (int sol = 0; sol < model->get(GRB_IntAttr_SolCount); ++sol) {
			model->set(GRB_IntParam_SolutionNumber, sol);
			if (model->get(GRB_DoubleAttr_PoolObjVal) >= 1e-6) {
				std::cout << "Viol " << model->get(GRB_DoubleAttr_PoolObjVal) << std::endl;
				std::unordered_set<int> S;
				for (int i : customers) {
					if (xi[i].get(GRB_DoubleAttr_Xn) >= 1 - 1e-6) {
						S.insert(i);
					}
				}
				violated_inequalities.emplace_back(BBV_RC_Inequality(S, S.size() - alpha.get(GRB_DoubleAttr_Xn), h));
				std::cout << "Alpha " << alpha.get(GRB_DoubleAttr_Xn) << std::endl;
			}
		}
	}
}