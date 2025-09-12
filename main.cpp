
#include <fstream>
#include <iostream>
#include <iomanip>
#include <thread>

#include "event_based.h"
#include "route_based.h"
#include "route_and_event_based.h"
#include "fmrsp.h"

#include "gurobi_c++.h"

const int THREADS = 12;
const int TIME_LIMIT = 3600;
const bool RELAX = false;
const bool SEPARATE_RCI = true;
const bool SILENT = true;
const bool PREPROCESS = true;


std::vector<std::string> A_INSTANCES = {
	"A-n24-k7", "A-n24-k8", "A-n25-k8", "A-n27-k8", "A-n27-k9", "A-n28-k9", "A-n29-k9", "A-n33-k10",
	"A-n33-k11", "A-n34-k11", "A-n36-k11", "A-n39-k13", "A-n40-k13", "A-n41-k13", "A-n45-k14", "A-n45-k15",
	"A-n46-k15", "A-n47-k15", "A-n48-k15", "A-n48-k16", "A-n51-k17", "A-n60-k19"
};
std::vector<std::string> B_INSTANCES = {
	"B-n23-k7", "B-n25-k8", "B-n26-k8", "B-n28-k9", "B-n29-k9", "B-n30-k10", "B-n32-k10", "B-n33-k10",
	"B-n33-k11", "B-n37-k12", "B-n38-k12", "B-n39-k12", "B-n42-k13", "B-n42-k14", "B-n47-k15", "B-n48-k15",
	"B-n49-k16", "B-n50-k16", "B-n51-k16", "B-n58-k19"
};
std::vector<std::string> E_INSTANCES = {
	"E-n16-k5", "E-n17-k5", "E-n22-k7", "E-n24-k8", "E-n38-k12", "E-n57-k18", "E-n75-k25"
};
std::vector<std::string> F_INSTANCES = { "F-n33-k11", "F-n54-k17", "F-n101-k33" };
std::vector<std::string> M_INSTANCES = { "M-n75-k25", "M-n90-k30", "M-n113-k37", "M-n150-k49" };
std::vector<std::string> P_INSTANCES = {
	"P-n12-k3", "P-n14-k4", "P-n15-k4", "P-n15-k5", "P-n16-k5", "P-n17-k5", "P-n30-k9", "P-n33-k11",
	"P-n37-k12", "P-n38-k12", "P-n41-k13", "P-n45-k14", "P-n48-k16", "P-n52-k17", "P-n57-k18", "P-n75-k25"
};
std::vector<std::string> X_SMALL_INSTANCES = {
	"X-n75-k25", "X-n79-k26", "X-n82-k27", "X-n86-k28", "X-n90-k29", "X-n93-k31", "X-n96-k32",
	"X-n100-k33", "X-n104-k34", "X-n107-k35", "X-n111-k36", "X-n114-k38", "X-n117-k39", "X-n121-k40",
	"X-n125-k41", "X-n129-k42", "X-n132-k43", "X-n135-k45", "X-n139-k46", "X-n142-k47", "X-n146-k48",
	"X-n150-k49", "X-n153-k50", "X-n156-k52", "X-n160-k53", "X-n164-k54", "X-n167-k55", "X-n171-k56",
	"X-n174-k58", "X-n177-k59", "X-n181-k60", "X-n185-k61", "X-n188-k62", "X-n192-k63", "X-n195-k65",
	"X-n199-k66"
};
std::vector<std::string> X_MEDIUM_INSTANCES = {
	"X-n202-k67", "X-n206-k68", "X-n210-k69", "X-n213-k70", "X-n216-k72", "X-n220-k73", "X-n223-k74",
	"X-n227-k75", "X-n231-k76", "X-n234-k78", "X-n237-k79", "X-n241-k80", "X-n245-k81", "X-n248-k82",
	"X-n252-k83", "X-n258-k85", "X-n263-k87", "X-n269-k89", "X-n275-k91", "X-n282-k93", "X-n288-k95",
	"X-n294-k98", "X-n300-k100", "X-n308-k102", "X-n315-k104", "X-n321-k107", "X-n329-k109", "X-n336-k112",
	"X-n344-k114", "X-n351-k117", "X-n360-k119", "X-n368-k122", "X-n376-k125", "X-n384-k128", "X-n393-k130"
};
std::vector<std::string> X_LARGE_INSTANCES = {
	"X-n402-k133", "X-n411-k136", "X-n420-k140", "X-n429-k143", "X-n439-k146", "X-n449-k149", "X-n459-k153",
	"X-n470-k156", "X-n480-k160", "X-n491-k163", "X-n502-k167", "X-n513-k171", "X-n525-k175", "X-n537-k178",
	"X-n549-k183", "X-n561-k187", "X-n574-k191", "X-n587-k195", "X-n600-k200", "X-n614-k204", "X-n627-k209",
	"X-n642-k213", "X-n657-k218", "X-n671-k223", "X-n687-k228", "X-n702-k233", "X-n717-k239", "X-n734-k244",
	"X-n750-k250"
};

std::vector<std::string> TEST_INSTANCES{
	"A-n51-k17", "A-n60-k19",
	"B-n50-k16", "B-n51-k16", "B-n58-k19",
	"E-n57-k18", "E-n75-k25",
	"F-n54-k17", "F-n101-k33",
	"M-n75-k25", "M-n90-k30", "M-n113-k37", "M-n150-k49",
	"P-n52-k17", "P-n57-k18", "P-n75-k25",
	"X-n75-k25", "X-n90-k29", "X-n100-k33", "X-n111-k36", "X-n121-k40", "X-n129-k42",  "X-n146-k48",
	"X-n153-k50", "X-n167-k55"
};

int main() {
	//{ "A-n60-k19", "B-n58-k19", "P-n75-k25", "X-n82-k27", "M-n90-k30", "F-n101-k33" }

	for (const std::string file : {"B-n50-k16"}) {
		auto start = std::chrono::high_resolution_clock::now();
		// Works only in visual studio
		//FILE* stream;
		//freopen_s(&stream, "Logs/REB_test.txt", "w", stdout);
		// Open the log file
		//std::ofstream logFile("Logs/2I_COCF_TTCF/LP_" + file + ".txt", std::ios::out);
		//if (!logFile) {
		//	std::cerr << "Failed to open log file." << std::endl;
		//	return 1;
		//}
		//// Redirect both std::cout and std::cerr
		//std::streambuf* originalCoutBuffer = std::cout.rdbuf();
		//std::streambuf* originalCerrBuffer = std::cerr.rdbuf();
		//std::cout.rdbuf(logFile.rdbuf());
		//std::cerr.rdbuf(logFile.rdbuf());
		//FMRSP inst = FMRSP(file, "2I", "COCF", "TTCF", THREADS, TIME_LIMIT);
		//inst.preprocess();
		//EB inst = EB(file, "TTCF", THREADS, TIME_LIMIT);
		REB inst = REB(file, THREADS, TIME_LIMIT);
		//inst.print_all_routes();

		if (RELAX) inst.solve_root_relaxation(SEPARATE_RCI, SILENT, PREPROCESS);
		else inst.solve_to_optimality(SEPARATE_RCI, SILENT, PREPROCESS);
		//int threads = std::thread::hardware_concurrency();
		//inst.model->getEnv().set(GRB_IntParam_Threads, 2);
		//std::cout << inst.model->get(GRB_IntParam_ThreadLimit) << std::endl;
		//inst.set_starting_edges();
		//inst.preprocess();
		//inst.create_model(relax, false);
		//inst.optimise_model();

		//inst.print_routes();

		//inst.gamma[1 * 0 + 5].set(GRB_DoubleAttr_LB, 1);
		// EB solution
		//for (int idx = 0; idx < inst.var_indices.size(); idx++) {
		//	if (inst.vars[idx].get(GRB_DoubleAttr_X) >= 0.5) {
		//		std::cout << inst.vertices[inst.edges[idx].first][0] << "->" << inst.vertices[inst.edges[idx].second][0] << std::endl;
		//		//std::cout << inst.vars[idx].get(GRB_DoubleAttr_X) << " " << inst.vars[idx].get(GRB_StringAttr_VarName) << std::endl;
		//	}
		//}
		// RB SOlution
		//for (const auto& idx : inst.var_indices_subset) {
		//	if (inst.vars[idx].get(GRB_DoubleAttr_X) <= 0.1) continue;
		//	const int k = inst.idx_to_route[idx].first;
		//	const int i = inst.idx_to_route[idx].second;
		//	std::cout << k << ": ";
		//	for (const int& j : inst.routes[k][i]) {
		//		std::cout << j << "->";
		//	}
		//	std::cout << "d" << std::endl;
		//}
		//REB Solution
		//for (const auto& idx : inst.var_indices_subset) {
		//	if (inst.vars[idx].get(GRB_DoubleAttr_X) <= 0.1) continue;
		//	int k = 0;
		//	while (idx >= inst.first_route_index[k + 1]) k++;
		//	std::cout << k << ": ";
		//	for (const int& j : inst.routes[idx]) {
		//		std::cout << inst.veh_vertices[k][j][0] << "->";
		//	}
		//	std::cout << "d" << std::endl;
		//}
		//for (int idx : inst.flow_indices) {
		//	if (inst.flow[idx].get(GRB_DoubleAttr_X) >= 0.5) std::cout << inst.flow[idx].get(GRB_DoubleAttr_X) << " " << inst.flow[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.rho_veh_indices) {
		//	std::cout << inst.rho_veh[idx].get(GRB_DoubleAttr_X) << " " << inst.rho_veh[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.rho_indices) {
		//	std::cout << inst.rho[idx].get(GRB_DoubleAttr_X) << " " << inst.rho[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.chi_indices) {
		//	if (inst.chi[idx].get(GRB_DoubleAttr_X) > 0) std::cout << inst.chi[idx].get(GRB_DoubleAttr_X) << " " << inst.chi[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.tau_indices) {
		//	if (inst.tau[idx].get(GRB_DoubleAttr_X) > 0) std::cout << inst.tau[idx].get(GRB_DoubleAttr_X) << " " << inst.tau[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.psi_indices) {
		//	if (inst.psi[idx].get(GRB_DoubleAttr_X) > 0) std::cout << inst.psi[idx].get(GRB_DoubleAttr_X) << " " << inst.psi[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.gamma_indices) {
		//	if (inst.gamma[idx].get(GRB_DoubleAttr_X) > 0) std::cout << inst.gamma[idx].get(GRB_DoubleAttr_X) << " " << inst.gamma[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.sigma_indices) {
		//	if (inst.sigma[idx].get(GRB_DoubleAttr_X) > 0) std::cout << inst.sigma[idx].get(GRB_DoubleAttr_X) << " " << inst.sigma[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.omega_indices) {
		//	if (inst.omega[idx].get(GRB_DoubleAttr_X) > 0) std::cout << inst.omega[idx].get(GRB_DoubleAttr_X) << " " << inst.omega[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (const int i : inst.dest_and_cust) {
		//	if (i > 0) continue;
		//	for (const int j : inst.dest_and_cust) {
		//		if (i >= j) continue;
		//		std::cout << i << "->" << j << ": " << inst.cust_dist[i][j] << std::endl;
		//	}
		//}
		//for (const int k : inst.vehicles) {
		//	for (const int j : inst.dest_and_cust) {
		//		std::cout << k << "->" << j << ": " << inst.veh_dist[k][j] << std::endl;
		//	}
		//}
		//GRBVar* vars = inst.model->getVars();
		//for (int idx = 0; idx < inst.model->get(GRB_IntAttr_NumVars); idx++) {
		//	 if (vars[idx].get(GRB_DoubleAttr_X) > 0) std::cout << vars[idx].get(GRB_DoubleAttr_X) << " " << vars[idx].get(GRB_StringAttr_VarName) << std::endl;
		//}
		//for (int idx : inst.var_indices_subset) {
		//	if (idx >= inst.no_dir_veh_vars) break;
		//	int k = floor(idx / inst.no_cust_d);
		//	int i = idx - k * inst.no_cust_d;
		//	if (inst.xi[i][k].get(GRB_DoubleAttr_X) >= 0.5) std::cout << inst.xi[i][k].get(GRB_DoubleAttr_X) << " " << inst.xi[i][k].get(GRB_StringAttr_VarName) << std::endl;
		//}
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(stop - start);
		// Prints
		std::cout << "---------------------------------- Results ----------------------------------" << std::endl;
		std::cout << "INSTANCE: " << file << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "LP Objective: " << inst.initial_lp << std::endl;
		if (SEPARATE_RCI) std::cout << std::fixed << std::setprecision(2) << "LP+RCI Objective: " << inst.strengthened_lp << std::endl;
		if (!RELAX) std::cout << std::fixed << std::setprecision(2) << "MIP Objective: " << inst.mip << std::endl;
		if (!RELAX) std::cout << std::fixed << std::setprecision(2) << "MIP Bound: " << inst.model->get(GRB_DoubleAttr_ObjBound) << std::endl;
		if (!RELAX) std::cout << std::fixed << std::setprecision(2) << "Gap: " << inst.model->get(GRB_DoubleAttr_MIPGap) * 100 << "%" << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "Total time: " << duration.count() << " seconds\n";
		std::cout << std::fixed << std::setprecision(2) << "	Graph Construction + Preprocessing Time: " << inst.graph_construction_time << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "	Model Construction Time: " << inst.model_construction_time << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "	LP Solve Time: " << inst.lp_solve_time << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "	LP+RCC Solve Time: " << inst.lp_rcc_solve_time << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "	MIP Solve Time: " << inst.mip_solve_time << std::endl;
		std::cout << std::fixed << std::setprecision(2) << "	Total LP Solve Time: " << inst.total_lp_solve_time << std::endl;
		std::cout << "# Variables: " << inst.model->get(GRB_IntAttr_NumVars) << std::endl;
		std::cout << "# Constraints: " << inst.model->get(GRB_IntAttr_NumConstrs) << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;

		//std::cout.rdbuf(originalCoutBuffer);
		//std::cerr.rdbuf(originalCerrBuffer);
		//logFile.close();
	}
}