#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>

struct Coords {
    double x;
    double y;
    
    Coords() = default;
    Coords(double x_val, double y_val) : x(x_val), y(y_val) {};
};

struct InstanceData {
    // General Data
    int no_veh = -1;
    int no_cust = -1;
    int time_ub = -1;

    // Vehicle Data
    std::vector<Coords> veh_coords;
    std::vector<int> capacity;
    std::vector<int> veh_occ;
    std::vector<int> veh_arr_time;
    std::vector<std::vector<int>> veh_dist; // Distance to each customer

    // Customer Data
    std::vector<Coords> cust_coords;
    std::vector<int> demand;
    std::vector<int> cust_arr_time;
    std::vector<std::vector<int>> cust_dist; // Upper triangular matrix of distances - customer zero is the destination
};

// Reads a .vrp file
void load_instance_data(const std::string& filename, InstanceData& data) {
    std::string instance_path = "Instances/" + filename + ".txt";
    std::ifstream file(instance_path);
    if (!file.is_open()) {
        file.open(instance_path);
    }

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << '\n';
        exit(1);
    }

    // Load general instance data
    std::string line;
    std::getline(file, line);
    std::getline(file, line);
    size_t pos = line.find(":");
    data.no_veh = std::stoi(line.substr(pos + 1));

    std::getline(file, line);
    pos = line.find(":");
    data.no_cust = std::stoi(line.substr(pos + 1));

    std::getline(file, line);
    pos = line.find(":");
    data.time_ub = std::stoi(line.substr(pos + 1));

    std::getline(file, line);
    std::getline(file, line);

    // Read vehicle data
    data.veh_coords.resize(data.no_veh);
    data.capacity.resize(data.no_veh);
    data.veh_occ.resize(data.no_veh);
    data.veh_arr_time.resize(data.no_veh);
    for (int k = 0; k < data.no_veh; k++) {
        int no;
        double x, y;
        file >> no >> x >> y >> data.capacity[k] >> data.veh_occ[k] >> data.veh_arr_time[k];
        data.veh_coords[k] = Coords(x, y);
    }

    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);

    // Read Customer Data        
    data.cust_coords.resize(data.no_cust + 1);
    data.demand.resize(data.no_cust + 1);
    data.cust_arr_time.resize(data.no_cust + 1);
    for (int i = 0; i <= data.no_cust; i++) {
        int no;
        double x, y;
        file >> no >> x >> y >> data.demand[i] >> data.cust_arr_time[i];
        data.cust_coords[i] = Coords(x, y);
    }

    file.close();

    // Calculate Distances
    data.veh_dist.resize(data.no_veh);
    for (int k = 0; k < data.no_veh; k++) {
        data.veh_dist[k].resize(data.no_cust + 1);
    }
    data.cust_dist.resize(data.no_cust + 1);
    for (int i = 0; i <= data.no_cust; i++) {
        data.cust_dist[i].resize(data.no_cust + 1);
    }

    for (int i = 0; i <= data.no_cust; i++) {
        for (int k = 0; k < data.no_veh; k++) {
            double dx = data.veh_coords[k].x - data.cust_coords[i].x;
            double dy = data.veh_coords[k].y - data.cust_coords[i].y;
            data.veh_dist[k][i] = int(std::round(std::sqrt(dx * dx + dy * dy)));
        }

        for (int j = i + 1; j <= data.no_cust; j++) {
            double dx = data.cust_coords[j].x - data.cust_coords[i].x;
            double dy = data.cust_coords[j].y - data.cust_coords[i].y;
            int dist = int(std::round(std::sqrt(dx * dx + dy * dy)));
            data.cust_dist[i][j] = dist;
            data.cust_dist[j][i] = dist;
        }
    }
}


// Loads instance data
InstanceData load_instance(const std::string& filename) {
    InstanceData data;
    std::cout << filename << "\n";

    auto start = std::chrono::high_resolution_clock::now(); // Start measuring loading time
    load_instance_data(filename, data);

    auto stop = std::chrono::high_resolution_clock::now(); // Stop measuring loading time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start); // Calculate loading time in milliseconds
    std::cout << "Loading time: " << duration.count() << " milliseconds\n";
    return data;
}
