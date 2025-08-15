#pragma once

#include <algorithm>
#include <coroutine>
#include <sstream>
#include <set>

struct DistanceLocationPair {
	int distance;
	int location;

	DistanceLocationPair(int d, int i) : distance(d), location(i) {}
};

std::string itos(int i) { std::stringstream s; s << i; return s.str(); }

std::set<std::pair<int, int>> edges_from_to(const std::set<std::pair<int, int>>& edges,
	const std::set<int>& set1,
	const std::set<int>& set2) {
	std::set<std::pair<int, int>> edges_set;
	for (const auto& edge : edges) {
		int i = edge.first;
		int j = edge.second;
		if (set1.find(i) != set1.end() && set2.find(j) != set2.end()) {
			edges_set.insert(edge);
		}
		else if (set2.find(i) != set2.end() && set1.find(j) != set1.end()) {
			edges_set.insert(edge);
		}
	}
	return edges_set;
}


std::set<std::pair<int, int>> edges_in(const std::set<std::pair<int, int>>& edges, const std::set<int> set) {
	std::set<std::pair<int, int>> edges_set;
	for (const auto& edge : edges) {
		int i = edge.first;
		int j = edge.second;
		if (set.find(i) != set.end() && set.find(j) != set.end()) {
			edges_set.insert(edge);
		}
	}
	return edges_set;
}
std::set<std::pair<int,int>> delta_edges(const std::set<std::pair<int, int>>& edges, const std::set<int> set, const std::set<int> vertices) {
	std::set<int> set_bar;
	std::set_difference(vertices.begin(), vertices.end(), set.begin(), set.end(), std::inserter(set_bar, set_bar.begin()));
	return edges_from_to(edges, set, set_bar);
}
