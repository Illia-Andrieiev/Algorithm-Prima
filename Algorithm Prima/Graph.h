#pragma once
#include <vector>
class Graph
{
	std::vector<std::vector<double>> matrix; //adjacency matrix
	std::vector<int> renamingList; // homomorphism of vertex numbers
	void setRenaimingList(const std::vector<int>& list);
	std::vector<int> getRenaimingList();
public:
	static Graph createRandomGraph(int n, double minWeight, double maxWeight);
	Graph();
	Graph(int n, double minWeight, double maxWeight);
	Graph(const std::vector<Graph>& connectComponentsList);
	std::vector<Graph> findConnectComponents();
	std::vector<Graph> findOstForest();

};

