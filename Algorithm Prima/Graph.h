#pragma once
#include <vector>
class Graph
{
	std::vector<std::vector<double>> matrix; //adjacency matrix
	std::vector<int> renamingList; // homomorphism of vertex numbers
	void setRenaimingList(const std::vector<int>& list);
	
	std::vector<int> connectionsForVertexes(int nomer, std::vector<bool>& checked);
	Graph findConnectionComponentForVertex(int nomer, std::vector<bool>& checked);
	bool isRenamingListContainsNomer(int nomer);
public:
	std::vector<int> getRenaimingList();
	static Graph createRandomGraph(int n, double minWeight, double maxWeight, double edgesPercent);
	static Graph createRandomGraph(int n, double minWeight, double maxWeight, unsigned seed, double edgesPercent);
	Graph();
	Graph(int n, double minWeight, double maxWeight, double edgesPercent=100);
	Graph(const std::vector<Graph>& connectComponentsList);
	void print();
	std::vector<Graph> findConnectComponents();
	std::vector<Graph> findOstForest();
};

