#pragma once
#include <vector>
#include <condition_variable>
#include <mutex>
#include<chrono>
class Graph
{
	std::mutex mtx;
	std::condition_variable cv;
	std::vector<std::vector<double>> matrix; //adjacency matrix
	std::vector<int> homomorphism; // homomorphism of vertex numbers
	void setRenaimingList(const std::vector<int>& list);
	bool isHomomorphismContainsNomer(const std::vector<int>& homomorphism, int nomer) const;
	std::vector<int> connectionsForVertexes(int nomer, std::vector<bool>& checked) const;
	Graph findConnectionComponentForVertex(int nomer, std::vector<bool>& checked) const;
	std::vector<int> fillPositionsHomomorphism(std::vector<int>& homomorphism, const std::vector<Graph>& components) const;
	void determineHomomrphizmForVertex(int curHomomorphism, int position, std::vector<int>& homomorphism) const;
	Graph findMinSpanningTree(const Graph& graph) const;
	void parallel_findConnectComponents(std::vector<Graph>& ConnectionComponents, bool& isReady);
public:
	std::vector<int> getRenaimingList() const;
	std::vector<std::vector<double>> getAdjacencyMatrix() const;
	static Graph createRandomGraph(int n, double minWeight, double maxWeight, double edgesPercent);
	static Graph createRandomGraph(int n, double minWeight, double maxWeight, unsigned seed, double edgesPercent);
	Graph();
	Graph(int n);
	Graph(const Graph& another);
	Graph(int n, double minWeight, double maxWeight, double edgesPercent=100);
	Graph(const std::vector<Graph>& componentsList);
	void printMatrix() const;
	void printEdges() const;
	void addEdge(int vertex1, int vertex2, double weight);
	void removeEdge(int vertex1, int vertex2);
	std::vector<Graph> findConnectComponents() const;
	std::vector<Graph> findMinSpanningForest() const;
	std::vector<Graph> parallel_findMinSpanningForest();
	int size() const;
	double findGraphWeight() const;
	//**********************************************
	static std::pair<std::chrono::duration<double>, std::chrono::duration<double>> sameGraphBenchMark(int n, double edgesPercent, unsigned amountOf‹easurements);
	static void differentGraphsBenchMark(int n, double minWeight, double maxWeight, double edgesPercent) ;
};
struct Edge {
	int v1 = 0;
	int v2 = 0;
	double weight = 0;
	Edge(int v1, int v2, double weight); 
};

