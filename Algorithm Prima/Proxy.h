#pragma once
#include<string>
#include "Graph.h"
#include<fstream>
/// Logger for Graph using Proxy pattern
class GraphProxyLogger : public GraphInterface {
	std::string log;
	GraphInterface& graph;
	void printInLog(std::string mes) const{
		std::fstream stream;
		stream.open(log,std::ios::app);
		stream << mes << '\n';
		stream.close();
	};
public:
	///  take graph and log name
	GraphProxyLogger(GraphInterface& graph, std::string log):graph(graph) {
		this->log = log;
	};
	std::vector<int> getHomomorphism() const override {
		printInLog("gethomomorphism started");
		std::vector<int> h = graph.getHomomorphism();
		printInLog("gethomomorphism ended");
		return h;
	} ;
	std::vector<std::vector<double>> getAdjacencyMatrix() const override {
		printInLog("getAdjacencyMatrix started");
		std::vector<std::vector<double>> h = graph.getAdjacencyMatrix();
		printInLog("getAdjacencyMatrix ended");
		return h;
	};
	void printMatrix() const override{
		printInLog("printMatrix started");
		graph.printMatrix();
		printInLog("printMatrix ended");
	};
	void printEdges() const override {
		printInLog("printEdges started");
		graph.printEdges();
		printInLog("printEdges ended");
	};
	void addEdge(int vertex1, int vertex2, double weight) override{
		printInLog("addEdge started");
		graph.addEdge(vertex1, vertex2, weight);
		printInLog("addEdge ended");
	};
	void removeEdge(int vertex1, int vertex2) override {
		printInLog("removeEdge started");
		graph.removeEdge(vertex1, vertex2);
		printInLog("removeEdge ended");
	};
	std::vector<Graph> findConnectComponents() const override {
		printInLog("findConnectComponents started");
		std::vector<Graph>  h = graph.findConnectComponents();
		printInLog("findConnectComponents ended");
		return h;
	};
	std::vector<Graph> findMinSpanningForest() const override {
		printInLog("findMinSpanningForest started");
		std::vector<Graph>  h = graph.findMinSpanningForest();
		printInLog("findMinSpanningForest ended");
		return h;
	};
	std::vector<Graph> parallel_findMinSpanningForest() override {
		printInLog("parallel_findMinSpanningForest started");
		std::vector<Graph>  h = graph.parallel_findMinSpanningForest();
		printInLog("parallel_findMinSpanningForest ended");
		return h;
	};
	int size() const override {
		printInLog("size started");
		int size = graph.size();
		printInLog("size ended");
		return size;
	};
	double findGraphWeight() const override {
		printInLog("findGraphWeight started");
		double h = graph.findGraphWeight();
		printInLog("findGraphWeight ended");
		return h;
	};
};

