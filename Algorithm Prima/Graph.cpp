#include "Graph.h"
#include<random>
#include <limits>
#include<iostream>
#include<list>
/// Create random graph with random generation seed
Graph Graph::createRandomGraph(int n, double minWeight, double maxWeight, double edgesPercent) {
    std::random_device rd;  
    return createRandomGraph(n, minWeight, maxWeight, rd(), edgesPercent);
}
/// Create graph whith generation seed
Graph Graph::createRandomGraph(int n, double minWeight, double maxWeight, unsigned seed, double edgesPercent) {
    if (n < 1) // Return default graph
        return Graph();
    if (minWeight > maxWeight) // swap diapasone
        std::swap(minWeight, maxWeight);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(minWeight, maxWeight);
    Graph rdGraph;
    rdGraph.matrix.resize(n);
    for (auto& row : rdGraph.matrix)
        row.resize(n, std::numeric_limits<double>::infinity());
    int edgesToCreate = static_cast<int>((n * (n - 1) / 2.0) * (edgesPercent / 100.0));
    std::vector<std::pair<int, int>> allEdges;
    for (int i = 0; i < n; i++) {
        rdGraph.matrix[i][i] = 0;
        rdGraph.renamingList.push_back(i);
        for (int j = i + 1; j < n; j++) {
            allEdges.push_back({ i, j });
        }
    }
    std::shuffle(allEdges.begin(), allEdges.end(), gen);
    for (int i = 0; i < edgesToCreate; i++) {
        double randNum = dis(gen);
        int u = allEdges[i].first, v = allEdges[i].second;
        rdGraph.matrix[u][v] = randNum;
        rdGraph.matrix[v][u] = randNum; // Symmetry about the main diagonal
    }
    return rdGraph;
}
void Graph::print() {
    for (int i = 0; i < matrix.size();i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout <<i<<"."<<j<<": " << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << "vertexes homomorphism: " << std::endl;
    for (int i = 0; i < matrix.size(); i++) {
        std::cout << i<<"->"<< renamingList[i]<<"; " << std::endl;
    }
}
Graph::Graph() {}
Graph::Graph(int n, double minWeight, double maxWeight,double edgesPercent) {
    Graph g = Graph::createRandomGraph(n, minWeight, maxWeight, edgesPercent);
    this->matrix = g.matrix;
    this->renamingList = g.renamingList;
}
Graph::Graph(const std::vector<Graph>& connectComponentsList) {

}

std::vector<Graph> Graph::findOstForest() {
    return std::vector<Graph>();
}
void Graph::setRenaimingList(const std::vector<int>& list) {
    renamingList = list;
}
std::vector<int> Graph::getRenaimingList() {
    return renamingList;
}
bool Graph::isRenamingListContainsNomer(int nomer) {
    for (auto& nom : renamingList) {
        if (nom = nomer)
            return true;
    }
    return false;
}std::vector<int> Graph::connectionsForVertexes(int nomer, std::vector<bool>& checked) {
    std::vector<int> connections;
    std::list<int> queue;
    queue.push_back(nomer);
    while (queue.size() > 0) {
        connections.push_back(queue.front());
        checked[queue.front()] = true;
        for (int i = 0; i < matrix.size(); i++) {
            if (matrix[queue.front()][i] != std::numeric_limits<double>::infinity() && !checked[i]) {
                queue.push_back(i);
                checked[i] = true;
            }
        }
        queue.pop_front();
    }
    return connections;
}
Graph Graph::findConnectionComponentForVertex(int nomer, std::vector<bool>& checked) {
    Graph connectionComponent;
    std::vector<int> connectionsForNomer = connectionsForVertexes(nomer, checked);
    size_t n = connectionsForNomer.size();
    connectionComponent.matrix.resize(n);
    for (auto& row : connectionComponent.matrix) {
        row.resize(n);
    }
    connectionComponent.renamingList = connectionsForNomer;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                connectionComponent.matrix[i][j] = 0;
            }
            else {
                connectionComponent.matrix[i][j] = matrix[connectionsForNomer[i]][connectionsForNomer[j]];
            }
        }
    }
    return connectionComponent;
}
std::vector<Graph> Graph::findConnectComponents() {
    std::vector<Graph> ConnectionComponents;
    std::vector<bool> checkedVertexes(matrix.size());
    for (int i = 0; i < matrix.size(); i++) {
        checkedVertexes[i] = false;
    }
    for (int i = 0; i < matrix.size(); i++) { // Checking all vertex
        if (!checkedVertexes[i]) { // If not connected with previous, find new component
            ConnectionComponents.push_back(findConnectionComponentForVertex(i, checkedVertexes));
        }
    }
    return ConnectionComponents;
}