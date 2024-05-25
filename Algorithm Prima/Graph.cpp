#include "Graph.h"
#include<random>
#include <limits>
#include<iostream>
#include<list>
#include <set>
#include <thread>
Edge::Edge(int v1, int v2, double weight) {
    if (v1 < 0)
        v1 = 0;
    if (v2 < 0)
        v2 = 0;
    this->v1 = v1;
    this->v2 = v2;
    this->weight = weight;
    if (v1 == v2 && v1 != std::numeric_limits<int>::infinity())
        this->weight = 0;
}
/// Copy constructor
Graph::Graph(const Graph& another) {
    this->matrix = another.matrix;
    this->homomorphism = another.homomorphism;
}
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
        rdGraph.homomorphism.push_back(i);
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
void Graph::printMatrix() const {
    for (int i = 0; i < matrix.size();i++) {
        for (int j = 0; j < matrix.size(); j++) {
            std::cout <<i<<"."<<j<<": " << matrix[i][j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << "vertexes homomorphism: " << std::endl;
    for (int i = 0; i < matrix.size(); i++) {
        std::cout << i<<"->"<< homomorphism[i]<<"; ";
    }
}
void Graph::printEdges() const {
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = i+1; j < matrix.size(); j++) {
            if (matrix[i][j] != std::numeric_limits<double>::infinity()) {
                std::cout << i << " <-->" << j << ": " << matrix[i][j] << std::endl;
            }
        }
    }
    std::cout << std::endl << "vertexes homomorphism: " << std::endl;
    for (int i = 0; i < matrix.size(); i++) {
        std::cout << i << "->" << homomorphism[i] << "; ";
    }
}
Graph::Graph() {}
Graph::Graph(int n, double minWeight, double maxWeight,double edgesPercent) {
    Graph g = Graph::createRandomGraph(n, minWeight, maxWeight, edgesPercent);
    this->matrix = g.matrix;
    this->homomorphism = g.homomorphism;
}
bool Graph::isHomomorphismContainsNomer(const std::vector<int>& homomorphism, int nomer) const{
    for (int i = 0; i < homomorphism.size(); i++) {
        if (homomorphism[i] == nomer)
            return true;
    }
    return false;
}
void Graph::determineHomomrphizmForVertex(int curHomomorphism, int position, std::vector<int>& homomorphism) const{
    if (curHomomorphism >= homomorphism.size())
        curHomomorphism = 0;
    while(isHomomorphismContainsNomer(homomorphism, curHomomorphism)) {
        curHomomorphism++;
        if (curHomomorphism >= homomorphism.size())
            curHomomorphism = 0;
    }
    homomorphism[position] = curHomomorphism;
}
std::vector<int> Graph::fillPositionsHomomorphism(std::vector<int>& homomorphism, const std::vector<Graph>& components) const{
    int n = 0;
    for (auto& component : components) {
        n += component.matrix.size();
    }
    homomorphism.clear();
    for (int i = 0; i < n; i++) {
        homomorphism.push_back(-1);
    }
    std::vector<int> componentsPositions;
    for (int i = 0; i < components.size(); i++) {
        if (i == 0) {
            componentsPositions.push_back(0);
        }
        else {
            componentsPositions.push_back(componentsPositions[i - 1] + components[i - 1].size());
        }
        for (int j = 0; j < components[i].size(); j++) {
            determineHomomrphizmForVertex(components[i].homomorphism[j], componentsPositions[i] + j, homomorphism);
        }
    }
    return componentsPositions;
}
Graph::Graph(const std::vector<Graph>& componentsList) {
    std::vector<int> positionsHomomorhphism;
    std::vector<int> componentsPositions = fillPositionsHomomorphism(positionsHomomorhphism, componentsList);
    matrix.resize(positionsHomomorhphism.size());
    for (auto& row : matrix)
        row.resize(positionsHomomorhphism.size(), std::numeric_limits<double>::infinity());
    for (int k = 0; k < componentsList.size(); k++) {
        for (int i = 0; i < componentsList[k].size(); i++) {
            for (int j = 0; j < componentsList[k].size(); j++) {
                matrix[positionsHomomorhphism[componentsPositions[k] + i]][positionsHomomorhphism[componentsPositions[k] + j]] =
                    componentsList[k].matrix[i][j];
            }
        }
    }
    for (int i = 0; i < positionsHomomorhphism.size(); i++)
        homomorphism.push_back(i);
}
/// Find minimum spanning tree for graph, using alghorithm Prima
Graph Graph::findMinSpanningTree(const Graph& graph) const{
    int n = graph.size();
    //init tree
    Graph spanningTree;
    spanningTree.matrix.resize(n);
    for (auto& row : spanningTree.matrix)
        row.resize(n, std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; i++)
        spanningTree.homomorphism.push_back(i);
    std::set<int> visited, unvisited; //visited unvisited vertexes
    std::vector<Edge> treeEdges; // tree edges
    //visit first vertex
    for (int v = 0; v < n; v++)
        unvisited.insert(v); 
    visited.insert(0);
    unvisited.erase(0);
    for (int i = 0; i < n; i++)
        spanningTree.matrix[i][i] = 0;
    // Initialize Finish -> Start main loop
    while (!unvisited.empty()) {
        Edge edge(std::numeric_limits<int>::infinity(), std::numeric_limits<int>::infinity(),
            std::numeric_limits<double>::infinity());// start edge with infinity params

        for (const auto& from : visited) {//for all visited vertexes
            for (int to = 0; to < n; to++) {  //choose minimum edge
                if (from != to) {
                    bool isUnvisitedVertex = unvisited.find(to) != unvisited.end();
                    bool isEdgeExists = graph.matrix[from][to] == std::numeric_limits<double>::infinity() ? false : true;

                    if (isEdgeExists && isUnvisitedVertex) {
                        if (edge.weight > graph.matrix[from][to])
                            edge = { from, to, graph.matrix[from][to] };
                    }
                }
            }
        }

        if (edge.weight != std::numeric_limits<double>::infinity()) {
            treeEdges.emplace_back(edge);
            visited.insert(edge.v2);
            unvisited.erase(edge.v2);
        }
        else {
            break;
        }
    }

    // Add edges in tree
    for (const auto& edge : treeEdges) {
        spanningTree.matrix[edge.v1][edge.v2] = edge.weight;
        spanningTree.matrix[edge.v2][edge.v1] = edge.weight;
    }

    return spanningTree;
}
std::vector<Graph> Graph::findMinSpanningForest() const {
    std::vector<Graph> connectionComponents = findConnectComponents();
    std::vector<Graph> minSpanningForest;
    for (const auto& component : connectionComponents)
        minSpanningForest.push_back(findMinSpanningTree(component));
    return minSpanningForest;
}
void Graph::setRenaimingList(const std::vector<int>& list) {
    homomorphism = list;
}
std::vector<int> Graph::getRenaimingList() const {
    return homomorphism;
}
std::vector<int> Graph::connectionsForVertexes(int nomer, std::vector<bool>& checked) const {
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
Graph Graph::findConnectionComponentForVertex(int nomer, std::vector<bool>& checked) const {
    Graph connectionComponent;
    std::vector<int> connectionsForNomer = connectionsForVertexes(nomer, checked);
    size_t n = connectionsForNomer.size();
    connectionComponent.matrix.resize(n);
    for (auto& row : connectionComponent.matrix) {
        row.resize(n);
    }
    connectionComponent.homomorphism = connectionsForNomer;
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
std::vector<Graph> Graph::findConnectComponents() const {
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
void Graph::parallel_findConnectComponents(std::vector<Graph>& ConnectionComponents, bool& isReady)  {
    isReady = false;
    std::vector<bool> checkedVertexes(matrix.size());
    for (int i = 0; i < matrix.size(); i++) {
        checkedVertexes[i] = false;
    }
    for (int i = 0; i < matrix.size(); i++) { // Checking all vertex
        if (!checkedVertexes[i]) { // If not connected with previous, find new component
            std::lock_guard<std::mutex> lock(this->mtx);
            ConnectionComponents.push_back(findConnectionComponentForVertex(i, checkedVertexes));
            cv.notify_one();
        }
    }
    isReady = true;
}

std::vector<Graph> Graph::parallel_findMinSpanningForest() {
    std::vector<Graph> connectionComponents;
    std::vector<Graph> minSpanningForest;
    bool isReady = false;
    std::thread t1(&Graph::parallel_findConnectComponents, this, std::ref(connectionComponents), std::ref(isReady));
    std::thread t2([&]() {
        std::unique_lock<std::mutex> lock(mtx);
        while (!isReady) {
            cv.wait(lock);
            while (!connectionComponents.empty()) {
                minSpanningForest.push_back(findMinSpanningTree(connectionComponents.back()));
                connectionComponents.pop_back();
            }
        }
        });
    t1.join();
    t2.join();
    return minSpanningForest;
}
int Graph::size() const{
    return matrix.size();
}
double Graph::findGraphWeight() const {
    double weight = 0;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            if (matrix[i][j] != std::numeric_limits<double>::infinity())
                weight += matrix[i][j];
        }
    }
    return weight;  
}