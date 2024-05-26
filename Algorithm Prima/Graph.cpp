#include "Graph.h"
#include<tbb/tbb.h>
#include <limits>
#include<iostream>
#include <set>
#include <thread>
#include<fstream>
///Constructor
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
Graph& Graph::operator =(const Graph& another) {
    matrix = another.matrix;
    homomorphism = another.homomorphism;
    return *this;
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
    for (auto& row : rdGraph.matrix) // init matrix as infinity
        row.resize(n, std::numeric_limits<double>::infinity());
    int edgesToCreate = static_cast<int>((n * (n - 1) / 2.0) * (edgesPercent / 100.0)); // how many edges create
    std::vector<std::pair<int, int>> allEdges;
    for (int i = 0; i < n; i++) {
        rdGraph.matrix[i][i] = 0; // elems on diagonal = 0
        rdGraph.homomorphism.push_back(i);
        for (int j = i + 1; j < n; j++) {
            allEdges.push_back({ i, j }); 
        }
    }
    std::shuffle(allEdges.begin(), allEdges.end(), gen); // shuffle edges
    for (int i = 0; i < edgesToCreate; i++) {
        double randNum = dis(gen);
        int u = allEdges[i].first, v = allEdges[i].second;
        rdGraph.matrix[u][v] = randNum;
        rdGraph.matrix[v][u] = randNum; // Symmetry about the main diagonal
    }
    return rdGraph;
}
/// Print adjacency matrix 
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
/// Print all graph edges 
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
/// Constructor
/*!
* \param[in] n number of vertexes in the graph
* \param[in] minWeight minimum possible edge weight.
* \param[in] maxWeight maximum possible edge weight.
* \param[in] edgesPercent percent of all edges that will be created
*/
Graph::Graph(int n, double minWeight, double maxWeight,double edgesPercent) {
    Graph g = Graph::createRandomGraph(n, minWeight, maxWeight, edgesPercent);
    this->matrix = g.matrix;
    this->homomorphism = g.homomorphism;
}
/// Constructor
/*!
* \param[in] componentsList list of components, that will merge into one
*/
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
/// Is homomorphism contains nomer
bool Graph::isHomomorphismContainsNomer(const std::vector<int>& homomorphism, int nomer) const{
    for (int i = 0; i < homomorphism.size(); i++) {
        if (homomorphism[i] == nomer)
            return true;
    }
    return false;
}
/// Determine homomorphism, that will be set for position vertex. curHomomorphism is default homomorphism 
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
/// determine general homomorphism for all components. Return vectors of positions components in homomorphism
std::vector<int> Graph::fillPositionsHomomorphism(std::vector<int>& homomorphism, const std::vector<Graph>& components) const{
    int n = 0;
    for (auto& component : components) {
        n += (int)component.matrix.size();
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

/// Find minimum spanning tree for graph, using alghorithm Prima. Works with O(n^2). 
Graph Graph::findMinSpanningTree(const Graph& graph) const{
    int n = graph.size();
    //init spanning tree
    Graph spanningTree;
    spanningTree.matrix.resize(n);
    for (auto& row : spanningTree.matrix)
        row.resize(n, std::numeric_limits<double>::infinity());
    spanningTree.homomorphism = graph.homomorphism;
    //                 Old algorithm 
    //std::set<int> visited, unvisited; //visited unvisited vertexes
    //std::vector<Edge> treeEdges; // tree edges
    ////visit first vertex
    //for (int v = 0; v < n; v++)
    //    unvisited.insert(v); 
    //visited.insert(0);
    //unvisited.erase(0);
    //for (int i = 0; i < n; i++)
    //    spanningTree.matrix[i][i] = 0;
    //// Initialize Finish -> Start main loop
    //while (!unvisited.empty()) {
    //    Edge edge(std::numeric_limits<int>::infinity(), std::numeric_limits<int>::infinity(),
    //        std::numeric_limits<double>::infinity());// start edge with infinity params

    //    for (const auto& from : visited) {//for all visited vertexes
    //        for (int to = 0; to < n; to++) {  //choose minimum edge
    //            if (from != to) {
    //                bool isUnvisitedVertex = unvisited.find(to) != unvisited.end();
    //                bool isEdgeExists = graph.matrix[from][to] == std::numeric_limits<double>::infinity() ? false : true;

    //                if (isEdgeExists && isUnvisitedVertex) {
    //                    if (edge.weight > graph.matrix[from][to])
    //                        edge = { from, to, graph.matrix[from][to] };
    //                }
    //            }
    //        }
    //    }
    //    if (edge.weight != std::numeric_limits<double>::infinity()) {
    //        treeEdges.emplace_back(edge);
    //        visited.insert(edge.v2);
    //        unvisited.erase(edge.v2);
    //    }
    //    else {
    //        break;
    //    }
    //}
    //// Add edges in tree
    //for (const auto& edge : treeEdges) {
    //    spanningTree.matrix[edge.v1][edge.v2] = edge.weight;
    //    spanningTree.matrix[edge.v2][edge.v1] = edge.weight;
    //}
    std::vector<bool> used(n); // vertexes that have been included in the minimum spanning tree.
    std::vector<double> minE(n, std::numeric_limits<double>::infinity()), selE(n, -1);  // minE[i] stores the weight of the edge connecting vertex i to the minimum spanning tree. selE[i] stores the selected edge that connects vertex i to the minimum spanning tree.
    minE[0] = 0; // The weight of the edge connecting the first vertex to the minimum spanning tree is set to 0.

    for (int i = 0; i < n; ++i) { 
        spanningTree.matrix[i][i] = 0; // The weight of the edge connecting a vertex to itself is always 0.
        int v = -1; // v will hold the index of the vertex to be added to the minimum spanning tree.
        for (int j = 0; j < n; ++j) // This loop finds the vertex with the smallest edge not yet included in the minimum spanning tree.
            if (!used[j] && (v == -1 || minE[j] < minE[v]))
                v = j;
        if (minE[v] == std::numeric_limits<double>::infinity()) { // If the smallest edge has infinite weight, then the graph is not connected, and the program exits.
            exit(0);
        }
        used[v] = true; // The vertex v is marked as used.
        if (selE[v] != -1) { // If there is an edge connecting vertex v to the minimum spanning tree, it is added to the spanning tree.
            spanningTree.addEdge(v, selE[v], graph.matrix[v][selE[v]]);
        }
        for (int to = 0; to < n; ++to) // This loop updates the weights of the edges connecting the vertices to the minimum spanning tree.
            if (graph.matrix[v][to] < minE[to]) {
                minE[to] = graph.matrix[v][to];
                selE[to] = v;
            }
    }
    return spanningTree;
}
/// This function finds the minimum spanning forest of a graph.
/*!
 * The function first finds the connected components of the graph. For each connected component, it finds the minimum spanning tree.
 * All these minimum spanning trees together form the minimum spanning forest of the graph.
 *
 * @return std::vector<Graph> - A vector of Graph objects, where each Graph object represents a minimum spanning tree of a connected component of the graph.
 */
std::vector<Graph> Graph::findMinSpanningForest() const {
    std::vector<Graph> connectionComponents = findConnectComponents(); ///< A vector of Graph objects, where each Graph object represents a connected component of the graph.
    std::vector<Graph> minSpanningForest; ///< A vector of Graph objects to store the minimum spanning trees of the connected components.
    for (const auto& component : connectionComponents) ///< For each connected component...
        minSpanningForest.push_back(findMinSpanningTree(component)); ///< ...find its minimum spanning tree and add it to the minimum spanning forest.
    return minSpanningForest; ///< Return the minimum spanning forest.
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
void Graph::findMinSpanningTreeThread(std::list<Graph>& componentsQueue, std::vector<Graph>& minSpanningForest, bool& isReady) {
    while (true) {

        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&] { return !componentsQueue.empty() || isReady; });
        if (!componentsQueue.empty()) {
            Graph component = componentsQueue.front();
            componentsQueue.pop_front();
            minSpanningForest.push_back(findMinSpanningTree(component));
        }
        if (isReady && componentsQueue.empty()) {
            return;
        }
    }
}
void Graph::parallel_findConnectComponents(std::list<Graph>& ConnectionComponents, bool& isReady) {
    isReady = false;
    std::vector<bool> checkedVertexes(matrix.size(), false);
    for (int i = 0; i < matrix.size(); i++) { // Checking all vertex
        if (!checkedVertexes[i]) { // If not connected with previous, find new component
            Graph component = findConnectionComponentForVertex(i, checkedVertexes);
            {
                std::lock_guard<std::mutex> lock(mtx);
                ConnectionComponents.push_back(component);
            }
            cv.notify_one();
        }
    }
    isReady = true;
    cv.notify_all(); // Notify all waiting threads that isReady is now true
}
std::vector<Graph> Graph::parallel_findMinSpanningForest() {
    std::list<Graph> componentsQueue; 
    std::vector<Graph> minSpanningForest;
    bool isReady = false;
    tbb::task_group g;

    g.run([&] { this->parallel_findConnectComponents(componentsQueue,isReady); });
    int threadsAmount = 2;
    for (int i = 0; i < threadsAmount; ++i) {
        g.run([&]() {findMinSpanningTreeThread(componentsQueue,minSpanningForest,isReady); });
    }

    g.wait();
    return minSpanningForest;
}
int Graph::size() const{
    return (int)matrix.size();
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
void Graph::addEdge(int vertex1, int vertex2, double weight) {
    if (vertex1 < 0 || vertex2 < 0 || vertex1 >= matrix.size()||vertex2 >=  matrix.size())
        return;
    if (vertex1 == vertex2)
        weight = 0;
    matrix[vertex2][vertex1] = weight;
    matrix[vertex1][vertex2] = weight;
}
void Graph::removeEdge(int vertex1, int vertex2) {
    if (vertex1 < 0 || vertex2 < 0 || vertex1 >= matrix.size() || vertex2 >= matrix.size())
        return;
    matrix[vertex2][vertex1] = std::numeric_limits<double>::infinity();
    matrix[vertex1][vertex2] = std::numeric_limits<double>::infinity();
}
std::vector<std::vector<double>> Graph::getAdjacencyMatrix() const {
    return matrix;
}
/*
    Benchmark
*/
std::pair<std::chrono::duration<double>, std::chrono::duration<double>> Graph::sameGraphBenchMark(int n,
    double edgesPercent, unsigned amountOfMeasurements) {
    std::random_device rd;
    return sameGraphBenchMark(n, edgesPercent, amountOfMeasurements, rd());
}
std::pair<std::chrono::duration<double>, std::chrono::duration<double>> Graph::sameGraphBenchMark(int n, double edgesPercent, unsigned amountOfMeasurements, unsigned seed) {
    int printEvery = 50;
    Graph graph = Graph::createRandomGraph(n, -1000, 1000,seed ,edgesPercent);
    int amount = amountOfMeasurements;
    // start counting
    auto startDefault = std::chrono::high_resolution_clock::now();
    while (amount > 0) {
        amount--;
        if (amount % printEvery == 0)
            std::cout << "default " << amount << std::endl;
        graph.findMinSpanningForest();
    }
    // stop counting
    auto endDefault = std::chrono::high_resolution_clock::now();
    // find time
    std::chrono::duration<double> resTimeDefault = endDefault - startDefault;
    amount = amountOfMeasurements;
    // start counting
    auto startParallel = std::chrono::high_resolution_clock::now();
    while (amount > 0) {
        amount--;
        if (amount % printEvery == 0)
            std::cout << "parallel " << amount << std::endl;
        graph.parallel_findMinSpanningForest();
    }
    // stop counting
    auto endParallel = std::chrono::high_resolution_clock::now();
    // find time
    std::chrono::duration<double> resTimeParallel = endParallel - startParallel;
    amount = amountOfMeasurements;
    // start counting
    auto startParallel_W = std::chrono::high_resolution_clock::now();
    while (amount > 0) {
        amount--;
        if (amount % printEvery == 0)
            std::cout << "parallel_without_Pool " << amount << std::endl;
        graph.p_f();
    }
    // stop counting
    auto endParallel_W = std::chrono::high_resolution_clock::now();
    // find time
    std::chrono::duration<double> resTimeParallel_W = endParallel_W - startParallel_W;
    // log
    std::cout << std::endl << "Graph vertexes: " << n << " edges percent: " << edgesPercent << " amount of measurements: " << amountOfMeasurements <<
        " default algorithm time: " << resTimeDefault.count() << " s" << " parallel algorithm time: " << resTimeParallel.count() << " s"
        << " parallel_without_pool algorithm time: " << resTimeParallel_W.count() << " s" << std::endl;
    std::fstream log;
    log.open("log_benchmark.txt", std::ios::app | std::ios::in | std::ios::out);
    log << "Graph vertexes: " << n << " edges percent: " << edgesPercent << " amount of measurements: " << amountOfMeasurements <<
        " default algorithm time: " << resTimeDefault.count() << " s" << " parallel algorithm time: " << resTimeParallel.count() << " s"
        << " parallel_without_pool algorithm time: " << resTimeParallel_W.count() << " s" << std::endl;
    log.close();
    return std::pair<std::chrono::duration<double>, std::chrono::duration<double>>(resTimeDefault, resTimeParallel);
}
std::pair<std::chrono::duration<double>, std::chrono::duration<double>> Graph::differentGraphsBenchMark(int n, double edgesPercent, unsigned amountOfMeasurements) {
    int printEvery = 50;
    std::chrono::duration<double> resTimeDefault(0), resTimeParallel(0), resTimeParallel_W(0);
    Graph graph = Graph::createRandomGraph(n, -1000, 1000, edgesPercent);
    int amount = amountOfMeasurements;
    while (amount > 0) {
        amount--;
        if (amount % printEvery == 0)
            std::cout << "default " << amount << std::endl;
        // start counting
        auto startDefault = std::chrono::high_resolution_clock::now();
        graph.findMinSpanningForest();
        // stop counting
        auto endDefault = std::chrono::high_resolution_clock::now();
        // find time
        resTimeDefault += (endDefault - startDefault);
        graph = Graph::createRandomGraph(n, -1000, 1000, edgesPercent);
    }
    amount = amountOfMeasurements;
    while (amount > 0) {
        amount--;
        if (amount % printEvery == 0)
            std::cout << "parallel " << amount << std::endl;
        // start counting
        auto startParallel = std::chrono::high_resolution_clock::now();
        graph.parallel_findMinSpanningForest();
        // stop counting
        auto endParallel = std::chrono::high_resolution_clock::now();
        // find time
        resTimeParallel += (endParallel - startParallel);
        graph = Graph::createRandomGraph(n, -1000, 1000, edgesPercent);
    }
    amount = amountOfMeasurements;
    while (amount > 0) {
        amount--;
        if (amount % printEvery == 0)
            std::cout << "parallel_without_Pool " << amount << std::endl;
        // start counting
        auto startParallel_W = std::chrono::high_resolution_clock::now();
        graph.p_f();
        // stop counting
        auto endParallel_W = std::chrono::high_resolution_clock::now();
        // find time
        resTimeParallel_W += (endParallel_W - startParallel_W);
        graph = Graph::createRandomGraph(n, -1000, 1000, edgesPercent);
    }
    // log
    std::cout << std::endl << "Graph vertexes: " << n << " edges percent: " << edgesPercent << " amount of measurements: " << amountOfMeasurements <<
        " default algorithm time: " << resTimeDefault.count() << " s" << " parallel algorithm time: " << resTimeParallel.count() << " s"
        << " parallel_without_pool algorithm time: " << resTimeParallel_W.count() << " s" << std::endl;
    std::fstream log;
    log.open("log_benchmark.txt", std::ios::app | std::ios::in | std::ios::out);
    log << "Graph vertexes: " << n << " edges percent: " << edgesPercent << " amount of measurements: " << amountOfMeasurements <<
        " default algorithm time: " << resTimeDefault.count() << " s" << " parallel algorithm time: " << resTimeParallel.count() << " s"
        << " parallel_without_pool algorithm time: " << resTimeParallel_W.count() << " s" << std::endl;
    log.close();
    return std::pair<std::chrono::duration<double>, std::chrono::duration<double>>(resTimeDefault, resTimeParallel);
}
/********************/
void Graph::p_c(std::vector<Graph>& ConnectionComponents, bool& isReady) {
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
    cv.notify_all();
}

std::vector<Graph> Graph::p_f() {
    std::vector<Graph> connectionComponents;
    std::vector<Graph> minSpanningForest;
    bool isReady = false;
    std::thread t1(&Graph::p_c, this, std::ref(connectionComponents), std::ref(isReady));
    std::unique_lock<std::mutex> lock(mtx);
    while (!isReady) {
        cv.wait(lock);
        while (!connectionComponents.empty()) {
            minSpanningForest.push_back(findMinSpanningTree(connectionComponents.back()));
            connectionComponents.pop_back();
        }
    }
    t1.join();
    return minSpanningForest;
}