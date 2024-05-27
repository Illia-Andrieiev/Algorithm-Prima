#include "Graph.h"
#include<tbb/tbb.h>
#include <limits>
#include<iostream>
#include <set>
#include <thread>
#include<fstream>
#include"libs/ut.hpp"
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
/// Constructor from matrix
Graph::Graph(const std::vector<std::vector<double>>& matrix) {
    int n = (int)matrix.size();
    this->matrix.resize(n);
    for (auto& row : this->matrix)
        row.resize(n);
    for (int i = 0; i < n; i++) {
        this->matrix[i][i] = 0;
        for (int j = i + 1; j < n; j++) {
            this->matrix[i][j] = matrix[i][j];
            this->matrix[j][i] = matrix[i][j];
        }
    }
}
Graph& Graph::operator =(const Graph& another) {
    matrix = another.matrix;
    homomorphism = another.homomorphism;
    return *this;
}
bool Graph::operator ==(const Graph& another) const{
    return matrix == another.matrix && homomorphism == another.homomorphism;
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

void Graph::setHomomorphism(const std::vector<int>& list) {
    homomorphism = list;
}
/// Get graphs homomorphism
std::vector<int> Graph::getHomomorphism() const {
    return homomorphism;
}
/// Find all vertexes, conected with 'nomer' vertex. in 'checked' vector set connected vertexes as true. 
std::vector<int> Graph::connectionsForVertexes(int nomer, std::vector<bool>& checked) const {
    std::vector<int> connections;
    std::list<int> queue; // queue for checking vertexes
    queue.push_back(nomer);
    while (queue.size() > 0) {
        // set connection with vertex
        connections.push_back(queue.front()); 
        checked[queue.front()] = true;
        for (int i = 0; i < matrix.size(); i++) {
            if (matrix[queue.front()][i] != std::numeric_limits<double>::infinity() && !checked[i]) { // if connected, add to queue
                queue.push_back(i);
                checked[i] = true;
            }
        }
        queue.pop_front(); // remove checked vertex
    }
    return connections;
}
/// This function finds the connected component of a graph for a given vertex. Works with O(n^2) 
/*!
 * The function takes as input a vertex and a reference to a vector of booleans that keeps track of which vertices have been checked.
 * It then finds all the vertices connected to the given vertex and creates a subgraph (connected component) with these vertices.
 * The function returns this subgraph.
 * @param nomer The vertex for which the connected component is to be found.
 * @param checked A reference to a vector of booleans that keeps track of which vertices have been checked.
 * @return Graph The connected component of the graph for the given vertex.
 */
Graph Graph::findConnectionComponentForVertex(int nomer, std::vector<bool>& checked) const {
    Graph connectionComponent; ///< The connected component of the graph for the given vertex.
    std::vector<int> connectionsForNomer = connectionsForVertexes(nomer, checked); ///< The vertexes connected to the given vertex.
    size_t n = connectionsForNomer.size(); ///< The number of vertexes connected to the given vertex.
    connectionComponent.matrix.resize(n); ///< Resize the adjacency matrix of the connected component to the number of connected vertices.
    for (auto& row : connectionComponent.matrix) { ///< For each row in the adjacency matrix...
        row.resize(n); ///< ...resize it to the number of connected vertices.
    }
    connectionComponent.homomorphism = connectionsForNomer; ///< The homomorphism of the connected component is the vertices connected to the given vertex.
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) { 
            if (i == j) { 
                connectionComponent.matrix[i][j] = 0; ///< the weight of the edge connecting the vertex to itself is 0.
            }
            else { 
                connectionComponent.matrix[i][j] = matrix[connectionsForNomer[i]][connectionsForNomer[j]]; ///< the weight of the edge connecting the two vertices is the same as in the original graph.
            }
        }
    }
    return connectionComponent; ///< Return the connected component.
}
///This function finds all the connected components of a graph. Works with O(n^2)
/*!
 * The function iterates over all the vertices of the graph. For each vertex, if it has not been checked (i.e., it is not part of a previously found connected component),
 * the function finds the connected component for that vertex and adds it to the list of connected components.
 *
 * @return std::vector<Graph> - A vector of Graph objects, where each Graph object represents a connected component of the graph.
 */
std::vector<Graph> Graph::findConnectComponents() const {
    std::vector<Graph> ConnectionComponents; ///< A vector of Graph objects to store the connected components of the graph.
    std::vector<bool> checkedVertexes(matrix.size()); ///< A vector of booleans to keep track of which vertices have been checked.
    for (int i = 0; i < matrix.size(); i++) { ///< Initialize all vertices as unchecked.
        checkedVertexes[i] = false;
    }
    for (int i = 0; i < matrix.size(); i++) {
        if (!checkedVertexes[i]) { ///< if the vertex has not been checked
            ConnectionComponents.push_back(findConnectionComponentForVertex(i, checkedVertexes)); ///< ...find the connected component for that vertex and add it to the list of connected components.
        }
    }
    return ConnectionComponents; ///< Return the list of connected components.
}
/// This function finds the minimum spanning tree of each connected component in a graph using multiple threads.
/*!
 * The function waits until there is a connected component in the queue or all connected components have been found.
 * If there is a connected component in the queue, it removes it from the queue, finds its minimum spanning tree, and adds it to the minimum spanning forest.
 * The function continues this process until all connected components have been found and the queue is empty.
 *
 * @param componentsQueue A reference to a queue of connected components.
 * @param minSpanningForest A reference to a vector of Graph objects to store the minimum spanning trees of the connected components.
 * @param isReady A reference to a boolean that indicates whether all connected components have been found.
 */
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
/// This function finds all the connected components of a graph using multiple threads.
/*!
 * The function iterates over all the vertices of the graph. For each vertex, if it has not been checked (i.e., it is not part of a previously found connected component),
 * the function finds the connected component for that vertex and adds it to the queue of connected components.
 * Once all connected components have been found, it notifies all waiting threads.
 *
 * @param ConnectionComponents A reference to a queue of connected components.
 * @param isReady A reference to a boolean that indicates whether all connected components have been found.
 */
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
/// This function finds the minimum spanning forest of a graph using multiple threads.
/*!
 * The function creates a task group and runs two tasks in parallel: one to find all the connected components of the graph, and the other to find the minimum spanning tree of each connected component.
 * Once all tasks have completed, the function returns the minimum spanning forest.
 *
 * @return std::vector<Graph> - A vector of Graph objects, where each Graph object represents a minimum spanning tree of a connected component of the graph.
 */
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
/// Return amount of graph`s vertexes
int Graph::size() const{
    return (int)matrix.size();
}
/// This function calculates the total weight of a graph.
/*!
 * The function iterates over all the edges in the adjacency matrix of the graph. If the weight of an edge is not infinity (i.e., the edge exists),
 * it adds the weight of the edge to the total weight of the graph. The function returns the total weight of the graph.
 *
 * @return double - The total weight of the graph.
 */
double Graph::findGraphWeight() const {
    double weight = 0; ///< The total weight of the graph.
    for (int i = 0; i < matrix.size(); i++) { ///< For each vertex in the graph...
        for (int j = i+1; j < matrix.size(); j++) { 
            if (matrix[i][j] != std::numeric_limits<double>::infinity()) ///< if the edge exists
                weight += matrix[i][j]; ///< åadd the weight of the edge to the total weight of the graph.
        }
    }
    return weight;
}
///  Add edge to graph
/*!
* if vertexes nomers out of range, method do not add edge. if vertex1 == vertex2 -> weight = 0. Also add edge for v2 -> v1.
* @param vertex1 nomer of first vertex in edge
* @param vertex2 nomer of second vertex in edge
* @param weight weight of edge
*/
void Graph::addEdge(int vertex1, int vertex2, double weight) {
    if (vertex1 < 0 || vertex2 < 0 || vertex1 >= matrix.size()||vertex2 >=  matrix.size())
        return;
    if (vertex1 == vertex2)
        weight = 0;
    matrix[vertex2][vertex1] = weight;
    matrix[vertex1][vertex2] = weight;
}
///  Remove edge from graph
/*!
* if vertexes nomers out of range, method do not remove edge. if vertex1 == vertex2 -> weight = 0. Also add edge for v2 -> v1.
* @param vertex1 nomer of first vertex in edge
* @param vertex2 nomer of second vertex in edge
*/
void Graph::removeEdge(int vertex1, int vertex2) {
    if (vertex1 < 0 || vertex2 < 0 || vertex1 >= matrix.size() || vertex2 >= matrix.size() || vertex1 == vertex2)
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

/// This function benchmarks the performance of the default and parallel implementations of the minimum spanning forest algorithm on the same graph.
/*!
 * The function creates a random graph and measures the time taken by the default and parallel implementations of the minimum spanning forest algorithm to process the graph.
 * The function repeats this process a specified number of times and returns the total time taken by each implementation.
 *
 * @param n The number of vertices in the graph.
 * @param edgesPercent The percentage of edges in the graph.
 * @param amountOfMeasurements The number of times the benchmarking process is repeated.
 * @return std::pair<std::chrono::duration<double>, std::chrono::duration<double>> - A pair of durations, where the first duration is the total time taken by the default implementation and the second duration is the total time taken by the parallel implementation.
 */
std::pair<std::chrono::duration<double>, std::chrono::duration<double>> Graph::sameGraphBenchMark(int n,
    double edgesPercent, unsigned amountOfMeasurements) {
    std::random_device rd;
    return sameGraphBenchMark(n, edgesPercent, amountOfMeasurements, rd());
}
/// This function benchmarks the performance of the default and parallel implementations of the minimum spanning forest algorithm on the same graph with a specified seed for the random number generator.
/*!
 * The function creates a random graph with a specified seed for the random number generator and measures the time taken by the default and parallel implementations of the minimum spanning forest algorithm to process the graph.
 * The function repeats this process a specified number of times and returns the total time taken by each implementation.
 *
 * @param n The number of vertices in the graph.
 * @param edgesPercent The percentage of edges in the graph.
 * @param amountOfMeasurements The number of times the benchmarking process is repeated.
 * @param seed The seed for the random number generator.
 * @return std::pair<std::chrono::duration<double>, std::chrono::duration<double>> - A pair of durations, where the first duration is the total time taken by the default implementation and the second duration is the total time taken by the parallel implementation.
 */
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
/// This function benchmarks the performance of the default and parallel implementations of the minimum spanning forest algorithm on different graphs.
/*!
 * The function creates a random graph and measures the time taken by the default and parallel implementations of the minimum spanning forest algorithm to process the graph.
 * The function repeats this process with a new random graph a specified number of times and returns the total time taken by each implementation.
 *
 * @param n The number of vertices in the graph.
 * @param edgesPercent The percentage of edges in the graph.
 * @param amountOfMeasurements The number of times the benchmarking process is repeated.
 * @return std::pair<std::chrono::duration<double>, std::chrono::duration<double>> - A pair of durations, where the first duration is the total time taken by the default implementation and the second duration is the total time taken by the parallel implementation.
 */
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


/// This function finds all the connected components of a graph.
/*!
 * The function iterates over all the vertices of the graph. For each vertex, if it has not been checked (i.e., it is not part of a previously found connected component),
 * the function finds the connected component for that vertex and adds it to the list of connected components.
 * Once all connected components have been found, it notifies all waiting threads.
 *
 * \param[out] ConnectionComponents A reference to a vector of Graph objects to store the connected components of the graph.
 * \param[in,out] isReady A reference to a boolean that indicates whether all connected components have been found.
 */
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
/// This function finds the minimum spanning forest of a graph using multiple threads.
/*!
 * The function creates a separate thread to find all the connected components of the graph. It then waits until a connected component is found,
 * finds its minimum spanning tree, and adds it to the minimum spanning forest.
 * The function continues this process until all connected components have been found.
 *
 * @return std::vector<Graph> - A vector of Graph objects, where each Graph object represents a minimum spanning tree of a connected component of the graph.
 */
std::vector<Graph> Graph::p_f() {
    std::vector<Graph> connectionComponents;
    std::vector<Graph> minSpanningForest;
    bool isReady = false;
    std::thread t1(&Graph::p_c, this, std::ref(connectionComponents), std::ref(isReady)); // find connection components
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
/*
                            Tests
*/

void Graph::test() {
    using namespace boost::ut;
    suite<"Find_Connect_Components_Test"> components_tests = [] {
        "test_findConnectComponents"_test = [] {
            // Arrange
            Graph g(5, -5, 5, 0);
            // Assuming you have a method to add edges to the graph
            g.addEdge(0, 1, 3);
            g.addEdge(1, 2, 5);
            g.addEdge(3, 4, 3);

            // Act
            auto components = g.findConnectComponents();

            // Assert
            expect(components.size() == 2 && components[0].size() == 3 && components[1].size() == 2) << fatal << "Number of connected components should be 2";
            std::vector<int> homomorphism1 = components[0].getHomomorphism();
            std::vector<int> homomorphism2 = components[1].getHomomorphism();
            expect(homomorphism1[0] == 0 && homomorphism1[1] == 1 && homomorphism1[2] == 2) << "Homomorphism do not correct!";
            expect(homomorphism2[0] == 3 && homomorphism2[1] == 4) << "Homomorphism do not correct!";
            };
        "test_findConnectComponentse_Symmetric "_test = [] {
            // Arrange
            Graph g(5, -5, 5, 0);
            // Assuming you have a method to add edges to the graph
            g.addEdge(0, 1, 3);
            g.addEdge(1, 2, 5);
            g.addEdge(3, 4, 3);

            // Act
            auto components = g.findConnectComponents();

            // Assert
            expect(components.size() == 2 && components[0].size() == 3 && components[1].size() == 2) << fatal << "Number of connected components should be 2";
            for (auto& component : components) {
                for (int i = 0; i < component.size(); i++) {
                    for (int j = i; j < component.size(); j++) {
                        expect(component.getAdjacencyMatrix()[i][j] == component.getAdjacencyMatrix()[j][i]);
                        if (i == j) {
                            expect(component.getAdjacencyMatrix()[i][i] == 0);
                        }
                    }
                }
            }
            for (int i = 0; i < components[0].size(); i++) {
                for (int j = i + 1; j < components[0].size(); j++) {
                    if (i == 0 && j == 1) {
                        expect(components[0].getAdjacencyMatrix()[i][j] == 3);
                    }
                    else {
                        if (i == 1 && j == 2) {
                            expect(components[0].getAdjacencyMatrix()[i][j] == 5);
                        }
                        else {
                            expect(components[0].getAdjacencyMatrix()[i][j] == std::numeric_limits<double>::infinity());
                        }
                    }
                }
            }
            };
        "test_findConnectComponents_homomorphismAndWeight"_test = [] {
            // Arrange
            Graph g(8, -5, 5, 0);
            // Assuming you have a method to add edges to the graph
            g.addEdge(0, 5, 3);
            g.addEdge(1, 2, 1.5);
            g.addEdge(1, 6, -3);
            g.addEdge(1, 7, -2);
            g.addEdge(3, 4, 5);
            // Act
            auto components = g.findConnectComponents();

            // Assert
            expect(components.size() == 3 && components[0].size() == 2 && components[1].size() == 4 && components[2].size() == 2) << fatal << "Number of connected components should be 3";

            std::vector<int> homomorphism1 = components[0].getHomomorphism();
            std::vector<int> homomorphism2 = components[1].getHomomorphism();
            std::vector<int> homomorphism3 = components[2].getHomomorphism();
            expect(homomorphism1[0] == 0 && homomorphism1[1] == 5) << "Homomorphism do not correct!";
            expect(homomorphism2[0] == 1 && homomorphism2[1] == 2 && homomorphism2[2] == 6 && homomorphism2[3] == 7) << "Homomorphism do not correct!";
            expect(homomorphism3[0] == 3 && homomorphism3[1] == 4) << "Homomorphism do not correct!";
            expect(components[0].findGraphWeight() == 3);
            expect(components[1].findGraphWeight() == -3.5);
            expect(components[2].findGraphWeight() == 5);
            };
        "test_findConnectComponents_single_vertex"_test = [] {
            // Arrange
            Graph g(1, -5, 5, 0);
            // Act
            auto components = g.findConnectComponents();

            // Assert
            expect(components.size() == 1 && components[0].size() == 1) << "Number of connected components should be 1 for a single vertex graph";
            };

        "test_findConnectComponents_disconnected_vertexes"_test = [] {
            // Arrange
            Graph g(2, -5, 5, 0);

            // Act
            auto components = g.findConnectComponents();

            // Assert
            expect(components.size() == 2 && components[0].size() == 1 && components[1].size() == 1) << "Number of connected components should be 2 for a graph with two disconnected vertexes with size 1";
            };

        "test_findConnectComponents_connected_vertexes"_test = [] {
            // Arrange
            Graph g(2, -5, 5, 0);
            g.addEdge(0, 1, 3);

            // Act
            auto components = g.findConnectComponents();

            // Assert
            expect(components.size() == 1 && components[0].size() == 2) << "Number of connected components should be 1 for a graph with two connected vertices";
            };
    };
 
    suite<"Create_Random_Graph_Test"> components_tests2 = [] {
        "test_createRandomGraph"_test = [] {
            int n = 5;
            double minWeight = 1.0;
            double maxWeight = 10.0;
            unsigned seed = 123;
            double edgesPercent = 50.0;

            Graph g = Graph::createRandomGraph(n, minWeight, maxWeight, seed, edgesPercent);

            expect( g.matrix.size() == n);

            int edgeCount = 0;
            for (const auto& row : g.matrix) {
                for (const auto& weight : row) {
                    if (weight != std::numeric_limits<double>::infinity() && weight != 0) {
                        edgeCount++;
                        expect( weight >= minWeight && weight <= maxWeight);
                    }
                }
            }
            edgeCount /= 2; // Since the graph is undirected

            int expectedEdges = static_cast<int>((n * (n - 1) / 2.0) * (edgesPercent / 100.0));
            expect( edgeCount == expectedEdges);

            for (int i = 0; i < n; i++) {
                expect( g.matrix[i][i] == 0);
                for (int j = i + 1; j < n; j++) {
                    expect( g.matrix[i][j] == g.matrix[j][i]);
                }
            }
            };

    };
    suite<"Constructor_Graph_Test"> components_tests3= [] {
        "test_constructorGraph"_test = [] {
            Graph g(8, -5, 5, 0);
            // Assuming you have a method to add edges to the graph
            g.addEdge(0, 5, 3);
            g.addEdge(1, 2, 1.5);
            g.addEdge(1, 6, -3);
            g.addEdge(1, 7, -2);
            g.addEdge(3, 4, 5);
            // Act
            auto components = g.findConnectComponents();
            Graph copyFromComponents(components);   
            expect(copyFromComponents == g);
            };
        };
    suite<"Find_Min_Spanning_Forest_Test"> components_tests4 = [] {
        "test_minSpanningForest_1"_test = [] {
            // Create a graph
            Graph g(5,-5,5,0);
            g.addEdge(0, 1, 1);
            g.addEdge(1, 2, 2);
            g.addEdge(2, 0, 3);
            g.addEdge(3, 4, 4);

            // Find the minimum spanning forest
            std::vector<Graph> forest = g.findMinSpanningForest();

            // Check the number of trees in the forest
            expect(forest.size() == 2)<< "Number of trees in the forest must be 2";

            // Check the total weight of each tree
            double totalWeight1 = forest[0].findGraphWeight();
            double totalWeight2 = forest[1].findGraphWeight();
            expect(totalWeight1 == 3)<< "Total weight of the first tree is not correct";
            expect(totalWeight2 == 4)<< "Total weight of the second tree is correct";
            };
        "test_minSpanningForest_2"_test = [] {
            // Create a graph
            Graph g(3, -5, 5, 0);
            g.addEdge(0, 1, 1);
            g.addEdge(1, 2, 2);
            g.addEdge(2, 0, 3);

            std::vector<Graph> forest = g.findMinSpanningForest();

            expect( forest.size() == 1);
            expect( forest[0].findGraphWeight() == 3);
            };
        "test_minSpanningForest_3"_test = [] {
            // Create a graph
            Graph g(7, -5, 5, 0);
          
            g.addEdge(0, 1, 1);
            g.addEdge(1, 2, 2);
            g.addEdge(2, 0, 3);
            g.addEdge(3, 4, 4);
            g.addEdge(5, 6, 5);

            std::vector<Graph> forest = g.findMinSpanningForest();

            expect( forest.size() == 3);
            expect(forest[0].findGraphWeight() == 3);
            expect(forest[1].findGraphWeight() == 4);
            expect(forest[2].findGraphWeight() == 5);
            };
        "test_minSpanningForest_4"_test = [] {
            // Create a graph
            Graph g(7, -5, 5, 0);
            g.addEdge(0, 1, 1);
            g.addEdge(1, 2, 2);
            g.addEdge(2, 0, 3);
            g.addEdge(3, 4, 4);
            g.addEdge(5, 6, 5);

            std::vector<Graph> forest = g.findMinSpanningForest();

            expect(forest.size() == 3);
            expect(forest[0].size() == 3);
            expect(forest[1].size() == 2);
            expect(forest[2].size() == 2);
            };
        "test_minSpanningForest_5"_test = [] {
            // Create a graph
            Graph g(7, -5, 5, 0);
            std::vector<Graph> forest = g.findMinSpanningForest();

            expect(forest.size() == 7);
            for (auto& tree : forest) {
                expect(tree.size() == 1);
            }
            };
        };

       
}