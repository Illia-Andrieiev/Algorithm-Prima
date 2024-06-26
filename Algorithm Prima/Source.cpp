#include <iostream>
#include <limits>
#include"Graph.h"
#include"Proxy.h"
int main(int argc, char* argv[]) {
	//double inf = std::numeric_limits<double>::infinity();
	//if (inf == std::numeric_limits<double>::infinity())
	//	inf = 5;
	//std::cout << "Infinity: " << inf << std::endl;
	//Graph g = Graph::createRandomGraph(10, -5, 5, 15);
	//g.printEdges();
	//std::cout << std::endl << "********************************************" << std::endl;
	//std::vector<Graph> compn = g.findConnectComponents();
	//for (auto& comp : compn) {
	//	comp.printEdges();
	//	std::cout<< std::endl<<"*******************************************" << std::endl;
	//}
	//Graph fromComponents(compn);
	////fromComponents.print();
	//std::cout << std::endl << "********************************************" << std::endl;
	//std::cout << std::endl << "********************************************" << std::endl;
	//std::vector<Graph> forest = g.parallel_findMinSpanningForest();
	//for (auto& tree : forest) {
	//	tree.printEdges();
	//	std::cout << std::endl << "*******************************************" << std::endl;
	//}
	//Graph::sameGraphBenchMark(2000,0.001, 1);
	Graph g(10, 1, 10, 20);
	GraphProxyLogger log(g, "loggerProxy.txt");
	log.findConnectComponents();
	log.getHomomorphism();
	//delete g;
	Graph::test();
	return 0;
}