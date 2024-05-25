#include <iostream>
#include <limits>
#include"Graph.h"
int main() {
	double inf = std::numeric_limits<double>::infinity();
	if (inf == std::numeric_limits<double>::infinity())
		inf = 5;
	std::cout << "Infinity: " << inf << std::endl;
	Graph g = Graph::createRandomGraph(5, -5, 5, 20);
	g.print();
	std::cout << std::endl << "********************************************" << std::endl;
	std::vector<Graph> compn = g.findConnectComponents();
	for (auto& comp : compn) {
		comp.print();
		std::cout<< std::endl<<"*******************************************" << std::endl;
	}
	compn.pop_back();
	Graph fromComponents(compn);
	fromComponents.print();
	
	return 0;
}