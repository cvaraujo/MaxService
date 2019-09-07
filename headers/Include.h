//
// Created by carlos on 05/03/19.
//

#ifndef MRP_INCLUDE_H
#define MRP_INCLUDE_H

#include <iostream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include "string"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include "boost/graph/adjacency_list.hpp"


ILOSTLBEGIN
using namespace std;
using namespace boost;

typedef adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, int>> BoostGraph;
typedef graph_traits<BoostGraph>::edge_descriptor EdgeDescriptor;
typedef graph_traits<BoostGraph>::vertex_descriptor VertexDescriptor;


#endif //MRP_INCLUDE_H
