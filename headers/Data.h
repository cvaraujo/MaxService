//
// Created by carlos on 05/03/19.
//

#ifndef MRP_DATA_H
#define MRP_DATA_H

#include "Include.h"
#include "Arc.h"


class Data {
    void graphAdapt();

    void insertEdge(int i, int j, int weight, bool isDelay);

public:
    vector<Arc *> arcs, nonDirectedArcs;
    vector<int> terminals, nonTerminals, DuS;
    int cntNodes, cntEdges, cntTerminals, root = 1;
    int paramDelay, paramJitter, paramDelayVariation, paramBandwidth;

    BoostGraph graphDelay, graphJitter;

    Data(const char *instance, const char *param);

    void showData();
};


#endif //MRP_DATA_H
