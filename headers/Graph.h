//
// Created by carlos on 27/05/19.
//

#ifndef MS_GRAPH_H
#define MS_GRAPH_H

#include "Arc.h"
#include "Include.h"

using namespace std;

class Graph {
    typedef graph_traits <BoostGraph>::vertex_descriptor vertex_descriptor;

    int n, m, paramDelay, paramJitter, paramVariation, paramBandwidth, root, bigMDelay = 0, bigMJitter = 0;

    BoostGraph graphShp;

    vector<vertex_descriptor> predecessors;
    vector<int> distance;

public:
    vector<vector<Arc *>> arcs;
    vector<int> terminals, nonTerminals, DuS, delayVector, jitterVector;
    vector<bool> removed;

    Graph(string instance, string param);

    void showGraph();

    int getN() const;

    void setN(int n);

    int getM() const;

    void setM(int m);

    int getParamDelay() const;

    void setParamDelay(int paramDelay);

    int getParamJitter() const;

    void setParamJitter(int paramJitter);

    int getParamVariation() const;

    void setParamVariation(int paramVariation);

    int getParamBandwidth() const;

    void setParamBandwidth(int paramBandwidth);

    int getRoot() const;

    void setRoot(int root);

    int getBigMDelay();

    int getBigMJitter();

    int getShpTerminal(int k);
};


#endif