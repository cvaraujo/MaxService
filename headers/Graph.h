//
// Created by carlos on 27/05/19.
//

#ifndef MS_GRAPH_H
#define MS_GRAPH_H

#include "Arc.h"
#include "Include.h"

using namespace std;
using namespace boost;

class Graph {
    int n, nRemoved, m, paramDelay, paramJitter, paramVariation, paramBandwidth, root, bigMDelay = 0, bigMJitter = 0;

    BoostGraph graphDelaySP;

    vector<VertexDescriptor> predecessors;
    vector<int> distanceDelay, distanceJitter;

public:
    vector<vector<Arc *>> arcs;
    vector<int> terminals, nonTerminals, DuS, delayVector, jitterVector;
    vector<bool> removed;
    vector<vector<vector<bool>>> removedF;
    vector<bool> notAttend;

    Graph(string instance, string param, string outputName);

    void graphReduction();

    void showGraph();

    int getN() const;

    int getNAfterRemove() const;

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