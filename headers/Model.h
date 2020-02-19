//
// Created by carlos on 06/03/19.
//

#ifndef MRP_MODEL_H
#define MRP_MODEL_H

#include "Include.h"
#include "Graph.h"

class Model {
    Graph *graph;
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    vector<vector<vector<GRBVar>>> f;
    vector<vector<GRBVar>> y;
    vector<GRBVar> z;

    int preprocessingTime;

    void preprocessing();

    void objectiveFunction();

    void rootFlow();

    void flowConservation();

    void terminalsFlow();

    void relXandY();

    void maxArcs();

    void limDelayAndJitter();

    void limVariation();

    void primeToTerminals();

    void nonTerminalsLeafs();

public:
    Model(Graph *graph, bool usePrep);

    void initialize();

    void initModel();

    void solve();

    void showSolution(string instance);

};


#endif //MRP_MODEL_H
