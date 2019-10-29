//
// Created by carlos on 06/03/19.
//

#ifndef MRP_MODEL_H
#define MRP_MODEL_H

#include "Include.h"
#include "Graph.h"

class Model {
//    typedef IloArray<IloNumVarArray> NumVarMatrix;
    Graph *graph;
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    vector<vector<vector<GRBVar>>> f;
    vector<vector<GRBVar>> y;
    vector<GRBVar> z;

    void objectiveFunction();

    void rootFlow();

    void flowConservation();

    void terminalsFlow();

    void relXandY();

    void maxArcs();

    void limDelayAndJitter();

    void limVariation();
public:
    Model(Graph *graph);

    void initialize();

    void initModel();

    void solve();

    void showSolution();

};


#endif //MRP_MODEL_H
