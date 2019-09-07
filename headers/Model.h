//
// Created by carlos on 06/03/19.
//

#ifndef MRP_MODEL_H
#define MRP_MODEL_H

#include "Include.h"
#include "Data.h"

class Model {
    typedef IloArray<IloNumVarArray> NumVarMatrix;
    Data *data;
    IloEnv env;
    IloModel model;
    IloCplex cplex;
    IloArray<NumVarMatrix> f;
    IloArray<IloNumVarArray> x;
    IloArray<IloNumVarArray> y;

    void objectiveFunction(double thetaC, double thetaP, double thetaD);

    void c2();

    void c3();

    void c4();

    void c5();

    void c6();

    void c7();

    void c8c9();

    void c10();

    void c11();

    void c12();

public:
    Model(Data *data);

    void initialize();

    void initModel(double thetaC, double thetaP, double thetaD);

    void solve();

    void showSolution(const char *input, const char *outputFile, double thetaC, double thetaP, double thetaD);

};


#endif //MRP_MODEL_H
