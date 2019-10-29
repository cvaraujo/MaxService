//
// Created by carlos on 06/03/19.
//

#include "../headers/Model.h"

Model::Model(Graph *graph) {
    if (graph != nullptr) {
        this->graph = graph;
        initialize();
    } else exit(EXIT_FAILURE);

}

void Model::initialize() {
    int o, d, n = graph->getN(), m = graph->getM();
    try {
        env.set("LogFile", "MS_mip.log");
        env.start();

        f = vector<vector<vector<GRBVar>>>(n, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));
        y = vector<vector<GRBVar>>(n, vector<GRBVar>(n));
        z = vector<GRBVar>(n);

        char name[20];
        for (auto *arc : graph->arcs) {
            o = arc->getO(), d = arc->getD();
            sprintf(name, "y%d%d", o, d);
            y[o][d] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
            for (int k: graph->DuS) {
                sprintf(name, "f%d%d%d", o, d, k);
                this->f[o][d][k] = model.addVar(0.0, 1.0, 0, GRB_CONTINUOUS, name);
            }
        }

        for (int i = 0; i < n; i++) {
            sprintf(name, "z%d", i);
            z[i] = model.addVar(0.0, 1.0, 0, GRB_BINARY, name);
        }
        model.update();
    } catch (GRBException &ex) {
        cout << ex.getMessage() << endl;
        cout << ex.getErrorCode() << endl;
        exit(EXIT_FAILURE);
    }
}

void Model::initModel() {
    cout << "Begin the model creation" << endl;
    objectiveFunction();
    rootFlow(), flowConservation(), terminalsFlow();
    relXandY(), maxArcs();
    limDelayAndJitter(), limVariation();
    cout << "All done!" << endl;
}

void Model::objectiveFunction() {
    GRBLinExpr objective;
    for (auto k : graph->terminals) objective += z[k];
    model.setObjective(objective, GRB_MINIMIZE);
    cout << "Objective Function was added successfully!" << endl;
}

void Model::rootFlow() {
    int o, d, root = graph->getRoot();
    for (auto k : graph->DuS) {
        GRBLinExpr flowExpr, rootExpr;
        for (auto *arc : graph->arcs) {
            o = arc->getO(), d = arc->getD();
            if (o == root) flowExpr += f[root][d][k];
            if (d == root) rootExpr += f[o][root][k];
        }
        model.addConstr((flowExpr - rootExpr) == 1, "root_flow");
    }
    model.update();
    cout << "Flow on root node" << endl;
}

void Model::flowConservation() {
    int o, d, root = graph->getRoot();
    for (auto k : graph->DuS) {
        for (int j = 0; j < graph->getN(); j++) {
            if (j != root && j != k) {
                GRBLinExpr flowIn, flowOut;
                for (auto *arc : graph->arcs) {
                    o = arc->getO(), d = arc->getD();
                    if (o == j) flowOut += f[j][d][k];
                    if (d == j) flowIn += f[o][j][k];
                }
                model.addConstr((flowIn - flowOut) == 0, "flow_conservation");
            }
        }
    }
    model.update();
    cout << "Flow conservation" << endl;
}

void Model::terminalsFlow() {
    int o, d;
    for (auto k : graph->DuS) {
        GRBLinExpr flowIn, flowOut;
        for (auto *arc : graph->arcs) {
            o = arc->getO(), d = arc->getD();
            if (o == k) flowOut += f[k][d][k];
            if (d == k) flowIn += f[o][k][k];
        }
        model.addConstr((flowOut - flowIn) == -1, "flow_on_terminals");
    }
    model.update();
    cout << "Flow on terminals" << endl;
}

void Model::relXandY() {
    int o, d;
    for (auto *arc : graph->arcs) {
        o = arc->getO(), d = arc->getD();
        for (auto k : graph->DuS) model.addConstr(f[o][d][k] <= y[o][d], "f_and_y_relation");
    }
    model.update();
    cout << "f and Y relation" << endl;
}

void Model::maxArcs() {
    GRBLinExpr totalArcs;
    for (auto *arc : graph->arcs)
        totalArcs += y[arc->getO()][arc->getD()];
    model.addConstr(totalArcs == graph->getN() - 1, "maximum_of_arcs");

    model.update();
    cout << "maximum of arcs in the tree" << endl;
}

void Model::limDelayAndJitter() {
    int o, d, paramDelay, paramJitter;
    for (auto k : graph->terminals) {
        GRBLinExpr limDelay, limJitter;
        for (auto *arc : graph->arcs) {
            o = arc->getO(), d = arc->getD();
            limDelay += arc->getDelay() * f[o][d][k];
            limJitter += arc->getJitter() * f[o][d][k];
        }
        paramDelay = graph->getParamDelay(), paramJitter = graph->getParamJitter();
        model.addConstr(limDelay <= (paramDelay + (graph->getBigMDelay() - paramDelay) * z[k]), "delay_limit");
        model.addConstr(limJitter <= (paramJitter + (graph->getBigMJitter() - paramJitter) * z[k]), "jitter_limit");
    }
    model.update();
    cout << "Delay and Jitter limits" << endl;
}

void Model::limVariation() {
    int o, d, bigMK, bigML;
    for (auto k : graph->terminals) {
        for (auto l : graph->terminals) {
            if (k != l) {
                GRBLinExpr delayVariation;
                for (auto *arc : graph->arcs) {
                    o = arc->getO(), d = arc->getD();
                    delayVariation += arc->getDelay() * (f[o][d][k] - f[o][d][l]);
                }
                bigMK = graph->getBigMDelay() -
                        min(graph->getShpTerminal(l) + graph->getParamVariation(), graph->getParamDelay());
                bigML = graph->getParamDelay() - graph->getParamVariation() - graph->getShpTerminal(l);
                model.addConstr(delayVariation <= graph->getParamVariation() + bigMK * z[k] + bigML * z[l],
                                "limit_of_variation_between_pairs");
            }
        }
    }
    model.update();
    cout << "Delay variation limits" << endl;
}

void Model::solve() {
    model.write("model.lp");
    try {
        model.optimize();
    } catch(GRBException &ex) {
        cout << ex.getMessage() << endl;
    }

}

void Model::showSolution() {
    cout << model.get(GRB_DoubleAttr_ObjVal) << endl;
    for (int i = 0; i < graph->getN(); i++)
        cout << z[i].get(GRB_DoubleAttr_X) << endl;

}