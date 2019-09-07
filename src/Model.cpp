//
// Created by carlos on 06/03/19.
//

#include "../headers/Model.h"

Model::Model(Data *data) {
    if (data != nullptr) {
        this->data = data;

        initialize();
    }
}

void Model::initialize() {
    try {
        char name[12];
        this->model = IloModel(env);
        this->cplex = IloCplex(model);

        this->f = IloArray<NumVarMatrix>(env, data->cntNodes);
        this->x = IloArray<IloNumVarArray>(env, data->cntNodes);
        this->y = IloArray<IloNumVarArray>(env, data->cntNodes);

        for (int i = 0; i < data->cntNodes; i++) {
            this->x[i] = IloNumVarArray(env, data->cntNodes);
            this->f[i] = NumVarMatrix(env, data->cntNodes);
            this->y[i] = IloNumVarArray(env, data->cntNodes);
            for (int j = 0; j < data->cntNodes; j++) {
                this->f[i][j] = IloNumVarArray(env, data->cntNodes);
            }
        }

        // Create the variable f^{ijk}
        for (auto *arc : data->arcs) {
            for (int k: data->DuS) {
                sprintf(name, "f%d%d%d", arc->getO(), arc->getD(), k);
                this->f[arc->getO()][arc->getD()][k] = IloNumVar(this->env, 0, 1, name);

                this->model.add(this->f[arc->getO()][arc->getD()][k]);
                this->model.add(IloConversion(env, this->f[arc->getO()][arc->getD()][k], ILOINT));
                // Linear Relaxation
                //this->model.add(IloConversion(env, this->f[arc->getO()][arc->getD()][k], ILOFLOAT));

            }
        }

        // Create the variable X^{ij}

        for (auto *arc : data->nonDirectedArcs) {
            sprintf(name, "x%d%d", arc->getO(), arc->getD());
            this->x[arc->getO()][arc->getD()] = IloNumVar(this->env, 0, 1, name);

            this->model.add(this->x[arc->getO()][arc->getD()]);
            this->model.add(IloConversion(env, this->x[arc->getO()][arc->getD()], ILOINT));
            // Linear Relaxation
            //this->model.add(IloConversion(env, this->x[arc->getO()][arc->getD()], ILOFLOAT));
        }

        // Create the variable Y_{ij}
        for (auto *arc : data->arcs) {
            sprintf(name, "y%d%d", arc->getO(), arc->getD());
            this->y[arc->getO()][arc->getD()] = IloNumVar(this->env, 0, 1, name);
            this->y[arc->getD()][arc->getO()] = IloNumVar(this->env, 0, 1, name);

            this->model.add(this->y[arc->getO()][arc->getD()]);
            this->model.add(IloConversion(env, this->y[arc->getO()][arc->getD()], ILOINT));
            //this->model.add(IloConversion(env, this->y[arc->getO()][arc->getD()], ILOFLOAT));

            this->model.add(this->y[arc->getD()][arc->getO()]);
            this->model.add(IloConversion(env, this->y[arc->getD()][arc->getO()], ILOINT));
            //this->model.add(IloConversion(env, this->y[arc->getD()][arc->getO()], ILOFLOAT));
        }


    } catch (IloException &ex) {
        cout << ex.getMessage() << endl;
        exit(EXIT_FAILURE);
    }

}

void Model::initModel(double thetaC, double thetaP, double thetaD) {
    cout << "Begin the model creation" << endl;
    cplex.setParam(IloCplex::Param::TimeLimit, 600);
    cplex.setParam(IloCplex::TreLim, 7000);
    cplex.setOut(env.getNullStream());
    objectiveFunction(thetaC, thetaP, thetaD);
    c2();
    c3();
    c4();
    c5();
    c6();
    c7();
    c8c9();
    c10();
    c11();
    c12();
    cout << "All done!" << endl;
}

void Model::objectiveFunction(double thetaC, double thetaP, double thetaD) {
    IloExpr objExpr(env);
    double mi;
    int i, j;
    for (auto *arc : data->nonDirectedArcs) {
        i = arc->getO(), j = arc->getD();
        if (i != 0 && j != 0) {
            mi = thetaC + (thetaP * arc->getDelay()) - (thetaD * arc->getEstimateLinkDuration());
            objExpr += this->x[i][j] * mi;
        }
    }
    IloObjective obj = IloObjective(env, objExpr, IloObjective::Minimize);
    this->model.add(obj);

    cout << "Objective Function was added successfully!" << endl;
}

void Model::c2() {
    for (int k : data->DuS) {
        IloExpr aux(env);
        IloExpr aux2(env);
        for (auto *arc : data->arcs) {
            if (arc->getO() == data->root) {
                aux += this->f[data->root][arc->getD()][k];
            } else if (arc->getD() == data->root) {
                aux2 += this->f[arc->getO()][data->root][k];
            }
        }
        this->model.add((aux - aux2) == 1);
        aux.end();
        aux2.end();
    }
    cout << "Constraint C2 was added successfully!" << endl;
}

void Model::c3() {
    for (int k : data->DuS) {
        for (int j = 0; j < data->cntNodes; j++) {
            if (j != k && j != data->root) {
                IloExpr aux(env);
                IloExpr aux2(env);
                for (auto *arc : data->arcs) {
                    if (arc->getD() == j)
                        aux += this->f[arc->getO()][j][k];
                    if (arc->getO() == j)
                        aux2 += this->f[j][arc->getD()][k];
                }
                this->model.add((aux - aux2) == 0);
                aux.end();
                aux2.end();
            }
        }
    }
    cout << "Constraint C3 was added successfully!" << endl;
}

void Model::c4() {
    int i, j;
    for (int k : data->DuS) {
        IloExpr aux(env);
        IloExpr aux2(env);
        for (auto *arc : data->arcs) {
            i = arc->getO(), j = arc->getD();
            if (i == k) aux += this->f[k][j][k];
            if (j == k) aux2 += this->f[i][k][k];
        }
        this->model.add((aux - aux2) == -1);
        aux.end();
        aux2.end();
    }

    cout << "Constraint C4 was added successfully!" << endl;
}

void Model::c5() {
    int i, j;
    for (auto *arc : data->arcs) {
        i = arc->getO(), j = arc->getD();
        for (int k : data->DuS) {
            model.add(this->f[i][j][k] <= this->y[i][j]);
        }
    }
    cout << "Constraint C5 was added successfully!" << endl;
}

void Model::c6() {
    IloExpr c6(env);
    for (auto *arc : data->arcs) {
        c6 += y[arc->getO()][arc->getD()];
    }
    model.add(c6 == (data->cntNodes - 1));

    cout << "Constraint C6 was added successfully!" << endl;
}

void Model::c7() {
    int i, j;
    for (auto arc : data->nonDirectedArcs) {
        i = arc->getO(), j = arc->getD();
        model.add(y[i][j] + y[j][i] == x[i][j]);
    }
    cout << "\nConstraint C7 was added successfully!" << endl;
}

void Model::c8c9() {
    for (int k : data->terminals) {
        IloExpr c8(env);
        IloExpr c9(env);
        for (auto *arc : data->arcs) {
            c8 += arc->getDelay() * this->f[arc->getO()][arc->getD()][k];
            c9 += arc->getJitter() * this->f[arc->getO()][arc->getD()][k];
        }
        this->model.add(c8 <= data->paramDelay);
        this->model.add(c9 <= data->paramJitter);

        c8.end();
        c9.end();
    }
    cout << "Constraints C8 and C9 was added successfully!" << endl;
}

void Model::c10() {
    for (auto k : data->terminals) {
        for (auto l : data->terminals) {
            if (k != l) {
                IloExpr c10(env);
                for (auto *arc : data->arcs) {
                    c10 += arc->getDelay() *
                           (this->f[arc->getO()][arc->getD()][k] - this->f[arc->getO()][arc->getD()][l]);
                }
                this->model.add(c10 <= data->paramDelayVariation);
                c10.end();
            }
        }
    }
    cout << "Constraint C10 was added successfully!" << endl;
}

void Model::c11() {
    for (int k : data->terminals) {
        this->model.add(this->f[data->root][0][k] == 0);
    }
    cout << "Constraint C11 was added successfully!" << endl;
}

void Model::c12() {
    for (int q : data->nonTerminals) {
        for (int e : data->DuS) {
            if (e != q)
                this->model.add(this->f[0][q][e] == 0);
        }
    }
    cout << "Constraint C12 was added successfully!" << endl;
}

void Model::solve() {
    this->cplex.exportModel("model.lp");
    this->cplex.solve();
}

void Model::showSolution(const char *input, const char *outputFile, double thetaC, double thetaP, double thetaD) {
    try {
        FILE *output;
        output = fopen(outputFile, "a");
        fprintf(output, "%s: ", input);
        double fo = 0;
        int i, j;
        for (auto *arc : data->nonDirectedArcs) {
            i = arc->getO(), j = arc->getD();
            if (i != 0 && j != 0 && cplex.getValue(this->x[i][j]) > 0.5) {
                printf("[%d, %d], ", i, j);
                fo += thetaC + arc->getDelay() * thetaP - arc->getEstimateLinkDuration() * thetaD;
            }
        }
        cout << "\n ---- CPLEX_FO = " << this->cplex.getObjValue() << " ---- SCALE_FO = " << fo
             << " ---- Time = "
             << this->cplex.getCplexTime() << " ---- GAP = " << this->cplex.getMIPRelativeGap()
             << " ---- LB = " << cplex.getBestObjValue() << endl;

        fprintf(output,
                " ---- CPLEX_FO = %g ---- FO = %lf ---- Time = %g ---- GAP = %g ---- LB = %g\n",
                cplex.getObjValue(), fo, cplex.getTime(), cplex.getMIPRelativeGap(), cplex.getBestObjValue());
        fclose(output);
    } catch (IloException &exception) {
        FILE *output;
        output = fopen(outputFile, "a");
        fprintf(output, "%s: ", input);
        fprintf(output, " ---- Time = %g ---- LB = %g\n", cplex.getTime(), cplex.getBestObjValue());
        fclose(output);
        cout << "No Feasible Solution " << exception.getMessage() << endl;
    }
}