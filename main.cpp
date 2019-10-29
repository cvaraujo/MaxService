#include "headers/Graph.h"
#include "headers/Model.h"

int main(int argc, const char *argv[]) {
    if (argc < 4) {
        auto *graph = new Graph("10/washington-50-10-3.txt",
                                "10/param-washington-50-10-3.txt");

        auto *model = new Model(graph);
//        getchar();
        model->initModel();
        model->solve();
        model->showSolution();
        return 0;
    }
//    } else {
//        auto *data = new Graph(argv[1], argv[2]);
//        auto *model = new Model(data);
//        //data->showData();
//        model->initialize();
//        model->initModel(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
//        model->solve();
//        model->showSolution(argv[1], "Results_model_integer.txt", atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
//    }

    return 0;
}
