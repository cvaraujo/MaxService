#include <sys/stat.h>
#include "headers/Graph.h"
#include "headers/Model.h"

int main(int argc, const char *argv[]) {
    if (argc < 4) {
        mkdir("results", 0777);
        auto *graph = new Graph("Washington-75/10/washington-75-10-4.txt",
                                "Washington-75/10/param-washington-75-10-4.txt", "results/result_teste.txt");
        auto *model = new Model(graph);
        model->initModel();
        model->solve();
        model->showSolution("results/result_teste.txt");
        return 0;
    } else {
        mkdir("results", 0777);
        auto *graph = new Graph(argv[1], argv[2], argv[3]);
        auto *model = new Model(graph);
        model->initModel();
        model->solve();
        model->showSolution(argv[3]);
    }

    return 0;
}
