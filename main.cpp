#include <sys/stat.h>
#include "headers/Graph.h"
#include "headers/Model.h"

int main(int argc, const char *argv[]) {
    if (argc < 4) {
        return 0;
    } else {
        mkdir("results", 0777);
        auto *graph = new Graph(argv[2], argv[3], argv[4]);
        string usePrep = "1";
        bool use = usePrep.compare(argv[1]);
        auto *model = new Model(graph, use);
        graph->showGraph();
        model->initModel();
        model->solve();
        model->showSolution(argv[4]);
    }

    return 0;
}
