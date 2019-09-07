#include "headers/Data.h"
#include "headers/Model.h"

int main(int argc, const char *argv[]) {
    if (argc < 4) {
        auto *data = new Data("Instances/30/washington-50-30-10.txt",
                              "Instances/30/param-washington-50-30-10.txt");
        auto *model = new Model(data);
        //data->showData();
        //getchar();
        model->initialize();
        model->initModel(60, 35, 5);
        model->solve();
        model->showSolution("test", "results.txt", 60, 35, 5);
        return 0;
    } else {
        auto *data = new Data(argv[1], argv[2]);
        auto *model = new Model(data);
        //data->showData();
        model->initialize();
        model->initModel(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
        model->solve();
        model->showSolution(argv[1], "Results_model_integer.txt", atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    }

    return 0;
}
