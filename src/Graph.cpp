//
// Created by carlos on 05/03/19.
//

#include "../headers/Data.h"

Data::Data(const char *instance, const char *param) {
    FILE *archiveInstance, *archiveParam;
    char *line = nullptr;
    size_t length = 0;
    char e;
    int v1, v2, vTerminal;
    double delay, jitter, bandwidth, estimateLinkDuration;

    archiveInstance = fopen(instance, "r");
    archiveParam = fopen(param, "r");

    cout << "Loading Archives" << endl;
    if (archiveInstance == nullptr || archiveParam == nullptr) {
        cout << "Error, Cannot open the file!" << endl;
        exit(EXIT_FAILURE);
    }

    getline(&line, &length, archiveParam);
    this->paramDelay = int(1e5 * atof(strtok(line, "Delay limit: ")));

    getline(&line, &length, archiveParam);
    this->paramJitter = int(1e6 * atof(strtok(line, "Jitter limit: ")));

    getline(&line, &length, archiveParam);
    this->paramDelayVariation = int(1e5 * atof(strtok(line, "Delay variation limit:  ")));

    getline(&line, &length, archiveParam);
    this->paramBandwidth = atoi(strtok(line, "Bandwidth limit: "));

    while ((getline(&line, &length, archiveInstance)) != -1) {
        if (strcasecmp("SECTION Graph\n", line) == 0) {
            getline(&line, &length, archiveInstance);
            this->cntNodes = atoi(strtok(line, "Nodes ")) + 1;
            getline(&line, &length, archiveInstance);
            this->cntEdges = atoi(strtok(line, "Edges "));
            getline(&line, &length, archiveInstance);

            this->arcs = vector<Arc *>();
            this->nonDirectedArcs = vector<Arc *>();
            this->graphDelay = BoostGraph(cntNodes);
            this->graphJitter = BoostGraph(cntNodes);
            for (int i = 0; i < this->cntEdges; i++) {
                fscanf(archiveInstance, "%c %d %d %lf %lf %lf %lf\n", &e, &v1, &v2, &delay, &jitter, &bandwidth,
                       &estimateLinkDuration);
                if (bandwidth >= this->paramBandwidth) {
                    Arc *arc = new Arc(v1, v2, int(1e5 * delay), int(1e6 * jitter), int(bandwidth),
                                       int(10 * estimateLinkDuration));
                    Arc *arcRev = new Arc(v2, v1, int(1e5 * delay), int(1e6 * jitter), int(bandwidth),
                                          int(10 * estimateLinkDuration));
                    arcs.push_back(arc), arcs.push_back(arcRev), nonDirectedArcs.push_back(arc);

                    //insertEdge(v1, v2, int(1e6 * delay), true);
                    //insertEdge(v1, v2, int(1e7 * jitter), false);
                }
            }

        } else if ((strcasecmp("SECTION Terminals\n", line) == 0)) {
            getline(&line, &length, archiveInstance);
            this->root = atoi(strtok(line, "Root "));
            getline(&line, &length, archiveInstance);
            this->cntTerminals = atoi(strtok(line, "Terminals "));


            this->terminals = vector<int>(this->cntTerminals);

            for (int i = 0; i < this->cntTerminals; ++i) {
                fscanf(archiveInstance, "%c %d\n", &e, &vTerminal);
                this->terminals[i] = vTerminal;
            }
        }
    }

    graphAdapt();
    cout << "All done!" << endl;
    fclose(archiveParam);
    fclose(archiveInstance);
}

void Data::insertEdge(int i, int j, int weight, bool isDelay) {
    if (isDelay) add_edge(i, j, weight, graphDelay), add_edge(j, i, weight, graphDelay);
    else add_edge(i, j, weight, graphJitter), add_edge(j, i, weight, graphJitter);
}

void Data::graphAdapt() {
    Arc *arc = new Arc(this->root, 0, 0, 0, 0, 0);
    this->nonDirectedArcs.push_back(arc);
    this->arcs.push_back(arc);
    DuS = vector<int>(terminals);

    bool isTerminal;
    for (int i = 1; i < this->cntNodes; ++i) {
        isTerminal = false;
        if (i == this->root) continue;
        for (int t : this->terminals)
            if (i == t) {
                isTerminal = true;
                break;
            }
        if (!isTerminal) this->nonTerminals.push_back(i);
    }

    DuS.insert(DuS.end(), nonTerminals.begin(), nonTerminals.end());


    for (int s : this->nonTerminals) {
        arc = new Arc(0, s, 0, 0, 0, 0);
        this->nonDirectedArcs.push_back(arc);
        this->arcs.push_back(arc);//, this->arcs.push_back(arcRev);
    }

    this->cntEdges = int(this->arcs.size());
}

void Data::showData() {
    cout << "Arcs" << endl;
    for (auto *arc : arcs) {
        cout << arc->getO() << " " << arc->getD() << ": " << arc->getDelay()
             << " " << arc->getJitter() << " " << arc->getBandwidth() << " " <<
             arc->getEstimateLinkDuration() << endl;
    }

    cout << "\n Param" << endl;
    cout << "Nodes: " << cntNodes << " Edges: " << cntEdges <<
         " CntTerminals: " << cntTerminals << " Root: " << root << endl;

    cout << "Delay: " << paramDelay << " Jitter: " << paramJitter <<
         " DelayVari.: " << paramDelayVariation << " Bandwidth: " << paramBandwidth << endl;

    cout << "Terminals" << endl;
    for (int i : terminals) {
        cout << "T: " << i << " ";
    }
    cout << "\nNonTerminals" << endl;
    for (int i : nonTerminals) {
        cout << "NT: " << i << " ";
    }

}