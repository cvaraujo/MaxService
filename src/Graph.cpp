//
// Created by carlos on 27/05/19.
//

#include "../headers/Graph.h"
#include "../headers/Include.h"

Graph::Graph(string instance, string param, string outputName) {
    int u, v;
    double delay, jitter, bandwidth, ldp, paramDelayToken, paramJitterToken, paramVariationToken, paramBandwidthToken;
    int delayInt, jitterInt;
    string token;
    ifstream fileGraph, fileParam;
    ofstream output;

    output.open(outputName);

    fileParam.open(param, fstream::in);

    while (!fileParam.eof()) {
        fileParam >> token;
        if (token == "Delay") {
            fileParam >> token;
            if (token == "variation") {
                fileParam >> token >> paramVariationToken;
                Graph::paramVariation = int(1e5 * paramVariationToken);
            } else {
                fileParam >> paramDelayToken;
                Graph::paramDelay = int(1e5 * paramDelayToken);
            }
        }
        if (token == "Jitter") {
            fileParam >> token >> paramJitterToken;
            Graph::paramJitter = int(1e6 * paramJitterToken);
        }
        if (token == "Bandwidth") {
            fileParam >> token >> paramBandwidthToken;
            Graph::paramBandwidth = int(paramBandwidthToken);
        }
    }

    fileGraph.open(instance, fstream::in);

    while (!fileGraph.eof()) {
        fileGraph >> token;
        if (token == "Nodes") {
            fileGraph >> n;
            output << n << "\n";
            graphShp = BoostGraph(n);
            arcs = vector<vector<Arc *>>(n, vector<Arc *>());
            removed = vector<bool>(n);
        }

        if (token == "Edges") {
            fileGraph >> m;
            output << m << "\n";
        }

        if (token == "E") {
            fileGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
            if (bandwidth >= paramBandwidth) {
                delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter), --u, --v;
                Arc *arc = new Arc(u, v, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                Arc *arcRev = new Arc(v, u, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                arcs[u].push_back(arc), arcs[v].push_back(arcRev);
                add_edge(u, v, delayInt, graphShp), add_edge(v, u, delayInt, graphShp);
            }
        }
        if (token == "Root") fileGraph >> root, root--;
        if (token == "T") fileGraph >> u, u--, terminals.push_back(u), DuS.push_back(u);
    }

    property_map<BoostGraph, edge_weight_t>::type weightMap = get(edge_weight, graphShp);
    predecessors = vector<vertex_descriptor>(n);
    distance = vector<int>(n);

    dijkstra_shortest_paths(graphShp, root, predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, graphShp))).distance_map(
            make_iterator_property_map(distance.begin(), get(vertex_index, graphShp))));

    int cntRemoved = n;
    for (int i = 0; i < n; i++) {
        removed[i] = distance[i] >= INT_MAX;
        if (removed[i]) cntRemoved--;
    }
    output << cntRemoved << "\n";

    bool isTerminal;
    for (int i = 0; i < n; ++i) {
        if (!removed[i]) {
            isTerminal = false;
            if (i != root) {
                for (auto t : terminals) {
                    if (i == t) {
                        isTerminal = true;
                        break;
                    }
                }
                if (!isTerminal) nonTerminals.push_back(i), DuS.push_back(i);
            }
        }
    }

    int o, d, cntRemovedEdge = m;
    for (o = 0; o < n; o++) {
        if (!removed[o]) {
            for (int j = 0; j < int(arcs[o].size()); j++) {
                auto *arc = arcs[o][j];
                d = arc->getD();
                if (removed[d]) {
                    arcs[o].erase(arcs[o].begin() + j);
                    cntRemovedEdge--;
                } else {
                    this->delayVector.push_back(arc->getDelay()), this->jitterVector.push_back(arc->getJitter());
                }
            }
        }
    }

    output << cntRemovedEdge << "\n";

    sort(delayVector.begin(), delayVector.end());
    sort(jitterVector.begin(), jitterVector.end());

    for (int i = 0; i < n; i++)
        bigMDelay += delayVector[i], bigMJitter += jitterVector[i];

//    cout << "BigM Delay = " << bigMDelay << endl;
//    cout << "BigM Jitter = " << bigMJitter << endl;

//    cout << "Distance to terminals" << endl;
//    for (int i = 0; i < n; i++) cout << "k: " << i << ", distance: " << distance[i] << endl;
//    output << "graph ends\n";
    output.close();
    cout << "Load graph successfully" << endl;
//    getchar();
}

void Graph::showGraph() {
    cout << "Arcs" << endl;
    for (int o = 0; o < n; o++) {
        if (!removed[o]) {
            for (auto *arc : arcs[o])
                cout << arc->getO() << " " << arc->getD() << ": " << arc->getDelay()
                     << " " << arc->getJitter() << " " << arc->getBandwidth() << " " <<
                     arc->getEstimateLinkDuration() << endl;
        }
    }

    cout << "\n Param" << endl;
    cout << "Nodes: " << n << " Edges: " << m <<
         " CntTerminals: " << int(terminals.size()) << " Root: " << root << endl;

    cout << "Delay: " << paramDelay << " Jitter: " << paramJitter <<
         " DelayVari.: " << paramVariation << " Bandwidth: " << paramBandwidth << endl;

    cout << "Terminals" << endl;
    for (int i : terminals) {
        cout << "T: " << i << " ";
    }
    cout << "\nNonTerminals" << endl;
    for (int i : nonTerminals) {
        cout << "NT: " << i << " ";
    }
}

int Graph::getBigMDelay() {
    return bigMDelay;
}

int Graph::getBigMJitter() {
    return bigMJitter;
}

int Graph::getShpTerminal(int k) {
    return distance[k];
}

int Graph::getN() const {
    return n;
}

void Graph::setN(int n) {
    Graph::n = n;
}

int Graph::getM() const {
    return m;
}

void Graph::setM(int m) {
    Graph::m = m;
}

int Graph::getParamDelay() const {
    return paramDelay;
}

void Graph::setParamDelay(int paramDelay) {
    Graph::paramDelay = paramDelay;
}

int Graph::getParamJitter() const {
    return paramJitter;
}

void Graph::setParamJitter(int paramJitter) {
    Graph::paramJitter = paramJitter;
}

int Graph::getParamVariation() const {
    return paramVariation;
}

void Graph::setParamVariation(int paramVariation) {
    Graph::paramVariation = paramVariation;
}

int Graph::getParamBandwidth() const {
    return paramBandwidth;
}

void Graph::setParamBandwidth(int paramBandwidth) {
    Graph::paramBandwidth = paramBandwidth;
}

int Graph::getRoot() const {
    return root;
}

void Graph::setRoot(int root) {
    Graph::root = root;
}

