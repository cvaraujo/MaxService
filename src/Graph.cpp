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
            nRemoved = n;
            graphDelaySP = BoostGraph(n);
            arcs = vector < vector < Arc * >> (n, vector<Arc *>());
            removed = vector<bool>(n);
            notAttend = vector<bool>(n);
            removedF = vector < vector < vector < bool >> > (n, vector < vector < bool >> (n, vector<bool>(n)));
        }

        if (token == "Edges") {
            fileGraph >> m;
            output << m << "\n";
        }

        if (token == "E") {
            fileGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
            u--, v--;
            if (bandwidth >= paramBandwidth) {
                delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter);
                Arc *arc = new Arc(u, v, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                Arc *arcRev = new Arc(v, u, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                delayVector.push_back(delayInt), jitterVector.push_back(jitterInt);
                arcs[u].push_back(arc), arcs[v].push_back(arcRev);
                add_edge(u, v, delayInt, graphDelaySP), add_edge(v, u, delayInt, graphDelaySP);
            }
        }
        if (token == "Root") fileGraph >> root, root--;
        if (token == "T") fileGraph >> u, u--, terminals.push_back(u), DuS.push_back(u);
    }

    property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, graphDelaySP);
    predecessors = vector<VertexDescriptor>(n);
    distanceDelay = vector<int>(n), distanceJitter = vector<int>(n);

    dijkstra_shortest_paths(graphDelaySP, root, predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, graphDelaySP))).distance_map(
            make_iterator_property_map(distanceDelay.begin(), get(vertex_index, graphDelaySP))));

    int cntRemoved = n;
    for (int i = 0; i < n; i++) {
        removed[i] = distanceDelay[i] > paramDelay;
        if (removed[i]) nRemoved--;
    }

    output << nRemoved << "\n";

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

    sort(delayVector.begin(), delayVector.end(), greater<int>());
    sort(jitterVector.begin(), jitterVector.end(), greater<int>());

    for (int i = 0; i < nRemoved-1; i++)
        bigMDelay += delayVector[i], bigMJitter += jitterVector[i];

    output.close();
    cout << "Load graph successfully" << endl;
}

void Graph::showGraph() {
    cout << "Arcs" << endl;
    for (int o = 0; o < n; o++) {
        for (auto *arc : arcs[o])
            cout << o + 1 << " " << arc->getD() + 1 << ": " << arc->getDelay()
                 << " " << arc->getJitter() << " " << arc->getBandwidth() << " " <<
                 arc->getEstimateLinkDuration() << endl;

    }

    cout << "\n Param" << endl;
    cout << "Nodes: " << n << " Edges: " << m <<
         " CntTerminals: " << int(terminals.size()) << " Root: " << root + 1 << endl;

    cout << "Delay: " << paramDelay << " Jitter: " << paramJitter <<
         " DelayVari.: " << paramVariation << " Bandwidth: " << paramBandwidth << endl;

    cout << "Terminals" << endl;
    for (int i : terminals) {
        cout << "T: " << i + 1 << " ";
    }
    cout << "\nNonTerminals" << endl;
    for (int i : nonTerminals) {
        cout << "NT: " << i + 1 << " ";
    }
}

int Graph::getBigMDelay() {
    return bigMDelay;
}

int Graph::getBigMJitter() {
    return bigMJitter;
}

int Graph::getShpTerminal(int k) {
    return distanceDelay[k];
}

int Graph::getN() const {
    return n;
}

void Graph::setN(int n) {
    Graph::n = n;
}

int Graph::getNAfterRemove() const {
    return nRemoved;
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

