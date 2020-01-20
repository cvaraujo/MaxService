//
// Created by carlos on 27/05/19.
//

#include "../headers/Graph.h"
#include "../headers/Include.h"

// data structures for shortest path problem with resource constraint
// ResourceContainer model
struct spp_spp_res_cont {
    spp_spp_res_cont(int c = 0, int r = 0) : cost(c), res(r) {}

    spp_spp_res_cont &operator=(const spp_spp_res_cont &other) {
        if (this == &other) return *this;
        this->~spp_spp_res_cont();
        new(this) spp_spp_res_cont(other);
        return *this;
    }

    int cost;
    int res;
};

bool operator==(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) {
    return (res_cont_1.cost == res_cont_2.cost && res_cont_1.res == res_cont_2.res);
}

bool operator<(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) {
    if (res_cont_1.cost > res_cont_2.cost) return false;
    if (res_cont_1.cost == res_cont_2.cost) return res_cont_1.res < res_cont_2.res;
    return true;
}

// ResourceExtensionFunction model
class ref_spprc {
public:
    inline bool operator()(const SPPRCGraph &g, spp_spp_res_cont &new_cont, const spp_spp_res_cont &old_cont,
                           graph_traits<SPPRCGraph>::edge_descriptor ed) const {

        const SPPRC_Graph_Arc &arc_prop = get(edge_bundle, g)[ed];
        const SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, g)[target(ed, g)];
        new_cont.cost = old_cont.cost + arc_prop.cost;
        int &i_res = new_cont.res;
        i_res = old_cont.res + arc_prop.res;
        return i_res <= vert_prop.con;
    }
};

// DominanceFunction model
class dominance_spptw {
public:
    inline bool operator()(const spp_spp_res_cont &res_cont_1, const spp_spp_res_cont &res_cont_2) const {
        // must be "<=" here!!!
        // must NOT be "<"!!!
        return res_cont_1.cost <= res_cont_2.cost && res_cont_1.res <= res_cont_2.res;
        // this is not a contradiction to the documentation
        // the documentation says:
        // "A label $l_1$ dominates a label $l_2$ if and only if both are resident
        // at the same vertex, and if, for each resource, the resource consumption
        // of $l_1$ is less than or equal to the resource consumption of $l_2$,
        // and if there is at least one resource where $l_1$ has a lower resource
        // consumption than $l_2$."
        // one can think of a new label with a resource consumption equal to that
        // of an old label as being dominated by that old label, because the new
        // one will have a higher number and is created at a later point in time,
        // so one can implicitly use the number or the creation time as a resource
        // for tie-breaking
    }
};
// end data structures for shortest path problem with time windows (spptw)

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
            fileGraph >> n, n++;
            output << n << "\n";
            graphDelaySP = BoostGraph(n);
            graphJitterSP = BoostGraph(n);
            arcs = vector<vector<Arc * >>(n, vector<Arc *>());
            removed = vector<bool>(n);
            notAttend = vector<bool>(n);
            removedF = vector<vector<vector<bool >>>(n, vector<vector<bool >>(n, vector<bool>(n)));
        }

        if (token == "Edges") {
            fileGraph >> m;
            output << m << "\n";
        }

        if (token == "E") {
            fileGraph >> u >> v >> delay >> jitter >> bandwidth >> ldp;
            if (bandwidth >= paramBandwidth) {
                delayInt = int(1e5 * delay), jitterInt = int(1e6 * jitter);
                Arc *arc = new Arc(u, v, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                Arc *arcRev = new Arc(v, u, delayInt, jitterInt, int(bandwidth), int(10 * ldp));
                delayVector.push_back(delayInt), jitterVector.push_back(jitterInt);
                arcs[u].push_back(arc), arcs[v].push_back(arcRev);
                add_edge(u, v, delayInt, graphDelaySP), add_edge(v, u, delayInt, graphDelaySP);
                add_edge(u, v, jitterInt, graphJitterSP), add_edge(v, u, jitterInt, graphJitterSP);
            }
        }
        if (token == "Root") fileGraph >> root;
        if (token == "T") fileGraph >> u, terminals.push_back(u), DuS.push_back(u);
    }

    property_map<BoostGraph, edge_weight_t>::type weightMapDelay = get(edge_weight, graphDelaySP);
    predecessors = vector<VertexDescriptor>(n);
    distanceDelay = vector<int>(n), distanceJitter = vector<int>(n);

    dijkstra_shortest_paths(graphDelaySP, root, predecessor_map(
            make_iterator_property_map(predecessors.begin(), get(vertex_index, graphDelaySP))).distance_map(
            make_iterator_property_map(distanceDelay.begin(), get(vertex_index, graphDelaySP))));

    int cntRemoved = n;
    for (int i = 1; i < n; i++) {
        removed[i] = distanceDelay[i] > paramDelay;
        if (removed[i]) cntRemoved--;
    }

    output << cntRemoved << "\n";

    arcs[root].push_back(new Arc(root, 0, paramDelay, paramJitter, 0, 0));

    bool isTerminal;
    for (int i = 1; i < n; ++i) {
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

    nonTerminals.push_back(0), DuS.push_back(0);
    for (auto s : DuS) {
        if (s == 0) continue;
        arcs[0].push_back(new Arc(0, s, paramDelay, paramJitter, 0, 0));
    }

    sort(delayVector.begin(), delayVector.end(), greater<int>());
    sort(jitterVector.begin(), jitterVector.end(), greater<int>());


    for (int i = 0; i < n - 2; i++)
        bigMDelay += delayVector[i], bigMJitter += jitterVector[i];

//    cout << "BigM Delay = " << bigMDelay << endl;
//    cout << "BigM Jitter = " << bigMJitter << endl;

//    cout << "Distance to terminals" << endl;
//    for (int i = 0; i < n; i++) cout << "k: " << i << ", distance: " << distance[i] << endl;
//    output << "graph ends\n";
    output.close();
    cout << "Load graph successfully" << endl;
}

void Graph::graphReduction() {
    int u, j, countEdges = 0, minSP, lbsi, prevDelay;
    vector<bool> alreadyVisited = vector<bool>(n);
    vector<int> jitterFromCShp = vector<int>(n), delayFromCShp = vector<int>(n), minSPVec = vector<int>(n);

    for (int i = 0; i < n; ++i) {
        if (!removed[i]) {
            dijkstra_shortest_paths(graphJitterSP, i, predecessor_map(
                    make_iterator_property_map(predecessors.begin(), get(vertex_index, graphJitterSP))).distance_map(
                    make_iterator_property_map(distanceJitter.begin(), get(vertex_index, graphJitterSP))));

            minSP = numeric_limits<int>::max();
            for (auto t : terminals)
                if (t != i && distanceJitter[t] < minSP)
                    minSP = distanceJitter[t];
            minSPVec[i] = minSP;

            add_vertex(SPPRC_Graph_Vert(i, paramJitter), graphDelay);
            add_vertex(SPPRC_Graph_Vert(i, paramDelay), graphJitter);
        }
    }

    for (u = 0; u < n; ++u) {
        if (!removed[u]) {
            for (auto arc : arcs[u]) {
                j = arc->getD();
                if (j == 0) continue;
                add_edge(u, j, SPPRC_Graph_Arc(countEdges, arc->getDelay(), arc->getJitter()), graphDelay);
                add_edge(u, j, SPPRC_Graph_Arc(countEdges++, arc->getJitter(), arc->getDelay()), graphJitter);
                add_edge(j, u, SPPRC_Graph_Arc(countEdges, arc->getDelay(), arc->getJitter()), graphDelay);
                add_edge(j, u, SPPRC_Graph_Arc(countEdges++, arc->getJitter(), arc->getDelay()), graphJitter);
            }
        }
    }

    vector<vector<graph_traits<SPPRCGraph>::edge_descriptor>> opt_solutions;
    vector<spp_spp_res_cont> pareto_opt;

    // CSHP root -> DuS
    for (auto j : DuS) {
        if (removed[j] || j == 0) continue;

        // Jitter
        r_c_shortest_paths
                (graphJitter,
                 get(&SPPRC_Graph_Vert::num, graphJitter),
                 get(&SPPRC_Graph_Arc::num, graphJitter),
                 root,
                 j,
                 opt_solutions,
                 pareto_opt,
                 spp_spp_res_cont(0, 0),
                 ref_spprc(),
                 dominance_spptw(),
                 allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                 default_r_c_shortest_paths_visitor());

        jitterFromCShp[j] = pareto_opt[0].cost;

        if (pareto_opt.empty() || jitterFromCShp[j] > paramJitter) notAttend[j] = true;

        // Delay
        SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];

        vert_prop.con = paramJitter - minSPVec[j];

        r_c_shortest_paths
                (graphDelay,
                 get(&SPPRC_Graph_Vert::num, graphDelay),
                 get(&SPPRC_Graph_Arc::num, graphDelay),
                 root,
                 j,
                 opt_solutions,
                 pareto_opt,
                 spp_spp_res_cont(0, 0),
                 ref_spprc(),
                 dominance_spptw(),
                 allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                 default_r_c_shortest_paths_visitor());


        delayFromCShp[j] = pareto_opt[0].cost;

        if (pareto_opt.empty() || delayFromCShp[j] > paramDelay) notAttend[j] = true;

        vert_prop.con = paramJitter;
    }

    // Moterated Vertex Elimination
    bool rem;
    for (auto i : nonTerminals) {
        if (i == 0) continue;
        if (!removed[i]) {
            rem = true;
            dijkstra_shortest_paths(graphJitterSP, i, predecessor_map(
                    make_iterator_property_map(predecessors.begin(), get(vertex_index, graphJitterSP))).distance_map(
                    make_iterator_property_map(distanceJitter.begin(), get(vertex_index, graphJitterSP))));

            for (auto t : terminals) {
                if (jitterFromCShp[i] + distanceJitter[t] <= paramJitter) {
                    rem = false;
                    break;
                }
            }
            if (rem) {
                cout << i << endl;
                notAttend[i] = true;
            }
        }
    }

    for (j = 1; j < n; ++j)
        if (notAttend[j]) removed[j] = true;

    for (j = 1; j < n; j++) {
        if (removed[j]) arcs[j].erase(arcs[j].begin(), arcs[j].end());
        for (int i = 0; i < int(arcs[j].size()); i++) {
            if (arcs[j][i]->getD() == 0) continue;
            if (removed[arcs[j][i]->getD()]) {
                arcs[j].erase(arcs[j].begin() + i);
                i--;
            }
        }
    }
    // =========== end ==============

    // Strong Arc Elimination

    // Terminals
    for (auto i : DuS) {
        if (removed[i] || i == 0) continue;
        lbsi = delayFromCShp[i];

        for (auto arc : arcs[i]) {
            j = arc->getD();

            SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
            vert_prop.con = paramJitter - jitterFromCShp[j];

            for (auto k : terminals) {
                if (k == i || k == j) continue;
                r_c_shortest_paths
                        (graphDelay,
                         get(&SPPRC_Graph_Vert::num, graphDelay),
                         get(&SPPRC_Graph_Arc::num, graphDelay),
                         k,
                         j,
                         opt_solutions,
                         pareto_opt,
                         spp_spp_res_cont(0, 0),
                         ref_spprc(),
                         dominance_spptw(),
                         allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                         default_r_c_shortest_paths_visitor());

                if (pareto_opt.empty() || lbsi + arc->getDelay() + pareto_opt[0].cost > paramDelay) {
//                    cout << "[" << i << ", " << j  << ", " << k << "]" << endl;
                    removedF[i][j][k] = true;
                }
            }
            vert_prop.con = paramJitter;
        }
    }

    getchar();

    /*
    for (int i = 0; i < n; ++i) {
        if (!removed[i]) {
            dijkstra_shortest_paths(graphJitterSP, i, predecessor_map(
                    make_iterator_property_map(predecessors.begin(), get(vertex_index, graphJitterSP))).distance_map(
                    make_iterator_property_map(distanceJitter.begin(), get(vertex_index, graphJitterSP))));

            minSP = numeric_limits<int>::max();
            for (auto t : terminals)
                if (t != i && distanceJitter[t] < minSP)
                    minSP = distanceJitter[t];
            minSPVec[i] = minSP;

            add_vertex(SPPRC_Graph_Vert(i, paramJitter), graphDelay);
            add_vertex(SPPRC_Graph_Vert(i, paramDelay), graphJitter);
        }
    }

    for (u = 0; u < n; ++u) {
        if (!removed[u]) {
            for (auto arc : arcs[u]) {
                j = arc->getD();
                if (j == 0) continue;
                add_edge(u, j, SPPRC_Graph_Arc(countEdges, arc->getDelay(), arc->getJitter()), graphDelay);
                add_edge(u, j, SPPRC_Graph_Arc(countEdges++, arc->getJitter(), arc->getDelay()), graphJitter);
                add_edge(j, u, SPPRC_Graph_Arc(countEdges, arc->getDelay(), arc->getJitter()), graphDelay);
                add_edge(j, u, SPPRC_Graph_Arc(countEdges++, arc->getJitter(), arc->getDelay()), graphJitter);
            }
        }
    }

    vector<vector<graph_traits<SPPRCGraph>::edge_descriptor>> opt_solutions;
    vector<spp_spp_res_cont> pareto_opt;

    // CSHP root -> nonterminals
    for (auto j : nonTerminals) {
        if (removed[j] || j == 0) continue;

        r_c_shortest_paths
                (graphJitter,
                 get(&SPPRC_Graph_Vert::num, graphJitter),
                 get(&SPPRC_Graph_Arc::num, graphJitter),
                 root,
                 j,
                 opt_solutions,
                 pareto_opt,
                 spp_spp_res_cont(0, 0),
                 ref_spprc(),
                 dominance_spptw(),
                 allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                 default_r_c_shortest_paths_visitor());

        jitterFromCShp[j] = pareto_opt[0].cost;

        if (pareto_opt.empty() || jitterFromCShp[j] > paramJitter) notAttend[j] = true;

        // Delay
        SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
        vert_prop.con = paramJitter - minSPVec[j];

        r_c_shortest_paths
                (graphDelay,
                 get(&SPPRC_Graph_Vert::num, graphDelay),
                 get(&SPPRC_Graph_Arc::num, graphDelay),
                 root,
                 j,
                 opt_solutions,
                 pareto_opt,
                 spp_spp_res_cont(0, 0),
                 ref_spprc(),
                 dominance_spptw(),
                 allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                 default_r_c_shortest_paths_visitor());


        delayFromCShp[j] = pareto_opt[0].cost;

        if (pareto_opt.empty() || delayFromCShp[j] > paramDelay) notAttend[j] = true;

        vert_prop.con = paramJitter;
    }

    // CSHP root -> terminals
    for (auto j : terminals) {
        r_c_shortest_paths
                (graphJitter,
                 get(&SPPRC_Graph_Vert::num, graphJitter),
                 get(&SPPRC_Graph_Arc::num, graphJitter),
                 root,
                 j,
                 opt_solutions,
                 pareto_opt,
                 spp_spp_res_cont(0, 0),
                 ref_spprc(),
                 dominance_spptw(),
                 allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                 default_r_c_shortest_paths_visitor());

        jitterFromCShp[j] = pareto_opt[0].cost;

        // Delay
        r_c_shortest_paths
                (graphDelay,
                 get(&SPPRC_Graph_Vert::num, graphDelay),
                 get(&SPPRC_Graph_Arc::num, graphDelay),
                 root,
                 j,
                 opt_solutions,
                 pareto_opt,
                 spp_spp_res_cont(0, 0),
                 ref_spprc(),
                 dominance_spptw(),
                 allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                 default_r_c_shortest_paths_visitor());

        delayFromCShp[j] = pareto_opt[0].cost;
    }

    // Terminals
    for (auto i : terminals) {
        lbsi = delayFromCShp[i];

        for (auto arc : arcs[i]) {
            j = arc->getD();

            SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
            vert_prop.con = paramJitter - jitterFromCShp[j];

            for (auto k : terminals) {
                if (k == i || k == j) continue;
                r_c_shortest_paths
                        (graphDelay,
                         get(&SPPRC_Graph_Vert::num, graphDelay),
                         get(&SPPRC_Graph_Arc::num, graphDelay),
                         k,
                         j,
                         opt_solutions,
                         pareto_opt,
                         spp_spp_res_cont(0, 0),
                         ref_spprc(),
                         dominance_spptw(),
                         allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                         default_r_c_shortest_paths_visitor());

                if (pareto_opt.empty() || lbsi + arc->getDelay() + pareto_opt[0].cost > paramDelay)
                    removedF[i][j][k] = true;
            }
            vert_prop.con = paramJitter;
        }
    }

    // Non-terminals
    for (auto i : nonTerminals) {
        if (removed[i] || i == 0) continue;

        lbsi = delayFromCShp[i];
        for (auto arc : arcs[i]) {
            j = arc->getD();
            if (j == root || j == 0) continue;

            SPPRC_Graph_Vert &vert_prop = get(vertex_bundle, graphDelay)[j];
            prevDelay = vert_prop.con;
            vert_prop.con = paramJitter - jitterFromCShp[j];

            for (auto k : terminals) {
                if (k == j) continue;

                r_c_shortest_paths
                        (graphDelay,
                         get(&SPPRC_Graph_Vert::num, graphDelay),
                         get(&SPPRC_Graph_Arc::num, graphDelay),
                         k,
                         j,
                         opt_solutions,
                         pareto_opt,
                         spp_spp_res_cont(0, 0),
                         ref_spprc(),
                         dominance_spptw(),
                         allocator<r_c_shortest_paths_label<SPPRCGraph, spp_spp_res_cont >>(),
                         default_r_c_shortest_paths_visitor());

                if (pareto_opt.empty() || lbsi + arc->getDelay() + pareto_opt[0].cost > paramDelay)
                    removedF[i][j][k] = true;
            }
            vert_prop.con = prevDelay;
        }
    }

//    for (j = 1; j < n; ++j)
//        if (notAttend[j]) removed[j] = true;

//    for (j = 1; j < n; j++) {
//        if (removed[j]) arcs[j].erase(arcs[j].begin(), arcs[j].end());
//        for (int i = 0; i < int(arcs[j].size()); i++) {
//            if (arcs[j][i]->getD() == 0) continue;
//            if (removed[arcs[j][i]->getD()]) {
//                arcs[j].erase(arcs[j].begin() + i);
//                i--;
//            }
//        }
//    }

    cout << "vertex to remove" << endl;
    for (j = 0; j < n; j++)
        if (removed[j]) {
            cout << j << ", ";
        }
    cout << endl;

//    for (j = 0; j < n; j++)
//        for (auto arc : arcs[j])
//            cout << j << " - " << arc->getD() << endl;

//    getchar();
 */
}

void Graph::showGraph() {
    cout << "Arcs" << endl;
    for (int o = 0; o < n; o++) {
        for (auto *arc : arcs[o])
            cout << o << " " << arc->getD() << ": " << arc->getDelay()
                 << " " << arc->getJitter() << " " << arc->getBandwidth() << " " <<
                 arc->getEstimateLinkDuration() << endl;

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
    return distanceDelay[k];
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

