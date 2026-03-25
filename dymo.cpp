#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include <queue>

using namespace std;

int update_count = 1;

struct PQEntry {
    int cand_size, neg_unmatched, u;
    bool operator>(const PQEntry& o) const {
#ifdef EXCLUDE_QUEUE
        return u > o.u;
#else
        if (cand_size != o.cand_size) return cand_size > o.cand_size;
        return neg_unmatched > o.neg_unmatched;
#endif
    }
};

struct HashPairII {
    size_t operator()(const pair<int, int>& p) const {
        return size_t(p.first) * 1000000007ULL + size_t(p.second);
    }
};
struct HashPairIIInt {
    size_t operator()(const pair<pair<int, int>, int>& p) const {
        return HashPairII()(p.first) * 1000000007ULL + size_t(p.second);
    }
};

int turns = 1;
struct Vertex {
    int label;
    unordered_set<int> neighbors;
};

vector<Vertex> G;
vector<Vertex> Q;
vector<pair<int, int>> updates;
int G_size, Q_size;
long long match_count;
double queue_time_ms = 0.0;

unordered_map<int, vector<int>> label_to_data_nodes;
unordered_map<int, vector<int>> label_to_query_nodes;

unordered_map<pair<int, int>, int, HashPairII> query_nlf;
unordered_map<pair<int, int>, int, HashPairII> data_nlf;
unordered_map<int, unordered_set<int>> nlf_view;
unordered_map<pair<pair<int, int>, int>, unordered_set<int>, HashPairIIInt> nlf_neighbor_view;
unordered_map<pair<int, int>, vector<pair<int, int>>, HashPairII> edge_inverted_index;
unordered_map<int, unordered_map<int, int>> query_nlf_by_u;

string g_path = "../dataset/dz_3/initial";
string q_path = "../dataset/dz_3/Q/8/q10";
string s_path = "../dataset/dz_3/s";
int max_updates = 0;

void inputQ() {
    ifstream qg(q_path);
    if (!qg) {
        cerr << "Failed to open query file: " << q_path << endl;
        exit(1);
    }

    char c;
    int id, id1, id2, label;
    Vertex u;

    while (qg >> c) {
        if (c == 'v') {
            qg >> id >> label;
            if (label == -1) label = 0;
            u.label = label;
            u.neighbors.clear();
            Q.push_back(u);
        } else {
            qg >> id1 >> id2;
            Q[id1].neighbors.insert(id2);
            Q[id2].neighbors.insert(id1);
            break;
        }
    }

    while (qg >> c >> id1 >> id2) {
        Q[id1].neighbors.insert(id2);
        Q[id2].neighbors.insert(id1);
    }

    qg.close();
    Q_size = Q.size();

    label_to_query_nodes.clear();
    for (int i = 0; i < Q_size; i++) {
        label_to_query_nodes[Q[i].label].push_back(i);
    }
}

void buildQueryNLF() {
    query_nlf.clear();
    query_nlf_by_u.clear();
    for (int u = 0; u < Q_size; u++) {
        for (int u_nei : Q[u].neighbors) {
            int L = Q[u_nei].label;
            query_nlf[make_pair(u, L)]++;
        }
    }
    for (const auto& pr : query_nlf) {
        int u = pr.first.first, L = pr.first.second, cnt = pr.second;
        query_nlf_by_u[u][L] = cnt;
    }
}

void buildDataNLF() {
    data_nlf.clear();
    unordered_set<pair<int, int>, HashPairII> valid_pairs;
    for (const auto& pr : edge_inverted_index) valid_pairs.insert(pr.first);
    for (int v = 0; v < G_size; v++) {
        int Lv = G[v].label;
        for (int v_nei : G[v].neighbors) {
            int L_nei = G[v_nei].label;
            int a = min(Lv, L_nei), b = max(Lv, L_nei);
            if (valid_pairs.count(make_pair(a, b))) data_nlf[make_pair(v, L_nei)]++;
        }
    }
}

bool satisfyNLF(int u, int v) {
    if (Q[u].label != G[v].label) return false;
    auto it_u = query_nlf_by_u.find(u);
    if (it_u == query_nlf_by_u.end()) return true;
    for (const auto& pr : it_u->second) {
        int L = pr.first, need = pr.second;
        int have = 0;
        auto it = data_nlf.find(make_pair(v, L));
        if (it != data_nlf.end()) have = it->second;
        if (have < need) return false;
    }
    return true;
}

void buildNlfView() {
    nlf_view.clear();
    for (int u = 0; u < Q_size; u++) {
        int Lu = Q[u].label;
        auto it = label_to_data_nodes.find(Lu);
        if (it == label_to_data_nodes.end()) continue;
        for (int v : it->second) {
            if (satisfyNLF(u, v)) nlf_view[u].insert(v);
        }
    }
}

void buildNlfNeighborView() {
    nlf_neighbor_view.clear();
    for (int u = 0; u < Q_size; u++) {
        for (int v : nlf_view[u]) {
            for (int u_nei : Q[u].neighbors) {
                int L = Q[u_nei].label;
                unordered_set<int> cand;
                for (int w : G[v].neighbors) {
                    if (G[w].label != L) continue;
                    if (!nlf_view[u_nei].count(w)) continue;
                    cand.insert(w);
                }
                if (!cand.empty()) nlf_neighbor_view[make_pair(make_pair(u, v), u_nei)] = cand;
            }
        }
    }
}

void buildEdgeInvertedIndex() {
    edge_inverted_index.clear();
    for (int u1 = 0; u1 < Q_size; u1++) {
        for (int u2 : Q[u1].neighbors) {
            if (u2 <= u1) continue;
            int L1 = Q[u1].label, L2 = Q[u2].label;
            int a = min(L1, L2), b = max(L1, L2);
            edge_inverted_index[make_pair(a, b)].push_back(make_pair(u1, u2));
        }
    }
}

void updateNlfOnEdgeAdd(int v1, int v2) {
    int L1 = G[v1].label, L2 = G[v2].label;
    int a = min(L1, L2), b = max(L1, L2);
    if (edge_inverted_index.find(make_pair(a, b)) == edge_inverted_index.end()) {
        return;
    }
    data_nlf[make_pair(v1, L2)]++;
    data_nlf[make_pair(v2, L1)]++;

    auto updateViewFor = [](int v) {
        int Lv = G[v].label;
        auto it = label_to_query_nodes.find(Lv);
        if (it == label_to_query_nodes.end()) return;
        for (int u : it->second) {
            if (satisfyNLF(u, v)) nlf_view[u].insert(v);
            else nlf_view[u].erase(v);
        }
    };
    updateViewFor(v1);
    updateViewFor(v2);

    for (int u : label_to_query_nodes[L1]) {
        if (!nlf_view[u].count(v1)) continue;
        for (int u_nei : Q[u].neighbors) {
            if (Q[u_nei].label != L2) continue;
            auto key1 = make_pair(make_pair(u, v1), u_nei);
            unordered_set<int> cand;
            for (int w : G[v1].neighbors) {
                if (G[w].label != L2) continue;
                if (!nlf_view[u_nei].count(w)) continue;
                cand.insert(w);
            }
            if (cand.empty()) nlf_neighbor_view.erase(key1);
            else nlf_neighbor_view[key1] = cand;
        }
    }
    for (int u : label_to_query_nodes[L2]) {
        if (!nlf_view[u].count(v2)) continue;
        for (int u_nei : Q[u].neighbors) {
            if (Q[u_nei].label != L1) continue;
            auto key2 = make_pair(make_pair(u, v2), u_nei);
            unordered_set<int> cand;
            for (int w : G[v2].neighbors) {
                if (G[w].label != L1) continue;
                if (!nlf_view[u_nei].count(w)) continue;
                cand.insert(w);
            }
            if (cand.empty()) nlf_neighbor_view.erase(key2);
            else nlf_neighbor_view[key2] = cand;
        }
    }
    auto updateNeighborViewForNeighborsOf = [](int v, int Lv) {
        for (int w : G[v].neighbors) {
            for (int u = 0; u < Q_size; u++) {
                if (!nlf_view[u].count(w)) continue;
                for (int u_nei : Q[u].neighbors) {
                    if (Q[u_nei].label != Lv) continue;
                    auto key = make_pair(make_pair(u, w), u_nei);
                    unordered_set<int> cand;
                    for (int x : G[w].neighbors) {
                        if (G[x].label != Lv) continue;
                        if (!nlf_view[u_nei].count(x)) continue;
                        cand.insert(x);
                    }
                    if (cand.empty()) nlf_neighbor_view.erase(key);
                    else nlf_neighbor_view[key] = cand;
                }
            }
        }
    };
    updateNeighborViewForNeighborsOf(v1, L1);
    updateNeighborViewForNeighborsOf(v2, L2);
}

void updateNlfOnEdgeDel(int v1, int v2) {
    int L1 = G[v1].label, L2 = G[v2].label;
    int a = min(L1, L2), b = max(L1, L2);
    if (edge_inverted_index.find(make_pair(a, b)) == edge_inverted_index.end()) {
        return;
    }
    data_nlf[make_pair(v1, L2)]--;
    data_nlf[make_pair(v2, L1)]--;
    if (data_nlf[make_pair(v1, L2)] == 0) data_nlf.erase(make_pair(v1, L2));
    if (data_nlf[make_pair(v2, L1)] == 0) data_nlf.erase(make_pair(v2, L1));

    auto updateViewFor = [](int v) {
        int Lv = G[v].label;
        auto it = label_to_query_nodes.find(Lv);
        if (it == label_to_query_nodes.end()) return;
        for (int u : it->second) {
            if (satisfyNLF(u, v)) nlf_view[u].insert(v);
            else nlf_view[u].erase(v);
        }
    };
    updateViewFor(v1);
    updateViewFor(v2);

    for (int u : label_to_query_nodes[L1]) {
        if (!nlf_view[u].count(v1)) continue;
        for (int u_nei : Q[u].neighbors) {
            if (Q[u_nei].label != L2) continue;
            auto key1 = make_pair(make_pair(u, v1), u_nei);
            unordered_set<int> cand;
            for (int w : G[v1].neighbors) {
                if (G[w].label != L2) continue;
                if (!nlf_view[u_nei].count(w)) continue;
                cand.insert(w);
            }
            if (cand.empty()) nlf_neighbor_view.erase(key1);
            else nlf_neighbor_view[key1] = cand;
        }
    }
    for (int u : label_to_query_nodes[L2]) {
        if (!nlf_view[u].count(v2)) continue;
        for (int u_nei : Q[u].neighbors) {
            if (Q[u_nei].label != L1) continue;
            auto key2 = make_pair(make_pair(u, v2), u_nei);
            unordered_set<int> cand;
            for (int w : G[v2].neighbors) {
                if (G[w].label != L1) continue;
                if (!nlf_view[u_nei].count(w)) continue;
                cand.insert(w);
            }
            if (cand.empty()) nlf_neighbor_view.erase(key2);
            else nlf_neighbor_view[key2] = cand;
        }
    }
    auto updateNeighborViewForNeighborsOfDel = [](int v, int Lv) {
        for (int w : G[v].neighbors) {
            for (int u = 0; u < Q_size; u++) {
                if (!nlf_view[u].count(w)) continue;
                for (int u_nei : Q[u].neighbors) {
                    if (Q[u_nei].label != Lv) continue;
                    auto key = make_pair(make_pair(u, w), u_nei);
                    unordered_set<int> cand;
                    for (int x : G[w].neighbors) {
                        if (G[x].label != Lv) continue;
                        if (!nlf_view[u_nei].count(x)) continue;
                        cand.insert(x);
                    }
                    if (cand.empty()) nlf_neighbor_view.erase(key);
                    else nlf_neighbor_view[key] = cand;
                }
            }
        }
    };
    updateNeighborViewForNeighborsOfDel(v1, L1);
    updateNeighborViewForNeighborsOfDel(v2, L2);
}

void inputG() {
    ifstream dg(g_path);
    if (!dg) {
        cerr << "Failed to open graph file: " << g_path << endl;
        exit(1);
    }

    char c;
    int id, id1, id2, label;
    Vertex v;

    label_to_data_nodes.clear();

    while (dg >> c) {
        if (c == 'v') {
            dg >> id >> label;
            if (label == -1) label = 0;
            v.label = label;
            v.neighbors.clear();
            G.push_back(v);
            
            label_to_data_nodes[label].push_back(id);
        } else {
            G_size = G.size();
            dg >> id1 >> id2;
            if (G[id1].neighbors.find(id2) != G[id1].neighbors.end()) {
                cerr << "Error: Graph contains duplicated neighbors!" << endl;
                continue;
            }
            if (id1 < G_size && id2 < G_size) {
                G[id1].neighbors.insert(id2);
                G[id2].neighbors.insert(id1);
            }
            break;
        }
    }

    while (dg >> c >> id1 >> id2) {
        if (id1 < G_size && id2 < G_size) {
            G[id1].neighbors.insert(id2);
            G[id2].neighbors.insert(id1);
        }
    }

    dg.close();
    if (G_size == 0) {
        G_size = G.size();
    }
}

void inputUpdates() {
    updates.clear();
    ifstream infile(s_path);
    if (!infile) {
        cerr << "Failed to open updates file: " << s_path << endl;
        exit(1);
    }

    char c;
    int v1, v2;
    int cnt = 0;

    while (infile >> c >> v1 >> v2) {
        if (max_updates != 0 && ++cnt > max_updates) {
            break;
        }
        updates.emplace_back(v1, v2);
    }

    infile.close();
}

static unordered_set<int> computeInitialCands(int u, const vector<int>& match,
                                              const unordered_set<int>& matched_u,
                                              const unordered_set<int>& matched_v) {
    vector<int> nei_in_matched;
    for (int i : matched_u)
        if (Q[i].neighbors.count(u)) nei_in_matched.push_back(i);
    unordered_set<int> cand;
    if (!nei_in_matched.empty()) {
        auto key0 = make_pair(make_pair(nei_in_matched[0], match[nei_in_matched[0]]), u);
        auto it0 = nlf_neighbor_view.find(key0);
        if (it0 != nlf_neighbor_view.end()) {
            cand = it0->second;
            for (size_t j = 1; j < nei_in_matched.size(); j++) {
                auto key = make_pair(make_pair(nei_in_matched[j], match[nei_in_matched[j]]), u);
                auto it = nlf_neighbor_view.find(key);
                if (it != nlf_neighbor_view.end()) {
                    unordered_set<int> inter;
                    for (int w : cand) if (it->second.count(w)) inter.insert(w);
                    cand = std::move(inter);
                }
                if (cand.empty()) break;
            }
        }
    }
    for (int w : matched_v) cand.erase(w);
    return cand;
}

static int countUnmatchedNeighbors(int u, const unordered_set<int>& matched_u) {
    int c = 0;
    for (int w : Q[u].neighbors) if (!matched_u.count(w)) c++;
    return c;
}

void expandBoundary(const unordered_map<int, unordered_set<int>>& u_cands,
                   unordered_set<int>& boundary,
                   vector<int>& match,
                   unordered_set<int>& matched_v) {
    if (boundary.empty()) {
        match_count++;
        return;
    }
    
    for (int u : boundary) {
        unordered_set<int> remaining = boundary;
        remaining.erase(u);
        
        for (int v : u_cands.at(u)) {
            if (matched_v.count(v)) continue;
            match[u] = v;
            matched_v.insert(v);
            
            expandBoundary(u_cands, remaining, match, matched_v);
            
            matched_v.erase(v);
            match[u] = -1;
        }
        return;
    }
}

void expandMatch(vector<int>& match, unordered_set<int>& matched_u, unordered_set<int>& matched_v,
                 unordered_map<int, unordered_set<int>>& u_cands,
                 unordered_set<int>& boundary,
                 priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>>& pq) {
    if (pq.empty()) return;
    int u;
    while (true) {
        if (pq.empty()) return;
        PQEntry e = pq.top();
        pq.pop();
        u = e.u;
        if (matched_u.count(u)) continue;
        if (!u_cands.count(u) || u_cands[u].empty()) return;
        break;
    }
    unordered_set<int> cands_u = u_cands[u];
    for (int v : cands_u) {
            if (matched_v.count(v)) continue;
            priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>> pq_copy = pq;
            match[u] = v;
            matched_u.insert(u);
            matched_v.insert(v);

            if (matched_u.size() == (size_t)Q_size) {
                match_count++;
            } else {
                vector<int> update_u;
                unordered_map<int, vector<int>> deleted_cands;
                bool abort_branch = false;

                auto uv_key = make_pair(u, v);
                for (int u_nei : Q[u].neighbors) {
                    if (matched_u.count(u_nei)) continue;
                    bool we_modified = false;
                    auto key_unei = make_pair(uv_key, u_nei);
                    if (!u_cands.count(u_nei)) {
                        auto it_cand = nlf_neighbor_view.find(key_unei);
                        unordered_set<int> cand = (it_cand != nlf_neighbor_view.end()) ? it_cand->second : unordered_set<int>{};
                        vector<int> deleted;
                        for (int w : cand) {
                            if (matched_v.count(w)) deleted.push_back(w);
                        }
                        for (int w : deleted) cand.erase(w);
                        if (cand.empty()) {
                            abort_branch = true;
                            break;
                        }
                        u_cands[u_nei] = std::move(cand);
                        we_modified = true;
                    } else {
                        auto it = nlf_neighbor_view.find(key_unei);
                        if (it == nlf_neighbor_view.end()) {
                            abort_branch = true;
                            break;
                        }
                        const unordered_set<int>& new_cands_ref = it->second;
                        unordered_set<int> u_nei_cands_copy = u_cands[u_nei];
                        for (int w : u_nei_cands_copy) {
                            if (matched_v.count(w) || !new_cands_ref.count(w)) {
                                deleted_cands[u_nei].push_back(w);
                            }
                        }
                        if (!deleted_cands[u_nei].empty()) {
                            for (int w : deleted_cands[u_nei])
                                u_cands[u_nei].erase(w);
                            we_modified = true;
                        }
                    }
                    if (u_cands[u_nei].empty()) {
                        if (we_modified) update_u.push_back(u_nei);
                        abort_branch = true;
                        break;
                    }
                    if (we_modified) update_u.push_back(u_nei);
                    int neg = -countUnmatchedNeighbors(u_nei, matched_u);
                    pq_copy.push(PQEntry{(int)u_cands[u_nei].size(), neg, u_nei});
                }
                if (!abort_branch){
                    if (matched_u.size() + boundary.size() == (size_t)Q_size) {
                        expandBoundary(u_cands, boundary, match, matched_v);
                    } else {
                        expandMatch(match, matched_u, matched_v, u_cands, boundary, pq_copy);
                    }
                } 

                for (int u_nei : update_u) {
                    if (deleted_cands.count(u_nei) && !deleted_cands[u_nei].empty()) {
                        for (int w : deleted_cands[u_nei]){
                            u_cands[u_nei].insert(w);
                        } 
                    } else {
                        u_cands.erase(u_nei);
                    }
                }
            }

            matched_v.erase(v);
            matched_u.erase(u);
            match[u] = -1;
        }
}

static bool buildInitialCandsAndQueue(int u1, int u2, int v1, int v2,
                                      vector<int>& match,
                                      unordered_set<int>& matched_u, unordered_set<int>& matched_v,
                                      unordered_map<int, unordered_set<int>>& u_cands,
                                      unordered_set<int>& boundary,
                                      priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>>& pq) {
    unordered_set<int> expandable;
    for (int u : Q[u1].neighbors) if (u != u2) expandable.insert(u);
    for (int u : Q[u2].neighbors) if (u != u1) expandable.insert(u);

    for (int u : expandable) {
        unordered_set<int> cand = computeInitialCands(u, match, matched_u, matched_v);
        u_cands[u] = std::move(cand);
        if (u_cands[u].empty()) return false;
        int neg = -countUnmatchedNeighbors(u, matched_u);
        pq.push(PQEntry{(int)u_cands[u].size(), neg, u});
    }
    return true;
}

void searchMatchesFromUpdate(int start_v1, int start_v2) {
    match_count = 0;
    queue_time_ms = 0.0;
    int L1 = G[start_v1].label, L2 = G[start_v2].label;
    int a = min(L1, L2), b = max(L1, L2);
    auto it_eii = edge_inverted_index.find(make_pair(a, b));
    if (it_eii == edge_inverted_index.end()) return;

    vector<int> match(Q_size, -1);
    unordered_set<int> matched_u, matched_v;
    unordered_map<int, unordered_set<int>> u_cands;
    unordered_set<int> boundary;
    priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>> pq;
    bool flag = false;
    for (const auto& e : it_eii->second) {
        int u1 = e.first, u2 = e.second;
        auto k1a = make_pair(make_pair(u1, start_v1), u2), k1b = make_pair(make_pair(u2, start_v2), u1);
        auto k2a = make_pair(make_pair(u1, start_v2), u2), k2b = make_pair(make_pair(u2, start_v1), u1);
        auto it1a = nlf_neighbor_view.find(k1a), it1b = nlf_neighbor_view.find(k1b);
        auto it2a = nlf_neighbor_view.find(k2a), it2b = nlf_neighbor_view.find(k2b);
        bool ok1 = (it1a != nlf_neighbor_view.end() && it1a->second.count(start_v2) &&
                    it1b != nlf_neighbor_view.end() && it1b->second.count(start_v1));
        bool ok2 = (it2a != nlf_neighbor_view.end() && it2a->second.count(start_v1) &&
                    it2b != nlf_neighbor_view.end() && it2b->second.count(start_v2));

        if (Q[u1].label == L1 && Q[u2].label == L2 && ok1) {
            flag = true;
            fill(match.begin(), match.end(), -1);
            matched_u.clear();
            matched_v.clear();
            u_cands.clear();
            boundary.clear();
            pq = priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>>();
            match[u1] = start_v1;
            match[u2] = start_v2;
            matched_u.insert(u1);
            matched_u.insert(u2);
            matched_v.insert(start_v1);
            matched_v.insert(start_v2);
            if (!buildInitialCandsAndQueue(u1, u2, start_v1, start_v2, match, matched_u, matched_v, u_cands, boundary, pq)) continue;
            expandMatch(match, matched_u, matched_v, u_cands, boundary, pq);
        }
        if (Q[u1].label == L2 && Q[u2].label == L1 && ok2) {
            flag = true;
            fill(match.begin(), match.end(), -1);
            matched_u.clear();
            matched_v.clear();
            u_cands.clear();
            boundary.clear();
            pq = priority_queue<PQEntry, vector<PQEntry>, greater<PQEntry>>();
            match[u1] = start_v2;
            match[u2] = start_v1;
            matched_u.insert(u1);
            matched_u.insert(u2);
            matched_v.insert(start_v2);
            matched_v.insert(start_v1);
            if (!buildInitialCandsAndQueue(u1, u2, start_v2, start_v1, match, matched_u, matched_v, u_cands, boundary, pq)) continue;
            expandMatch(match, matched_u, matched_v, u_cands, boundary, pq);
        }
    }
    if (flag) turns++;
}

void updateAndMatching() {
    long long add_matches = 0;
    long long del_matches = 0;
    double update_time = 0.0;
    double search_time = 0.0;
    clock_t time;

    for (const auto& update : updates) {
        int v1 = update.first;
        int v2 = update.second;
        int current_matches = 0;

        if (v1 < 0) {
            v1 = -v1 - 1;
            v2 = -v2 - 1;

            if (v1 >= G_size || v2 >= G_size) {
                update_count++;
                continue;
            }

            if (!G[v1].neighbors.count(v2)) {
                update_count++;
                continue;
            }

            time = clock();
            searchMatchesFromUpdate(v1, v2);
            current_matches = match_count;
            del_matches += current_matches;
            search_time += double(clock() - time) * 1000 / CLOCKS_PER_SEC;

            time = clock();
            G[v1].neighbors.erase(v2);
            G[v2].neighbors.erase(v1);
            updateNlfOnEdgeDel(v1, v2);
            update_time += double(clock() - time) * 1000 / CLOCKS_PER_SEC;
        } else {
            if (v1 >= G_size || v2 >= G_size) {
                update_count++;
                continue;
            }

            if (G[v1].neighbors.count(v2)) {
                update_count++;
                continue;
            }

            time = clock();
            G[v1].neighbors.insert(v2);
            G[v2].neighbors.insert(v1);
            updateNlfOnEdgeAdd(v1, v2);
            update_time += double(clock() - time) * 1000 / CLOCKS_PER_SEC;

            time = clock();
            searchMatchesFromUpdate(v1, v2);
            current_matches = match_count;
            add_matches += current_matches;
            search_time += double(clock() - time) * 1000 / CLOCKS_PER_SEC;

        }

        update_count++;
    }

    cerr << "added matches: " << add_matches << endl;
    cerr << "deleted matches: " << del_matches << endl;
    cerr << "updated matches: " << add_matches + del_matches << endl;
    cerr << "update time: " << update_time << "ms" << endl;
    cerr << "search time: " << search_time << "ms" << endl;
}

int main(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        if (string(argv[i]) == "-d") {
            g_path = argv[i + 1];
        } else if (string(argv[i]) == "-q") {
            q_path = argv[i + 1];
        } else if (string(argv[i]) == "-c") {
            max_updates = atoi(argv[i + 1]);
        } else if (string(argv[i]) == "-s") {
            s_path = argv[i + 1];
        }
    }

    clock_t total_start = clock();
    clock_t time = clock();

    cerr << "========== Start inputting. ==========" << endl;
    clock_t t_q = clock();
    inputQ();
    cerr << "  inputQ cost (ms): " << double(clock() - t_q) * 1000 / CLOCKS_PER_SEC << endl;
    clock_t t_g = clock();
    inputG();
    cerr << "  inputG cost (ms): " << double(clock() - t_g) * 1000 / CLOCKS_PER_SEC << endl;
    clock_t t_s = clock();
    inputUpdates();
    cerr << "  inputUpdates cost (ms): " << double(clock() - t_s) * 1000 / CLOCKS_PER_SEC << endl;
    cerr << "Inputting cost (ms): " << double(clock() - time) * 1000 / CLOCKS_PER_SEC << endl;

    cerr << "========== Start NLF index building. ==========" << endl;
    clock_t nlf_start = clock();
    clock_t t0 = clock();
    buildQueryNLF();
    cerr << "  buildQueryNLF cost (ms): " << double(clock() - t0) * 1000 / CLOCKS_PER_SEC << endl;
    t0 = clock();
    buildEdgeInvertedIndex();
    cerr << "  buildEdgeInvertedIndex cost (ms): " << double(clock() - t0) * 1000 / CLOCKS_PER_SEC << endl;
    t0 = clock();
    buildDataNLF();
    cerr << "  buildDataNLF cost (ms): " << double(clock() - t0) * 1000 / CLOCKS_PER_SEC << endl;
    t0 = clock();
    buildNlfView();
    cerr << "  buildNlfView cost (ms): " << double(clock() - t0) * 1000 / CLOCKS_PER_SEC << endl;
    t0 = clock();
    buildNlfNeighborView();
    cerr << "  buildNlfNeighborView cost (ms): " << double(clock() - t0) * 1000 / CLOCKS_PER_SEC << endl;
    cerr << "NLF index building total cost (ms): " << double(clock() - nlf_start) * 1000 / CLOCKS_PER_SEC << endl;
    cerr << "========== End NLF index building. ==========" << endl;
    cerr << "The number of vertices in Q: " << Q_size << endl;
    cerr << "The number of vertices in G: " << G_size << endl;
    cerr << "========== End inputting. ==========" << endl;

    time = clock();
    cerr << "========== Start updating. ==========" << endl;
    updateAndMatching();
    cerr << "Updating totally cost (ms): " << double(clock() - time) * 1000 / CLOCKS_PER_SEC << endl;
    cerr << "========== End updating. ==========" << endl;

    cerr << "Total cost (ms): " << double(clock() - total_start) * 1000 / CLOCKS_PER_SEC << endl;
    cerr << "========== DONE!! ==========" << endl;
    cerr << "The number of turns: " << turns << endl;

    return 0;
}
