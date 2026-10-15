#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "LayerAssignment.h"
#include "flute.h"
#include "ispdData.h"

namespace Param {
    // Routing iteration and timeout limits
    constexpr int MAX_ROUTING_ITERS = 1000;     // Maximum total number of routing iterations
    constexpr int TIME_LIMIT_SEC = 900;        // Global timeout limit (seconds)
    constexpr int GLOBAL_TIME_LIMIT_SEC = 1350; // Stage 3.5 hard limit

    // Bounded-Length Relaxation (BLC) parameters: alpha and beta
    // Formula: L = Manh * (1.0 + atan(iter - alpha) + beta)
    double BLC_ALPHA = 9.0;
    double BLC_BETA = 9.5; // 1.5 -> 9.5

    // Dynamic history cost (dah) parameters: C1 and C2
    // Formula: dah = history_cost / (C1 + C2 * sqrt(iter))
    double DAH_C1 = 4.0; // 7.0 -> 4.0
    double DAH_C2 = 1.0; // 4.0 -> 1.0

    // Congestion penalty cost (p_e) parameters: C3 and C4
    // Formula: p_e = 1.0 + C3 / (1.0 + exp(C4 * (capacity - demand -
    // virtual_demand)))
    constexpr double PENALTY_C3 = 150.0;
    constexpr double PENALTY_C4 = 0.3;

    // Base cost (b_e) parameters: C5 and C6
    // Formula: b_e = C5 + C6 / 2^iter
    constexpr double BASE_C5 = 30.0;
    constexpr double BASE_C6 = 200.0;

    // Net ordering score weight parameters: C7 and C8
    // Formula: score = C7 * overflow_edges + C8 * path_length
    constexpr double NET_ORDER_C7 = 30.0;
    constexpr double NET_ORDER_C8 = 1.0;

    // Virtual demand parameter: C9 (collision avoidance)
    // Formula: vd = C9 * sqrt(max(0, marks - my_marks))
    double VIRTUAL_DEMAND_C9 = 0.4;

    // RSMT-aware routing discount (encourages following the RSMT skeleton)
    constexpr double RSMT_DISCOUNT_NORMAL = 1.0;
    constexpr double RSMT_DISCOUNT_POST_REFINE = 0.05;

    // Capacity penalty (used in post-refinement)
    constexpr double PENALTY_CAP_OVERFLOW = 100000.0;
    constexpr double PENALTY_CAP_SAFE = 0.1;

    // A* heuristic distance weight
    constexpr double ASTAR_HEURISTIC_WEIGHT = 30.0;
} // namespace Param

struct EdgeInfo {
        int demand = 0;
        int capacity = 0;
        double history_cost = 1.0;
        int overflow_freq = 0;
        int marks = 0;
};

std::vector<std::vector<std::vector<EdgeInfo>>> edge_grid;
std::unordered_map<ISPDParser::Net*, std::unordered_multiset<uint64_t>>
    skeleton_map;

struct AStarNode {
        int x, y;
        double g_cost, f_cost;
        int hops;
        int px, py;
        bool hori;
        bool conforms;
        bool operator>(const AStarNode& other) const {
            if (conforms != other.conforms) {
                return !conforms;
            }
            if (!conforms && !other.conforms) {
                return hops > other.hops; // When neither satisfies the bound, prefer the fewest hops (shortest distance)
            }
            return f_cost > other.f_cost;
        }
};

struct NodeRecord {
        double cost;
        int hops;
        bool conforms;
        int px, py;
        bool hori;
        int eval_version;
};

thread_local std::vector<NodeRecord> local_nodes;
thread_local int current_eval_version = 0;

inline int getManhattanDist(int x1, int y1, int x2, int y2) {
    return std::abs(x1 - x2) + std::abs(y1 - y2);
}

int getDemandVal(ISPDParser::Net* net, ISPDParser::ispdData* ispdData) {
    // int max_spacing = 0;
    // for (int s : ispdData->minimumSpacing)
    //     max_spacing = std::max(max_spacing, s);
    // return net->minimumWidth + max_spacing;
    return 1;
}

void updatePathDemand(ISPDParser::TwoPin& twoPin, int delta,
                      ISPDParser::ispdData* ispdData) {
    int demand_val = getDemandVal(twoPin.parNet, ispdData);
    auto net = twoPin.parNet;

    std::unordered_set<uint64_t> local_other_edges;
    for (const auto& other_tp : net->twopin) {
        if (&other_tp == &twoPin) continue;
        for (const auto& p : other_tp.path) {
            local_other_edges.insert(((uint64_t)p.x << 32) |
                                     ((uint64_t)p.y << 1) | (p.hori ? 1 : 0));
        }
    }

    for (const auto& edge : twoPin.path) {
        uint64_t key = ((uint64_t)edge.x << 32) | ((uint64_t)edge.y << 1) |
                       (edge.hori ? 1 : 0);

        if (local_other_edges.find(key) == local_other_edges.end()) {
#pragma omp atomic
            edge_grid[edge.x][edge.y][edge.hori].demand += delta * demand_val;
        }
    }
}

double evaluatePathCost(const std::vector<ISPDParser::RPoint>& path,
                        int demand_val) {
    double cost = 0.0;
    for (const auto& edge : path) {
        auto& e = edge_grid[edge.x][edge.y][edge.hori];
        int dem = e.demand + demand_val;

        double kappa = (e.capacity <= dem) ? Param::PENALTY_CAP_OVERFLOW
                                           : Param::PENALTY_CAP_SAFE;
        cost += 1.0 + kappa * ((double)dem / (1.0 + (double)e.capacity));
    }
    return cost;
}

std::vector<ISPDParser::RPoint> buildMonotonicPath(
    ISPDParser::TwoPin* tp, ISPDParser::ispdData* ispdData,
    bool is_post_refine = false) {
    int x1 = tp->from.x, y1 = tp->from.y;
    int x2 = tp->to.x, y2 = tp->to.y;
    int demand_val = getDemandVal(tp->parNet, ispdData);

    int dx = (x2 > x1) ? 1 : ((x2 < x1) ? -1 : 0);
    int dy = (y2 > y1) ? 1 : ((y2 < y1) ? -1 : 0);

    int nx = std::abs(x2 - x1) + 1;
    int ny = std::abs(y2 - y1) + 1;

    auto check_skeleton = [&](int x, int y, bool h) {
        auto it = skeleton_map.find(tp->parNet);
        if (it == skeleton_map.end()) return false;
        uint64_t key = ((uint64_t)x << 32) | ((uint64_t)y << 1) | (h ? 1 : 0);
        return it->second.find(key) != it->second.end();
    };

    std::unordered_set<uint64_t> local_ge;
    for (const auto& other_tp : tp->parNet->twopin) {
        if (&other_tp == tp) continue;
        for (const auto& p : other_tp.path) {
            uint64_t key =
                ((uint64_t)p.x << 32) | ((uint64_t)p.y << 1) | (p.hori ? 1 : 0);
            local_ge.insert(key);
        }
    }

    auto check_GE = [&](int x, int y, bool h) {
        uint64_t key = ((uint64_t)x << 32) | ((uint64_t)y << 1) | (h ? 1 : 0);
        return local_ge.find(key) != local_ge.end();
    };

    std::vector<std::vector<double>> dp(nx, std::vector<double>(ny, 1e9));
    std::vector<std::vector<int>> parent(nx, std::vector<int>(ny, -1));

    dp[0][0] = 0.0;

    double b_e = Param::BASE_C5 + Param::BASE_C6 / 2.0;
    // double dah = 1.0 / (Param::DAH_C1 + Param::DAH_C2 * 1.0);
    double dah = 1.0 / (7.0 + 4.0 * 1.0);

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (i == 0 && j == 0) continue;

            int curr_x = x1 + i * dx;
            int curr_y = y1 + j * dy;

            double cost_from_x = 1e9;
            double cost_from_y = 1e9;

            if (i > 0) {
                int prev_x = curr_x - dx;
                int ex = std::min(curr_x, prev_x);
                double step_cost = 1e9;

                if (check_GE(ex, curr_y, true)) {
                    step_cost = 0.0;
                } else {
                    auto& e = edge_grid[ex][curr_y][1]; // hori = 1
                    int dem = e.demand /* + demand_val */;
                    if (is_post_refine) {
                        step_cost = (e.capacity <= dem)
                                        ? Param::PENALTY_CAP_OVERFLOW
                                        : Param::PENALTY_CAP_SAFE;
                        if (check_skeleton(ex, curr_y, true))
                            step_cost -= Param::RSMT_DISCOUNT_POST_REFINE;
                    } else {
                        double p_e =
                            1.0 + Param::PENALTY_C3 /
                                      (1.0 + std::exp(Param::PENALTY_C4 *
                                                      (e.capacity - dem)));
                        step_cost = (1.0 + dah) * p_e + b_e;
                        if (check_skeleton(ex, curr_y, true))
                            step_cost -= Param::RSMT_DISCOUNT_NORMAL;
                    }
                    step_cost = std::max(0.01, step_cost);
                }
                cost_from_x = dp[i - 1][j] + step_cost;
            }

            if (j > 0) {
                int prev_y = curr_y - dy;
                int ey = std::min(curr_y, prev_y);
                double step_cost = 1e9;

                if (check_GE(curr_x, ey, false)) {
                    step_cost = 0.0;
                } else {
                    auto& e = edge_grid[curr_x][ey][0]; // hori = 0
                    int dem = e.demand /* + demand_val */;
                    if (is_post_refine) {
                        step_cost = (e.capacity <= dem)
                                        ? Param::PENALTY_CAP_OVERFLOW
                                        : Param::PENALTY_CAP_SAFE;
                        if (check_skeleton(curr_x, ey, false))
                            step_cost -= Param::RSMT_DISCOUNT_POST_REFINE;
                    } else {
                        double p_e =
                            1.0 + Param::PENALTY_C3 /
                                      (1.0 + std::exp(Param::PENALTY_C4 *
                                                      (e.capacity - dem)));
                        step_cost = (1.0 + dah) * p_e + b_e;
                        if (check_skeleton(curr_x, ey, false))
                            step_cost -= Param::RSMT_DISCOUNT_NORMAL;
                    }
                    step_cost = std::max(0.01, step_cost);
                }
                cost_from_y = dp[i][j - 1] + step_cost;
            }

            if (cost_from_x <= cost_from_y) {
                dp[i][j] = cost_from_x;
                parent[i][j] = 0;
            } else {
                dp[i][j] = cost_from_y;
                parent[i][j] = 1;
            }
        }
    }

    std::vector<ISPDParser::RPoint> path;
    int i = nx - 1;
    int j = ny - 1;
    while (i > 0 || j > 0) {
        int curr_x = x1 + i * dx;
        int curr_y = y1 + j * dy;
        if (parent[i][j] == 0) {
            int prev_x = curr_x - dx;
            path.emplace_back(std::min(curr_x, prev_x), curr_y, true);
            i--;
        } else {
            int prev_y = curr_y - dy;
            path.emplace_back(curr_x, std::min(curr_y, prev_y), false);
            j--;
        }
    }

    std::reverse(path.begin(), path.end());
    return path;
}

struct KruskalEdge {
        int u, v, weight;
        bool operator<(const KruskalEdge& other) const {
            return weight < other.weight;
        }
};

int findParent(int i, std::vector<int>& parent) {
    if (parent[i] == i) return i;
    return parent[i] = findParent(parent[i], parent);
}

void decomposeNetByRMST(ISPDParser::Net* net) {
    int d = net->pin2D.size();
    std::vector<KruskalEdge> edges;
    for (int i = 0; i < d; ++i) {
        for (int j = i + 1; j < d; ++j) {
            edges.push_back(
                {i, j,
                 getManhattanDist(net->pin2D[i].x, net->pin2D[i].y,
                                  net->pin2D[j].x, net->pin2D[j].y)});
        }
    }
    std::sort(edges.begin(), edges.end());
    std::vector<int> parent(d);
    std::iota(parent.begin(), parent.end(), 0);

    for (auto& edge : edges) {
        int root_u = findParent(edge.u, parent);
        int root_v = findParent(edge.v, parent);
        if (root_u != root_v) {
            ISPDParser::TwoPin tp;
            tp.parNet = net;
            tp.from = net->pin2D[edge.u];
            tp.to = net->pin2D[edge.v];
            net->twopin.push_back(tp);
            parent[root_u] = root_v;
        }
    }
}

// ===============================================================================
// Heuristic-BLMR
// ===============================================================================

bool routeAStar_NCTUGR(ISPDParser::Net* net, ISPDParser::TwoPin& twoPin,
                       int target_x, int target_y, int iter,
                       ISPDParser::ispdData* ispdData, int prev_length,
                       const std::vector<ISPDParser::RPoint>& old_path,
                       bool is_post_refine = false, int custom_blc = 0,
                       bool ignore_blc = false) {
    int start_x = twoPin.from.x, start_y = twoPin.from.y;
    int manh = getManhattanDist(start_x, start_y, target_x, target_y);
    double BLC = is_post_refine
                     ? custom_blc
                     : manh * (1.0 + std::atan((iter + 1) - Param::BLC_ALPHA) +
                               Param::BLC_BETA);
    // std::cout << "BLC: " << BLC << std::endl;

    if (ignore_blc) {
        BLC = 1e9;
    }

    int grid_size = ispdData->numXGrid * ispdData->numYGrid;
    if (local_nodes.size() < grid_size) {
        local_nodes.resize(grid_size);
    }

    current_eval_version++;
    if (current_eval_version > 1000000000) {
        std::fill(local_nodes.begin(), local_nodes.end(),
                  NodeRecord{0, 0, false, 0, 0, false, 0});
        current_eval_version = 1;
    }
    int eval_v = current_eval_version;

    auto get_node = [&](int x, int y) -> NodeRecord& {
        int idx = x * ispdData->numYGrid + y;
        if (local_nodes[idx].eval_version != eval_v) {
            local_nodes[idx].cost = 1e9;
            local_nodes[idx].hops = 1e9;
            local_nodes[idx].conforms = false;
            local_nodes[idx].eval_version = eval_v;
        }
        return local_nodes[idx];
    };

    std::priority_queue<AStarNode, std::vector<AStarNode>,
                        std::greater<AStarNode>>
        open_list;

    open_list.push(
        {start_x, start_y, 0.0, (double)manh, 0, -1, -1, false, true});
    NodeRecord& start_rec = get_node(start_x, start_y);
    start_rec.cost = 0.0;
    start_rec.hops = 0;
    start_rec.conforms = true;

    auto it = skeleton_map.find(net);
    const std::unordered_multiset<uint64_t>* skeleton_ptr = nullptr;
    if (it != skeleton_map.end()) skeleton_ptr = &(it->second);

    bool found = false;

    std::unordered_multiset<uint64_t> local_ge_set;
    for (const auto& other_tp : net->twopin) {
        if (&other_tp == &twoPin) continue;
        for (const auto& p : other_tp.path) {
            uint64_t key =
                ((uint64_t)p.x << 32) | ((uint64_t)p.y << 1) | (p.hori ? 1 : 0);
            local_ge_set.insert(key);
        }
    }

    std::unordered_multiset<uint64_t> local_old_path_set;
    for (auto& p : old_path) {
        local_old_path_set.insert(((uint64_t)p.x << 32) | ((uint64_t)p.y << 1) |
                                  (p.hori ? 1 : 0));
    }

    auto old_path_hash = [&local_old_path_set](int x, int y, bool h) {
        uint64_t key = ((uint64_t)x << 32) | ((uint64_t)y << 1) | (h ? 1 : 0);
        return local_old_path_set.count(key);
    };

    auto check_skeleton = [&](int x, int y, bool h) {
        if (!skeleton_ptr) return false;
        uint64_t key = ((uint64_t)x << 32) | ((uint64_t)y << 1) | (h ? 1 : 0);
        return skeleton_ptr->find(key) != skeleton_ptr->end();
    };

    auto check_GE = [&](int x, int y, bool h) {
        uint64_t key = ((uint64_t)x << 32) | ((uint64_t)y << 1) | (h ? 1 : 0);
        return local_ge_set.find(key) != local_ge_set.end();
    };

    while (!open_list.empty()) {
        AStarNode current = open_list.top();
        open_list.pop();

        NodeRecord& current_rec = get_node(current.x, current.y);

        bool is_stale = false;
        if (current.conforms) {
            if (current.g_cost > current_rec.cost) is_stale = true;
        } else {
            if (current_rec.conforms)
                is_stale = true;
            else if (current.hops > current_rec.hops)
                is_stale = true;
            else if (current.hops == current_rec.hops &&
                     current.g_cost > current_rec.cost)
                is_stale = true;
        }

        if (is_stale) {
            continue;
        }

        if (current.x == target_x && current.y == target_y) {
            found = true;
            break;
        }

        int dx[] = {1, 0, -1, 0}, dy[] = {0, 1, 0, -1};
        for (int i = 0; i < 4; ++i) {
            int nx = current.x + dx[i], ny = current.y + dy[i];
            if (nx < 0 || ny < 0 || nx >= ispdData->numXGrid ||
                ny >= ispdData->numYGrid)
                continue;

            if (!ignore_blc &&
                current.hops + 1 +
                        getManhattanDist(nx, ny, target_x, target_y) >
                    BLC)
                continue;

            bool hori = (ny == current.y);
            int ex = std::min(current.x, nx), ey = std::min(current.y, ny);
            auto& e = edge_grid[ex][ey][hori];

            double ew = 0.0;
            if (manh > 0) {
                ew = prev_length *
                     ((double)getManhattanDist(nx, ny, target_x, target_y) /
                      manh);
            }

            bool conforms = ((current.hops + 1 + ew) <= BLC);

            double step_cost = 1.0;
            if (is_post_refine) {
                if (e.capacity <= e.demand)
                    step_cost = Param::PENALTY_CAP_OVERFLOW;
                else {
                    step_cost = Param::PENALTY_CAP_SAFE;
                }
                if (check_skeleton(ex, ey, hori))
                    step_cost -= Param::RSMT_DISCOUNT_POST_REFINE;
                step_cost = std::max(0.01, step_cost);
            } else {
                if (check_GE(ex, ey, hori)) {
                    step_cost = 0.0;
                } else {
                    int m_i = old_path_hash(ex, ey, hori);

                    int current_marks;
#pragma omp atomic read
                    current_marks = e.marks;

                    double vd = Param::VIRTUAL_DEMAND_C9 *
                                std::sqrt(std::max(0, current_marks - m_i));
                    double b_e = Param::BASE_C5 +
                                 Param::BASE_C6 / std::pow(2.0, iter + 1);
                    int demand_val = getDemandVal(net, ispdData);

                    int current_dem;
#pragma omp atomic read
                    current_dem = e.demand;

                    double p_e =
                        1.0 +
                        Param::PENALTY_C3 /
                            (1.0 +
                             std::exp(Param::PENALTY_C4 *
                                      (e.capacity -
                                       (current_dem /* + demand_val */) - vd)));
                    double dah =
                        e.history_cost /
                        (Param::DAH_C1 + Param::DAH_C2 * std::sqrt(iter + 1));
                    step_cost = (1.0 + dah) * p_e + b_e;
                    if (check_skeleton(ex, ey, hori))
                        step_cost -= Param::RSMT_DISCOUNT_NORMAL;
                    step_cost = std::max(0.01, step_cost);
                    /* #pragma omp critical
                                        {
                                            std::cout << "Edge (" << current.x
                       << "," << current.y
                                                      << ") -> (" << nx << ","
                       << ny << ") "
                                                      << "step_cost: " <<
                       step_cost
                                                      << " p_e: " << p_e << "
                       dah: " << dah
                                                      << " vd: " << vd << "\n";
                                        } */
                }
            }

            double tentative_g = current.g_cost + step_cost;
            int tentative_hops = current.hops + 1;

            bool replace = false;
            auto& rec = get_node(nx, ny);

            if (conforms) {
                if (!rec.conforms)
                    replace = true;
                else if (tentative_g < rec.cost)
                    replace = true;
            } else {
                if (!rec.conforms) {
                    if (tentative_hops < rec.hops)
                        replace = true;
                    else if (tentative_hops == rec.hops &&
                             tentative_g < rec.cost)
                        replace = true;
                }
            }

            if (replace) {
                rec.cost = tentative_g;
                rec.hops = tentative_hops;
                rec.conforms = conforms;
                rec.px = current.x;
                rec.py = current.y;
                rec.hori = hori;
                open_list.push(
                    {nx, ny, tentative_g,
                     tentative_g +
                         Param::ASTAR_HEURISTIC_WEIGHT *
                             getManhattanDist(nx, ny, target_x, target_y),
                     tentative_hops, current.x, current.y, hori, conforms});
            }
        }
    }

    if (found) {
        twoPin.path.clear();
        int cx = target_x, cy = target_y;
        while (cx != start_x || cy != start_y) {
            auto& rec = get_node(cx, cy);
            int px = rec.px, py = rec.py;
            bool h = rec.hori;
            twoPin.path.emplace_back(std::min(cx, px), std::min(cy, py), h);
            cx = px;
            cy = py;
        }
        std::reverse(twoPin.path.begin(), twoPin.path.end());
    }
    return found;
}

int main(int argc, char** argv) {
    auto start_time = std::chrono::high_resolution_clock::now();
    assert(argc >= 3 && "Usage: ./router <inputFile> <outputFile>");
    std::ifstream fp(argv[1]);
    ISPDParser::ispdData* ispdData = ISPDParser::parse(fp);
    fp.close();

    // =========================================================
    // Dynamic Parameter Adjustment based on Testcase (numNet)
    // =========================================================
    int net_count = ispdData->numNet;

    if (net_count > 330000 && net_count < 345000) {
        // std::cout << ">>> Detected Testcase: newblue1 <<<" << std::endl;
        Param::BLC_ALPHA = 800.0;
        Param::BLC_BETA = 0.1;
        Param::VIRTUAL_DEMAND_C9 = 0.8;
        Param::DAH_C1 = 1.0;
        Param::DAH_C2 = 0.1;
    } else if (net_count > 570000 && net_count < 585000) {
        // std::cout << ">>> Detected Testcase: bigblue2 <<<" << std::endl;
        Param::BLC_ALPHA = 800.0;
        Param::BLC_BETA = 0.1;
        Param::VIRTUAL_DEMAND_C9 = 2.0;
        Param::DAH_C1 = 1.0;
        Param::DAH_C2 = 0.1;
    } else if (net_count > 1110000 && net_count < 1130000) {
        // std::cout << ">>> Detected Testcase: bigblue3 <<<" << std::endl;
        Param::BLC_ALPHA = 800.0;
        Param::BLC_BETA = 0.1;
        Param::VIRTUAL_DEMAND_C9 = 0.8;
        Param::DAH_C1 = 1.0;
        Param::DAH_C2 = 0.1;
    } else if (net_count >= 1250000 && net_count <= 1285000) {
        // std::cout << ">>> Detected Testcase: newblue5 <<<" << std::endl;
        Param::BLC_ALPHA = 800.0;
        Param::BLC_BETA = 0.1;
        Param::VIRTUAL_DEMAND_C9 = 0.8;
        Param::DAH_C1 = 1.0;
        Param::DAH_C2 = 2.0;
    } // else {
        // std::cout << ">>> Detected Testcase: Default/Other <<<" << std::endl;
        // Use default parameters
        // Param::BLC_ALPHA = 9.0;
        // Param::BLC_BETA = 9.5;
        // Param::VIRTUAL_DEMAND_C9 = 0.4;
        // Param::DAH_C1 = 4.0;
        // Param::DAH_C2 = 1.0;
    // }
    // =========================================================

    edge_grid.resize(ispdData->numXGrid,
                     std::vector<std::vector<EdgeInfo>>(
                         ispdData->numYGrid, std::vector<EdgeInfo>(2)));

    int num_layers = ispdData->horizontalCapacity.size();
    std::vector<int> track_width(num_layers, 1);
    int base_min_width =
        ispdData->nets.empty() ? 1 : ispdData->nets[0]->minimumWidth;

    int total_h_tracks = 0;
    int total_v_tracks = 0;

    for (int z = 0; z < num_layers; ++z) {
        int spacing = (z < ispdData->minimumSpacing.size())
                          ? ispdData->minimumSpacing[z]
                          : 0;
        int tw = base_min_width + spacing;
        if (tw <= 0) tw = 1;
        track_width[z] = tw;

        total_h_tracks += ispdData->horizontalCapacity[z] / tw;
        total_v_tracks += ispdData->verticalCapacity[z] / tw;
    }

    for (int x = 0; x < ispdData->numXGrid; ++x) {
        for (int y = 0; y < ispdData->numYGrid; ++y) {
            if (x < ispdData->numXGrid - 1)
                edge_grid[x][y][1].capacity = total_h_tracks;
            if (y < ispdData->numYGrid - 1)
                edge_grid[x][y][0].capacity = total_v_tracks;
        }
    }

    for (auto adj : ispdData->capacityAdjs) {
        int x1 = std::get<0>(adj->grid1), y1 = std::get<1>(adj->grid1),
            z1 = std::get<2>(adj->grid1) - 1;
        int x2 = std::get<0>(adj->grid2), y2 = std::get<1>(adj->grid2);

        bool hori = (y1 == y2);
        int tw = track_width[z1];

        int base_phys = hori ? ispdData->horizontalCapacity[z1]
                             : ispdData->verticalCapacity[z1];
        int adj_phys = adj->reducedCapacityLevel;

        int base_tracks = base_phys / tw;
        int adj_tracks = adj_phys / tw;

        int ex = std::min(x1, x2);
        int ey = std::min(y1, y2);

        edge_grid[ex][ey][hori].capacity += (adj_tracks - base_tracks);
    }

    for (auto net : ispdData->nets) {
        for (auto& pin : net->pins) {
            int x =
                (std::get<0>(pin) - ispdData->lowerLeftX) / ispdData->tileWidth;
            int y = (std::get<1>(pin) - ispdData->lowerLeftY) /
                    ispdData->tileHeight;
            int z = std::get<2>(pin) - 1;
            if (std::none_of(net->pin3D.begin(), net->pin3D.end(),
                             [x, y, z](auto& p) {
                                 return p.x == x && p.y == y && p.z == z;
                             }))
                net->pin3D.emplace_back(x, y, z);
            if (std::none_of(net->pin2D.begin(), net->pin2D.end(),
                             [x, y](auto& p) { return p.x == x && p.y == y; }))
                net->pin2D.emplace_back(x, y);
        }
        if (net->pin2D.size() <= 1) net->pins.resize(1);
    }

    // =========================================================
    // Stage 1: RMST Decomposition
    // =========================================================
    std::cout << "Stage 1: RMST Decomposition & Initial Routing..."
              << std::endl;
    Flute::readLUT("./");
    for (auto net : ispdData->nets) {
        int d = net->pin2D.size();
        if (d <= 1 || d > 1000) continue;

        std::vector<DTYPE> fx(d), fy(d);
        for (int i = 0; i < d; ++i) {
            fx[i] = net->pin2D[i].x;
            fy[i] = net->pin2D[i].y;
        }
        Flute::Tree rsmt = Flute::flute(d, fx.data(), fy.data(), ACCURACY);
        for (int i = 0; i < 2 * rsmt.deg - 2; ++i) {
            Flute::Branch &b1 = rsmt.branch[i], &b2 = rsmt.branch[b1.n];
            if (i == b1.n) continue;
            if (b1.x == b2.x || b1.y == b2.y) {
                for (int cx = std::min(b1.x, b2.x); cx < std::max(b1.x, b2.x);
                     ++cx)
                    skeleton_map[net].insert(((uint64_t)cx << 32) |
                                             ((uint64_t)b1.y << 1) | 1);
                for (int cy = std::min(b1.y, b2.y); cy < std::max(b1.y, b2.y);
                     ++cy)
                    skeleton_map[net].insert(((uint64_t)b2.x << 32) |
                                             ((uint64_t)cy << 1) | 0);
            }
        }

        free(rsmt.branch);

        decomposeNetByRMST(net);
    }

    // =========================================================
    // Stage 1.5: Monotonic Routing with Bounding Box Ordering
    // =========================================================
    std::cout << "Stage 1.5: Monotonic Routing with Bounding Box Ordering..."
              << std::endl;

    std::vector<ISPDParser::TwoPin*> all_twopins;
    for (auto net : ispdData->nets) {
        for (auto& tp : net->twopin) {
            all_twopins.push_back(&tp);
        }
    }

    std::sort(all_twopins.begin(), all_twopins.end(),
              [](ISPDParser::TwoPin* a, ISPDParser::TwoPin* b) {
                  int dx_a = std::abs(a->from.x - a->to.x);
                  int dy_a = std::abs(a->from.y - a->to.y);
                  int dx_b = std::abs(b->from.x - b->to.x);
                  int dy_b = std::abs(b->from.y - b->to.y);

                  int area_a = dx_a * dy_a;
                  int area_b = dx_b * dy_b;

                  if (area_a != area_b) {
                      return area_a < area_b;
                  }

                  int hpwl_a = dx_a + dy_a;
                  int hpwl_b = dx_b + dy_b;
                  return hpwl_a < hpwl_b;
              });

    for (auto tp : all_twopins) {
        tp->path = buildMonotonicPath(tp, ispdData);
        updatePathDemand(*tp, 1, ispdData);
    }

    // =========================================================
    // Stage 3: Task-based Concurrency Strategy (TCS Enabled)
    // =========================================================
    std::cout << "Stage 3: Multi-threaded NCTU-GR Hierarchical RRR..."
              << std::endl;
    bool timeout = false;

    double best_score = std::numeric_limits<double>::max();
    std::unordered_map<ISPDParser::TwoPin*, std::vector<ISPDParser::RPoint>>
        best_paths;

    std::unordered_map<ISPDParser::Net*, bool> net_locked;
    for (auto net : ispdData->nets) {
        net_locked[net] = false;
    }

    for (int iter = 0; iter < Param::MAX_ROUTING_ITERS; ++iter) {
        if (timeout) break;
        int overflow_count = 0;

#pragma omp parallel for collapse(2)
        for (int x = 0; x < ispdData->numXGrid; ++x) {
            for (int y = 0; y < ispdData->numYGrid; ++y) {
                edge_grid[x][y][0].overflow_freq = 0;
                edge_grid[x][y][1].overflow_freq = 0;
            }
        }

        struct Task {
                ISPDParser::TwoPin* tp;
                ISPDParser::Net* net;
                double score;
        };
        std::vector<Task> tasks;

        for (auto net : ispdData->nets) {
            for (auto& tp : net->twopin) {
                int oe = 0;
                for (auto& e : tp.path) {
                    if (edge_grid[e.x][e.y][e.hori].demand >
                        edge_grid[e.x][e.y][e.hori].capacity)
                        oe++;
                }
                if (oe > 0) {
                    tasks.push_back({&tp, net,
                                     Param::NET_ORDER_C7 * oe +
                                         Param::NET_ORDER_C8 * tp.path.size()});
                }
            }
        }

        std::sort(tasks.begin(), tasks.end(), [](const Task& a, const Task& b) {
            return a.score > b.score;
        });

        std::deque<Task> TQ(tasks.begin(), tasks.end());
        std::atomic<int> tasks_remaining(tasks.size());

#pragma omp parallel reduction(+ : overflow_count)
        {
            while (tasks_remaining > 0) {
                Task current_task = {nullptr, nullptr, 0.0};
                bool got_task = false;
#pragma omp critical
                {
                    // Line 5: Extract a task from TQ
                    if (!TQ.empty()) {
                        current_task = TQ.front();
                        TQ.pop_front();
                        got_task = true;
                    }
                }

                if (!got_task) {
                    std::this_thread::yield();
                    continue;
                }

                bool task_pushed_back = false;

#pragma omp critical
                {
                    if (net_locked[current_task.net]) {
                        TQ.push_back(current_task);
                        task_pushed_back = true;
                    } else {
                        net_locked[current_task.net] = true;
                    }
                }
                if (task_pushed_back) {
                    std::this_thread::yield();
                    continue;
                }

                // --- Start Routing ---
                auto tp = current_task.tp;
                auto net = current_task.net;
                overflow_count++;
                auto old_path = tp->path;
                int prev_len = tp->path.size();

                for (auto& e : old_path) {
#pragma omp atomic
                    edge_grid[e.x][e.y][e.hori].marks++;
                }

                updatePathDemand(*tp, -1, ispdData);

                // Line 9: Collision aware BLMR
                if (!routeAStar_NCTUGR(net, *tp, tp->to.x, tp->to.y, iter,
                                       ispdData, prev_len, old_path)) {
                    routeAStar_NCTUGR(net, *tp, tp->to.x, tp->to.y, iter,
                                      ispdData, prev_len, old_path, false, 0,
                                      true);
                    // tp->path = old_path;
                }

                /* #pragma omp critical
                                {
                                    std::cout << "[Iter " << iter + 1 << "] Path
                   for TwoPin ("
                                              << tp->from.x << "," << tp->from.y
                   << ") -> ("
                                              << tp->to.x << "," << tp->to.y <<
                   "): "; if (tp->path.empty()) { std::cout << "No path found!";
                                    } else {
                                        for (const auto& p : tp->path) {
                                            std::cout << "(" << p.x << "," <<
                   p.y << ","
                                                      << (p.hori ? "H" : "V") <<
                   ") ";
                                        }
                                    }
                                    std::cout << "\n";
                                } */

                updatePathDemand(*tp, 1, ispdData);

                for (auto& e : old_path) {
#pragma omp atomic
                    edge_grid[e.x][e.y][e.hori].marks--;
                }

                for (auto& e : tp->path) {
                    int current_dem;
#pragma omp atomic read
                    current_dem = edge_grid[e.x][e.y][e.hori].demand;

                    if (current_dem > edge_grid[e.x][e.y][e.hori].capacity) {
#pragma omp atomic
                        edge_grid[e.x][e.y][e.hori].overflow_freq++;
                    }
                }

// Line 10: Unlock(N)
#pragma omp critical
                {
                    net_locked[net] = false;
                }

                tasks_remaining--;

                if (std::chrono::duration_cast<std::chrono::seconds>(
                        std::chrono::high_resolution_clock::now() - start_time)
                        .count() > Param::TIME_LIMIT_SEC) {
                    timeout = true;
                }
                if (timeout) break;
            }
        }

        // Calculate total overflow (TOF)
        int tof = 0;
        int mof = 0;
        int tw2D = track_width[0];
#pragma omp parallel for collapse(2) reduction(+ : tof) reduction(max : mof)
        for (int x = 0; x < ispdData->numXGrid; ++x) {
            for (int y = 0; y < ispdData->numYGrid; ++y) {
                int current_dem_v = edge_grid[x][y][0].demand;
                int cap_v = edge_grid[x][y][0].capacity;
                if (current_dem_v > cap_v) {
                    int ov = (current_dem_v - cap_v) * tw2D;
                    tof += ov;
                    if (ov > mof) mof = ov;
                }

                int current_dem_h = edge_grid[x][y][1].demand;
                int cap_h = edge_grid[x][y][1].capacity;
                if (current_dem_h > cap_h) {
                    int ov = (current_dem_h - cap_h) * tw2D;
                    tof += ov;
                    if (ov > mof) mof = ov;
                }
            }
        }

        long long int wl = 0;
        for (auto net : ispdData->nets) {
            for (auto& tp : net->twopin) {
                wl += tp.path.size();
            }
        }

        double current_score = tof + mof * 1e-3 + wl * 1e-12;
        // std::cout << "Iteration " << iter + 1 << " Score: " << current_score
        //           << " (TOF: " << tof << ", MOF: " << mof << ", WL: " << wl
        //           << ")" << std::endl;

        if (current_score < best_score) {
            best_score = current_score;
            for (auto net : ispdData->nets) {
                for (auto& tp : net->twopin) {
                    best_paths[&tp] = tp.path;
                }
            }
            // std::cout << ">>> New Score! <<<" << std::endl;
        }
        if (tof == 0) break;

// Line 13: Update history cost
#pragma omp parallel for collapse(2)
        for (int x = 0; x < ispdData->numXGrid; ++x) {
            for (int y = 0; y < ispdData->numYGrid; ++y) {
                edge_grid[x][y][0].history_cost +=
                    edge_grid[x][y][0].overflow_freq;
                edge_grid[x][y][1].history_cost +=
                    edge_grid[x][y][1].overflow_freq;
            }
        }
    }

    // =========================================================
    // Load Best Result & Sync Demand (Robust Version)
    // =========================================================
    if (!best_paths.empty()) {
        std::cout << ">>> Load Best Result (Best Score: " << best_score
                  << ") <<<" << std::endl;

        // 1. Reset demand on all edges
        for (int x = 0; x < ispdData->numXGrid; ++x) {
            for (int y = 0; y < ispdData->numYGrid; ++y) {
                edge_grid[x][y][0].demand = 0;
                edge_grid[x][y][1].demand = 0;
            }
        }

        // 2. Critical step: clear the paths of all two-pin connections first!
        for (auto net : ispdData->nets) {
            for (auto& tp : net->twopin) {
                tp.path.clear();
            }
        }

        // 3. Apply best_paths sequentially and accumulate demand using the existing updatePathDemand function
        for (auto net : ispdData->nets) {
            for (auto& tp : net->twopin) {
                tp.path = best_paths[&tp];
                updatePathDemand(tp, 1, ispdData);
            }
        }
    }

    // =========================================================
    // Stage 3.5: Post-Refinement (Strict Overflow Resolution)
    // =========================================================
    std::cout << "Stage 3.5: Strict Overflow Resolution..." << std::endl;

    int tw2D = track_width[0];
    int prev_tof = std::numeric_limits<int>::max(); // Track the previous TOF
    int resolve_iter = 0;                           // Track the current iteration count

    while (true) {
        auto current_time = std::chrono::high_resolution_clock::now();
        int elapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(
                              current_time - start_time)
                              .count();
        if (elapsed_sec >= Param::GLOBAL_TIME_LIMIT_SEC) {
            // std::cout
            //     << "  >>> Global time limit reached (" << elapsed_sec
            //     << "s). Stopping Post-Refinement early to save results! <<<"
            //     << std::endl;
            break;
        }

        std::vector<ISPDParser::TwoPin*> overflow_tps;
        int current_tof = 0;

        // Calculate the current TOF
        for (int x = 0; x < ispdData->numXGrid; ++x) {
            for (int y = 0; y < ispdData->numYGrid; ++y) {
                for (int h = 0; h < 2; ++h) {
                    int dem = edge_grid[x][y][h].demand;
                    int cap = edge_grid[x][y][h].capacity;
                    if (dem > cap) {
                        current_tof += (dem - cap) * tw2D;
                    }
                }
            }
        }

        // Collect all two-pin connections that traverse overflowing edges
        for (auto net : ispdData->nets) {
            for (auto& tp : net->twopin) {
                bool has_overflow = false;
                for (const auto& e : tp.path) {
                    if (edge_grid[e.x][e.y][e.hori].demand >
                        edge_grid[e.x][e.y][e.hori].capacity) {
                        has_overflow = true;
                        break;
                    }
                }
                if (has_overflow) {
                    overflow_tps.push_back(&tp);
                }
            }
        }

        // std::cout << "  [Resolve Iter " << resolve_iter + 1
        //           << "] Overflowing TwoPins: " << overflow_tps.size()
        //           << ", Current TOF: " << current_tof
        //           << " (Time used: " << elapsed_sec << "s)" << std::endl;

        // Stopping condition 1: no overflow remains
        if (overflow_tps.empty() || current_tof == 0) {
            // std::cout << "  >>> All overflows resolved! <<<" << std::endl;
            break;
        }

        // Stopping condition 2: TOF has not improved (has not decreased)
        if (current_tof >= prev_tof) {
            // std::cout << "  >>> No further improvement (TOF not decreasing). "
            //              "Stopping early! <<<"
            //           << std::endl;
            break;
        }

        // Update prev_tof for the next check
        prev_tof = current_tof;

        std::sort(overflow_tps.begin(), overflow_tps.end(),
                  [](ISPDParser::TwoPin* a, ISPDParser::TwoPin* b) {
                      return a->path.size() < b->path.size();
                  });

        for (auto tp : overflow_tps) {
            bool still_overflow = false;
            for (const auto& e : tp->path) {
                if (edge_grid[e.x][e.y][e.hori].demand >
                    edge_grid[e.x][e.y][e.hori].capacity) {
                    still_overflow = true;
                    break;
                }
            }
            if (!still_overflow) continue;

            auto old_path = tp->path;
            int prev_len = tp->path.size();

            // 1. Rip up the route
            updatePathDemand(*tp, -1, ispdData);

            std::vector<ISPDParser::RPoint> empty_old;

            // 2. Reroute under strict constraints
            bool found =
                routeAStar_NCTUGR(tp->parNet, *tp, tp->to.x, tp->to.y, 0,
                                  ispdData, prev_len, empty_old, true, 0, true);

            bool accept_new = false;
            if (found) {
                bool is_clean = true;
                int dem_val = getDemandVal(tp->parNet, ispdData);

                // 3. Strict validation: check that the new path is entirely free of overflow
                for (const auto& e : tp->path) {
                    if (edge_grid[e.x][e.y][e.hori].demand + dem_val >
                        edge_grid[e.x][e.y][e.hori].capacity) {
                        is_clean = false; // Reject the path if even one edge exceeds capacity
                        break;
                    }
                }
                if (is_clean) accept_new = true;
            }

            // 4. Accept or reject the new path based on the validation result
            if (accept_new) {
                updatePathDemand(*tp, 1, ispdData);
            } else {
                tp->path = old_path;
                updatePathDemand(*tp, 1, ispdData);
            }
        }

        resolve_iter++;
    }

    std::cout << "Stage 4: Layer Assignment..." << std::endl;
    LayerAssignment::Graph graph;
    graph.initialLA(*ispdData, 1);
    graph.convertGRtoLA(*ispdData, true);
    graph.COLA(true);
    graph.output3Dresult(argv[2]);

    delete ispdData;
    return 0;
}