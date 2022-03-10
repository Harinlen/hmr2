#include <algorithm>
#include <queue>
#include <thread>
#include <cmath>

#include "hmr_global.h"
#include "hmr_cgd_parse.h"
#include "hmr_ui.h"

#include "partition.h"

typedef std::unordered_map<int, double> EDGE_MAP;

typedef struct NODE_INFO
{
    size_t trust_pos;
    EDGE_MAP edge_map;
    CONTIG_NODE *node;
} NODE_INFO;

inline double edge_weight(const EDGE_MAP &map, int node_id)
{
    const auto finder = map.find(node_id);
    return finder == map.end() ? 0.0 : finder->second;
}

typedef struct TRUST_EDGE
{
    int start;
    int end;
    double weight;
} TRUST_EDGE;

bool compare_trust_edge(const TRUST_EDGE &lhs, const TRUST_EDGE &rhs)
{
    return lhs.weight > rhs.weight;
}

template <typename T>
std::vector<T> derivative(const std::vector<T> &fx)
{
    std::vector<T> result;
    result.reserve(fx.size() - 1);
    //Loop and calculate the position.
    for(size_t i=1; i<fx.size(); ++i)
    {
        result.push_back(hAbs(fx[i] - fx[i-1]));
    }
    return result;
}

template <typename T>
size_t maximum_pos(const std::vector<T> &fx)
{
    return std::distance(fx.begin(), std::max_element(fx.begin(), fx.end()));
}

typedef union EDGE_INFO
{
    typedef struct EDGE
    {
        int32_t start;
        int32_t end;
    } EDGE;
    EDGE edge;
    uint64_t data;
} EDGE_INFO;

uint64_t pack_edge_data(int32_t x, int32_t y)
{
    EDGE_INFO info;
    if(x < y)
    {
        info.edge = EDGE_INFO::EDGE {x, y};
    }
    else
    {
        info.edge = EDGE_INFO::EDGE {y, x};
    }
    return info.data;
}

void unpack_edge_data(uint64_t data, int32_t &x, int32_t &y)
{
    EDGE_INFO info;
    info.data = data;
    x = info.edge.start;
    y = info.edge.end;
}

typedef struct EDGE_VOTE_RESULT
{
    uint64_t edge;
    size_t votes;
} EDGE_VOTE_RESULT;

typedef std::unordered_set<int32_t> NODE_ID_SET;
typedef std::unordered_map<uint64_t, NODE_ID_SET> EDGE_VOTER;

void vote_edge(EDGE_VOTER &voter, int32_t x, int32_t y, int32_t host_node)
{
    uint64_t edge_key = pack_edge_data(x, y);
    auto edge_finder = voter.find(edge_key);
    if(edge_finder == voter.end())
    {
        std::unordered_set<int32_t> node_set;
        node_set.insert(host_node);
        voter.insert(std::make_pair(edge_key, node_set));
    }
    else
    {
        edge_finder->second.insert(host_node);
    }
}

typedef struct NODE_SET_MARK
{
    NODE_ID_SET node_ids;
    int subsets;
} NODE_SET_MARK;

bool is_subset(const NODE_ID_SET &x, const NODE_ID_SET &y)
{
    //Likely.
    if(y.size() <= x.size())
    {
        for(int32_t y_id: y)
        {
            if(!hInSet(y_id, x))
            {
                return false;
            }
        }
        return true;
    }
    return false;
}

std::vector<size_t> find_belongs(const std::vector<NODE_SET_MARK> &node_sets, const NODE_ID_SET &id_set)
{
    std::list<size_t> belong_ids;
    for(size_t i=0; i<node_sets.size(); ++i)
    {
        if(is_subset(node_sets[i].node_ids, id_set))
        {
            belong_ids.push_back(i);
        }
    }
    return std::vector<size_t>(belong_ids.begin(), belong_ids.end());
}

typedef struct NODE_SET_RELATION
{
    int set_a;
    int set_b;
    NODE_ID_SET inter;
} NODE_SET_RELATION;

template <typename T>
bool has_intersection(const std::unordered_set<T> &x, const std::unordered_set<T> &y)
{
    for(int y_id: y)
    {
        if(hInSet(y_id, x))
        {
            return true;
        }
    }
    return false;
}

NODE_ID_SET set_intersection(const NODE_ID_SET &x, const NODE_ID_SET &y)
{
    std::list<int> intersection_ids;
    for(int y_id: y)
    {
        if(hInSet(y_id, x))
        {
            intersection_ids.push_back(y_id);
        }
    }
    return NODE_ID_SET(intersection_ids.begin(), intersection_ids.end());
}

void extract_set(NODE_ID_SET &x, NODE_ID_SET &y, NODE_ID_SET &intersection)
{
    std::list<int> intersection_ids;
    for(int y_id: y)
    {
        if(hInSet(y_id, x))
        {
            intersection_ids.push_back(y_id);
        }
    }
    if(!intersection_ids.empty())
    {
        //Time the modify x and y.
        for(int node_id: intersection_ids)
        {
            x.erase(node_id);
            y.erase(node_id);
        }
        intersection = NODE_ID_SET(intersection_ids.begin(), intersection_ids.end());
    }
}

typedef struct NODE_SET_RELATION_GROUP
{
    NODE_ID_SET inter;
    std::unordered_set<int> set_ids;
    int relation_count;
} NODE_SET_RELATION_GROUP;

bool has_cross_ids(const std::vector<NODE_ID_SET> &node_id_sets)
{
    //Loop check each pair.
    for(size_t i=0; i<node_id_sets.size(); ++i)
    {
        for(size_t j=i+1; j<node_id_sets.size(); ++j)
        {
            if(has_intersection(node_id_sets[i], node_id_sets[j]))
            {
                return true;
            }
        }
    }
    return false;
}

double group_relation(const NODE_ID_SET &lhs, const NODE_ID_SET &rhs, NODE_INFO *node_info, size_t edge_limit)
{
    //Loop for all edges in smaller group.
    std::priority_queue<double, std::vector<double>, std::less<double> > weight_queue;
    //Find the first matched ids from the big group.
    for(int s_id: rhs)
    {
        //The edges are already sorted descending, the first appeared edge is
        //the best match edge.
        const auto &edges = node_info[s_id].edge_map;
        for(const auto &edge: edges)
        {
            if(hInSet(edge.first, lhs))
            {
                if(weight_queue.size() < edge_limit)
                {
                    //When less than heap size, just insert.
                    weight_queue.push(edge.second);
                }
                else
                {
                    //When reach the limitation, compare with the minimum value
                    //of the heap.
                    if(edge.second > weight_queue.top())
                    {
                        weight_queue.pop();
                        weight_queue.push(edge.second);
                    }
                }
            }
        }
    }
    //Sum up the max weighted edges.
    double relation = 0.0;
    while(!weight_queue.empty())
    {
        relation += weight_queue.top();
        weight_queue.pop();
    }
    return relation;
}

typedef std::unordered_set<int> PIECES_BELONGS;

typedef struct PIECES_GROUP
{
    PIECES_BELONGS belongs;
    std::unordered_set<int> pieces_ids;
} PIECES_GROUP;

typedef struct HANA_GROUP
{
    NODE_ID_SET ids;
    size_t length;
} HANA_GROUP;

typedef struct UNUSED_NODE
{
    int node_id;
    double best_edge;
} UNUSED_NODE;

typedef struct GROUP_PREDICT
{
    int group_id;
    double mark;
} GROUP_PREDICT;

typedef std::vector<GROUP_PREDICT> GROUP_PREDICTION;

GROUP_PREDICTION node_predict(int node_id, CONTIG_NODE *nodes, NODE_INFO *node_infos,
                              const std::vector<HANA_GROUP> &core_groups)
{
    //Loop and find the marking limit.
    size_t edge_limit = core_groups[0].ids.size();
    for(size_t i=1; i<core_groups.size(); ++i)
    {
        edge_limit = hMin(edge_limit, core_groups[i].ids.size());
    }
    //Loop and generate the weight array for all core groups.
    GROUP_PREDICTION results;
    results.reserve(core_groups.size());
    size_t *edge_counter = new size_t[core_groups.size()];
    for(size_t i=0; i<core_groups.size(); ++i)
    {
        edge_counter[i] = 0;
        results.push_back(GROUP_PREDICT {static_cast<int>(i), -1.0});
    }
    //Loop and calculate the limit.
    const auto &edges = nodes[node_id].links;
    size_t trust_pos = node_infos[node_id].trust_pos;
    for(size_t edge_id=0; edge_id < edges.size(); ++edge_id)
    {
        const auto &edge = edges[edge_id];
        for(size_t i=0; i<core_groups.size(); ++i)
        {
            if(edge_counter[i] < edge_limit && hInSet(edge.id, core_groups[i].ids))
            {
                double edge_w = edge.weight;
                //Cut the weight if the id is greater than trust pos.
                if(edge_id > trust_pos)
                {
                    edge_w /= 2.0;
                }
                //Initial the result.
                if(results[i].mark < 0)
                {
                    results[i].mark = 0.0;
                }
                results[i].mark += edge.weight;
                ++edge_counter[i];
                break;
            }
        }
        size_t full_count = 0;
        for(size_t i=0; i<core_groups.size(); ++i)
        {
            if(edge_counter[i] == edge_limit)
            {
                ++full_count;
            }
        }
        if(full_count == core_groups.size())
        {
            break;
        }
    }
    delete[] edge_counter;
    //Sort the result.
    std::sort(results.begin(), results.end(), [](const GROUP_PREDICT &lhs, const GROUP_PREDICT &rhs) {
        return lhs.mark > rhs.mark;
    });
    return results;
}

void group_hanamaru(CONTIG_NODE *nodes, size_t node_size, int no_of_group, int threads,
                    std::vector<std::vector<int> > &groups, std::list<int> &unknown_ids)
{
    time_print("Generating node informations...");
    NODE_INFO *node_infos = new NODE_INFO[node_size];
    std::unordered_set<int> best_nodes;
    {
        size_t trust_edge_size = (node_size * node_size) >> 2;
        std::priority_queue<TRUST_EDGE, std::vector<TRUST_EDGE>, decltype(&compare_trust_edge)> best_edge_queue(compare_trust_edge);
        //Initial the node information.
        for(size_t i=0; i<node_size; ++i)
        {
            //Initial the parent.
            node_infos[i].node = &(nodes[i]);
            //Create the edge map.
            auto &edges = nodes[i].links;
            std::sort(edges.begin(), edges.end(), [](const CONTIG_EDGE &lhs, const CONTIG_EDGE &rhs) {
                return lhs.weight > rhs.weight;
            });
            //Construct the edge map.
            auto &edge_map = node_infos[i].edge_map;
            for(const auto &edge: edges)
            {
                edge_map.insert(std::make_pair(edge.id, edge.weight));
                //Append the edge information to best edge.
                TRUST_EDGE edge_info {static_cast<int>(i), edge.id, edge.weight};
                if(best_edge_queue.size() < trust_edge_size)
                {
                    //Directly add the edge to the queue.
                    best_edge_queue.push(edge_info);
                }
                else
                {
                    if(best_edge_queue.top().weight < edge.weight)
                    {
                        //Insert it.
                        best_edge_queue.pop();
                        best_edge_queue.push(edge_info);
                    }
                }
            }
            //Calculate the trust edge size.
            node_infos[i].trust_pos = 0;
            if(edges.size() > 2)
            {
                //Calculate the deviation.
                std::vector<double> edge_weights;
                edge_weights.reserve(edges.size());
                for(const auto &edge: edges)
                {
                    //Calculate the difference of the left and right side.
                    edge_weights.push_back(edge.weight);
                }
                //Calculate the 1st derivative result.
                std::vector<double> median_diff_deri = derivative(edge_weights);
                //Find the maximum position of the derivcative.
                node_infos[i].trust_pos = maximum_pos(median_diff_deri) + 1;
            }
            else
            {
                node_infos[i].trust_pos = 1;
            }
        }
        //Gather the trust edges.
        size_t edge_size = best_edge_queue.size(), edge_size_1 = edge_size - 1;
        std::vector<TRUST_EDGE> best_edges;
        best_edges.resize(edge_size);
        for(size_t i=0; i<edge_size; ++i)
        {
            best_edges[edge_size_1 - i] = best_edge_queue.top();
            best_edge_queue.pop();
        }
        //Pick up the best nodes.
        size_t expected_nodes = node_size >> 1, edge_id = 0;
        while(best_nodes.size() < expected_nodes && edge_id < edge_size)
        {
            const auto &edge = best_edges[edge_id];
            best_nodes.insert(edge.start);
            best_nodes.insert(edge.end);
            ++edge_id;
        }
    }

    // -- HANA Stage --
    time_print("HANA Stage...");
    std::vector<NODE_SET_MARK> node_set_marks;
    {
        std::vector<EDGE_VOTE_RESULT> edge_vote_results;
        EDGE_VOTER edge_voter;
        //Vote the best edges from the best nodes.
        size_t vote_edge_range = no_of_group * 2;
        time_print("GCN for edges, window size: %zu", vote_edge_range);
        for(int node_id: best_nodes)
        {
            const auto &edge_list = nodes[node_id].links;
            size_t edge_list_range = hMin(vote_edge_range, edge_list.size());
            for(size_t i=0; i<edge_list_range; ++i)
            {
                if(hInSet(edge_list[i].id, best_nodes))
                {
                    vote_edge(edge_voter, node_id, edge_list[i].id, node_id);
                }
            }
            for(size_t i=0; i<edge_list_range-1; ++i)
            {
                const int32_t i_id = edge_list[i].id;
                if(!hInSet(i_id, best_nodes))
                {
                    continue;
                }
                for(size_t j=i+1; j<edge_list_range; ++j)
                {
                    if(hInSet(edge_list[j].id, best_nodes))
                    {
                        vote_edge(edge_voter, i_id, edge_list[j].id, node_id);
                    }
                }
            }
        }
        //Sort the edge voting result.
        edge_vote_results.reserve(edge_voter.size());
        for(auto edge_vote: edge_voter)
        {
            edge_vote_results.push_back(EDGE_VOTE_RESULT {edge_vote.first, edge_vote.second.size()});
        }
        std::sort(edge_vote_results.begin(), edge_vote_results.end(), [](const EDGE_VOTE_RESULT &lhs, const EDGE_VOTE_RESULT &rhs) {
            return lhs.votes > rhs.votes;
        });
        //Extract the expected node id group.
        time_print("Extracting node set from the trusted edges...");
        node_set_marks.reserve(edge_vote_results.size());
        for(const auto &edge_result: edge_vote_results)
        {
            const auto node_set = edge_voter.find(edge_result.edge)->second;
            std::vector<size_t> subset_ids = find_belongs(node_set_marks, node_set);
            if(subset_ids.empty())
            {
                //Add this node set to the node set marks.
                node_set_marks.push_back(NODE_SET_MARK {node_set, 0});
            }
            else
            {
                for(size_t set_id: subset_ids)
                {
                    ++node_set_marks[set_id].subsets;
                }
            }
        }
        //Find out the subset marks.
        node_set_marks.reserve(node_set_marks.size());
        std::sort(node_set_marks.begin(), node_set_marks.end(), [](const NODE_SET_MARK &lhs, const NODE_SET_MARK &rhs) {
            return lhs.subsets > rhs.subsets;
        });
    }
    //Calculate the intersection of the node sets.
    time_print("Finding the intersection of the node sets...");
    std::vector<NODE_SET_RELATION> node_set_relations;
    node_set_relations.reserve(node_set_marks.size() * node_set_marks.size());
    for(size_t i=0; i<node_set_marks.size()-1; ++i)
    {
        const NODE_ID_SET &i_set = node_set_marks[i].node_ids;
        for(size_t j=i+1; j<node_set_marks.size(); ++j)
        {
            NODE_ID_SET inter = set_intersection(i_set, node_set_marks[j].node_ids);
            if(inter.size() > 1)
            {
                node_set_relations.push_back(
                            NODE_SET_RELATION {static_cast<int>(i),
                                               static_cast<int>(j), inter});
            }
        }
    }
    node_set_relations.reserve(node_set_relations.size());
    std::sort(node_set_relations.begin(), node_set_relations.end(),
              [](const NODE_SET_RELATION &lhs, const NODE_SET_RELATION &rhs) {
        return lhs.inter.size() > rhs.inter.size();
    });
    //Find out the largest subset of the relations.
    time_print("Merging the intersection of the node sets...");
    std::vector<NODE_SET_RELATION_GROUP> relation_groups;
    relation_groups.reserve(node_set_relations.size());
    for(const auto &relation: node_set_relations)
    {
        //Add relation to relation groups.
        bool relation_added = false;
        for(auto &relation_group: relation_groups)
        {
            if(is_subset(relation_group.inter, relation.inter))
            {
                //Add the set id to relation group.
                relation_group.set_ids.insert(relation.set_a);
                relation_group.set_ids.insert(relation.set_b);
                ++relation_group.relation_count;
                relation_added = true;
                break;
            }
        }
        if(!relation_added)
        {
            std::unordered_set<int> set_ids;
            set_ids.insert(relation.set_a);
            set_ids.insert(relation.set_b);
            relation_groups.push_back(NODE_SET_RELATION_GROUP {relation.inter, set_ids, 1});
        }
    }
    std::sort(relation_groups.begin(), relation_groups.end(),
              [](const NODE_SET_RELATION_GROUP &lhs, const NODE_SET_RELATION_GROUP &rhs) {
        return lhs.relation_count > rhs.relation_count;
    });
    time_print("%zu node sets generaeted.", relation_groups.size());
    //Create the node sets by merging the most common relation groups.
    time_print("Find the greatest no-intersection union...");
    std::vector<NODE_ID_SET> node_set_guess;
    {
        std::vector<NODE_ID_SET> relation_guess;
        relation_guess.reserve(relation_groups.size());
        for(auto relation_group: relation_groups)
        {
            NODE_ID_SET node_ids;
            for(int set_id: relation_group.set_ids)
            {
                const auto &node_id_set = node_set_marks[set_id].node_ids;
                node_ids.insert(node_id_set.begin(), node_id_set.end());
            }
            relation_guess.push_back(node_ids);
        }
        std::sort(relation_guess.begin(), relation_guess.end(),
                  [](const NODE_ID_SET &lhs, const NODE_ID_SET &rhs) {
            return lhs.size() > rhs.size();
        });
        //Filter the unique set.
        node_set_guess.reserve(relation_groups.size());
        for(auto node_set: relation_guess)
        {
            bool not_exist = true;
            for(auto existed_node_set: node_set_guess)
            {
                if(is_subset(existed_node_set, node_set))
                {
                    not_exist = false;
                    break;
                }
            }
            if(not_exist)
            {
                node_set_guess.push_back(node_set);
            }
        }
    }
    //Deep copy the node set guess.
    std::vector<NODE_ID_SET> node_set_pieces;
    node_set_pieces.reserve(node_set_guess.size());
    for(const auto &node_set: node_set_guess)
    {
        node_set_pieces.push_back(NODE_ID_SET(node_set.begin(), node_set.end()));
    }
    //Keep finding the common parts of the guessed core sets.
    while(has_cross_ids(node_set_pieces))
    {
        for(size_t i=0; i<node_set_pieces.size()-1; ++i)
        {
            auto &i_set = node_set_pieces[i];
            bool group_break = false;
            for(size_t j=i+1; j<node_set_pieces.size(); ++j)
            {
                //Update the new common sets.
                auto i_j_intersection = set_intersection(i_set, node_set_pieces[j]);
                if(!i_j_intersection.empty())
                {
                    //Break i, j into three set.
                    auto &j_set = node_set_pieces[j];
                    NODE_ID_SET i_j_common;
                    extract_set(i_set, j_set, i_j_common);
                    if(!i_j_common.empty())
                    {
                        if(j_set.empty())
                        {
                            node_set_pieces.erase(node_set_pieces.begin() + j);
                        }
                        if(i_set.empty())
                        {
                            node_set_pieces.erase(node_set_pieces.begin() + i);
                        }
                        node_set_pieces.push_back(i_j_common);
                        std::sort(node_set_pieces.begin(), node_set_pieces.end(), [](const NODE_ID_SET &lhs, const NODE_ID_SET &rhs) {
                            return lhs.size() > rhs.size();
                        });
                        group_break = true;
                        break;
                    }
                }
            }
            if(group_break)
            {
                break;
            }
        }
    }
    //Filtered the uncertained nodes.
    size_t uncertain_edge = 0;
    for(size_t i=0; i<node_set_pieces.size(); ++i)
    {
        if(node_set_pieces[i].size() == 1)
        {
            uncertain_edge = i;
            break;
        }
    }
    if(uncertain_edge > 0)
    {
        node_set_pieces.resize(uncertain_edge);
    }
    time_print("%zu group pieces generated.", node_set_pieces.size());
    //Now merged groups into target groups, always merge the minimum size into the larger groups.
    time_print("Merging into core groups...");
    std::vector<PIECES_BELONGS> pieces_belongs;
    pieces_belongs.reserve(node_set_pieces.size());
    for(size_t i=0; i<node_set_pieces.size(); ++i)
    {
        const auto &node_set_piece = node_set_pieces[i];
        std::unordered_set<int> belongs;
        for(size_t j=0; j<node_set_guess.size(); ++j)
        {
            if(is_subset(node_set_guess[j], node_set_piece))
            {
                belongs.insert(static_cast<int>(j));
            }
        }
        pieces_belongs.push_back(belongs);
    }
    //Merge the pieces into bigger groups.
    std::vector<PIECES_GROUP> pieces_groups;
    pieces_groups.reserve(pieces_belongs.size());
    for(size_t i=0; i<pieces_belongs.size(); ++i)
    {
        const auto &belongs = pieces_belongs[i];
        //Loop and check whether it has intersect with the nodes.
        bool group_new = true;
        for(size_t j=0; j<pieces_groups.size(); ++j)
        {
            auto &pieces_group = pieces_groups[j];
            if(has_intersection(pieces_group.belongs, belongs))
            {
                //Merge the belongs into groups.
                pieces_group.belongs.insert(belongs.begin(), belongs.end());
                pieces_group.pieces_ids.insert(static_cast<int>(i));
                group_new = false;
                break;
            }
        }
        //Check shall we insert a new group.
        if(group_new)
        {
            std::unordered_set<int> pieces_ids;
            pieces_ids.insert(static_cast<int>(i));
            pieces_groups.push_back(PIECES_GROUP {belongs, pieces_ids});
        }
    }
    pieces_groups.reserve(pieces_groups.size());
    //Merge the pieces into groups.
    std::vector<HANA_GROUP> core_groups;
    core_groups.reserve(pieces_groups.size());
    for(size_t i=0; i<pieces_groups.size(); ++i)
    {
        NODE_ID_SET core_group;
        for(int piece_id: pieces_groups[i].pieces_ids)
        {
            const auto &node_id_piece = node_set_pieces[piece_id];
            core_group.insert(node_id_piece.begin(), node_id_piece.end());
        }
        size_t length = 0;
        for(int node_id: core_group)
        {
            length += nodes[node_id].length;
        }
        core_groups.push_back(HANA_GROUP {core_group, length});
    }
    time_print("%zu core groups generated.", pieces_groups.size());
    //Keep merge to the no of group.
    time_print("Merged to target groups.", pieces_groups.size());
    while(core_groups.size() > static_cast<size_t>(no_of_group))
    {
        //Find the limitation.
        size_t relation_limit = core_groups[0].ids.size();
        for(size_t i=1; i<core_groups.size(); ++i)
        {
            if(core_groups[i].ids.size() < relation_limit)
            {
                relation_limit = core_groups[i].ids.size();
            }
        }
        relation_limit <<= 1;
        //Extract the top group.
        size_t best_i = 0, best_j = 0;
        double max_relation = -1.0;
        for(size_t i=0; i<core_groups.size()-1; ++i)
        {
            const auto &set_i = core_groups[i];
            for(size_t j=i+1; j<core_groups.size(); ++j)
            {
                double relation_i_j = group_relation(set_i.ids, core_groups[j].ids, node_infos, relation_limit) / static_cast<double>(hMin(set_i.length, core_groups[j].length));
                if(relation_i_j > max_relation)
                {
                    best_i = i;
                    best_j = j;
                    max_relation = relation_i_j;
                }
            }
        }
        //Merge group i and group j.
        core_groups[best_i].ids.insert(core_groups[best_j].ids.begin(), core_groups[best_j].ids.end());
        core_groups.erase(core_groups.begin() + best_j);
    }
    time_print("%zu groups generated.", core_groups.size());
    time_print("MARU Stage...");
    //Find all the rest of the nodes.
    std::unordered_set<int> used_nodes;
    for(size_t i=0; i<core_groups.size(); ++i)
    {
        used_nodes.insert(core_groups[i].ids.begin(), core_groups[i].ids.end());
    }
    std::unordered_set<int> unused_nodes;
    unused_nodes.reserve(node_size - used_nodes.size());
    for(size_t i=0; i<node_size; ++i)
    {
        if(!hInSet(static_cast<int>(i), used_nodes))
        {
            unused_nodes.insert(static_cast<int>(i));
        }
    }
    //Goes from top to bottom, find out the best matching group.
    while(!unused_nodes.empty())
    {
        int best_node_id = -1, best_group_id = -1;
        double best_weight = -1.0;
        //Predict all the nodes.
        for(const int node_id: unused_nodes)
        {
            //Predict the node belonging result.
            GROUP_PREDICTION predicts = node_predict(node_id, nodes, node_infos, core_groups);
            //Pick the best matching group.
            const auto &best_predict = predicts[0];
            //Assign the best prediction.
            if(best_predict.mark > best_weight)
            {
                best_node_id = node_id;
                best_group_id = best_predict.group_id;
                best_weight = best_predict.mark;
            }
        }
        //Pick the top one as the best node id.
        if(best_weight < 0.0)
        {
            //Bad things happens.
            break;
        }
        //Merged our choices.
        core_groups[best_group_id].ids.insert(best_node_id);
        core_groups[best_group_id].length += nodes[best_node_id].length;
        //Remove the node id from the set.
        unused_nodes.erase(best_node_id);
    }
    //Dump the data to output parameters.
    groups = std::vector<std::vector<int> >();
    groups.reserve(core_groups.size());
    for(auto core_group: core_groups)
    {
        groups.push_back(std::vector<int>(core_group.ids.begin(), core_group.ids.end()));
    }
    //Generaete unknown ids.
    unknown_ids = std::list<int>(unused_nodes.begin(), unused_nodes.end());
    //Free the memory.
    delete[] node_infos;
}
