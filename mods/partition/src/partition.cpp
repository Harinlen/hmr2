#include <cassert>
#include <queue>
#include <thread>

#include "hmr_bin_file.h"
#include "hmr_contig_graph.h"
#include "hmr_global.h"
#include "hmr_ui.h"

#include "partition.h"

typedef struct CONTIG_INFO
{
    size_t trust_pos;
    CONTIG_EDGE_MAP edge_map;
    const HMR_CONTIG* contig;
} CONTIG_INFO;

typedef std::unordered_map<uint64_t, CONTIG_ID_SET> EDGE_VOTER;

typedef struct EDGE_VOTERS
{
    HMR_EDGE edge;
    CONTIG_ID_SET supporters;
} EDGE_VOTERS;

typedef struct KERNEL_CANDIDATE
{
    CONTIG_ID_SET ids;
    int32_t subsets;
} KERNEL_CANDIDATE;

typedef std::vector<KERNEL_CANDIDATE> KERNEL_CANDIDATES;

typedef struct KERNEL_SET_INTERSECTION
{
    size_t set_a, set_b;
    CONTIG_ID_SET ids;
} KERNEL_SET_INTERSECTION;

typedef std::unordered_set<size_t> KERNEL_PIECE_BELONGS;

typedef struct KERNEL_SET_INTERSECTION_MARK
{
    CONTIG_ID_SET contig_ids;
    KERNEL_PIECE_BELONGS set_ids;
    int32_t relation_count;
} KERNEL_SET_INTERSECTION_MARK;

typedef struct TRUST_EDGE
{
    int start;
    int end;
    double weight;
} TRUST_EDGE;

bool has_same_belongs(const KERNEL_PIECE_BELONGS& x, const KERNEL_PIECE_BELONGS& y)
{
    for (size_t x_id : x)
    {
        if (hInSet(x_id, y))
        {
            return true;
        }
    }
    return false;
}

typedef struct KERNEL_PIECE_MERGE
{
    std::vector<size_t> kernel_piece_ids;
    KERNEL_PIECE_BELONGS belongs;
} KERNEL_PIECE_MERGE;

typedef struct HANA_GROUP
{
    CONTIG_ID_SET contig_ids;
    int64_t length;
} HANA_GROUP;

bool trust_edge_comp(const TRUST_EDGE& lhs, const TRUST_EDGE& rhs)
{
    return lhs.weight > rhs.weight;
}

bool contig_edge_comp(const CONTIG_EDGE& lhs, const CONTIG_EDGE& rhs)
{
    return lhs.weight > rhs.weight;
}

bool voter_comp(const EDGE_VOTERS& lhs, const EDGE_VOTERS& rhs)
{
    return lhs.supporters.size() > rhs.supporters.size();
}

bool kernel_set_comp(const KERNEL_CANDIDATE& lhs, const KERNEL_CANDIDATE& rhs)
{
    return lhs.subsets > rhs.subsets;
}

bool kernel_inter_comp(const KERNEL_SET_INTERSECTION& lhs, const KERNEL_SET_INTERSECTION& rhs)
{
    return lhs.ids.size() > rhs.ids.size();
}

bool kernel_inter_mark_comp(const KERNEL_SET_INTERSECTION_MARK& lhs, const KERNEL_SET_INTERSECTION_MARK& rhs)
{
    return lhs.relation_count > rhs.relation_count;
}

bool contig_id_set_comp(const CONTIG_ID_SET& lhs, const CONTIG_ID_SET& rhs)
{
    return lhs.size() > rhs.size();
}

bool is_subset(const CONTIG_ID_SET& parent, const CONTIG_ID_SET& child)
{
    if (child.size() > parent.size())
    {
        return false;
    }
    //Check each item one by one.
    for (const int32_t y : child)
    {
        if (!hInSet(y, parent))
        {
            return false;
        }
    }
    return true;
}

CONTIG_ID_SET get_intersection(const CONTIG_ID_SET& x, const CONTIG_ID_SET& y)
{
    std::list<int32_t> ids;
    for (int32_t y_id : y)
    {
        if (hInSet(y_id, x))
        {
            ids.push_back(y_id);
        }
    }
    return CONTIG_ID_SET(ids.begin(), ids.end());
}

void get_differ(CONTIG_ID_SET& x, CONTIG_ID_SET& y, const CONTIG_ID_SET& intersection)
{
    std::list<int32_t> nx, ny;
    for (int32_t x_id : x)
    {
        if (!hInSet(x_id, intersection))
        {
            nx.push_back(x_id);
        }
    }
    for (int32_t y_id : y)
    {
        if (!hInSet(y_id, intersection))
        {
            ny.push_back(y_id);
        }
    }
    x = CONTIG_ID_SET(nx.begin(), nx.end());
    y = CONTIG_ID_SET(ny.begin(), ny.end());
}

void partition_load_edges(const char* filepath, CONTIG_EDGES& edges)
{
    //Load the file.
    FILE* edge_file;
    if (!bin_open(filepath, &edge_file, "rb"))
    {
        time_error(-1, "Failed to read edge file %s", filepath);
    }
    //Read the number of edges.
    HMR_EDGE_WEIGHT edge_weight{};
    size_t edge_count;
    fread(&edge_count, sizeof(size_t), 1, edge_file);
    for (size_t i = 0; i < edge_count; ++i)
    {
        fread(&edge_weight, sizeof(HMR_EDGE_WEIGHT), 1, edge_file);
        //Increase the record.
        edges[edge_weight.edge.pos.start].push_back(CONTIG_EDGE{ edge_weight.edge.pos.end, edge_weight.weight });
        edges[edge_weight.edge.pos.end].push_back(CONTIG_EDGE{ edge_weight.edge.pos.start, edge_weight.weight });
    }
    //Read and parse the edges.
    fclose(edge_file);
}

void vote_edge(EDGE_VOTER& voter, int32_t x, int32_t y, int32_t host_node)
{
    uint64_t voting_key = hmr_graph_edge(x, y).data;
    auto edge_finder = voter.find(voting_key);
    if (edge_finder == voter.end())
    {
        std::unordered_set<int32_t> node_set;
        node_set.insert(host_node);
        voter.insert(std::make_pair(voting_key, node_set));
    }
    else
    {
        edge_finder->second.insert(host_node);
    }
}

std::vector<size_t> find_candidate_belongs(const CONTIG_ID_SET& id_set, const KERNEL_CANDIDATES& kernel_sets)
{
    std::list<size_t> result;
    for (size_t i = 0; i < kernel_sets.size(); ++i)
    {
        //Find out whether the id set is the subset of the current set.
        if (is_subset(kernel_sets[i].ids, id_set))
        {
            result.push_back(i);
        }
    }
    return std::vector<size_t>(result.begin(), result.end());
}

double core_group_relation(const CONTIG_ID_SET& lhs, const CONTIG_ID_SET& rhs, const CONTIG_INFO* contig_info, size_t edge_limit)
{
    //Loop for all edges in smaller group.
    std::priority_queue<double, std::vector<double>, std::less<double> > weight_heap;
    for (int32_t r_id : rhs)
    {
        const auto& edge_map = contig_info[r_id].edge_map;
        //size_t edge_counter = 0;
        for (int32_t l_id: lhs)
        {
            auto l_id_finder = edge_map.find(l_id);
            if (l_id_finder != edge_map.end())
            {
                if (weight_heap.size() < edge_limit)
                {
                    weight_heap.push(l_id_finder->second);
                }
                else
                {
                    if (l_id_finder->second > weight_heap.top())
                    {
                        weight_heap.pop();
                        weight_heap.push(l_id_finder->second);
                    }
                }
            }
        }
    }
    //Sum up all the weights in the heap.
    double relation = 0.0;
    while (!weight_heap.empty())
    {
        relation += weight_heap.top();
        weight_heap.pop();
    }
    return relation;
}

std::vector<HANA_GROUP> partition_hana(size_t gcn_window, const HMR_CONTIGS& contigs, const CONTIG_EDGES& edges, const CONTIG_INFO* contig_info, const CONTIG_ID_SET& best_contigs, size_t num_of_group)
{
    time_print("%zu - HANA stage start...", gcn_window);
    KERNEL_CANDIDATES kernel_candidate_sets;
    {
        std::vector<EDGE_VOTERS> voter_ids;
        {
            //Vote the edges when both sides are in the best contig set.
            time_print("%zu - Voting edges...", gcn_window);
            EDGE_VOTER edge_voter;
            for (int32_t contig_id : best_contigs)
            {
                const auto& contig_edges = edges[contig_id];
                size_t edge_boundary = hMin(gcn_window, contig_edges.size());
                for (size_t i = 0; i < edge_boundary; ++i)
                {
                    const auto& edge = contig_edges[i];
                    if (hInSet(edge.id, best_contigs))
                    {
                        vote_edge(edge_voter, contig_id, edge.id, contig_id);
                    }
                }
                for (size_t i = 0; i < edge_boundary - 1; ++i)
                {
                    const int32_t i_id = contig_edges[i].id;
                    if (hInSet(i_id, best_contigs))
                    {
                        for (size_t j = i + 1; j < edge_boundary; ++j)
                        {
                            const int32_t j_id = contig_edges[j].id;
                            if (hInSet(j_id, best_contigs))
                            {
                                vote_edge(edge_voter, i_id, j_id, contig_id);
                            }
                        }
                    }
                }
            }
            //If one edge is supported by many voters, these voters should come from the same group.
            time_print("%zu - Collecting contig voting sets...", gcn_window);
            //Gathering the voter groups based on the number of contigs.
            voter_ids.reserve(edge_voter.size());
            for (const auto& edge_vote_result : edge_voter)
            {
                HMR_EDGE vote_edge{};
                vote_edge.data = edge_vote_result.first;
                voter_ids.push_back(EDGE_VOTERS{ vote_edge, edge_vote_result.second });
            }
            std::sort(voter_ids.begin(), voter_ids.end(), voter_comp);
            time_print("%zu - %zu set(s) collected", gcn_window, voter_ids.size());
        }
        //Extract the trust edges node sets.
        time_print("%zu - Extract kernel candidate contig sets...", gcn_window);
        kernel_candidate_sets.reserve(voter_ids.size());
        for (const auto& edge_voter_info : voter_ids)
        {
            const auto& contig_set = edge_voter_info.supporters;
            auto parent_ids = find_candidate_belongs(contig_set, kernel_candidate_sets);
            if (parent_ids.empty())
            {
                //Add a new record in the kernel sets.
                kernel_candidate_sets.push_back(KERNEL_CANDIDATE{ contig_set, 0 });
            }
            else
            {
                for (const size_t set_id : parent_ids)
                {
                    ++kernel_candidate_sets[set_id].subsets;
                }
            }
        }
        kernel_candidate_sets.reserve(kernel_candidate_sets.size());
        std::sort(kernel_candidate_sets.begin(), kernel_candidate_sets.end(), kernel_set_comp);
        time_print("%zu - %zu kernel candidate set(s) generated.", gcn_window, kernel_candidate_sets.size());
    }
    std::vector<HANA_GROUP> core_groups;
    {
        //Find out the intersection of the kernel candidate sets.
        time_print("%zu - Calculating the intersection of the candidates...", gcn_window);
        std::vector<KERNEL_SET_INTERSECTION> candidate_relations;
        candidate_relations.reserve(hSquare(kernel_candidate_sets.size()));
        for (size_t i = 0; i < kernel_candidate_sets.size() - 1; ++i)
        {
            const auto& i_set = kernel_candidate_sets[i].ids;
            for (size_t j = i + 1; j < kernel_candidate_sets.size(); ++j)
            {
                auto intersection = get_intersection(i_set, kernel_candidate_sets[j].ids);
                //Only care about the intersection with more than 1 nodes.
                if (intersection.size() > 1)
                {
                    candidate_relations.push_back(KERNEL_SET_INTERSECTION{ i, j, intersection });
                }
            }
        }
        //Check relation size.
        if (candidate_relations.empty())
        {
            //Directly construct the core groups from candidate sets.
            time_print("%zu - no intersection set(s) found, exit.", gcn_window);
            return core_groups;
        }
        candidate_relations.reserve(candidate_relations.size());
        std::sort(candidate_relations.begin(), candidate_relations.end(), kernel_inter_comp);
        time_print("%zu - %zu intersection set(s) found.", gcn_window, candidate_relations.size());
        //Combine the intersection sets, make sure there is no subsets.
        time_print("%zu - Merging intersections...", gcn_window);
        std::vector<KERNEL_SET_INTERSECTION_MARK> intersection_marks;
        intersection_marks.reserve(candidate_relations.size());
        for (const auto& relation : candidate_relations)
        {
            //Search inside the intersection marks.
            bool find_parent = false;
            for (auto& intersection_mark : intersection_marks)
            {
                if (is_subset(intersection_mark.contig_ids, relation.ids))
                {
                    //Add the set id to relation group.
                    intersection_mark.set_ids.insert(relation.set_a);
                    intersection_mark.set_ids.insert(relation.set_b);
                    ++intersection_mark.relation_count;
                    find_parent = true;
                    break;
                }
            }
            if (!find_parent)
            {
                std::unordered_set<size_t> set_ids;
                set_ids.insert(relation.set_a);
                set_ids.insert(relation.set_b);
                intersection_marks.push_back(KERNEL_SET_INTERSECTION_MARK{ relation.ids, set_ids, 1 });
            }
        }
        intersection_marks.reserve(intersection_marks.size());
        std::sort(intersection_marks.begin(), intersection_marks.end(), kernel_inter_mark_comp);
        time_print("%zu - %zu intersection(s) filtered.", gcn_window, intersection_marks.size());
        //Merge the kernel pieces based on the set ids.
        time_print("%zu - Merging intersections based on the relations...", gcn_window);
        bool merge_performed = false;
        do
        {
            merge_performed = false;
            //Loop and find whether we can find intersections between the set ids.
            for (size_t i = 0; i < intersection_marks.size() - 1; ++i)
            {
                auto& i_inter = intersection_marks[i];
                const auto& i_set = intersection_marks[i].set_ids;
                //Check whether we have 
                for (size_t j = i + 1; j < intersection_marks.size(); ++j)
                {
                    auto& j_inter = intersection_marks[j];
                    //When they have any intersection, then merge it.
                    if (has_same_belongs(i_set, j_inter.set_ids))
                    {
                        //Merge j to i.
                        i_inter.contig_ids.insert(j_inter.contig_ids.begin(), j_inter.contig_ids.end());
                        i_inter.set_ids.insert(j_inter.set_ids.begin(), j_inter.set_ids.end());
                        //Erase the j.
                        intersection_marks.erase(intersection_marks.begin() + j);
                        //Update the perform position.
                        merge_performed = true;
                        break;
                    }
                }
                if (merge_performed)
                {
                    break;
                }
            }
        } while (merge_performed);
        //Keep merging the groups into target size.
        core_groups.reserve(intersection_marks.size());
        for (size_t i = 0; i < intersection_marks.size(); ++i)
        {
            //Calculate the length of the contigs.
            core_groups.push_back(HANA_GROUP{ intersection_marks[i].contig_ids, 0 });
        }
        time_print("%zu - %zu intersection(s) left.", gcn_window, core_groups.size());
    }
    //Keep merging to target groups.
    time_print("%zu - Merged to target group...", gcn_window);
    while (core_groups.size() > num_of_group)
    {
        //Find out the minimum size of the groups as the limitation.
        size_t edge_limit = core_groups[0].contig_ids.size();
        for (size_t i = 1; i < core_groups.size(); ++i)
        {
            if (core_groups[i].contig_ids.size() < edge_limit)
            {
                edge_limit = core_groups[i].contig_ids.size();
            }
        }
        //Find out the best matched related groups.
        size_t group_i = 0, group_j = 0;
        double max_relation = -1.0;
        for (size_t i = 0; i < core_groups.size() - 1; ++i)
        {
            const auto& i_ids = core_groups[i].contig_ids;
            for (size_t j = i + 1; j < core_groups.size(); ++j)
            {
                double i_j_relation = core_group_relation(i_ids, core_groups[j].contig_ids, contig_info, edge_limit) / static_cast<double>(hMin(core_groups[i].contig_ids.size(), core_groups[j].contig_ids.size()));
                if (i_j_relation > max_relation)
                {
                    group_i = i; group_j = j;
                    max_relation = i_j_relation;
                }
            }
        }
        //Merge group i and group j.
        core_groups[group_i].contig_ids.insert(core_groups[group_j].contig_ids.begin(), core_groups[group_j].contig_ids.end());
        core_groups.erase(core_groups.begin() + group_j);
    }
    time_print("%zu - HANA stage complete.", gcn_window);
    return core_groups;
}

typedef struct GROUP_PREDICT
{
    int group_id;
    double mark;
} GROUP_PREDICT;

typedef std::vector<GROUP_PREDICT> GROUP_PREDICTION;

GROUP_PREDICTION contig_predict(int32_t contig_id, size_t edge_limit, const CONTIG_EDGES& edges, const CONTIG_INFO* contig_info, const std::vector<HANA_GROUP> &core_groups)
{
    //Loop and generate the weight array for all core groups.
    GROUP_PREDICTION results;
    results.reserve(core_groups.size());
    size_t* edge_counter = new size_t[core_groups.size()];
    for (size_t i = 0; i < core_groups.size(); ++i)
    {
        edge_counter[i] = 0;
        results.push_back(GROUP_PREDICT{ static_cast<int>(i), -1.0 });
    }
    //Calculate the best matching edges.
    const auto& contig_edges = edges[contig_id];
    size_t trust_pos = contig_info[contig_id].trust_pos;
    for (size_t edge_id = 0; edge_id < contig_edges.size(); ++edge_id)
    {
        const auto& edge = contig_edges[edge_id];
        //Find if this edge belongs to any one of the group.
        for (size_t i = 0; i < core_groups.size(); ++i)
        {
            if (edge_counter[i] < edge_limit && hInSet(edge.id, core_groups[i].contig_ids))
            {
                double edge_w = edge.weight;
                //Cut the weight if the id is greater than trust pos.
                if (edge_id > trust_pos)
                {
                    edge_w /= 2.0;
                }
                //Initial the result.
                if (results[i].mark < 0)
                {
                    results[i].mark = 0.0;
                }
                results[i].mark += edge.weight;
                ++edge_counter[i];
                break;
            }
        }
        //Check whether all the edge counter reaches the limitation.
        size_t full_count = 0;
        for (size_t i = 0; i < core_groups.size(); ++i)
        {
            if (edge_counter[i] == edge_limit)
            {
                ++full_count;
            }
        }
        if (full_count == core_groups.size())
        {
            break;
        }
    }
    delete[] edge_counter;
    //Sort the result.
    std::sort(results.begin(), results.end(), [](const GROUP_PREDICT& lhs, const GROUP_PREDICT& rhs) {
        return lhs.mark > rhs.mark;
        });
    return results;
}

double partition_mark(const std::vector<HANA_GROUP>& core_groups)
{
    //Calculate the standard derivation of the nodes.
    double length_ave = 0.0, num_of_groups = static_cast<double>(core_groups.size());
    for (const auto &core_group : core_groups)
    {
        length_ave += static_cast<double>(core_group.length);
    }
    length_ave /= num_of_groups;
    //Calculate the variance.
    double length_var = 0.0;
    for (const auto& core_group : core_groups)
    {
        length_var += hSquare(static_cast<double>(core_group.length) - length_ave);
    }
    length_var /= num_of_groups;
    return sqrt(length_var);
}

typedef struct HANAMARU_PARAM
{
    size_t gcn_window;
    const HMR_CONTIGS& contigs;
    const CONTIG_EDGES& edges;
    const CONTIG_INFO* contig_info;
    const CONTIG_ID_SET& best_contigs;
    const size_t num_of_group;
    std::vector<CONTIG_ID_SET>* result;
    double* result_mark;
} HANAMARU_PARAM;

void partition_hanamaru(const HANAMARU_PARAM &param)
{
    //Unpack the parameters.
    const size_t& gcn_window = param.gcn_window;
    const HMR_CONTIGS& contigs = param.contigs;
    const CONTIG_EDGES& edges = param.edges;
    const CONTIG_INFO* contig_info = param.contig_info;
    const CONTIG_ID_SET& best_contigs = param.best_contigs;
    const size_t num_of_group = param.num_of_group;
    std::vector<CONTIG_ID_SET>* result = param.result;
    double* result_mark = param.result_mark;
    // -- HANA stage --
    auto core_groups = partition_hana(gcn_window, contigs, edges, contig_info, best_contigs, num_of_group);
    if (core_groups.empty())
    {
        *result_mark = -1.0;
        return;
    }
    // -- MARU stage --
    time_print("%zu - MARU stage start...", gcn_window);
    //Find all the rest of the nodes.
    CONTIG_ID_SET used_nodes, unused_nodes;
    for (const auto &core_group : core_groups)
    {
        used_nodes.insert(core_group.contig_ids.begin(), core_group.contig_ids.end());
    }
    //Loop and construct unused nodes.
    for (int32_t i = 0, num_of_contigs = static_cast<int32_t>(contigs.size()); i < num_of_contigs; ++i)
    {
        if (!hInSet(i, used_nodes))
        {
            unused_nodes.insert(i);
        }
    }
    time_print("%zu - %zu contig(s) need to be classified.", gcn_window, unused_nodes.size());
    while (!unused_nodes.empty())
    {
        int best_node_id = -1, best_group_id = -1;
        double best_weight = -1.0;
        //Calculate the group limitation.
        size_t edge_limit = core_groups[0].contig_ids.size();
        for (size_t i = 1; i < core_groups.size(); ++i)
        {
            edge_limit = hMin(edge_limit, core_groups[i].contig_ids.size());
        }
        //Predict all the nodes.
        for (const int32_t contig_id : unused_nodes)
        {
            //Predict the contig belongs.
            auto predicts = contig_predict(contig_id, edge_limit, edges, contig_info, core_groups);
            //Pick the best matching group.
            const auto& best_predict = predicts[0];
            //Assign the best prediction.
            if (best_predict.mark > best_weight)
            {
                best_node_id = contig_id;
                best_group_id = best_predict.group_id;
                best_weight = best_predict.mark;
            }
        }
        //Pick the top one as the best node id.
        if (best_weight < 0.0)
        {
            //Bad things happens.
            break;
        }
        //Merged our choices.
        core_groups[best_group_id].contig_ids.insert(best_node_id);
        core_groups[best_group_id].length += contigs[best_node_id].length;
        //Remove the node id from the set.
        unused_nodes.erase(best_node_id);
    }
    time_print("%zu - MARU stage complete, %zu node(s) are undefined.", gcn_window, unused_nodes.size());
    //Calculate the mark of the core groups.
    (*result).reserve(core_groups.size());
    for (auto& core_group : core_groups)
    {
        size_t length = 0;
        for (const int32_t contig_id : core_group.contig_ids)
        {
            length += contigs[contig_id].length;
        }
        core_group.length = length;
        (*result).push_back(core_group.contig_ids);
    }
    //Save the result.
    *result_mark = partition_mark(core_groups);
}

std::vector<CONTIG_ID_SET> partition_run(const HMR_CONTIGS& contigs, CONTIG_EDGES& edges, size_t num_of_group, const int32_t& threads)
{
    //If the group is 1, no need to seperate.
    if (num_of_group < 2)
    {
        //Just result a full id set.
        std::vector<CONTIG_ID_SET> result;
        CONTIG_ID_SET full_id_set;
        for (int32_t i = 0, i_max = static_cast<int32_t>(contigs.size()); i < i_max; ++i)
        {
            full_id_set.insert(i);
        }
        result.push_back(full_id_set);
        return result;
    }
    //Build the node information.
    size_t contig_size = contigs.size(), max_trust_pos = 0;
    CONTIG_INFO* contig_info = new CONTIG_INFO[contig_size];
    time_print("Searching for high-quality contigs relations...");
    CONTIG_ID_SET best_contigs;
    {
        // 1/4 total edges (expected) would contains half of the nodes.
        //It actually contains more than this, because the graph is sparse.
        size_t trust_edge_size = (contig_size * contig_size) >> 2;
        std::priority_queue<TRUST_EDGE, std::vector<TRUST_EDGE>, decltype(&trust_edge_comp)> trust_edge_heap(trust_edge_comp);
        //Build the contig information, and build the priority queue at the same time.
        for (size_t i = 0; i < contig_size; ++i)
        {
            //Update the contig map.
            contig_info[i].contig = &contigs[i];
            //Sort the edges.
            auto& contig_edges = edges[i];
            std::sort(contig_edges.begin(), contig_edges.end(), contig_edge_comp);
            //Build the edge map, finding trust edges.
            auto& contig_edge_map = contig_info[i].edge_map;
            for (const auto& edge: contig_edges)
            {
                //Insert to edge map.
                contig_edge_map.insert(std::make_pair(edge.id, edge.weight));
                //Append the edge information to trust edges.
                TRUST_EDGE candidate{ static_cast<int32_t>(i), edge.id, edge.weight };
                if (trust_edge_heap.size() < trust_edge_size)
                {
                    trust_edge_heap.push(candidate);
                }
                else
                {
                    if (trust_edge_heap.top().weight < edge.weight)
                    {
                        trust_edge_heap.pop();
                        trust_edge_heap.push(candidate);
                    }
                }
            }
            //Find out the trust position.
            contig_info[i].trust_pos = 0;
            if (contig_edges.size() > 2)
            {
                //Calculate the 1st deviation.
                std::vector<double> weight_d;
                weight_d.reserve(contig_edges.size());
                for (size_t i = 1; i < contig_edges.size(); ++i)
                {
                    weight_d.push_back(contig_edges[i].weight - contig_edges[i - 1].weight);
                }
                //Find the maximum position of the derivcative.
                contig_info[i].trust_pos = std::distance(weight_d.begin(), std::max_element(weight_d.begin(), weight_d.end())) + 1;
            }
            else
            {
                contig_info[i].trust_pos = 1;
            }
            max_trust_pos = hMax(max_trust_pos, contig_info[i].trust_pos);
        }
        //Reverse the heap into vector.
        size_t edge_size = trust_edge_heap.size(), edge_size_1 = edge_size - 1;
        std::vector<TRUST_EDGE> trust_edges;
        trust_edges.resize(trust_edge_heap.size());
        for (size_t i = 0; i < edge_size; ++i)
        {
            trust_edges[i] = trust_edge_heap.top();
            trust_edge_heap.pop();
        }
        time_print("%zu edge(s) gathered.", trust_edges.size());
        //Gathering the best nodes.
        time_print("Gathering high-quality contigs...");
        size_t expected_contig_size = contig_size >> 1, edge_id = 0;
        while (best_contigs.size() < expected_contig_size && edge_id < edge_size)
        {
            best_contigs.insert(trust_edges[edge_id].start);
            best_contigs.insert(trust_edges[edge_id].end);
            ++edge_id;
        }
        time_print("%zu contig(s) selected for kernel building.", best_contigs.size());
    }
    //Runs Hana-Maru algorithm for multiple times, find out the best voting edge range.
    size_t window_min = hMin(num_of_group, static_cast<size_t>(3)), 
        window_max = hMax(num_of_group, max_trust_pos);
    bool bouncing_detected = false;
    double result_mark = -1.0, *thread_marks = new double[threads];
    std::vector<CONTIG_ID_SET> result, *thread_results = new std::vector<CONTIG_ID_SET>[threads];
    std::thread* gcn_worker = new std::thread[threads];
    HANAMARU_PARAM param{ 0, contigs, edges, contig_info, best_contigs, num_of_group, NULL, NULL };
    for (size_t gcn_window = window_min; gcn_window <= window_max; gcn_window += threads)
    {
        //Parallel calculate the result.
        size_t threads_used = (gcn_window + threads > window_max) ? (window_max - gcn_window + 1) : threads;
        //Start up thread to calculate the result.
        for (size_t i = 0; i < threads_used; ++i)
        {
            param.gcn_window = gcn_window + i;
            param.result = &thread_results[i];
            param.result_mark = &thread_marks[i];
            gcn_worker[i] = std::thread(partition_hanamaru, param);
        }
        for (size_t i = 0; i < threads_used; ++i)
        {
            gcn_worker[i].join();
            if (result_mark < 0.0)
            {
                //Trust the result.
                result_mark = thread_marks[i];
                result = thread_results[i];
            }
            else
            {
                if (!bouncing_detected)
                {
                    //Check whether the result mark is decending.
                    if (thread_marks[i] < result_mark)
                    {
                        result_mark = thread_marks[i];
                        result = thread_results[i];
                    }
                    else
                    {
                        //We found the bouncing, exit running.
                        bouncing_detected = true;
                    }
                }
            }
        }
        if (bouncing_detected)
        {
            break;
        }
    }
    delete[] gcn_worker;
    delete[] thread_results;
    delete[] thread_marks;
    //Give back the result.
    return result;
}
