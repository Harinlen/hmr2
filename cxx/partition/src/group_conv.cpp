#include <algorithm>
#include <thread>
#include <cmath>
#include <vector>

#include "hmr_cgd_parse.h"
#include "hmr_ui.h"

#include "group_conv.h"

typedef struct THREAD_BLOCK
{
    int idx;
    int total;
} THREAD_BLOCK;

typedef struct NODE_INFO
{
    int id;
    std::vector<CONTIG_EDGE> trust_edges;
    std::unordered_map<int, double> edge;
} NODE_INFO;

typedef struct NODE_MARK
{
    NODE_INFO *info;
    double mark;
} NODE_MARK;

typedef struct HANA_GROUP
{
    std::vector<int> nodes;
    double mark;
} HANA_GROUP;

typedef struct MARU_GROUP
{
    std::set<int> nodes;
    std::unordered_map<int, double> trust_edges;
    double mark;
    std::unordered_map<MARU_GROUP *, double> relation;
} MARU_GROUP;

typedef struct MARU_GROUP_RELATION
{
    MARU_GROUP *left;
    MARU_GROUP *right;
    double mark;
} MARU_GROUP_RELATION;

typedef struct MARU_NODE_RELATION
{
    int group_id;
    double mark;
} MARU_NODE_RELATION;

template <typename T>
constexpr inline const T &hMin(const T &a, const T &b) { return (a < b) ? a : b; }

template <typename T>
constexpr inline const T hAbs(const T &x) { return (x < 0) ? (-1 * x) : x; }

void group_hana_trust_edge(const THREAD_BLOCK &block, CONTIG_NODE *nodes, size_t node_size, NODE_INFO *node_infos)
{
    //Decide the start and end based on thread ID.
    size_t range = (node_size + block.total - 1) / block.total;
    for(size_t i=range*block.idx, i_max=hMin(i+range, node_size); i<i_max; ++i)
    {
        //Sort the edges first.
        std::sort(nodes[i].links.begin(), nodes[i].links.end(), [](const auto &lhs, const auto &rhs){
            return lhs.weight > rhs.weight;
        });
        //Calculate the sum of middle value of the two side.
        const auto &node = nodes[i];
        const auto &edges = node.links;
        std::unordered_map<int, double> edge;
        for(size_t j=0; j<edges.size(); ++j)
        {
            edge.insert(std::make_pair(edges[j].id, edges[j].weight));
        }
        std::vector<CONTIG_EDGE> trust_edges;
        if(edges.size() > 2)
        {
            std::vector<double> mid_sum;
            mid_sum.reserve(edges.size());
            for(size_t pos=1; pos<edges.size()-1; ++pos)
            {
                //Calculate the total loss of the left hand side and right hand size.
                double mid_left = (edges[0].weight + edges[pos-1].weight) / 2.0,
                        mid_right = (edges[pos].weight + edges[edges.size()-1].weight) / 2.0;
                mid_sum.push_back(mid_left + mid_right);
            }
            //Calculate the derivative of the mid-sum.
            std::vector<double> deri_mid;
            deri_mid.reserve(mid_sum.size());
            for(size_t j=1; j<mid_sum.size(); ++j)
            {
                deri_mid.push_back(hAbs(mid_sum[j] - mid_sum[j-1]));
            }
            //Find the maximum position of the dericative.
            size_t trust_pos = std::distance(deri_mid.begin(), std::max_element(deri_mid.begin(), deri_mid.end())) + 1;
            //Prepare the trust edges.
            trust_edges.reserve(trust_pos);
            for(size_t j=0; j<trust_pos; ++j)
            {
                trust_edges.push_back(edges[j]);
            }
        }
        node_infos[i] = NODE_INFO {static_cast<int>(i), trust_edges, edge};
    }
}

double group_maru_relation(MARU_GROUP *m_group, MARU_GROUP *n_group, NODE_INFO *node_infos)
{
    double mark = 0.0;
    for(int m: m_group->nodes)
    {
        for(int n: n_group->nodes)
        {
            if(m < n)
            {
                const auto &edge_finder = node_infos[m].edge.find(n);
                if(edge_finder != node_infos[m].edge.end())
                {
                    mark += edge_finder->second;
                }
            }
            else
            {
                const auto &edge_finder = node_infos[n].edge.find(m);
                if(edge_finder != node_infos[n].edge.end())
                {
                    mark += edge_finder->second;
                }
            }
        }
    }
    return mark;
}

void group_maru_marks(const THREAD_BLOCK &block, MARU_GROUP **groups, size_t group_size, NODE_INFO *node_infos)
{
    //Decide the start and end based on thread ID.
    size_t range = (group_size + block.total - 1) / block.total;
    for(size_t i=range*block.idx, i_max=hMin(i+range, group_size); i<i_max; ++i)
    {
        //Calculate the best 4 mark of the group.
        double mark = 0.0;
        //Create a vector for the nodes.
        std::vector<int> node_vec(groups[i]->nodes.begin(), groups[i]->nodes.end());
        for(size_t m=0; m<node_vec.size(); ++m)
        {
            const auto &edges = node_infos[i].edge;
            for(size_t n=0; n<node_vec.size(); ++n)
            {
                const auto &edge_finder = edges.find(node_vec[n]);
                if(edge_finder != edges.end())
                {
                    mark += edge_finder->second;
                }
            }
        }
        groups[i]->mark = mark;
        //Calculate the relation mark.
        std::unordered_map<MARU_GROUP *, double> relation;
        relation.reserve(group_size - 1);
        for(size_t j=0; j<group_size; ++j)
        {
            //Ignore the relation to itself.
            if(j==i)
            {
                continue;
            }
            //Insert the mark to relation.
            relation.insert(std::make_pair(groups[j], group_maru_relation(groups[i], groups[j], node_infos)));
        }
        groups[i]->relation = relation;
    }
}

void group_apply_merge(const THREAD_BLOCK &block, MARU_GROUP** groups, size_t group_size, NODE_INFO *node_infos,
                       MARU_GROUP *group_left, MARU_GROUP *group_right)
{
    //Decide the start and end based on thread ID.
    size_t range = (group_size + block.total - 1) / block.total;
    for(size_t i=range*block.idx, i_max=hMin(i+range, group_size); i<i_max; ++i)
    {
        //Ignore the group right.
        if(group_right == groups[i])
        {
            continue;
        }
        //For all the groups, remove the group right record.
        groups[i]->relation.erase(group_right);
        //Update the value of the group left mark.
        if(group_left != groups[i])
        {
            auto left_finder = groups[i]->relation.find(group_left);
            left_finder->second = group_maru_relation(groups[i], group_left, node_infos);
        }
    }
}

double group_maru_mark(const std::vector<int> &right_nodes,
                       const std::vector<int> &left_nodes,
                       NODE_INFO *node_infos)
{
    double edge_1st = 0.0, edge_2nd = 0.0;
    for(auto n: right_nodes)
    {
        const auto &right_edges = node_infos[n].edge;
        for(auto m: left_nodes)
        {
            //Loop and find the edges.
            const auto &left_result = right_edges.find(m);
            if(left_result != right_edges.end())
            {
                if(left_result->second > edge_1st)
                {
                    edge_2nd = edge_1st;
                    edge_1st = left_result->second;
                }
                else if(left_result->second > edge_2nd)
                {
                    edge_2nd = left_result->second;
                }
                continue;
            }
            //Find in dst edges.
            const auto &left_edges = node_infos[m].edge;
            const auto &right_result = left_edges.find(n);
            if(right_result != left_edges.end())
            {
                if(right_result->second > edge_1st)
                {
                    edge_2nd = edge_1st;
                    edge_1st = right_result->second;
                }
                else if(right_result->second > edge_2nd)
                {
                    edge_2nd = right_result->second;
                }
            }
        }
    }
    return edge_1st + edge_2nd;
}

double group_maru_mark_4(const std::vector<int> &right_nodes,
                         const std::vector<int> &left_nodes,
                         NODE_INFO *node_infos)
{
    double edge_1 = 0.0, edge_2 = 0.0, edge_3 = 0.0, edge_4 = 0.0;
    for(auto n: right_nodes)
    {
        const auto &right_edges = node_infos[n].edge;
        for(auto m: left_nodes)
        {
            //Loop and find the edges.
            {
                const auto &left_result = right_edges.find(m);
                if(left_result != right_edges.end())
                {
                    if(left_result->second > edge_1)
                    {
                        edge_4 = edge_3;
                        edge_3 = edge_2;
                        edge_2 = edge_1;
                        edge_1 = left_result->second;
                    }
                    else if(left_result->second > edge_2)
                    {
                        edge_4 = edge_3;
                        edge_3 = edge_2;
                        edge_2 = left_result->second;
                    }
                    else if(left_result->second > edge_3)
                    {
                        edge_4 = edge_3;
                        edge_3 = left_result->second;
                    }
                    else if(left_result->second > edge_4)
                    {
                        edge_4 = left_result->second;
                    }
                    continue;
                }
            }
            //Find in dst edges.
            {
                const auto &left_edges = node_infos[m].edge;
                const auto &right_result = left_edges.find(n);
                if(right_result != left_edges.end())
                {
                    if(right_result->second > edge_1)
                    {
                        edge_4 = edge_3;
                        edge_3 = edge_2;
                        edge_2 = edge_1;
                        edge_1 = right_result->second;
                    }
                    else if(right_result->second > edge_2)
                    {
                        edge_4 = edge_3;
                        edge_3 = edge_2;
                        edge_2 = right_result->second;
                    }
                    else if(right_result->second > edge_3)
                    {
                        edge_4 = edge_3;
                        edge_3 = right_result->second;
                    }
                    else if(right_result->second > edge_4)
                    {
                        edge_4 = right_result->second;
                    }
                }
            }
        }
    }
    return edge_1 + edge_2 + edge_3 + edge_4;
}

//MARU_GROUP *group_maru_best_total(MARU_GROUP *group_src, NODE_INFO *nodes)
//{
//    (void)nodes;
//    //Based on the relation.
//    const std::unordered_map<MARU_GROUP *, double> &relation = group_src->relation;
//    MARU_GROUP *group = NULL;
//    double max_mark = 0.0;
//    //Find the maximum marks.
//    for(auto i: relation)
//    {
//        if(group == NULL || max_mark < i.second)
//        {
//            group = i.first;
//            max_mark = i.second;
//        }
//    }
//    return group;
//}

void group_maru_select_smallest(std::vector<MARU_GROUP *> &maru_group, NODE_INFO *node_infos, MARU_GROUP **group_left, MARU_GROUP **group_right)
{
    //Merged the most group from the least nodes.
    std::sort(maru_group.begin(), maru_group.end(), [](MARU_GROUP *lhs, MARU_GROUP *rhs) {
        return lhs->nodes.size() == rhs->nodes.size() ? lhs->mark > rhs->mark : lhs->nodes.size() < rhs->nodes.size();
    });
    MARU_GROUP *right = maru_group[0];
    MARU_GROUP *left = NULL;
    double max_mark = 0.0;
    //Use the two best relations border to the other.
    std::vector<int> right_nodes(right->nodes.begin(), right->nodes.end());
    for(auto i: right->relation)
    {
        //Calculate the group mark.
        MARU_GROUP *group_dst = i.first;
        double group_mark = group_maru_mark_4(right_nodes, std::vector<int>(group_dst->nodes.begin(), group_dst->nodes.end()), node_infos);
        if(group_mark > max_mark)
        {
            left = group_dst;
            max_mark = group_mark;
        }
    }
    //Set the value.
    *group_left = left;
    *group_right = right;
}

typedef struct MARU_GROUP_MARK
{
    MARU_GROUP *right;
    MARU_GROUP *left;
    double mark;
} MARU_GROUP_MARK;

void group_maru_select_all(std::vector<MARU_GROUP *> &maru_group, NODE_INFO *node_infos, MARU_GROUP **group_left, MARU_GROUP **group_right)
{
    //Merged the most group from the least nodes.
    std::sort(maru_group.begin(), maru_group.end(), [](MARU_GROUP *lhs, MARU_GROUP *rhs) {
        return lhs->nodes.size() == rhs->nodes.size() ? lhs->mark > rhs->mark : lhs->nodes.size() < rhs->nodes.size();
    });
    MARU_GROUP *right = maru_group[0];
    MARU_GROUP *left = NULL;
    double max_mark = 0.0;
    //Use the two best relations border to the other.
    std::vector<int> right_nodes(right->nodes.begin(), right->nodes.end());
    for(auto i: right->relation)
    {
        //Calculate the group mark.
        MARU_GROUP *group_dst = i.first;
        double group_mark = group_maru_mark_4(right_nodes, std::vector<int>(group_dst->nodes.begin(), group_dst->nodes.end()), node_infos);
        if(group_mark > max_mark)
        {
            left = group_dst;
            max_mark = group_mark;
        }
    }
    std::vector<MARU_GROUP_MARK> group_marks;
    group_marks.reserve(maru_group.size() * maru_group.size());
    for(size_t i=0; i<maru_group.size()-1; ++i)
    {
        MARU_GROUP *right = maru_group[i];
        std::vector<int> right_nodes(right->nodes.begin(), right->nodes.end());
        for(size_t j=i+1; j<maru_group.size(); ++j)
        {
            group_marks.push_back(MARU_GROUP_MARK
            {right,
             maru_group[j],
             group_maru_mark(right_nodes, std::vector<int>(maru_group[j]->nodes.begin(), maru_group[j]->nodes.end()), node_infos)});
        }
    }
    //Sort the group marks.
    std::sort(group_marks.begin(), group_marks.end(), [](const MARU_GROUP_MARK &lhs, const MARU_GROUP_MARK &rhs) {
        return lhs.mark > rhs.mark;
    });

    //Set the value.
    *group_left = left;
    *group_right = right;

}

typedef void (*MARU_GROUP_SELECT)(std::vector<MARU_GROUP *> &, NODE_INFO *, MARU_GROUP **, MARU_GROUP **);

void group_maru_reduce(MARU_GROUP_SELECT select_method, std::set<MARU_GROUP *> &group_set, std::vector<MARU_GROUP *> &maru_group, NODE_INFO *node_infos, int target_groups, int threads)
{
    while(static_cast<int>(group_set.size()) > target_groups)
    {
        //Select the merge groups.
        MARU_GROUP *group_left, *group_right;
        //Merged the most group from the least nodes.
        select_method(maru_group, node_infos, &group_left, &group_right);
        //Merge the group right to group left.
        group_left->nodes.insert(group_right->nodes.begin(), group_right->nodes.end());
        group_left->mark += group_right->mark + group_left->relation.find(group_right)->second;
        //Update all relations.
        {
            int update_threads = hMin(threads, static_cast<int>(group_set.size()));
            std::thread *update_workers = new std::thread[update_threads];
            for(int i=0; i<update_threads; ++i)
            {
                update_workers[i] = std::thread(
                            group_apply_merge, THREAD_BLOCK {i, update_threads},
                            maru_group.data(), maru_group.size(), node_infos,
                            group_left, group_right);
            }
            for(int i=0; i<update_threads; ++i)
            {
                update_workers[i].join();
            }
            delete[] update_workers;
        }
        //Remove the group right from the set and recover the memory.
        group_set.erase(group_right);
        delete group_right;
        //Update the maru group.
        maru_group = std::vector<MARU_GROUP *>(group_set.begin(), group_set.end());
    }
}

std::vector<std::set<int> > group_hanamaru(CONTIG_NODE *nodes, size_t node_size, int groups, int threads)
{
    time_print("Find trusted edges...");
    //Generate the node information.
    NODE_INFO *node_infos = new NODE_INFO[node_size];
    {
        std::thread *node_info_workers = new std::thread[threads];
        for(int i=0; i<threads; ++i)
        {
            node_info_workers[i] = std::thread(
                        group_hana_trust_edge,
                        THREAD_BLOCK{i, threads}, nodes, node_size, node_infos);
        }
        for(int i=0; i<threads; ++i)
        {
            node_info_workers[i].join();
        }
        delete[] node_info_workers;
    }
    time_print("Ranking nodes for HANA stage...");
    //Normalized node ranking.
    int mark_range = static_cast<int>(node_size) / groups / 2;
    if(mark_range < 1)
    {
        mark_range = 1;
    }
    int *marks = new int[mark_range];
    if(mark_range < 10)
    {
        //Directly use log marks.
        int mark = 1;
        for(int i=mark_range-1; i>-1; --i)
        {
            marks[i] = mark;
            mark <<= 1;
        }
    }
    else
    {
        //Use linear for all the other.
        int mark = 1;
        for(int i=mark_range-1; i>9; --i)
        {
            marks[i] = mark;
            ++mark;
        }
        for(int i=9; i>-1; --i)
        {
            marks[i] = mark;
            mark <<= 1;
        }
    }
    NODE_MARK *node_marks = new NODE_MARK[node_size];
    for(size_t i=0; i<node_size; ++i)
    {
        node_marks[i] = NODE_MARK{ &node_infos[i], 0 };
    }
    for(size_t i=0; i<node_size; ++i)
    {
        const std::vector<CONTIG_EDGE> &edges = nodes[i].links;
        //Only mark for the ranges.
        for(int j=0, j_max=hMin(mark_range, static_cast<int>(edges.size()));
            j<j_max; ++j)
        {
            node_marks[edges[j].id].mark += marks[j];
        }
    }
    //Sort the node marks.
    std::sort(node_marks, node_marks + node_size,
              [](const auto &lhs, const auto &rhs) {
        return lhs.mark > rhs.mark;
    });
    // -- HANA Stage --
    //Top level nodes should be the central part of the chromosome.
    time_print("HANA stage...");
    std::vector<HANA_GROUP> core_groups;
    {
        std::list<HANA_GROUP> core_group_list;
        for(size_t i=0; i<node_size; ++i)
        {
            const auto &node_mark = node_marks[i];
            if(node_mark.info->id == -1)
            {
                continue;
            }
            //Use the trust edges to grow the edge group.
            std::set<int> group_nodes;
            std::list<NODE_INFO *> info_queue;
            info_queue.push_back(node_mark.info);
            while(!info_queue.empty())
            {
                //Pop the info from the head.
                NODE_INFO *info = info_queue.front();
                info_queue.pop_front();
                //Check whether the node is already used.
                if(info->id == -1)
                {
                    continue;
                }
                //Push the node and related info.
                group_nodes.insert(info->id);
                //Reset the node info id.
                info->id = -1;
                //Loop for all the trust edges, add the info.
                for(auto &edge: info->trust_edges)
                {
                    //Push the node info to queue.
                    info_queue.push_back(&node_infos[edge.id]);
                }
            }
            //Calculate the group marks.
            std::vector<int> node_ids(group_nodes.begin(), group_nodes.end());
            double group_mark = 0.0;
            if(node_ids.size() > 1)
            {
                for(size_t n=0; n<node_ids.size()-1; ++n)
                {
                    const auto &edges = node_infos[node_ids[n]].edge;
                    for(size_t m=n+1; m<node_ids.size(); ++m)
                    {
                        auto edge_finder = edges.find(node_ids[m]);
                        if(edge_finder != edges.end())
                        {
                            group_mark += edge_finder->second;
                        }
                    }
                }
            }
            core_group_list.push_back(HANA_GROUP {node_ids, group_mark});
        }
        core_groups = std::vector<HANA_GROUP>(core_group_list.begin(),
                                              core_group_list.end());
    }
    std::sort(core_groups.begin(), core_groups.end(),
              [](const auto &lhs, const auto &rhs) {
        return lhs.nodes.size() > rhs.nodes.size();
    });
    time_print_size("HANA stage complete, find %zu groups.", core_groups.size());

    // -- MARU Stage --
    time_print("MARU stage...");
    //Find the groups which is contains more than 1 node.
    size_t valid_groups = 0;
    for(size_t i=0; i<core_groups.size(); ++i)
    {
        if(core_groups[i].nodes.size() < 2)
        {
            valid_groups = i;
            break;
        }
    }
    time_print_size("%zu groups are marked as valid groups.", valid_groups);
    //Combine those groups under the limits.
    std::vector<MARU_GROUP *> node_belongs;
    node_belongs.resize(node_size);
    for(size_t i=0; i<node_size; ++i)
    {
        node_belongs[i] = NULL;
    }
    //Based on the trusted edges, try to merge the top level groups.
    {
        time_print("Merging top level groups with trusted edges...");
        std::list<MARU_GROUP *> maru_queue;
        for(size_t i=0; i<valid_groups; ++i)
        {
            const auto &node_list = core_groups[i].nodes;
            MARU_GROUP *group = new MARU_GROUP();
            group->nodes = std::set<int>(node_list.begin(), node_list.end());
            const auto &node_set = group->nodes;
            for(auto node_id: node_list)
            {
                //Loop for all the edge.
                for(const auto &edge: node_infos[node_id].trust_edges)
                {
                    if(node_set.find(edge.id) == node_set.end())
                    {
                        group->trust_edges.insert(
                                    std::make_pair(edge.id, edge.weight));
                    }
                }
                node_belongs[node_id] = group;
            }
            maru_queue.push_back(group);
        }
        //Try to merge the group based on trust edges.
        while(!maru_queue.empty())
        {
            MARU_GROUP *group = maru_queue.front();
            maru_queue.pop_front();
            //Check whether the group could be merged with existed groups.
            if(group->trust_edges.empty())
            {
                continue;
            }
            for(auto i: group->trust_edges)
            {
                MARU_GROUP *target_group = node_belongs[i.first];
                if(target_group != NULL && target_group != group)
                {
                    //Change all the target node group to current group.
                    for(auto j: target_group->nodes)
                    {
                        node_belongs[j] = group;
                    }
                    //Combine the set.
                    group->nodes.insert(target_group->nodes.begin(),
                                        target_group->nodes.end());
                    const auto &node_set = group->nodes;
                    //Construct the trust edges.
                    std::unordered_map<int, double> trust_edges;
                    for(auto j: group->trust_edges)
                    {
                        if(node_set.find(j.first) == node_set.end())
                        {
                            trust_edges.insert(j);
                        }
                    }
                    for(auto j: target_group->trust_edges)
                    {
                        if(node_set.find(j.first) == node_set.end())
                        {
                            trust_edges.insert(j);
                        }
                    }
                    group->trust_edges = trust_edges;
                    //Merge target group to group.
                    maru_queue.remove(target_group);
                    delete target_group;
                    //Push the group back to queue.
                    maru_queue.push_front(group);
                    break;
                }
            }
        }
    }
    //Keep merging until the group reaches the group limitations.
    std::set<MARU_GROUP *> group_set;
    for(size_t i=0; i<node_size; ++i)
    {
        if(node_belongs[i] != NULL)
        {
            group_set.insert(node_belongs[i]);
        }
    }
    time_print_size("%zu groups left after merge.", group_set.size());
    std::vector<MARU_GROUP *> maru_group(group_set.begin(), group_set.end());
    //Now calculate the relationship of the groups.
    {
        time_print("Building relation matrix of the groups...");
        std::thread *mark_workers = new std::thread[threads];
        for(int i=0; i<threads; ++i)
        {
            mark_workers[i] = std::thread(group_maru_marks,
                                          THREAD_BLOCK{i, threads}, maru_group.data(), maru_group.size(), node_infos);
        }
        for(int i=0; i<threads; ++i)
        {
            mark_workers[i].join();
        }
        delete[] mark_workers;
    }
    //Loop and merge the group until the expected groups.
    time_print("Merging groups based on length...");
    group_maru_reduce(group_maru_select_smallest, group_set, maru_group, node_infos, groups << 1, threads);
    time_print_size("%zu core groups left after merge.", group_set.size());
    time_print("Merging groups based on relations...");
    group_maru_reduce(group_maru_select_all, group_set, maru_group, node_infos, groups, threads);
    time_print_size("%zu core groups left after merge.", group_set.size());
    //Calculate the marks of the rest nodes to these core groups.
    time_print("Merging individual nodes into core groups...");
    MARU_NODE_RELATION *relation = new MARU_NODE_RELATION[groups];
    for(size_t i=valid_groups; i<core_groups.size(); ++i)
    {
        const int node_id = core_groups[i].nodes[0];
        const auto &node_edges = node_infos[node_id].edge;
        //Calculate the marks to all the core groups.
        for(int j=0; j<groups; ++j)
        {
            double target_mark = 0.0;
            for(int target_id: maru_group[j]->nodes)
            {
                const auto &target_result = node_edges.find(target_id);
                if(target_result != node_edges.end())
                {
                    target_mark += target_result->second;
                }
                else
                {
                    const auto &target_edges = node_infos[target_id].edge;
                    const auto &node_result = target_edges.find(node_id);
                    if(node_result != target_edges.end())
                    {
                        target_mark += target_result->second;
                    }
                }
            }
            relation[j] = MARU_NODE_RELATION {j, target_mark};
        }
        //Sort the marks.
        std::sort(relation, relation+groups,
                  [](const auto &lhs, const auto &rhs){
            return lhs.mark > rhs.mark;
        });
        //Chose the top group to merge into.
        maru_group[relation[0].group_id]->nodes.insert(node_id);
    }
    time_print_size("Clustering complete, %zu groups generated.", maru_group.size());
    //Construct the maru group sets.
    std::vector<std::set<int> > cluster_result;
    cluster_result.reserve(groups);
    for(auto i: maru_group)
    {
        cluster_result.push_back(i->nodes);
        delete i;
    }
    delete[] node_marks;
    delete[] marks;
    delete[] node_infos;
    //Provide the cluster result.
    return cluster_result;
}
