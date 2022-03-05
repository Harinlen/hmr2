#include <algorithm>
#include <queue>
#include <thread>
#include <cmath>

#include "hmr_global.h"
#include "hmr_cgd_parse.h"
#include "hmr_ui.h"

#include "partition.h"

typedef struct HANA_GROUP
{
    std::unordered_set<int> ids;
    size_t length;
} HANA_GROUP;
typedef std::unordered_map<int, double> EDGE_MAP;

typedef struct NODE_INFO
{
    size_t trust_pos;
    EDGE_MAP edge_map;
    HANA_GROUP *parent;
    CONTIG_NODE *node;
} NODE_INFO;

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

inline bool in_group(int node_id, HANA_GROUP *group)
{
    return hInSet(node_id, group->ids);
}

inline double edge_weight(const EDGE_MAP &map, int node_id)
{
    const auto finder = map.find(node_id);
    return finder == map.end() ? 0.0 : finder->second;
}

void fit_linear(const std::vector<double> &x, const std::vector<double> &y, double &k, double &b)
{
    size_t n = x.size();
    //If we only have two points, we directly solve it.
    if(n == 2)
    {
        k = (y[0] - y[1]) / (x[0] - x[1]);
        b = y[0] - k * x[0];
        return;
    }
    double x_ave = 0.0, y_ave = 0.0, xy_sum = 0.0, x2_sum = 0.0;
    for(size_t i=0; i<n; ++i)
    {
        x_ave += x[i];
        y_ave += y[i];
        xy_sum += x[i] * y[i];
        x2_sum += x[i] * x[i];
    }
    x_ave /= static_cast<double>(n);
    y_ave /= static_cast<double>(n);
    //Calculate k and b with least squares method.
    k = (xy_sum - n * x_ave * y_ave) / (x2_sum - n * x_ave * x_ave);
    b = y_ave - k * x_ave;
}

double linear_loss(const std::vector<double> &x, const std::vector<double> &y, const double &k, const double &b)
{
    size_t n = x.size();
    double loss = 0.0, A = k, B = -1, C = b;
    double bottom = sqrt(A * A + B * B);
    for(size_t i=0; i<n; ++i)
    {
        loss += hAbs(A * x[i] + B * y[i] + C) / bottom;
    }
    return loss;
}

void generate_mark(int **marks, int mark_range)
{
    int *mark_arr = new int[mark_range];
    int mark = 1;
    //The last 10 marks are exponential, others are linear.
    if(mark_range < 10)
    {
        for(int i=mark_range - 1; i>-1; --i)
        {
            mark_arr[i] = mark;
            mark <<= 1;
        }
    }
    else
    {
        for(int i=mark_range - 1; i>9; --i)
        {
            mark_arr[i] = mark;
            ++mark;
        }
        for(int i=9; i>-1; --i)
        {
            mark_arr[i] = mark;
            mark <<= 1;
        }
    }
    *marks = mark_arr;
}

void node_merge(HANA_GROUP *dst, int node_id, NODE_INFO *node_infos)
{
    //Update the node information.
    node_infos[node_id].parent = dst;
    //Update the node group.
    dst->ids.insert(node_id);
    dst->length += node_infos[node_id].node->length;
}

void group_merge(HANA_GROUP *dst, HANA_GROUP *src, NODE_INFO *node_infos)
{
    //Update all the parent of source group to target group.
    dst->ids.insert(src->ids.begin(), src->ids.end());
    for(int src_id: src->ids)
    {
        node_infos[src_id].parent = dst;
    }
    dst->length += src->length;
}

double group_relation(HANA_GROUP *lhs, HANA_GROUP *rhs, NODE_INFO *node_info, size_t edge_limit)
{
    HANA_GROUP *small = lhs, *big = rhs;
    if(rhs->ids.size() < lhs->ids.size())
    {
        small = rhs;
        big = lhs;
    }
    //Loop for all edges in smaller group.
    std::priority_queue<double, std::vector<double>, std::less<double> > weight_queue;
    //Find the first matched ids from the big group.
    for(int s_id: small->ids)
    {
        //The edges are already sorted descending, the first appeared edge is
        //the best match edge.
        const auto &edges = node_info[s_id].edge_map;
        for(const auto &edge: edges)
        {
            if(in_group(edge.first, big))
            {
                //When less than heap size, just insert.
                if(weight_queue.size() < edge_limit)
                {
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

void group_compare(size_t i, size_t j, const std::vector<HANA_GROUP *> &groups, NODE_INFO *node_info, size_t edge_limit,
                   size_t &max_i, size_t &max_j, double &max_relation)
{
    double relation = group_relation(groups[i], groups[j], node_info, edge_limit);
    if(relation > max_relation)
    {
        max_i = i;
        max_j = j;
        max_relation = relation;
    }
}

void thread_group_max_relation(const THREAD_BLOCK &block,
                               const std::vector<HANA_GROUP *> &core_groups, NODE_INFO *node_info, size_t edge_limit,
                               size_t *result_i, size_t *result_j, double *result_relation)
{
    //Calculate the start nodes.
    size_t thread_range = (core_groups.size() + block.total - 1) / block.total;
    size_t node_start = thread_range * block.idx,
            node_end = hMin(node_start + thread_range, core_groups.size());
    size_t max_i = 0, max_j = 0;
    double max_relation = 0.0;
    //Complete the thread workload.
    for(size_t i=node_start; i<node_end; ++i)
    {
        for(size_t j=i+1; j<core_groups.size(); ++j)
        {
            group_compare(i, j, core_groups, node_info, edge_limit, max_i, max_j, max_relation);
        }
    }
    //Set the answer of the part.
    *result_i = max_i;
    *result_j = max_j;
    *result_relation = max_relation;
}

typedef struct GROUP_PREDICT
{
    int group_id;
    double mark;
} GROUP_PREDICT;

bool group_predict_compare(const GROUP_PREDICT &lhs, const GROUP_PREDICT &rhs)
{
    return lhs.mark > rhs.mark;
}

std::vector<GROUP_PREDICT> node_predict(
        int node_id, CONTIG_NODE *nodes, NODE_INFO *node_infos,
        const std::vector<HANA_GROUP *> &core_groups)
{
    //Find the edge limits.
    size_t edge_limit = 0;
    for(const HANA_GROUP *group: core_groups)
    {
        if(edge_limit == 0 || group->ids.size())
        {
            edge_limit = group->ids.size();
        }
    }
    //Loop and generate the weight array for all core groups.
    std::vector<GROUP_PREDICT> result;
    result.reserve(core_groups.size());
    size_t *edge_counter = new size_t[core_groups.size()];
    for(size_t i=0; i<core_groups.size(); ++i)
    {
        edge_counter[i] = 0;
        result.push_back(GROUP_PREDICT {static_cast<int>(i), 0.0});
    }
    //From the biggest edges to minimum groups.
    const auto &edges = nodes[node_id].links;
    size_t trust_pos = node_infos[node_id].trust_pos;
    for(size_t edge_id=0; edge_id < edges.size(); ++edge_id)
    {
        const auto &edge = edges[edge_id];
        for(size_t i=0; i<core_groups.size(); ++i)
        {
            if(edge_counter[i] < edge_limit && in_group(edge.id, core_groups[i]))
            {
                double edge_w = edge.weight;
                if(edge_id > trust_pos)
                {
                    edge_w /= 2.0;
                }
                result[i].mark += edge.weight;
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
    std::sort(result.begin(), result.end(), group_predict_compare);
    //Generate the group predict result.
    return result;
}

bool dual_linear_detect(const std::vector<double> &marks)
{
    if(marks.size() < 2)
    {
        return true;
    }
    double best_left_k, best_left_b, best_right_k, best_right_b, best_loss = -1.0;
    //Generate the x and y groups.
    for(size_t i=1, i_max = hMin(4, static_cast<int>(marks.size()>>1)); i<=i_max; ++i)
    {
        size_t left_bound = i + 1;
        //Fit the left side.
        double left_k, left_b, right_k, right_b;
        std::vector<double> left_x, left_y;
        for(size_t j=0; j<left_bound; ++j)
        {
            left_x.push_back(static_cast<double>(j));
            left_y.push_back(marks[j]);
        }
        fit_linear(left_x, left_y, left_k, left_b);
        std::vector<double> right_x, right_y;
        for(size_t j=i; j<marks.size(); ++j)
        {
            right_x.push_back(static_cast<double>(j));
            right_y.push_back(marks[j]);
        }
        fit_linear(right_x, right_y, right_k, right_b);
        //Calculate the total loss.
        double total_loss = linear_loss(left_x, left_y, left_k, left_b) +
                linear_loss(right_x, right_y, right_k, right_b);
        if(best_loss < 0.0 || total_loss < best_loss)
        {
            best_left_k = left_k;
            best_left_b = left_b;
            best_right_k = right_k;
            best_right_b = right_b;
            best_loss = total_loss;
        }
    }
    //Find the cross point X.
    double x = (best_right_b - best_left_b) / (best_left_k - best_right_k);
//    printf("%lf\n", x);
//    for(double k: marks)
//    {
//        printf("%lf\t", k);
//    }
//    printf("\n");
    //It should close to 1.0. Which is the first point.
    return x < 1.5;
}

typedef struct NODE_MARK
{
    int id;
    size_t mark;
} NODE_MARK;

void group_hanamaru(CONTIG_NODE *nodes, size_t node_size, int no_of_group, int threads,
                    std::vector<std::vector<int> > &groups, std::list<int> &unknown_ids)
{
    time_print("Find trusted edges...");
    NODE_INFO *node_infos = new NODE_INFO[node_size];
    //Initial the node information.
    for(size_t i=0; i<node_size; ++i)
    {
        //Initial the parent.
        node_infos[i].parent = NULL;
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
    //Normalized node ranking.
    time_print("Ranking nodes for HANA stage...");
    int mark_range = static_cast<int>(node_size) / no_of_group / 2;
    int *marks;
    generate_mark(&marks, mark_range);
    NODE_MARK *node_marks = new NODE_MARK[node_size];
    for(size_t i=0; i<node_size; ++i)
    {
        node_marks[i] = NODE_MARK {static_cast<int>(i), 0};
    }
    for(size_t i=0; i<node_size; ++i)
    {
        const auto &edges = nodes[i].links;
        int node_mark_range = hMin(mark_range, static_cast<int>(edges.size()));
        for(int j=0; j<node_mark_range; ++j)
        {
            node_marks[edges[j].id].mark += marks[j];
        }
    }
    //Sort the nodes as mark.
    std::sort(node_marks, node_marks+node_size,
              [](const NODE_MARK &lhs, const NODE_MARK &rhs) {
        return lhs.mark > rhs.mark;
    });


    // -- HANA Stage --
    //Top level nodes should be the central part of the chromosome.
    time_print("HANA stage start...");
    std::vector<HANA_GROUP *> chr_groups;
    {
        std::unordered_set<int> hana_nodes;
        //Hyper-parameter: hana_size
        //Determines the central part.
        size_t hana_size = node_size / 2;
//        for(size_t i=0; i<hana_size; ++i)
//        {
//            int node_id = node_marks[i].id;
//            printf("%s\t", nodes[node_id].name);
//            auto &es = nodes[node_id].links;
//            for(int j=0; j<8; ++j)
//            {
//                auto &e = es[j];
//                printf("%s %lf", nodes[e.id].name, e.weight);
//                if(node_infos[node_id].trust_pos == j+1)
//                {
//                    printf("|");
//                }
//                printf("\t");
//            }
//            printf("\n");
//        }
//        exit(-1);
        time_print("Classify best ranking %zu node(s)...", hana_size);
        for(size_t i=0; i<hana_size; ++i)
        {
            hana_nodes.insert(node_marks[i].id);

        }
        //Based on the trust edges, build the core groups.
        //Classify the hana nodes.
        size_t group_counter = 0;
        for(size_t i=0; i<hana_size; ++i)
        {
            const int node_id = node_marks[i].id;
            //Only check the nodes whose parent is NULL.
            if(node_infos[node_id].parent)
            {
                continue;
            }
            //Create the group for the node.
            HANA_GROUP *group = new HANA_GROUP();
            group->length = 0;
            ++group_counter;
            //Search all the trust edges.
            std::list<int> search_id_queue;
            search_id_queue.push_back(node_id);
            while(!search_id_queue.empty())
            {
                //Pop the search id.
                int head_id = search_id_queue.front();
                search_id_queue.pop_front();
                //Check whether this node is already added to the group.
                if(in_group(head_id, group))
                {
                    continue;
                }
                //Push head id into group.
                node_merge(group, head_id, node_infos);
                //Update the node info.
                auto &head_info = node_infos[head_id];
                auto &head_edges = nodes[head_id].links;
                //Insert the trust edge node ids into group.
                for(size_t j=0; j<head_info.trust_pos; ++j)
                {
                    const int edge_id = head_edges[j].id;
                    //The target node is HANA node.
                    if(node_infos[edge_id].parent == NULL && //The node is not in any group.
//                            hInSet(edge_id, hana_nodes) &&   //The node must be a trusted nodes.
                            !in_group(edge_id, group))       //The node is not in the current group.
                    {
                        if(!hInSet(edge_id, hana_nodes))
                        {
                            hana_nodes.insert(edge_id);
                        }
                        search_id_queue.push_back(edge_id);
                    }
                }
            }
        }
        time_print("First round complete, %zu group(s) generated.", group_counter);
        time_print("Keep merging trust edges...");
        //Loop until all the trust edges are combined.
        bool is_merged = true;
        while(is_merged)
        {
            is_merged = false;
            for(size_t i=0; i<hana_size; ++i)
            {
                //Check from the start.
                const int node_id = node_marks[i].id;
                const auto &edges = nodes[node_id].links;
                auto &node_info = node_infos[node_id];
                //Check for all the trust edges.
                for(size_t j=0; j<node_info.trust_pos; ++j)
                {
                    const int edge_id = edges[j].id;
                    if(hInSet(edge_id, hana_nodes) && node_infos[edge_id].parent != node_info.parent)
                    {
                        //Merge two group into a single group.
                        group_merge(node_info.parent, node_infos[edge_id].parent, node_infos);
                        is_merged = true;
                        break;
                    }
                }
                if(is_merged)
                {
                    break;
                }
            }
        }
        //Summarize the groups and nodes.
        std::unordered_set<HANA_GROUP *> core_groups;
        for(size_t i=0; i<hana_size; ++i)
        {
            core_groups.insert(node_infos[node_marks[i].id].parent);
        }
        time_print("%zu group(s) generated.", core_groups.size());
        // ---- Definitely correct above ----

        // ---- Experiment codes ----
//        for(auto t: core_groups)
//        {
//            printf("-> %zu", t->ids.size());
//            size_t s=0;
//            for(auto nid: t->ids)
//            {
//                s += nodes[nid].length;
//            }
//            printf(" size: %zu\n", s);
//            for(auto nid: t->ids)
//            {
//                printf("%s\n", nodes[nid].name);
//            }
//            printf("--------\n");
//        }
//        for(int node_id: singles)
//        {
//            printf("%s\t", nodes[node_id].name);
//            auto &es = nodes[node_id].links;
//            for(int j=0; j<8; ++j)
//            {
//                auto &e = es[j];
//                printf("%s %lf %zu", nodes[e.id].name, e.weight, node_infos[e.id].parent);
//                if(node_infos[node_id].trust_pos == j+1)
//                {
//                    printf("|");
//                }
//                printf("\t");
//            }
//            printf("\n");
//        }
//        exit(-1);
        //Keep merge the groups until to the number of groups.
        time_print("Generating core groups...");
        size_t *max_is = new size_t[threads], *max_js = new size_t[threads];
        double *max_relations = new double[threads];
        //Generate workers.
        std::thread *workers = new std::thread[threads];
        while(core_groups.size() > static_cast<size_t>(no_of_group))
        {
            std::vector<HANA_GROUP *> core_group_list(core_groups.begin(), core_groups.end());
            //Use the size of minimum nodes group as the limitation.
            size_t edge_limit = 0;
            for(size_t i=0; i<core_group_list.size(); ++i)
            {
                if(edge_limit == 0 || core_group_list[i]->ids.size() < edge_limit)
                {
                    edge_limit = core_group_list[i]->ids.size();
                    //Reach the least minimum, ignore all the others.
                    if(edge_limit == 2)
                    {
                        break;
                    }
                }
            }
//            edge_limit <<= 1;

            //Find the best matched group.
            size_t max_i, max_j;
            double max_relation = -1.0;
            {
//                std::sort(core_group_list.begin(), core_group_list.end(), [](HANA_GROUP *left, HANA_GROUP *right)
//                {
//                    return left->ids.size() < right->ids.size();
//                });
                std::sort(core_group_list.begin(), core_group_list.end(), [](HANA_GROUP *left, HANA_GROUP *right)
                {
                    return left->length < right->length;
                });

                for(size_t j=1; j<core_group_list.size(); ++j)
                {
                    group_compare(0, j, core_group_list, node_infos, edge_limit,
                                  max_i, max_j, max_relation);
                }
            }

//            {
//                printf("Best matching!!\n");
                // -- Single thread --
//                {
//                    for(size_t i=0; i<core_group_list.size() - 1; ++i)
//                    {
//                        for(size_t j=i+1; j<core_group_list.size(); ++j)
//                        {
//                            group_compare(i, j, core_group_list, node_infos, edge_limit,
//                                          max_i, max_j, max_relation);
//                        }
//                    }
//                }
//                // -- Single thread end --
//                //            // -- Multi thread --
//                //            {
//                //                //Prepare the result area.

//                //                for(int i=0; i<threads; ++i)
//                //                {
//                //                    max_relations[i] = 0.0;
//                //                    workers[i] = std::thread(thread_group_max_relation, THREAD_BLOCK{i, threads}, core_group_list, node_infos, edge_limit,
//                //                                             &max_is[i], &max_js[i], &max_relations[i]);
//                //                }
//                //                for(int i=0; i<threads; ++i)
//                //                {
//                //                    workers[i].join();
//                //                    if(max_relations[i] > max_relation)
//                //                    {
//                //                        max_i = max_is[i];
//                //                        max_j = max_js[i];
//                //                        max_relation = max_relations[i];
//                //                    }
//                //                }
//                //            }
//            }
            // -- Multi thread end --
            //Merge group[max_j] to group[max_i].
//            printf("Merge: (limits: %zu, weight: %lf) nodes: %zu + %zu = %zu\tlength: %zu + %zu = %zu)\n",
//                   edge_limit, max_relation,
//                   core_group_list[max_j]->ids.size(), core_group_list[max_i]->ids.size(), core_group_list[max_j]->ids.size()+core_group_list[max_i]->ids.size(),
//                   core_group_list[max_j]->length, core_group_list[max_i]->length, core_group_list[max_j]->length+core_group_list[max_i]->length);
//            for(int nid: core_group_list[max_j]->ids)
//            {
//                printf("%s\t%d\n", nodes[nid].name, nid);
//            }
//            printf("to:\n");
//            for(int nid: core_group_list[max_i]->ids)
//            {
//                printf("%s\t%d\n", nodes[nid].name, nid);
//            }
            core_groups.erase(core_group_list[max_j]);
            group_merge(core_group_list[max_i], core_group_list[max_j], node_infos);
//            printf("==============\n");
//            time_print_size("%zu groups remain.", core_groups.size());
        }
        delete[] workers;
        delete[] max_is;
        delete[] max_js;
        delete[] max_relations;
        //Merge the rest of the trust nodes to core groups.
        chr_groups = std::vector<HANA_GROUP *>(core_groups.begin(), core_groups.end());
        time_print("%zu core groups generated.", chr_groups.size());
    }

    // -- MARU Stage --
    time_print("Classifying individual nodes...");
    for(size_t i=0; i<node_size; ++i)
    {
        if(node_infos[i].parent == NULL)
        {
            //Predict the group of node it should be.
            std::vector<GROUP_PREDICT> predict = node_predict(static_cast<int>(i), nodes, node_infos, chr_groups);
            //Check whether the prediction result is correct.
            std::vector<double> marks;
            marks.reserve(predict.size());
            for(const auto &result: predict)
            {
                marks.push_back(result.mark);
            }
            if(dual_linear_detect(marks))
            {
                //Apply the predicted merge action.
                const auto &action = predict[0];
                node_merge(chr_groups[action.group_id], i, node_infos);
            }
            else
            {
                //This node is detected as defected node.
                time_print("Ambiguous contig '%s' is marked as unknown.", nodes[i].name);
                unknown_ids.push_back(static_cast<int>(i));
            }
        }
    }
    time_print("%zu group(s) are generated.", chr_groups.size());
    //Generate the group ids.
    groups.reserve(chr_groups.size());
    for(HANA_GROUP *group: chr_groups)
    {
        //Add node ids to group.
        std::vector<int> node_ids(group->ids.begin(), group->ids.end());
        groups.push_back(node_ids);
        time_print("\tGroup %2zu: %zu contig(s).", groups.size(), node_ids.size());
        //Delete the group.
        delete group;
    }
    //Free the memory.
    delete[] node_marks;
    delete[] marks;
    delete[] node_infos;
}
