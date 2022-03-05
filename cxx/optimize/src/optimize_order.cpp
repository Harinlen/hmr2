//Harerun ID: 826917359
#include <algorithm>
#include <set>
#include <unordered_set>

#include "hmr_cgd_type.h"
#include "hmr_global.h"

#include "optimize_order.h"

#define INVALID     (-1)
#define NO_WEIGHT   (-2)

typedef struct EDGE
{
    int32_t a;
    int32_t b;
} EDGE;

typedef union EDGE_INFO
{
    EDGE ids;
    uint64_t data;
} EDGE_INFO;

typedef struct TRUST_EDGE
{
    EDGE_INFO edge;
    double weight;
} TRUST_EDGE;

bool contig_edge_compare(const CONTIG_EDGE &lhs, const CONTIG_EDGE &rhs)
{
    return lhs.weight > rhs.weight;
}

template <typename T>
inline bool in_set(const std::unordered_set<T> &set, const T &element)
{
    return set.find(element) != set.end();
}

void fit_linear(const std::vector<double> &x, const std::vector<double> &y,
                double &k, double &b)
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

inline double linear_loss(double x, double y, double k, double b)
{
    double A = k, B = -1, C = b;
    return hAbs(A * x + B * y + C) / sqrt(A * A + B * B);
}

inline TRUST_EDGE to_trust_edge(int32_t src_id, const CONTIG_EDGE &edge)
{
    //Smaller id always comes first.
    int32_t small_id = src_id, big_id = edge.id;
    if(small_id > big_id)
    {
        small_id = edge.id;
        big_id = src_id;
    }
    //Construct the trust edge.
    return TRUST_EDGE { EDGE_INFO { EDGE {small_id, big_id} }, edge.weight};
}

inline int32_t node_order(int32_t dst_id, const ORDER_EDGE_MAP &src_map)
{
    auto id_finder = src_map.find(dst_id);
    return id_finder == src_map.end() ? NO_WEIGHT : id_finder->second.order;
}

inline double node_weight(int32_t dst_id, const ORDER_EDGE_MAP &src_map)
{
    auto id_finder = src_map.find(dst_id);
    return id_finder == src_map.end() ? 0.0 : id_finder->second.weight;
}

inline double edge_penalty(int32_t src_id, int32_t dst_id, ORDER_NODE_INFO *node_infos)
{
    return hMin(node_infos[src_id].penalty, node_infos[dst_id].penalty);
}

typedef struct SEQ_WEIGHT_RANK
{
    int32_t pos;
    double weight;
} SEQ_WEIGHT_RANK;

bool weight_rank_compare(const SEQ_WEIGHT_RANK &lhs, const SEQ_WEIGHT_RANK &rhs)
{
    return lhs.weight > rhs.weight;
}

void merge_node(int32_t node_id, SEQUENCE *seq, ORDER_NODE_INFO *node_infos)
{
    ORDER_NODE_INFO &node_info = node_infos[node_id];
    const auto &node_edge_map = node_info.edge_map;
    //Find out the most related position of node in the sequence.
    std::vector<SEQ_WEIGHT_RANK> weight_seq;
    weight_seq.reserve(seq->size());
    int32_t pos = 0;
    for(int32_t seq_id: *seq)
    {
        weight_seq.push_back(SEQ_WEIGHT_RANK {pos, node_weight(seq_id, node_edge_map) * edge_penalty(node_id, seq_id, node_infos)});
        ++pos;
    }
    std::stable_sort(weight_seq.begin(), weight_seq.end(), weight_rank_compare);
    //Depends on the position, update the sequence and node info.
    node_info.parent = seq;
    int32_t position = weight_seq[0].pos;
    //Check whether we can insert into this position.
    int32_t position_id = (*seq)[position];
    ORDER_NODE_INFO &position_info = node_infos[position_id];
    //Check the position.
    if(position == 0)
    {
        // Update the node infos.
        node_info.prev_order = INVALID;
        node_info.next_order = node_order(position_id, node_edge_map);
        position_info.prev_order = node_order(node_id, position_info.edge_map);
        // Go left-most.
        seq->insert(seq->begin(), node_id);
        return;
    }
    if(position == static_cast<int32_t>(weight_seq.size() - 1))
    {
        // Update the node infos.
        node_info.prev_order = node_order(position_id, node_edge_map);
        node_info.next_order = INVALID;
        position_info.next_order = node_order(node_id, position_info.edge_map);
        // Go right-most.
        seq->push_back(node_id);
        return;
    }
    //Check should we insert the node here.
    // Compare the current order with left and right.
    int32_t left_id = (*seq)[position - 1], right_id = (*seq)[position + 1];
    //Decide which side should be pushed.
    if(node_order(left_id, node_edge_map) < node_order(right_id, node_edge_map))
    {
        //Update the node info itself.
        ORDER_NODE_INFO &left_info = node_infos[left_id];
        node_info.prev_order = node_order(left_id, node_edge_map);
        node_info.next_order = node_order(position, node_edge_map);
        //Update the left and position info.
        left_info.next_order = node_order(node_id, left_info.edge_map);
        position_info.prev_order = node_order(node_id, position_info.edge_map);
        //Insert to left side and position.
        seq->insert(seq->begin() + position, node_id);
    }
    else
    {
        //Update the node info itself.
        ORDER_NODE_INFO &right_info = node_infos[right_id];
        node_info.prev_order = node_order(position, node_edge_map);
        node_info.next_order = node_order(right_id, node_edge_map);
        //Update the position and right info.
        position_info.next_order = node_order(node_id, position_info.edge_map);
        right_info.prev_order = node_order(node_id, right_info.edge_map);
        //Insert to position and right side.
        seq->insert(seq->begin() + position + 1, node_id);
    }
}

#define HEAD_TO_HEAD    (0)
#define HEAD_TO_TAIL    (1)
#define TAIL_TO_HEAD    (2)
#define TAIL_TO_TAIL    (3)

SEQUENCE *merge_sequence(SEQUENCE *a, SEQUENCE *b, int mode, ORDER_NODE_INFO *node_infos)
{
    switch(mode)
    {
    //Head to head:
    //  ida ---
    //  idb ---
    case HEAD_TO_HEAD:
    {
        //Swap all the seq a node prev and next.
        for(int32_t id: *a)
        {
            std::swap(node_infos[id].prev_order, node_infos[id].next_order);
        }
        //Update the parent of seq b.
        for(int32_t id: *b)
        {
            node_infos[id].parent = a;
        }
        //Reverse seq_a, append seq_b to seq_a.
        std::reverse(a->begin(), a->end());
        a->insert(a->end(), b->begin(), b->end());
        delete b;
        return a;
    }
    //Tail to tail:
    //  --- ida
    //  --- idb
    case TAIL_TO_TAIL:
    {
        //Swap all the seq b node prev and next, and update the parent.
        for(int32_t id: *b)
        {
            auto &node_info = node_infos[id];
            std::swap(node_info.prev_order, node_info.next_order);
            node_info.parent = a;
        }
        //Reverse seq_b, append seq_b to seq_a.
        std::reverse(b->begin(), b->end());
        a->insert(a->end(), b->begin(), b->end());
        delete b;
        return a;
    }
    //Tail to Head.
    //  --- ida
    //  idb ---
    case TAIL_TO_HEAD:
    {
        //Update all the seq b parent.
        for(int32_t id: *b)
        {
            node_infos[id].parent = a;
        }
        //Append seq_b to seq_a.
        a->insert(a->end(), b->begin(), b->end());
        delete b;
        return a;
    }
    //Head to tail.
    //  ida ---
    //  --- idb
    case HEAD_TO_TAIL:
    {
        //Update all the seq a parent.
        for(int32_t id: *a)
        {
            node_infos[id].parent = b;
        }
        //Append seq_a to seq_b.
        b->insert(b->end(), a->begin(), a->end());
        delete a;
        return b;
    }
    default:
        return NULL;
    }
}

typedef struct SEQ_MERGE_INFO
{
    SEQUENCE *a;
    SEQUENCE *b;
    int32_t mode;
    double marks;
} SEQ_MERGE_INFO;

bool merge_info_compare(const SEQ_MERGE_INFO &lhs, const SEQ_MERGE_INFO &rhs)
{
    return lhs.marks > rhs.marks;
}

#define MAX_EDGES   (4)

SEQ_MERGE_INFO sequence_relation(SEQUENCE *a, SEQUENCE *b, ORDER_NODE_INFO *node_infos)
{
    typedef struct WEIGHT_EDGE
    {
        size_t s_pos;
        size_t b_pos;
        double weight;
    } WEIGHT_EDGE;
    //Find the best 4 edges of both sequences.
    WEIGHT_EDGE weight[MAX_EDGES];
    for(int i=0; i<MAX_EDGES; ++i)
    {
        weight[i].weight = 0.0;
    }
    SEQUENCE *small = a, *big = b;
    if(b->size() < a->size())
    {
        small = b;
        big = a;
    }
    for(size_t i=0; i<(*small).size(); ++i)
    {
        int32_t s_id = (*small)[i];
        //Extract the node weight map.
        const ORDER_EDGE_MAP &s_map = node_infos[s_id].edge_map;
        for(size_t j=0; j<(*big).size(); ++j)
        {
            int32_t b_id = (*big)[j];
            double b_weight = node_weight(b_id, s_map) * edge_penalty(s_id, b_id, node_infos);
            if(b_weight > weight[MAX_EDGES - 1].weight)
            {
                //Time to insert the data.
                for(int32_t k=0; k<MAX_EDGES; ++k)
                {
                    if(b_weight > weight[k].weight)
                    {
                        for(int32_t p=MAX_EDGES-1; p>k; --p)
                        {
                            weight[p] = weight[p-1];
                        }
                        weight[k] = WEIGHT_EDGE {i, j, b_weight};
                        break;
                    }
                }
            }
        }
    }
    //Calculate the sum of the total weights and decide the mode.
    double weight_sum = 0.0;
    size_t s_head_pos = 0, b_head_pos = 0, s_tail_pos = 0, b_tail_pos = 0,
            s_len = small->size() - 1, b_len = big->size() - 1;
    for(int i=0; i<MAX_EDGES; ++i)
    {
        weight_sum += weight[i].weight;
        //Calculate the pos.
        s_head_pos += weight[i].s_pos;
        b_head_pos += weight[i].b_pos;
        s_tail_pos += s_len - weight[i].s_pos;
        b_tail_pos += b_len - weight[i].b_pos;
    }
    int mode = 0;
    if(s_head_pos < s_tail_pos)
    {
        if(b_head_pos < b_tail_pos)
        {
            mode = HEAD_TO_HEAD;
        }
        else
        {
            mode = HEAD_TO_TAIL;
        }
    }
    else
    {
        if(b_head_pos < b_tail_pos)
        {
            mode = TAIL_TO_HEAD;
        }
        else
        {
            mode = TAIL_TO_TAIL;
        }
    }
    return SEQ_MERGE_INFO {small, big, mode, weight_sum};
}

SEQUENCE optimize_order(const std::vector<int32_t> &node_ids, CONTIG_NODE *nodes,
                        ORDER_NODE_INFO *node_infos)
{
    std::unordered_set<int> node_set(node_ids.begin(), node_ids.end());
    std::vector<TRUST_EDGE> trust_edges;
    {
        std::unordered_set<uint64_t> existed_edges;
        //Prepare the node information.
        for(const int32_t id: node_ids)
        {
            ORDER_NODE_INFO &node_info = node_infos[id];
            node_info.parent = NULL;
            node_info.prev_order = INVALID;
            node_info.next_order = INVALID;
            //Sort the edge of the nodes.
            auto &edges = nodes[id].links;
            std::sort(edges.begin(), edges.end(), contig_edge_compare);
            //Calculate the penalty of the node.
            size_t hit_edges = 0,
                    penalty_range = hMin(node_ids.size(), edges.size());
            for(size_t i=0; i<penalty_range; ++i)
            {
                if(in_set(node_set, edges[i].id))
                {
                    ++hit_edges;
                }
            }
            node_info.penalty = static_cast<double>(hit_edges) / static_cast<double>(penalty_range);
            //Find the trust edges and build the edge map.
            std::vector<CONTIG_EDGE> node_edges;
            node_edges.reserve(edges.size());
            int32_t order = 1;
            for(const auto &edge: edges)
            {
                if(in_set(node_set, edge.id))
                {
                    //Record the edge.
                    node_edges.push_back(edge);
                    //Create the edge map.
                    node_info.edge_map.insert(
                                std::make_pair(edge.id, ORDER_EDGE {edge.weight, order}));
                    ++order;
                }
            }
            //Calculate the 1st derivative.
            std::vector<double> deri_1st;
            deri_1st.reserve(node_edges.size());
            deri_1st.push_back(0.0);
            for(size_t i=1; i<node_edges.size(); ++i)
            {
                deri_1st.push_back(hAbs(node_edges[i].weight - node_edges[i-1].weight));
            }
            size_t trust_pos = std::distance(
                        deri_1st.begin(),
                        std::max_element(deri_1st.begin(), deri_1st.end()));
            node_edges.resize(trust_pos);
            //Save as the trust edges.
            trust_edges.reserve(trust_edges.size() + node_edges.size());
            for(const auto &edge: node_edges)
            {
                TRUST_EDGE trust_edge = to_trust_edge(id, edge);

                //Check edge existance.
                if(!in_set(existed_edges, trust_edge.edge.data))
                {
                    trust_edges.push_back(trust_edge);
                    //Record the edge in set.
                    existed_edges.insert(trust_edge.edge.data);
                }
            }
        }
    }
    //Sort the trust edges.
    std::stable_sort(trust_edges.begin(), trust_edges.end(),
                     [](const TRUST_EDGE &lhs, const TRUST_EDGE &rhs){
        return lhs.weight > rhs.weight;
    });
    //First half of the turst edges are good edges, train a linear equation and
    //calculation the maximum loss.
    {
        size_t best_range = trust_edges.size() >> 1;
        std::vector<double> edge_x, edge_y;
        edge_x.reserve(best_range);
        edge_y.reserve(best_range);
        for(size_t i=0; i<best_range; ++i)
        {
            edge_x.push_back(static_cast<double>(i));
            edge_y.push_back(trust_edges[i].weight);
        }
        double trust_k, trust_b, max_loss = -1.0;
        fit_linear(edge_x, edge_y, trust_k, trust_b);
        for(size_t i=0; i<best_range; ++i)
        {
            max_loss = hMax(linear_loss(static_cast<double>(i), trust_edges[i].weight, trust_k, trust_b),
                            max_loss);
        }
        //Find out the trust edges exceed the maximum loss.
        std::list<size_t> invalid_ids;
        for(size_t i=best_range; i<trust_edges.size(); ++i)
        {
            double loss = linear_loss(static_cast<double>(i), trust_edges[i].weight, trust_k, trust_b);
            if(loss > max_loss)
            {
                invalid_ids.push_front(i);
            }
        }
        //Remove the invalid edges.
        if(invalid_ids.size() > 0)
        {
            for(size_t id: invalid_ids)
            {
                trust_edges.erase(trust_edges.begin() + id);
            }
        }
    }
    //Apply the combination of the optimize trusted edge.
    for(const auto &trust_edge: trust_edges)
    {
        int32_t a_id = trust_edge.edge.ids.a, b_id = trust_edge.edge.ids.b;
        ORDER_NODE_INFO &a_info = node_infos[a_id],
                &b_info = node_infos[b_id];
        //Check the group information.
        if(a_info.parent == NULL)
        {
            if(b_info.parent == NULL)
            {
                //Create a new group which contains these two nodes.
                SEQUENCE *seq = new SEQUENCE();
                seq->reserve(2);
                seq->push_back(a_id);
                seq->push_back(b_id);
                a_info.next_order = node_order(b_id, a_info.edge_map);
                b_info.prev_order = node_order(a_id, b_info.edge_map);
                a_info.parent = seq;
                b_info.parent = seq;
            }
            else
            {
                //Insert the node a to sequence of node b.
                merge_node(a_id, b_info.parent, node_infos);
            }
        }
        else
        {
            if(b_info.parent == NULL)
            {
                //Insert the node b to sequence of node a.
                merge_node(b_id, a_info.parent, node_infos);
            }
            else
            {
                //Check whether both node are already inside a group.
                if(a_info.parent != b_info.parent)
                {
                    //Check the position of these two node.
                    SEQUENCE *seq_a = a_info.parent, *seq_b = b_info.parent;
                    int32_t a_head = seq_a->front(), b_head = seq_b->front(),
                            a_tail = seq_a->back(), b_tail = seq_b->back();
                    if(a_head == a_id)
                    {
                        if(b_head == b_id)
                        {
                            merge_sequence(seq_a, seq_b, HEAD_TO_HEAD, node_infos);
                        }
                        else if(b_tail == b_id)
                        {
                            merge_sequence(seq_a, seq_b, HEAD_TO_TAIL, node_infos);
                        }
                    }
                    else if(a_tail == a_id)
                    {
                        if(b_head == b_id)
                        {
                            merge_sequence(seq_a, seq_b, TAIL_TO_HEAD, node_infos);
                        }
                        else if(b_tail == b_id)
                        {
                            merge_sequence(seq_a, seq_b, TAIL_TO_TAIL, node_infos);
                        }
                    }
                    //For all the other cases, just ignore the edge.
                }
            }
        }
    }
    //Pick out all the sequence.
    std::unordered_set<SEQUENCE *> sequences;
    std::vector<int32_t> individuals;
    individuals.reserve(node_ids.size());
    for(const int32_t id: node_ids)
    {
        SEQUENCE *seq = node_infos[id].parent;
        if(seq == NULL)
        {
            individuals.push_back(id);
        }
        else
        {
            sequences.insert(seq);
        }
    }
    //Keep merging all the sequences into a single sequences.
    while(sequences.size() > 1)
    {
        std::vector<SEQUENCE *> seq_vec(sequences.begin(), sequences.end());
        std::vector<SEQ_MERGE_INFO> merge_info;
        merge_info.reserve(seq_vec.size() * seq_vec.size());
        //Calculate the relationship.
        for(size_t i=0; i<seq_vec.size()-1; ++i)
        {
            for(size_t j=i+1; j<seq_vec.size(); ++j)
            {
                merge_info.push_back(sequence_relation(seq_vec[i], seq_vec[j], node_infos));
            }
        }
        //Sort the merge information.
        std::sort(merge_info.begin(), merge_info.end(), merge_info_compare);
        //Pick out the top level merge option.
        const SEQ_MERGE_INFO &action = merge_info[0];
        //We merge those two groups.
        SEQUENCE *seq_remain = merge_sequence(action.a, action.b, action.mode, node_infos);
        if(seq_remain == action.a)
        {
            //Remove seq b.
            sequences.erase(action.b);
        }
        else
        {
            sequences.erase(action.a);
        }
    }
    //Now we pick out the sequences.
    SEQUENCE *group_seq = *sequences.begin();
    //Keep merging the individual nodes into the sequence.
    for(const int32_t id: individuals)
    {
        merge_node(id, group_seq, node_infos);
    }
    SEQUENCE result = SEQUENCE(group_seq->begin(), group_seq->end());
    delete group_seq;
    return result;
}
