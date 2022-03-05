#ifndef HMR_POLAR_TYPE_H
#define HMR_POLAR_TYPE_H

#include <vector>
#include <list>

typedef struct NODE_POLAR NODE_POLAR;

typedef struct READ_PAIR
{
    int32_t pos;
    NODE_POLAR *pair;
    int32_t pair_pos;
} READ_PAIR;

typedef struct NODE_POLAR
{
    int node_id, seq_pos, length;
    bool plus;
    std::list<READ_PAIR> reads;
} NODE_POLAR;

typedef std::vector<NODE_POLAR> POLAR_INFO;

#endif // HMR_POLAR_TYPE_H
