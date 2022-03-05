#include <cstdlib>
#include <list>

#include "optimize_direction.h"

void optimize_direction(POLAR_INFO &polar_infos)
{
    //Use greedy to determine the direction.
    for(size_t i=0; i<polar_infos.size(); ++i)
    {
        int my_seq = static_cast<int>(i);
        int my_length = polar_infos[i].length;
        //Find whether the position is on left or right.
        std::list<READ_PAIR> on_left, on_right;
        for(auto read: polar_infos[i].reads)
        {
            int pair_seq_pos = read.pair->seq_pos;
            if(pair_seq_pos < my_seq)
            {
                on_left.push_back(read);
            }
            else
            {
                on_right.push_back(read);
            }
        }
        size_t plus_distance = 0, minus_distance = 0;
        for(const auto &pair: on_left)
        {
            plus_distance += pair.pos;
            minus_distance += my_length - pair.pos;
        }
        for(const auto &pair: on_right)
        {
            plus_distance += my_length - pair.pos;
            minus_distance += pair.pos;
        }
        //Decide the direction.
        polar_infos[i].plus = plus_distance > minus_distance;
    }
}
