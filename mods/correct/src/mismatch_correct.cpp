#include <cstring>
#include <queue>
#include <unordered_map>

#include "hmr_bin_file.h"
#include "hmr_ui.h"

#include "mismatch_correct.h"

typedef std::unordered_map<int32_t, double> SCORE_DB;

inline double round5(double value)
{
    return static_cast<double>(static_cast<int64_t>(value * 100000.0)) / 100000.0;
}

double sat_level(const HIC_DB& hic_db, double percent)
{
    //Find the non-self related ranges.
    std::vector<int32_t> hic_impact_queue;
    hic_impact_queue.reserve(hic_db.size());
    std::priority_queue<int32_t, std::vector<int32_t>, std::greater<int32_t> > hic_impact(std::greater<int32_t>(), std::move(hic_impact_queue));
    for (const auto &range_info : hic_db)
    {
        POS_PAIR pair {};
        pair.data = range_info.first;
        if (pair.pos.a != pair.pos.b)
        {
            hic_impact.push(range_info.second);
        }
    }
    //Check the size of the relation map.
    if (hic_impact.empty())
    {
        return -1.0;
    }
    if (hic_impact.size() == 1)
    {
        return static_cast<double>(hic_impact.top());
    }
    //Calculate the expected position.
    double pos = static_cast<double>(hic_impact.size() + 1) * percent;
    if (pos < 1.0)
    {
        return static_cast<double>(hic_impact.top());
    }
    std::vector<int32_t> hic_impacts;
    hic_impacts.reserve(hic_impact.size());
    while (!hic_impact.empty()) { hic_impacts.push_back(hic_impact.top()); hic_impact.pop(); }
    if (pos >= static_cast<double>(hic_impacts.size()))
    {
        return static_cast<double>(hic_impacts.back());
    }
    int32_t pos_d = static_cast<int32_t>(pos);
    double pos_dpv = static_cast<double>(hic_impacts[pos_d - 1]), pos_dv = static_cast<double>(hic_impacts[pos_d]);
    return pos_dpv + (pos - static_cast<double>(pos_d)) * (pos_dv - pos_dpv);
}

void precompute_dep_score(const HIC_DB& hic_db, int32_t bin_size, int32_t dep_size, int32_t sat_level, SCORE_DB &score_db, int32_t &min_score, int32_t &max_score)
{
    //Initial the min and max.
    min_score = 0;
    max_score = 0;
    //Loop for all the ranges in hic_db.
    for (const auto& range_info : hic_db)
    {
        POS_PAIR pair{};
        pair.data = range_info.first;
        int32_t s = pair.pos.a, e = pair.pos.b;
        //Seems like a marking scheme?
        if (e - s <= dep_size)
        {
            double se_count = static_cast<double>(range_info.second);
            if (se_count >= sat_level)
            {
                se_count = sat_level;
            }
            for (int32_t i = s + bin_size; i < e; i += bin_size)
            {
                auto score_finder = score_db.find(i);
                if (score_finder == score_db.end())
                {
                    score_db.insert(std::make_pair(i, se_count));
                    min_score = (min_score == 0) ? i : hMin(i, min_score);
                    max_score = hMax(i, max_score);
                }
                else
                {
                    (*score_finder).second += se_count;
                }
            }
        }
    }
}

inline SCORE_DB get_sub_score_db(const SCORE_DB& score_db, int32_t min_pos, int32_t max_pos, int32_t bin_size, int32_t dep_size)
{
    SCORE_DB sub_db;
    for (int32_t i = min_pos + dep_size - (bin_size * 2), max_i = (max_pos - dep_size + 3 * bin_size); i < max_i; i += bin_size)
    {
        auto score_finder = score_db.find(i);
        if (score_finder == score_db.end())
        {
            sub_db.insert(std::make_pair(i, 0));
        }
        else
        {
            sub_db.insert(std::make_pair(i, (*score_finder).second));
        }
    }
    return sub_db;
}

HMR_PFOR_FUNC(mismatch_calc, HIC_DB* wide_db, HIC_DB* narrow_db, double percent, double sens, int32_t dep, int32_t wide, int32_t narrow, RANGE_LIST* mismatches)
{
    //When the contig database is empty, just go back.
    if (wide_db[idx].empty())
    {
        return;
    }
    //Calculate the sat level.
    double sat_wide = round5(sat_level(wide_db[idx], percent));
    double dep_f = static_cast<double>(dep), wide_f = static_cast<double>(wide), narrow_f = static_cast<double>(narrow);
    //If sat is -1, we don't have to calculate the mismatch array.
    RANGE_LIST wide_mismatch, narrow_mismatch;
    if (sat_wide != -1)
    {
        //Get the dep score.
        SCORE_DB wide_score; 
        int32_t wide_min_pos, wide_max_pos;
        double threshold = sens * sat_wide * 0.5 * dep_f / wide_f * (dep_f / wide_f - 1.0);
        precompute_dep_score(wide_db[idx], wide, dep, sat_wide, wide_score, wide_min_pos, wide_max_pos);
        if (!wide_score.empty())
        {
            wide_score = get_sub_score_db(wide_score, wide_min_pos, wide_max_pos, wide, dep);
            if (!wide_score.empty())
            {
                //Construct the wide mismatch.
                std::vector<int32_t> wide_poses_queue;
                wide_poses_queue.reserve(wide_score.size());
                std::priority_queue<int32_t, std::vector<int32_t>, std::greater<int32_t> > wide_poses(std::greater<int32_t>(), std::move(wide_poses_queue));
                for (const auto& wide_score_map : wide_score)
                {
                    wide_poses.push(wide_score_map.first);
                }
                //Pop the data out.
                bool is_a = true;
                POS_PAIR pair {};
                int32_t wide_pos;
                while (!wide_poses.empty())
                {
                    wide_pos = wide_poses.top();
                    //Check the top mark.
                    if ((*wide_score.find(wide_pos)).second < threshold)
                    {
                        if(is_a)
                        {
                            //Only available for position a.
                            pair.pos.a = wide_pos;
                            is_a = false;
                        }
                    }
                    else
                    {
                        if (!is_a)
                        {
                            //Only available for position b.
                            pair.pos.b = wide_pos;
                            wide_mismatch.push_back(pair);
                            is_a = true;
                        }
                    }
                    //Pop the top.
                    wide_poses.pop();
                }
                //Check the position.
                if (!is_a)
                {
                    pair.pos.b = wide_pos + wide;
                    wide_mismatch.push_back(pair);
                }
            }
        }
    }
    //When we have the wide mismatch, we calcualte the narrow mismatch
    if (!wide_mismatch.empty())
    {
        double sat_narrow = sat_level(narrow_db[idx], percent);
        //Get the dep score.
        SCORE_DB narrow_score;
        int32_t narrow_min_pos, narrow_max_pos;
        precompute_dep_score(narrow_db[idx], narrow, wide, sat_narrow, narrow_score, narrow_min_pos, narrow_max_pos);
        if (!narrow_score.empty())
        {
            narrow_score = std::move(get_sub_score_db(narrow_score, narrow_min_pos, narrow_max_pos, narrow, wide));
        }
        //If no narrow score, then use the wide mismatch.
        if (narrow_score.empty())
        {
            narrow_mismatch = std::move(wide_mismatch);
        }
        else
        {
            //Merge the narrow score into wide mismatch, get the narrow mismatch.
            int32_t idx_wide = 0, wide_length = static_cast<int32_t>(wide_mismatch.size());
            double min_val = 0.0;
            RANGE_LIST tmp_list;
            std::vector<int32_t> narrow_poses_queue;
            narrow_poses_queue.reserve(narrow_score.size());
            std::priority_queue<int32_t, std::vector<int32_t>, std::greater<int32_t> > narrow_poses(std::greater<int32_t>(), std::move(narrow_poses_queue));
            for (const auto& narrow_score_map : narrow_score)
            {
                narrow_poses.push(narrow_score_map.first);
            }
            while (!narrow_poses.empty())
            {
                int32_t pos = narrow_poses.top();
                narrow_poses.pop();
                //Check the index wide reaches the limit.
                if (idx_wide >= wide_length)
                {
                    break;
                }
                double pos_narrow_score = (*narrow_score.find(pos)).second;
                const auto& wide_mismatch_idx = wide_mismatch[idx_wide];
                if (pos <= wide_mismatch_idx.pos.a)
                {
                    min_val = pos_narrow_score;
                }
                else
                {
                    if (pos_narrow_score < min_val)
                    {
                        min_val = pos_narrow_score;
                    }
                }
                if (pos + narrow <= wide_mismatch_idx.pos.a)
                {
                    continue;
                }
                if (pos >= wide_mismatch_idx.pos.b)
                {
                    for (int32_t i = wide_mismatch_idx.pos.a; i < wide_mismatch_idx.pos.b; i += narrow)
                    {
                        auto narrow_i_finder = narrow_score.find(i);
                        if (narrow_i_finder != narrow_score.end() && (*narrow_i_finder).second == min_val)
                        {
                            tmp_list.push_back(POS_PAIR{ {i, i + narrow} });
                        }
                    }
                    ++idx_wide;
                }
            }
            if (idx_wide < wide_length)
            {
                const auto& wide_mismatch_idx = wide_mismatch[idx_wide];
                for (int32_t i = wide_mismatch_idx.pos.a; i < wide_mismatch_idx.pos.b; i += narrow)
                {
                    auto narrow_i_finder = narrow_score.find(i);
                    if (narrow_i_finder != narrow_score.end() && (*narrow_i_finder).second == min_val)
                    {
                        tmp_list.push_back(POS_PAIR{ {i, i + narrow} });
                    }
                }
            }
            if (tmp_list.empty())
            {
                narrow_mismatch = std::move(wide_mismatch);
            }
            else
            {
                //Construct the narrow mismatch.
                int32_t last_e = 0;
                POS_PAIR narrow_pair {};
                for (const POS_PAIR& pair : tmp_list)
                {
                    const int32_t s = pair.pos.a, e = pair.pos.b;
                    if (last_e == 0)
                    {
                        narrow_pair.pos.a = s;
                    }
                    else
                    {
                        if (s != last_e)
                        {
                            narrow_pair.pos.b = last_e;
                            narrow_mismatch.push_back(narrow_pair);
                            narrow_pair.pos.a = s;
                        }
                    }
                    last_e = e;
                }
                narrow_pair.pos.b = last_e;
                narrow_mismatch.push_back(narrow_pair);
            }
        }
    }
    //Set the mismatch result.
    mismatches[idx] = std::move(narrow_mismatch);
}

void mismatch_correct_open(const char* filepath, MISMATCH_CORRECTING* correct_file)
{
    //Try to open the file for written, use binary type to open it.
    if (!bin_open(filepath, &correct_file->fp, "wb"))
    {
        time_error(-1, "Failed to open corrected FASTA file: %s", filepath);
    }
}

void mismatch_corrected(int32_t index, char* seq_name, size_t seq_name_size, char* seq, size_t seq_size, void* user)
{
    MISMATCH_CORRECTING* correct_file = reinterpret_cast<MISMATCH_CORRECTING*>(user);
    //Check whether the contig is splited.
    FILE* fp = correct_file->fp;
    const RANGE_LIST& idx_range = correct_file->mismatches[index];
    if (idx_range.empty())
    {
        //Just write the name and sequence.
        fwrite(">", 1, 1, fp);
        fwrite(seq_name, 1, seq_name_size, fp);
        fwrite("\n", 1, 1, fp);
        fwrite(seq, 1, seq_size, fp);
        fwrite("\n", 1, 1, fp);
    }
    else
    {
        //Split the sequence based on mismatch information.
        char name_buf[1024];
        int32_t base = 0;
        for (const auto& edge : idx_range)
        {
            int32_t s = edge.pos.a - 1, e = edge.pos.b - 1;
            //Name: >name_base_s
            fwrite(">", 1, 1, fp);
            fwrite(seq_name, 1, seq_name_size, fp);
#ifdef _MSC_VER
            sprintf_s(name_buf, 1023, "_%d_%d\n", base + 1, s);
#else
            sprintf(name_buf, "_%d_%d\n", base + 1, s);
#endif
            fwrite(name_buf, 1, strlen(name_buf), fp);
            fwrite(seq + base, 1, s - base, fp);
            fwrite("\n", 1, 1, fp);
            //Name: >name_s_e
            fwrite(">", 1, 1, fp);
            fwrite(seq_name, 1, seq_name_size, fp);
#ifdef _MSC_VER
            sprintf_s(name_buf, 1023, "_%d_%d\n", s + 1, e);
#else
            sprintf(name_buf, "_%d_%d\n", s + 1, e);
#endif
            fwrite(name_buf, 1, strlen(name_buf), fp);
            fwrite(seq + s, 1, e - s, fp);
            fwrite("\n", 1, 1, fp);
            //Update the base.
            base = e;
        }
        //Check whether the e reaches the end.
        if (base < seq_size)
        {
            //Name: >name_base_seqsize
            fwrite(">", 1, 1, fp);
            fwrite(seq_name, 1, seq_name_size, fp);
#ifdef _MSC_VER
            sprintf_s(name_buf, 1023, "_%d_%zu\n", base, seq_size);
#else
            sprintf(name_buf, "_%d_%zu\n", base, seq_size);
#endif
            fwrite(name_buf, 1, strlen(name_buf), fp);
            fwrite(seq + base, 1, seq_size - base, fp);
            fwrite("\n", 1, 1, fp);
        }
    }
    //Recover the memory.
    free(seq);
    free(seq_name);
}
