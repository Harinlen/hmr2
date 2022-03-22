#ifndef MISMATCH_CORRECT_H
#define MISMATCH_CORRECT_H

#include <vector>
#include "hmr_parallel.h"

#include "mapping_correct_type.h"

typedef std::vector<POS_PAIR> RANGE_LIST;

HMR_PFOR_FUNC(mismatch_calc, HIC_DB* wide_db, HIC_DB* narrow_db, double percent, double sens, int32_t dep, int32_t wide, int32_t narrow, RANGE_LIST* mismatches);

typedef struct MISMATCH_CORRECTING
{
    RANGE_LIST* mismatches;
    FILE* fp;
} MISMATCH_CORRECTING;

void mismatch_correct_open(const char* filepath, MISMATCH_CORRECTING* correct_file);
void mismatch_corrected(int32_t index, char* seq_name, size_t seq_name_size, char* seq, size_t seq_size, void* user);

#endif // MISMATCH_CORRECT_H