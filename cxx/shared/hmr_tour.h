#ifndef HMR_TOUR_H
#define HMR_TOUR_H

#include "hmr_polar_type.h"

typedef struct CONTIG_NODE CONTIG_NODE;

void hmr_tour_dump(const char *filepath, CONTIG_NODE *nodes, const POLAR_INFO &polar_infos);

#endif // HMR_TOUR_H
