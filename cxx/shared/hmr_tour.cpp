#include <vector>

#include "hmr_cgd_type.h"
#include "hmr_ui.h"

#include "hmr_tour.h"

void hmr_tour_dump(const char *filepath, CONTIG_NODE *nodes, const POLAR_INFO &polar_infos)
{
#ifdef _MSC_VER
    FILE *tour_file = NULL;
    fopen_s(&tour_file, filepath, "w");
#else
    FILE *tour_file = fopen(filepath, "w");
#endif
    if(tour_file == NULL)
    {
        time_error(-1, "Failed to write tour node file %s", filepath);
    }
    //Loop and write the node and plus minus.
    for(size_t i=0; i<polar_infos.size(); ++i)
    {
        const auto &node_polar = polar_infos[i];
        if(i != 0)
        {
            fprintf(tour_file, " ");
        }
        //Write the node name with plus and minus.
        fprintf(tour_file, "%s%c", nodes[node_polar.node_id].name, node_polar.plus ? '+' : '-');
    }
    fclose(tour_file);
}
