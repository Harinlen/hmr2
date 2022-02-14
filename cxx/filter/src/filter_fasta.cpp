#include "hmr_args_type.h"
#include "filter_fasta_type.h"

#include "filter_fasta.h"

extern HMR_ARGS opts;

static std::mutex search_mutex;

void filter_fasta_push_work(char *seq_name, size_t seq_name_length, char *seq, size_t seq_size, void *user)
{
    //Recast the user into filter user.
    FASTA_FILTER_USER *filter_user = reinterpret_cast<FASTA_FILTER_USER *>(user);
    //Now push the task to workers.
    filter_user->workers->push_task(FASTA_FILTER_WORK_ARGS {seq_name, seq_name_length, seq, seq_size, &(filter_user->enzyme_poses)});
}

void filter_fasta_search_enzyme(const FASTA_FILTER_WORK_ARGS &args)
{
    //Use strstr to search all the enzyme poses.
    std::list<size_t> poses;
    char *inseq_pos = strstr(args.seq, opts.nuc_seq);
    while(inseq_pos != NULL)
    {
        //Save the record.
        poses.push_back(inseq_pos - args.seq);
        //Keep searching.
        inseq_pos = strstr(inseq_pos + opts.nuc_seq_length, opts.nuc_seq);
    }
    //Check whether the poses are emtpy.
    if(!poses.empty())
    {
        //Construct the range list.
        ENZYME_RANGE *pos_ranges = static_cast<ENZYME_RANGE *>(malloc(sizeof(ENZYME_RANGE) * poses.size()));
        size_t counter = 0;
        for(auto i : poses)
        {
            if(args.seq_size < 1000) { pos_ranges[counter] = ENZYME_RANGE {0, args.seq_size}; }
            //Now it must be larger than 1000, check the position.
            else if(i < 500) { pos_ranges[counter] = ENZYME_RANGE {0, i + 500}; }
            else if(i + opts.nuc_seq_length > args.seq_size) { pos_ranges[counter] = ENZYME_RANGE {i - 500, args.seq_size}; }
            else { pos_ranges[counter] = ENZYME_RANGE {i - 500, i + 500}; }
            ++counter;
        }
        //Fetch the modify lock, update the result.
        std::unique_lock<std::mutex> search_lock(search_mutex);
        auto enzyme_poses = args.enzyme_poses;
        size_t ref_id = enzyme_poses->seq_total;
        //Append the name.
        enzyme_poses->seq_names = static_cast<char **>(realloc(enzyme_poses->seq_names, sizeof(char *) * (1 + ref_id)));
        enzyme_poses->seq_names[ref_id] = static_cast<char *>(malloc(args.seq_name_length + 1));
        strncpy(enzyme_poses->seq_names[ref_id], args.seq_name, args.seq_name_length + 1);
        //Record the sequence length.
        enzyme_poses->seq_length = static_cast<size_t *>(realloc(enzyme_poses->seq_length, sizeof(size_t) * (1 + ref_id)));
        enzyme_poses->seq_length[ref_id] = args.seq_size;
        //Record the range size.
        enzyme_poses->enz_range_size = static_cast<size_t *>(realloc(enzyme_poses->enz_range_size, sizeof(size_t) * (1 + ref_id)));
        enzyme_poses->enz_range_size[ref_id] = counter;
        //Append the range.
        enzyme_poses->enz_ranges = static_cast<ENZYME_RANGE **>(realloc(enzyme_poses->enz_ranges, sizeof(ENZYME_RANGE *) * (1 + ref_id)));
        enzyme_poses->enz_ranges[ref_id] = pos_ranges;
        //Increase the total.
        ++enzyme_poses->seq_total;
    }
}

void filter_fasta_dump_enzyme_count(const char *filepath, const FASTA_ENZYME_POSES &poses)
{
#ifdef _MSC_VER
    FILE *enz_count_file = NULL;
    fopen_s(&enz_count_file, filepath, "w");
#else
    FILE *enz_count_file = fopen(filepath, "w");
#endif
    //Loop and write the file.
    fprintf(enz_count_file, "%zu\n", poses.seq_total);
    for(size_t i=0; i<poses.seq_total; ++i)
    {
        fprintf(enz_count_file, "%s\t%zu\t%zu\n",
                poses.seq_names[i],
                poses.enz_range_size[i],
                poses.seq_length[i]);
    }
    fclose(enz_count_file);
}
