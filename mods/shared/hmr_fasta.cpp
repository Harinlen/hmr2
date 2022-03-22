#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>

#include "hmr_text_file.h"
#include "hmr_bin_file.h"
#include "hmr_path.h"
#include "hmr_ui.h"

#include "hmr_fasta.h"

typedef struct FASTA_PARSE
{
    char *seq_name, *seq_data;
    size_t seq_name_len, seq_data_len;
    FASTA_PROC user_parser;
    void *user;
} FASTA_PARSE;

inline void fasta_yield_line(char *line, size_t line_size, FASTA_PARSE &args, int32_t &index)
{
    if(line_size == 0 || line == NULL)
    {
        return;
    }
    //Check whether the line start is with '>'.
    if(line[0] == '>')
    {
        //Check the last is empty or not.
        if(args.seq_name != NULL)
        {
            if(args.seq_data != NULL)
            {
                //Yield the line parser.
                args.user_parser(index, args.seq_name, args.seq_name_len, args.seq_data, args.seq_data_len, args.user);
                ++index;
            }
            else
            {
                //Clear the sequence name.
                free(args.seq_name);
            }
        }
        //Set the sequence name.
        args.seq_name = static_cast<char *>(malloc(line_size));
        assert(NULL != args.seq_name);
        args.seq_name_len = line_size-1;
        memcpy(args.seq_name, line+1, args.seq_name_len);
        args.seq_name[args.seq_name_len] = '\0';
        //Reset the sequence.
        args.seq_data = NULL;
        args.seq_data_len = 0;
    }
    else
    {
        //Check the sequence.
        size_t seq_extend_len = args.seq_data_len + line_size;
        args.seq_data = static_cast<char *>(realloc(args.seq_data, seq_extend_len + 1));
        //Copy the data.
        assert(args.seq_data + args.seq_data_len);
        memcpy(args.seq_data + args.seq_data_len, line, line_size);
        args.seq_data_len = seq_extend_len;
        args.seq_data[args.seq_data_len] = '\0';
    }
}

void hmr_fasta_read(const char *filepath, FASTA_PROC parser, void *user)
{
    FASTA_PARSE args;
    TEXT_LINE_HANDLE line_handle;
    if(!text_open_read_line(filepath, &line_handle))
    {
        time_error(-1, "Failed to read FASTA file %s", filepath);
    }
    //Loop and detect line.
    char *line = NULL;
    size_t len = 0, line_length = 0;
    int32_t index = 0;
    ssize_t line_size = 0;
    args.seq_name = NULL;
    args.seq_data = NULL;
    args.seq_name_len = 0;
    args.seq_data_len = 0;
    args.user_parser = parser;
    args.user = user;
    while((line_size = (line_handle.parser(&line, &len, &line_handle.buf, line_handle.file_handle))) != -1)
    {
        //Trimmed line.
        line_length = static_cast<size_t>(line_size);
        trimmed_right(line, line_length);
        //Yield the line, call the function.
        fasta_yield_line(line, line_length, args, index);
    }
    //At the end of the line, yield the last result.
    parser(index, args.seq_name, args.seq_name_len, args.seq_data, args.seq_data_len, user);
    //Close the file.
    text_close_read_line(&line_handle);
}

std::string hmr_fasta_path_contig(const char* filepath)
{
    return path_basename(filepath) + ".hmr_contig";
}

bool hmr_fasta_load_contig(const char* filepath, std::vector<HMR_CONTIG>& contigs)
{
    //Open the contig input file to read the data.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "rb"))
    {
        return false;
    }
    //Read the length of the contigs.
    size_t contig_sizes = 0;
    fread(&contig_sizes, sizeof(size_t), 1, contig_file);
    contigs = std::vector<HMR_CONTIG>();
    contigs.reserve(contig_sizes);
    for (size_t i = 0; i < contig_sizes; ++i)
    {
        //Read the contig from the file.
        HMR_CONTIG contig;
        fread(&contig.seq_name_length, sizeof(size_t), 1, contig_file);
        fread(&contig.seq_name, sizeof(char), contig.seq_name_length, contig_file);
        fread(&contig.seq_length, sizeof(size_t), 1, contig_file);
        contigs.push_back(contig);
    }
    fclose(contig_file);
    return true;
}

bool hmr_fasta_save_contig(const char* filepath, const std::vector<HMR_CONTIG>& contigs)
{
    //Open the contig output file to write the data.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "wb"))
    {
        return false;
    }
    //Write the length of the contigs.
    size_t contig_sizes = contigs.size();
    fwrite(&contig_sizes, sizeof(size_t), 1, contig_file);
    for (const auto& contig : contigs)
    {
        //Write the contig information.
        fwrite(&contig.seq_name_length, sizeof(size_t), 1, contig_file);
        fwrite(&contig.seq_name, sizeof(char), contig.seq_name_length, contig_file);
        fwrite(&contig.seq_length, sizeof(size_t), 1, contig_file);
    }
    fclose(contig_file);
    return true;
}
