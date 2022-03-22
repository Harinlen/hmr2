#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "hmr_text_file.h"
#include "hmr_ui.h"
#include "hmr_fasta_type.h"
#include "hmr_fastq.h"

inline bool fetch_line(TEXT_GETLINE parser, void* handle, char** line, size_t* length)
{
    ssize_t result;
    while ((result = parser(line, length, handle)) != -1)
    {
        if (result > 0)
        {
            *line = static_cast<char*>(realloc(*line, result));
            *length = result;
            return true;
        }
    }
    return false;
}

inline bool fetch_name(TEXT_GETLINE parser, void* handle, char** line, size_t* length)
{
    //Loop until get a line starts with '@'.
    while (fetch_line(parser, handle, line, length))
    {
        //Check whether the line starts with '@'.
        if ((*line)[0] == '@')
        {
            return true;
        }
    }
    //If we goes here, then it means we are failed to fetch a line.
    return false;
}

bool fastq_fetch_block(TEXT_GETLINE parser, void *handle, FASTQ_READ *read)
{
    char *block_name = NULL;
    size_t block_name_length = 0;
    /*
    * A FASTQ file containing a single sequence might look like this:
    * 
    *   @SEQ_ID
    *   GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
    *   +
    *   !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
    * 
    */
    if (!fetch_name(parser, handle, &block_name, &block_name_length))
    {
        //Failed to fetch the block.
        return false;
    }
    //Read the sequence of the read.
    char* seq = NULL, * split = NULL, * quality = NULL;
    size_t seq_length = 0, split_length = 0, quality_length = 0;
    //Fetch the block lines.
    if (!fetch_line(parser, handle, &seq, &seq_length)) { return false; }
    if (!fetch_line(parser, handle, &split, &split_length) || split[0] != '+') { return false; }
    if (!fetch_line(parser, handle, &quality, &quality_length)) { return false; }
    //Trim the lines.
    trimmed_right(block_name, block_name_length);
    trimmed_right(seq, seq_length);
    trimmed_right(split, split_length);
    trimmed_right(quality, quality_length);
    //Then length of quality should equal to sequence.
    if (seq_length != quality_length) { return false; }
    //Assign the data to read result.
    *read = FASTQ_READ{block_name_length, seq_length, split_length, block_name, seq, quality, split};
    return true;
}

void hmr_fastq_pair_read(const char* left_filepath, const char* right_filepath, FASTQ_PAIR_PROC parser, void* user)
{
    TEXT_GETLINE left_parser, right_parser;
    void* left_handle, * right_handle;
    //Read the left and right FASTQ file.
    if (!text_open(left_filepath, &left_parser, &left_handle))
    {
        time_error(-1, "Failed to read FASTQ file %s", left_filepath);
    }
    if (!text_open(right_filepath, &right_parser, &right_handle))
    {
        time_error(-1, "Failed to read FASTQ file %s", right_filepath);
    }
    //Loop and fetch FASTQ block from left side file, and then fetch from the right side.
    FASTQ_READ left_read, right_read;
    bool keep_parsing = true;
    while (keep_parsing)
    {
        int left_result = fastq_fetch_block(left_parser, left_handle, &left_read),
            right_result = fastq_fetch_block(right_parser, right_handle, &right_read);
        if (left_result && right_result)
        {
            //Call the parser to handle the read pair.
            parser(FASTQ_READ_PAIR{ left_read, right_read }, user);
        }
        else
        {
            //When both reaches the end, mission complete, finished.
            if ((!left_result) && (!right_result))
            {
                keep_parsing = false;
                continue;
            }
            //Check which side of the file is not paired.
            time_error(-1, "Failed to fetch read from %s", left_result ? right_filepath : left_filepath);
        }
    }
    //Close the file.
    text_close(left_parser, left_handle);
    text_close(right_parser, right_handle);
}

void hmr_fastq_read_free(const FASTQ_READ& read)
{
    free(read.name);
    free(read.seq);
    free(read.split);
    free(read.quality);
}
