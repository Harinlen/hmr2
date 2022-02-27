#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>

#include "hmr_ui.h"

#include "hmr_fasta.h"

#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;

ssize_t getline(char **buf, size_t *bufsiz, FILE *fp)
{
    char *ptr, *eptr;
    //Pre-alloc memory for buffer.
    if (*buf == NULL || *bufsiz == 0)
    {
        //4K for initial.
        *bufsiz = 4096;
        if ((*buf = static_cast<char *>(malloc(*bufsiz))) == NULL)
        {
            return -1;
        }
    }
    for (ptr = *buf, eptr = *buf + *bufsiz;;)
    {
        //Fetch a char.
        int c = fgetc(fp);
        //EOF checking.
        if(c == -1)
        {
            if (feof(fp))
            {
                return ptr == *buf ? -1 : ptr - *buf;
            }
            return -1;
        }
        *ptr++ = c;
        //Reach the delimiter.
        if (c == '\n')
        {
            *ptr = '\0';
            return ptr - *buf;
        }
        //Enlarge the buffer by doubled.
        if (ptr + 2 >= eptr)
        {
            char *nbuf;
            size_t nbufsiz = *bufsiz * 2;
            ssize_t d = ptr - *buf;
            if ((nbuf = static_cast<char *>(realloc(*buf, nbufsiz))) == NULL)
            {
                return -1;
            }
            *buf = nbuf;
            *bufsiz = nbufsiz;
            eptr = nbuf + nbufsiz;
            ptr = nbuf + d;
        }
    }
}
#endif

static char *seq_name = NULL, *seq_data = NULL;
static size_t seq_name_len, seq_data_len;

inline int is_space(char c)
{
    return c=='\t' || c=='\n' || c=='\v' || c=='\f' || c=='\r' || c==' ';
}

inline void trimmed_right(char *buf, size_t &size)
{
    while(size > 0 && is_space(buf[size-1]))
    {
        buf[size--] = '\0';
    }
}

inline void fasta_yield_line(char *line, size_t line_size, FASTA_SEQ_PROC parser, void *user)
{
    if(line_size == 0 || line == NULL)
    {
        return;
    }
    //Check whether the line start is with '>'.
    if(line[0] == '>')
    {
        //Check the last is empty or not.
        if(seq_name != NULL)
        {
            if(seq_data != NULL)
            {
                //Yield the line parser.
                parser(seq_name, seq_name_len, seq_data, seq_data_len, user);
            }
            else
            {
                //Clear the sequence name.
                free(seq_name);
            }
        }
        //Set the sequence name.
        seq_name = static_cast<char *>(malloc(line_size));
        seq_name_len = line_size-1;
        memcpy(seq_name, line+1, seq_name_len);
        seq_name[seq_name_len] = '\0';
        //Reset the sequence.
        seq_data = NULL;
        seq_data_len = 0;
    }
    else
    {
        //Check the sequence.
        size_t seq_extend_len = seq_data_len + line_size;
        seq_data = static_cast<char *>(realloc(seq_data, seq_extend_len + 1));
        //Copy the data.
        memcpy(seq_data+seq_data_len, line, line_size);
        seq_data_len = seq_extend_len;
        seq_data[seq_data_len] = '\0';
    }
}

void fasta_parser(const char *filepath, FASTA_SEQ_PROC parser, void *user)
{
    //Open and load the file.
#ifdef _MSC_VER
    FILE *fasta_file = NULL;
    fopen_s(&fasta_file, filepath, "r");
#else
    FILE *fasta_file = fopen(filepath, "r");
#endif
    if(!fasta_file)
    {
        time_error_str(1, "Failed to open fasta file %s", filepath);
    }
    //Loop and detect line.
    char *line = NULL;
    size_t len = 0, line_length = 0;
    ssize_t line_size = 0;
    while((line_size = (getline(&line, &len, fasta_file))) != -1)
    {
        //Trimmed line.
        line_length = static_cast<size_t>(line_size);
        trimmed_right(line, line_length);
        //Yield the line, call the function.
        fasta_yield_line(line, line_length, parser, user);
    }
    //At the end of the line, yield the last result.
    parser(seq_name, seq_name_len, seq_data, seq_data_len, user);
    //Reset the name and length to the initial state.
    seq_name = NULL;
    seq_name_len = 0;
    seq_data = NULL;
    seq_data_len = 0;
    fclose(fasta_file);
}
