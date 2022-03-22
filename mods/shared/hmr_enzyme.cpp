#include <unordered_map>

#include "hmr_ui.h"
#include "hmr_bin_file.h"

#include "hmr_enzyme.h"

std::unordered_map<std::string, const char*> known_enzyme_alias
{
    {std::string("HINDIII"), "AAGCTT"},
    {std::string("HIND3"), "AAGCTT"},
    {std::string("NCOI"), "CCATGG"},
    {std::string("DPN1"), "GATC"},
    {std::string("MBOI"), "GATC"},
};

void hmr_enzyme_formalize(char* enzyme, const char** nuc_seq, int* nuc_seq_size)
{
    size_t length = strlen(enzyme);
    //Convert the original char in upper case letter.
    for (char* s = enzyme, *e = enzyme + length; s < e; ++s)
    {
        if ((*s) >= 'a' && (*s) <= 'z')
        {
            (*s) = 'A' + ((*s) - 'a');
        }
    }
    //Check whether it is an known alias name.
    auto known_finder = known_enzyme_alias.find(std::string(enzyme));
    if (known_finder != known_enzyme_alias.end())
    {
        (*nuc_seq) = known_finder->second;
        (*nuc_seq_size) = static_cast<int>(strlen(*nuc_seq));
        return;
    }
    //Or else we have to check the enzyme is valid or not.
    for (char* s = enzyme, *e = enzyme + length; s < e; ++s)
    {
        //Check invalid nuc.
        if ((*s) != 'A' && (*s) != 'C' && (*s) != 'T' && (*s) != 'G')
        {
            time_error(-1, "Invalid nucleotide base '%c' found in enzyme '%s'.", *s, enzyme);
        }
    }
    *nuc_seq = enzyme;
    *nuc_seq_size = static_cast<int>(length);
}

bool hmr_enzyme_load_range(const char* filepath, std::vector<CONTIG_ENZYME_RANGE>& contig_ranges)
{
    FILE* enzyme_range_file;
    if (!bin_open(filepath, &enzyme_range_file, "rb"))
    {
        return false;
    }
    //Read the contig size.
    size_t contig_size;
    fread(&contig_size, sizeof(size_t), 1, enzyme_range_file);
    contig_ranges = std::vector<CONTIG_ENZYME_RANGE>();
    contig_ranges.reserve(contig_size);
    for (size_t i = 0; i < contig_size; ++i)
    {
        CONTIG_ENZYME_RANGE range;
        fread(&range.range_size, sizeof(size_t), 1, enzyme_range_file);
        range.ranges = static_cast<ENZYME_RANGE*>(malloc(sizeof(ENZYME_RANGE) * range.range_size));
        for (size_t j = 0; j < range.range_size; ++j)
        {
            fread(&range.ranges[j].start, sizeof(size_t), 1, enzyme_range_file);
            fread(&range.ranges[j].end, sizeof(size_t), 1, enzyme_range_file);
        }
    }
    fclose(enzyme_range_file);
    return true;
}

bool hmr_enzyme_save_range(const char* filepath, const std::vector<CONTIG_ENZYME_RANGE>& enzyme_ranges)
{
    FILE* enzyme_range_file;
    if (!bin_open(filepath, &enzyme_range_file, "wb"))
    {
        return false;
    }
    //Write the contig size.
    size_t contig_sizes = enzyme_ranges.size();
    fwrite(&contig_sizes, sizeof(size_t), 1, enzyme_range_file);
    for (const auto& contig : enzyme_ranges)
    {
        //Write the contig information.
        fwrite(&contig.range_size, sizeof(size_t), 1, enzyme_range_file);
        for (size_t i=0; i<contig.range_size; ++i)
        {
            const auto& range = contig.ranges[i];
            fwrite(&range.start, sizeof(size_t), 1, enzyme_range_file);
            fwrite(&range.end, sizeof(size_t), 1, enzyme_range_file);
        }
    }
    fclose(enzyme_range_file);
    return true;
}
