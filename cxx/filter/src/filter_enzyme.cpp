#include <cstring>
#include <unordered_map>

#include "hmr_ui.h"

#include "filter_enzyme.h"

#ifdef _MSC_VER
#define sprintf     sprintf_s
#endif

std::unordered_map<std::string, const char *> known_enzyme_alias
{
    {std::string("HINDIII"), "AAGCTT"},
    {std::string("HIND3"), "AAGCTT"},
    {std::string("NCOI"), "CCATGG"},
    {std::string("DPN1"), "GATC"},
    {std::string("MBOI"), "GATC"},
};

void filter_enzyme_formalize(char *enzyme, const char **nuc_seq, int *nuc_seq_size)
{
    size_t length = strlen(enzyme);
    //Convert the original char in upper case letter.
    for(char *s=enzyme, *e = enzyme + length; s<e; ++s)
    {
        if((*s) >= 'a' && (*s) <= 'z')
        {
            (*s) = 'A' + ((*s) - 'a');
        }
    }
    //Check whether it is an known alias name.
    auto known_finder = known_enzyme_alias.find(std::string(enzyme));
    if(known_finder != known_enzyme_alias.end())
    {
        (*nuc_seq) = known_finder->second;
        (*nuc_seq_size) = static_cast<int>(strlen(*nuc_seq));
        return;
    }
    //Or else we have to check the enzyme.
    for(char *s = enzyme, *e = enzyme + length; s<e; ++s)
    {
        //Check invalid nuc.
        if((*s) != 'A' && (*s) != 'C' && (*s) != 'T' && (*s) != 'G')
        {
            char buf[1024];
            sprintf(buf, "Invalid nucleotide base '%c' found in enzyme.", *s);
            time_print(buf);
            time_error(-1, "Make sure your enzyme should be only composed of 'A', 'C', 'T' and 'G'.");
        }
    }
    *nuc_seq = enzyme;
    *nuc_seq_size = static_cast<int>(length);
}
