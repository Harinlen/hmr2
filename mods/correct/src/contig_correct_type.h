#ifndef CONTIG_CORRECT_TYPE_H
#define CONTIG_CORRECT_TYPE_H

#include <cstdint>
#include <string>

#include <unordered_map>

typedef struct CONTIG_INFO
{
    int32_t id;
    size_t length;
} CONTIG_INFO;

typedef std::unordered_map<std::string, CONTIG_INFO> CONTIG_MAP;

#endif // CONTIG_CORRECT_TYPE_H