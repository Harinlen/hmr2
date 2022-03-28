#include <unordered_map>

#include "hmr_path.h"
#include "hmr_ui.h"
#include "hmr_text_file.h"

#include "allele.h"

typedef std::unordered_map<std::string, int32_t> CONTIG_NAME_ID;

CONTIG_NAME_ID build_contig_name_id(const HMR_CONTIGS& contigs)
{
    CONTIG_NAME_ID contig_ids;
    for (int32_t i = 0, i_max = static_cast<int32_t>(contigs.size()); i < i_max; ++i)
    {
        contig_ids.insert(std::make_pair(std::string(contigs[i].name), i));
    }
    return contig_ids;
}

inline int32_t contig_id(const CONTIG_NAME_ID& id_map, const std::string &name)
{
    auto id_finder = id_map.find(name);
    return id_finder == id_map.end() ? -1 : id_finder->second;
}

ALLELE_TABLE allele_ctg_table_load(const char* filepath, const HMR_CONTIGS& contigs, int32_t allele_groups)
{
    auto contig_id_map = build_contig_name_id(contigs);
    //Read this file line by line.
    TEXT_LINE_HANDLE table_handle;
    if (!text_open_read_line(filepath, &table_handle))
    {
        time_error(-1, "Failed to read allele contig table file %s", filepath);
    }
    char* line = NULL;
    size_t len = 0, line_length = 0;
    ssize_t line_size = 0;
    std::list<ALLELE_IDS> allele_rules;
    while ((line_size = (table_handle.parser(&line, &len, &table_handle.buf, table_handle.file_handle))) != -1)
    {
        //Get the line data.
        line_length = static_cast<size_t>(line_size);
        trimmed_right(line, line_length);
        if (line_length == 0)
        {
            continue;
        }
        //Loop and find all the '\t' inside the table.
        std::vector<std::string> columns;
        size_t last_pos = 0;
        for (size_t i = 0; i < line_length; ++i)
        {
            if (line[i] == '\t')
            {
                columns.push_back(std::string(line + last_pos, i - last_pos));
                last_pos = i + 1;
            }
        }
        if (last_pos != line_length - 1)
        {
            columns.push_back(std::string(line + last_pos, line_length - last_pos));
        }
        //We have to ignore the first two columns.
        if (columns.size() < 4)
        {
            continue;
        }
        //The column sizes must be less than allele groups.
        if (static_cast<int32_t>(columns.size()) > allele_groups + 2)
        {
            time_error(-1, "Allele contig line has more contigs (%zu) than allele groups (%d) in file %s", columns.size() - 2, allele_groups, filepath);
        }
        //Finding each name from the index column.
        ALLELE_IDS line_ids;
        line_ids.reserve(columns.size() - 2);
        for (auto i = columns.begin() + 2; i != columns.end(); ++i)
        {
            line_ids.push_back(contig_id(contig_id_map, *i));
        }
        allele_rules.push_back(line_ids);
    }
    //Close the file.
    text_close_read_line(&table_handle);
    return ALLELE_TABLE(allele_rules.begin(), allele_rules.end());
}

ALLELE_TABLE allele_load(const char* filepath, const HMR_CONTIGS& contigs, int32_t allele_groups)
{
    //Check the suffix of the filepath.
    std::string suffix = path_suffix(filepath);
    if (suffix == ".table")
    {
        //Read as .ctg.table format.
        return allele_ctg_table_load(filepath, contigs, allele_groups);
    }
    time_error(-1, "Unknown suffix for allele table file %s", suffix.data());
}
