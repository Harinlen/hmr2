#include <cstring>

#include "hmr_bin_file.h"

#include "hmr_path.h"

std::string path_suffix(const char *filepath, size_t length)
{
    //When the length is 0, auto detect the length of the path.
    if (length == 0)
    {
        length = strlen(filepath);
    }
    std::string filepath_str(filepath, length);
    //Find the dot from the last one.
    auto dot_pos = filepath_str.find_last_of('.');
    if(dot_pos == std::string::npos)
    {
        return std::string();
    }
    //Or else, extract the last string.
    return filepath_str.substr(dot_pos);
}

std::string path_basename_core(const char* filepath, size_t length)
{
    //When the length is 0, auto detect the length of the path.
    if (length == 0)
    {
        length = strlen(filepath);
    }
    std::string filepath_str(filepath, length);
    //Find the dot from the last one.
    auto dot_pos = filepath_str.find_last_of('.');
    if (dot_pos == std::string::npos)
    {
        return std::string();
    }
    //Extract the first part of the string.
    return filepath_str.substr(0, dot_pos);
}

std::string path_basename(const char* filepath, size_t length)
{
    //If the suffix is gz, remove the gz.
    if (path_suffix(filepath, length) == ".gz")
    {
        return path_basename_core(filepath, strlen(filepath) - 3);
    }
    //Just get the basename.
    return path_basename(filepath, length);
}

bool path_can_read(const char* filepath)
{
    //Try to open the file to have a test.
    FILE* file_test;
    bool result = bin_open(filepath, &file_test, "rb");
    //If we open it successfully, close it.
    if (result)
    {
        fclose(file_test);
    }
    return result;
}
