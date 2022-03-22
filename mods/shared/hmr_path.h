#ifndef HMR_PATH_H
#define HMR_PATH_H

#include <string>

std::string path_suffix(const char *filepath, size_t length = 0);
std::string path_basename(const char* filepath, size_t length = 0);

bool path_can_read(const char* filepath);

#endif // HMR_PATH_H
