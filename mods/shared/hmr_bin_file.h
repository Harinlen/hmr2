#ifndef HMR_BIN_FILE_H
#define HMR_BIN_FILE_H

#include <cstdio>

bool bin_open(const char* filepath, FILE** file, char const* mode);

#endif // HMR_BIN_FILE_H