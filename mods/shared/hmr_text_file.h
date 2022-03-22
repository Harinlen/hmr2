#ifndef HMR_TEXT_FILE_H
#define HMR_TEXT_FILE_H

#include <string>

inline int is_space(char c)
{
    return c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r' || c == ' ';
}

inline void trimmed_right(char* buf, size_t& size)
{
    while (size > 0 && is_space(buf[size - 1]))
    {
        buf[--size] = '\0';
    }
}

#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

typedef struct TEXT_LINE_BUF
{
    char* buf;
    size_t buf_size, reserved, offset;
} TEXT_LINE_BUF;

typedef ssize_t(*TEXT_GETLINE)(char** line, size_t* line_size, TEXT_LINE_BUF *buf, void* file_handle);

typedef struct TEXT_LINE_HANDLE
{
    TEXT_GETLINE parser;
    void* file_handle;
    TEXT_LINE_BUF buf;
} TEXT_LINE_HANDLE;

std::string text_open(const char *filepath, void** handle);
bool text_open_line(const char* filepath, TEXT_LINE_HANDLE *handle);
void text_close_line(TEXT_LINE_HANDLE* handle);

#endif // HMR_TEXT_FILE_H
