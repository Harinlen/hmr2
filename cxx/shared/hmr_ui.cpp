#include <cstdio>
#include <cstdlib>
#include <ctime>
#ifdef _MSC_VER
#include <Windows.h>
#endif

#include "hmr_ui.h"

#ifdef _MSC_VER
HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
#endif

void time_print(const char *str)
{
    time_t now = time(0);
#ifdef _MSC_VER
    struct tm tstruct;
    localtime_s(&tstruct, &now);
    SetConsoleTextAttribute(hStdOut, 2);
    printf("[%02d:%02d:%02d]",
           tstruct.tm_hour, tstruct.tm_min, tstruct.tm_sec);
    SetConsoleTextAttribute(hStdOut, 7);
    printf(" %s\n", str);
#else
    struct tm *tstruct = localtime(&now);
    printf("\033[32m[%02d:%02d:%02d]\033[0m %s\n",
           tstruct->tm_hour, tstruct->tm_min, tstruct->tm_sec, str);
#endif
}

void time_error(int exitCode, const char *str)
{
    //Print the data.
    time_print(str);
    //Exit the program.
    exit(exitCode);
}

#define TIME_BUF_ERROR  { \
    char buf[1024]; \
    sprintf(buf, fmt_str, value); \
    time_error(exitCode, buf); \
}

#define TIME_BUF_PRINT  { \
    char buf[1024]; \
    sprintf(buf, fmt_str, value); \
    time_print(buf); \
}

void time_print_int(const char *fmt_str, int value)
TIME_BUF_PRINT

void time_print_float(const char *fmt_str, float value)
TIME_BUF_PRINT

void time_print_str(const char *fmt_str, const char *value)
TIME_BUF_PRINT

void time_print_size(const char *fmt_str, size_t value)
TIME_BUF_PRINT

void time_error_int(int exitCode, const char *fmt_str, int value)
TIME_BUF_ERROR

void time_error_str(int exitCode, const char *fmt_str, const char *value)
TIME_BUF_ERROR

void time_error_size(int exitCode, const char *fmt_str, size_t value)
TIME_BUF_ERROR
