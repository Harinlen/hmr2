#ifndef HMR_UI_H
#define HMR_UI_H

/*!
 * \brief Print the text with the current time at the beginning.
 * \param str The text to be printed.
 */
void time_print(const char *str);

/*!
 * \brief Exit the program with an error code, and print an error info.
 * \param exitCode The exit code to be set.
 * \param str The text to be printed.
 */
void time_error(int exitCode, const char *str);

void time_print_int(const char *fmt_str, int value);

void time_print_float(const char *fmt_str, float value);

void time_print_str(const char *fmt_str, const char *value);

void time_print_size(const char *fmt_str, size_t value);

void time_error_int(int exitCode, const char *fmt_str, int value);

void time_error_str(int exitCode, const char *fmt_str, const char *value);

void time_error_size(int exitCode, const char *fmt_str, size_t value);

#endif // HMR_UI_H
