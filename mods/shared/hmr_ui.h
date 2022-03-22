#ifndef HMR_UI_H
#define HMR_UI_H

/*!
 * \brief Print the text with the current time at the beginning.
 * \param str The text to be printed.
 */
void time_print(const char *str, ...);

/*!
 * \brief Exit the program with an error code, and print an error info.
 * \param exitCode The exit code to be set.
 * \param str The text to be printed.
 */
void time_error(int exitCode, const char *str, ...);

#endif // HMR_UI_H
