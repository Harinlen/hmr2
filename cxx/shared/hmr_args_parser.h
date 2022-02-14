#ifndef HMR_ARGS_PARSER_H
#define HMR_ARGS_PARSER_H

#include <unordered_map>
#include <functional>
#include <set>
#include <string>

typedef struct HMR_ARG_INFO
{
    std::string val_name;
    std::set<std::string> alias;
    std::string description;
    std::function<void(char *)> set_val;
} HMR_ARG_INFO;

typedef std::unordered_map<std::string, HMR_ARG_INFO> HMR_ARG_PARSER;

typedef struct HMR_ARG_COND
{
    std::function<bool(void)> check;
    std::string error_str;
} HMR_ARG_COND;
typedef std::vector<HMR_ARG_COND> HMR_CONSTRAINTS;

#endif // HMR_ARGS_PARSER_H
