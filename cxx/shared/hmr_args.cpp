#include <cstdio>
#include <cstring>

#include "hmr_args_type.h"
#include "hmr_args_parser.h"

#include "hmr_args.h"

extern HMR_ARG_PARSER args_parser;
extern HMR_CONSTRAINTS args_constrains;

void show_argument(const char *arg, const HMR_ARG_INFO &arg_info)
{
    char arg_buf[1024], arg_item[1024];
    const char *val_name = arg_info.val_name.c_str();
    sprintf(arg_buf, "%s %s", arg, val_name);
    for(auto i=arg_info.alias.begin(); i!=arg_info.alias.end(); ++i)
    {
        sprintf(arg_item, ", %s %s", (*i).c_str(), val_name);
#ifdef _MSC_VER
        strcat_s(arg_buf, 1024, arg_item);
#else
        strcat(arg_buf, arg_item);
#endif
    }
    printf("  %s", arg_buf);
    size_t arg_len = strlen(arg_buf);
    if(arg_len < 21)
    {
        //Fill to 22 spacing.
        printf("%*s", static_cast<int>(20 - arg_len), "");
    }
    else
    {
        printf("\n%*s  ", 20, "");
    }
    printf("  %s\n", arg_info.description.c_str());
}

void help_exit(char *prog_name)
{
    printf("usage: %s [-h]", prog_name);
    for(auto i=args_parser.begin(); i!=args_parser.end(); ++i)
    {
        auto arg_info = i->second;
        //Print the header and value name.
        printf(" %s %s", i->first.c_str(), arg_info.val_name.c_str());
    }
    printf("\n");
    printf("optional arguments:\n"
           "  -h, --help            Show this help message and exit\n");
    for(auto i=args_parser.begin(); i!=args_parser.end(); ++i)
    {
        //Print the header and value name.
        show_argument(i->first.c_str(), i->second);
    }
    exit(0);
}

void parse_arguments(int argc, char *argv[])
{
    //Loop and check all the arguments.
    for(int i=1; i<argc; i+=2)
    {
        std::string arg(argv[i]);
        //Check the help arguments.
        if(std::string("-h") == arg || std::string("--help") == arg) { help_exit(argv[0]); }
        //When find the argument, run the program.
        for(auto j=args_parser.begin(); j!=args_parser.end(); ++j)
        {
            //Call the set function.
            auto arg_info = j->second;
            if(j->first == arg || arg_info.alias.find(arg) != arg_info.alias.end())
            {
                //Set the value.
                arg_info.set_val(argv[i+1]);
            }
        }
    }
    //Check whether the arguments meet the conditions.
    for(auto i=args_constrains.begin(); i!=args_constrains.end(); ++i)
    {
        //When the condition is not meet.
        if(!(*i).check())
        {
            printf("%s\n", (*i).error_str.c_str());
            help_exit(argv[0]);
        }
    }
}
