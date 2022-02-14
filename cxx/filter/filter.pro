TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = \
    ../shared
    src

HEADERS += \
    ../shared/hmr_args.h \
    ../shared/hmr_args_parser.h \
    ../shared/hmr_fasta.h \
    ../shared/hmr_thread_pool.h \
    ../shared/hmr_ui.h \
    src/filter_enzyme.h \
    src/filter_fasta.h \
    src/filter_fasta_type.h \
    src/hmr_args_type.h

SOURCES += \
    ../shared/hmr_args.cpp \
    ../shared/hmr_fasta.cpp \
    ../shared/hmr_ui.cpp \
    src/filter_enzyme.cpp \
    src/filter_fasta.cpp \
    src/hmr_args_type.cpp \
    src/main.cpp
