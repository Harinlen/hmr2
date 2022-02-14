TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH = ../shared

HEADERS += \
    ../shared/hmr_fasta.h \
    ../shared/hmr_ui.h

SOURCES += \
    ../shared/hmr_fasta.cpp \
    ../shared/hmr_ui.cpp \
    src/main.cpp
