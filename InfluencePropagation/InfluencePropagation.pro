TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    mygraphlib.cpp \
    mycolor.cpp \
    geompack.cpp \
    myutils.cpp \
    columngenerator.cpp \
    glcipinstance.cpp

HEADERS += \
    main.h \
    mygraphlib.h \
    mycolor.h \
    geompack.h \
    myutils.h \
    easyscip.h \
    columngenerator.h \
    star.h \
    glcipinstance.h
