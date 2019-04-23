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
    glcipinstance.cpp \
    glcipsolution.cpp \
    covmodel.cpp \
    arcmodel.cpp \
    cyclecutsgenerator.cpp \
    graphviewer.cpp

HEADERS += \
    main.h \
    mygraphlib.h \
    mycolor.h \
    geompack.h \
    myutils.h \
    easyscip.h \
    columngenerator.h \
    star.h \
    glcipinstance.h \
    glcipsolution.h \
    glcipsolution.h \
    covmodel.h \
    arcmodel.h \
    cyclecutsgenerator.h \
    graphviewer.h
