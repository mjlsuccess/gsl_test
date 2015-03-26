#-------------------------------------------------
#
# Project created by QtCreator 2015-03-23T09:23:49
#
#-------------------------------------------------

QT       += core gui opengl

QT += widgets

TARGET = gsl_test
TEMPLATE = app

INCLUDEPATH +=/home/pku-m/SDK/GLViewer/include

LIBS += /usr/local/lib/libgsl.a
LIBS += /usr/local/lib/libgslcblas.a

LIBS += /home/pku-m/SDK/GLViewer/lib/libGLViewer_Debug.a

LIBS += -lGLU

SOURCES += main.cpp\
        mainwindow.cpp \
    frenettraj.cpp \
    frenet2cartesian.cpp

HEADERS  += mainwindow.h \
    frenettraj.h \
    frenet2cartesian.h

FORMS    += mainwindow.ui

DISTFILES += \
    test.py \
    ../../下载/down.png
