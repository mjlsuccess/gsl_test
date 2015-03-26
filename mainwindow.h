#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <glviewer.h>
#include <frenettraj.h>
#include <GL/glu.h>
#include <qtimer.h>
#include <QImage>
#include "frenet2cartesian.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
public:
    Frenet2Cartesian frenet2Cart_test;
    FrenetTraj traj;
    GLViewer *glview;
    GLuint shownum;
    GLuint displaylistbase;
    QImage image;
public:
    void latTest();
    void velTest();
    void combineTest();
public slots:
    void on_timer();
    void on_plotOnline();
private:
    Ui::MainWindow *ui;
    QTimer timer;
    GLubyte *pixeldata;
};

#endif // MAINWINDOW_H
