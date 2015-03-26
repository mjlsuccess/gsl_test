#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <cstdio>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    glview = new GLViewer;

    shownum = 5;
    glview->makeCurrent();
    glview->clearDisplayList();
    displaylistbase = glGenLists(shownum);
    for(unsigned int i=0; i<shownum; i++)
        glview->addDisplayList(displaylistbase+i);

    ui->scrollArea->setWidget(glview);

    timer.start(100);
    //connect(&timer, SIGNAL(timeout()), this, SLOT(on_timer()));
    connect(ui->plotOnlineButton, SIGNAL(clicked()), this, SLOT(on_plotOnline()));

    //load image
    image.load("down.png");
    pixeldata = new GLubyte [image.byteCount()];
    pixeldata = image.bits();



    printf("%d %d\n", image.byteCount(), image.bytesPerLine());
}

MainWindow::~MainWindow()
{
    glview->clearDisplayList();
    if(glview->parent() ==NULL)
        delete glview;
    delete ui;

    delete pixeldata;
}

int cmp(const void *a, const void *b)
{
    Shift *aa = (Shift*)a;
    Shift *bb = (Shift*)b;
    return (aa->cost - bb->cost);
}

void MainWindow::on_plotOnline()
{
    glNewList(displaylistbase,GL_COMPILE_AND_EXECUTE);
   // latTest();
    //velTest();
   // glNewList(displaylistbase,GL_COMPILE_AND_EXECUTE);
    //glDrawPixels(image.width(), image.height(), GL_BGRA, GL_UNSIGNED_BYTE, pixeldata);
    //combineTest();
    Point2D startPt;
//    startPt.x = 2;//-2 -> circle
//    startPt.y = 2;//0 -> circle
    startPt.x = 60;// cos
    startPt.y = 10;//
    frenet2Cart_test.generateOfflineTraj();
    frenet2Cart_test.generateOnlineTraj(startPt, traj);

    glPushMatrix();
    glPointSize(10);
    glColor3f(1,0,0);
    glBegin(GL_POINTS);
    glVertex3f(startPt.x, startPt.y, 0);
    glEnd();
    glPopMatrix();

    glEndList();
}

void MainWindow::on_timer()
{

}
void MainWindow::latTest()
{

    double t0 = 0;
    traj.lat(t0, 2, 0.1, -0.1);
    traj.plotLat(t0, traj.Clat);

//    double gap=1;
//    double d0[3];
//    double *dt=traj.Clat[0].quintic;
//    gsl_poly_eval_derivs(dt,6,gap,d0,3); //计算Clat[0]在1秒处的位置，速度，加速度
//    traj.lat(t0+gap,d0[0],d0[1],d0[2]);
//    qsort(traj.Clat,traj.dNum*traj.tNum,sizeof(Shift),cmp);
//    glColor3f(1,0,0);
//    glLineWidth(10.0);
//    traj.plotBestLat(t0+gap,traj.Clat[0].end_t,traj.Clat[0].quintic);
//    glColor3f(1,0.5,0.5);
//    traj.plotLat(t0+gap,traj.Clat);

}

void MainWindow::velTest()
{
   // glNewList(displaylistbase,GL_COMPILE_AND_EXECUTE);

    int flag=1;//0为轨迹，1为速度

    double t0=0;

    traj.vel(t0,0,15,-0.1);
    qsort(traj.Clon,traj.vNum*traj.tNum,sizeof(Shift),cmp);

    glColor3f(0.5,0.5,0.8);
    traj.plotVel(t0,traj.Clon,flag);

//    glColor3f(0,1,0);
//    glLineWidth(4.0);
//    traj.plotBestVel(t0,traj.Clon[0].end_t,traj.Clon[0].quintic,flag);
    //draw picture
    glDrawPixels(image.width(), image.height(), GL_BGRA, GL_UNSIGNED_BYTE, pixeldata);
   // glEndList();
}

void MainWindow::combineTest()
{
    double t0, d0, d0_1, d0_2, s0, s0_1, s0_2;
    t0 = 0;
    d0 = 2;
    d0_1 = 0.1;
    d0_2 = 0.04;

    s0 = 0;
    s0_1 = 2;
    s0_2 = 1;

    traj.combine(t0, d0, d0_1, d0_2, s0, s0_1, s0_2);

    glColor3f(0.5,0.5,0.8);
    traj.plotCombine(t0, traj.traj);
}
