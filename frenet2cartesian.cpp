#include "frenet2cartesian.h"
#include <qmath.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_heapsort.h>

#include <qopengl.h>

#define HAVE_INLINE
Frenet2Cartesian::Frenet2Cartesian()
{
    //generateOfflineTraj();
}

void Frenet2Cartesian::generateOfflineTraj()
{
    //center (R, 0), x = R + Rcos(theta), y = Rsin(theta)
//    double R = 50;
//    double resolution = 5800;
//    int i;
//    double theta = 0;
//    Point2D point;
//    offlineTraj.clear();

//    glPushMatrix();
//    glLineWidth(4);
//    glColor3f(1,0,0);
//    glBegin(GL_LINE_STRIP);

//    for(i=0; i<=resolution; i++)
//    {
//        theta =i * M_PI/resolution;
//        point.x = R + R*cos(theta);
//        point.y = R*sin(theta);
//        offlineTraj.push_front(point);
//        glVertex3f(point.x, point.y, 0);
//    }

//    glEnd();
//    glPopMatrix();
//    for(i=0; i<offlineTraj.size(); i++)
//        offlineTraj[i].index = i;

    ////////////////////////////////////
    /// \brief test cos - function
    //////////
    double T = 60;
    double resolution = 4800;
    int i;
    Point2D point;
    offlineTraj.clear();

    glPushMatrix();
    glLineWidth(4);
    glColor3f(1,0,0);
    glBegin(GL_LINE_STRIP);

    for(i=0; i<=resolution; i++)
    {
        point.x  = i*T/resolution;
        point.y = 20*cos((2*M_PI/T)*point.x) + 10;
        offlineTraj.push_front(point);
        glVertex3f(point.x, point.y, 0);
    }
    glEnd();
    glPopMatrix();

    for(i=0; i<offlineTraj.size(); i++)
        offlineTraj[i].index = i;
}

inline double dist(const Point2D& p0,const Point2D& p1)
{
    return sqrt((p0.x-p1.x)*(p0.x-p1.x) + (p0.y-p1.y)*(p0.y - p1.y));
}

Point2D Frenet2Cartesian::findNearestPoint(const Point2D &startPt)
{
    int size = offlineTraj.size();
    double distance = INT_MAX, tmp;
    int index = 0;
    for (int i=0; i<size; i++)
    {
        tmp = dist(startPt, offlineTraj[i]);
        if(distance > tmp)
        {
            distance = tmp;
            index = i;
        }
    }
    return offlineTraj[index];
}
//给定startPt 和 弧长，计算弧线的终点, startPt在弧线上
Point2D Frenet2Cartesian::findArcEndPoint(const Point2D &startPt, double arcLength)
{
    int size = offlineTraj.size();
    int i = startPt.index >= (size-1) ? (size -1)  : startPt.index;
    double distSum = 0;
    int index = i;
    for(; i<(size-1); i++ )
    {
        distSum += dist(offlineTraj[i], offlineTraj[i+1]);
        if(distSum >= arcLength)
        {
            index = i+1;
            break;
        }
    }

    return offlineTraj[index];
}

void Frenet2Cartesian::generateOnlineTraj(const Point2D &startPt, FrenetTraj &frenetTraj)
{
    double t0, d0, d0_1, d0_2, s0, s0_1, s0_2;
    //find nearest point on offline trajectory
    Point2D nearPt = findNearestPoint(startPt);
    Point2D nearPt_1;
    d0 = dist(startPt, nearPt);
    //判断方向
    nearPt_1 = nearPt.index >= (offlineTraj.size() - 1) ? offlineTraj[nearPt.index -1] : offlineTraj[nearPt.index+1];
    double k = (nearPt.y - nearPt_1.y)/(nearPt.x - nearPt_1.x);
    double direct = nearPt.y+k*(startPt.x - nearPt.x) - startPt.y; //判断初始点在曲线的上方还是下方
    if(direct > 0) //在下方
        d0 = -d0;

    d0_1 = 0.6;
    d0_2 = 0.1;

    s0 = 0;
    s0_1 = 1;
    s0_2 = 0.2;

    t0 = 0;

//    frenetTraj.lat(t0, d0, d0_1, d0_2);
//    frenetTraj.vel(t0, s0, s0_1, s0_2);
    frenetTraj.combine(t0, d0, d0_1, d0_2, s0, s0_1, s0_2);

    double x,y, ratio, theta;
    Traj* tj = frenetTraj.traj;

    int num = frenetTraj.tNum*frenetTraj.dNum*frenetTraj.vNum;

    glPushMatrix();
    glLineWidth(1);
    glColor3f(0,1,0);
    for (int i=0;i<num;i++)
    {
        glBegin(GL_LINE_STRIP);
        for(double t=t0;t<tj[i].end_t+t0;t+=0.01)
        {
            double arcLength = gsl_poly_eval(tj[i].Lon.quintic, 6, t);
            double latLength = gsl_poly_eval(tj[i].Lat.quintic, 6, t);

            //find the arc end point
            Point2D arcEndPt = findArcEndPoint(nearPt, arcLength);
            //calc the tangent of arcEndPoint
            Point2D nextPt ;
            if(arcEndPt.index >= (offlineTraj.size()-1))
                nextPt = offlineTraj[arcEndPt.index -1];
            else
                nextPt = offlineTraj[arcEndPt.index+1];

            ratio = (arcEndPt.y - nextPt.y)/(arcEndPt.x - nextPt.x + 1e-22);
            theta = atan(ratio);

            x = arcEndPt.x - latLength*sin(theta);
            y = arcEndPt.y + latLength*cos(theta);

            glVertex3f(x, y, 0);
        }
        glEnd();
    }

    glPopMatrix();
}

















