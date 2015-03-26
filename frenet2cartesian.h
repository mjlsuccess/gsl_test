#ifndef FRENET2CARTESIAN
#define FRENET2CARTESIAN
#include <qvector.h>
#include "frenettraj.h"
struct Point2D
{
    double x, y;
    long index;
};

class Frenet2Cartesian
{
public:
    Frenet2Cartesian();
    //~Frenet2Cartesian();


    void transform(const FrenetTraj & frenetTraj);
    Point2D findNearestPoint(const Point2D& startPt);
    Point2D findArcEndPoint(const Point2D& startPt, double arcLength);
    void generateOnlineTraj(const Point2D& startPt, FrenetTraj& frenetTraj);

     void generateOfflineTraj();
private:
    QVector<Point2D> offlineTraj;
    QVector< QVector<Point2D> > onlineTrajs;
};

#endif // FRENET2CARTESIAN

