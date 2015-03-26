#ifndef FRENETTRAJ_H
#define FRENETTRAJ_H
struct Shift
{
    double quintic[6];
    double cost;
    double end_offset;
    double end_t;
    bool flag; //速度或加速度是否超过最大值的标记，是否可用
};

struct Traj
{
    Shift Lat;
    Shift Lon;
    double cost;
    double end_t;
    bool flag;
};

class FrenetTraj
{
public:
    FrenetTraj();
    ~FrenetTraj();
public:
    double *tt, *dd, *vv, *ss;            //保存离散后的时间，横向位移，速度，切向位移
    int tNum, dNum, vNum, sNum;//离散的数目
    double targetSpeed; //目标速度
    double target_d; //目标中线，此处需要修改，应该是目标点
    Shift *Clat;//横向所有轨迹，以时间为变量
    Shift *Clon;//侧向所有轨迹，以时间为变量
    Traj *traj;//
    //
    double MaxLatSpeed;
    double MaxLatAcc;
    double MaxLonSpeed;
    double MaxLonAcc;
public:
    //parameters，用于计算cost
    double para_lat_kt; //时间
    double para_lat_kj; //jerk
    double para_lat_kd; //和目标lateral距离
    double para_lat_k; // 整个lateral的权重

    double para_lon_kt;
    double para_lon_kj;
    double para_lon_ks;
    double para_lon_k;

    //Lateral movement
    void lat(double t0, double d0, double d0_1, double d0_2);
    void plotLat(double t0, const Shift clat[]);
    void plotBestLat(double t0, double ti, const double q[]);

    void lon(double t0, double s0, double s0_1, double s0_2);
    void plotLon(double t0, const Shift clon[]);

    //Longitudinal movement for Velocity Keeping
    void vel(double t0,double s0,double s0_1,double s0_2);
    void plotVel(double t0,const Shift clon[],int flag);
    void plotBestVel(double t0,double ti,const double q[],int flag);

    //Combine movement
    void combine(double t0,double d0,double d0_1,double d0_2,double s0,double s0_1,double s0_2);
    void plotCombine(double t0,const Traj tj[]);

private:
    double * linspace(double begin,double end,int num);
    //double poly_eval(const double c[], const int len, const double x);
};

#endif // FRENETTRAJ_H
