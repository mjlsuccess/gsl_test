#include "frenettraj.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_heapsort.h>

#include <qopengl.h>

#define HAVE_INLINE

FrenetTraj::FrenetTraj()
{
    tNum=3;
    dNum=5;//21
    sNum=11;//21
    vNum=5;//21
    tt=linspace(2,4,tNum);
    dd=linspace(-4,8,dNum);
    ss=linspace(-0.6,0.6,sNum);
    vv=linspace(0,34,vNum);

    Clat=new Shift[tNum*dNum];
    Clon=new Shift[tNum*vNum];
    traj=new Traj[tNum*dNum*vNum];

    targetSpeed = 12;
    target_d=0;

    para_lat_kj=1;
    para_lat_kt=1;
    para_lat_kd=20;
    para_lat_k=1;

    para_lon_kj=1;
    para_lon_kt=1;
    para_lon_ks=1;
    para_lon_k=50;

    //
    MaxLonSpeed=40;
    MaxLonAcc=10;
    MaxLatSpeed=10;
    MaxLatAcc=10;
}

FrenetTraj::~FrenetTraj()
{
    delete []tt;
    delete []dd;
    delete []vv;
    delete []ss;
    delete []Clat;
    delete []Clon;
    delete []traj;
}
//求解横向位移d(t)，起点状态(d0, d0_1, d0_2)，即起点处的位置，速度和加速度
// 终点状态(T, dt, dt_1, dt_2)，即2终点处的位置，速度和加速度，这些可以构成6个方程组，来求解5次方程
// 在终点时我们希望轨迹能平行与中线，所以在终点处横向速度为0，加速度为0
//最后我们通过离散时间T和横向的位移来获得轨迹集
void FrenetTraj::lat(double t0,double d0,double d0_1,double d0_2)
{
    int i,k;
    int h=0;
    double T,di;
    gsl_vector *x = gsl_vector_alloc (6);
    gsl_permutation * p = gsl_permutation_alloc (6);
    for (k=0;k<tNum;k++)    //离散时间
    {
        T=tt[k];
        for (i=0;i<dNum;i++)  //离散横向位移
        {
            di=dd[i];
            double A[]=             //
            {1,		t0,				gsl_pow_int(t0,2),		gsl_pow_int(t0,3),			gsl_pow_int(t0,4),				gsl_pow_int(t0,5),
            0,		1,				2*t0,							3*gsl_pow_int(t0,2),		4*gsl_pow_int(t0,3),			5*gsl_pow_int(t0,4),
            0,		0,				2,								6*t0,								12*gsl_pow_int(t0,2),			20*gsl_pow_int(t0,3),
            1,		T+t0,		gsl_pow_int(T+t0,2),	gsl_pow_int(T+t0,3),		gsl_pow_int(T+t0,4),			gsl_pow_int(T+t0,5),
            0,		1,				2*(T+t0),					3*gsl_pow_int(T+t0,2),	4*gsl_pow_int(T+t0,3),		5*gsl_pow_int(T+t0,4),
            0,		0,				2,								6*(T+t0),						12*gsl_pow_int(T+t0,2),		20*gsl_pow_int(T+t0,3)};
            double B[]={d0,		d0_1,		d0_2,		di,		0,		0};
            gsl_matrix_view m = gsl_matrix_view_array (A, 6, 6);
            gsl_vector_view b = gsl_vector_view_array (B, 6);

            //求解线性方程 AX = B
            int s;
            gsl_linalg_LU_decomp (&m.matrix, p, &s);
            gsl_linalg_LU_solve(&m.matrix,p,&b.vector,x);
            double c[6];
            for (int w=0;w<6;w++)
            {//保存求解的5次方程系数
                c[w]=gsl_vector_get(x,w);
                Clat[h].quintic[w]=c[w];
            }
            //去掉速度和加速度超过阈值的轨迹
            double res[3];
            bool flag=1;
            for (double j=t0;j<t0+T;j+=0.1)
            {
                gsl_poly_eval_derivs(c,6,j,res,3);
                if ((abs(res[1])>MaxLatSpeed)||(abs(res[2])>MaxLatAcc))
                {
                    flag=0;
                    break;
                }
              }
            if (flag==0)
            {
                Clat[h].flag=0;
            }
            //计算cost，jerk 是先计算quintic的三阶导数，然后平方，再积分得到Jd系数
            double Jd[]={0,		36*c[3]*c[3],		144*c[3]*c[4],		192*c[4]*c[4]+240*c[3]*c[5],	720*c[4]*c[5],		720*c[5]*c[5]};
            double cost=para_lat_kj*(gsl_poly_eval(Jd,6,T+t0)-gsl_poly_eval(Jd,6,t0))+para_lat_kt*T+para_lat_kd*(di-target_d)*(di-target_d); // jerk + Time + diff(lateral)
            Clat[h].cost=cost;
            Clat[h].end_t=T;
            Clat[h].end_offset=di;
            h++;
        }
    }
    gsl_permutation_free (p);
    gsl_vector_free (x);
}

void FrenetTraj::plotLat(double t0,const Shift clat[])
{
    glPushMatrix();
    glTranslatef(-1.5,-0.5,0);
    glScalef(0.5,0.5,0.5);
    glLineWidth(2.0);
    for (int i=1;i<tNum*dNum;i++)
    {
        glBegin(GL_LINE_STRIP);
        for(double t=t0;t<clat[i].end_t+t0;t+=0.01)
            glVertex3f(t,gsl_poly_eval(clat[i].quintic,6,t),0.0);
        glEnd();
    }
    glPopMatrix();
}
void FrenetTraj::plotBestLat(double t0,double ti,const double q[])
{
    glPushMatrix();
    glTranslatef(-1.5,-0.5,0);
    glScalef(0.5,0.5,0.5);
    glBegin(GL_LINE_STRIP);
    for(double t=t0;t<ti+t0;t+=0.01)
        glVertex3f(t,gsl_poly_eval(q,6,t),0.0);
    glEnd();
    glPopMatrix();
}

//和函数vel的功能相同
void FrenetTraj::lon(double t0, double s0, double s0_1, double s0_2)
{
    int i,k;
    int h=0;
    double T,vi;
    gsl_vector *x = gsl_vector_alloc (6);
    gsl_permutation * p = gsl_permutation_alloc (6);
    for (k=0; k<tNum; k++)
    {
        T = tt[k];
        for(i=0; i<vNum; i++)
        {
            double A[] =
            {1,		t0,				gsl_pow_int(t0,2),		gsl_pow_int(t0,3),			gsl_pow_int(t0,4),				gsl_pow_int(t0,5),
            0,		1,				2*t0,							3*gsl_pow_int(t0,2),		4*gsl_pow_int(t0,3),			5*gsl_pow_int(t0,4),
            0,		0,				2,								6*t0,								12*gsl_pow_int(t0,2),			20*gsl_pow_int(t0,3),
            0,		1,				2*(T+t0),					3*gsl_pow_int(T+t0,2),	4*gsl_pow_int(T+t0,3),		5*gsl_pow_int(T+t0,4),
            0,		1,				2*(T+t0+0.2),					3*gsl_pow_int(T+t0+0.2,2),	4*gsl_pow_int(T+t0+0.2,3),		5*gsl_pow_int(T+t0+0.2,4),
            0,		0,				2,								6*(T+t0),						12*gsl_pow_int(T+t0,2),		20*gsl_pow_int(T+t0,3)};

            vi=vv[i];
            double B[]={s0,		s0_1,		s0_2,		vi,		vi,  0}; //认为在终点附近速度已经保持不变了
            gsl_matrix_view m = gsl_matrix_view_array (A, 6, 6);
            gsl_vector_view b = gsl_vector_view_array (B, 6);
            //求解线性方程 AX = B
            int s;
            gsl_linalg_LU_decomp (&m.matrix, p, &s);
            gsl_linalg_LU_solve(&m.matrix,p,&b.vector,x);
            double c[6];
            for (int w=0;w<6;w++)
            {//保存求解的5次方程系数
                c[w]=gsl_vector_get(x,w);
                Clon[h].quintic[w]=c[w];
            }

            //去掉速度和加速度超过阈值的轨迹
            double res[3];
            bool flag=1;
            for (double j=t0;j<t0+T;j+=0.1)
            {
                gsl_poly_eval_derivs(c,6,j,res,3);
                if ((abs(res[1])>MaxLonSpeed)||(abs(res[2])>MaxLonAcc))
                {
                    flag=0;
                    break;
                }
            }
            if (flag==0)
            {
                Clon[h].flag=0;
            }
            //
            double Js[]={0,		36*c[3]*c[3],		144*c[3]*c[4],		192*c[4]*c[4]+240*c[3]*c[5],	720*c[4]*c[5],		720*c[5]*c[5]};
            double cost=para_lon_kj*(gsl_poly_eval(Js,6,T+t0)-gsl_poly_eval(Js,6,t0))+para_lon_kt*T+para_lon_ks*(targetSpeed-vi)*(targetSpeed-vi);
            Clon[h].cost=cost;
            Clon[h].end_t=T;
            Clon[h].end_offset=vi;
            h++;
        }
    }
    gsl_permutation_free (p);
    gsl_vector_free (x);
}

void FrenetTraj::plotLon(double t0, const Shift clon[])
{
    glPushMatrix();
    glTranslatef(-15,-2.5,0);
    glScalef(0.25,0.25,0.25);
    glLineWidth(2.0);
    for (int i=1;i<tNum*vNum;i++)
    {
        glBegin(GL_LINE_STRIP);
        for(double t=t0;t<clon[i].end_t+t0;t+=0.01)
            glVertex3f(t,gsl_poly_eval(clon[i].quintic,6,t),0.0);
        glEnd();
    }
    glPopMatrix();
}

//dt的计算中因为车必须在车道内行驶，所以可以直接离散终点的位置
//s(t)是无限延伸的，无法对s(t)离散化，如果想得到一组与纵向位移无关的轨迹，一个合理的方法是对终点处的速度离散话
//另外，在终点附近车以恒定速度行驶是我们可接受的理想状态，所以起点状态(s0, s0_1, s0_2),终点(T, st_1, 0)
//
void FrenetTraj::vel(double t0,double s0,double s0_1,double s0_2)
{
    int i,k;
    int h=0;
    double T,vi;
    gsl_vector *x = gsl_vector_alloc (5);
    gsl_permutation * p = gsl_permutation_alloc (5);
    for (k=0;k<tNum;k++)//离散时间
    {
        T=tt[k];
        for (i=0;i<vNum;i++)//离散终点速度
        {//
            double A[]=
            {1,		t0,				gsl_pow_int(t0,2),		gsl_pow_int(t0,3),			gsl_pow_int(t0,4),
            0,		1,				2*t0,							3*gsl_pow_int(t0,2),		4*gsl_pow_int(t0,3),
            1,		T+t0,		gsl_pow_int(T+t0,2),	gsl_pow_int(T+t0,3),		gsl_pow_int(T+t0,4),
            0,		1,				2*(T+t0),					3*gsl_pow_int(T+t0,2),	4*gsl_pow_int(T+t0,3),
            1,		T+t0+0.2,		gsl_pow_int(T+t0+0.2,2),	gsl_pow_int(T+t0+0.2,3),		gsl_pow_int(T+t0+0.2,4)}; //认为在终点附近速度已经保持不变了
            gsl_matrix_view m = gsl_matrix_view_array (A, 5, 5);
            vi=vv[i];
            double B[]={s0_1,		s0_2,		vi,		0,		vi}; //认为在终点附近速度已经保持不变了
            gsl_vector_view b = gsl_vector_view_array (B, 5);

            int s;
            gsl_linalg_LU_decomp (&m.matrix, p, &s);
            gsl_linalg_LU_solve(&m.matrix,p,&b.vector,x);
            //此时得到的解为速度的描述
            double c[6];
            c[0]=0;
            //积分得到纵向位移的描述
            for (int w=0;w<5;w++)
            {
                c[w+1]=gsl_vector_get(x,w)/(w+1);
                Clon[h].quintic[w+1]=c[w+1];
            }
            //计算s0，特别注意！！！！
            //之前把gsl_poly_eval(c,6,t0)写成了gsl_poly_eval(c,5,t0)，找错误找了一天！！！
            double s0_=s0-gsl_poly_eval(c,6,t0);
            c[0]=s0_;
            Clon[h].quintic[0]=c[0];

            //去掉速度和加速度超过阈值的轨迹
            double res[3];
            bool flag=1;
            for (double j=t0;j<t0+T;j+=0.1)
            {
                gsl_poly_eval_derivs(c,6,j,res,3);
                if ((abs(res[1])>MaxLonSpeed)||(abs(res[2])>MaxLonAcc))
                {
                    flag=0;
                    break;
                }
            }
            if (flag==0)
            {
                Clon[h].flag=0;
            }
            //
            double Js[]={0,		36*c[3]*c[3],		144*c[3]*c[4],		192*c[4]*c[4]+240*c[3]*c[5],	720*c[4]*c[5],		720*c[5]*c[5]};
            double cost=para_lon_kj*(gsl_poly_eval(Js,6,T+t0)-gsl_poly_eval(Js,6,t0))+para_lon_kt*T+para_lon_ks*(targetSpeed-vi)*(targetSpeed-vi);
            Clon[h].cost=cost;
            Clon[h].end_t=T;
            Clon[h].end_offset=vi;
            h++;
        }
    }
    gsl_permutation_free (p);
    gsl_vector_free (x);
}
void FrenetTraj::plotVel(double t0,const Shift clon[],int flag)
{
    glPushMatrix();
    if(flag==0) //画轨迹
    {
        glTranslatef(-1.5,-2.3,0);
        glScalef(0.25,0.25,0.5);

    }
    else if (flag==1) //画速度
    {
        glTranslatef(-1.5,-2.0,0);
        glScalef(0.5,0.5,0.5);
    }
    glLineWidth(2.0);
    double y[2];
    for (int i=1;i<tNum*vNum;i++)
    {
        glBegin(GL_LINE_STRIP);
        for(double t=t0;t<clon[i].end_t+t0;t+=0.01)
        {
            gsl_poly_eval_derivs(clon[i].quintic,6,t,y,2);
            glVertex3f(t,y[flag],0.0);
        }
        glEnd();
    }
    glPopMatrix();
}


void FrenetTraj::plotBestVel(double t0,double ti,const double q[],int flag)
{
    glPushMatrix();
    if(flag==0)
    {
        glTranslatef(-1.5,-2.3,0);
        glScalef(0.25,0.25,0.5);
    }
    else if (flag==1)
    {
        glTranslatef(-1.5,-2.0,0);
        glScalef(0.5,0.5,0.5);
    }
    double y[2];
    glBegin(GL_LINE_STRIP);
    for(double t=t0;t<ti+t0;t+=0.01)
    {
        gsl_poly_eval_derivs(q,6,t,y,2);
        glVertex3f(t,y[flag],0.0);
    }
    glEnd();
    glPopMatrix();
}

void FrenetTraj::combine(double t0,double d0,double d0_1,double d0_2,double s0,double s0_1,double s0_2)
{
    int i,j,k;
    lat(t0,d0,d0_1,d0_2);
    vel(t0,s0,s0_1,s0_2);
    int h=0;
    for (k=0;k<tNum;k++)
    {
        for (i=0;i<dNum;i++)
        {
            for (j=0;j<vNum;j++)
            {
                traj[h].flag=1;
                if ((Clat[k*dNum+i].flag==0)||(Clon[k*vNum+j].flag==0))
                {
                    traj[h].flag=0;
                }
                traj[h].Lat=Clat[k*dNum+i];
                traj[h].Lon=Clon[k*vNum+j];
                traj[h].cost=para_lat_k*traj[h].Lat.cost+para_lon_k*traj[h].Lon.cost;
                traj[h].end_t=tt[k];
                h++;
            }
        }
    }
}

void FrenetTraj::plotCombine(double t0,const Traj tj[])
{
    double r=200;
    glPushMatrix();
    glTranslatef(40.0,-2.2,0);
    glScalef(0.2,0.2,0.2);
    glLineWidth(1.0);
    double x,y;
    for (int i=1;i<tNum*dNum*vNum;i++)
    {
        glBegin(GL_LINE_STRIP);
        for(double t=t0;t<tj[i].end_t+t0;t+=0.01)
        {
            x=cos(M_PI-gsl_poly_eval(tj[i].Lon.quintic,6,t)/r)*(r+gsl_poly_eval(tj[i].Lat.quintic,6,t));
            y=sin(M_PI-gsl_poly_eval(tj[i].Lon.quintic,6,t)/r)*(r+gsl_poly_eval(tj[i].Lat.quintic,6,t));
            glVertex3f(x,y,0.0);
        }
        glEnd();
    }
    glPopMatrix();
}


double *  FrenetTraj::linspace(double begin,double end,int num)
{
    if (num==0)
    {
        return 0;
    }
    double *temp=new double[num];
    double gap;
    gap=(end-begin)/(num-1);
    for (int i=0;i<num;i++)
    {
        temp[i]=begin+i*gap;
    }
    return temp;
}
