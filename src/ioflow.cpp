#ifndef H_IOFLOW
#define H_IOFLOW

#include "Point.h"
#include "Particle.h"

int Particle::ioflow (double DP, double geometry[])
{
    Point P (r, DP);

    double x = P.x;
    double y = P.y;
    double z = P.z;

    Point O;
    Point X (1,0,0);
    Point Y (0,1,0);
    Point Z (0,0,1);

    //////////////////////////////////////////
    //// Grating modeling
    //// Variaveis (double, Point, etc..)
    double h0 = 0.240;   /// C4 h0 = 0.300, C3 h0 = 0.24, C2 h0 = 0.20, C1 h0 = 0.171
    double l0 = 0.300;   
    double h1 = 0.120;
    double l1 = 0.505;
    double h2 = 0.150;
    double l2 = 0.195;
    double h = 0.450;
    double l = 1.000;
    double b0 = 0.180;
    double a0 = 0.150;

    /// Points

    ///Regiões (Region)

    ///Operações

    ///Return padrão (caso nenhuma operação seja feita)
    /*if (y > h)
    {
        return -1;
    }

    if ((id == 2 || id == 3) && (x >= l0+l1 && x <= l0+l1+3*DP) && (y > 0) && (z > 0 && z <= DP))
    {
        return id+2;
    }

    if ((id == 2 || id == 3) && (x >= l0+l1+(l2-b0) && x <= l0+l1+l2) && (y >= h2-3*DP && y <= h2) && (z > 0 && z <= DP))
    {
        return id+4;
    }

    if ((id == 2 || id == 3) && (x >= l0+l1+l2 && x <= l0+l1+l2+3*DP) && (y >= h2-3*DP && y <= h2+a0) && (z > 0 && z <= DP))
    {
        return id+6;
    }

    if ((id == 2 || id == 3) && (x >= l0+l1+3*DP && x < l0+l1+(l2-b0)) && (y >= h2-3*DP && y <= h2) && (z > 0 && z <= DP))
    {
        return id+8;
    }*/

    /*if (x < 0 || x > 0.2)
    {
        return -1;
    }*/

    if(x > 0.075 && x < 0.125 && y > 0.1 && z > 0.0 && z < DP)
    {
        return id + 2;
    }

    if((x < 0 || x > 0.2 || y > 0) && (id == 2 || id == 3) && z > 0.0 && z < DP)
    {
        return id + 8;
    }

    if (z < 0 || z >= DP)
    {
        return -1;
    }

    return id;
}
#endif