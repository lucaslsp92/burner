#ifndef H_IOSIDE
#define H_IOSIDE

#include "Point.h"
#include "Particle.h"

int Particle::ioside (double DP, double geometry[])
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
    //// Poiseuille flow 2D
    //// Variaveis (double, Point, etc..)
    double R = 0.1;
    double r = 0.05;
    double H = 0.3;
    double L = 0.35;

    Point A(0.0, 0.0, 0.0);

    ///Regi�es (Region)
    Region cone = P.transformation(A, Y).cone(R, r, H);
    Region table = P.transformation(A, Y).cylinder(L, 2.0*H);

    ///Opera��es
    if(table)
    {
        if (cone && (z>0.0 && z<=DP)) 
            return 0;
        else
            return -1;
    }    

    ///Return padr�o (constr�i a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif