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

    /////////////////////////////////////////////
    /// Fixture drain recieving the flush of a WC
    /// Variaveis (double, Point, etc..)
    /// Tank wall dimensions (meters)
    double L = 10.0;        // pipe length
    double D = 0.078;       // pipe diameter
    double R = 1.5*D/2;     // elbow revolution radius
    double H = 0.01;        // inflow pipe length

    /// Reference points
    Point inflowPipeBeginn (0,R+D/2,0);
    Point elbowBeginn (R,R+D/2,0);
    Point pipeBeginn (R,D/2,0);

    /// Regiões (Region)
    /// Geometries and geometric transformations from file Point.h
    Region inflowPipe  = P.transformation(inflowPipeBeginn,Y).cylinder(D/2, H);
    Region elbow = P.transformation(elbowBeginn).toroid(R, D/2);
    Region pipe  = P.transformation(pipeBeginn,X).cylinder(D/2, L-R);

    ///Operações
    if (inflowPipe || (elbow && x<=R && y<=R+D/2) || pipe)
    {
    	return -1;                
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall

}

#endif