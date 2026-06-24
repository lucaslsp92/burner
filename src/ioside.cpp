#ifndef H_IOSIDE
#define H_IOSIDE

#include "Point.h"
#include "Particle.h"

CStlFile *stlProa = new CStlFile("data/CASCO-180-C01-ASCII.stl");

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

    ///Variaveis (double, Point, etc..)
    double tankLength = 1500.0;
    double tankHeight = 292.0;
    double waterHeight = 250.0;
    double beachLength = 200.0;

    double boxLenght = 311.0;
    double boxHeight = 33.5;
    double boxX = 544.5;
    double boxY = 236.5;

    Point base(1200.0, 0, 0);
	Point normal(0.7808688094, -0.6246950476, 0.0);
	Point wm(-20.0, 0.0, 0.0);
	Point boxCorner(boxX, boxY, -5.0*DP);

    ///Regiões (Region)
    Region tank = P.parallelepiped(tankLength, tankHeight, 10*DP);
    Region water = P.parallelepiped(tankLength, waterHeight, DP);
    Region wavemaker = P.transformation(wm).parallelepiped(20.0-5*DP, tankHeight, 10.0*DP);
    Region beach = P.plane(base, normal);
    Region proa = P.stlmodel(stlProa);
    Region box = P.transformation(boxCorner).parallelepiped(boxLenght, boxHeight, 10*DP);

    ///Operações
    if(proa)
        return 6;
    /*if(box)
        return 6;*/
    else if(water && beach)
        return 0;
    else if(tank && beach)
        return -1;

   return 2;
}
#endif