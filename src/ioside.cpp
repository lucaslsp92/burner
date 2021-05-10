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
    //// Dam break tank 2D with grating elements
    //// Variaveis (double, Point, etc..)
    bool grating = false;                       // grating flag
    double H = 0.3;                             // water column
    double W = 0.4;                             // water length
    double L = 1.6;                             // tank length
    double th = 0.75;                           // tank height
    double Hb = 0.01;                           // water height tank bottom
    double Wb = 1.6;                            // water length tank bottom
    double gs = 0.010;                          // grating spacing
    double gw = 0.004;                          // grating width
    double gx = 0.90;                           // grating position X
    double theta = 45*M_PI/180;                 // grating angle
    double gl = 0.12/cos(theta);                // grating length
    double gsy = (gs+gw)/cos(theta);            // grating spacing Y

    Point GE1 (gx,gs+0*gsy,0.0);
    Point GE2 (gx,gs+1*gsy,0.0);
    Point GE3 (gx,gs+2*gsy,0.0);
    Point GE4 (gx,gs+3*gsy,0.0);
    Point GE5 (gx,gs+4*gsy,0.0);
    Point GE6 (gx,gs+5*gsy,0.0);
    Point GE7 (gx,gs+6*gsy,0.0);
    Point GE8 (gx,gs+7*gsy,0.0);
    Point GE9 (gx,gs+8*gsy,0.0);
    Point GE10 (gx,gs+9*gsy,0.0);
    Point GE11 (gx,gs+10*gsy,0.0);
    Point GE12 (gx,gs+11*gsy,0.0);
    Point GE13 (gx,gs+12*gsy,0.0);
    Point GE14 (gx,gs+13*gsy,0.0);
    Point GE15 (gx,gs+14*gsy,0.0);
    Point GE16 (gx,gs+15*gsy,0.0);
    Point GE17 (gx,gs+16*gsy,0.0);
    Point GE18 (gx,gs+17*gsy,0.0);
    Point GE19 (gx,gs+18*gsy,0.0);
    Point GE20 (gx,gs+19*gsy,0.0);
    Point GE21 (gx,gs+20*gsy,0.0);
    Point GE22 (gx,gs+21*gsy,0.0);
    Point GE23 (gx,gs+22*gsy,0.0);
    Point GE24 (gx,gs+23*gsy,0.0);
    Point GE25 (gx,gs+24*gsy,0.0);
    Point GE26 (gx,gs+25*gsy,0.0);
    Point GE27 (gx,gs+26*gsy,0.0);
    Point GE28 (gx,gs+27*gsy,0.0);
    Point GE29 (gx,gs+28*gsy,0.0);
    Point GE30 (gx,gs+29*gsy,0.0);
    Point GE31 (gx,gs+30*gsy,0.0);
    Point GE32 (gx,gs+31*gsy,0.0);

    ///Regiões (Region)
    Region fluid = P.rectangleXY(W, H, DP);
    Region bottomFluid = P.rectangleXY(Wb, Hb, DP);
    Region tank  = P.rectangleXY(L, th, 10*DP);

    Region gratingElement1 = P.transformation(GE1,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement2 = P.transformation(GE2,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement3 = P.transformation(GE3,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement4 = P.transformation(GE4,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement5 = P.transformation(GE5,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement6 = P.transformation(GE6,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement7 = P.transformation(GE7,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement8 = P.transformation(GE8,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement9 = P.transformation(GE9,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement10 = P.transformation(GE10,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement11 = P.transformation(GE11,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement12 = P.transformation(GE12,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement13 = P.transformation(GE13,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement14 = P.transformation(GE14,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement15 = P.transformation(GE15,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement16 = P.transformation(GE16,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement17 = P.transformation(GE17,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement18 = P.transformation(GE18,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement19 = P.transformation(GE19,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement20 = P.transformation(GE20,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement21 = P.transformation(GE21,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement22 = P.transformation(GE22,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement23 = P.transformation(GE23,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement24 = P.transformation(GE24,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement25 = P.transformation(GE25,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement26 = P.transformation(GE26,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement27 = P.transformation(GE27,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement28 = P.transformation(GE28,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement29 = P.transformation(GE29,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement30 = P.transformation(GE30,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement31 = P.transformation(GE31,Z,0).rectangleXY(gl, gw, 10*DP);
    Region gratingElement32 = P.transformation(GE32,Z,0).rectangleXY(gl, gw, 10*DP);

    ///Operações
    if (tank)
    {
        if (fluid || bottomFluid) 
            return 0;

        if (grating)
        {
            if (gratingElement1 || gratingElement2 || gratingElement3 || gratingElement4 || gratingElement5 || 
            gratingElement6 || gratingElement7 || gratingElement8 || gratingElement9 || gratingElement10||
            gratingElement11|| gratingElement12|| gratingElement13|| gratingElement14|| gratingElement15||
            gratingElement16|| gratingElement17|| gratingElement18|| gratingElement19|| gratingElement20||
            gratingElement21|| gratingElement22|| gratingElement23|| gratingElement24|| gratingElement25||
            gratingElement26|| gratingElement27|| gratingElement28|| gratingElement29|| gratingElement30||
            gratingElement31|| gratingElement32)
            return 4;
        }

    	return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif