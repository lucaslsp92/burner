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
    bool grating = true;                       // grating flag
    double H = 0.3;                             // water column
    double W = 0.4;                             // water length
    double L = 1.6;                             // tank length
    double th = 0.75;                           // tank height
    double Hb = 0.01;                           // water height tank bottom
    double Wb = 1.6;                            // water length tank bottom
    double gs = 0.010;                          // grating spacing
    double gb = 0.010;                          // space between tank bottom and first grating element
    double gw = 0.005;                          // grating width
    double gx = 0.90;                           // grating position X
    double glx = 0.12;                          // grating length X
    double theta = 45*M_PI/180;                 // grating angle
    double gl = glx/cos(theta);                 // grating length
    double gsy = (gs+gw)/cos(theta);            // grating spacing Y
    double ch = gw/cos(theta)-0.5*DP;           // chamfer length
    double gm;                                  // move granting downward Y
    if (DP == 0.001)
        gm = 0.013767767; 
    else if (DP == 0.0005)
        gm = 0.013767767;
    else if (DP == 0.0002)
        gm = 0.013767767;

    Point BF2 (gx+0.015,0.0,0.0);

    Point GE1 (gx,gb+0*gsy,0.0);
    Point GE2 (gx,gb+1*gsy,0.0);
    Point GE3 (gx,gb+2*gsy,0.0);
    Point GE4 (gx,gb+3*gsy,0.0);
    Point GE5 (gx,gb+4*gsy,0.0);
    Point GE6 (gx,gb+5*gsy,0.0);
    Point GE7 (gx,gb+6*gsy,0.0);
    Point GE8 (gx,gb+7*gsy,0.0);
    Point GE9 (gx,gb+8*gsy,0.0);
    Point GE10 (gx,gb+9*gsy,0.0);
    Point GE11 (gx,gb+10*gsy,0.0);
    Point GE12 (gx,gb+11*gsy,0.0);
    Point GE13 (gx,gb+12*gsy,0.0);
    Point GE14 (gx,gb+13*gsy,0.0);
    Point GE15 (gx,gb+14*gsy,0.0);
    Point GE16 (gx,gb+15*gsy,0.0);
    Point GE17 (gx,gb+16*gsy,0.0);
    Point GE18 (gx,gb+17*gsy,0.0);
    Point GE19 (gx,gb+18*gsy,0.0);
    Point GE20 (gx,gb+19*gsy,0.0);
    Point GE21 (gx,gb+20*gsy,0.0);
    Point GE22 (gx,gb+21*gsy,0.0);
    Point GE23 (gx,gb+22*gsy,0.0);
    Point GE24 (gx,gb+23*gsy,0.0);
    Point GE25 (gx,gb+24*gsy,0.0);
    Point GE26 (gx,gb+25*gsy,0.0);
    Point GE27 (gx,gb+26*gsy,0.0);
    Point GE28 (gx,gb+27*gsy,0.0);
    Point GE29 (gx,gb+28*gsy,0.0);
    Point GE30 (gx,gb+29*gsy,0.0);

    ///Regiões (Region)
    Region fluid = P.rectangleXY(W, H, DP);
    Region bottomFluid1 = P.rectangleXY(gx-0.002, Hb, DP);
    Region bottomFluid2 = P.transformation(BF2,Z,0).rectangleXY(Wb-gx-0.015, Hb, DP);
    Region tank  = P.rectangleXY(L, th, 10*DP);
    Region chamfer = P.transformation(GE1,Z,theta*180/M_PI).rectangleXYCenter(ch, ch, 10*DP);

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

    ///Operações
    if (tank)
    {
        if (fluid || bottomFluid1 || bottomFluid2) 
            return 0;

        if (grating)
        {
            if ((gratingElement1 && !(chamfer)) || gratingElement2 || gratingElement3 || gratingElement4 || gratingElement5 || 
            gratingElement6 || gratingElement7 || gratingElement8 || gratingElement9 || gratingElement10||
            gratingElement11|| gratingElement12|| gratingElement13|| gratingElement14|| gratingElement15||
            gratingElement16|| gratingElement17|| gratingElement18|| gratingElement19|| gratingElement20||
            gratingElement21|| gratingElement22|| gratingElement23|| gratingElement24|| gratingElement25||
            gratingElement26|| gratingElement27|| gratingElement28|| gratingElement29|| gratingElement30)
            return 4;
        }

    	return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif