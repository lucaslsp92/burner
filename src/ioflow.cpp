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
    //// Dam break tank 2D with grating elements
    bool grating = false;                       // grating flag
    double gs = 0.010;                          // grating spacing
    double gw = 0.004;                          // grating width
    double gx = 0.90;                           // grating position X
    double theta = 45*M_PI/180;                 // grating angle
    double gl = 0.12/cos(theta);                // grating length
    double gsy = (gs+gw)/cos(theta);            // grating spacing Y

    Point IGE1 (gx+DP,gs+0*gsy+DP,0.0);
    Point IGE2 (gx+DP,gs+1*gsy+DP,0.0);
    Point IGE3 (gx+DP,gs+2*gsy+DP,0.0);
    Point IGE4 (gx+DP,gs+3*gsy+DP,0.0);
    Point IGE5 (gx+DP,gs+4*gsy+DP,0.0);
    Point IGE6 (gx+DP,gs+5*gsy+DP,0.0);
    Point IGE7 (gx+DP,gs+6*gsy+DP,0.0);
    Point IGE8 (gx+DP,gs+7*gsy+DP,0.0);
    Point IGE9 (gx+DP,gs+8*gsy+DP,0.0);
    Point IGE10 (gx+DP,gs+9*gsy+DP,0.0);
    Point IGE11 (gx+DP,gs+10*gsy+DP,0.0);
    Point IGE12 (gx+DP,gs+11*gsy+DP,0.0);
    Point IGE13 (gx+DP,gs+12*gsy+DP,0.0);
    Point IGE14 (gx+DP,gs+13*gsy+DP,0.0);
    Point IGE15 (gx+DP,gs+14*gsy+DP,0.0);
    Point IGE16 (gx+DP,gs+15*gsy+DP,0.0);
    Point IGE17 (gx+DP,gs+16*gsy+DP,0.0);
    Point IGE18 (gx+DP,gs+17*gsy+DP,0.0);
    Point IGE19 (gx+DP,gs+18*gsy+DP,0.0);
    Point IGE20 (gx+DP,gs+19*gsy+DP,0.0);
    Point IGE21 (gx+DP,gs+20*gsy+DP,0.0);
    Point IGE22 (gx+DP,gs+21*gsy+DP,0.0);
    Point IGE23 (gx+DP,gs+22*gsy+DP,0.0);
    Point IGE24 (gx+DP,gs+23*gsy+DP,0.0);
    Point IGE25 (gx+DP,gs+24*gsy+DP,0.0);
    Point IGE26 (gx+DP,gs+25*gsy+DP,0.0);
    Point IGE27 (gx+DP,gs+26*gsy+DP,0.0);
    Point IGE28 (gx+DP,gs+27*gsy+DP,0.0);
    Point IGE29 (gx+DP,gs+28*gsy+DP,0.0);
    Point IGE30 (gx+DP,gs+29*gsy+DP,0.0);
    Point IGE31 (gx+DP,gs+30*gsy+DP,0.0);
    Point IGE32 (gx+DP,gs+31*gsy+DP,0.0);

    Point VGE1 (gx+3*DP,gs+0*gsy+3*DP,0.0);
    Point VGE2 (gx+3*DP,gs+1*gsy+3*DP,0.0);
    Point VGE3 (gx+3*DP,gs+2*gsy+3*DP,0.0);
    Point VGE4 (gx+3*DP,gs+3*gsy+3*DP,0.0);
    Point VGE5 (gx+3*DP,gs+4*gsy+3*DP,0.0);
    Point VGE6 (gx+3*DP,gs+5*gsy+3*DP,0.0);
    Point VGE7 (gx+3*DP,gs+6*gsy+3*DP,0.0);
    Point VGE8 (gx+3*DP,gs+7*gsy+3*DP,0.0);
    Point VGE9 (gx+3*DP,gs+8*gsy+3*DP,0.0);
    Point VGE10 (gx+3*DP,gs+9*gsy+3*DP,0.0);
    Point VGE11 (gx+3*DP,gs+10*gsy+3*DP,0.0);
    Point VGE12 (gx+3*DP,gs+11*gsy+3*DP,0.0);
    Point VGE13 (gx+3*DP,gs+12*gsy+3*DP,0.0);
    Point VGE14 (gx+3*DP,gs+13*gsy+3*DP,0.0);
    Point VGE15 (gx+3*DP,gs+14*gsy+3*DP,0.0);
    Point VGE16 (gx+3*DP,gs+15*gsy+3*DP,0.0);
    Point VGE17 (gx+3*DP,gs+16*gsy+3*DP,0.0);
    Point VGE18 (gx+3*DP,gs+17*gsy+3*DP,0.0);
    Point VGE19 (gx+3*DP,gs+18*gsy+3*DP,0.0);
    Point VGE20 (gx+3*DP,gs+19*gsy+3*DP,0.0);
    Point VGE21 (gx+3*DP,gs+20*gsy+3*DP,0.0);
    Point VGE22 (gx+3*DP,gs+21*gsy+3*DP,0.0);
    Point VGE23 (gx+3*DP,gs+22*gsy+3*DP,0.0);
    Point VGE24 (gx+3*DP,gs+23*gsy+3*DP,0.0);
    Point VGE25 (gx+3*DP,gs+24*gsy+3*DP,0.0);
    Point VGE26 (gx+3*DP,gs+25*gsy+3*DP,0.0);
    Point VGE27 (gx+3*DP,gs+26*gsy+3*DP,0.0);
    Point VGE28 (gx+3*DP,gs+27*gsy+3*DP,0.0);
    Point VGE29 (gx+3*DP,gs+28*gsy+3*DP,0.0);
    Point VGE30 (gx+3*DP,gs+29*gsy+3*DP,0.0);
    Point VGE31 (gx+3*DP,gs+30*gsy+3*DP,0.0);
    Point VGE32 (gx+3*DP,gs+31*gsy+3*DP,0.0);

    Region innerGratingElement1 = P.transformation(IGE1,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement2 = P.transformation(IGE2,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement3 = P.transformation(IGE3,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement4 = P.transformation(IGE4,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement5 = P.transformation(IGE5,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement6 = P.transformation(IGE6,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement7 = P.transformation(IGE7,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement8 = P.transformation(IGE8,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement9 = P.transformation(IGE9,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement10 = P.transformation(IGE10,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement11 = P.transformation(IGE11,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement12 = P.transformation(IGE12,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement13 = P.transformation(IGE13,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement14 = P.transformation(IGE14,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement15 = P.transformation(IGE15,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement16 = P.transformation(IGE16,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement17 = P.transformation(IGE17,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement18 = P.transformation(IGE18,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement19 = P.transformation(IGE19,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement20 = P.transformation(IGE20,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement21 = P.transformation(IGE21,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement22 = P.transformation(IGE22,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement23 = P.transformation(IGE23,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement24 = P.transformation(IGE24,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement25 = P.transformation(IGE25,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement26 = P.transformation(IGE26,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement27 = P.transformation(IGE27,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement28 = P.transformation(IGE28,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement29 = P.transformation(IGE29,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement30 = P.transformation(IGE30,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement31 = P.transformation(IGE31,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);
    Region innerGratingElement32 = P.transformation(IGE32,Z,0).rectangleXY(gl-2*DP, gw-2*DP, 10*DP);

    Region voidGratingElement1 = P.transformation(VGE1,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement2 = P.transformation(VGE2,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement3 = P.transformation(VGE3,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement4 = P.transformation(VGE4,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement5 = P.transformation(VGE5,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement6 = P.transformation(VGE6,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement7 = P.transformation(VGE7,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement8 = P.transformation(VGE8,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement9 = P.transformation(VGE9,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement10 = P.transformation(VGE10,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement11 = P.transformation(VGE11,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement12 = P.transformation(VGE12,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement13 = P.transformation(VGE13,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement14 = P.transformation(VGE14,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement15 = P.transformation(VGE15,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement16 = P.transformation(VGE16,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement17 = P.transformation(VGE17,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement18 = P.transformation(VGE18,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement19 = P.transformation(VGE19,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement20 = P.transformation(VGE20,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement21 = P.transformation(VGE21,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement22 = P.transformation(VGE22,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement23 = P.transformation(VGE23,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement24 = P.transformation(VGE24,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement25 = P.transformation(VGE25,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement26 = P.transformation(VGE26,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement27 = P.transformation(VGE27,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement28 = P.transformation(VGE28,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement29 = P.transformation(VGE29,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement30 = P.transformation(VGE30,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement31 = P.transformation(VGE31,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);
    Region voidGratingElement32 = P.transformation(VGE32,Z,0).rectangleXY(gl-6*DP, gw-6*DP, 10*DP);

    if (z < 0.0 || z >= DP)
    {
        return -1;
    }

    if (grating)
    {
        if (voidGratingElement1 || voidGratingElement2 || voidGratingElement3 || voidGratingElement4 || voidGratingElement5 ||
        voidGratingElement6 || voidGratingElement7 || voidGratingElement8 || voidGratingElement9 || voidGratingElement10||
        voidGratingElement11|| voidGratingElement12|| voidGratingElement13|| voidGratingElement14|| voidGratingElement15||
        voidGratingElement16|| voidGratingElement17|| voidGratingElement18|| voidGratingElement19|| voidGratingElement20||
        voidGratingElement21|| voidGratingElement22|| voidGratingElement23|| voidGratingElement24|| voidGratingElement25||
        voidGratingElement26|| voidGratingElement27|| voidGratingElement28|| voidGratingElement29|| voidGratingElement30||
        voidGratingElement31|| voidGratingElement32)
        {
            return -1;
        }

        else if (innerGratingElement1 || innerGratingElement2 || innerGratingElement3 || innerGratingElement4 || innerGratingElement5 ||
            innerGratingElement6 || innerGratingElement7 || innerGratingElement8 || innerGratingElement9 || innerGratingElement10||
            innerGratingElement11|| innerGratingElement12|| innerGratingElement13|| innerGratingElement14|| innerGratingElement15||
            innerGratingElement16|| innerGratingElement17|| innerGratingElement18|| innerGratingElement19|| innerGratingElement20||
            innerGratingElement21|| innerGratingElement22|| innerGratingElement23|| innerGratingElement24|| innerGratingElement25||
            innerGratingElement26|| innerGratingElement27|| innerGratingElement28|| innerGratingElement29|| innerGratingElement30||
            innerGratingElement31|| innerGratingElement32)
        {
            return id+1;
        }
    }
    

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif