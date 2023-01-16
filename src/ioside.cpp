#ifndef H_IOSIDE
#define H_IOSIDE

#include "Point.h"
#include "Particle.h"
#include <algorithm>

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
    //// Grating modelin
    //// Variaveis (double, Point, etc..)
    double h0 = 0.240;   /// C4 h0 = 0.300, C3 h0 = 0.24, C2 h0 = 0.20, C1 h0 = 0.171
    double lb = 0.075;   /// lb = 0.040, lb = 0.075
    double l0 = 0.300;
    double h1 = 0.120;
    double l1 = 0.505;
    double h2 = 0.150;
    double l2 = 0.195;
    double h = 0.450;
    double l = 1.000;
    double b0 = 0.180;    
    double a0 = 0.150;
    double gratingWidth = 6*DP;
    //double gratingLength = 26*DP;               /// DP = 0.000234375
    //double gratingLength = 14*DP;               /// DP = 0.00046875
    //double gratingSpacing = 0.008855;           /// DP = 0.000234375
    double gratingSpacing = 0.008855;           /// DP = 0.00046875
    double gratingSpacingX = 0.00889;
    double firstGratingHeight = h2 + 0.008985 + gratingWidth/2;
    double firstGratingX = l0 + l1 - 0.008985 - gratingWidth/2;

    double boxHeight = 0.1;
    double boxWidth = 0.2;
    double firstGratingElementX = 0.06557077;
    //double firstGratingElementY = 0.01;
    double firstGratingElementY = 0.005;        /// DP = 0.00005859375
    //double gratingSpace = 0.00839961;         /// DP = 0.00046875
    //double gratingSpace = 0.007230004;          /// DP = 0.000234375
    //double gratingSpace = 0.006644004;          /// DP = 0.0001171875
    //double gratingSpace = 0.006270001;          /// DP = 0.00005859375
    double gratingSpace = 0.006144;               /// DP = 0.000029296875
    //double gratingLength = 14*DP;               /// DP = 0.00046875
    //double gratingLength = 27*DP;               /// DP = 0.000234375
    //double gratingLength = 52*DP;               /// DP = 0.0001171875
    //double gratingLength = 103*DP;               /// DP = 0.00005859375
    double gratingLength = 206*DP;               /// DP = 0.000029296875

    /// Points
    Point TNK(0,0,-5*DP);
    Point VD(l0+l1,0,-5*DP);
    Point GR0(l0+l1+gratingLength/2,firstGratingHeight+0*gratingSpacing,0);
    Point GR1(l0+l1+gratingLength/2,firstGratingHeight+1*gratingSpacing,0);
    Point GR2(l0+l1+gratingLength/2,firstGratingHeight+2*gratingSpacing,0);
    Point GR3(l0+l1+gratingLength/2,firstGratingHeight+3*gratingSpacing,0);
    Point GR4(l0+l1+gratingLength/2,firstGratingHeight+4*gratingSpacing,0);
    Point GR5(l0+l1+gratingLength/2,firstGratingHeight+5*gratingSpacing,0);
    Point GR6(l0+l1+gratingLength/2,firstGratingHeight+6*gratingSpacing,0);
    Point GR7(l0+l1+gratingLength/2,firstGratingHeight+7*gratingSpacing,0);
    Point GR8(l0+l1+gratingLength/2,firstGratingHeight+8*gratingSpacing,0);
    Point GR10(firstGratingX-0*gratingSpacingX,h2-gratingLength/2,0);
    Point GR11(firstGratingX-1*gratingSpacingX,h2-gratingLength/2,0);
    Point GR12(firstGratingX-2*gratingSpacingX,h2-gratingLength/2,0);
    Point GR13(firstGratingX-3*gratingSpacingX,h2-gratingLength/2,0);
    Point GR14(firstGratingX-4*gratingSpacingX,h2-gratingLength/2,0);
    Point GR15(firstGratingX-5*gratingSpacingX,h2-gratingLength/2,0);
    Point GR16(firstGratingX-6*gratingSpacingX,h2-gratingLength/2,0);
    Point GR17(firstGratingX-7*gratingSpacingX,h2-gratingLength/2,0);
    Point GR18(firstGratingX-8*gratingSpacingX,h2-gratingLength/2,0);
    Point VS(l0+l1+gratingWidth/2,h2+lb/2,0);
    Point HS(l0+l1-lb/2,h2-gratingWidth/2,0);

    Point BX(0,0,-5*DP);
    Point GR20(firstGratingElementX+0*gratingSpace,firstGratingElementY,0);
    Point GR21(firstGratingElementX+1*gratingSpace,firstGratingElementY,0);
    Point GR22(firstGratingElementX+2*gratingSpace,firstGratingElementY,0);
    Point GR23(firstGratingElementX+3*gratingSpace,firstGratingElementY,0);
    Point GR24(firstGratingElementX+4*gratingSpace,firstGratingElementY,0);
    Point GR25(firstGratingElementX+5*gratingSpace,firstGratingElementY,0);
    Point GR26(firstGratingElementX+6*gratingSpace,firstGratingElementY,0);
    Point GR27(firstGratingElementX+7*gratingSpace,firstGratingElementY,0);
    Point GR28(firstGratingElementX+8*gratingSpace,firstGratingElementY,0);
    
    ///Regiões (Region)
    Region tank = P.transformation(TNK).rectangleXY(l,h,10*DP);
    Region water0 = P.rectangleXY(l0,h0,DP);
    Region water1 = P.rectangleXY(l0+l1,h1,DP);
    Region ship = P.transformation(VD).rectangleXY(l2,h2,10*DP);
    Region grating0 = P.transformation(GR0).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating1 = P.transformation(GR1).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating2 = P.transformation(GR2).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating3 = P.transformation(GR3).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating4 = P.transformation(GR4).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating5 = P.transformation(GR5).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating6 = P.transformation(GR6).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating7 = P.transformation(GR7).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating8 = P.transformation(GR8).rectangleXYCenter(gratingLength,gratingWidth,10*DP);
    Region grating10 = P.transformation(GR10).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating11 = P.transformation(GR11).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating12 = P.transformation(GR12).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating13 = P.transformation(GR13).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating14 = P.transformation(GR14).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating15 = P.transformation(GR15).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating16 = P.transformation(GR16).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating17 = P.transformation(GR17).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating18 = P.transformation(GR18).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region verticalShield = P.transformation(VS).rectangleXYCenter(gratingWidth,lb,10*DP);
    Region horizontalShield = P.transformation(HS).rectangleXYCenter(lb,gratingWidth,10*DP);

    Region box = P.transformation(BX).rectangleXY(boxWidth,boxHeight,10*DP);
    Region grating20 = P.transformation(GR20).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating21 = P.transformation(GR21).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating22 = P.transformation(GR22).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating23 = P.transformation(GR23).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating24 = P.transformation(GR24).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating25 = P.transformation(GR25).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating26 = P.transformation(GR26).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating27 = P.transformation(GR27).rectangleXYCenter(gratingWidth,gratingLength,10*DP);
    Region grating28 = P.transformation(GR28).rectangleXYCenter(gratingWidth,gratingLength,10*DP);

    ///Operações
    /*if(tank && !(ship))
    {
        if (water0 || water1)
            return 0;
        else if (grating0 || grating1 || grating2 || grating3 || grating4 || grating5 || grating6 || grating7 || grating8)
            return 12;
        //else if (grating10 || grating11 || grating12 || grating13 || grating14 || grating15 || grating16 || grating17 || grating18)
            //return 12;
        //else if (verticalShield)
            //return 12;
        //else if (horizontalShield)
            //return 12;
        else
            return -1;
    }*/

    if(box)
    {
        if (grating20 || grating21 || grating22 || grating23 || grating24 || grating25 || grating26  || grating27  || grating28)
            return 6;
        else
            return -1;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif