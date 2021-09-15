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
    //// Moonpool sloshing
    //// Variaveis (double, Point, etc..)
    int WBid = 4;                   // id wave breaker  
    double L = 38.4;                // tank length
    double H = 19.7;              // tank height
    double W = 2.0;               // tank width
    double SL = 12.8;               // step length
    double SH = 7.3;                // step height
    double SC = 2.912043956;        // step chamfer length
    double SA = 15.945395774;       // stpe chamfer angle
    double WL = 0.8;                // wave breaker length
    //double WH = 0.25+0.8;           // wave breaker height
    double WH = 6*DP;               // wave breaker height
    double WR = 0.8;                // wave breaker radius
    double WF = 0.5;                // wave breaker front
    double WC = 0.360555128;        // wave breaker chamfer length
    double alfa = 56.309932333;     // wave breaker chamfer angle left
    double WBP = 10.5;              // wave breaker bottom position
    double WTP = 13.7;              // wave breaker top position
    double FH = 11.5;               // fluid height
    double DFL = 1.13137085;        // deflector length
    double DFA = 45;                // deflector angle
    double WBHD = 0.2;              // wave breaker hole diameter
    double WBHP = 0.4;              // wave breaker hole x and z position

    /// Points
    Point SCP(SL-0.8,0,0);                  // step chamfer position

    Point WBL (0,WBP-WH/2,0);               // wave breaker position bottom left
    Point WBLR (WR,WBP-WH/2,-W/2);             // wave breaker radius bottom left
    Point WBLC (0.5,WBP,-W/2);                 // wave breaker chamfer bottom left
    Point WBLF (WL-6*DP,WBP-WF/2,0);        // wave breaker front bottom left

    Point WTL (0,WTP-WH/2,0);                   // wave breaker position top left
    Point WTLR (WR,WTP-WH/2,-W/2);             // wave breaker radius top left
    Point WTLC (0.5,WTP,-W/2);                 // wave breaker chamfer top left
    Point WTLF (WL-6*DP,WTP-WF/2,0);            // wave breaker front top left
    Point WTLH (WBHP, WTP-WH/2, WBHP);           // wave breaker hole top left

    Point DFP (0,WTP-WL-3*DP,0);            // deflector position 

    ///Regiões (Region)
    Region tank = P.rectangleXY(L, H, W);
    Region fluid = P.rectangleXY(L, FH, W);

    Region step = P.rectangleXY(SL, SH, W);
    Region stepChamfer = P.transformation(SCP, Z, SA).rectangleXY(SC, SC, W);

    Region waveBreakerBL = P.transformation(WBL).rectangleXY(WL, WH, W);
    Region waveBreakerBLR = P.transformation(WBLR).cylinder(WR, W);
    Region waveBreakerBLC = P.transformation(WBLC, Z, alfa).rectangleXY(WC, WC, W);
    Region waveBreakerBLF = P.transformation(WBLF).rectangleXY(6*DP, WF, W);

    Region waveBreakerTL = P.transformation(WTL).rectangleXY(WL, WH, W);
    Region waveBreakerTLR = P.transformation(WTLR).cylinder(WR, W);
    Region waveBreakerTLC = P.transformation(WTLC, Z, alfa).rectangleXY(WC, WC, W);
    Region waveBreakerTLF = P.transformation(WTLF).rectangleXY(6*DP, WF, W);
    Region waveBreakerTLH = P.transformation(WTLH, Y).cylinder(WBHD/2,WH);

    Region deflector = P.transformation(DFP, Z, -DFA).rectangleXY(DFL, DFL, W);

    ///Operações
    if (tank && !(step && !(stepChamfer))) //&& 
      //!(waveBreakerBL /*&& !(waveBreakerBLR)*/ || (waveBreakerBLF) /*||  (waveBreakerBLC && x <= WL && y >= WBP)*/) &&     /// wave breaker bottom left
      //!((waveBreakerTL /*&& !(waveBreakerTLH)*/) /*&& !(waveBreakerTLR)*/ /*|| (waveBreakerTLF)*/ /*||  (waveBreakerTLC && x <= WL && y >= WTP)*/ || (deflector && x >= 0.0 && y <= WTP)))  /// wave breaker top left
    {
        if(fluid)
            return 0;
        else
            return -1; 
    }         

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif