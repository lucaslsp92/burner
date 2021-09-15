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
    double beta = 213.690067667;    // wave breaker chamfer angle right
    double WBP = 10.5;              // wave breaker bottom position
    double WTP = 13.7;              // wave breaker top position
    double FH = 11.5;               // fluid height
    double DFL = 1.13137085;        // deflector length
    double DFA = 45;                // deflector angle
    double WBHD = 0.2;              // wave breaker hole diameter
    double WBHP = 0.4;              // wave breaker hole x and z position

    /// Points
    Point SCP(SL-0.8,0,-W/2);                  // step chamfer position

    Point WBL (0,WBP-WH/2,0);                  // wave breaker position bottom left
    Point WBLR (WR,WBP-WH/2,-W/2);             // wave breaker radius bottom left
    Point WBLC (0.5,WBP,-W/2);                 // wave breaker chamfer bottom left
    Point WBLF (WL-6*DP,WBP-WF/2,0);           // wave breaker front bottom left

    Point WTL (0,WTP-WH/2,0);                  // wave breaker position top left
    Point WTLR (WR,WTP-WH/2,-W/2);             // wave breaker radius top left
    Point WTLC (0.5,WTP,-W/2);                 // wave breaker chamfer top left
    Point WTLF (WL-6*DP,WTP-WF/2,0);           // wave breaker front top left
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

    if (y > H)
    {
        return -1;
    }

    /*if (z <= 0.0 || z > DP)
    {
        return -1;
    }*/
    
    if (x < 0.0 && y <= H)                                            /// left wall
    {
        return id+2;
    }

    //if (x > 0.0 && x <= WL && y > WTP-WL-WF && y < WTP && z >= 0.0 && z <= W /*&& ((waveBreakerTL && !(waveBreakerTLH)) || waveBreakerTLF)*/)        /// wave breaker top left (bottom half)
    /*{
        return id+4;
    }*/

    //if (x > 0.0 && x <= WL && y > WTP && y < WTP+WF && z >= 0.0 && z <= W && ((waveBreakerTL /*&& !(waveBreakerTLH)*/) || waveBreakerTLF))        /// wave breaker top left (top half)
    /*{
        return id+6;
    }*/

    /*if (x > 0.0 && x <= WL && y > WBP-WF && y < WBP && z >= 0.0 && z <= W && (waveBreakerBL || waveBreakerBLF))        /// wave breaker bottom left (bottom half)
    {
        return id+8;
    } 

    if (x > 0.0 && x <= WL && y > WBP && y < WBP+WF && z >= 0.0 && z <= W && (waveBreakerBL || waveBreakerBLF))        /// wave breaker bottom left (top half)
    {
        return id+10;
    }*/

    /*if (y > H/2)
    {
        return id+12;
    }*/  

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif