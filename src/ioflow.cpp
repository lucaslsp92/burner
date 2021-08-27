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
    double H = 19.7;                // tank height
    double W = 10*DP;               // tank width
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
    double FH = 12.0;               // fluid height

    /// Points
    Point SCP(SL-0.8,0,-W/2);                  // step chamfer position

    Point WBL (0,WBP-WH,0);                    // wave breaker position bottom left
    Point WBLR (WR,WBP-WH,-W/2);               // wave breaker radius bottom left
    Point WBLC (0.5,WBP,-W/2);                 // wave breaker chamfer bottom left
    Point WBLF (WL-4*DP,WBP-WF/2,0);           // wave breaker front bottom left

    Point WTL (0,WTP-WH,0);                    // wave breaker position top left
    Point WTLR (WR,WTP-WH,-W/2);               // wave breaker radius top left
    Point WTLC (0.5,WTP,-W/2);                 // wave breaker chamfer top left
    Point WTLF (WL-4*DP,WTP-WF/2,0);           // wave breaker front top left

    Point WBR (L-WL,WBP-WH,-W/2);              // wave breaker position bottom right
    Point WBRR (L-WR,WBP-WH,-W/2);             // wave breaker radius bottom right
    Point WBRC (L-0.5,WBP,-W/2);               // wave breaker chamfer bottom right
    Point WBRF (L-WL-4*DP,WBP-WF/2,-W/2);      // wave breaker front bottom right

    Point WTR (L-WL,WTP-WH,-W/2);              // wave breaker position top right
    Point WTRR (L-WR,WTP-WH,-W/2);             // wave breaker radius top right
    Point WTRC (L-0.5,WTP,-W/2);               // wave breaker chamfer top right
    Point WTRF (L-WL-4*DP,WTP-WF/2,-W/2);      // wave breaker front top right

    ///Regiões (Region)
    Region tank = P.rectangleXY(L, H, W);
    Region fluid = P.rectangleXY(L, FH, DP);

    Region step = P.rectangleXY(SL, SH, W);
    Region stepChamfer = P.transformation(SCP, Z, SA).rectangleXY(SC, SC, W);

    Region waveBreakerBL = P.transformation(WBL).rectangleXY(WL, WH, DP);
    Region waveBreakerBLR = P.transformation(WBLR).cylinder(WR, W);
    Region waveBreakerBLC = P.transformation(WBLC, Z, alfa).rectangleXY(WC, WC, W);
    Region waveBreakerBLF = P.transformation(WBLF).rectangleXY(6*DP, WF, DP);

    Region waveBreakerTL = P.transformation(WTL).rectangleXY(WL, WH, DP);
    Region waveBreakerTLR = P.transformation(WTLR).cylinder(WR, W);
    Region waveBreakerTLC = P.transformation(WTLC, Z, alfa).rectangleXY(WC, WC, W);
    Region waveBreakerTLF = P.transformation(WTLF).rectangleXY(6*DP, WF, DP);

    Region waveBreakerBR = P.transformation(WBR).rectangleXY(WL, WH, W);
    Region waveBreakerBRR = P.transformation(WBRR).cylinder(WR, W);
    Region waveBreakerBRC = P.transformation(WBRC, Z, beta).rectangleXY(WC, WC, W);
    Region waveBreakerBRF = P.transformation(WBRF).rectangleXY(6*DP, WF, W);

    Region waveBreakerTR = P.transformation(WTR).rectangleXY(WL, WH, W);
    Region waveBreakerTRR = P.transformation(WTRR).cylinder(WR, W);
    Region waveBreakerTRC = P.transformation(WTRC, Z, beta).rectangleXY(WC, WC, W);
    Region waveBreakerTRF = P.transformation(WTRF).rectangleXY(6*DP, WF, W);

    if ((waveBreakerBL || (waveBreakerBLF)) ||      /// wave breaker bottom left
        (waveBreakerTL || (waveBreakerTLF)))        /// wave breaker top left
    {
        return id+2;
    }

    if (y > H)
    {
        return -1;
    }

    if (z <= 0.0 || z > DP)
    {
        return -1;
    }

    ///Return padrão (caso nenhuma operação seja feita)
    return id;

}

#endif