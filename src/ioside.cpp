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
    int dim = 2;
    double L = 38.4;                // tank length
    double H = 2*19.7;              // tank height
    //double H = 19.7;                // tank height
    double W = 10*DP;               // tank width
    double WH = 6*DP;               // wave breaker height
    //double H = 2*19.7;              // tank height
    //double W = 0.8;                 // tank width
    //double H = 19.7;                // tank height
    //double W = 2.0;                 // tank width
    //double WH = 4*DP;               // wave breaker height
    double SL = 12.8;               // step length
    double SH = 7.3;                // step height
    double SC = 2.912043956;        // step chamfer length
    double SA = 15.945395774;       // stpe chamfer angle
    double WL = 0.8;                // wave breaker length
    double WF = 0.5;                // wave breaker front
    double WBP = 10.5;              // wave breaker bottom position
    double WTP = 13.7;              // wave breaker top position
    double FH = 11.5;               // fluid height
    double DFL = 1.13137085;        // deflector length
    double DFA = 45;                // deflector angle
    double WBHD = 0.2;              // wave breaker hole diameter
    double WBHP = 0.4;              // wave breaker hole x and z position
    double VDH = 1.0;               // division height
    double VDP = 3.2;               // division position
    double BTL = 26.4;              // bottom tank length
    double BTH = 1.6;               // bottom tank height

    /// Points
    Point SCP2D(SL-0.8,0,-W/2);                 // step chamfer position
    Point WBL2D(0,WBP-WH/2,-W/2);               // wave breaker position bottom 
    Point WBLF2D(WL-WH,WBP-WF/2,-W/2);          // wave breaker front bottom 
    Point WTL2D(0,WTP-WH/2,-W/2);               // wave breaker position top 
    Point WTLF2D(WL-WH,WTP-WF/2,-W/2);          // wave breaker front top 
    Point DFP2D(0,WTP-WL-WH/2,-W/2);            // deflector position 
    Point TDP2D(SL-WH-VDP,SH,-W/2);                 // tank division position
    Point BTT2D(SL-0.8,-BTH,-W/2);              // bottom tank position
    Point EXT2D(SL-0.8,-BTH*2,-W/2);            // bottom tank extension

    Point SCP3D(SL-0.8,0,0);                 // step chamfer position
    Point WBL3D(0,WBP-WH/2,0);               // wave breaker position bottom 
    Point WBLF3D(WL-WH,WBP-WF/2,0);          // wave breaker front bottom 
    Point WTL3D(0,WTP-WH/2,0);               // wave breaker position top 
    Point WTLF3D(WL-WH,WTP-WF/2,0);          // wave breaker front top 
    Point WTLH3D(WBHP, WTP-WH/2, WBHP);      // wave breaker hole top 
    Point DFP3D(0,WTP-WL-WH/2,0);            // deflector position 

    ///Regiões (Region)
    Region tank2D = P.rectangleXY(L, H, W);
    Region fluid2D = P.rectangleXY(L, FH, DP);
    Region step2D = P.rectangleXY(SL, SH, W);
    Region stepChamfer2D = P.transformation(SCP2D, Z, SA).rectangleXY(SC, SC, W);
    Region waveBreakerBL2D = P.transformation(WBL2D).rectangleXY(WL, WH, W);
    Region waveBreakerBLF2D = P.transformation(WBLF2D).rectangleXY(WH, WF, W);
    Region waveBreakerTL2D = P.transformation(WTL2D).rectangleXY(WL, WH, W);
    Region waveBreakerTLF2D = P.transformation(WTLF2D).rectangleXY(WH, WF, W);
    Region deflector2D = P.transformation(DFP2D, Z, -DFA).rectangleXY(DFL, DFL, W);
    Region division2D = P.transformation(TDP2D).rectangleXY(WH, VDH, W);
    Region bottomTank2D = P.transformation(BTT2D).rectangleXY(BTL, BTH, DP);
    Region extension2D = P.transformation(EXT2D).rectangleXY(BTL, BTH-3*DP, W);

    Region tank3D = P.rectangleXY(L, H, W);
    Region fluid3D = P.rectangleXY(L, FH, W);
    Region step3D = P.rectangleXY(SL, SH, W);
    Region stepChamfer3D = P.transformation(SCP3D, Z, SA).rectangleXY(SC, SC, W);
    Region waveBreakerBL3D = P.transformation(WBL3D).rectangleXY(WL, WH, W);
    Region waveBreakerBLF3D = P.transformation(WBLF3D).rectangleXY(WH, WF, W);
    Region waveBreakerTL3D = P.transformation(WTL3D).rectangleXY(WL, WH, W);
    Region waveBreakerTLF3D = P.transformation(WTLF3D).rectangleXY(WH, WF, W);
    Region waveBreakerTLH3D = P.transformation(WTLH3D, Y).cylinder(WBHD/2,WH);
    Region deflector3D = P.transformation(DFP3D, Z, -DFA).rectangleXY(DFL, DFL, W);
    Region bottomTank3D = P.transformation(BTT2D).rectangleXY(BTL, BTH, W);

    ///Operações
    if (dim == 2)
    {
        if ((tank2D && !(step2D && !(stepChamfer3D)) //&& 
          //!(waveBreakerBL2D || (waveBreakerBLF2D)) &&                                                      /// wave breaker bottom 
          /*!((waveBreakerTL2D) || (waveBreakerTLF2D) || (deflector2D && x >= 0.0 && y <= WTP)) &&*/           /// wave breaker top 
          /*!(division2D)*/) || bottomTank3D || extension2D)
        {
            if(fluid2D || bottomTank2D)
                return 0;
            else
                return -1; 
        }  
    }
    else if (dim == 3)
    {
        if (tank3D && !(step3D && !(stepChamfer3D)))// && 
          //!(waveBreakerBL3D || (waveBreakerBLF3D)) &&                                                                           /// wave breaker bottom 
          //!((waveBreakerTL3D && !(waveBreakerTLH3D)) || (waveBreakerTLF3D) /*|| (deflector3D && x >= 0.0 && y <= WTP)*/))           /// wave breaker top 
        {
            if(fluid3D)
                return 0;
            else
                return -1; 
        }  
    }           

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif