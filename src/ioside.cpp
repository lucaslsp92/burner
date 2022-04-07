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
    //// Moonpool sloshing
    //// Variaveis (double, Point, etc..)
    double cresthillHeight = 0.028;
    double crestHillDistance = 9.0*cresthillHeight;
    double channelHeight = 3.035*cresthillHeight;

    /// Points
    
    ///Regiões (Region)
    Region water = P.rectangleXY(crestHillDistance,channelHeight,DP);

    ///Operações
    if(water)
    {
        /// left hill
        if(x*1000>=0 && x*1000<9 && y*1000<std::min(28.0+0.006775070969851*(x*1000)*(x*1000)-0.0021245277758*(x*1000)*(x*1000)*(x*1000),28.0))
            return 2;
        else if(x*1000>=9 && x*1000<14 && y*1000<25.07355893131+0.9754803562315*(x*1000)-0.1016116352781*(x*1000)*(x*1000)+0.001889794677820*(x*1000)*(x*1000)*(x*1000))
            return 2;
        else if(x*1000>=14 && x*1000<20 && y*1000<25.79601052357+0.8206693007457*(x*1000)-0.09055370274339*(x*1000)*(x*1000)+0.001626510569859*(x*1000)*(x*1000)*(x*1000))
            return 2;
        else if(x*1000>=20 && x*1000<30 && y*1000<40.46435022819-1.379581654948*(x*1000)+0.01945884504128*(x*1000)*(x*1000)-0.0002070318932190*(x*1000)*(x*1000)*(x*1000))
            return 2;
        else if(x*1000>=30 && x*1000<40 && y*1000<17.92461334664+0.8743920332081*(x*1000)-0.05567361123058*(x*1000)*(x*1000)+0.0006277731764683*(x*1000)*(x*1000)*(x*1000))
            return 2;
        else if(x*1000>=40 && x*1000<54 && y*1000<std::max(56.39011190988-2.010520359035*(x*1000)+0.01644919857549*(x*1000)*(x*1000)+0.00002674976141766*(x*1000)*(x*1000)*(x*1000),0.0))
            return 2;

        /// right hill
        else if((crestHillDistance-x)*1000>=0 && (crestHillDistance-x)*1000<9 && y*1000<std::min(28.0+0.006775070969851*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)-0.0021245277758*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000),28.0))
            return 2;
        else if((crestHillDistance-x)*1000>=9 && (crestHillDistance-x)*1000<14 && y*1000<25.07355893131+0.9754803562315*((crestHillDistance-x)*1000)-0.1016116352781*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)+0.001889794677820*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000))
            return 2;
        else if((crestHillDistance-x)*1000>=14 && (crestHillDistance-x)*1000<20 && y*1000<25.79601052357+0.8206693007457*((crestHillDistance-x)*1000)-0.09055370274339*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)+0.001626510569859*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000))
            return 2;
        else if((crestHillDistance-x)*1000>=20 && (crestHillDistance-x)*1000<30 && y*1000<40.46435022819-1.379581654948*((crestHillDistance-x)*1000)+0.01945884504128*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)-0.0002070318932190*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000))
            return 2;
        else if((crestHillDistance-x)*1000>=30 && (crestHillDistance-x)*1000<40 && y*1000<17.92461334664+0.8743920332081*((crestHillDistance-x)*1000)-0.05567361123058*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)+0.0006277731764683*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000))
            return 2;
        else if((crestHillDistance-x)*1000>=40 && (crestHillDistance-x)*1000<54 && y*1000<std::max(56.39011190988-2.010520359035*((crestHillDistance-x)*1000)+0.01644919857549*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)+0.00002674976141766*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000)*((crestHillDistance-x)*1000),0.0))
            return 2;
        return 0;
    }

    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall
}
#endif