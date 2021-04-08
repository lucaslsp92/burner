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
    //// Dam-break 2D with internal solid //// Created only on plane XY
    ///Variaveis (double, Point, etc..)
    // Tank wall dimensions (meters)
    double Wx = 2.0;
    double Wy = 1.0;
    double Wz = DP;         // Distance between particles

    // Fluid column dimensions (meters)
    double Fx = 0.3*Wx;     // Column length 30% of tank length
    double Fy = 0.4*Wy;     // Column height 40% of tank height
    double Fz = DP;         // Distance between particles

    // Solid dimensions (meters)
    double Sx = 0.5;
    double Sy = 0.4;
    double Sz = DP;         // Distance between particles

    // Solid position (meters)
    Point A (0.25*Wx,0.4*Wy,0);
    // Solid rotation around Z axis
    Point Rotation (0,0,1);
    // Angle of rotation (degrees)
    double theta = 30;

    ///Regiões (Region)
    // Geometries and geometric transformations from file Point.h
    Region wall  = P.rectangleXY(Wx, Wy, Wz);
	Region fluid = P.rectangleXY(Fx, Fy, Fz);
	// Parallelepiped translated for point A and rotated by theta degrees around Z axis
	// For 2D geometries, Burner don't create internal or external dummies for internal solids.
	// The user have to create the empty region and the intern or external dummies.
	Region solid    = P.transformation(A, Rotation, theta).rectangleXYCenter(Sx, Sy, Sz);            // Solid external walls
	Region solidint = P.transformation(A, Rotation, theta).rectangleXYCenter(Sx-2*DP, Sy-2*DP, Sz);  // Solid internal walls
	Region solidemp = P.transformation(A, Rotation, theta).rectangleXYCenter(Sx-6*DP, Sy-6*DP, Sz);  // Solid empty region

    ///Operações
    if (wall) // External tank wall
    {
        if (solidemp) return -1;    // Internal solid - empty region
        if (solidint) return 5;     // Internal solid - internal walls (dummies)
        if (solid) return 4;        // Internal solid - external walls
    	if (fluid) return 0;        // Fluid Column

    	return -1;                  // Empty region
    }
    //// Dam-break 2D with internal solid //// Created only on plane XY
    //////////////////////////////////////////

/*
    //////////////////////////////////////////
    //// Dam-break 3D with internal solid ////
    ///Variaveis (double, Point, etc..)
    // Tank wall dimensions (meters)
    double Wx = 2.0;
    double Wy = 1.0;
    double Wz = 1.0;

    // Fluid column dimensions (meters)
    double Fx = 0.3*Wx;     // Column length 30% of tank length
    double Fy = Wy;
    double Fz = 0.4*Wz;     // Column height 40% of tank height

    // Solid dimensions (meters)
    double Sx = 0.5;
    double Sy = 0.4;
    double Sz = 0.3;

    // Solid position (meters)
    Point A (0.5*Wx,0.3*Wy,0.3*Wz);
    // Solid rotation around Z axis
    Point Rotation (0,0,1);
    // Angle of rotation (degrees)
    double theta = 30;

    ///Regiões (Region)
    // Geometries and geometric transformations from file Point.h
    Region wall = P.parallelepiped (Wx, Wy, Wz);
	Region fluid = P.parallelepiped (Fx, Fy, Fz);
	// Parallelepiped translated for point A and rotated by 30 degrees around Z axis
	// For 3D geometries, Burner create internal or external dummies for internal solids.
	Region solid = P.transformation(A, Rotation, theta).parallelepiped(Sx, Sy, Sz);

    ///Operações
    if (wall) // External tank wall
    {
        if (solid) return 4;    // Internal solid
    	if (fluid) return 0;    // Fluid Column

    	return -1;              // Empty region
    }
    //// Dam-break 3D with internal solid ////
    //////////////////////////////////////////
*/
    ///Return padrão (constrói a parede externa)
    // DO NOT CHANGE HERE !!!
    return 2; // External wall

}

#endif

