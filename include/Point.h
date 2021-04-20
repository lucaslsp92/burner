#ifndef H_POINT
#define H_POINT

#include "define.h"

struct Point
{
    double x;
    double y;
    double z;

    Point (double X = 0, double Y = 0, double Z = 0)
    {
        this->x = X;
        this->y = Y;
        this->z = Z;
    }

    Point (long r[3], double DP)
    {
        this->x = (r[0]+POSITIONX*0.5)*DP;
        this->y = (r[1]+POSITIONY*0.5)*DP;
        this->z = (r[2]+POSITIONZ*0.5)*DP;
    }

    Point (double r[3], double DP)
    {
        this->x = (r[0]+POSITIONX*0.5)*DP;
        this->y = (r[1]+POSITIONY*0.5)*DP;
        this->z = (r[2]+POSITIONZ*0.5)*DP;
    }

    Point (Point P, double X, double Y, double Z = 0)
    {
        this->x = P.x + X;
        this->y = P.y + Y;
        this->z = P.z + Z;
    }

    Point transformation (Point O, Point Vector = Point (0,0,1), double rotation = 0)
    {
        Point P;

        double senTeta, cosTeta, senFi, cosFi;
        double projXY, norma;

        norma = sqrt(pow(Vector.x, 2) + pow(Vector.y, 2) + pow(Vector.z, 2));
        if (fabs(norma) <= EPS) {
            printf ("Error: ||Vector|| = 0!\n");
        }
        projXY = sqrt(pow(Vector.x, 2) + pow(Vector.y, 2));
        if (!(fabs(projXY) <= EPS)) {

            senTeta = (Vector.y)/projXY;
            cosTeta = (Vector.x)/projXY;

            senFi = (Vector.z)/norma;
            cosFi = projXY/norma;

            P.x = ((x-O.x)*cosTeta + (y-O.y)*senTeta)*senFi - (z-O.z)*cosFi;
            P.y = ((y-O.y)*cosTeta - (x-O.x)*senTeta);
            P.z = (z-O.z)*senFi + ((x-O.x)*cosTeta + (y-O.y)*senTeta)*cosFi;

        } else {
            P.x = x - O.x;
            P.y = y - O.y;
            if (Vector.z > 0) P.z = z - O.z;
            else P.z = O.z - z;
        }

        double alpha;

        if (P.x == 0) {
            if (P.y > 0) alpha = PI/2;
            else if (P.y < 0) alpha = 3*PI/2;
        }
        else if (P.x > 0)
            if (P.y > 0) alpha = atan(P.y/P.x);
            else alpha = atan(P.y/P.x) + 2*PI;
        else if (P.x < 0) alpha = atan(P.y/P.x) + PI;

        if (rotation != 0 && rotation != 180) {
            double radius = sqrt(pow(P.x,2)+pow(P.y,2));
            double angle = alpha + rotation*PI/180;
            P.x = radius*cos(angle);
            P.y = radius*sin(angle);
        } else if (rotation == 180) {
            P.x = -P.x;
            P.y = -P.y;
        }

        return P;

    }

    Region region2D (double DP)
    {
        if (z > -DP && z < DP)
            return true;

        return false;
    }

    Region sphere (double R)
    {
        if (x*x + y*y + z*z < R*R)
            return true;

        return false;
    }

    Region toroid (double R, double r)
    {
        if (pow(R-sqrt(x*x + y*y),2) + z*z < r*r) return true;

        return false;
    }

    Region helix (double R, double r, double d)
    {
        double alpha;

        if (x == 0) {
            if (y > 0) alpha = PI/2;
            else if (y < 0) alpha = 3*PI/2;
        }
        else if (x > 0)
            if (y > 0) alpha = atan(y/x);
            else alpha = atan(y/x) + 2*PI;
        else if (x < 0) alpha = atan(y/x) + PI;

        alpha = alpha*180/PI;

        if (pow(R-sqrt(x*x + y*y),2) + pow(z - (alpha/360)*(d),2) < r*r) return true;

        return false;

    }

    Region cylinder (double R, double H)
    {
        if (z > 0 && z < H)
            if (x*x + y*y < R*R)
                return true;

        return false;
    }

    Region extrudeEllipse (double a, double b, double R, double H)
    {
        if (z > 0 && z < H)
            if (x*x/(a*a) + y*y/(b*b) < R*R)
                return true;

        return false;
    }

    Region heart (double H)
    {
        if (z > 0 && z < H)
            if (pow (x*x + y*y -1,3) < x*x*y*y*y)
                return true;

        return false;
    }

    Region quadrica (double A, double B, double C, double D, double E, double F, double G, double H, double I, double J, double R)
    {
        if (x*x + y*y + z*z < R*R)
            if (A*x*x + B*y*y + C*z*z + D*x*y + E*y*z + F*z*x + G*x + H*y + I*z + J < 0)
                return true;

        return false;
    }

    Region ellipsoid (double a, double b, double c, double R)
    {
        if (x*x/(a*a) + y*y/(b*b) + z*z/(c*c) < R*R)
            return true;

        return false;
    }

    Region paraboloid (double a, double b, double H)
    {
        if (x*x/(a*a) + y*y/(b*b) < z && z < H)
            return true;

        return false;
    }

    Region rectangleXYCenter (double X, double Y, double Z)
    {
        if (x > -X/2 && x < X/2)
            if (y > -Y/2 && y < Y/2)
                if (z > -Z/2 && z < Z/2)
                    return true;

        return false;
    }

    Region rectangleXY (double X, double Y, double Z)
    {
        if (x > 0 && x < X)
            if (y > 0 && y < Y)
                if (z > 0 && z < Z)
                    return true;

        return false;
    }

    Region parallelepipedCenter (double X, double Y, double Z)
    {
        if (x > -X/2 && x < X/2)
            if (y > -Y/2 && y < Y/2)
                if (z > -Z/2 && z < Z/2)
                    return true;

        return false;
    }

    Region parallelepiped (double X, double Y, double Z)
    {
        if (x > 0 && x < X)
            if (y > 0 && y < Y)
                if (z > 0 && z < Z)
                    return true;

        return false;
    }

    Region ring (double R, double r, double H)
    {
        if (z > 0 && z < H)
            if (x*x + y*y < R*R && x*x + y*y > r*r)

        return false;
    }

    Region wigley (double L, double B, double d)
    {
        if (z < 0 && z > -d)
            if (fabs(x) < L/2)
                if (fabs(y) < 0.5*B*(1-pow(2*x/L,2))*(1-pow(z/d,2)))
                    return true;

        return false;
    }

    Region wigley (double L, double B, double d, double a2, double a4, double alpha)
    {

        double X = x*2/L;
        double Y = y*2/B;
        double Z = -z/d;

        if (Z > 0 && Z < 1)
			if (fabs(X) < 1)
				if (fabs(Y) < (1 - Z*Z)*(1 - X*X)*(1 + a2*X*X + a4*pow(X,4)) + alpha*Z*Z*(1 - pow(Z,8))*pow(1-X*X, 4))
					return true;

        return false;
    }

    Region podalyro (double a, double b, double c, double h)
    {
        if (z > 0&& z <= h)
            if (pow(x-(z/h*(a-c)+c),2)/pow(z/h*(a-c)+c,2) + y*y/pow(z/h*(b-c)+c,2) < 1)
                return true;

        return false;
    }

    Region podalyro2 (double a, double b, double r, double h, double d)
    {
        double da = d*sqrt(h*h + (a-r)*(a-r))/h;
        double db = d*sqrt(h*h + (b-r)*(b-r))/h;

        if (z > 0 && z < h)
            if (x*x/pow(z/h*(a-r)+r,2) + y*y/pow(z/h*(b-r)+r,2) < 1)
                if (x*x/pow(z/h*(a-r)+r - da,2) + y*y/pow(z/h*(b-r)+r - db,2) > 1)
                    return true;


        return false;
    }

	Region parshall (double W, double A, double B, double C, double D, double E, double F, double G, double K, double N)
	{
		if (x > 0 && x <= B)
			if (fabs(y) < 0.5*D + x * (W-D)/(2*B))
				if (z > 0 && z < E)
					return true;

		if (x > B && x <= B + F)
			if (fabs(y) < 0.5*W)
				if (z > -N/F*(x-B) && z < E)
					return true;

		if (x > B + F && x < B + F + G)
			if (fabs(y) < (C-W)/(2*G) * (x-(B+F)) + W/2)
				if (z > (-K+N)/G*(x - (B+F)) - N && z < E)
					return true;

		return false;

	}

	double plane (Point O, Point Vec)
	{
	    return Vec.x * ( x - O.x ) + Vec.y * ( y - O.y ) + Vec.z * ( z - O.z );
	}

};

#endif

