#ifndef CSTLFILE_H
#define	CSTLFILE_H

#include <string.h>
#include <cstdio>
#include <iostream>

#include "stlfile.h"


class CStlFile
{
private:

    StlFile stlFile;

public:
    double xmax;
    double ymax;
    double zmax;
    double xmin;
    double ymin;
    double zmin;
    double P1x;
    double P1y;
    double P1z;
    double P2x;
    double P2y;
    double P2z;
    double P3x;
    double P3y;
    double P3z;


    StlFile::Facet *facets;
    StlFile::Stats stats;


    struct Temp
    {
    double x;
    double y;
    double z;

    Temp(double X = 0, double Y = 0, double Z = 0)
    {
            this->x = X;
            this->y = Y;
            this->z = Z;
    }

    Temp(double r[3])
    {
            this->x = r[0];
            this->y = r[1];
            this->z = r[2];
    }


    inline Temp cross(Temp B)
    {
            return Temp(this->y * B.z - this->z * B.y, - this->x * B.z + this->z * B.x, this->x * B.y - this->y * B.x);
    };

    inline double area()
    {
        return double(this->x*this->x + this->y*this->y + this->z*this->z);
    };

    };

    CStlFile(std::string fileName)
    {
        stlFile.setFormat(StlFile::ASCII); //ASCII
        stlFile.open(fileName.c_str());

        facets = stlFile.getFacets();
        stats = stlFile.getStats();
        int j = 0;
        int k = 0;
        int l = 0;

        double A1 = 0;
        double A2 = A1;
        double A3 = A1;
        for (int i = 0; i < stats.numFacets; i ++)
        {

            Temp q(facets[i].vector[0].x - facets[i].vector[1].x,facets[i].vector[0].y - facets[i].vector[1].y, facets[i].vector[0].z - facets[i].vector[1].z);
            Temp w(facets[i].vector[2].x - facets[i].vector[1].x,facets[i].vector[2].y - facets[i].vector[1].y, facets[i].vector[2].z - facets[i].vector[1].z);

            double A0 = (q.cross(w)).area();

            if (A0>A1)
            {
                A1 = A0;
                j=i;
            }
        }

        for (int i = 0; i < stats.numFacets; i ++)
        {
            if (i!=j)
            {
                Temp q(facets[i].vector[0].x - facets[i].vector[1].x,facets[i].vector[0].y - facets[i].vector[1].y, facets[i].vector[0].z - facets[i].vector[1].z);
                Temp w(facets[i].vector[2].x - facets[i].vector[1].x,facets[i].vector[2].y - facets[i].vector[1].y, facets[i].vector[2].z - facets[i].vector[1].z);

                double A0 = (q.cross(w)).area();

                if (A0>A2)
                {
                    A2 = A0;
                    k=i;
                }
            }
        }

        for (int i = 0; i < stats.numFacets; i ++)
        {
            if (i!=j and i!=k)
            {
                Temp q(facets[i].vector[0].x - facets[i].vector[1].x,facets[i].vector[0].y - facets[i].vector[1].y, facets[i].vector[0].z - facets[i].vector[1].z);
                Temp w(facets[i].vector[2].x - facets[i].vector[1].x,facets[i].vector[2].y - facets[i].vector[1].y, facets[i].vector[2].z - facets[i].vector[1].z);
                double A0 = (q.cross(w)).area();

                if (A0>A3)
                {
                    A3 = A0;
                    l=i;
                }
            }

        }

         P1x = (facets[j].vector[0].x +facets[j].vector[1].x +facets[j].vector[2].x)/3;
         P1y = (facets[j].vector[0].y +facets[j].vector[1].y +facets[j].vector[2].y)/3;
         P1z = (facets[j].vector[0].z +facets[j].vector[1].z +facets[j].vector[2].z)/3;

         P2x = (facets[k].vector[0].x +facets[k].vector[1].x +facets[l].vector[2].x)/3;
         P2y = (facets[k].vector[0].y +facets[k].vector[1].y +facets[l].vector[2].y)/3;
         P2z = (facets[k].vector[0].z +facets[k].vector[1].z +facets[l].vector[2].z)/3;

         P3x = (facets[l].vector[0].x +facets[l].vector[1].x +facets[l].vector[2].x)/3;
         P3y = (facets[l].vector[0].y +facets[l].vector[1].y +facets[l].vector[2].y)/3;
         P3z = (facets[l].vector[0].z +facets[l].vector[1].z +facets[l].vector[2].z)/3;


        //std::cout << "\n"<<j << "\t"<< k <<"\t"<< l << std::endl;
        //std::cout << "\n"<< A1 << "\t"<< A2 <<"\t"<< A3 << std::endl;
        //std::cout << "\n"<< P1x << "\t \t"<< P1y <<"\t\t"<< P1z << std::endl;
        //std::cout << "\n"<< P2x << "\t \t"<< P2y <<"\t\t"<< P2z << std::endl;
        //std::cout << "\n"<< P3x << "\t \t"<< P3y <<"\t\t"<< P3z << std::endl;

        xmax = stats.max.x;
        xmin = stats.min.x;
        ymin = stats.min.y;
        ymax = stats.max.y;
        zmin = stats.min.z;
        zmax = stats.max.z;
        Display();
    }

    void Display()
    {
        std::cout  << "\n CStlFile: numFacets = " << stats.numFacets << "\n"
                << "\t" << xmin
                << "\t      " << xmax
                << "\t" << ymin
                << "\t      " << ymax
                << "\t" << zmin
                << "\t      " << zmax
                 << "\n" << std::flush;
    }
};


#endif	/* CSTLFILE_H */

