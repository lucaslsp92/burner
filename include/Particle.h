#ifndef H_PARTICLE
#define H_PARTICLE

#include "define.h"
#include <stdlib.h>
#include <time.h>

struct Particle
{
    int id;
    long r[3];

    Particle ()
    {
        id = -1;
    }

    ~Particle ()
    {

    }

    int isID (int ID, int listID[NUMMAXID], int numID)
    {
        for (int n = 0; n < numID; n++)
            if (ID == listID[n]) return 1;

        return 0;
    }

    void setCell (Particle zero)
    {

        id = 0;
        for (int n = 0; n < 3; n++)
            r[n] = zero.r[n];

    }

    void setParticle (int ID, Particle zero, int R[])
    {

        id = ID;
        for (int n = 0; n < 3; n++)
            r[n] = zero.r[n] + R[n];

    }

    int ioside (double DP, double geometry[]);

    int ioflow (double DP, double geometry[]);

    void getMatrix (int Matrix[CELLDIAMETER+2*CELLBORDER][CELLDIAMETER+2*CELLBORDER][CELLDIAMETER+2*CELLBORDER], double DP, double geometry[])
    {
        Particle zero;

        for (int d = 0; d < 3; d++)
            zero.r[d] = r[d] - CELLBORDER;

#pragma omp parallel for

        for (int i = 0;  i < CELLDIAMETER+2*CELLBORDER; i++)
            for (int j = 0; j < CELLDIAMETER+2*CELLBORDER; j++)
                for (int k = 0; k < CELLDIAMETER+2*CELLBORDER; k++) {
                    Particle temp;
                    int R[3] = {i, j, k};
                    temp.setParticle(0, zero, R);
                    Matrix[R[0]][R[1]][R[2]] = temp.ioside(DP, geometry);
                }

    }

    long getParticles (FILE* file, Particle particle[], long contParticles, int Matrix[CELLDIAMETER+2*CELLBORDER][CELLDIAMETER+2*CELLBORDER][CELLDIAMETER+2*CELLBORDER], int fluidID[], int numFluidID, int solidID[], int numSolidID, double DP, double geometry[])
    {

        long contAux = contParticles;

        int R[3];

        for (R[0] = 0; R[0] < CELLDIAMETER; R[0]++)
            for (R[1] = 0; R[1] < CELLDIAMETER; R[1]++)
                for (R[2] = 0; R[2] < CELLDIAMETER; R[2]++)
                {
                    if (isID (Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER], fluidID, numFluidID))
                    {
                        particle[contParticles].setParticle (Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER], *this, R);
                        Point X (particle[contParticles].r , DP);
                        fprintf (file, "%d\t%lf\t%lf\t%lf\n", particle[contParticles].id, X.x, X.y, X.z);
                        contParticles++;
                    } else if (isID (Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER], solidID, numSolidID))
                    {

                        int flag = 0;
                        for (int I = R[0]+CELLBORDER-1; I <= R[0]+CELLBORDER+1 && flag == 0; I++)
                            for (int J = R[1]+CELLBORDER-1; J <= R[1]+CELLBORDER+1 && flag == 0; J++)
                                for (int K = R[2]+CELLBORDER-1; K <= R[2]+CELLBORDER+1 && flag == 0; K++)
                                    if (Matrix[I][J][K] != Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER])
                                        flag = 1;

                        if (flag == 1) {
                            particle[contParticles].setParticle (Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER], *this, R);
                            particle[contParticles].id = particle[contParticles].ioflow(DP, geometry);
                            if (particle[contParticles].id >= 0) {
                                Point X (particle[contParticles].r , DP);
                                fprintf (file, "%d\t%lf\t%lf\t%lf\n", particle[contParticles].id, X.x, X.y, X.z);
                                contParticles ++;
                            }
                        } else {
                            for (int I = R[0]+CELLBORDER-3; I <= R[0]+CELLBORDER+3 && flag == 0; I++)
                                for (int J = R[1]+CELLBORDER-3; J <= R[1]+CELLBORDER+3 && flag == 0; J++)
                                    for (int K = R[2]+CELLBORDER-3; K <= R[2]+CELLBORDER+3 && flag == 0; K++)
                                        if (Matrix[I][J][K] != Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER])
                                            flag = 1;

                            if (flag == 1) {
                                particle[contParticles].setParticle (Matrix[R[0]+CELLBORDER][R[1]+CELLBORDER][R[2]+CELLBORDER]+1, *this, R);
                                particle[contParticles].id = particle[contParticles].ioflow(DP, geometry);
                                if (particle[contParticles].id >= 0) {
                                    Point X (particle[contParticles].r , DP);
                                    fprintf (file, "%d\t%lf\t%lf\t%lf\n", particle[contParticles].id, X.x, X.y, X.z);
                                    contParticles ++;
                                }
                            }

                        }
                    }
                }

        return contParticles-contAux;

    }

    int cellNeighbor (Particle particlezero, int cube, int contCellNeighbor)
    {

        int R[3];

        int diameter = 2*cube;

        R[0] = int(floor(contCellNeighbor/pow(diameter+1,2)));
        R[1] = int(floor((contCellNeighbor-R[0]*pow(diameter+1,2))/(diameter+1)));
        R[2] = int(contCellNeighbor-R[0]*pow(diameter+1,2)-R[1]*(diameter+1));

        for (int i = 0; i < 3; i++)
            r[i] = particlezero.r[i] + (R[i]-cube)*CELLDIAMETER;

        if (fabs(R[0]-cube)+fabs(R[1]-cube)+fabs(R[2]-cube) == 1) return 1;
        else return 0;

    }

    int comparesCells (Particle tempCell)
    {
        if (tempCell.r[0] == r[0])
            if (tempCell.r[1] == r[1])
                if (tempCell.r[2] == r[2])
                    return 1;

        return 0;
    }

    long burnCells (FILE *file, Particle particle[], double DP, double geometry[], int fluidID[], int numFluidID, int solidID[], int numSolidID)
    {

        struct Particle *cell = ( struct Particle* ) malloc ( NUMMAXCELLS * sizeof(struct Particle));

        printf ("\nInitialization of cells: ");
        if (cell) {
            printf ("Done!\n\n");
        }
        else
        {
            printf ("Failure!");
            system("PAUSE");
        }

        cell = new Particle[NUMMAXCELLS];

        int contCells = 0;
        cell[contCells].setCell(*this);
        contCells++;

        long contParticles = 0;

        int badCells = 0;

        int Matrix[CELLDIAMETER+2*CELLBORDER][CELLDIAMETER+2*CELLBORDER][CELLDIAMETER+2*CELLBORDER];

        for (long n = 0; n < contCells; n++)
        {
            cell[n].getMatrix(Matrix, DP, geometry);
            cell[n].id = cell[n].getParticles(file, particle, contParticles, Matrix, fluidID, numFluidID, solidID, numSolidID, DP, geometry);
            contParticles += cell[n].id;
            if (cell[n].id == 0)
            {
                badCells++;
            }
            else
            {
                for (int contCellNeighbor = 0; contCellNeighbor < 27; contCellNeighbor++)
                {
                    Particle cellTemp;
                    if(cellTemp.cellNeighbor(cell[n], SIMPLECUBE, contCellNeighbor))
                    {
                        int compares, flag;

                        for (compares = 0, flag = 0; compares < contCells && flag == 0; compares++)
                            flag += cellTemp.comparesCells(cell[compares]);

                        if (!flag)
                        {
                            cell[contCells].setCell(cellTemp);
                            contCells++;
                        }
                    }
                }
            }

            printf ("\rcontParticles = %ld", contParticles);
        }

        printf ("\n\nbadCells = %d\n", badCells);
        printf ("goodCells = %d\n\n", contCells-badCells);

        printf ("\a");

        return contParticles;

    }

    void print_Paraview (FILE *paraview_xyz, long cont, Particle particle[], double DP, double geometry[], int numgeometry, int decimal)
    {
        double theta = 45*M_PI/180;                 // grating angle
        double gx = 0.90;                           // grating position X
        double gs = 0.010;                          // grating spacing
        double gb = 0.010;                          // space between tank bottom and first grating element
        double gw = 0.005;                          // grating width
        double gsy = (gs+gw)/cos(theta);            // grating spacing Y
        double gm = 0.013767767;                    // move granting downward Y
        double gy;

        fprintf (paraview_xyz, "%ld\n", cont);

        for (int i = 0; i < numgeometry; i++)
            fprintf (paraview_xyz, "geometry[%d] = %9lf   ", i, geometry[i]);

        for (long i = 0; i < cont; i++)
        {
            Point X (particle[i].r , DP);
            if (particle[i].id < 4)
            {
                if (decimal == 2)
                    fprintf (paraview_xyz, "\n%d %9.2lf %9.2lf %9.2lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 3)
                    fprintf (paraview_xyz, "\n%d %9.3lf %9.3lf %9.3lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 4)
                    fprintf (paraview_xyz, "\n%d %9.4lf %9.4lf %9.4lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 5)
                    fprintf (paraview_xyz, "\n%d %9.5lf %9.5lf %9.5lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 6)
                    fprintf (paraview_xyz, "\n%d %9.6lf %9.6lf %9.6lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 7)
                    fprintf (paraview_xyz, "\n%d %9.7lf %9.7lf %9.7lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 8)
                    fprintf (paraview_xyz, "\n%d %9.8lf %9.8lf %9.8lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 9)
                    fprintf (paraview_xyz, "\n%d %9.9lf %9.9lf %9.9lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 10)
                    fprintf (paraview_xyz, "\n%d %9.10lf %9.10lf %9.10lf", particle[i].id, X.x, X.y, 0.0);
            }
            else
            {
                if (X.y < gb+gw+0*gsy)
                    gy = gb+0.5*gw+0*gsy;                  // grating position Y
                else if (X.y < gb+gw+1*gsy)
                    gy = gb+0.5*gw+1*gsy;                  // grating position Y
                else if (X.y < gb+gw+2*gsy)
                    gy = gb+0.5*gw+2*gsy;                  // grating position Y
                else if (X.y < gb+gw+3*gsy)
                    gy = gb+0.5*gw+3*gsy;                  // grating position Y
                else if (X.y < gb+gw+4*gsy)
                    gy = gb+0.5*gw+4*gsy;                  // grating position Y
                else if (X.y < gb+gw+5*gsy)
                    gy = gb+0.5*gw+5*gsy;                  // grating position Y
                else if (X.y < gb+gw+6*gsy)
                    gy = gb+0.5*gw+6*gsy;                  // grating position Y
                else if (X.y < gb+gw+7*gsy)
                    gy = gb+0.5*gw+7*gsy;                  // grating position Y
                else if (X.y < gb+gw+8*gsy)
                    gy = gb+0.5*gw+8*gsy;                  // grating position Y
                else if (X.y < gb+gw+9*gsy)
                    gy = gb+0.5*gw+9*gsy;                  // grating position Y
                else if (X.y < gb+gw+10*gsy)
                    gy = gb+0.5*gw+10*gsy;                 // grating position Y
                else if (X.y < gb+gw+11*gsy)
                    gy = gb+0.5*gw+11*gsy;                 // grating position Y
                else if (X.y < gb+gw+12*gsy)
                    gy = gb+0.5*gw+12*gsy;                 // grating position Y
                else if (X.y < gb+gw+13*gsy)
                    gy = gb+0.5*gw+13*gsy;                 // grating position Y
                else if (X.y < gb+gw+14*gsy)
                    gy = gb+0.5*gw+14*gsy;                 // grating position Y
                else if (X.y < gb+gw+15*gsy)
                    gy = gb+0.5*gw+15*gsy;                 // grating position Y
                else if (X.y < gb+gw+16*gsy)
                    gy = gb+0.5*gw+16*gsy;                 // grating position Y
                else if (X.y < gb+gw+17*gsy)
                    gy = gb+0.5*gw+17*gsy;                 // grating position Y
                else if (X.y < gb+gw+18*gsy)
                    gy = gb+0.5*gw+18*gsy;                 // grating position Y
                else if (X.y < gb+gw+19*gsy)
                    gy = gb+0.5*gw+19*gsy;                 // grating position Y
                else if (X.y < gb+gw+20*gsy)
                    gy = gb+0.5*gw+20*gsy;                 // grating position Y
                else if (X.y < gb+gw+21*gsy)
                    gy = gb+0.5*gw+21*gsy;                 // grating position Y
                else if (X.y < gb+gw+22*gsy)
                    gy = gb+0.5*gw+22*gsy;                 // grating position Y
                else if (X.y < gb+gw+23*gsy)
                    gy = gb+0.5*gw+23*gsy;                 // grating position Y
                else if (X.y < gb+gw+24*gsy)
                    gy = gb+0.5*gw+24*gsy;                 // grating position Y
                else if (X.y < gb+gw+25*gsy)
                    gy = gb+0.5*gw+25*gsy;                 // grating position Y
                else if (X.y < gb+gw+26*gsy)
                    gy = gb+0.5*gw+26*gsy;                 // grating position Y
                else if (X.y < gb+gw+27*gsy)
                    gy = gb+0.5*gw+27*gsy;                 // grating position Y
                else if (X.y < gb+gw+28*gsy)
                    gy = gb+0.5*gw+28*gsy;                 // grating position Y
                else if (X.y < gb+gw+29*gsy)
                    gy = gb+0.5*gw+29*gsy;                 // grating position Y
                else if (X.y < gb+gw+30*gsy)
                    gy = gb+0.5*gw+30*gsy;                 // grating position Y

                if (decimal == 2)
                    fprintf (paraview_xyz, "\n%d %9.2lf %9.2lf %9.2lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 3)
                    fprintf (paraview_xyz, "\n%d %9.3lf %9.3lf %9.3lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 4)
                    fprintf (paraview_xyz, "\n%d %9.4lf %9.4lf %9.4lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 5)
                    fprintf (paraview_xyz, "\n%d %9.5lf %9.5lf %9.5lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 6)
                    fprintf (paraview_xyz, "\n%d %9.6lf %9.6lf %9.6lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 7)
                    fprintf (paraview_xyz, "\n%d %9.7lf %9.7lf %9.7lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 8)
                    fprintf (paraview_xyz, "\n%d %9.8lf %9.8lf %9.8lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 9)
                    fprintf (paraview_xyz, "\n%d %9.9lf %9.9lf %9.9lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 10)
                    fprintf (paraview_xyz, "\n%d %9.10lf %9.10lf %9.10lf", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
            }
            
        }

        printf ("paraview.xyz: Done!\n");

        fclose (paraview_xyz);

        printf ("\a\a\a");
    }

    void print_Paraview_2D (FILE *paraview_2D_xyz, long cont, Particle particle[], double DP, double geometry[], int numgeometry, long section, char d, int decimal)
    {
        long cont2D = 0;
        int dim = d-120;

        for (long i = 0; i < cont; i++)
            if (particle[i].r[dim] == section)
                cont2D++;

        fprintf (paraview_2D_xyz, "%ld\n", cont2D);

        for (int i = 0; i < numgeometry; i++)
            fprintf (paraview_2D_xyz, "geometry[%d] = %9lf   ", i, geometry[i]);

        for (long i = 0; i < cont; i++)
            if (particle[i].r[dim] == section)
            {
                Point X (particle[i].r , DP);
                if (decimal == 2)
                    fprintf (paraview_2D_xyz, "\n%d %9.2lf %9.2lf %9.2lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 3)
                    fprintf (paraview_2D_xyz, "\n%d %9.3lf %9.3lf %9.3lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 4)
                    fprintf (paraview_2D_xyz, "\n%d %9.4lf %9.4lf %9.4lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 5)
                    fprintf (paraview_2D_xyz, "\n%d %9.5lf %9.5lf %9.5lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 6)
                    fprintf (paraview_2D_xyz, "\n%d %9.6lf %9.6lf %9.6lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 7)
                    fprintf (paraview_2D_xyz, "\n%d %9.7lf %9.7lf %9.7lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 8)
                    fprintf (paraview_2D_xyz, "\n%d %9.8lf %9.8lf %9.8lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 9)
                    fprintf (paraview_2D_xyz, "\n%d %9.9lf %9.9lf %9.9lf", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 10)
                    fprintf (paraview_2D_xyz, "\n%d %9.10lf %9.10lf %9.10lf", particle[i].id, X.x, X.y, 0.0);
            }

        printf ("paraview_2D.xyz: Done!\n");
    }

    void print_MPS (FILE *mps_grid, long cont, Particle particle[], double DP, int decimal)
    {
        double theta = 45*M_PI/180;                 // grating angle
        double gx = 0.90;                           // grating position X
        double gs = 0.010;                          // grating spacing
        double gb = 0.010;                          // space between tank bottom and first grating element
        double gw = 0.005;                          // grating width
        double gsy = (gs+gw)/cos(theta);            // grating spacing Y
        double gm = 0.013767767;                    // move granting downward Y
        double gy;

        fprintf (mps_grid, "%d\n%ld", 0, cont);

        for (long i = 0; i < cont; i++)
        {
            Point X (particle[i].r , DP);
            if (particle[i].id < 4)
            {
                if (decimal == 2)
                    fprintf (mps_grid, "\n%d %9.2lf %9.2lf %9.2lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 3)
                    fprintf (mps_grid, "\n%d %9.3lf %9.3lf %9.3lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 4)
                    fprintf (mps_grid, "\n%d %9.4lf %9.4lf %9.4lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 5)
                    fprintf (mps_grid, "\n%d %9.5lf %9.5lf %9.5lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 6)
                    fprintf (mps_grid, "\n%d %9.6lf %9.6lf %9.6lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 7)
                    fprintf (mps_grid, "\n%d %9.7lf %9.7lf %9.7lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 8)
                    fprintf (mps_grid, "\n%d %9.8lf %9.8lf %9.8lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 9)
                    fprintf (mps_grid, "\n%d %9.9lf %9.9lf %9.9lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 10)
                    fprintf (mps_grid, "\n%d %9.10lf %9.10lf %9.10lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
            }
            else
            {
                if (X.y < gb+gw+0*gsy)
                    gy = gb+0.5*gw+0*gsy;                  // grating position Y
                else if (X.y < gb+gw+1*gsy)
                    gy = gb+0.5*gw+1*gsy;                  // grating position Y
                else if (X.y < gb+gw+2*gsy)
                    gy = gb+0.5*gw+2*gsy;                  // grating position Y
                else if (X.y < gb+gw+3*gsy)
                    gy = gb+0.5*gw+3*gsy;                  // grating position Y
                else if (X.y < gb+gw+4*gsy)
                    gy = gb+0.5*gw+4*gsy;                  // grating position Y
                else if (X.y < gb+gw+5*gsy)
                    gy = gb+0.5*gw+5*gsy;                  // grating position Y
                else if (X.y < gb+gw+6*gsy)
                    gy = gb+0.5*gw+6*gsy;                  // grating position Y
                else if (X.y < gb+gw+7*gsy)
                    gy = gb+0.5*gw+7*gsy;                  // grating position Y
                else if (X.y < gb+gw+8*gsy)
                    gy = gb+0.5*gw+8*gsy;                  // grating position Y
                else if (X.y < gb+gw+9*gsy)
                    gy = gb+0.5*gw+9*gsy;                  // grating position Y
                else if (X.y < gb+gw+10*gsy)
                    gy = gb+0.5*gw+10*gsy;                 // grating position Y
                else if (X.y < gb+gw+11*gsy)
                    gy = gb+0.5*gw+11*gsy;                 // grating position Y
                else if (X.y < gb+gw+12*gsy)
                    gy = gb+0.5*gw+12*gsy;                 // grating position Y
                else if (X.y < gb+gw+13*gsy)
                    gy = gb+0.5*gw+13*gsy;                 // grating position Y
                else if (X.y < gb+gw+14*gsy)
                    gy = gb+0.5*gw+14*gsy;                 // grating position Y
                else if (X.y < gb+gw+15*gsy)
                    gy = gb+0.5*gw+15*gsy;                 // grating position Y
                else if (X.y < gb+gw+16*gsy)
                    gy = gb+0.5*gw+16*gsy;                 // grating position Y
                else if (X.y < gb+gw+17*gsy)
                    gy = gb+0.5*gw+17*gsy;                 // grating position Y
                else if (X.y < gb+gw+18*gsy)
                    gy = gb+0.5*gw+18*gsy;                 // grating position Y
                else if (X.y < gb+gw+19*gsy)
                    gy = gb+0.5*gw+19*gsy;                 // grating position Y
                else if (X.y < gb+gw+20*gsy)
                    gy = gb+0.5*gw+20*gsy;                 // grating position Y
                else if (X.y < gb+gw+21*gsy)
                    gy = gb+0.5*gw+21*gsy;                 // grating position Y
                else if (X.y < gb+gw+22*gsy)
                    gy = gb+0.5*gw+22*gsy;                 // grating position Y
                else if (X.y < gb+gw+23*gsy)
                    gy = gb+0.5*gw+23*gsy;                 // grating position Y
                else if (X.y < gb+gw+24*gsy)
                    gy = gb+0.5*gw+24*gsy;                 // grating position Y
                else if (X.y < gb+gw+25*gsy)
                    gy = gb+0.5*gw+25*gsy;                 // grating position Y
                else if (X.y < gb+gw+26*gsy)
                    gy = gb+0.5*gw+26*gsy;                 // grating position Y
                else if (X.y < gb+gw+27*gsy)
                    gy = gb+0.5*gw+27*gsy;                 // grating position Y
                else if (X.y < gb+gw+28*gsy)
                    gy = gb+0.5*gw+28*gsy;                 // grating position Y
                else if (X.y < gb+gw+29*gsy)
                    gy = gb+0.5*gw+29*gsy;                 // grating position Y
                else if (X.y < gb+gw+30*gsy)
                    gy = gb+0.5*gw+30*gsy;                 // grating position Y

                if (decimal == 2)
                    fprintf (mps_grid, "\n%d %9.2lf %9.2lf %9.2lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 3)
                    fprintf (mps_grid, "\n%d %9.3lf %9.3lf %9.3lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 4)
                    fprintf (mps_grid, "\n%d %9.4lf %9.4lf %9.4lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 5)
                    fprintf (mps_grid, "\n%d %9.5lf %9.5lf %9.5lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 6)
                    fprintf (mps_grid, "\n%d %9.6lf %9.6lf %9.6lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 7)
                    fprintf (mps_grid, "\n%d %9.7lf %9.7lf %9.7lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 8)
                    fprintf (mps_grid, "\n%d %9.8lf %9.8lf %9.8lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 9)
                    fprintf (mps_grid, "\n%d %9.9lf %9.9lf %9.9lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
                else if (decimal == 10)
                    fprintf (mps_grid, "\n%d %9.10lf %9.10lf %9.10lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0);
            }
        }
        printf ("MPS.grid: Done!\n");

        fclose (mps_grid);

        printf ("\a\a\a");
    }

    void print_MPS_2D (FILE *grid_2D, long cont, Particle particle[], double DP, int section, char d, int decimal)
    {
        long cont2D = 0;
        int dim = d-120;

        for (long i = 0; i < cont; i++)
            if (particle[i].r[dim] == section)
                cont2D++;

        fprintf (grid_2D, "%d\n%ld", 0, cont2D);

        for (long i = 0; i < cont; i++)
            if (particle[i].r[dim] == section)
            {
                Point X (particle[i].r , DP);
                if (decimal == 2)
                    fprintf (grid_2D, "\n%d %9.2lf %9.2lf %9.2lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 3)
                    fprintf (grid_2D, "\n%d %9.3lf %9.3lf %9.3lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 4)
                    fprintf (grid_2D, "\n%d %9.4lf %9.4lf %9.4lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 5)
                    fprintf (grid_2D, "\n%d %9.5lf %9.5lf %9.5lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 6)
                    fprintf (grid_2D, "\n%d %9.6lf %9.6lf %9.6lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 7)
                    fprintf (grid_2D, "\n%d %9.7lf %9.7lf %9.7lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 8)
                    fprintf (grid_2D, "\n%d %9.8lf %9.8lf %9.8lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 9)
                    fprintf (grid_2D, "\n%d %9.9lf %9.9lf %9.9lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
                else if (decimal == 10)
                    fprintf (grid_2D, "\n%d %9.10lf %9.10lf %9.10lf      0.0      0.0      0.0      0.0      0.0", particle[i].id, X.x, X.y, 0.0);
            }

        printf ("MPS_2D.grid: Done!\n");
    }

    void print_dMPS (FILE *dmps_grid, long cont, Particle particle[], double DP, int decimal)
    {
        double theta = 45*M_PI/180;                 // grating angle
        double gx = 0.90;                           // grating position X
        double gs = 0.010;                          // grating spacing
        double gb = 0.010;                          // space between tank bottom and first grating element
        double gw = 0.005;                          // grating width
        double gsy = (gs+gw)/cos(theta);            // grating spacing Y
        double gm = 0.013767767;                    // move granting downward Y
        double gy;


        int rmin[3];
        int rmax[3];
        int delta[3];
        double zero = 0.0;

        for (int i = 0; i < 3; i++)
        {
            rmin[i] = particle[0].r[i];
            rmax[i] = particle[0].r[i];
        }

        for (long p = 0; p < cont; p++)
            for (int i = 0; i < 3; i++)
            {
                if (rmin[i] > particle[p].r[i]) rmin[i] = particle[p].r[i];
                if (rmax[i] < particle[p].r[i]) rmax[i] = particle[p].r[i];
            }

        for (int i = 0; i < 3; i++)
            delta[i] = rmax[i] - rmin[i];

        char axis;

        if (delta[0] >= delta [1] && delta[0] >= delta[2]) axis = 'x';
        else if (delta[1] >= delta [2] && delta[1] >= delta[0]) axis = 'y';
        else if (delta[2] >= delta [0] && delta[2] >= delta[1]) axis = 'z';

        if (decimal == 2)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.2lf %10.2lf %10.2lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.2lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.2lf %10.2lf %10.2lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.2lf %10.2lf %10.2lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 3)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.3lf %10.3lf %10.3lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.3lf %10.3lf %10.3lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.3lf %10.3lf %10.3lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        
                    }
        }
        else if (decimal == 4)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.4lf %10.4lf %10.4lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.4lf %10.4lf %10.4lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.4lf %10.4lf %10.4lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 5)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.5lf %10.5lf %10.5lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.5lf %10.5lf %10.5lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.5lf %10.5lf %10.5lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 6)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.6lf %10.6lf %10.6lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.6lf %10.6lf %10.6lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.6lf %10.6lf %10.6lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 7)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.7lf %10.7lf %10.7lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.7lf %10.7lf %10.7lf %10.7lf %10.7lf %10.7lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.7lf %10.7lf %10.7lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.7lf %10.7lf %10.7lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 8)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.8lf %10.8lf %10.8lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.8lf %10.8lf %10.8lf %10.8lf %10.8lf %10.8lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.8lf %10.8lf %10.8lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.8lf %10.8lf %10.8lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 9)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.9lf %10.9lf %10.9lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.9lf %10.9lf %10.9lf %10.9lf %10.9lf %10.9lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.9lf %10.9lf %10.9lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.9lf %10.9lf %10.9lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }
        else if (decimal == 10)
        {
            fprintf(dmps_grid, "%9d %10d %10d %10d %10.10lf %10.10lf %10.10lf %10.1lf %10.1lf %10.1lf\n",
                0, (axis == 'x'), (axis == 'y'), (axis == 'z'), zero, zero, zero, zero, zero, zero);
            fprintf(dmps_grid, "%9ld %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.10lf %10.1lf %10.1lf %10.1lf\n", cont,
                 (rmin[0]+POSITIONX*0.5)*DP, (rmax[0]+POSITIONX*0.5)*DP,
                 (rmin[1]+POSITIONY*0.5)*DP, (rmax[1]+POSITIONY*0.5)*DP,
                 zero, zero, zero, zero, zero);

            for (long r = rmin[axis-120]; r <= rmax[axis-120]; r++)
                for (long i = 0; i < cont; i++)
                    if (particle[i].r[axis-120] == r)
                    {
                        Point X (particle[i].r , DP);
                        if (particle[i].id < 4)
                        {
                            fprintf (dmps_grid, "%9d %10.10lf %10.10lf %10.10lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                            particle[i].id, X.x, X.y, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                        else
                        {
                            if (X.y < gb+gw+0*gsy)
                                gy = gb+0.5*gw+0*gsy;                  // grating position Y
                            else if (X.y < gb+gw+1*gsy)
                                gy = gb+0.5*gw+1*gsy;                  // grating position Y
                            else if (X.y < gb+gw+2*gsy)
                                gy = gb+0.5*gw+2*gsy;                  // grating position Y
                            else if (X.y < gb+gw+3*gsy)
                                gy = gb+0.5*gw+3*gsy;                  // grating position Y
                            else if (X.y < gb+gw+4*gsy)
                                gy = gb+0.5*gw+4*gsy;                  // grating position Y
                            else if (X.y < gb+gw+5*gsy)
                                gy = gb+0.5*gw+5*gsy;                  // grating position Y
                            else if (X.y < gb+gw+6*gsy)
                                gy = gb+0.5*gw+6*gsy;                  // grating position Y
                            else if (X.y < gb+gw+7*gsy)
                                gy = gb+0.5*gw+7*gsy;                  // grating position Y
                            else if (X.y < gb+gw+8*gsy)
                                gy = gb+0.5*gw+8*gsy;                  // grating position Y
                            else if (X.y < gb+gw+9*gsy)
                                gy = gb+0.5*gw+9*gsy;                  // grating position Y
                            else if (X.y < gb+gw+10*gsy)
                                gy = gb+0.5*gw+10*gsy;                 // grating position Y
                            else if (X.y < gb+gw+11*gsy)
                                gy = gb+0.5*gw+11*gsy;                 // grating position Y
                            else if (X.y < gb+gw+12*gsy)
                                gy = gb+0.5*gw+12*gsy;                 // grating position Y
                            else if (X.y < gb+gw+13*gsy)
                                gy = gb+0.5*gw+13*gsy;                 // grating position Y
                            else if (X.y < gb+gw+14*gsy)
                                gy = gb+0.5*gw+14*gsy;                 // grating position Y
                            else if (X.y < gb+gw+15*gsy)
                                gy = gb+0.5*gw+15*gsy;                 // grating position Y
                            else if (X.y < gb+gw+16*gsy)
                                gy = gb+0.5*gw+16*gsy;                 // grating position Y
                            else if (X.y < gb+gw+17*gsy)
                                gy = gb+0.5*gw+17*gsy;                 // grating position Y
                            else if (X.y < gb+gw+18*gsy)
                                gy = gb+0.5*gw+18*gsy;                 // grating position Y
                            else if (X.y < gb+gw+19*gsy)
                                gy = gb+0.5*gw+19*gsy;                 // grating position Y
                            else if (X.y < gb+gw+20*gsy)
                                gy = gb+0.5*gw+20*gsy;                 // grating position Y
                            else if (X.y < gb+gw+21*gsy)
                                gy = gb+0.5*gw+21*gsy;                 // grating position Y
                            else if (X.y < gb+gw+22*gsy)
                                gy = gb+0.5*gw+22*gsy;                 // grating position Y
                            else if (X.y < gb+gw+23*gsy)
                                gy = gb+0.5*gw+23*gsy;                 // grating position Y
                            else if (X.y < gb+gw+24*gsy)
                                gy = gb+0.5*gw+24*gsy;                 // grating position Y
                            else if (X.y < gb+gw+25*gsy)
                                gy = gb+0.5*gw+25*gsy;                 // grating position Y
                            else if (X.y < gb+gw+26*gsy)
                                gy = gb+0.5*gw+26*gsy;                 // grating position Y
                            else if (X.y < gb+gw+27*gsy)
                                gy = gb+0.5*gw+27*gsy;                 // grating position Y
                            else if (X.y < gb+gw+28*gsy)
                                gy = gb+0.5*gw+28*gsy;                 // grating position Y
                            else if (X.y < gb+gw+29*gsy)
                                gy = gb+0.5*gw+29*gsy;                 // grating position Y
                            else if (X.y < gb+gw+30*gsy)
                                gy = gb+0.5*gw+30*gsy;                 // grating position Y

                            fprintf (dmps_grid, "%9d %10.10lf %10.10lf %10.10lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
                                particle[i].id, (X.x-gx)*cos(theta)+gx-(X.y-gy)*sin(theta), (X.y-gy)*cos(theta)+gy+(X.x-gx)*sin(theta)-gm, 0.0, zero, zero, zero, zero, zero, zero);
                        }
                    }
        }

        printf ("dMPS.grid: Done!\n\n");

        fclose(dmps_grid);

        printf ("\a\a\a");

    }

    void print_alt (FILE *alt, long cont, Particle particle[], double DP)
    {
        for (long i = 0; i < cont; i++)
        {
            Point X (particle[i].r , DP);
            if (particle[i].id == 2 && 
                X.x >= 0.5-0.5*DP && X.x < 0.5+0.5*DP &&
                X.y >= 0.0-1.0*DP && X.y < 0.0+0.0*DP)
                fprintf (alt, "\n<condition id=\"%ld\" y=\"0\" z=\"0\" x=\"0\" rz=\"0\" ry=\"0\" rx=\"0\" />", i );
        }

        printf ("alt.grid: Done!\n");

        fclose (alt);
    }

};

#endif

