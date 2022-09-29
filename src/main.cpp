#include "define.h"
#include "Point.h"
#include "Particle.h"

int main() {

    time_t start, end;

    //createFiles(paraview_xyz, mps_grid, dmps_grid, input, temp);

    FILE *paraview_xyz, *paraview_2D_xyz, *mps_grid, *mps_2D_grid, *dmps_grid, *input, *temp;
    int mode;
    long section;
    char d;
    char file_name[128];
    char file_name_temp[128];

    {
		//printf ("Opening default input: burner.txt ");
        printf ("Enter the name of the case: ");
        scanf ("%s", file_name);

        if (strcmp (file_name,"def") == 0)
            sprintf (file_name, "data/burner");

        printf("\n\n");

        time(&start);

        strcpy(file_name_temp, file_name);
        strcat(file_name_temp, ".txt");
        input = fopen(file_name_temp, "r");
        if (input == NULL) {
            printf("The file %s does not exist or an error occurred while opening.\n\n\n", file_name_temp);
            system("PAUSE");
            return 0;
        }
        printf("Initialization of file %s: Done!\n\n", file_name_temp);
        fscanf(input, "%d", &mode);
        printf("mode = %d\n\n", mode);

        if (mode % 2 == 0) {

            //fscanf (input, "\n%c", &d);
            d = 'z';
            printf("The 2D cases are have to be sketched in XY\n", d);
            //fscanf (input, "%ld", &section);
            section = 0;
            //printf ("section = %ld\n\n", section);

            //printf ("2D case in section %ld of the axis %c.\n\n", section, d);

            strcpy(file_name_temp, file_name);
            strcat(file_name_temp, ".xyz");
            paraview_2D_xyz = fopen(file_name_temp, "w");
            printf("Initialization of file %s: Done!\n", file_name_temp);

            strcpy(file_name_temp, file_name);
            strcat(file_name_temp, "_mps.grid");
            mps_2D_grid = fopen(file_name_temp, "w");
            printf("Initialization of file %s: Done!\n", file_name_temp);

        }


        if (mode % 3 == 0) {

            strcpy(file_name_temp, file_name);
            strcat(file_name_temp, ".xyz");
            paraview_xyz = fopen(file_name_temp, "w");
            printf("Initialization of file %s: Done!\n", file_name_temp);

            strcpy(file_name_temp, file_name);
            strcat(file_name_temp, "_mps.grid");
            mps_grid = fopen(file_name_temp, "w");
            printf("Initialization of file %s: Done!\n", file_name_temp);

            if (DMPS == 1) {
                strcpy(file_name_temp, file_name);
                strcat(file_name_temp, "_dmps.grid");
                dmps_grid = fopen(file_name_temp, "w");
                printf("Initialization of file %s: Done!\n\n", file_name_temp);
            }

        }

        strcpy(file_name_temp, file_name);
        strcat(file_name_temp, "_temp.xyz");
        temp = fopen(file_name_temp, "w");
        printf("Initialization of file %s: Done!\n\n", file_name_temp);

        printf("Files initialization: Done!\n\n");

        //createFiles(paraview_xyz, mps_grid, dmps_grid, input, temp);

    }

    //readGeometryFile ();

    double DP;
    int numGeometry;
    double geometry[NUMMAXGEOMETRY];
    int numFluidID;
    int fluidID[NUMMAXID];
    int numSolidID;
    int solidID[NUMMAXID];
    double r[3];
    Particle zeroPoint;

    {
        fscanf(input, "%lf", &DP);
        printf("\nDP = %lf m\n", DP);

        fscanf(input, "%d", &numGeometry);
        printf("\nnumGeometry = %d\n", numGeometry);
        for (int n = 0; n < numGeometry; n++) {
            fscanf(input, "%lf", &geometry[n]);
            printf("geometry[%d] = %lf\n", n + 1, geometry[n]);
        }

        fscanf(input, "%d", &numFluidID);
        printf("\nnumFluidID = %d\n", numFluidID);
        for (int n = 0; n < numFluidID; n++) {
            fscanf(input, "%d", &fluidID[n]);
            printf("fluidID[%d] = %d\n", n + 1, fluidID[n]);
        }

        fscanf(input, "%d", &numSolidID);
        printf("\nnumSolidID = %d\n", numSolidID);
        for (int n = 0; n < numSolidID; n++) {
            fscanf(input, "%d", &solidID[n]);
            printf("solidID[%d] = %d\n", n + 1, solidID[n]);
        }

        for (int n = 0; n < 3; n++) {
            fscanf(input, "%lf", &r[n]);
            zeroPoint.r[n] = floor(r[n] / DP) - CELLRADIUS;
        }

        printf("\nzeroPoint = ( %lf , %lf , %lf )\n", r[0], r[1], r[2]);

    }

    int decimal = 0;
    double num = DP/2;
    while (num >= 0.0000001)
    {
        num = num * 10;
        decimal = decimal + 1;
        num = num - int(num);
    }

    printf ("\nMinimum decimal places: %d ", decimal);
    printf ("\nEnter with decimal places (between 2 and 10): ");
    scanf ("%d", &decimal);

    if (decimal < 2)
        decimal = 2;
    else if (decimal > 10)
        decimal = 10;

    //initializationParticles ();

    Particle *particle = new Particle[NUMMAXPARTICLES];

    //initializationParticles ();

    long contParticles = zeroPoint.burnCells(temp, particle, DP, geometry, fluidID, numFluidID, solidID, numSolidID);

    Particle print;
    if (DIMENSION == 2) {
        print.print_Paraview_2D(paraview_xyz, contParticles, particle, DP, geometry, numGeometry, decimal);
        print.print_MPS_2D(mps_grid, contParticles, particle, DP, decimal);
        if (DMPS == 1) 
            print.print_dMPS_2D(dmps_grid, contParticles, particle, DP, decimal);
    }

    if (DIMENSION == 3) {
        print.print_Paraview_3D(paraview_xyz, contParticles, particle, DP, geometry, numGeometry, decimal);
        print.print_MPS_3D(mps_grid, contParticles, particle, DP, decimal);
        if (DMPS == 1) 
            print.print_dMPS_3D(dmps_grid, contParticles, particle, DP, decimal);
    }

    if (ALT == 1) { //condition, random position, random velocity
        strcat(file_name, "_alt.grid");
        FILE *alt = fopen(file_name, "w");
        printf("\nInitialization of file %s: Done!\n", file_name);
        print.print_alt(alt, contParticles, particle, DP);
    }

    time(&end);
    printf("\ntotal time = %.0lf s\n", difftime(end, start));

    fclose(temp);
    remove(file_name_temp);

    return 0;

}
