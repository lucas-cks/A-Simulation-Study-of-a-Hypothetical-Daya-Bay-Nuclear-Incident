//3D
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//Grid parameters
#define NX 200
#define NY 200
#define NZ 50
#define LX 100000.0   // 100 km
#define LY 100000.0   // 100 km
#define LZ 5000.0     // 5 km
#define DX (LX/(NX-1))
#define DY (LY/(NY-1))
#define DZ (LZ/(NZ-1))

//Physical parameters
#define z_speed 0.01     //Vertical velocity m/s
#define K_x 100.0       //Horizontal Diffusion Constant m^2/s
#define K_y 100.0       //Horizontal Diffusion Constant m^2/s
#define K_z 10.0        //Vertical Diffusion Constant m^2/s
#define lamda 1.0e-6    //Decay Constant s^-1
#define USTAR 0.3       // Friction velocity m/s
#define Z0 0.1          // Roughness length m
#define KAPPA 0.4       // von Kármán constant

//Time parameters
#define DT 10.0
#define NT 3600 //total setps (10 hours)

#define SOURCE_CONC 5.0 //Source concentration at the inflow boundary
#define DEPOSITION_RATE 0.001 // Deposition rate at the ground

double wind_profile(double z) {
    return (USTAR / KAPPA) * log((z + Z0) / Z0);
}

int main() {
    
    double *C = malloc(NX * NY * NZ * sizeof(double));
    double *C_new = malloc(NX * NY * NZ * sizeof(double));
    
    // initialise 
    for (int i = 0; i < NX * NY * NZ; i++) C[i] = 0.0;
    int plane = NX * NY; // number of points in one horizontal plane
    
    // stability check
    double max_u = wind_profile(LZ);  
    double cfl_x = max_u * DT / DX;
    double diff_stable_x = K_x * DT / (DX*DX);
    double diff_stable_y = K_y * DT / (DY*DY);
    double diff_stable_z = K_z * DT / (DZ*DZ);
    printf("Stability\n");
    printf("CFL_x = %.4f (<=1)\n", cfl_x);
    printf("Diff_x = %.4f (<=0.5)\n", diff_stable_x);
    printf("Diff_y = %.4f (<=0.5)\n", diff_stable_y);
    printf("Diff_z = %.4f (<=0.5)\n", diff_stable_z);
    if (cfl_x > 1.0 || diff_stable_x > 0.5 || diff_stable_y > 0.5 || diff_stable_z > 0.5)
        printf("WARNING: Unstable, Reduce DT.\n");
    else
        printf("Stable.\n");
        
    // output
    FILE *f = fopen("output_3d.bin", "wb");

    // Time loop
    for (int n = 0; n <= NT; n++) {
        // Output 
        if (n % 500 == 0 || n == NT) {
            fwrite(C, sizeof(double), NX * NY * NZ, f);
        }

        // Update concentration field
        for (int k = 1; k < NZ - 1; k++) {
            double z = k * DZ;
            double u_x = wind_profile(z);

            for (int j = 1; j < NY - 1; j++) {
                for (int i = 1; i < NX - 1; i++) {
                    int idx = i + j * NX + k * plane;

                    // x direction (advection + diffusion)
                    double advec_x = -u_x * (C[idx] - C[idx - 1]) / DX;
                    double diff_x = K_x * (C[idx + 1] - 2 * C[idx] + C[idx - 1]) / (DX * DX);

                    // y direction (diffusion only, assuming no lateral wind)
                    double diff_y = K_y * (C[idx + NX] - 2 * C[idx] + C[idx - NX]) / (DY * DY);

                    // z direction (vertical velocity + diffusion)
                    double advec_z = -z_speed * (C[idx] - C[idx - plane]) / DZ; 
                    double diff_z = K_z * (C[idx + plane] - 2 * C[idx] + C[idx - plane]) / (DZ * DZ);

                    double decay = -lamda * C[idx];

                    C_new[idx] = C[idx] + DT * (advec_x + diff_x + diff_y + advec_z + diff_z + decay);
                }
            }
        }

        // 3D BC

        // 1. Inflow (X=0) and point source
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NY; j++) {
                C_new[0 + j * NX + k * plane] = 0.0;
            }
        }
        // (X=0, Y=NY/2, Z=2km)
        int leak_y_idx = NY / 2;
        int leak_z_idx = 20; 
        C_new[0 + leak_y_idx * NX + leak_z_idx * plane] = SOURCE_CONC;

        // 2. Outflow (X=LX, Y=0, Y=LY) - Zero Gradient
        for (int k = 0; k < NZ; k++) {
            for (int j = 0; j < NY; j++) {
                C_new[(NX-1) + j * NX + k * plane] = C_new[(NX-2) + j * NX + k * plane];
            }
            for (int i = 0; i < NX; i++) {
                C_new[i + 0 * NX + k * plane] = C_new[i + 1 * NX + k * plane];
                C_new[i + (NY-1) * NX + k * plane] = C_new[i + (NY-2) * NX + k * plane];
            }
        }

        // 3. Top (Z=LZ)
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                C_new[i + j * NX + (NZ-1) * plane] = C_new[i + j * NX + (NZ-2) * plane];
            }
        }

        // 4. Ground (Z=0) - Deposition
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                int idx_bot = i + j * NX + 0 * plane;
                int idx_up  = i + j * NX + 1 * plane;
                C_new[idx_bot] = C_new[idx_up] * (1.0 - DEPOSITION_RATE);
            }
        }

        memcpy(C, C_new, NX * NY * NZ * sizeof(double));
    }
    
    fclose(f);
    free(C); free(C_new);
    return 0;
}
