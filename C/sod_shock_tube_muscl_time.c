// This program solves the Sod Shock Tube problem using the Steger-Warming scheme and MUSCL formalism.
// Further, this version of the code stores the solution at every time-step that can later be used to visualize the evolution over time.
// To compile: gcc -Wall -o sod_shock_tube_muscl_time sod_shock_tube_muscl_time.c
// To run: ./sod_shock_tube_muscl_time
// Author: Saron Manikam Bhoopathy
// Email: saron.bhoopathy@gmail.com

#include<math.h>
#include<stdlib.h>
#include<stdio.h>

// Global constants
const int nx = 1001;        // Number of grid points
const int nt = 12000;        // Number of time steps
const double g = 1.4;       // Specific heat ratio for air
const double nu = 0.1;      // Safety factor

// Function prototypes
void initialConditions(double uold[][3], double xmin, double xmax, double dx);
void computeFlux(double uold[][3], double unew[][3]);
void computeMUSCL(double x[nx], double xr[nx], double xl[nx]);
void writeSolutionToFile(double **pressure, double **density, double **temperature, double **velocity, double xmin, double xmax);
double computeTimeStep(double uold[][3], double dx);
double computeDensity(double u[]);
double computeVelocity(double u[]);
double computeTemperature(double u[]);
double computePressure(double u[]);

// Main function
int main(int argc, char **argv) {

  // Initialize solution vectors and flux array
  double uold[nx][3];   // Solution vector at previous time step
  double unew[nx][3];   // Solution vector at current time step
  double flux[nx][3];   // Flux array

  // Initialize solution vectors with time-step
  //double pressureSolution[nx][nt];
  //double densitySolution[nx][nt];
  //double temperatureSolution[nx][nt];
  //double velocitySolution[nx][nt];

  double **pressureSolution = (double **)malloc(nx * sizeof(double *));
  for (int i = 0; i < nx; i++) {
    pressureSolution[i] = (double *)malloc(nt * sizeof(double));
  }

  double **densitySolution = (double **)malloc(nx * sizeof(double *));
  for (int i = 0; i < nx; i++) {
    densitySolution[i] = (double *)malloc(nt * sizeof(double));
  }

  double **temperatureSolution = (double **)malloc(nx * sizeof(double *));
  for (int i = 0; i < nx; i++) {
    temperatureSolution[i] = (double *)malloc(nt * sizeof(double));
  }

  double **velocitySolution = (double **)malloc(nx * sizeof(double *));
  for (int i = 0; i < nx; i++) {
    velocitySolution[i] = (double *)malloc(nt * sizeof(double));
  }

  // Initialize parameters and grid size
  double xmin = -5.0;
  double xmax = 5.0;
  double dx = (xmax - xmin) / (nx - 1.0);
  double dt;

  // Apply initial conditions
  initialConditions(uold, xmin, xmax, dx);

  // Loop through each step in time
  for (int i = 0; i < nt; i++) {
    // Compute time-step
    dt = computeTimeStep(uold, dx);
    // Compute flux
    computeFlux(uold, flux);
    // Loop through each step in space
    for (int j = 1; j < nx; j++) {
      // Loop through the three solution parameters: rho, momentum, and energy
      for (int k = 0; k < 3; k++) {
        unew[j][k] = uold[j][k] - (dt / dx) * (flux[j][k] - flux[j - 1][k]); 
        uold[j][k] = unew[j][k];
      }
    }
    
    // Apply reflective boundary conditions
    // Left boundary (j = 0)
    for (int k = 0; k < 3; k++) {
      unew[0][k] = uold[1][k]; // Reflect the solution from the first interior cell
    }

    // Right boundary (j = nx-1)
    for (int k = 0; k < 3; k++) {
      unew[nx-1][k] = uold[nx-2][k]; // Reflect the solution from the second-to-last interior cell
    }

    // Update uold with new solution
    for (int j = 0; j < nx; j++) {
      for (int k = 0; k < 3; k++) {
        uold[j][k] = unew[j][k];
      }
      pressureSolution[j][i] = computePressure(uold[j]);
      densitySolution[j][i] = computeDensity(uold[j]);
      temperatureSolution[j][i] = computeTemperature(uold[j]);
      velocitySolution[j][i] = computeVelocity(uold[j]);
    }
  }

  // Write solution to output file
  writeSolutionToFile(pressureSolution, densitySolution, temperatureSolution, velocitySolution, xmin, xmax);

  free(pressureSolution);
  free(densitySolution);
  free(temperatureSolution);
  free(velocitySolution);

  printf("Simulation completed.\n");

  return 0;
}

// Function to apply initial conditions
void initialConditions(double uold[][3], double xmin, double xmax, double dx) {
  double x[nx + 1];
  x[0] = xmin;

  double rho[nx];
  double p[nx];
  double u[nx];

  for (int i = 0; i < nx; i++) {
    if(x[i] <= 0) {
      rho[i] = 1;
      p[i] = 1;
      u[i] = 0;
    } else if (x[i] > 0) {
      rho[i] = 0.125;
      p[i] = 0.1;
      u[i] = 0;
    }
    uold[i][0] = rho[i];
    uold[i][1] = rho[i] * u[i];
    uold[i][2] = p[i] / (g - 1.0) + 0.5 * (rho[i] * pow(u[i], 2.0));
    x[i + 1] = x[i] + dx; 
  }
}

// Function to compute time step
double computeTimeStep(double uold[][3], double dx) {
  double max_u = 0.0;
  double max_T = 0.0;

  double p[nx];
  double T[nx];
  double u[nx];

  for (int i = 0; i < nx; i++) {
    p[i] = computePressure(uold[i]);
    T[i] = computeTemperature(uold[i]);
    u[i] = computeVelocity(uold[i]);
    double abs_u = fabs(u[i]);
    if (abs_u > max_u) {
      max_u = abs_u;
    }
    if (T[i] > max_T) {
      max_T = T[i];
    }
  }
  double dt = (nu * dx) / (max_u + sqrt(g * max_T));
  return dt;
}

// Function to compute density
double computeDensity(double u[]) {
  double rho = u[0];
  return rho;
}

// Function to compute velocity
double computeVelocity(double u[]) {
  double velocity = u[1] / u[0];
  return velocity;
}

// Function to compute temperature
double computeTemperature(double u[]) {
  double rho = computeDensity(u);
  double p = computePressure(u);
  double T = p / rho;
  return T;
}

// Function to compute pressure
double computePressure(double u[]) {
  double velocity = computeVelocity(u);
  double rho = computeDensity(u);
  double p = (g - 1.0) * (u[2] - 0.5 * rho * pow(velocity, 2.0));
  return p;
}

// Extrapolate to cell boundaries using the MUSCL scheme
void computeMUSCL(double x[nx], double xr[nx], double xl[nx]) {
  double kappa = 1.0 / 3.0;
  // double kappa = -1.0;
  double phiord = 1.0;
  double smalleps_extrap = 1.0e-8;
  double d1, d2, s;

  // Extrapolate at left-cell boundary
  for (int i = 1; i < nx; i++) {
    d1 = x[i] - x[i - 1];
    d2 = x[i + 1] - x[i];
    s = (2.0 * d2 * d1 + smalleps_extrap) / (pow(d2, 2.0) + pow(d1, 2.0) + smalleps_extrap);
    xl[i] = x[i] + phiord * 0.25 * s * ((1.0 - kappa * s) * d1 + (1.0 + kappa * s) * d2);
  }

  // Extrapolate at right-cell boundary
  for (int i = 0; i < nx; i++) {
    d1 = x[i + 1] - x[i];
    d2 = x[i + 2] - x[i + 1];
    s = (2.0 * d2 * d1 + smalleps_extrap) / (pow(d2, 2.0) + pow(d1, 2.0) + smalleps_extrap);
    xr[i] = x[i + 1] - phiord * 0.25 * s * ((1.0 + kappa * s) * d1 + (1.0 - kappa * s) * d2);
  }

  // Boundary conditions
  xl[0] = x[0];
  xr[nx - 2] = x[nx - 1];
}

// Compute fluxes using Steger-Warming scheme
void computeFlux(double uold[][3], double flux[][3]) {
  double eps_fix = 0.001;
  double aa, uu, rr;
  double lambda1, lambda2, lambda3;
  double E1, E2, E3;
  double Tl[nx], Tr[nx], ul[nx], ur[nx], rhol[nx], rhor[nx];
  double T[nx], u[nx], rho[nx];

  // Compute temperature, velocity, and density
  for (int i = 0; i < nx; i++) {
    T[i] = computeTemperature(uold[i]);
    u[i] = computeVelocity(uold[i]);
    rho[i] = computeDensity(uold[i]);
  }

  // Extrapolate temperature, velocity, and density at right- and left-cell boundaries using MUSCL scheme
  computeMUSCL(T, Tr, Tl);
  computeMUSCL(u, ur, ul);
  computeMUSCL(rho, rhor, rhol);

  // F+ flux
  for (int i = 0; i < nx; i++) {
    aa = sqrt(g * Tl[i]);
    uu = ul[i];
    rr = rhol[i];

    lambda1 = uu;
    lambda2 = uu + aa;
    lambda3 = uu - aa;

    lambda1 = 0.5 * (lambda1 + sqrt(pow(lambda1, 2.0) + pow(eps_fix, 2.0)));
    lambda2 = 0.5 * (lambda2 + sqrt(pow(lambda2, 2.0) + pow(eps_fix, 2.0)));
    lambda3 = 0.5 * (lambda3 + sqrt(pow(lambda3, 2.0) + pow(eps_fix, 2.0)));

    E1 = rr * (2.0 * (g - 1.0) * lambda1 + lambda2 + lambda3) / (2.0 * g);
    E2 = rr * (2.0 * (g - 1.0) * uu * lambda1 + (uu + aa) * lambda2 + (uu - aa) * lambda3) / (2.0 * g);
    E3 = rr * (2.0 * pow(uu, 2.0) * (g - 1.0) * (g - 1.0) * lambda1 +
    (2.0 * pow(aa, 2.0) + 2.0 * (g - 1.0) * aa * uu + (g - 1.0) * pow(uu, 2.0)) * lambda2 +
    (2.0 * pow(aa, 2.0) - 2.0 * (g - 1.0) * aa * uu + (g - 1.0) * pow(uu, 2.0)) * lambda3) /
    (4.0 * g * (g - 1.0));

    flux[i][0] = E1;
    flux[i][1] = E2;
    flux[i][2] = E3;
  }

  // F- flux
  for (int i = 0; i < nx; i++) {
    aa = sqrt(g * Tr[i]);
    uu = ur[i];
    rr = rhor[i];

    lambda1 = uu;
    lambda2 = uu + aa;
    lambda3 = uu - aa;

    lambda1 = 0.5 * (lambda1 - sqrt(pow(lambda1, 2.0) + pow(eps_fix, 2.0)));
    lambda2 = 0.5 * (lambda2 - sqrt(pow(lambda2, 2.0) + pow(eps_fix, 2.0)));
    lambda3 = 0.5 * (lambda3 - sqrt(pow(lambda3, 2.0) + pow(eps_fix, 2.0)));

    E1 = rr * (2.0 * (g - 1.0) * lambda1 + lambda2 + lambda3) / (2.0 * g);
    E2 = rr * (2.0 * (g - 1.0) * uu * lambda1 + (uu + aa) * lambda2 + (uu - aa) * lambda3) / (2.0 * g);
    E3 = rr * (2.0 * pow(uu, 2.0) * (g - 1.0) * (g - 1.0) * lambda1 +
    (2.0 * pow(aa, 2.0) + 2.0 * (g - 1.0) * aa * uu + (g - 1.0) * pow(uu, 2.0)) * lambda2 +
    (2.0 * pow(aa, 2.0) - 2.0 * (g - 1.0) * aa * uu + (g - 1.0) * pow(uu, 2.0)) * lambda3) /
    (4.0 * g * (g - 1.0));

    flux[i][0] = flux[i][0] + E1;
    flux[i][1] = flux[i][1] + E2;
    flux[i][2] = flux[i][2] + E3;
  }
}

// Function to write solution to file
void writeSolutionToFile(double **pressure, double **density, double **temperature, double **velocity, double xmin, double xmax) {
  FILE *fpDensity;
  FILE *fpVelocity;
  FILE *fpPressure;
  FILE *fpTemperature;

  fpDensity = fopen("density.txt", "w");
  fpVelocity = fopen("velocity.txt", "w");
  fpPressure = fopen("pressure.txt", "w");
  fpTemperature = fopen("temperature.txt", "w");

  for (int i = 0; i < nt; i++) {
    for (int j = 0; j < nx; j++) {
      fprintf(fpPressure, "%.6f\t", pressure[j][i]);
      fprintf(fpDensity, "%.6f\t", density[j][i]);
      fprintf(fpTemperature, "%.6f\t", temperature[j][i]);
      fprintf(fpVelocity, "%.6f\t", velocity[j][i]);
    }
    fprintf(fpPressure, "\n");
    fprintf(fpDensity, "\n");
    fprintf(fpTemperature, "\n");
    fprintf(fpVelocity, "\n");
  }

  fclose(fpPressure);
  fclose(fpDensity);
  fclose(fpTemperature);
  fclose(fpVelocity);
}
