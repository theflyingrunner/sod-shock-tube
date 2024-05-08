// This program solves the Sod Shock Tube problem using the Steger-Warming scheme with a first-order extrapolation at the cell faces.
// To compile: gcc -Wall -o sod_shock_tube sod_shock_tube.c
// To run: ./sod_shock_tube 
// Author: Saron Bhoopathy
// Email: saron.bhoopathy@gmail.com

#include<math.h>
#include<stdlib.h>
#include<stdio.h>

// Global constants
const int nx = 1001;        // Number of grid points
const int nt = 3000;        // Number of time steps
const double g = 1.4;       // Specific heat ratio for air
const double nu = 0.1;      // Safety factor

// Function prototypes
void initialConditions(double uold[][3], double xmin, double xmax, double dx);
void computeFlux(double uold[][3], double unew[][3]);
void writeSolutionToFile(double uold[][3], double flux[][3], double xmin, double xmax);
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
    }
  }

  // Write solution to output file
  writeSolutionToFile(uold, flux, xmin, xmax);

  printf("Simulation completed. Solution written to 'solution.txt'.\n");

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

// Compute fluxes using Steger-Warming scheme
void computeFlux(double uold[][3], double flux[][3]) {
  double eps_fix = 0.001;
  double aa, uu, rr;
  double lambda1, lambda2, lambda3;
  double E1, E2, E3;

  // F+ flux
  for (int i = 0; i < nx; i++) {
    aa = sqrt(g * computeTemperature(uold[i]));
    uu = computeVelocity(uold[i]);
    rr = computeDensity(uold[i]);

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
    aa = sqrt(g * computeTemperature(uold[i + 1]));
    uu = computeVelocity(uold[i + 1]);
    rr = computeDensity(uold[i + 1]);

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
void writeSolutionToFile(double uold[][3], double flux[][3], double xmin, double xmax) {
  FILE *fp;
  fp = fopen("solution.txt", "w");

  double dx = (xmax - xmin) / (nx - 1);
  double x = xmin;

  for (int i = 0; i < nx; i++) {
    double p = computePressure(uold[i]);
    double rho = computeDensity(uold[i]);
    double T = computeTemperature(uold[i]);
    double u = computeVelocity(uold[i]);

    fprintf(fp, "%.6f %.6f %.6f %.6f %.6f\n", x, p, rho, T, u);
    x += dx;
  }

  fclose(fp);
}
