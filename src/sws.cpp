//============================================================================
// Name        : sws.cpp
// Author      : Tom McCormack
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cmath>
#include <cstdlib>
#include <cstdio>

#define M_DIM 100

void AddArray(double (*A)[M_DIM], double (*B)[M_DIM], int x, int y, double t);
void ddy(double (*V)[M_DIM], double (*KV)[M_DIM], double (*MOD)[M_DIM], int length, double inc);
void ddx(double (*U)[M_DIM], double (*KU)[M_DIM], double (*MOD)[M_DIM], int length, double inc);
void rk4(double (*U)[M_DIM], double (*V)[M_DIM], double (*H)[M_DIM],
         double (*KU)[M_DIM], double (*KV)[M_DIM], double (*KH)[M_DIM], int x_l,
int y_l, double x_inc, double y_inc, double (*G)[M_DIM]);

int main() {

  // Assigning variables for operating space
  double x_min = 0;
  double x_max = 100;
  double x_inc = 1;
  double y_min = 0;
  double y_max = 100;
  double y_inc = 1;
  double t_min = 0;
  double t_max = 100;
  double t_inc = 0.01;

  // Calculating variable to determine the number of grid points in x, y, and t
  int x_grid_points = (x_max - x_min) / x_inc;
  int y_grid_points = (y_max - y_min) / y_inc;
  int t_grid_points = (t_max - t_min) / t_inc;

  // Creating arrays for operating space
  double x[x_grid_points];
  double y[y_grid_points];
  double t[t_grid_points];

  // Assigning dimensions values to the x, y, and t matrices
  for (int i = 0; i < x_grid_points; i++) {
    for (int j = 0; j < y_grid_points; j++) {
      x[i] = i * x_inc;
      y[j] = j * y_inc;
    }
  }

  // Creating Matrices to store values of u (x-velocity), v (y-velocity), and h (height), for each grid point of x and y
  double u[x_grid_points][y_grid_points];
  double v[x_grid_points][y_grid_points];
  double h[x_grid_points][y_grid_points];
  double g[x_grid_points][y_grid_points];

  // Assigning initial values to the u, v, and h matrices
  for (int i = 0; i < x_grid_points; i++) {
    for (int j = 0; j < y_grid_points; j++) {
      u[i][j] = 0;
      v[i][j] = 0;
      h[i][j] = 1 + 0.5 * exp(-(pow(x[i] - 30, 2) + pow(y[j] - 30, 2)) / 25);
      g[i][j] = 9.81;
    }
  }

  // Matrices for the storage of modifiers to the u, h, and v matrices, to ensure they don't get modified before their time
  double KU1[x_grid_points][y_grid_points];
  double KV1[x_grid_points][y_grid_points];
  double KH1[x_grid_points][y_grid_points];
  double KU2[x_grid_points][y_grid_points];
  double KV2[x_grid_points][y_grid_points];
  double KH2[x_grid_points][y_grid_points];
  double KU3[x_grid_points][y_grid_points];
  double KV3[x_grid_points][y_grid_points];
  double KH3[x_grid_points][y_grid_points];
  double KU4[x_grid_points][y_grid_points];
  double KV4[x_grid_points][y_grid_points];
  double KH4[x_grid_points][y_grid_points];

  FILE * pFile;

  pFile = fopen("SWS.csv", "w+");

  // Time stepping loop begins here, where h, u, and v shall be over-written for each time step to save space

  for (int i = 0; i < t_grid_points; i++) {

    for (int i = 0; i < x_grid_points; i++) {
      for (int j = 0; j < y_grid_points; j++) {
        KU1[i][j] = 0;
        KV1[i][j] = 0;
        KH1[i][j] = 0;
        KU2[i][j] = 0;
        KV2[i][j] = 0;
        KH2[i][j] = 0;
        KU3[i][j] = 0;
        KV3[i][j] = 0;
        KH3[i][j] = 0;
        KU4[i][j] = 0;
        KV4[i][j] = 0;
        KH4[i][j] = 0;
      }
    }

    rk4(u, v, h, KU1, KV1, KH1, x_grid_points, y_grid_points, t_inc, g);

    rk4(AddArray(u, KU1, x_grid_points, y_grid_points, t_inc / 2),
        AddArray(v, KV1, x_grid_points, y_grid_points, t_inc / 2),
        AddArray(h, KH1, x_grid_points, y_grid_points, t_inc / 2), KU2, KV2,
        KH2, x_grid_points, y_grid_points, t_inc, g);

    rk4(AddArray(u, KU2, x_grid_points, y_grid_points, t_inc / 2),
        AddArray(v, KV2, x_grid_points, y_grid_points, t_inc / 2),
        AddArray(h, KH2, x_grid_points, y_grid_points, t_inc / 2), KU3, KV3,
        KH3, x_grid_points, y_grid_points, t_inc, g);

    rk4(AddArray(u, KU3, x_grid_points, y_grid_points, t_inc),
        AddArray(v, KV3, x_grid_points, y_grid_points, t_inc),
        AddArray(h, KH3, x_grid_points, y_grid_points, t_inc), KU4, KV4, KH4,
        x_grid_points, y_grid_points, t_inc, g);

    //  Updating u for one time step
    AddArray(u, KU1, x_grid_points, y_grid_points, t_inc / 6);

    AddArray(u, KU2, x_grid_points, y_grid_points, t_inc / 3);

    AddArray(u, KU3, x_grid_points, y_grid_points, t_inc / 3);

    AddArray(u, KU4, x_grid_points, y_grid_points, t_inc / 6);

    // Updating V for one time step
    AddArray(v, KV1, x_grid_points, y_grid_points, t_inc / 6);

    AddArray(v, KV2, x_grid_points, y_grid_points, t_inc / 3);

    AddArray(v, KV3, x_grid_points, y_grid_points, t_inc / 3);

    AddArray(v, KV4, x_grid_points, y_grid_points, t_inc / 6);

    // Updating h for one time step
    AddArray(h, KH1, x_grid_points, y_grid_points, t_inc / 6);

    AddArray(h, KH2, x_grid_points, y_grid_points, t_inc / 3);

    AddArray(h, KH3, x_grid_points, y_grid_points, t_inc / 3);

    AddArray(h, KH4, x_grid_points, y_grid_points, t_inc / 6);

    for (int i = 0; i < x_grid_points; i++) {
      for (int j = 0; j < y_grid_points; j++) {
        fprintf(pFile, "%f ,", h[i][j]);
      }
    }
  fprintf(pFile, "\n %d n\", i);

    }
  }

  void rk4(double (*U)[M_DIM], double (*V)[M_DIM], double (*H)[M_DIM],
      double (*KU)[M_DIM], double (*KV)[M_DIM], double (*KH)[M_DIM], int x_l,
      int y_l, double x_inc, double y_inc, double (*G)[M_DIM]) {

    double k_temp[x_l][y_l];
    for (int i = 0; i < x_l; i++) {
      for (int j = 0; j < y_l; j++) {
        k_temp[i][j] = 0;
      }
    }

// Creating a temporary matrix to store values of various spatial derivative, and the adding them to the KU matrix

    ddx(U, k_temp, U);

    AddArray(KU, k_temp, x_l, y_l, -1);

    ddy(U, k_temp, V);

    AddArray(KU, k_temp, x_l, y_l, -1);

    ddx(H, k_temp, G);

    AddArray(KU, k_temp, x_l, y_l, -1);

// Creating a temporary matrix to store values of various spatial derivative, and the adding them to the KV matrix

    ddx(V, k_temp, U);

    AddArray(KV, k_temp, x_l, y_l, -1);

    ddy(V, k_temp, V);

    AddArray(KV, k_temp, x_l, y_l, -1);

    ddy(H, k_temp, G);

    AddArray(KV, k_temp, x_l, y_l, -1);

// Creating a temporary matrix to store values of various spatial derivative, and the adding them to the KH matrix

    ddx(H, k_temp, U);

    AddArray(KH, k_temp, x_l, y_l, -1);

    ddy(H, k_temp, V);

    AddArray(KH, k_temp, x_l, y_l, -1);

    ddx(U, k_temp, H);

    AddArray(KH, k_temp, x_l, y_l, -1);

    ddy(V, k_temp, H);

    AddArray(KH, k_temp, x_l, y_l, -1);

  }

  void ddx(double (*U)[M_DIM], double (*KU)[M_DIM], double (*MOD)[M_DIM],
      int length, double inc) {

    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        if (i == 1) {
          KU(i, j) = (MOD(i, j) / (12 * inc))
              * (U(length - 1, j) - 8 * U(length, j) + 8 * U(i + 1, j)
                  - U(i + 2, j));
        } else if (i == 2) {
          KU(i, j) = (MOD(i, j) / (12 * inc))
              * (U(length, j) - 8 * U(1, j) + 8 * U(3, j) - U(4, j));
        } else if (i == length - 1) {
          KU(i, j) = (MOD(i, j) / (12 * inc))
              * (U(length - 3, j) - 8 * U(length - 2, j) + 8 * U(length, j)
                  - U(1, j));
        } else if (i == length) {
          KU(i, j) =
              (MOD(i, j) / (12 * inc))
                  * (U(length - 2, j) - 8 * U(length - 1, j) + 8 * U(1, j)
                      - U(2, j));
        } else {
          KU(i, j) = (MOD(i, j) / (12 * inc))
              * (U(i - 2, j) - 8 * U(i - 1, j) + 8 * U(i + 1, j) - U(i + 2, j));
        }
      }
    }
  }

  void ddy(double (*V)[M_DIM], double (*KV)[M_DIM], double (*MOD)[M_DIM],
      int length, double inc) {

    for (int i = 0; i < length; i++) {
      for (int j = 0; j < length; j++) {
        if (j == 1) {
          KV(i, j) = (MOD(i, j) / (12 * inc))
              * (V(i, length - 1) - 8 * V(i, length) + 8 * V(i, 2) - V(i, 3));
        } else if (j == 2) {
          KV(i, j) = (MOD(i, j) / (12 * inc))
              * (V(i, length) - 8 * V(i, 1) + 8 * V(i, 3) - V(i, 3));
        } else if (j == length - 1) {
          KV(i, j) = (MOD(i, j) / (12 * inc))
              * (V(i, length - 3) - 8 * V(i, length - 2) + 8 * V(i, length)
                  - V(i, 1));
        } else if (j == length) {
          KV(i, j) =
              (MOD(i, j) / (12 * inc))
                  * (V(i, length - 2) - 8 * V(i, length - 1) + 8 * V(i, 1)
                      - V(i, 2));
        } else {
          KV(i, j) = (MOD(i, j) / (12 * inc))
              * (V(i, j - 2) - 8 * V(i, j - 1) + 8 * V(i, j + 1) - V(i, j + 2));
        }
      }
    }
  }

  void AddArray(double (*A)[M_DIM], double (*B)[M_DIM], int x, int y,
      double t) {

    for (int i = 0; i < x; i++) {
      for (int j = 0; j < y; j++) {
        A[i][j] += t * B[i][j];
      }
    }
  }
