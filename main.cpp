#include <iostream>

int main() {
long double x = -1.5e8;  // Starting position
long double vx = 8949.0; // Starting velocity
long double dt = 0.01;   // Time step
const long double G  = 6.67430e-11;    // Gravitational Constant
const long double MU = 3.986004418e14; // Earth's Gravitational Parameter
const long double R_E = 6378137.0;    // Earth's Radius in meters
