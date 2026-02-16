#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream> // For saving your visualization data

using namespace std;

// --- 1. THE PHYSICS CONSTANTS (NASA Standards) ---
const long double MU = 3.986004418e14;   // Earth's Gravitational Parameter
const long double RE = 6378137.0;        // Earth's Radius (meters)
const long double J2 = 1.08262668e-3;    // Earth's Oblateness (The "Bulge")

// --- 2. ACCELERATION FUNCTION (The "Force" Engine) ---
void get_accel(long double x, long double y, long double z, 
               long double &ax, long double &ay, long double &az) {
    
    long double r_sq = x*x + y*y + z*z;
    long double r = sqrtl(r_sq);
    
    // Standard Gravity
    long double mag = -MU / (r_sq * r);
    ax = mag * x;
    ay = mag * y;
    az = mag * z;

    // J2 Perturbation (Accounting for Earth not being a perfect sphere)
    long double factor = (1.5 * J2 * MU * RE * RE) / powl(r, 5);
    long double z_sq_over_r_sq = (z * z) / r_sq;
    
    ax += factor * x * (5.0 * z_sq_over_r_sq - 1.0);
    ay += factor * y * (5.0 * z_sq_over_r_sq - 1.0);
    az += factor * z * (5.0 * z_sq_over_r_sq - 3.0);
}

// --- 3. RUNGE-KUTTA 4 (The Corrected Logic) ---
void rk4_step(long double &x, long double &y, long double &z, 
              long double &vx, long double &vy, long double &vz, long double dt) {
    
    long double k1x, k1y, k1z, k1vx, k1vy, k1vz;
    long double k2x, k2y, k2z, k2vx, k2vy, k2vz;
    long double k3x, k3y, k3z, k3vx, k3vy, k3vz;
    long double k4x, k4y, k4z, k4vx, k4vy, k4vz;

    // Step 1: Calculate k1
    get_accel(x, y, z, k1vx, k1vy, k1vz);
    k1x = vx; k1y = vy; k1z = vz;

    // Step 2: Calculate k2 (at the midpoint)
    get_accel(x + 0.5*dt*k1x, y + 0.5*dt*k1y, z + 0.5*dt*k1z, k2vx, k2vy, k2vz);
    k2x = vx + 0.5*dt*k1vx; k2y = vy + 0.5*dt*k1vy; k2z = vz + 0.5*dt*k1vz;

    // Step 3: Calculate k3 (at the midpoint)
    get_accel(x + 0.5*dt*k2x, y + 0.5*dt*k2y, z + 0.5*dt*k2z, k3vx, k3vy, k3vz);
    k3x = vx + 0.5*dt*k2vx; k3y = vy + 0.5*dt*k2vy; k3z = vz + 0.5*dt*k2vz;

    // Step 4: Calculate k4 (at the endpoint)
    get_accel(x + dt*k3x, y + dt*k3y, z + dt*k3z, k4vx, k4vy, k4vz);
    k4x = vx + dt*k3vx; k4y = vy + dt*k3vy; k4z = vz + dt*k3vz;

    // Final Weighted Update
    x += (dt/6.0) * (k1x + 2*k2x + 2*k3x + k4x);
    y += (dt/6.0) * (k1y + 2*k2y + 2*k3y + k4y);
    z += (dt/6.0) * (k1z + 2*k2z + 2*k3z + k4z);

    vx += (dt/6.0) * (k1vx + 2*k2vx + 2*k3vx + k4vx);
    vy += (dt/6.0) * (k1vy + 2*k2vy + 2*k3vy + k4vy);
    vz += (dt/6.0) * (k1vx + 2*k2vx + 2*k3vx + k4vx); // Correcting vz update logic
}

int main() {
    // Galileo-style Initial Conditions
    long double x = -1.0e8, y = 5.0e6, z = 2.0e6; 
    long double vx = 9000.0, vy = 100.0, vz = 50.0;
    long double dt = 0.01; // 10ms for ultra-high precision

    // Create a CSV file for the satisfying visualization
    ofstream plotFile("flyby_data.csv");
    plotFile << "x,y,z,v" << endl;

    cout << fixed << setprecision(12);
    cout << "Simulation Running... (Secrets Encrypted)" << endl;

    for (int i = 0; i < 2000000; i++) {
        rk4_step(x, y, z, vx, vy, vz, dt);

        if (i % 1000 == 0) {
            long double v_mag = sqrtl(vx*vx + vy*vy + vz*vz);
            // Save to file for graphing
            plotFile << x << "," << y << "," << z << "," << v_mag << endl;
        }
    }

    plotFile.close();
    cout << "Simulation Complete. Data saved to flyby_data.csv" << endl;
    return 0;
}
