#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

// --- 1. GLOBAL CONSTANTS (NASA Standards) ---
const long double MU = 3.986004418e14;   // Earth's Gravitational Parameter
const long double RE = 6378137.0;        // Earth's Equatorial Radius (meters)
const long double J2 = 1.08262668e-3;    // Earth's Oblateness Coefficient

// --- 2. THE PHYSICS ENGINE (State Derivative) ---
// This calculates acceleration including J2 Perturbation
void get_acceleration(long double x, long double y, long double z, long double &ax, long double &ay, long double &az) {
    long double r_sq = x*x + y*y + z*z;
    long double r = sqrtl(r_sq);
    
    // Standard Newtonian Gravity
    long double mag = -MU / (r_sq * r);
    ax = mag * x;
    ay = mag * y;
    az = mag * z;

    // J2 Perturbation (The "Bulge" correction)
    long double factor = (1.5 * J2 * MU * RE * RE) / powl(r, 5);
    long double z_sq_r_sq = (z * z) / r_sq;

    ax += factor * x * (5.0 * z_sq_r_sq - 1.0);
    ay += factor * y * (5.0 * z_sq_r_sq - 1.0);
    az += factor * z * (5.0 * z_sq_r_sq - 3.0);
}

// --- 3. RUNGE-KUTTA 4th ORDER INTEGRATOR ---
// This moves the ship forward in time with extreme accuracy
void rk4_step(long double &x, long double &y, long double &z, 
              long double &vx, long double &vy, long double &vz, long double dt) {
    long double k1v[3], k2v[3], k3v[3], k4v[3];
    long double k1r[3], k2r[3], k3r[3], k4r[3];

    // k1
    get_acceleration(x, y, z, k1v[0], k1v[1], k1v[2]);
    k1r[0] = vx; k1r[1] = vy; k1r[2] = vz;

    // k2
    get_acceleration(x + 0.5*dt*k1r[0], y + 0.5*dt*k1r[1], z + 0.5*dt*k1r[2], k2v[0], k2v[1], k2v[2]);
    k2r[0] = vx + 0.5*dt*k1v[0]; k2r[1] = vy + 0.5*dt*k1v[1]; k2r[2] = vz + 0.5*dt*k1v[2];

    // k3
    get_acceleration(x + 0.5*dt*k2r[0], y + 0.5*dt*k2r[1], z + 0.5*dt*k2r[2], k3v[0], k3v[1], k3v[2]);
    k3r[0] = vx + 0.5*dt*k2v[0]; k3r[1] = vy + 0.5*dt*k2v[1]; k3r[2] = vz + 0.5*dt*k2v[2];

    // k4
    get_acceleration(x + dt*k3r[0], y + dt*k3r[1], z + dt*k3r[2], k4v[0], k4v[1], k4v[2]);
    k4r[0] = vx + dt*k3v[0]; k4r[1] = vy + dt*k3v[1]; k4r[2] = vz + dt*k3v[2];

    // Update state using weighted average
    x += (dt/6.0) * (k1r[0] + 2*k2r[0] + 2*k3r[0] + k4r[0]);
    y += (dt/6.0) * (k1r[1] + 2*k2r[1] + 2*k3r[1] + k4r[1]);
    z += (dt/6.0) * (k1r[2] + 2*k2r[2] + 2*k3r[2] + k4r[2]);

    vx += (dt/6.0) * (k1v[0] + 2*k2v[0] + 2*k3v[0] + k4v[0]);
    vy += (dt/6.0) * (k1v[1] + 2*k2v[1] + 2*k3v[1] + k4v[1]);
    vz += (dt/6.0) * (k1v[2] + 2*k2v[2] + 2*k3v[2] + k4v[2]);
}

int main() {
   
    long double x = -8.0e7, y = 3.0e7, z = 1.0e7; 
    long double vx = 8000.0, vy = -2000.0, vz = -1000.0;
    long double dt = 0.1;

    cout << fixed << setprecision(12);
    cout << "Step, Time, Velocity(m/s), Distance(km)" << endl;

    for (int i = 0; i < 1000000; i++) {
        rk4_step(x, y, z, vx, vy, vz, dt);

        if (i % 10000 == 0) {
            long double v_mag = sqrtl(vx*vx + vy*vy + vz*vz);
            long double r_mag = sqrtl(x*x + y*y + z*z) / 1000.0;
            cout << i << ", " << i*dt << ", " << v_mag << ", " << r_mag << endl;
        }
    }
    return 0;
}
;
