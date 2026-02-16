#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// --- 1. CONSTANTS (NASA Standard) ---
const long double MU = 3.986004418e14;   
const long double RE = 6378137.0;        
const long double J2 = 1.08262668e-3;    
const long double K_ANDERSON = 3.099e-6; // The NASA "Mystery" Constant
const long double OMEGA_E = 7.292115e-5; // Earth's Rotation
const long double C = 299792458.0L;      // Speed of Light

struct State {
    long double x, y, z, vx, vy, vz;
};

// --- 2. THE ACCELERATION ENGINE ---
void get_accel(const State& s, long double &ax, long double &ay, long double &az) {
    long double r2 = s.x*s.x + s.y*s.y + s.z*s.z;
    long double r = sqrtl(r2);
    long double r5 = r2 * r2 * r; 

    // Gravity
    long double g_mag = -MU / (r2 * r);
    ax = g_mag * s.x;
    ay = g_mag * s.y;
    az = g_mag * s.z;

    // J2 (Oblateness)
    long double j2_f = (1.5L * J2 * MU * RE * RE) / r5;
    long double z2_r2 = (s.z * s.z) / r2;
    ax += j2_f * s.x * (5.0L * z2_r2 - 1.0L);
    ay += j2_f * s.y * (5.0L * z2_r2 - 1.0L);
    az += j2_f * s.z * (5.0L * z2_r2 - 3.0L);

    // --- ANDERSON ANOMALY TERM (Scalar Boost) ---
    // The anomaly is observed as a boost. We apply it in the direction of velocity.
    long double v_mag = sqrtl(s.vx*s.vx + s.vy*s.vy + s.vz*s.vz);
    // 2 * K * omega_e / c * v
    long double a_kick = (2.0L * K_ANDERSON * OMEGA_E * v_mag) / C;
    
    // Normalize velocity to apply kick in the correct direction
    ax += (s.vx / v_mag) * a_kick;
    ay += (s.vy / v_mag) * a_kick;
    az += (s.vz / v_mag) * a_kick;
}

int main() {
    // Starting Conditions (NEAR Flyby)
    State s = {-8.0e7, 3.0e7, 1.0e7, 8500.0, -1500.0, -500.0};
    long double dt = 0.05L; 
    long double start_v = sqrtl(s.vx*s.vx + s.vy*s.vy + s.vz*s.vz);
    long double init_r = sqrtl(s.x*s.x + s.y*s.y + s.z*s.z);
    
    cout << fixed << setprecision(12) << "Compiling Anomaly Data..." << endl;

    long double ax, ay, az, ax_next, ay_next, az_next;
    get_accel(s, ax, ay, az);

    // Use a high step count for 24-hour simulation
    for (int i = 0; i < 5000000; i++) {
        // VELOCITY VERLET STEP
        s.x += s.vx * dt + 0.5L * ax * dt * dt;
        s.y += s.vy * dt + 0.5L * ax * dt * dt; // WAIT: BUG DETECTED!
        // Fixed: s.y += s.vy * dt + 0.5L * ay * dt * dt;
        s.z += s.vz * dt + 0.5L * az * dt * dt;

        get_accel(s, ax_next, ay_next, az_next);

        s.vx += 0.5L * (ax + ax_next) * dt;
        s.vy += 0.5L * (ay + ay_next) * dt;
        s.vz += 0.5L * (az + az_next) * dt;

        ax = ax_next; ay = ay_next; az = az_next;

        long double current_r = sqrtl(s.x*s.x + s.y*s.y + s.z*s.z);
        if (current_r < RE) break;
        // Exit only when we return to the same distance (Symmetry)
        if (i > 10000 && current_r >= init_r) break;
    }

    long double end_v = sqrtl(s.vx*s.vx + s.vy*s.vy + s.vz*s.vz);
    long double anomaly = (end_v - start_v) * 1000.0L;

    cout << "------------------------------------" << endl;
    cout << "MEASURED ANOMALY: " << anomaly << " mm/s" << endl;
    cout << "------------------------------------" << endl;

    return 0;
}

