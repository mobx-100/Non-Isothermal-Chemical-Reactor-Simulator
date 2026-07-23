#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

// ==========================================
// CLASS 1: THE REACTOR MODEL
// ==========================================
class CSTRModel {
private:
    // Reactor Parameters
    double F;       // Volumetric flow rate (m^3/s)
    double V;       // Reactor volume (m^3)
    double C_A0;    // Inlet concentration of A (kmol/m^3)
    double T0;      // Inlet temperature (K)

    // Kinetic & Thermodynamic Parameters
    double k0;      // Pre-exponential factor (1/s)
    double E_a;     // Activation energy (J/mol)
    double R;       // Universal gas constant (J/(mol K))
    double dH;      // Heat of reaction (J/kmol) - negative for exothermic
    double rho;     // Density of mixture (kg/m^3)
    double Cp;      // Heat capacity (J/(kg K))

    // Cooling Jacket Parameters
    double UA;      // Overall heat transfer coeff * Area (W/K)
    double Tc;      // Cooling jacket temperature (K)

public:
    // Constructor to initialize parameters
    CSTRModel(double f, double v, double ca0, double t0, double k_0, double ea,
              double dh, double r, double c_p, double ua, double tc) {
        F = f; V = v; C_A0 = ca0; T0 = t0; k0 = k_0; E_a = ea; R = 8.314;
        dH = dh; rho = r; Cp = c_p; UA = ua; Tc = tc;
              }

              // Method to calculate derivatives (dCa/dt and dT/dt)
              vector<double> getDerivatives(double C_A, double T) {
                  // Calculate rate constant k using Arrhenius equation
                  double k = k0 * exp(-E_a / (R * T));

                  // Mass balance ODE: dCa/dt
                  double dCadt = (F / V) * (C_A0 - C_A) - k * C_A;

                  // Energy balance ODE: dT/dt
                  double heat_generated = (-dH / (rho * Cp)) * k * C_A;
                  double heat_removed = (UA / (rho * Cp * V)) * (T - Tc);
                  double sensible_heat = (F / V) * (T0 - T);

                  double dTdt = sensible_heat + heat_generated - heat_removed;

                  return {dCadt, dTdt}; // Return as vector {dCa/dt, dT/dt}
              }
};

// ==========================================
// CLASS 2: THE RK4 NUMERICAL SOLVER
// ==========================================
class RK4Solver {
public:
    void solve(CSTRModel& reactor, double initial_Ca, double initial_T, double dt, double t_final, string filename) {
        ofstream outFile(filename);
        if (!outFile.is_open()) {
            cerr << "Error: Could not create CSV file!" << endl;
            return;
        }

        // Write CSV Headers
        outFile << "Time(s),Concentration(kmol/m3),Temperature(K)\n";

        double t = 0.0;
        double Ca = initial_Ca;
        double T = initial_T;

        // Simulation Loop using Runge-Kutta 4th Order
        while (t <= t_final) {
            // Log current state to CSV
            outFile << fixed << setprecision(4) << t << "," << Ca << "," << T << "\n";

            // RK4 Step 1
            vector<double> k1 = reactor.getDerivatives(Ca, T);

            // RK4 Step 2
            vector<double> k2 = reactor.getDerivatives(Ca + 0.5 * dt * k1[0], T + 0.5 * dt * k1[1]);

            // RK4 Step 3
            vector<double> k3 = reactor.getDerivatives(Ca + 0.5 * dt * k2[0], T + 0.5 * dt * k2[1]);

            // RK4 Step 4
            vector<double> k4 = reactor.getDerivatives(Ca + dt * k3[0], T + dt * k3[1]);

            // Update state variables
            Ca = Ca + (dt / 6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
            T = T + (dt / 6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);

            t += dt; // Advance time
        }

        outFile.close();
        cout << "Simulation complete. Data exported to " << filename << endl;
    }
};

// ==========================================
// MAIN FUNCTION
// ==========================================
int main() {
    // 1. Define typical chemical reactor parameters
    double F = 0.1;           // Flow rate (m^3/s)
    double V = 2.0;           // Volume (m^3)
    double C_A0 = 1.0;        // Inlet concentration (kmol/m^3)
    double T0 = 350.0;        // Inlet temp (K)
    double k0 = 7.2e10;       // Arrhenius pre-exponential (1/s)
    double E_a = 7.275e4;     // Activation energy (J/mol)
    double dH = -5.0e7;       // Heat of reaction (J/kmol) - Exothermic
    double rho = 1000.0;      // Density (kg/m^3)
    double Cp = 4184.0;       // Heat capacity (J/kg/K) (Water-like)
    double UA = 5.0e4;        // Heat transfer coeff (W/K)
    double Tc = 300.0;        // Cooling jacket temp (K)

    // 2. Instantiate the Reactor Model
    CSTRModel myReactor(F, V, C_A0, T0, k0, E_a, dH, rho, Cp, UA, Tc);

    // 3. Define Simulation Conditions
    double initial_Ca = 0.0;  // Start with empty reactor (Concentration = 0)
    double initial_T = 300.0; // Start at room temperature (300 K)
    double dt = 0.1;          // Time step (seconds)
    double t_final = 200.0;   // Run simulation for 200 seconds

    // 4. Instantiate Solver and Run
    RK4Solver solver;
    cout << "Starting CSTR Dynamic Simulation..." << endl;
    solver.solve(myReactor, initial_Ca, initial_T, dt, t_final, "reactor_simulation_results.csv");

    return 0;
}
