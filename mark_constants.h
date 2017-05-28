/**
 * Definition of constants in the code.
 * All definition are in CAP letters.
 */

// Universal constants:

#define PI 3.1415926
#define KB 8.617e-05       // [ev/K] Boltzmann's constant.
  
// Material properties:

#define DENSITY 6.30705e+22   // [atoms/cm^3] Atomic density for W.
#define ALATT 3.165e-08     // [cm] Lattice parameter for W.
#define DISLOCATION 1e+10   // [cm^-2] Dislocation density.
#define ODS_R 2.5e-07       // [cm] ODS-particle radius.
#define ODS_DENSITY 2.6e+17 // [cm^-3] ODS-particle density.
#define GRAIN_SIZE 0.01    // [cm] Typical grain size.
#define NU0 6.1e+12           // [Hz] Attempt frequency.
#define C_DENSITY 10        // [appm] C-atom density
#define GAMMA 1.0           // Fraction of surface emission.
#define TDE 90              // [eV] Threshold displacement energy for W.

// Run parameters:

#define PKA                 // (ELECTRON, NEUTRON, ION, PKA) Irradiation type.
#define TOTAL_TIME 1.0e+08  // [s] Total simulated time.
#define VOLUME 1.0e-12     // [cm^3] System volume.
#define TEMPERATURE 300.0  // [K] System temperature.
#define RATIO_HE 0.0      // [appm/dpa] He-to-dpa ratio.
#define RATIO_H 0.0
#define DPA_RATE 5.8e-08   // [dpa/s] Damage rate.
#define CHANNELS 2         // Irradiation channels used (1:W, 2:He, 3:H,...).
#define PSTEPS 50000        // Print data every so many.
#define TSTEPS 30000000    // Run these many steps.
#define EXP10 3.0           // Number of significant figures per key field.

// This is specific to Mark's code:
#define PKA_RATE 1.1e+15 // [pka/s/cm^3] for W
#define HE_RATE 6.7145e-8 // [appm He/s] for W
#define H_RATE 3.0957e-7 // [appm H/s] for W

// Auxiliary definitions: 

//#define RESTART            // Do restart.
/* #define RATE_DUMP          // Dump rate spectrum every PSTEPS. */
/* #define DEBUG              // Check cascade damage. */
