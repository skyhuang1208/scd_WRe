/**
 * Definition of constants in the code.
 * All definition are in CAP letters.
 */

// Universal constants:

#define PI 3.1415926
#define KB 8.617e-05       // [ev/K] Boltzmann's constant.
  
// Material properties:
#define DENSITY 6.30705e+22 // [atoms/cm^3] Atomic density for W.
#define ALATT 3.16e-08      // [cm] Lattice parameter for W.
#define DISLOCATION 1e+10   // [cm^-2] Dislocation density.
#define ODS_R 2.5e-07       // [cm] ODS-particle radius.    !!!<OFF>!!!
#define ODS_DENSITY 2.6e+17 // [cm^-3] ODS-particle density.!!!<OFF>!!!
#define GRAIN_SIZE 0.01     // [cm] Typical grain size.     !!!<OFF>!!!
#define NU0I 1.50e+12       // [Hz] Attempt frequency of ITL
#define NU0V 6.46e+12       // [Hz] Attempt frequency of VCC
#define C_DENSITY 10        // [appm] C-atom density        !!!<OFF>!!!(probably)
#define GAMMA 1.0           // Fraction of surface emission.
#define TDE 90              // [eV] Threshold displacement energy for W.

// Run parameters:
#define VOLUME 1.0e-15      // [cm^3] System volume.
#define TEMPERATURE 1073.0  // [K] System temperature.
#define CHANNELS 2          // Irradiation channels used (1:W, 2:He, 3:H,...).
#define EXP10 4.0           // Number of significant figures per key field.

// Running time pars:
#define PSTEPS     1000000  // Print data every so many.
//#define TSTEPS  1000000000  // Run these many steps.
//#define TOTAL_TIME 4.8e+07  // [s] Total simulated time.
#define TOTAL_DPA      2.1  // [dpa] Total simulated dpas (DPA/rate)

// All rates 
#define PKA                 	// Irradiation type.
#define DPA_RATE 1.88742463e-07	// [dpa/s] Damage rate. <if PKA, only for cal other rates>
#define PKA_RATE        5.61e15	// [pka/s/cm^3] for W
#define RE_RATE1   181.65606738	// [appm/dpa] Re transmutation rate 
#define RE_RATE0   521.58147823	// rate= RE_RATE1 * dpa + RE_RATE0

#define RATIO_HE 0.0        // [appm/dpa] He-to-dpa ratio.  !!!<OFF>!!!
#define RATIO_H 0           //                              !!!<OFF>!!!

#define TRACK          // output Re contribution to "trackRE.out"

// Auxiliary definitions: 
/* #define RESTART            // Do restart. */
/* #define RATE_DUMP          // Dump rate spectrum every PSTEPS. */
/* #define DEBUG              // Check cascade damage. */
