#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <malloc.h>


/**
 * Include uthash utilities.
 */

#include "uthash.h"

/**
 * Generic definitions.
 */

#define RETURN_SUCCESS EXIT_SUCCESS
#define RETURN_FAILURE EXIT_FAILURE
#define LEVELS 2 // 1: ion; 2: ion+He; 3: ion+He+H; 4: ion+He+H+?;
#define SIZE 40000 // Overdimensioned size of file "cpdf.dat"

typedef long long int int64;

/**
 * Object attributes.
 */

struct object_t{
  int64                key;                /* object id (hash key)                    */
  int                  attributes[LEVELS]; /* number of attributes                    */ 
  /* Level 0: number of defects (ndef, if <0 vacancies, if >0 SIAs)
     Level 1: number of He atoms.
     Level 2: number of H atoms.
     Level 3: number of C atoms.
     etc...
  */
  int                  number;             /* number of objects with these attributes */
  unsigned int         dimensionality;     /* 1 or 3, depending on defects            */
  double               diff;               /* diffusion coefficient                   */
  double               binding[LEVELS];    /* binding energy term                     */
  UT_hash_handle       hh1, hh2;           /* hash for all (1) and mobile (2) objects */
};

/**
 * Cumulative probability distribution function from TRIM.
 */

struct cpdf_t{
  int    size;
  double energy[SIZE];
  double cumul[SIZE];
};

/**
 * Declare structures:
 */

extern struct RateMatrix* rateMatrix;
extern struct ConstantValues constantValues;
extern struct object_t *all_objects, *m_objects;
extern struct cpdf_t cpdf;

/**
 * Prototypes:
 */

inline int signof (int64 a);
inline double zero (int a);

extern double *compute_rates_objects (const struct object_t *all,
				      const struct object_t *mobile,
				      int                    n1,
				      int                    n2,
				      double                *tot_rate);

/*modified!*/
extern struct TempMaterial* process_event(struct object_t *my_object,
                                          int64            k2,
                                          int              n2);

extern void compute_sinks (double *s);

extern double compute_diff_coeff (const int attr[LEVELS]);

extern double *compute_bind_term (const int attr[LEVELS]);

/*modified!*/
extern struct TempMaterial* get_insertion(const int channel);

/*modified!*/
extern struct object_t* add_object_hash(int64         new_object_id,
                                        int           new_object_attributes[LEVELS],
                                        int           new_object_number,
                                        unsigned int  new_object_dimensionality,
                                        double        new_object_diff,
                                        double        new_object_binding[LEVELS]);

extern void delete_object (struct object_t *P_object);

extern struct object_t *find_object_hash (int64 object_id);

extern struct object_t *find_object_hash2 (int64 object_id);

extern int long long key_sort (const struct object_t *a, 
			       const struct object_t *b);

extern double dimension_term (double                rab,
			      const struct object_t *a, 
			      const struct object_t *b);

extern int **generate_ion_damage (double energy,
				  int    *ndef);

extern int *get_attributes (int64 key);

extern int mergesort (double *input, 
		      int     size);

extern void merge_helper(double *input, 
			 int     left, 
			 int     right, 
			 double *scratch);

extern int **generate_pka_damage (double energy,
				  int    *ndef);
