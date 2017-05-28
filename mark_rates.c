/**
 * Compute rates and other utilities.	
 */

#include "types.h"
#include "constants.h"
#include "rvgs.h"

#define AVG_ION_EN 1.4522e+06 // (from TRIM) Average ion energy (in eV) expended on damage from 10.5-MeV Fe.
#define AVG_NEUTRON_EN 40.6 // (from SPECTER) Total damage energy in keV produced by a neutron in ITER.

/**
 * Obtain sign of an integer:
 */

static inline int signof (int64 a) { return (a < 0) ? -1 : 1; }

/**
 * Obtain zero factor:
 */

static inline double zero (int a) { return (abs(a) > 1) ? 1.0 : 0.0; }

/**
 * This function samples energies from the cumulative PKA-energy distribution.  
 **/

static inline double sample_pka_energy (double xi) { 
  
  int i, index=0;
  double slope;
  
  for (i = 1; i < cpdf.size; i++){
    if ( (xi > cpdf.cumul[i-1]) && (xi < cpdf.cumul[i]) ){
      index = i;
      break;
    }
  }
  
  // Interpolation:
  slope = (cpdf.energy[index] - cpdf.energy[index - 1])/
    (cpdf.cumul[index] - cpdf.cumul[index - 1]);
  
  return slope*(xi - cpdf.cumul[index - 1]) + cpdf.energy[index - 1];
  
}


/**
 * Compute object rates. 
 */

double *compute_rates_objects (const struct object_t *all, 
			       const struct object_t *mobile,
			       int n1,
			       int n2,
			       double *tot_rate){
  
  const struct object_t *P_1;
  const struct object_t *P_2;
  int N;
  int i=0, j=0, k;
  int ndef;
  double R_1, R_2, R12;
  double R_1e, concentration;
  double agg_rate=0;
  double damage[CHANNELS]={0.0};
  // For bcc lattice.
  const double avol=ALATT*ALATT*ALATT/2; 
  const double jumpd=sqrt(3.0)*ALATT/2.0; 
  
  // Create total rate vector of size (hh1)*(hh2).
  N = n1*(n2 + 1 + LEVELS); 
  /* number of mobiles + 1 (1st-rate: absorption) + 
     + LEVELS (1st-rate: emissions) */
  double *rate = malloc(sizeof(*rate)* (N + CHANNELS));
  
  // Compute sink strengths;
  double *ss = malloc(2*sizeof(*ss));
  if (ss == NULL){
    fprintf (stderr, "Could not allocate memory for ss[2].\n");
    exit(EXIT_FAILURE);
  }
  compute_sinks(ss);
  
  // Define damage rates:
#ifdef ELECTRON
#define NRT 1;  /* Number of number of Frenkel pairs per electron.*/  
  damage[0] = DPA_RATE*DENSITY*VOLUME/NRT;
#elif NEUTRON
#define NRT 361;  /* Number of number of Frenkel pairs per neutron (SPECTER).*/
  damage[0] = DPA_RATE*DENSITY*VOLUME/NRT; /* damage[0] is the ion or neutron insertion rate. */
  if (CHANNELS>1) damage[1] = RATIO_HE*1.0e-06*DPA_RATE*DENSITY*VOLUME;    // damage[1] is the He insertion rate.
  if (CHANNELS>2) damage[2] = RATIO_H*1.0e-06*DPA_RATE*DENSITY*VOLUME;    // damage[2] is the H insertion rate.
#elif ION
#define NRT 14523; /* Number of vacancies per incident 10.5-MeV ion (TRIM). */
  damage[0] = DPA_RATE*DENSITY*VOLUME/NRT;
  if (CHANNELS>1) damage[1] = RATIO_HE*1.0e-06*DPA_RATE*DENSITY*VOLUME;    // damage[1] is the He insertion rate.
  if (CHANNELS>2) damage[2] = RATIO_H*1.0e-06*DPA_RATE*DENSITY*VOLUME;    // damage[2] is the H insertion rate.
#elif defined PKA
  damage[0] = PKA_RATE*VOLUME; /* Insert PKAs directly according to some rate in units of [pka/s] */
  if (CHANNELS>1) damage[1] = HE_RATE*1.0e-06*DENSITY*VOLUME;    // damage[1] is the He insertion rate.
  if (CHANNELS>2) damage[2] = H_RATE*1.0e-06*DENSITY*VOLUME;    // damage[2] is the H insertion rate.
#endif

  if (CHANNELS>3) 
    for (k=3; k<CHANNELS; k++)
      damage[k] = 0.0;   
  
  for (P_1 = all; P_1 != NULL; P_1 = P_1->hh1.next){
    
    ndef = P_1->attributes[0];
    
    // 0-th order rates (as many as CHANNELS) are inserted at the end of rate[].
    
    /* ( The pow(base, exp) function returns base raised to the exp power. 
     *   There's a domain error if base is zero and exp is less than or equal to zero. 
     *   There's also a domain error if base is negative and exp is not an integer.) 
     */

    if (ndef <= 0){
      R_1 = zero(ndef)*pow(3.0*fabs((double) ndef)*avol/4.0/PI, 0.333333333333333333) + jumpd;
      if (ndef != 0) 
	R_1e = pow(3.0*(fabs((double) ndef) - 1)*avol/4.0/PI, 0.333333333333333333) + jumpd;
      else 
	R_1e = jumpd;	
    } 
    else if (ndef > 0){
      R_1 = zero(ndef)*sqrt((double)ndef*avol/jumpd/PI) + jumpd;
      R_1e = zero(ndef)*sqrt(((double)ndef - 1)*avol/jumpd/PI) + jumpd; 
      // This ensures that point defects have dissociation rate =0, so no explicit condition is needed.
    }
    
    // Determine the sink strength:
    double sink_strength = (ndef < 0) ? ss[0] : ss[1];
    // sink_strength = 0; // For testing
    
    // Compute rate for 1-st order process: absorption to sinks.
    rate[i] = P_1->number*P_1->diff*sink_strength;
    agg_rate += rate[i];
    i++;
    
    // Compute rate for 1-st order process: emission of monomer.
    for (j=0; j<LEVELS; j++){
/*       rate[i+j] =  */
/* 	( jumpd/(jumpd + R_1e) )*( (4.0*PI*R_1e*R_1e)/(GAMMA*ALATT*ALATT) )*P_1->binding[j]*P_1->number; */
      if (P_1->attributes[j] != 0){
	int attr[LEVELS]={0};
	attr[j] = signof(P_1->attributes[j]);
	rate[i+j] = (4.0*PI*R_1e/avol)*compute_diff_coeff(attr)*P_1->binding[j]*((double) P_1->number);  
      }
      else 
	rate[i+j] = 0.0;
      
      agg_rate += rate[i+j];
    }    
    
    i += LEVELS;
    
    for (P_2 = mobile; P_2 != NULL; P_2 = P_2->hh2.next){
      
      if (P_1->key != P_2->key){
	concentration = (double) P_1->number * (double) P_2->number / VOLUME;
      }
      else{
	assert(P_1->number == P_2->number);
	concentration = (double) P_1->number*( (double) P_1->number - 1.0 )/2.0/VOLUME; // Take into account 1/2 of rate.
      }
      
      ndef = P_2->attributes[0];	
      
      if (ndef > 0) // Loops.
	R_2 = zero(ndef)*sqrt((double)ndef*avol/jumpd/PI) + jumpd;
      else // 3D objects.
	R_2 = zero(ndef)*pow(3.0*fabs((double) ndef)*avol/4.0/PI, 0.333333333333333333) + jumpd;
      
      R12 = R_1 + R_2; // Reaction radius.

      // compute rate for 2-nd order process: reactions between two species.
      rate[i] = 4.0*PI*concentration*R12*dimension_term(R12, P_1, P_2);     
            
      agg_rate += rate[i];
      i++;
      
    }  
  }
  
  // Append 0-th order insertion:
  assert(i == N);
  for (k=0; k<CHANNELS; k++){
    rate[N+k] = damage[k];
    agg_rate += damage[k];
  }
  
  *tot_rate = agg_rate;
  
  if (ss != NULL)
    free(ss);
  else {
    printf("Error freeing 'ss'\n");
    exit(EXIT_FAILURE);
  }
  
  return rate;
  
} 

/**
 * Execute event and return object undergoing it. 
 */

void process_event (struct object_t *my_object, 
		    int64            n1,
		    int              n2) {
  
  int64 mono_key=0;
  int64 cluster_key=0;
  int64 new_key=0;
  int i, j=1+LEVELS;
  int64 k2;
  int *attr, *clstr_attr;
  int dim, att[LEVELS]={0};
  struct object_t *new_object;
  struct object_t *other_object;
  struct object_t *P_object;
  
  if (n2 == 0) {  
//////////  Absorption of object at sink. //////////
    
    my_object->number--;
    
  }
  else if ( (n2 > 0) && (n2 <= LEVELS) ) { 
////////// Emission of monomer. //////////
    
    // Monomer:
    mono_key = signof(my_object->attributes[n2 - 1])*
      ((int64) pow( 10.0, (double) EXP10*(LEVELS - n2) )); /* signof() provides the correct sign */
    attr = get_attributes(mono_key);
    cluster_key = (int64) signof(my_object->attributes[0])*( abs(my_object->key) - abs(mono_key) );
    clstr_attr = get_attributes(cluster_key);
    
    if (find_object_hash(mono_key) == NULL){ /* Create new object for monomer. */
      add_object_hash(mono_key, attr, 1, 3, compute_diff_coeff(attr),
		      compute_bind_term(attr)); 
    }
    else{ /* Update object */
      new_object = find_object_hash(mono_key);
      new_object->number++;
    }	
    
    // Cluster:
    if (find_object_hash(cluster_key) == NULL){ /* Create new object for cluster. */
      dim = (clstr_attr[0] > 4) ? 1 : 3; /* Only SI clusters N>4 diffuse one-dimensionally. */      
      add_object_hash(cluster_key, clstr_attr, 1, dim, compute_diff_coeff(clstr_attr),
		      compute_bind_term(clstr_attr)); 
    }
    else{ /* Update object */
      new_object = find_object_hash(cluster_key);
      new_object->number++;
    }	
    
    // Update parent cluster:
    my_object->number--;
    
    if (attr != NULL && clstr_attr != NULL){
      free(attr);
      free(clstr_attr);
    }
    else {
      printf("Error freeing pointers for monomer emission\n");
      exit(EXIT_FAILURE);
    }
    
  }
  else if ( n2 > LEVELS ) { 
////////// Reactions: //////////
    
    for (P_object = m_objects; P_object != NULL; P_object = P_object->hh2.next){
      if (j == n2){
	k2 = P_object->key;
	break;    
      }
      j++;
    }

    other_object = find_object_hash(k2);
    /* I now have both reactants. */

    for (i=0; i<LEVELS; i++){
      att[i] = my_object->attributes[i] + other_object->attributes[i];
      cluster_key += labs(att[i])*((int64) pow( 10.0, (double) EXP10*(LEVELS - 1 - i) ));
    }
    
    cluster_key *= signof(att[0]);
    
    if (find_object_hash(cluster_key) == NULL){ /* Create new object for cluster. */
      dim = (att[0] > 4) ? 1 : 3; /* Only SI clusters N>4 diffuse one-dimensionally. */
      add_object_hash(cluster_key, att, 1, dim, compute_diff_coeff(att),
		      compute_bind_term(att)); 
    }
    else{ /* Update object */
      new_object = find_object_hash(cluster_key);
      new_object->number++;
    }	
    
    // Update reactants:
    my_object->number--;
    other_object->number--;
    
  }
  
  return;
}

/**
 * Compute sink strength due to dislocations, precipitates, and grain boundaries. 
 */

void compute_sinks (double *s) {
  
  /* The total sink strength for all defects are stored in the array s. 
     [0] for vacancies; [1] for SIAs; */

  double Zdv=1.0, Zdi=1.1;
  double Zodsv=0.0, Zodsi=0.0;
  double Zgbv=1.0, Zgbi=1.0;
  double Sd, Sods, Sgbv, Sgbi;
    
  // Dislocation sink strength: 
  Sd = DISLOCATION;
  // ODS-particle sink strength:
  Sods = 4*PI*ODS_R*ODS_DENSITY; 
  // Grain boundary sink strength:
  Sgbv = 6*sqrt(Zdv*Sd + Zodsv*Sods)/GRAIN_SIZE;
  Sgbi = 6*sqrt(Zdi*Sd + Zodsi*Sods)/GRAIN_SIZE;

  s[0] = Zdv*Sd + Zgbv*Sgbv + Zodsv*Sods;
  s[1] = Zdi*Sd + Zgbi*Sgbi + Zodsi*Sods;
  
}

/**
 * Compute object diffusion coefficient based on its attributes. 
 * FeCr data:
 */

double compute_diff_coeff (const int attr[LEVELS]) {

  const double jumpd=sqrt(3.0)*ALATT/2.0;
  const double fi = 0.9, fv = 0.7; // Diffusion correlationm factors.
  const double gi = 0.5, gv = 0.125; // Geometric factor for diffusion.
  double prefactor=0, energy_m=0;
  double diff=0;
  int check_all = 0;
  int check_He = 0;
  int check_H = 0;
  int check_C = 0;
  int i;
  
  for (i = 1; i < LEVELS; i++) {
    check_all |= attr[i];  
    if (i>=2) check_He |= attr[i];
    if (i>=3) check_H |= attr[i];
  }

  // Pure defect clusters:
  if (!check_all){

    if (attr[0]>0) { // SIAs
      if (abs(attr[0])==1){ // 1I
        prefactor = gi*jumpd*jumpd*fi*NU0;
        energy_m = 0.013;
      } else if (abs(attr[0])>1){ // >1I
        prefactor = gi*jumpd*jumpd*fi*NU0*pow(fabs(attr[0]),-0.5);
        energy_m = 0.013;
      }
    } else if (attr[0]<0) { // Vacancies.
      if (abs(attr[0])==1){ // 1V
        prefactor = gv*jumpd*jumpd*fv*NU0;
        energy_m = 1.66;
      } else if (abs(attr[0])>1){ // >1V
        prefactor = gv*jumpd*jumpd*fv*NU0*pow(0.001, fabs(attr[0]) - 1.0);
        energy_m = 1.66;
      }
    }       
    
  }
  
  else if (!check_He){ 
    
    if (attr[0]!=0) { // SIA-He and V-He clusters immobile.
      prefactor = 0.0;
    } else
      if (attr[1]==1){ // He1.
        prefactor = gi*jumpd*jumpd*NU0;
        energy_m = 0.01;
      } else  if (attr[1]==2){ // He2.
        prefactor = gi*jumpd*jumpd*NU0*0.01;
        energy_m = 0.03;
      } else  if (attr[1]==3){ // He3.
        prefactor = gi*jumpd*jumpd*NU0*0.01;
        energy_m = 0.05;
      } else  if (attr[1]>3){ // He>3.
        prefactor = gi*jumpd*jumpd*NU0*0.01;
        energy_m = 0.06;
      } else prefactor = 0.0;

  }
  
  else if (!check_H){ 
    
    ;
    
  }
  /* All data from [CS Becquart et al., J Nucl Mater 403 (2010) 75] */
  
  diff = prefactor*exp(-energy_m/KB/TEMPERATURE);  
  
  return diff;
  
}

/**
 * Compute object binding terms based on its attributes. 
 */

double *compute_bind_term(const int attr[LEVELS]) {

  double *bind;
  double energy_d[LEVELS]={0.0};
  double energy_b=0.0, energy_m=0.0;
  double attfreq=1.0;
  double efi=9.96, emi=0.013; // Ab initio migration and formation energies of V and SIA in pure W.
  double efv=3.23, emv=1.66;
  double eb2v=-0.1, eb2i=2.12, eb2he=1.03;
  double efhe=4.0, emhe = 0.01;
  int check_all=0, check_He=0, check_H=0;
  int i;
  
  bind = malloc(sizeof(*bind)*LEVELS);
  for (i = 0; i < LEVELS; i++) {
    bind[i] = 0;
    if (i>=1) check_all |= attr[i];  
    if (i>=2) check_He |= attr[i];
    if (i>=3) check_H |= attr[i];
  }
  
  // Pure defect clusters:
  if (!check_all){
    
    if (attr[0]>0) { // SIAs
      energy_m = emi;
      if (abs(attr[0])==1){ // 1I
        attfreq = 0.0;
      } else if (abs(attr[0])==2){ // 2I
        energy_b = 2.12;
      } else if (abs(attr[0])==3){ // 3I
        energy_b = 3.02;
      } else if (abs(attr[0])==4){ // 4I
        energy_b = 3.60;
      } else if (abs(attr[0])==5){ // 5I
        energy_b = 3.98;
      } else if (abs(attr[0])==6){ // 6I
        energy_b = 4.27;
      } else if (abs(attr[0])==7){ // 7I
        energy_b = 5.39;
      } else if (abs(attr[0])>7) // > 7I
        energy_b = efi + (eb2i - efi)*
          ( pow(fabs((double)attr[0]),0.6666667) - pow((fabs((double)attr[0]) - 1.0),0.6666667) )/0.5847;      
    } else if (attr[0]<0) { // Vacancies.
      energy_m = emv;
      if (abs(attr[0])==1){ // 1V
        attfreq = 0.0;
      } else if (abs(attr[0])==2){ // 2V
        energy_b = eb2v;
      } else if (abs(attr[0])==3){ // 3V
        energy_b = 0.04;
      } else if (abs(attr[0])==4){ // 4V
        energy_b = 0.64;
      } else if (abs(attr[0])==5){ // 5V
        energy_b = 0.72;
      } else if (abs(attr[0])==6){ // 6V
        energy_b = 0.89;
      } else if (abs(attr[0])==7){ // 7V
        energy_b = 0.72;
      } else if (abs(attr[0])==8){ // 8V
        energy_b = 0.88;
      } else if (abs(attr[0])>8) // > 8V
        energy_b = efv + (eb2v - efv)*
          ( pow(fabs((double)attr[0]),0.6666667) - pow((fabs((double)attr[0]) - 1.0),0.6666667) )/0.5874;  
    }
    
    energy_d[0] = energy_b; // + energy_m
    bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);

  }
  // He-defect clusters:
  else if (!check_He){
    assert (attr[1] != 0);
    
    if (attr[0]<0) { // He-V clusters:
      double ratio = fabs( ((double) attr[1])/((double) attr[0]) );
      printf("%dV - %dHe\n", abs(attr[0]), abs(attr[1]));      
      assert(ratio > 0);
      if (abs(attr[0])==1 && abs(attr[1]==1)){ // 1V-1He
        energy_d[0] = 4.6;
        energy_d[1] = energy_d[0];
      } else {
        energy_d[0] = 2.4 + 3.5*log10(ratio) + 1.7*log10(ratio)*log10(ratio);
        // binding energy of V to cluster.
        energy_d[1] = 4.6 - 1.1*log10(ratio) - 0.3*log10(ratio)*log10(ratio);
        // binding energy of He to cluster.
      }
      bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
      bind[1] = attfreq*exp(-energy_d[1]/KB/TEMPERATURE);      
    }
    else if (attr[0]>0) // He-SIA clusters.
      attfreq = 0.0; // No dissociation between He and SIA clusters.
    else if (attr[0]==0){ // pure He clusters.
      if (attr[1]==1) { // He1.
        attfreq = 0;
      } else if (attr[1]==2){ // He2.
        energy_b = 1.03;
      } else if (attr[1]==3){ // He3.
        energy_b = 1.36;
      } else if (attr[1]==4){ // He4.
        energy_b = 1.52;
      } else { // He>4.
        energy_b = efhe + (eb2he - efhe)*
          ( pow(fabs((double)attr[0]),0.6666667) - pow((fabs((double)attr[0]) - 1.0),0.6666667) )/0.5874;;
      }
      energy_d[1] = energy_b + emhe;
      bind[1] = attfreq*exp(-energy_d[1]/KB/TEMPERATURE);      
    }
  }
  
  // H-defect clusters:
  else if (!check_H)
    attfreq = 0.0;
  
  return bind;
}

/**
 * Compute 0-th order (insertion rate). 
 */


void get_insertion (const int channel) {
  
  struct object_t *my_object;
  
  int *attr;
  int i,j,k;
  int check_all=0;
  int64 e_key=0;
  int e_defs=1;
  
  int **damage;
  double tot_energy;
  double ion_energy=0.0, neutron_energy=0.0;
  const int types=2; // Only vacancies and SIAs in ion cascades.
  int ndef=0;
  
  switch (channel) {
    
  case 1: // Insert damage:
    
#ifdef ELECTRON // Insert Frenkel pairs:
    
    e_key = (int64) pow( 10.0, (double) EXP10*(LEVELS - 1) ); /* Key for SIA. */
    if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
      attr = get_attributes(e_key);
      add_object_hash(e_key, attr, e_defs, 3, compute_diff_coeff(attr),
		      compute_bind_term(attr));  /* Create new object for SIAs. */
      // Free pointers:
      if (attr != NULL)
	free(attr);
      else {
	printf("Error freeing 'attr in electron insert'\n");
	exit(EXIT_FAILURE);
      }
    }
    else{ /* Update object */      
      my_object = find_object_hash(e_key);
      my_object->number += e_defs;
    }
    
    e_key *= -1; /* Key for vacancy. */
    if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
      attr = get_attributes(e_key);
      add_object_hash(e_key, attr, e_defs, 3, compute_diff_coeff(attr),
		      compute_bind_term(attr));  /* Create new object for vacancies. */
      // Free pointers:
      if (attr != NULL)
	free(attr);
      else {
	printf("Error freeing 'attr in electron insert'\n");
	exit(EXIT_FAILURE);
      }
    }
    else{ /* Update object */
      my_object = find_object_hash(e_key);
      my_object->number += e_defs;
    }	
    
#elif defined NEUTRON

    //    printf("Neutron insertion\n");
    tot_energy = (double) Poisson(AVG_NEUTRON_EN);
    //    printf("Neutron damage energy [keV] = %le\n", tot_energy);
    while (neutron_energy < tot_energy){ // Damage energy from inserted neutron (from SPECTER).
      double pka_energy = 0.0;
      while (pka_energy < 0.62){ // Limit to produce a stable Frenkel pair in keV (from Troev et al (2011)).
        double xi1 = drand48();
        neutron_energy += pka_energy;
        pka_energy = sample_pka_energy(xi1);
      }
      damage = generate_pka_damage( pka_energy, &ndef );
#ifdef DEBUG
      printf ("PKA energy=%f [eV]\n", pka_energy);
#endif      
      for (j = 0; j < types; j++){
	int sign = (j == 0) ? -1 : 1;
	for(i = 0; i < ndef; i++){
	  if (damage[j][i] != 0){
	    e_key = (int64) sign*(i + 1)*(pow( 10.0, (double) EXP10*(LEVELS - 1))); /* Key for cluster. */
#ifdef DEBUG
	    if (sign < 0) printf("%d V%d cluster(s)\n", damage[j][i], i+1);
	    if (sign > 0) printf("%d I%d cluster(s)\n", damage[j][i], i+1);
	    printf("key %lld\n", e_key);
#endif
	    if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
	      attr = get_attributes(e_key);
	      int dim = (attr[0] > 4) ? 1 : 3; 
	      add_object_hash(e_key, attr, damage[j][i], dim, compute_diff_coeff(attr),
			      compute_bind_term(attr));  /* Create new object for cluster. */
	      // Free pointers:
	      if (attr != NULL)
		free(attr);
	      else {
		printf("Error freeing 'attr' in ion damage insert\n");
		exit(EXIT_FAILURE);
	      }
	    }
	    else{ /* Update object */      
	      my_object = find_object_hash(e_key);
	      my_object->number += damage[j][i];
	    } // else
	  } // if damage
	} // for i
      } // for j
      
      // Free pointers:
      if (damage != NULL){
	for (j = 0; j < types; j++)
	  free(damage[j]);
	free(damage); 
      }
      else {
	printf("Error freeing 'damage'\n");
	exit(EXIT_FAILURE);
      }
      
      neutron_energy += pka_energy;
    }    

#elif defined ION
    
    printf("Ion insertion\n");
    tot_energy = (double) Poisson(AVG_ION_EN);
    printf("Ion energy [eV] = %le\n", tot_energy);
    while (ion_energy < tot_energy){ // Energy per incident 24.2-MeV ion expended on recoils.
                                     //  (from TRIM calculation)
      double pka_energy = 0.0;
      while (pka_energy < 333.3334){ // Limit to produce a stable Frenkel pair.
	double xi1 = drand48();
	ion_energy += pka_energy;
	pka_energy = sample_pka_energy(xi1);
      }
      damage = generate_ion_damage( pka_energy, &ndef );
#ifdef DEBUG
      printf ("PKA energy=%f [eV]\n", pka_energy);
#endif      
      for (j = 0; j < types; j++){
	int sign = (j == 0) ? -1 : 1;
	for(i = 0; i < ndef; i++){
	  if (damage[j][i] != 0){
	    e_key = (int64) sign*(i + 1)*(pow( 10.0, (double) EXP10*(LEVELS - 1))); /* Key for cluster. */
#ifdef DEBUG
	    if (sign < 0) printf("%d V%d cluster(s)\n", damage[j][i], i+1);
	    if (sign > 0) printf("%d I%d cluster(s)\n", damage[j][i], i+1);
	    printf("key %lld\n", e_key);
#endif
	    if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
	      attr = get_attributes(e_key);
	      int dim = (attr[0] > 4) ? 1 : 3; 
	      add_object_hash(e_key, attr, damage[j][i], dim, compute_diff_coeff(attr),
			      compute_bind_term(attr));  /* Create new object for cluster. */
	      // Free pointers:
	      if (attr != NULL)
		free(attr);
	      else {
		printf("Error freeing 'attr' in ion damage insert\n");
		exit(EXIT_FAILURE);
	      }
	    }
	    else{ /* Update object */      
	      my_object = find_object_hash(e_key);
	      my_object->number += damage[j][i];
	    } // else
	  } // if damage
	} // for i
      } // for j
      
      // Free pointers:
      if (damage != NULL){
	for (j = 0; j < types; j++)
	  free(damage[j]);
	free(damage); 
      }
      else {
	printf("Error freeing 'damage'\n");
	exit(EXIT_FAILURE);
      }
      
      ion_energy += pka_energy;
    }    
    
#elif defined PKA

    printf("PKA insertion\n");      
    double pka_energy = 0.0;
    
    double xi1 = drand48();
    pka_energy = sample_pka_energy(xi1);
#ifdef DEBUG
    printf ("PKA energy=%f [eV]\n", pka_energy);
#endif
    if (pka_energy > 620) // Limit to produce a stable Frenkel pair in eV (from Troev et al (2011)).
      damage = generate_pka_damage( pka_energy, &ndef );
    else
      break;
    
    for (j = 0; j < types; j++){
      int sign = (j == 0) ? -1 : 1;
      for(i = 0; i < ndef; i++){
	if (damage[j][i] != 0){
	  e_key = (int64) sign*(i + 1)*(pow( 10.0, (double) EXP10*(LEVELS - 1))); /* Key for cluster. */
#ifdef DEBUG
	  if (sign < 0) printf("%d V%d cluster(s)\n", damage[j][i], i+1);
	  if (sign > 0) printf("%d I%d cluster(s)\n", damage[j][i], i+1);
	  printf("key %lld\n", e_key);
#endif
	  if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
	    attr = get_attributes(e_key);
	    int dim = (attr[0] > 4) ? 1 : 3; 
	    add_object_hash(e_key, attr, damage[j][i], dim, compute_diff_coeff(attr),
			    compute_bind_term(attr));  /* Create new object for cluster. */
	    // Free pointers:
	    if (attr != NULL)
	      free(attr);
	    else {
	      printf("Error freeing 'attr' in ion damage insert\n");
	      exit(EXIT_FAILURE);
	    }
	  }
	  else{ /* Update object */      
	    my_object = find_object_hash(e_key);
	    my_object->number += damage[j][i];
	  } // else
	} // if damage
      } // for i
    } // for j
    
      // Free pointers:
    if (damage != NULL){
      for (j = 0; j < types; j++)
	free(damage[j]);
      free(damage); 
    }
    else {
      printf("Error freeing 'damage'\n");
      exit(EXIT_FAILURE);
    }
    
    tot_energy += pka_energy;
    
#endif
    
    break;  
    
  case 2: // Insert He:
        
    if (drand48() < 0.0){ // Fraction of 1.7-MeV He that ends up as interstitial after insertion (P Erhart).   
      e_key = -1*( (int64) pow( 10.0, (double) EXP10*(LEVELS - 1) ) + (int64) pow( 10.0, (double) EXP10*(LEVELS - 2) ) ); 
      /* Key for He substitutional. */
      printf("s-He insertion\n");
    }
    else {
      e_key = (int64) pow( 10.0, (double) EXP10*(LEVELS - 2) ); 
      /* Key for He interstitial. */
      printf("i-He insertion\n");
    }  
    
    if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
      attr = get_attributes(e_key);
      add_object_hash(e_key, attr, 1, 3, compute_diff_coeff(attr),
		      compute_bind_term(attr));  /* Create new object for He ions. */
      // Free pointers:
      if (attr != NULL)
	free(attr);
      else {
	printf("Error freeing 'attr in He insert'\n");
	exit(EXIT_FAILURE);
      }
    }
    else{ /* Update object */
      my_object = find_object_hash(e_key);
      my_object->number++;
    }	
    
    break;
    
  case 3: // Insert H:
        
    if (drand48() < 0.0){
      e_key = -1*( (int64) pow( 10.0, (double) EXP10*(LEVELS - 1) ) + (int64) pow( 10.0, (double) EXP10*(LEVELS - 3) ) ); 
      /* Key for H substitutional. */
      printf("s-H insertion\n");
    }
    else {
      e_key = (int64) pow( 10.0, (double) EXP10*(LEVELS - 3) ); 
      /* Key for H interstitial. */
      printf("i-H insertion\n");
    }
    
    if (find_object_hash(e_key) == NULL){ /* If object does not exist. */
      attr = get_attributes(e_key);
      add_object_hash(e_key, attr, 1, 3, compute_diff_coeff(attr),
		      compute_bind_term(attr));  /* Create new object for He ions. */
      // Free pointers:
      if (attr != NULL)
	free(attr);
      else {
	printf("Error freeing 'attr in H insert'\n");
	exit(EXIT_FAILURE);
      }
    }
    else{ /* Update object */
       my_object = find_object_hash(e_key);
      my_object->number++;
    }	
    
    break;
  }
  
  return;
  
  
}

double dimension_term (double                 rab,
		       const struct object_t *a, 
		       const struct object_t *b) {

  double term;
  /*   int dimsum = a->dimensionality + b->dimensionality; */
  int dimsum = 6; 
  double alpha_a = -log(PI*PI*pow(rab,3.0)/VOLUME/a->number);
  double alpha_b = -log(PI*PI*pow(rab,3.0)/VOLUME/b->number);

   switch (dimsum) {
   case 6: // 3D + 3D
     term = a->diff + b->diff;     
     break;
   case 4: // 3D + 1D
     if (a->dimensionality == 1 && b->dimensionality == 3) 
       term = a->diff*(b->number/VOLUME)*(2.0*PI*pow(rab,3.0)) + b->diff;
     else 
       term = b->diff*(b->number/VOLUME)*(2.0*PI*pow(rab,3.0)) + a->diff;     
     break;
   case 2: // 1D + 1D
     term = a->diff/alpha_b + b->diff/alpha_a;
     break;
   } 
  
  return term;

}

/**
 * This function generates a group of SIA- and V-type defects representing 
 * damage produced by a single cascade of a given energy
 * (automatic cascade generation from [L Malerba, JNM 351 (2006) 28]).
 */

int **generate_ion_damage (double energy,
			   int    *ndef) {

  const double fcli = 0.55; // Fraction of interstitials in clusters
  const double fclv = 0.25; // Fraction of vacancies in clusters, both from: [L Malerba, JNM 351 (2006) 28].
  const double fmd = 0.65; // kMC escape probabilty 
                           // from [Soneda and Diaz de la Rubia, Phil Mag A 78 (1998) 995].
  const int types = 2; // Only vacancies and SIAs in a cascade.
  
  double nrt, eta;
  double psia, pv;
  int n = 0;
  int k, j, i;
  int **damage;
  
  while (n < 1) { // Repeat until at least one stable Frenkel pair produced.
    eta = exp(-3.57*energy) + 0.3; // from: [L Malerba, JNM 351 (2006) 28].
    nrt = eta*0.8*energy/TDE/2;
    
    n = rint( Poisson(nrt) ); // Given the average number of defects, sample  from some statistical 
                              // distribution to give some variability (Poisson distribution only. 
  }    

  // Dynamically allocate first dimension of pointer array:
  damage = malloc(types * sizeof(int *)); // [0]: vacancy defects; [1]: SIA defects.
  if(damage == NULL){
    fprintf(stderr, "out of memory\n");
    exit(EXIT_FAILURE);
  }
  // Dynamically allocate second dimension of pointer array:
  for(i=0; i<types; i++){ 
    damage[i] = (int*) malloc(n * sizeof(int));
    if (damage[i] == NULL){
      fprintf(stderr, "out of memory\n");
      exit(EXIT_FAILURE);
    }
  }
  // This allocation overdimensions damage[][].

  // Initialize array:
  for(j=0; j<types; j++)
    for(i=0; i<n; i++)
      damage[j][i]=0;
  
  /* damage[2][n] is an integer array containing the counts (from 0 to the total number of defects)
     of defects and clusters ([0]: vacancy; [1]: SIA) produced by the cascade
     cascade. For example, damage[1][3]=2 means that two I3 clusters were formed. */
  
  *ndef = n;
  if (n > 1) {
    
    n *= fmd;   
    psia = 1.0 - pow( (1.0 - fcli), 1.0/( (double) types * (double) n - 1.0) ); 
    // Parameter for the binomial sampling of SIA-clusters.
    pv = 1.0 - pow( (1.0 - fclv), 1.0/( (double) types * (double) n - 1.0) ); 
    // Parameter for the binomial sampling of V-clusters.
    
    // Randomly sample SIA-clusters from the Binomial distribution:
    int nv = 0;
    while (nv < n){
      k = Binomial( (n - 1), pv) + 1;
      assert(k!=0);
      damage[0][k - 1]++;
      nv += k;
    }
    int nsia = 0;
    while (nsia < n){
      k = Binomial( (n - 1), psia) + 1;
      assert(k!=0);
      damage[1][k - 1]++;
      nsia += k;
    }
    
    // Complete with remaining point defects:
    if (nv < nsia) damage[0][0] += (nsia - nv);
    else if (nsia < nv) damage[1][0] += (nv - nsia);
  }
  else {       
    double xi = drand48();
    if (xi<fmd){
      damage[0][0] = 1;
      damage[1][0] = 1;
    }
  }

  return damage;
  
}


/**
 * This function generates a group of SIA- and V-type defects representing 
 * damage produced by a single cascade of a given energy
 * (automatic cascade generation from Troev et al., Fikar et al.).
 */

int **generate_pka_damage (double energy,
			   int    *ndef) {

  const double fcli = 0.5; // Fraction of interstitials in clusters
  const double fclv = 0.2; // Fraction of vacancies in clusters.
  // Both from Fikar et al (2009), Troev et al (2011).
  const double fmd = 1.0; // kMC escape probabilty.
  const int types = 2; // Only vacancies and SIAs in a cascade.

  double nrt, eta;
  double psia, pv;
  int n = 0;
  int i, j, k = 0;
  int **damage;
  
  while (n < 1) // Repeat until at least one stable Frenkel pair produced.
    n = rint(1.49*pow(energy,0.82)); // Relation from Fikar et al (2009).
  
  // Dynamically allocate first dimension of pointer array:
  damage = malloc(types * sizeof(int *)); // [0]: vacancy defects; [1]: SIA defects.
  if(damage == NULL){
    fprintf(stderr, "out of memory\n");
    exit(EXIT_FAILURE);
  }
  // Dynamically allocate second dimension of pointer array:
  for(i=0; i<types; i++){
    damage[i] = (int*) malloc(n * sizeof(int));
    if (damage[i] == NULL){
      fprintf(stderr, "out of memory\n");
      exit(EXIT_FAILURE);
    }
  }
  // This allocation hopefully overdimensions damage[][].
  // Initialize array:
  for(j=0; j<types; j++)
    for(i=0; i<n; i++)
      damage[j][i]=0;
  
  /* damage[2][n] is an integer array containing the counts (from 0 to the total number of defects)
     of defects and clusters ([0]: vacancy; [1]: SIA) produced by the cascade
     cascade. For example, damage[1][3]=2 means that two I3 clusters were formed. */
  
  *ndef = n;
  if (n > 1) {
    
    n *= fmd;
    psia = 1.0 - pow( (1 - fcli), 1.0/( (double) types * (double) n - 1.0) );
    // Parameter for the binomial sampling of SIA-clusters.
    pv = 1.0 - pow( (1 - fclv), 1.0/( (double) types * (double) n - 1.0) );
    // Parameter for the binomial sampling of V-clusters.
    
    // Randomly sample SIA-clusters from the Binomial distribution:
    int nv = 0;
    while (nv <= n){
      k = Binomial((n - 1), pv) + 1;
      assert(k!=0);
      damage[0][k - 1]++;
      nv += k;
    }
    int nsia = 0;
    while (nsia <= n){
      k = Binomial((n - 1), psia) + 1;
      assert(k!=0);
      damage[1][k - 1]++;
      nsia += k;
    }
    
    // Complete with remaining point defects:
    if (nv < nsia) damage[0][0] += (nsia - nv);
    if (nsia < nv) damage[1][0] += (nv - nsia);
  } else {
    double xi = drand48();
    if (xi<fmd){
      damage[0][0] = 1;
      damage[1][0] = 1;
    }
  }
  
  return damage;
  
}

/***
 * Get attributes from key. 
 */

int *get_attributes (int64 key) {
  
  int i;
  int64 okey = abs(key);
  int *attr;

  attr = malloc(sizeof(*attr)*LEVELS);
  
  for (i=0; i<LEVELS; i++){
    attr[i] = (double) okey / pow( 10.0, (double) EXP10*(LEVELS - 1 - i));
    okey -= ((int64) attr[i]) * ((int64) pow( 10.0, (double) EXP10*(LEVELS - 1 - i))); 
  }
  attr[0] *= signof(key);
  
  return attr;
}
