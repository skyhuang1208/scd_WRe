/**
 * Stochastic cluster dynamics implementation with Gillespie algorithm.
 * J Marian and VV Bulatov, Apr 2010.
 **/	

#include "types.h"
#include "constants.h"

#define is_mobile(x) (((struct object_t*)x)->diff > 0) 

// Initialize structures:
struct object_t *all_objects = NULL, *m_objects = NULL;
struct cpdf_t cpdf = {0, {0},{0}};

int 
main(int argc, char *argv[])
{
  
  /**
   * (0) Declarations. 
   */  
  
  // Seed:
  srand48( (unsigned) time(NULL) );
  //srand48( 223 );
    
  int i=0, j=0, k=0;
  int i_step = 0;
  int process;
  int num_all=0, num_m=0;
  int n1=0, n2=0;
  int64 k1;
  int *attr;
  int nis=0, nvs=0, nhes=0, nhs=0;

  double adv_time = 0.0;
  double xi1, xi2, dt;
  double *rate;
  double tot_rate = 0.0;
  double integer;

  struct object_t *event_object, *P_object, *P_tmp;

  const char TIME[] = "time.out";
  const char SPECIES[] = "species.out";
  FILE *fp, *fo;
  
  /* Open output files:*/
  fp = fopen(TIME, "w");
  fo = fopen(SPECIES, "w");
  if (fp == NULL) {
    fprintf(stderr,"Error: Unable to open %s\n", TIME);
    exit(EXIT_FAILURE);
  }
  /*   fprintf(fp, "# SIAs SIAc SIA-He SIA-H SIA-He-H / Vs Vc V-He V-H V-He-H \n"); */
  fprintf(fp, "# time / SIAs SIAc SIA-He / Vs Vc V-He \n");
  if (fo == NULL) {
    fprintf(stderr,"Error: Unable to open %s\n", SPECIES);
    exit(EXIT_FAILURE);
  }

#if defined (ION) || defined (NEUTRON) || defined (PKA)

  const char CUM_PDF[] = "cpdf.dat";
  FILE *fc;
  
  /* Open file with PKA cumulative distribution function: */
  fc = fopen(CUM_PDF, "r");
  if (fc == NULL) {
    fprintf(stderr,"Error: Unable to open %s\n", CUM_PDF);
    exit(EXIT_FAILURE);
  }
  
  /* Read file: */
  while (!feof(fc)){
    if (fscanf(fc, "%le %le", &cpdf.energy[i], &cpdf.cumul[i]) != 2)
      break;
    i++;
  }
  cpdf.size = i;
  printf("lines in cpdf file: %d\n", cpdf.size);
  fclose(fc);
#endif

#ifdef RESTART // June 2010.
  const char RESTART_FILE[] = "restart.out";
  FILE *fr;
  int dim;
  int64 object;
  int number;
  int dimensionality;
  
  /* Open restart file: */
  fr = fopen(RESTART_FILE, "r");
  if (fr == NULL) {
    fprintf(stderr, "Error: Unable to open %s\n", RESTART_FILE);
    exit(EXIT_FAILURE);
  }
  
  /* Read file: */
  fscanf(fr, "i_step idim = %d %d\n", &i_step, &dim);
  fscanf(fr, "num_all = %d\n", &num_all);
  fscanf(fr, "Aggregate time = %le\n", &adv_time);
  for (j=0; j<num_all; j++){
    fscanf(fr, "object %lld, number %d\n", &object, &number);
    printf("Read: object %lld, number %d\n", object, number);
    // Get attributes and dimensionality:
    attr = get_attributes(object);
    dimensionality = (attr[0] > 4) ? 1 : 3; /* Only SI clusters N>4 diffuse one-dimensionally. */
    add_object_hash(object, attr, number, dimensionality, 
		    compute_diff_coeff(attr), compute_bind_term(attr));
  }
  if (attr != NULL)
    free(attr);
  else {
    printf("Error freeing 'attr in restart'\n");
    exit(EXIT_FAILURE);
  };

  fclose(fr);
#endif

#ifdef DEBUG // June 2010.
  const char RATE_FILE[] = "system_rates.out";
  FILE *frs;
  
  /* Open rate file: */
  frs = fopen(RATE_FILE, "w");
  if (frs == NULL) {
    fprintf(stderr, "Error: Unable to open %s\n", RATE_FILE);
    exit(EXIT_FAILURE);
  }
#endif

  const char SINK_FILE[] = "sinks.out";
  FILE *fs;
  
  /* Open rate file: */
  fs = fopen(SINK_FILE, "w");
  if (fs == NULL) {
    fprintf(stderr, "Error: Unable to open %s\n", SINK_FILE);
    exit(EXIT_FAILURE);
  }
  /***********************************************************
   * (1) Execute calculation loop.
   */
  
  while (adv_time < TOTAL_TIME){

    /**
     * (3) Compute and add rates.
     */
    
    num_all = HASH_CNT(hh1, all_objects);
    
    // Separate mobile objects into separate hash:
    if (num_all != 0){
      HASH_SELECT(hh2, m_objects, hh1, all_objects, is_mobile); 
      num_m = HASH_CNT(hh2, m_objects);
    }
    
    // Compute rates:
    rate = 
      compute_rates_objects(all_objects, m_objects, num_all, num_m, &tot_rate);
    int idim = (num_m + 1 + LEVELS)*num_all + CHANNELS;
    if (i_step%PSTEPS == 0){
      fprintf(fo, "i_step idim = %d %d\n", i_step, idim);
      fprintf(fo, "num_all = %d\n", num_all);
#ifdef RATE_DUMP
      fprintf(frs, "i_step idim = %d %d\n", i_step, idim);
      fprintf(frs, "num_all = %d\n", num_all); 
      for (k = 0; k < idim; k++){
	if (rate[k]!=0) fprintf(frs, "%le\n", rate[k]);	
      }
#endif
    }

    /**
     * (4) Sample frequency line and choose event.
     */

    xi1 = drand48();
    dt = - log(xi1)/tot_rate;
    xi2 = drand48();
    
    double acc_rate = 0.0;
    j=0;
    
    // Sort 'rate' for improved efficiency.
    /*     mergesort(rate, idim); */

    for (process = 0; process < idim; process++){

      acc_rate += rate[process];
      if (acc_rate > xi2*tot_rate){ // rate[process] cannot be zero if this 'if' is satisfied.
	modf( ((double) process/((double)(num_m + 1 + LEVELS)) ), &integer);
	/* The function modf(double num, double *i) splits num into its integer 
	 * and fraction parts. It returns the fractional part and loads the integer part into i.
	 */ 
	if ((int) integer >=  num_all){	     
	  get_insertion(CHANNELS - idim + process);
	  break;
	}
	else {
	  n1 = (int) integer;
	  n2 = process % (num_m + 1 + LEVELS);
	  assert(n2 + n1*(num_m + 1 + LEVELS) == process);
	  for (P_object = all_objects; P_object != NULL; P_object = P_object->hh1.next){
	    if (j == n1){
	      k1 = P_object->key;
	      break;
	    }
	    j++;
	  }
	  
	  /* The object that will undergo a transformation (k1) has now been
	     chosen, as well,as which transformation (n2). */
	  
      
	  /**
	   * (5) Execute event.
	   */
	  
	  // Extract primary species involved:
	  event_object = find_object_hash(k1);
	  
	  // Identify transformation.
	  process_event(event_object, k1, n2);
	  
	  // Tally sink losses.
	  if (n2 == 0){
	    attr = get_attributes(event_object->key);
	    if (attr[0]>0) nis += abs(attr[0]);
	    else if (attr[0]<0) nvs += abs(attr[0]);
	    if (LEVELS>1 && attr[1]!=0) nhes += abs(attr[1]);
	    if (LEVELS>2 && attr[2]!=0) nhs += abs(attr[2]);
	    if (attr != NULL)
	      free(attr);
	    else {
	      printf("Error freeing 'attr in restart'\n");
	      exit(EXIT_FAILURE);
	    }
	  } 
	  if (i_step%PSTEPS == 0){
	    if (CHANNELS==1) fprintf(fs, "%e %d %d\n", adv_time, nis, nvs);   
	    else if (CHANNELS==2) fprintf(fs, "%e %d %d %d\n", adv_time, nis, nvs, nhes);
	    else if (CHANNELS==3) fprintf(fs, "%e %d %d %d %d\n", adv_time, nis, nvs, nhes, nhs);
	  }
	  
	  break;
	  
	} // if-else
      } // if
      
    } // for
 
    /* next: */
     if (rate != NULL)
       free(rate);
     else {
       printf("Error freeing 'rate'\n");
       exit(EXIT_FAILURE);
     }
     
    /**
     * (6) Update time and structures.
     */
    
    adv_time += dt;

    // Delete hash with mobile objects.
    HASH_CLEAR(hh2, m_objects);

    /**
     * (7) Analyze and print data.
     */
    
    // Sort hash by key.
    HASH_SRT(hh1, all_objects, key_sort);
    
    if (i_step%PSTEPS == 0){
      fprintf(fo, "Aggregate time = %e\n", adv_time); 
#ifdef RATE_DUMP
      fprintf(frs, "Aggregate time = %e\n", adv_time); 
#endif
    }
    
    int sia=0, v=0, siac=0, vc=0, siahe=0, vhe=0, siah=0, vh=0, siaheh=0, vheh=0;
    
    // Delete 'zero' objects from overall hash.

    P_object = all_objects;
    while (P_object != NULL){
      
      int check_zero = 0;
      P_tmp = P_object->hh1.next; // Temporary pointer to avoid freeing a hash element in use.
   
      for (i = 0; i < LEVELS; i++)
	check_zero |= P_object->attributes[i];
    
      if ( (!check_zero) || (P_object->number==0) || (P_object->key==0) ){
	/* 	delete_object(P_object); */
	HASH_DELETE(hh1, all_objects, P_object); 
	free(P_object); 
      }
      else {
	if (i_step%PSTEPS == 0){ 
	  fprintf(fo, "object %lld, number %d\n", P_object->key, P_object->number);
	  
	  if (P_object->attributes[0]>0){
	    if (P_object->attributes[0]==1)
	      sia += P_object->number; // Single SIAs.
	    else {
	      siac += P_object->number; // Total number of SIA clusters.	
	      if (P_object->attributes[1]!=0)
		siahe += P_object->number; // SIA clusters with He. 	
	      /* 	      if (P_object->attributes[1]!=0 && P_object->attributes[2]==0) */
	      /* 		siahe += P_object->number; // SIA clusters with He. */
	      /* 	      else if (P_object->attributes[1]==0 && P_object->attributes[2]!=0) */
	      /* 		siah += P_object->number; // SIA clusters with H. */
	      /* 	      else if (P_object->attributes[1]!=0 && P_object->attributes[2]!=0) */
	      /* 		siaheh += P_object->number; // SIA clusters with He and H. */
	    }
	  }
	  else if (P_object->attributes[0]<0){
	    if (P_object->attributes[0]==-1)
	      v += P_object->number; // Single Vs.
	    else {
	      vc += P_object->number; // Total number of V clusters.	      
	      if (P_object->attributes[1]!=0) 
		vhe += P_object->number; // SIA clusters with He. 
	      /* 	      if (P_object->attributes[1]!=0 && P_object->attributes[2]==0) */
	      /* 		vhe += P_object->number; // SIA clusters with He. */
	      /* 	      else if (P_object->attributes[1]==0 && P_object->attributes[2]!=0) */
	      /* 		vh += P_object->number; // V clusters with H. */
	      /* 	      else if (P_object->attributes[1]!=0 && P_object->attributes[2]!=0) */
	      /* 		vheh += P_object->number; // V clusters with He and H. */
	    }
	  }
	}
      }
      P_object = P_tmp;
    }

    if (i_step%PSTEPS == 0){
      fprintf(fp, "%le %le %le %le %le %le %le\n", adv_time, 
	      ((double)sia)/VOLUME, ((double)siac)/VOLUME, ((double)siahe)/VOLUME, ((double)v)/VOLUME, ((double)vc)/VOLUME, ((double)vhe)/VOLUME);
      /*       fprintf(fp, "%le %le %le %le %le %le %le %le %le %le %le\n", adv_time, */
      /*     	      ((double)sia)/VOLUME, ((double)siac)/VOLUME, ((double)siahe)/VOLUME, ((double)siah)/VOLUME, ((double)siaheh)/VOLUME, */
      /*     	      ((double)v)/VOLUME, ((double)vc)/VOLUME, ((double)vhe)/VOLUME, ((double)vh)/VOLUME, ((double)vheh)/VOLUME); */
#if defined ION || defined NEUTRON
      printf("******** Iteration = %d **** Dose = %e [dpa]\n", i_step, adv_time*DPA_RATE);
#elif defined PKA
      printf("******** Iteration = %d **** Dose = %e [pka]\n", i_step, adv_time*PKA_RATE*VOLUME);
#endif
    }
    
    /**
     * (8) Go to (2) unless desired time is reached.
     ************************************************************/
    i_step++;
    
  }
  
  /**
   * Exit.
   */
  
  fclose(fp);
  fclose(fo);
#ifdef RATE_DUMP
  fclose(frs);
#endif

  exit(EXIT_SUCCESS);
  return(EXIT_SUCCESS);
}
