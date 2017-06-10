/**
* Stochastic cluster dynamics implementation with Gillespie algorithm.
* J Marian and VV Bulatov, Apr 2010.
**/

#include "types.h"
#include "constants.h"
#include "RateMatrix.h"
#include "time.h"

#define is_mobile(x) (((struct object_t*)x)->diff > 0)

// Initialize structures:
struct object_t *all_objects = NULL, *m_objects = NULL;
struct cpdf_t cpdf = { 0,{ 0 },{ 0 } };

/*from "RateMatrix,h"*/
struct RateMatrix* rateMatrix;
struct ConstantValues constantValues;

int64 i_step = 0;

/*void displayAllObject() {
    FILE* fp = fopen("AllObject.txt", "wt+");
    fprintf(fp, "All Object Keys:\t");
    struct object_t* object;
    for (object = all_objects; object != NULL; object = (struct object_t*)(object->hh1.next)) {
        fprintf(fp, "%d\t", object->key);
    }
    fclose(fp);
}

void displayRateDifference(struct RateMatrix* rateMatrix) {
    FILE *fp = fopen("RateDifference.txt", "at+");

    int lIndex, cZFIndex, cSIndex;
    int lLength = rateMatrix->lineLength;
    int cZFLength = LEVELS + 1, cSLength = rateMatrix->secondLength;
    struct Material* lineHead = rateMatrix->lineHead;

    double realRate = 0.0;
    double* damage = rateMatrix->damage;
    for (lIndex = 0; lIndex < lLength; ++lIndex) {
        struct LinkedArrayNode* line = lineHead[lIndex].line;
        for (cZFIndex = 0; cZFIndex < cZFLength; ++cZFIndex) {
            realRate += line->array[cZFIndex];
        }
        line = line->next;
        for (cSIndex = 0; cSIndex < cSLength; ++cSIndex) {
            realRate += line->array[cSIndex];
        }
    }
    realRate += damage[0] + damage[1];
    fprintf(fp, "%lf\n", realRate - rateMatrix->totalRate);
    fclose(fp);
}*/

int
main(int argc, char *argv[])
{

	/**
	* (0) Declarations.
	*/

	// Seed:
	srand48((unsigned)time(NULL));
	//srand48( 223 );

	int i = 0, j = 0, k = 0;
	int process;
	int num_all = 0, num_m = 0;
	int n1 = 0, n2 = 0;
	int64 k1;
	int *attr;
	int nis = 0, nvs = 0, nhes = 0, nhs = 0;
	int lIndex, cIndex;


	double adv_time = 0.0;
	//double xi1, xi2,
	double dt;
	//double *rate;
	//double tot_rate = 0.0;
	//double integer;

	struct object_t *event_object = NULL, *P_object, *P_tmp;

	const char TIME[] = "time.out";
	const char SPECIES[] = "species.out";
    clock_t start, end; /*start and end times*/
    double elapsed; // elapsed CPU time in seconds
	FILE *fp, *fo, *ft, *fdl; // SKY: add fdl(diff length)

	/* Open output files:*/
	fp = fopen(TIME, "w");
	fo = fopen(SPECIES, "w");
    ft = fopen("timeconsumption.out", "wt+"); // SKY: .txt => .out
    fdl= fopen("dlengthMAX.out", "at+");      // SKY
	if (fp == NULL) {
		fprintf(stderr, "Error: Unable to open %s\n", TIME);
		exit(EXIT_FAILURE);
	}
	/*   fprintf(fp, "# SIAs SIAc SIA-He SIA-H SIA-He-H / Vs Vc V-He V-H V-He-H \n"); */
	fprintf(fp, "# time / SIAs SIAc SIA-He / Vs Vc V-He \n");
	if (fo == NULL) {
		fprintf(stderr, "Error: Unable to open %s\n", SPECIES);
		exit(EXIT_FAILURE);
	}
    fprintf(ft, "time consumption\n");

#if defined (ION) || defined (NEUTRON) || defined (PKA)

	const char CUM_PDF[] = "cpdf.dat";
	FILE *fc;

	/* Open file with PKA cumulative distribution function: */
	fc = fopen(CUM_PDF, "r");
	if (fc == NULL) {
		fprintf(stderr, "Error: Unable to open %s\n", CUM_PDF);
		exit(EXIT_FAILURE);
	}

	/* Read file: */
	while (!feof(fc)) {
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
	fscanf(fr, "i_step idim = %lld %d\n", &i_step, &dim);
	fscanf(fr, "num_all = %d\n", &num_all);
	fscanf(fr, "Aggregate time = %le\n", &adv_time);
	for (j = 0; j<num_all; j++) {
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

    /*initialization a rateMatrix, also compute damages and constant values*/
	rateMatrix = initialization();

	/***********************************************************
	* (1) Execute calculation loop.
	*/
    start = clock();

//////////////// SKY: choose one simulation time criterior
//	while (i_step <= TSTEPS) {
//	while (adv_time <= TOTAL_TIME) {
    while (adv_time <= TOTAL_DPA*1.0/DPA_RATE) {
//////////////// SKY: CHOOSE END

		int i, length = rateMatrix->lineLength;
		/**
		* (2) Compute and add rates.
		*/

		num_all = HASH_CNT(hh1, all_objects);
        // Separate mobile objects into separate hash:
		if (num_all != 0)
        {
			HASH_SELECT(hh2, m_objects, hh1, all_objects, is_mobile);
			num_m = HASH_CNT(hh2, m_objects);
		}

		addRemoveInRateMatrix(rateMatrix, m_objects);
		int idim = (num_m + 1 + LEVELS)*num_all + CHANNELS;

		if (i_step%PSTEPS == 0) {
			fprintf(fo, "i_step idim = %lld %d\n", i_step, idim);
			fprintf(fo, "num_all = %d\n", num_all);
#ifdef RATE_DUMP
			fprintf(frs, "i_step idim = %lld %d\n", i_step, idim);
			fprintf(frs, "num_all = %d\n", num_all);
			for (k = 0; k < idim; k++) {
				if (rate[k] != 0) fprintf(frs, "%le\n", rate[k]);
			}
#endif
		}

		 /**
		 * (4) Sample frequency line and choose event.
		 */

//        display(rateMatrix); // SKY

        computeDamages(rateMatrix, adv_time*DPA_RATE); // SKY: transmutation rate varies with dpa
        calibrateTotalRate(rateMatrix); // SKY
		dt = -log(drand48()) / rateMatrix->totalRate;
//        display(rateMatrix);
        selectReaction(rateMatrix, &lIndex, &cIndex);
		if (lIndex == -1)
		{			// damage
			struct TempMaterial* newMaterialList = get_insertion(cIndex);
			int count = 0;
			struct TempMaterial* newMaterial;
			for (newMaterial = newMaterialList; newMaterial != NULL; newMaterial = newMaterial->next)
			{
				++count;
			}
			rateMatrix->newMaterialList = newMaterialList;
			rateMatrix->newMaterialLength = count;
		}
		else
		{						// Matrix
			struct Material* material = rateMatrix->lineHead + lIndex;
			struct Reaction* reaction = rateMatrix->columnHead + cIndex;
			event_object = find_object_hash(material->object->key);
            if(event_object==NULL){
                int aaaa;
                ++aaaa;
            }
			int64 k2 = reaction->object == NULL ? 0 : reaction->object->key;

			struct TempMaterial* newMaterialList;                                                     // SKY
			if(i_step%PSTEPS==0) newMaterialList = process_event(event_object, k2, cIndex, adv_time); // SKY: add advtime
            else                 newMaterialList = process_event(event_object, k2, cIndex, -1.0);     // SKY: dont print
			
            int count = 0;
			struct TempMaterial* newMaterial;
			for (newMaterial = newMaterialList; newMaterial != NULL; newMaterial = newMaterial->next)
			{
				++count;
			}
			rateMatrix->newMaterialList = newMaterialList;
			rateMatrix->newMaterialLength = count;

			// Tally sink losses.
            if(cIndex==0)
            {
                attr = get_attributes(event_object->key);
				if (attr[0]>0) nis += abs(attr[0]);
				else if (attr[0]<0) nvs += abs(attr[0]);
				if (LEVELS>1 && attr[1] != 0) nhes += abs(attr[1]);
				if (LEVELS>2 && attr[2] != 0) nhs += abs(attr[2]);
				if (attr != NULL)
					free(attr);
				else {
					printf("Error freeing 'attr in restart'\n");
					exit(EXIT_FAILURE);
				}
			}
			if (i_step%PSTEPS == 0) {
				if (CHANNELS == 1) fprintf(fs, "%e %d %d\n", adv_time, nis, nvs);
				else if (CHANNELS == 2) fprintf(fs, "%e %d %d %d\n", adv_time, nis, nvs, nhes);
				else if (CHANNELS == 3) fprintf(fs, "%e %d %d %d %d\n", adv_time, nis, nvs, nhes, nhs);
            
                // SKY: --- check no diff out of box ---
                double lMAX= calculateDiffusionLength(rateMatrix); 
                fprintf(fdl, "%e\n", lMAX);
                // SKY: ---
			}
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

		if (i_step%PSTEPS == 0) {
			fprintf(fo, "Aggregate time = %e\n", adv_time);
#ifdef RATE_DUMP
			fprintf(frs, "Aggregate time = %e\n", adv_time);
#endif
		}

		int sia = 0, v = 0, siac = 0, vc = 0, siahe = 0, vhe = 0, siah = 0, vh = 0, siaheh = 0, vheh = 0;

		// Delete 'zero' objects from overall hash.

		P_object = all_objects;
		while (P_object != NULL) {

			int check_zero = 0;
			P_tmp = P_object->hh1.next; // Temporary pointer to avoid freeing a hash element in use.

			for (i = 0; i < LEVELS; i++)
				check_zero |= P_object->attributes[i];

			if ((!check_zero) || (P_object->number == 0) || (P_object->key == 0)) {
				/* 	delete_object(P_object); */
				HASH_DELETE(hh1, all_objects, P_object);
				deleteLine(rateMatrix, P_object->key);
                free(P_object);
			}
			else {
				if (i_step%PSTEPS == 0) {
					fprintf(fo, "object %lld, number %d\n", P_object->key, P_object->number);

					if (P_object->attributes[0]>0) {
						if (P_object->attributes[0] == 1)
							sia += P_object->number; // Single SIAs.
						else {
							siac += P_object->number; // Total number of SIA clusters.
							                               /*           if (P_object->attributes[1] != 0)
                                                           /*      siahe += P_object->number; // SIA clusters with He. */
                                                           /*           if (P_object->attributes[1]!=0 && P_object->attributes[2]==0) */
                                                           /*      siahe += P_object->number; // SIA clusters with He. */
                            if (P_object->attributes[1]==0 && P_object->attributes[2]!=0)
                                   siah += P_object->number; // SIA clusters with H.
                                                           /*    else if (P_object->attributes[1]!=0 && P_object->attributes[2]!=0) */
														   /* 		siaheh += P_object->number; // SIA clusters with He and H. */
						}
					}
					else if (P_object->attributes[0]<0) {
						if (P_object->attributes[0] == -1)
							v += P_object->number; // Single Vs.
						else {
							vc += P_object->number; // Total number of V clusters.
                                                         /*            if (P_object->attributes[1] != 0)
								                         /*        vhe += P_object->number; // SIA clusters with He. */
														 /* 	      if (P_object->attributes[1]!=0 && P_object->attributes[2]==0) */
														 /* 		vhe += P_object->number; // SIA clusters with He. */
                            if (P_object->attributes[1]==0 && P_object->attributes[2]!=0)
                                   vh += P_object->number; // V clusters with H.
														 /* 	      else if (P_object->attributes[1]!=0 && P_object->attributes[2]!=0) */
														 /* 		vheh += P_object->number; // V clusters with He and H. */
						}
					}
				}
			}
			P_object = P_tmp;
		}
		if (i_step%PSTEPS == 0) {
			fprintf(fp, "%le %le %le %le %le %le %le\n", adv_time,
				((double)sia) / VOLUME, ((double)siac) / VOLUME, ((double)siahe) / VOLUME, ((double)v) / VOLUME, ((double)vc) / VOLUME, ((double)vhe) / VOLUME);
			/*       fprintf(fp, "%le %le %le %le %le %le %le %le %le %le %le\n", adv_time, */
			/*     	      ((double)sia)/VOLUME, ((double)siac)/VOLUME, ((double)siahe)/VOLUME, ((double)siah)/VOLUME, ((double)siaheh)/VOLUME, */
			/*     	      ((double)v)/VOLUME, ((double)vc)/VOLUME, ((double)vhe)/VOLUME, ((double)vh)/VOLUME, ((double)vheh)/VOLUME); */
			printf("******** Iteration = %lld **** Dose = %e [dpa]. (Total, d0, d1 rates) %e %e %e\n", i_step, adv_time*DPA_RATE, rateMatrix->totalRate, rateMatrix->damage[0], rateMatrix->damage[1]);
            end = clock();
            elapsed = (end - start)/(double)CLOCKS_PER_SEC;
            fprintf(ft, "%lld %lf\n",i_step,elapsed);
		}

		/**
		* (8) Go to (2) unless desired time is reached.
		************************************************************/
		
        fflush(stdout); // SKY
        fflush(fo); // SKY
        fflush(fs); // SKY

        i_step++;

	}

	/**
	* Exit.
	*/

	fclose(fp);
	fclose(fo);
    fclose(ft);
    fclose(fdl); // SKY
#ifdef RATE_DUMP
	fclose(frs);
#endif

	exit(EXIT_SUCCESS);
	return(EXIT_SUCCESS);
}
