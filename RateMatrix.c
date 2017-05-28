#include "RateMatrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>



// Initialization Functions
inline double computeDamages(struct RateMatrix* rateMatrix, double dpa) {
    // SKY: if input DPA rates, we can compute PKA using NRT
    // SKY: if input KPA rates, directly use it
    // Define damage rates:
    #ifdef ELECTRON
    #define NRT 1;  /* Number of number of Frenkel pairs per electron.*/  
    rateMatrix->damage[0] = DPA_RATE*DENSITY*VOLUME/NRT;
    #elif defined NEUTRON
    #define NRT 361;  /* Number of number of Frenkel pairs per neutron (SPECTER).*/
    rateMatrix->damage[0] = DPA_RATE*DENSITY*VOLUME/NRT; /* damage[0] is the ion or neutron insertion rate. */
    #elif defined ION 
    #define NRT 14523; /* Number of vacancies per incident 10.5-MeV ion (TRIM). */
    rateMatrix->damage[0] = DPA_RATE*DENSITY*VOLUME/NRT;
    #elif defined PKA
    rateMatrix->damage[0] = PKA_RATE*VOLUME; /* Insert PKAs directly according to some rate in units of [pka/s] */
    #endif
    if(dpa==0) printf("damage[0]=  %le\n", rateMatrix->damage[0]);
    //////////// SKY: END
	
    double agg_rate = rateMatrix->damage[0];
    double concentration_H = 2.68e+5, flux_H = 4.00e+16;
	// damage[0] is the ion insertion rate.
	if (CHANNELS > 1) {
        double rh = RATIO_HE;
        double dr = DPA_RATE;
        double de = DENSITY;
        double vl = VOLUME;
		rateMatrix->damage[1] = (RE_RATE1*dpa+RE_RATE0)*DPA_RATE*1.0e-06*DENSITY*VOLUME; // SKY: Re transmutation rate
        // SKY: NOTICE!! RE_RATE is a function of dpa [unit: appm/dpa]

		//rateMatrix->damage[1] = RATIO_HE*1.0e-06*DPA_RATE*DENSITY*VOLUME;    // damage[1] is the He insertion rate.
        if(dpa==0) printf("damage[1]=  %le\n", rateMatrix->damage[1]); // SKY
		agg_rate += rateMatrix->damage[1];
	}
	if (CHANNELS > 2) {
		rateMatrix->damage[2] = concentration_H * flux_H * VOLUME;    // damage[2] is the H insertion rate.
		agg_rate += rateMatrix->damage[2];
        printf("damage[2]=  %le\n", rateMatrix->damage[2]);
	}
	if (CHANNELS>3) {
		int k;
		for (k = 3; k<CHANNELS; k++)
			rateMatrix->damage[k] = 0.0;
	}
	return agg_rate;
}

inline struct RateMatrix* mallocRateMatrix() {
	struct RateMatrix* rateMatrix = (struct RateMatrix*)malloc(sizeof(struct RateMatrix));
	rateMatrix->lineMaxLength = 0;
	rateMatrix->lineLength = 0;
	rateMatrix->secondMaxLength = 0;
	rateMatrix->secondLength = 0;
	rateMatrix->lineHead = NULL;
	rateMatrix->columnHead = (struct Reaction*)malloc(sizeof(struct Reaction) * (LEVELS + 1));
	memset(rateMatrix->columnHead, 0, sizeof(struct Reaction) * (LEVELS + 1));

	rateMatrix->newMaterialList = NULL;
	rateMatrix->newMaterialLength = 0;
	rateMatrix->newMobileList = NULL;
	rateMatrix->newMobileLength = 0;

	rateMatrix->totalRate = computeDamages(rateMatrix, 0.0);

	return rateMatrix;
}

inline void initConstantValues() {
	constantValues.avol = ALATT*ALATT*ALATT / 2.0;
	constantValues.jumpd = sqrt(3.0)*ALATT / 2.0;
	compute_sinks(constantValues.ss);
}

struct RateMatrix* initialization() {
	struct RateMatrix* rateMatrix = mallocRateMatrix();
	initConstantValues();
	return rateMatrix;
}
// Initialization End

// Calculate Functions
void computeR1R1e(struct RateMatrix* rateMatrix, struct Material* material) {
    int ndef = material->object->attributes[0];
    double avol = constantValues.avol, jumpd = constantValues.jumpd;
    if (ndef <= 0) {
        material->r1 = zero(ndef)*pow(3.0*fabs((double)ndef)*avol / 4.0 / PI, 0.333333333333333333) + jumpd;
        if (ndef != 0)
            material->r1e = pow(3.0*(fabs((double)ndef) - 1)*avol / 4.0 / PI, 0.333333333333333333) + jumpd;
        else
            material->r1e = jumpd;
    }
    else if (ndef > 0) {
        /*
         double z = zero(ndef);
         double bs1 = (double)ndef*avol / jumpd / PI;
         double bs1e = ((double)ndef - 1)*avol / jumpd / PI;
         double s1 = sqrt(bs1);
         double s1e = sqrt(bs1e);
         */
        material->r1 = zero(ndef)*sqrt((double)ndef*avol / jumpd / PI) + jumpd;
        material->r1e = zero(ndef)*sqrt(((double)ndef - 1)*avol / jumpd / PI) + jumpd;
        // This ensures that point defects have dissociation rate =0, so no explicit condition is needed.
    }
}

inline double computeZeroReaction(const struct object_t* P_1, const double* ss) { // sink
    // Determine the sink strength:
    double sink_strength = (P_1->attributes[0] < 0) ? ss[0] : ss[1];
    // sink_strength = 0; // For testing

    // Compute rate for 1-st order process: absorption to sinks.
    return P_1->number*P_1->diff*sink_strength;
}

inline double computeOneReaction(struct RateMatrix* rateMatrix, struct Material* material, int index) { // dissociation
    struct object_t* P_1 = material->object;
    if (P_1->attributes[index] != 0) {
        int attr[LEVELS] = { 0 };
        attr[index] = signof(P_1->attributes[index]);
        
        if(P_1->attributes[0]>=4 && P_1->attributes[1]!=0) attr[index] = P_1->attributes[index]; // SKY: SIA-Re detach
        
        return (4.0 * PI * material->r1e / constantValues.avol)*compute_diff_coeff(attr)*P_1->binding[index] * ((double)P_1->number);
    }
    else
        return 0;
}

inline double computeTwoReaction(struct RateMatrix* rateMatrix, struct Material* material, const struct object_t* P_2) { // reaction
    struct  object_t* P_1 = material->object;
    double avol = constantValues.avol, jumpd = constantValues.jumpd;
    double concentration;
    if (P_1->key != P_2->key) {
        concentration = (double)P_1->number * (double)P_2->number / VOLUME;
    }
    else {
        assert(P_1->number == P_2->number);
        concentration = (double)P_1->number*((double)P_1->number - 1.0) / 2.0 / VOLUME; // Take into account 1/2 of rate.
    }

    int ndef = P_2->attributes[0];
    double R_2, R12;
    if (ndef > 0) // Loops.
        R_2 = zero(ndef)*sqrt((double)ndef*avol / jumpd / PI) + jumpd;
    else // 3D objects.
        R_2 = zero(ndef)*pow(3.0*fabs((double)ndef)*avol / 4.0 / PI, 0.333333333333333333) + jumpd;

    R12 = material->r1 + R_2; // struct Reaction radius.

    // compute rate for 2-nd order process: reactions between two species.
//    return 4.0*PI*concentration*R12*dimension_term(R12, P_1, P_2);
    return 4.0*PI*concentration*R12*dimension_term(R12, P_1, P_2); // SKY
}
// calculate End

// Add Remove Functions
void removeSecondArrayColumn(struct RateMatrix* rateMatrix, int deleteIndex) {
	int lIndex, cIndex, swapSecondIndex;
	int secondLength = rateMatrix->secondLength, deleteSecondIndex = deleteIndex - LEVELS - 1;
	int columnLength = LEVELS + 1 + secondLength;
	int lineLength = rateMatrix->lineLength;
    struct Reaction* columnHead = rateMatrix->columnHead;
	struct Material* lineHead = rateMatrix->lineHead;
	// Find the last column of available reaction column(not marked with 'd'), but the available column must be after the cuting one
	// If find, store its index
	for (cIndex = columnLength - 1; cIndex > deleteIndex; --cIndex) {
		if (columnHead[cIndex].status == 's') {
			break;
		}
	}

	if (cIndex > deleteIndex) { // If find an available column then
		// Swap the deleting column and the available column
		// But minus the rates of the removing column from rateMatrix->totalRate before swapping
		swapSecondIndex = cIndex - LEVELS - 1;
		for (lIndex = 0; lIndex < lineLength; ++lIndex) {
			double* secondArray = lineHead[lIndex].line->next->array;
			// Minus
			rateMatrix->totalRate -= secondArray[deleteSecondIndex];
			// Swapping
			secondArray[deleteSecondIndex] = secondArray[swapSecondIndex];
		}

		// Swap the column head.
		columnHead[deleteIndex] = columnHead[cIndex];

		// Set the stored index as the new length of second array
		rateMatrix->secondLength = swapSecondIndex;
	} else {
		// Or, only Minus the rates of the removing column from rateMatrix->totalRate
		for (lIndex = 0; lIndex < lineLength; ++lIndex) {
			rateMatrix->totalRate -= lineHead[lIndex].line->next->array[deleteSecondIndex];
		}

		// And then, check whether the deleting column is the last line according to the record.
		if (deleteSecondIndex < secondLength) {
			// If yes, minus columnLength by 1
			rateMatrix->secondLength = deleteSecondIndex;
		} // Or, do nothing.
	}
}

inline void recordColumnHead(struct RateMatrix* rateMatrix, int requestSecondLength) {
	int secondMaxLength = rateMatrix->secondMaxLength;
	int requestColumnLength = LEVELS + 1 + requestSecondLength;
	int columnLength = LEVELS + 1 + rateMatrix->secondLength;
	struct Reaction* columnHead = rateMatrix->columnHead;
	int cIndex;

	if (secondMaxLength < requestSecondLength) {
		struct Reaction* newColumnHead = (struct Reaction*)malloc(sizeof(struct Reaction) * requestColumnLength);
		memset(newColumnHead, 0, sizeof(struct Reaction) * requestColumnLength);

		for (cIndex = LEVELS + 1; cIndex < columnLength; ++cIndex) {
			newColumnHead[cIndex] = columnHead[cIndex];
		}
		free(columnHead);
		rateMatrix->columnHead = newColumnHead;
		columnHead = newColumnHead;
	}

	struct TempMaterial* newMobileObject;
	cIndex = columnLength;
	for (newMobileObject = rateMatrix->newMobileList; newMobileObject != NULL; newMobileObject = newMobileObject->next, ++cIndex) {
		columnHead[cIndex].object = newMobileObject->newObject;
	}
}

inline void calculateSecondColumnAndEnlargeIfNeed(struct RateMatrix* rateMatrix) {
	int requestSecondLength = rateMatrix->secondLength + rateMatrix->newMobileLength;
	recordColumnHead(rateMatrix, requestSecondLength);

	struct Material* lineHead = rateMatrix->lineHead;
	int lineLength = rateMatrix->lineLength, lineMaxLength = rateMatrix->lineMaxLength;
	int secondByteLength = sizeof(double) * requestSecondLength, secondLength = rateMatrix->secondLength;
	int lIndex, cSecondIndex;
	double addRate = 0.0, rate;
	for (lIndex = 0; lIndex < lineMaxLength; ++lIndex) {
		struct Material* material = lineHead + lIndex;
		struct LinkedArrayNode* secondArrayNode = material->line->next;
		double* secondArray = secondArrayNode->array;
		if (rateMatrix->secondMaxLength < requestSecondLength) {
			double* newSecondArray = (double*)malloc(secondByteLength);
			memset(newSecondArray, 0, secondByteLength);

			if (lIndex < lineLength) {
				for (cSecondIndex = 0; cSecondIndex < secondLength; ++cSecondIndex) {
					newSecondArray[cSecondIndex] = secondArray[cSecondIndex];
				}
			}

			free(secondArray);
			secondArrayNode->array = newSecondArray;
			secondArray = newSecondArray;
		}

		if (lIndex < lineLength) {
			struct Reaction* columnHead = rateMatrix->columnHead;
			for (cSecondIndex = secondLength; cSecondIndex < requestSecondLength; ++cSecondIndex) {
				struct object_t* mobileObject = columnHead[cSecondIndex + LEVELS + 1].object;
				rate = computeTwoReaction(rateMatrix, material, mobileObject);
				addRate += rate;
				secondArray[cSecondIndex] = rate;
			}
		}
	}
	rateMatrix->totalRate += addRate;
	rateMatrix->secondLength = requestSecondLength;
	rateMatrix->secondMaxLength = requestSecondLength;
}

inline void enlargeLineHead(struct RateMatrix* rateMatrix, int requestLineLength) {
	struct Material* lineHead = rateMatrix->lineHead;
	struct Material* newLineHead = (struct Material*)malloc(sizeof(struct Material) * requestLineLength);
	int lIndex = 0, lineMaxLength = rateMatrix->lineMaxLength;

	memset(newLineHead, 0, sizeof(struct Material) * requestLineLength);
	for (lIndex = 0; lIndex < lineMaxLength; ++lIndex) {
		newLineHead[lIndex] = lineHead[lIndex];
	}

	free(lineHead);
	rateMatrix->lineHead = newLineHead;
	rateMatrix->lineMaxLength = requestLineLength;
}

struct Material* addLine(struct RateMatrix* rateMatrix, struct object_t* newObject) {
	struct Material* newMaterial = rateMatrix->lineHead + rateMatrix->lineLength;
	int secondMaxLength = rateMatrix->secondMaxLength;
	newMaterial->object = newObject;

	// If the line attribute is NULL, meaning it is a new line, allocate the storage of this line.
	if (newMaterial->line == NULL) {
		struct LinkedArrayNode* zeroFirstArray = (struct LinkedArrayNode*)malloc(sizeof(struct LinkedArrayNode));
		struct LinkedArrayNode* secondArray = (struct LinkedArrayNode*)malloc(sizeof(struct LinkedArrayNode));

		zeroFirstArray->array = (double*)malloc(sizeof(double) * (LEVELS + 1));
		secondArray->array = (double*)malloc(sizeof(double) * secondMaxLength);

		zeroFirstArray->next = secondArray;
		secondArray->next = NULL;

		newMaterial->line = zeroFirstArray;
	}
	//Or(else), this line is abandoned by a delete action before, while the memory of the line still exist.

	++rateMatrix->lineLength;
	return newMaterial;
}

inline void computeLine(struct RateMatrix* rateMatrix, struct Material* material, double* ss, char status) {
	struct LinkedArrayNode* zeroFirstRates = material->line;
	struct LinkedArrayNode* secondRates = zeroFirstRates->next;
	double rate, addRate = 0.0;

	double* arrayzf = zeroFirstRates->array, *arrays = secondRates->array;

	// 0th Reaction

    computeR1R1e(rateMatrix, material);

    rate = computeZeroReaction(material->object, constantValues.ss);
	if (status == 'a') {
		addRate += rate;
	} else {
		addRate += rate - arrayzf[0];
	}
	arrayzf[0] = rate;

	// 1st Reaction
	int cFirstIndex;
	for (cFirstIndex = 0; cFirstIndex < LEVELS; ++cFirstIndex) {
		rate = computeOneReaction(rateMatrix, material, cFirstIndex);
		if (status == 'a') {
			addRate += rate;
		} else {
			addRate += rate - arrayzf[cFirstIndex + 1];
		}
		arrayzf[cFirstIndex + 1] = rate;
	}

	// 2nd Reaction
	int cSecondIndex, secondLength = rateMatrix->secondLength;
	struct Reaction* columnHead = rateMatrix->columnHead;
	for (cSecondIndex = 0; cSecondIndex < secondLength; ++cSecondIndex) {
		rate = computeTwoReaction(rateMatrix, material, columnHead[LEVELS + 1 + cSecondIndex].object);
		if (status == 'a') {
			addRate += rate;
		} else {
			addRate += rate - arrays[cSecondIndex];
		}
		arrays[cSecondIndex] = rate;
	}

	rateMatrix->totalRate += addRate;
}

void addRemoveInRateMatrix(struct RateMatrix* rateMatrix, struct object_t* mobileObjects) {
	// 0 Declare variable
	struct Reaction* columnHead = rateMatrix->columnHead;
	int cIndex, columnLength = LEVELS + 1 + rateMatrix->secondLength;

	// 1 Compare mobileObjects with columnTitle, detect which remain, to delete or to add.
	for (cIndex = LEVELS + 1; cIndex < columnLength; ++cIndex) {
		columnHead[cIndex].status = 'd';
	}

	struct object_t* mobileObject;
	for (mobileObject = mobileObjects; mobileObject != NULL; mobileObject = (struct object_t*)mobileObject->hh2.next) {
		for (cIndex = LEVELS + 1; cIndex < columnLength; ++cIndex) {
			if (mobileObject->key == columnHead[cIndex].object->key) {
				columnHead[cIndex].status = 's';
				break;
			}
		}
		if (cIndex == columnLength) {
			struct TempMaterial* newMobileObject = (struct TempMaterial*)malloc(sizeof(struct TempMaterial));
			newMobileObject->newObject = mobileObject;
			newMobileObject->next = rateMatrix->newMobileList;
			rateMatrix->newMobileList = newMobileObject;
			++rateMatrix->newMobileLength;
		}
	}

	// 2 Delete those status marked by 'd'.
	for (cIndex = LEVELS + 1; cIndex < columnLength; ++cIndex) {
		if (columnHead[cIndex].status == 'd') {
			removeSecondArrayColumn(rateMatrix, cIndex);
		}
	}

	// 3 Calculate rates of new columns for the exist lines, and increase the storage of columns and column heads if necessary.
	calculateSecondColumnAndEnlargeIfNeed(rateMatrix);

	// 4 Enlarge the storage of line heads.
	int requestLineLength = rateMatrix->lineLength + rateMatrix->newMaterialLength;
	if (rateMatrix->lineMaxLength < requestLineLength) {
		enlargeLineHead(rateMatrix, requestLineLength);
	}

	// 5 Allocate the storage of new lines if necessary and calculate rates of new lines.
	struct TempMaterial* newMaterial;
	for (newMaterial = rateMatrix->newMaterialList; newMaterial != NULL; newMaterial = newMaterial->next) {
		struct Material* material = addLine(rateMatrix, newMaterial->newObject);
		computeLine(rateMatrix, material, constantValues.ss, 'a');
	}

	struct TempMaterial* tempDelete, *tempNext;
	for (tempDelete = rateMatrix->newMaterialList; tempDelete != NULL; tempDelete = tempNext) {
		tempNext = tempDelete->next;
		free(tempDelete);
	}
	rateMatrix->newMaterialList = NULL;
	rateMatrix->newMaterialLength = 0;
	for (tempDelete = rateMatrix->newMobileList; tempDelete != NULL; tempDelete = tempNext) {
		tempNext = tempDelete->next;
		free(tempDelete);
	}
	rateMatrix->newMobileList = NULL;
	rateMatrix->newMobileLength = 0;
}
// Add Remove End

// Delete Line Function
void deleteLine(struct RateMatrix* rateMatrix, int64 key) {
	struct Material* lineHead = rateMatrix->lineHead;
	int lineLength = rateMatrix->lineLength;
	int lIndex;
	for (lIndex = 0; lIndex < lineLength; ++lIndex) {
		if (lineHead[lIndex].object->key == key) {
			break;
		}
	}
	if (lIndex == lineLength) {
		struct TempMaterial* deleteTemp, *prevTemp = NULL;
		deleteTemp = rateMatrix->newMaterialList;

        int count= 0;
        while (deleteTemp != NULL) {
            count ++;
            if (deleteTemp->newObject->key == key){
                break;
            }
            prevTemp = deleteTemp;
            deleteTemp = deleteTemp->next;
            
            if(count>10000){ // SKY: cant break the loop, exit with error
                printf("(deleteLine) key not expected. something's wrong");
                exit(1);
            }                // SKY
		}

		if (deleteTemp != NULL) {
			if (prevTemp == NULL) {
				rateMatrix->newMaterialList = deleteTemp->next;
			} else {
				prevTemp->next = deleteTemp->next;
			}
			free(deleteTemp);
			--rateMatrix->newMaterialLength;
		}
		return;
	}

	double* zeroFirstArray = lineHead[lIndex].line->array;
	double* secondArray = lineHead[lIndex].line->next->array;
	double minusRate = 0.0;
	int zeroFirstLength = LEVELS + 1, secondLength = rateMatrix->secondLength;
	int cZFIndex, cSIndex;
	for (cZFIndex = 0; cZFIndex < zeroFirstLength; ++cZFIndex) {
		minusRate += zeroFirstArray[cZFIndex];
	}
	for (cSIndex = 0; cSIndex < secondLength; ++cSIndex) {
		minusRate += secondArray[cSIndex];
	}
	rateMatrix->totalRate -= minusRate;

	struct Material temp = lineHead[lIndex];
	lineHead[lIndex].object = NULL;
	lineHead[lIndex] = lineHead[lineLength - 1];
	lineHead[lineLength - 1] = temp;
	--rateMatrix->lineLength;
}
// Delete Line End

// Update
void computeSecondArrayColumn(struct RateMatrix* rateMatrix, int secondIndex) {
	struct Material* lineHead = rateMatrix->lineHead;
	double* secondArray;
	int lIndex, lineLength = rateMatrix->lineLength;
	double changeRate = 0.0, rate;
	struct object_t* mobileObject = rateMatrix->columnHead[LEVELS + 1 + secondIndex].object;

	for (lIndex = 0; lIndex < lineLength; ++lIndex) {
		secondArray = lineHead[lIndex].line->next->array;
		rate = computeTwoReaction(rateMatrix, lineHead + lIndex, mobileObject);
		changeRate += rate - secondArray[secondIndex];
		secondArray[secondIndex] = rate;
	}
	rateMatrix->totalRate += changeRate;
}

void update(struct RateMatrix* rateMatrix, int64 key) {
	struct Material* lineHead = rateMatrix->lineHead;
	struct Reaction* columnHead = rateMatrix->columnHead;
	int lIndex, cIndex;
	int lineLength = rateMatrix->lineLength, columnLength = rateMatrix->secondLength + LEVELS + 1;

	for (lIndex = 0; lIndex < lineLength; ++lIndex) {
		if (lineHead[lIndex].object->key == key) {
			computeLine(rateMatrix, lineHead + lIndex, constantValues.ss, 'u');
			break;
		}
	}

	for (cIndex = LEVELS + 1; cIndex < columnLength; ++cIndex) {
		if (columnHead[cIndex].object->key == key) {
			computeSecondArrayColumn(rateMatrix, cIndex - LEVELS - 1);
			break;
		}
	}
}
// Update End

// Select Reaction
void selectReaction(struct RateMatrix* rateMatrix, int* lIndex, int* cIndex) {
    double xi2=drand48();
	double randRate = rateMatrix->totalRate * xi2;
    double addRate = 0.0, addRate2;
	struct Material* lineHead = rateMatrix->lineHead;
	int lineLength = rateMatrix->lineLength, zeroFirstLength = LEVELS + 1, secondLength = rateMatrix->secondLength;
	int tlIndex, tcZFIndex, tcSIndex;
	double* zeroFirstArray, *secondArray;
	double* damage = rateMatrix->damage;

	for (tlIndex = 0; tlIndex < lineLength; ++tlIndex) {
		struct LinkedArrayNode* line = lineHead[tlIndex].line;
		zeroFirstArray = line->array;
		for (tcZFIndex = 0; tcZFIndex < zeroFirstLength; ++tcZFIndex) {
			addRate += zeroFirstArray[tcZFIndex];
			if (addRate > randRate) {
				*lIndex = tlIndex;
				*cIndex = tcZFIndex;
				return;
			}
		}

		secondArray = line->next->array;
		for (tcSIndex = 0; tcSIndex < secondLength; ++tcSIndex) {
			addRate += secondArray[tcSIndex];
			if (addRate > randRate) {
				*lIndex = tlIndex;
				*cIndex = tcSIndex + zeroFirstLength;
				return;
			}
		}
	}
    addRate2 = addRate;
	for (tcZFIndex = 0; tcZFIndex < CHANNELS; ++tcZFIndex) {
		addRate += damage[tcZFIndex];
		if (addRate > randRate) {
			*lIndex = -1;
			*cIndex = tcZFIndex ;
			return;
		}
	}

    *lIndex = -1;
    *cIndex = CHANNELS - 1;
    printf("Warning: cant select a rate!! (Toatlrate) %e\n", rateMatrix->totalRate); // SKY
    calibrateTotalRate(rateMatrix);
}
// Select End

/*
void display(struct RateMatrix* rateMatrix) { // SKY MODIFIED
    FILE *fp = fopen("rmresult.txt", "at+");

    int lIndex, cIndex, cZFIndex, cSIndex;
    int lLength = rateMatrix->lineLength;
    int cZFLength = LEVELS + 1, cSLength = rateMatrix->secondLength, cLength = cZFLength + cSLength;
    struct Material* lineHead = rateMatrix->lineHead;
    struct Reaction* columnHead = rateMatrix->columnHead;

    fprintf(fp, "Column Head:\t");

    if (columnHead != NULL) {
        for (cIndex = LEVELS + 1; cIndex < cLength; ++cIndex) {
            fprintf(fp, "%d ", columnHead[cIndex].object->key);
        }
    }
    fprintf(fp, "\n");

    fprintf(fp, "Line Head:  \t");
    if (lineHead != NULL) {
        for (lIndex = 0; lIndex < lLength; ++lIndex) {
            fprintf(fp, "%d ", lineHead[lIndex].object->key);
        }
    }
    fprintf(fp, "\n");
    
    fprintf(fp, "In Object:  \t");
    struct object_t* object;
    for (object = all_objects; object != NULL; object = (struct object_t*)(object->hh1.next)) {
        fprintf(fp, "%d ", object->key);
    }
    fprintf(fp, "\n\n");

    fprintf(fp, "RateMatrix:\n");
    double rateInMatrix = 0.0;
    for (lIndex = 0; lIndex < lLength; ++lIndex) {
        struct LinkedArrayNode* line = lineHead[lIndex].line;
        for (cZFIndex = 0; cZFIndex < cZFLength; ++cZFIndex) {
            fprintf(fp, "%lf\t", line->array[cZFIndex]);
            rateInMatrix += line->array[cZFIndex];
        }
        line = line->next;
        for (cSIndex = 0; cSIndex < cSLength; ++cSIndex) {
            fprintf(fp, "%lf\t", line->array[cSIndex]);
            rateInMatrix += line->array[cSIndex];
        }
        fprintf(fp, "\n\n");
    }

    fprintf(fp, "MatrixRate: %lf\n", rateInMatrix);

    fprintf(fp, "TotalRate: %lf\n", rateMatrix->totalRate);

    fprintf(fp, "----------------------------\n\n", rateMatrix->totalRate);

    fclose(fp);
}
*/

// Calibrate TotalRate
void calibrateTotalRate(struct RateMatrix* rateMatrix) {
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
    int i;
    for(i=0; i<CHANNELS; i++){
        realRate += damage[i];
    }
    if(abs(realRate - rateMatrix->totalRate)>1) printf("Warning: rates inconsist (realRate, totalRate): %.5e %.5e\n", realRate, rateMatrix->totalRate); // SKY
    rateMatrix->totalRate = realRate;
}
// Calibrate End

// SKY: make sure diffusion length<VOL ---
double calculateDiffusionLength(struct RateMatrix* rateMatrix){
    double lMAX= 0; // SKY
    int lIndex, cIndex, cZFIndex, cSIndex;
    int lLength = rateMatrix->lineLength;
    int cZFLength = LEVELS + 1, cSLength = rateMatrix->secondLength, cLength = cZFLength + cSLength;
    struct Material* lineHead = rateMatrix->lineHead;
    struct Reaction* columnHead = rateMatrix->columnHead;
    for (lIndex = 0; lIndex < lLength; ++lIndex) {
        double RI = 0;
        struct LinkedArrayNode* line = lineHead[lIndex].line;
        int number  = lineHead[lIndex].object->number;
        double diff = lineHead[lIndex].object->diff;
        for (cZFIndex = 0; cZFIndex < cZFLength; ++cZFIndex) {
            RI += line->array[cZFIndex];
        }
        line = line->next;
        for (cSIndex = 0; cSIndex < cSLength; ++cSIndex) {
            RI += number * (line->array[cSIndex]);
        }
        double length = sqrt(diff/RI);
        if(length >lMAX) lMAX= length;
        double lcubic = length*length*length;
        if(lcubic>VOLUME) printf("Warning: out of box (l^3, VOL): %e %e\n", lcubic, VOLUME); // SKY
    }

    return lMAX;
}
// SKY
