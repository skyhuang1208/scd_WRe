#pragma once

#include "constants.h"
#include <stdlib.h>
#include <time.h>
#include "types.h"

//#define NRT 361

extern int64 i_step;

//basic structure declaration
struct ConstantValues {
	double avol, jumpd;
	double ss[2];
};

struct LinkedArrayNode {
	double* array;
	struct LinkedArrayNode* next;
};

struct Material {
    double r1, r1e;
	struct object_t* object;
	struct LinkedArrayNode* line;
};

struct Reaction {
	struct object_t* object;
	char status;
};

struct TempMaterial {
	struct object_t* newObject;
	struct TempMaterial* next;
};

struct RateMatrix {
	struct Material* lineHead;
	int lineLength;
	int lineMaxLength;

	struct Reaction* columnHead;
	int secondLength;
	int secondMaxLength;

	struct TempMaterial* newMaterialList;
	int newMaterialLength;
	struct TempMaterial* newMobileList;
	int newMobileLength;

	double damage[CHANNELS];
	double totalRate;
};
//declaration end

//major operations that appeared in main.c
struct RateMatrix* initialization();

void addRemoveInRateMatrix(struct RateMatrix* rateMatrix, struct object_t* mobileObject);

void deleteLine(struct RateMatrix* rateMatrix, int64 key);

void update(struct RateMatrix* rateMatrix, int64 key);

void selectReaction(struct RateMatrix* rateMatrix, int* lIndex, int* cIndex);

//void display(struct RateMatrix* rateMatrix);

void calibrateTotalRate(struct RateMatrix* rateMatrix);
//end


//auxiliary functions that used in major oprations

/*1. functions for initialization()*/
inline double computeDamages(struct RateMatrix* rateMatrix, double dpa);

inline struct RateMatrix* mallocRateMatrix();

inline void initConstantValues();
/* 1.end */

/*2. functions for calculate 0/1/2 reaction*/
void computeR1R1e(struct RateMatrix* rateMatrix, struct Material* material);

inline double computeZeroReaction(const struct object_t* P_1, const double* ss);

inline double computeOneReaction(struct RateMatrix* rateMatrix, struct Material* material, int index);

inline double computeTwoReaction(struct RateMatrix* rateMatrix, struct Material* material, const struct object_t* P_2);

/*3.functions for addRemoveInRateMatrix()*/
void removeSecondArrayColumn(struct RateMatrix* rateMatrix, int deleteIndex);

inline void recordColumnHead(struct RateMatrix* rateMatrix, int requestSecondLength);

inline void calculateSecondColumnAndEnlargeIfNeed(struct RateMatrix* rateMatrix);

inline void enlargeLineHead(struct RateMatrix* rateMatrix, int requestLineLength);

struct Material* addLine(struct RateMatrix* rateMatrix, struct object_t* newObject);

inline void computeLine(struct RateMatrix* rateMatrix, struct Material* material, double* ss, char status);

/* functions for update()*/
void computeSecondArrayColumn(struct RateMatrix* rateMatrix, int secondIndex);

double calculateDiffusionLength(struct RateMatrix* rateMatrix);
