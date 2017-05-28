/**
 * Hash table operations.
 */

#include "types.h"

#define max(x, y) (((double) x > (double) y) ? x : y)

/* Add object to hash */

/*add_object_hash modified!!!*/
struct object_t* add_object_hash (int64         new_object_id,
                                  int           *new_object_attributes,
                                  int           new_object_number,
                                  unsigned int  new_object_dimensionality,
                                  double        new_object_diff,
                                  double        *new_object_binding) {
    
    struct object_t *P_object;
    int i_level;
    P_object = (struct object_t*)malloc(sizeof(*P_object));
    
    P_object->key = new_object_id;
    for (i_level=0; i_level<LEVELS; i_level++)
        P_object->attributes[i_level] = new_object_attributes[i_level];
    P_object->number = new_object_number;
    P_object->dimensionality = new_object_dimensionality;
    P_object->diff = new_object_diff;
    for (i_level=0; i_level<LEVELS; i_level++)
        P_object->binding[i_level] = new_object_binding[i_level];
    
    HASH_ADD(hh1, all_objects, key, sizeof(int64), P_object);
    // key: name of key field.
    
    if (new_object_binding != NULL)
        free(new_object_binding);
    else {
        printf("Error freeing 'new_object_binding'\n");
        exit(EXIT_FAILURE);
    }
    
    /*   if (P_object != NULL) */
    /*     free(P_object); */
    return P_object;
    
}
/* Delete object from hash */

void delete_object (struct object_t *P_object) {

  HASH_DELETE(hh1, all_objects, P_object);
  free(P_object);

}

/* Find an object in hash */

struct object_t *find_object_hash (int64 object_id) {
  
  struct object_t *P_object;
  HASH_FIND(hh1, all_objects, &object_id, sizeof(int64), P_object);
  return P_object;
}

struct object_t *find_object_hash2 (int64 object_id) {
  
  struct object_t *P_object;
  HASH_FIND(hh2, m_objects, &object_id, sizeof(int64), P_object);
  return P_object;
}

/* Sort objects by key */

int64 key_sort (const struct object_t *a, 
		const struct object_t *b) {
  return (b->key - a->key);
}

/**
 * Routines for sorting rate vector in descending order. 
 */ 

int mergesort (double *input, 
	       int     size) {
  
  double *scratch = (int *)malloc(size * sizeof(int));
  
  if(scratch != NULL){
    merge_helper(input, 0, size, scratch);
    free(scratch);
    return 1;
  }
  else
    return 0;
}

/* 'left' is the index of the leftmost element of the subarray; 'right' is one
   past the index of the rightmost element. */

void merge_helper(double *input, 
		  int     left, 
		  int     right, 
		  double *scratch) {
  
  /* base case: one element */
  if(right == left + 1)
    return;
  else{
    int i = 0;
    int length = right - left;
    int midpoint_distance = length/2;
    /* l and r are to the positions in the left and right subarrays */
    int l = left, r = left + midpoint_distance;
    
    /* sort each subarray */
    merge_helper(input, left, left + midpoint_distance, scratch);
    merge_helper(input, left + midpoint_distance, right, scratch);
    
    /* merge the arrays together using scratch for temporary storage */ 
    for(i = 0; i < length; i++){
      /* Check to see if any elements remain in the left array; if so,
       * we check if there are any elements left in the right array; if
       * so, we compare them.  Otherwise, we know that the merge must
       * use take the element from the left array */
      if(l < left + midpoint_distance && 
	 (r == right || max(input[l], input[r]) == input[l])){
	scratch[i] = input[l];
	l++;
      }
      else{
	scratch[i] = input[r];
	r++;
      }
    }
    /* Copy the sorted subarray back to the input */
    for(i = left; i < right; i++) 
      input[i] = scratch[i - left];
 
  }
}
