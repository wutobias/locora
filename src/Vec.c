#include <math.h>
#include "Vec.h"

void crossp(double *vec1, double *vec2, double *result){

    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

}

double dotp(double *vec1, double *vec2){

    double prod;

    prod  = vec1[0] * vec2[0];
    prod += vec1[1] * vec2[1];
    prod += vec1[2] * vec2[2];

    return prod;

}

void norm3(double *vec1){

    double norm;

    norm = dotp(vec1, vec1);
    norm = sqrt(norm);

    vec1[0] /= norm;
    vec1[1] /= norm;
    vec1[2] /= norm;

}