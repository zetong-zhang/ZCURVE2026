#ifndef ENCODING
#define ENCODING

#include <omp.h>
#include <cfloat>
#include <cmath>
#include <assert.h>
#include <iostream>
#include "BioStruct.hpp"

const int DIM  = 765;
const int SDIM = 189;

namespace encoding {
    void     mono_trans(const char *, int, double *);
    void     di_trans(const char *, int, double *);
    void     tri_trans(const char *, int, double *);
    void     quart_trans(const char *, int, double *);
    void     encode_orfs(bio::orf_array &, double *);
    void     encode_orfs(bio::record_array &, double *);
    double *std_scale(double *, int, int, double *, double *);
    double *std_trans(double *, int, int, float *, float *);
}

#endif