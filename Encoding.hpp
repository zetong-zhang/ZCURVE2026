/**
 * @brief       Encoding functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.0.6-SNAPSHOT
 * @date        2025-11-30
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef ENCODING
#define ENCODING

#include <omp.h>
#include <cfloat>
#include <cmath>
#include <assert.h>
#include <iostream>
#include "BioStruct.hpp"

#define DIM_A 765

namespace encoding {
    /**
     * @brief           Calculate 1-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     mono_trans(const char *seq, int len, double *params);
    /**
     * @brief           Calculate 2-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     di_trans  (const char *seq, int len, double *params);
    /**
     * @brief           Calculate 3-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     tri_trans (const char *seq, int len, double *params);
    /**
     * @brief           Calculate 4-mer Z-curve params for a given sequence.
     * @param seq       The input sequence.
     * @param len       The length of the input sequence.
     * @param params    The output Z-curve parameters.
     */
    void     quart_trans(const char *seq, int len, double *params);
    /**
     * @brief       Encode a given sequence into Z-curve params.
     * @param seq   The input sequence.
     * @param len   The length of the input sequence.
     * @param data  The output Z-curve parameters.
     */
    void     encode     (const char *seq, int len, double *data, int n_trans);
    /**
     * @brief       Encode ORFs into Z-curve params.
     * @param orfs  The input ORF array.
     * @param data  The output Z-curve params.
     */
    void     encode_orfs(bio::orf_array &orfs, double *data);
    /**
     * @brief           Standardize Z-curve params.
     * @param data      The input Z-curve params.
     * @param n         The number of samples.
     * @param dim       The dimension of Z-curve params.
     * @param means     The means of Z-curve params.
     * @param stds      The standard deviations of Z-curve params.
     * @return          The standardized Z-curve params.
     */
    double * std_scale(double *data, int n, int dim, double *means, double *stds);
    /**
     * @brief           Transform Z-curve params to standardized data.
     * @param data      The input Z-curve params.
     * @param n         The number of samples.
     * @param dim       The dimension of Z-curve params.
     * @param means     The means of Z-curve params.
     * @param stds      The standard deviations of Z-curve params.
     * @return          The transformed Z-curve params.
     */
    double * std_trans(double *data, int n, int dim, float *means, float *stds);
    // TODO
    double   get_slope(double *params, int len);
    // TODO
    double   x_prime_curve(char *seq, int len, double *params);
    // TODO
    double   y_prime_curve(char *seq, int len, double *params);
    // TODO
    double   z_prime_curve(char *seq, int len, double *params);
}

#endif