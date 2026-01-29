/**
 * @brief       Encoding functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.1.0
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
#define DIM_S 189

/* 
 * Map for converting ASCII chars into one-hot vectors
 *
 * A = [1, 0, 0, 0] G = [0, 1, 0, 0]
 * C = [0, 0, 1, 0] T = [0, 0, 0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */

extern double ONE_HOT[][4];

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
     * @brief       Encode sequences into Z-curve params.
     * @param orfs  The input sequence.
     * @param data  The output Z-curve params.
     */
    void     encode_seqs(pch_array &seqs, int len, double *data, int n_trans);
    /**
     * @brief       Encode ORFs into Z-curve params.
     * @param orfs  The input ORF array.
     * @param data  The output Z-curve params.
     */
    void     encode_orfs(bio::orf_array &orfs, double *data, int n_trans);
    /**
     * @brief           Transform Z-curve params to standardized data.
     * @param data      The input Z-curve params.
     * @param n         The number of samples.
     * @param dim       The dimension of Z-curve params.
     * @param means     The means of Z-curve params.
     * @param stds      The standard deviations of Z-curve params.
     * @return          The transformed Z-curve params.
     */
    double * std_scale(double *data, int n, int dim, const float *means, const float *stds);
    /**
     * @brief           Transform Z-curve params to min-max scaled data.
     * @param data      The input Z-curve params.
     * @param n         The number of samples.
     * @param dim       The dimension of Z-curve params.
     * @param mins      The minimums of Z-curve params.
     * @param maxs      The maximums of Z-curve params.
     * @return          The transformed Z-curve params.
     */
    double * minmax_scale(double *data, int n, int dim, double *mins, double *maxs);
    // TODO
    double   get_slope(double *params, int len);
    // TODO
    double   x_prime_curve(char *seq, int len, double *params);
    // TODO
    double   y_prime_curve(char *seq, int len, double *params);
    // TODO
    double   z_prime_curve(char *seq, int len, double *params);
    // TODO
    void     z_curve(char *seq, int len, double *params);
    /**
     * Find island regions by a sliding window calculating PCC
     * 
     * @param values y-values of 2D curves.
     * @param length Length of y-value array.
     * @param locs   Locations of found islands.
     * @param window Window size of the algorithm.
     * @param min_pcc Threshold of PCC.
     * 
     * @return count of islands.
     */
    int      find_island(double *values, int length, int window, double min_pcc, bio::region *root);
    /**
     * @brief           Apply Gaussian smoothing in-place.
     * @param params    The input/output parameter array.
     * @param len       The length of the parameter array.
     * @param sigma     The standard deviation of the Gaussian kernel.
     */
    void     gaussian_smooth(double *params, int len, int sigma);
}

#endif