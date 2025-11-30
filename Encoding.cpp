#include "Encoding.hpp"
static const int    X = 0, Y = 1, Z = 2, PHASE = 3;
static const double V = 1.0F/2, W = 1.0F/3, Q = 1.0F/4;
static int dims[] = { 9, 36, 144, 576 };
/* 
 * Map for converting ASCII chars into one-hot vectors
 *
 * A = [1, 0, 0, 0] G = [0, 1, 0, 0]
 * C = [0, 0, 1, 0] T = [0, 0, 0, 1]
 * 
 * Degenerate symbols are calculated by probabilities
 */
static double ONE_HOT[][4] = 
{   
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0},  
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0},
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0}, 
    {0, 0, V, V}, {1, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {1, 0, 0, 0}, {0, W, W, W}, {0, 0, 1, 0}, {W, W, 0, W},
    {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 1, 0, 0}, {W, 0, W, W}, 
    {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, V, 0, V}, {0, 0, 0, 0}, 
    {V, 0, V, 0}, {Q, Q, Q, Q}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
    {0, 0, 0, 0}, {V, V, 0, 0}, {0, V, V, 0}, {0, 0, 0, 1},
    {0, 0, 0, 1}, {W, W, W, 0}, {V, 0, 0, V}, {0, 0, 0, 0},
    {0, 0, V, V}, {1, 0, 0, 0}
};

/* 
 * Map for converting ASCII chars into Z-curve coordinates
 *
 * A = [+1, +1, +1]  G = [+1, -1, -1]
 * C = [-1, +1, -1]  T = [-1, -1, +1]
 * 
 * Degenerate symbols are calculated by weighted vector sum
 */
static double Z_COORD[][3] = {
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+1, +1, +1}, {-W, -W, -W}, {-1, +1, -1}, {+W, -W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+1, -1, -1}, {-W, +W, +W},
    {+0, +0, +0}, {+0, +0, +0}, {+0, -1, +0}, {+0, +0, +0},
    {+0, +1, +0}, {+0, +0, +0}, {+0, +0, +0}, {+0, +0, +0},
    {+0, +0, +0}, {+1, +0, +0}, {+0, +0, -1}, {-1, -1, +1},
    {-1, -1, +1}, {+W, +W, -W}, {+0, +0, +1}, {+0, +0, +0},
    {+1, +0, +0}, {+1, +1, +1}
};

void encoding::mono_trans(const char *seq, int len, double *params) {
    double counts[PHASE][3] = {{0.0f}};
        
    for (int i = 0; i < len; i ++) {
        int p = i % PHASE;

        counts[p][X] += Z_COORD[seq[i]][X];
        counts[p][Y] += Z_COORD[seq[i]][Y];
        counts[p][Z] += Z_COORD[seq[i]][Z];
    }
        
    for (int p = 0; p < PHASE; p ++, params += 3) {
        params[X] = counts[p][X] / len * PHASE;
        params[Y] = counts[p][Y] / len * PHASE;
        params[Z] = counts[p][Z] / len * PHASE;
    }
}

void encoding::di_trans(const char *seq, int len, double *params) {
    double counts[PHASE][4][3] = {{{0.0f}}}, denom;
    int sublen = len - 1;

    for (int i = 0; i < sublen; i++) {
        int p = i % PHASE;
        for (int b = 0; b < 4; b ++) {
            counts[p][b][X] += ONE_HOT[seq[i]][b] * Z_COORD[seq[i + 1]][X];
            counts[p][b][Y] += ONE_HOT[seq[i]][b] * Z_COORD[seq[i + 1]][Y];
            counts[p][b][Z] += ONE_HOT[seq[i]][b] * Z_COORD[seq[i + 1]][Z];
        }
    }
        
    for (int p = 0; p < PHASE; p ++)
        for (int b = 0; b < 4; b ++, params += 3) {
            params[X] = counts[p][b][X] / len * PHASE;
            params[Y] = counts[p][b][Y] / len * PHASE;
            params[Z] = counts[p][b][Z] / len * PHASE;
        }
}

void encoding::tri_trans(const char *seq, int len, double *params) {
    double counts[PHASE][4][4][3] = {{{{0.0f}}}}, denom;
    int sublen = len - 2;

    for (int i = 0; i < sublen; i ++) {
        int p = i % PHASE;
        for (int s = 0; s < 4; s ++)
        for (int b = 0; b < 4; b ++) {
            counts[p][s][b][X] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                  Z_COORD[seq[i + 2]][X];
            counts[p][s][b][Y] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                  Z_COORD[seq[i + 2]][Y];
            counts[p][s][b][Z] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                  Z_COORD[seq[i + 2]][Z];
        }
    }

    for (int p = 0; p < PHASE; p ++)
    for (int s = 0; s < 4; s ++)
    for (int b = 0; b < 4; b ++, params += 3) {
        params[X] = counts[p][s][b][X] / len * PHASE;
        params[Y] = counts[p][s][b][Y] / len * PHASE;
        params[Z] = counts[p][s][b][Z] / len * PHASE;
    }
}

void encoding::quart_trans(const char *seq, int len, double *params) {
    double counts[PHASE][4][4][4][3] = {{{{{0.0f}}}}};
    int sublen = len - 3;

    for (int i = 0; i < sublen; i ++) {
        int p = i % PHASE;
        for (int s = 0; s < 4; s ++)
        for (int b = 0; b < 4; b ++)
        for (int t = 0; t < 4; t ++) {
            counts[p][s][b][t][X] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                     ONE_HOT[seq[i + 2]][t] * Z_COORD[seq[i + 3]][X];
            counts[p][s][b][t][Y] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                     ONE_HOT[seq[i + 2]][t] * Z_COORD[seq[i + 3]][Y];
            counts[p][s][b][t][Z] += ONE_HOT[seq[i]][s] * ONE_HOT[seq[i + 1]][b] * 
                                     ONE_HOT[seq[i + 2]][t] * Z_COORD[seq[i + 3]][Z];
        }
    }

    for (int p = 0; p < PHASE; p ++)
    for (int s = 0; s < 4; s ++)
    for (int b = 0; b < 4; b ++)
    for (int t = 0; t < 4; t ++, params += 3) {
        params[X] = counts[p][s][b][t][X] / len * PHASE;
        params[Y] = counts[p][s][b][t][Y] / len * PHASE;
        params[Z] = counts[p][s][b][t][Z] / len * PHASE;
    }
}

void encoding::encode_orfs(bio::orf_array &orfs, double *data) {
    static void (*k_trans[4])(const char *, int, double *) = {
        mono_trans, di_trans, tri_trans, quart_trans
    };
    const int count = (int) orfs.size();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < count; i ++) {
        double *head = data + (i * DIM);
        for (int j = 0; j < 4; j ++) {
            (*k_trans[j])(orfs[i].pstr, orfs[i].len, head);
            head += dims[j];
        }
    }
}

void encoding::encode_orfs(bio::record_array &samples, double *data) {
    static void (*k_trans[4])(const char *, int, double *) = {
        mono_trans, di_trans, tri_trans, quart_trans
    };
    const int count = (int) samples.size();
#ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
#endif
    for (int i = 0; i < count; i ++) {
        double *head = data + (i * DIM);
        for (int j = 0; j < 4; j ++) {
            const char *pstr = samples[i].sequence.c_str();
            int seqlen = (int) samples[i].sequence.length();
            (*k_trans[j])(pstr, seqlen, head);
            head += dims[j];
        }
    }
}

double *encoding::std_scale(
    double *data, 
    int n, int dim,
    double *means,
    double *stds
) {
    double *cache = new double[n*dim];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            means[j] += data[i*dim+j];
        }
    }
    for (int j = 0; j < dim; j++) means[j] /= n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
            double diff = data[i*dim+j] - means[j];
            stds[j] += diff * diff;
        }
    }
    for (int j = 0; j < dim; j++) {
        stds[j] = std::sqrt(stds[j] / n);
        // 改这行代码时冯鑫师姐在唱: "忍不住化身一条固执的鱼~"
        if (stds[j] < 1E-5) stds[j] = 1.0;
    }
    for (int i = 0; i < n; i++) {
        int j;
        double *p_node = cache + i*dim;
        for (j = 0; j < dim; j++, p_node ++) {
            *p_node = (data[i*dim+j]-means[j])/stds[j];
        }
    }
    return cache;
}

double *encoding::std_trans(
    double *data, 
    int n, int dim,
    float *means,
    float *stds
) {
    double *cache = new double[n*dim];
    for (int i = 0; i < n; i++) {
        int j;
        double *p = cache + i*dim;
        for (j = 0; j < dim; j++, p ++) {
            *p = (data[i*dim+j]-(double)means[j])/(double)stds[j];
        }
    }
    return cache;
}