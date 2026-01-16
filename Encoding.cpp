#include "Encoding.hpp"
/* const value for calculate PCC */
#define PCC_CONST 3.4641016151
// offset constants
static const int    X = 0, Y = 1, Z = 2, PHASE = 3;
// float constants
static const double V = 1.0F/2, W = 1.0F/3, Q = 1.0F/4;
// dimension constants
static int dims[] = { 9, 36, 144, 576 };

double ONE_HOT[][4] = 
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

void encoding::encode(const char *seq, int len, double *data, int n_trans) {
    static void (*k_trans[4])(const char *, int, double *) = {
        mono_trans, di_trans, tri_trans, quart_trans
    };
    for (int j = 0; j < n_trans; j ++) {
        (*k_trans[j])(seq, len, data);
        data += dims[j];
    }
}

void encoding::encode_orfs(bio::orf_array &orfs, double *data, int n_trans) {
    const int count = (int) orfs.size();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < count; i ++) {
        double *head = data + (i * DIM_A);
        encode(orfs[i].pstr, orfs[i].len, head, n_trans);
    }
}

double *encoding::std_trans(
    double *data, 
    int n, int dim,
    const float *means,
    const float *stds
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

double encoding::get_slope(double *params, int len) {
    double ySum, xySum, kp;

    for (int i = 0, ySum = 0, xySum = 0; i < len; i ++)
        ySum += params[i], xySum += params[i] * i;

    kp = (xySum / (len - 1) - ySum / 2) / len / (len + 1) * 12;

    return kp;
}

double encoding::x_prime_curve(char *seq, int len, double *params) {
    double *cache = new double[len];
    double counts = 0.0;
    double ySum, xySum, kp;
    int i;

    for (i = 0; i < len; i ++) {
        counts += Z_COORD[seq[i]][X];
        cache[i] = counts;
    }

    kp = get_slope(cache, len);
    
    for (i = 0; i < len; i++) {
        cache[i] -= kp * i;
        params[i] += cache[i];
    }

    delete[] cache;
    return kp;
}

double encoding::y_prime_curve(char *seq, int len, double *params) {
    double *cache = new double[len];
    double counts = 0.0;
    double ySum, xySum, kp;
    int i;

    for (i = 0; i < len; i ++) {
        counts -= Z_COORD[seq[i]][Y];
        cache[i] = counts;
    }

    kp = get_slope(cache, len);
    
    for (i = 0; i < len; i++) {
        cache[i] -= kp * i;
        params[i] += cache[i];
    }
    
    delete[] cache;
    return kp;
}

double encoding::z_prime_curve(char *seq, int len, double *params) {
    double *cache = new double[len];
    double counts = 0.0;
    double ySum, xySum, kp;
    int i;

    for (i = 0; i < len; i ++) {
        counts -= Z_COORD[seq[i]][Z];
        cache[i] = counts;
    }

    kp = get_slope(cache, len);
    
    for (i = 0; i < len; i++) {
        cache[i] -= kp * i;
        params[i] += cache[i];
    }
    
    delete[] cache;
    return kp;
}

void encoding::z_curve(char *seq, int len, double *params) {
    static double (*prime_curve[3])(char *, int, double *) = {
        y_prime_curve, z_prime_curve, x_prime_curve
    };
    double *p_data = params;
    for (int i = 0; i < 3; i ++) {
        (*prime_curve[i])(seq, len, p_data);
        p_data += len;
    }
}

int encoding::find_island(
    double *values, 
    int length,
    int window,
    double min_pcc,
    bio::region *root
) {
    // Window size-related constants
    const double lxx = sqrt(window);
    const int maxIndex = window - 1;
    // Values for PCC calculation
    double xySum = 0.0, ySum = 0.0, y2Sum = 0.0;

    // Initialization of sliding window
    for (int i = 0; i < window; i ++) {
        xySum += i * (double) values[i];
        ySum += (double) values[i];
        y2Sum += (double) values[i] * values[i]; 
    }

    // Initialization of PCC calculating
    double lxy = xySum / maxIndex - ySum / 2;
    double lyy = sqrt(y2Sum - ySum * ySum / window);
    double pcc = PCC_CONST * lxy / lyy / lxx;
    
    bool recording = pcc > min_pcc; // Recording switch
    // The start point of a new island and count of islands
    int start = 0, count = 0;

    bio::region *nextNode = root;
    const int stop = length - window;
    for (int winStart = 0, winEnd = window; winStart < stop; winStart ++, winEnd ++) {
        // Update the sum values of the next sliding window
        ySum += values[winEnd] - values[winStart];
        xySum += window * values[winEnd] - ySum;
        y2Sum += values[winEnd] * values[winEnd] - values[winStart] * values[winStart];

        // Calculate the PCC value of the next sliding window
        lxy = xySum / maxIndex - ySum / 2;
        lyy = sqrt(y2Sum - ySum * ySum / window);
        pcc = PCC_CONST * lxy / lyy / lxx;

        if (recording && pcc <= min_pcc) {
            nextNode->next = new bio::region(start, winEnd + 1);
            nextNode = nextNode->next;
            recording = false;
            count ++;
        } else if (!recording && pcc > min_pcc) {
            start = winStart + 1;
            recording = start > nextNode->end;
            if (!recording) nextNode->end = winEnd + 1;
        }
    }

    // Finally operation
    if (recording) {
        nextNode->next = new bio::region(start, length);
        count ++;
    }

    // Refine
    nextNode = root;
    int minPoint, maxPoint, endPoint;
    float minValue, maxValue;
    while(nextNode->next) {
        nextNode = nextNode->next;
        endPoint = nextNode->end;
        minPoint = nextNode->start, minValue = values[minPoint];
        maxPoint = nextNode->end, maxValue = values[maxPoint];
        
        for (int i = minPoint + 1; i < endPoint; i ++) {
            if (values[i] > maxValue)
                maxPoint = i, maxValue = values[i];
            else if (values[i] < minValue)
                minPoint = i, minValue = values[i];
        }
        
        nextNode->start = minPoint;
        nextNode->end = maxPoint;
    }

    return count;
}