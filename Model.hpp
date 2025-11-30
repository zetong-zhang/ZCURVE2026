#ifndef ZCURVE_MODEL
#define ZCURVE_MODEL

// The number of mlp models in models.bin.gz
#define N_MODELS  55
// The number of neural neurons
#define N_HIDDEN  200
// The number of total params in a mlp model
#define N_PARAMS  154931

#include <iostream>
#include <fstream>
#include "Encoding.hpp"
#include "svm.h"

namespace model {
    bool init_models(const fs::path);
    void mlp_predict(int, double *, int, double *);
    bool train_predict(double *, int, double *, double *);
}

#endif