/**
 * @brief       Model functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.1.0
 * @date        2025-11-30
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef ZCURVE_MODEL
#define ZCURVE_MODEL

// The number of mlp models in meta.bin
#define N_MODELS  61
// The number of neural neurons
#define N_HIDDEN  200
// The number of total params in a mlp model
#define N_PARAMS  154931
// The upper limit for threshold of seed ORFs
#define UP_PROBA  0.9

#include <iostream>
#include <fstream>
#include <random>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include "Encoding.hpp"
#include "svm.h"

namespace model {
    /**
     * @brief               Initialize the mlp models from embedded source.
     * @return              True if the models are successfully initialized.
     */
    bool init_models();
    /**
     * @brief           Predict the output of a mlp model.
     * @param model_id  The index of the mlp model.
     * @param data      The input data.
     * @param size      The size of the input data.
     * @param probas    The output probabilities.
     */
    void mlp_predict(int index, double *data, int size, double *probas);
    /**
     * @brief               Train and predict the output of a SVM model.
     * @param params        The SVM parameters.
     * @param size          The size of the input data.
     * @param init_score    The initial scores.
     * @param score         The output scores.
     * @return              True if the training and prediction are successful.
     */
    bool train_predict(double *params, int size, double *init_score, double *score);
    bool GS_Finder(bio::orf_array &seeds, int flanking, std::mt19937_64 &ran_eng, bool QUIET);
}

#endif