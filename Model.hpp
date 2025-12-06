/**
 * @brief       Model functions for Z-curve.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.0.5-SNAPSHOT
 * @date        2025-11-30
 * @license     GNU GPLv3
 * @contact     ylin@tju.edu.cn | fgao@tju.edu.cn
 */
#ifndef ZCURVE_MODEL
#define ZCURVE_MODEL

// The number of mlp models in meta.bin
#define N_MODELS  60
// The number of neural neurons
#define N_HIDDEN  200
// The number of total params in a mlp model
#define N_PARAMS  154931

#include <iostream>
#include <fstream>
#include "Encoding.hpp"
#include "svm.h"

namespace model {
    /**
     * @brief               Initialize the mlp models from a binary file.
     * @param model_path    The path to the binary file.
     * @return              True if the models are successfully initialized.
     */
    bool init_models(const fs::path model_path);
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
}

#endif