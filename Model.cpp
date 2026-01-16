#include "Model.hpp"
#include<chrono>
extern const unsigned char _binary_meta_bin_start[];
extern const unsigned char _binary_meta_bin_end[];
#define FLOAT_ZERO 1E-10
// Minimum set size for SVM training.
static const int MIN_SET_SIZE = 5;
/* SVM parameters. */
static svm_parameter param = {
    C_SVC,   /* svm_type     */
    RBF   ,  /* kernel_type  */
    3,       /* degree       */
    1/DIM_S,   /* gamma        */
    0.0,     /* coef0        */
    200,     /* cache_size   */
    1e-3,    /* eps          */
    10.0,    /* C            */
    0,       /* nr_weight    */
    nullptr, /* weight_label */
    nullptr, /* weight       */
    0.0,     /* nu           */
    0.0,     /* p            */
    true,    /* shrinking    */
    false    /* probability  */
};
/* Model parameters. */
const float *MODELS[N_MODELS];
/* Start codon frequency (ATG:GTG:XTG) */
double START_FREQ[] = { 0.78, 0.14, 0.08 };
/**
 * @brief           Calculate dot product of two vectors using AVX instructions.
 * @param x         The first vector.
 * @param y         The second vector.
 * @param dim       The dimension of the vectors.
 * @return          The dot product of the two vectors.
 */
inline double dot_df_avx(const double* x, const float* y, const int dim) {
#ifdef __AVX__
    // ---- AVX version ----
    __m256d sum_vec = _mm256_setzero_pd();
    int d = 0;
    for (; d + 4 <= dim; d += 4) {
        __m256d xv = _mm256_loadu_pd(x + d);
        __m128  yf = _mm_loadu_ps(y + d);
        __m256d yv = _mm256_cvtps_pd(yf);
        __m256d prod = _mm256_mul_pd(xv, yv);
        sum_vec = _mm256_add_pd(sum_vec, prod);
    }
    double buf[4];
    _mm256_storeu_pd(buf, sum_vec);
    double sum = buf[0] + buf[1] + buf[2] + buf[3];
    for (; d < dim; ++d)
        sum += x[d] * y[d];
    return sum;
#else
    // ---- Non-AVX fallback ----
    double sum = 0.0;
    for (int d = 0; d < dim; ++d)
        sum += x[d] * y[d];
    return sum;
#endif
}

inline double dot_dd_avx(const double *x, const double *y, int dim)
{
#ifdef __AVX__
    __m256d sum = _mm256_setzero_pd();
    int i = 0;

    #ifdef __FMA__
    for (; i + 3 < dim; i += 4) {
        __m256d vx = _mm256_loadu_pd(x + i);
        __m256d vy = _mm256_loadu_pd(y + i);
        sum = _mm256_fmadd_pd(vx, vy, sum);
    }
    #else

    for (; i + 3 < dim; i += 4) {
        __m256d vx = _mm256_loadu_pd(x + i);
        __m256d vy = _mm256_loadu_pd(y + i);
        sum = _mm256_add_pd(sum, _mm256_mul_pd(vx, vy));
    }
    #endif

    __m128d lo = _mm256_castpd256_pd128(sum);
    __m128d hi = _mm256_extractf128_pd(sum, 1);
    __m128d s2 = _mm_add_pd(lo, hi);
    double result = _mm_cvtsd_f64(_mm_add_pd(s2, _mm_unpackhi_pd(s2, s2)));
    for (; i < dim; ++i)
        result += x[i] * y[i];

    return result;

#else
    double result = 0.0;
    for (int i = 0; i < dim; ++i)
        result += x[i] * y[i];
    return result;

#endif
}

static void normalize(double *v, int n) {
    double norm = std::sqrt(dot_dd_avx(v, v, n));
    if (norm < FLOAT_ZERO) return;
    for (int i = 0; i < n; ++i) v[i] /= norm;
}

bool model::init_models() {
    const size_t expected_size = 
        N_MODELS * N_PARAMS * sizeof(float);

    const size_t actual_size =
        _binary_meta_bin_end - _binary_meta_bin_start;

    if (actual_size != expected_size || 
        reinterpret_cast<uintptr_t>(_binary_meta_bin_start) %
        alignof(float) != 0) 
    {
        std::cerr << "Error: failed to load embedded model data\n";
        return false;
    }

    const float* base = reinterpret_cast<const float*>(_binary_meta_bin_start);

    for (int i = 0; i < N_MODELS; ++i) MODELS[i] = base + i * N_PARAMS;

    return true;
}

void model::mlp_predict(int index, double *data, int size, double *probas) {
    static double tmp = 3.0;
    if (!size) return;
    /* scale data */
    const float *means = MODELS[index], *stds = means + DIM_A;
    double *scaled = encoding::std_trans(data, size, DIM_A, means, stds);
    /* mlp process */
    const float *model = stds + DIM_A; // first hidden w
    const float *hid_bs = model + N_HIDDEN*DIM_A;  // first hidden b
    const float *out_ws = hid_bs + N_HIDDEN;  // output w
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < size; i ++) {
        double *x = scaled + i*DIM_A;
        // Hidden Layer (765, ) -> (200, )
        double hid_out[N_HIDDEN] = {0.0};
        for (int j = 0; j < N_HIDDEN; j ++) {
            const float *hid_w = model + j*DIM_A;
            // hid_out = hid_w * x + hid_b
            hid_out[j] = dot_df_avx(x, hid_w, DIM_A);
            // ReLU activation function
            hid_out[j] = std::max(0.0, hid_out[j]+(double)hid_bs[j]);
        }
        // Output Layer (200, ) -> (1, )
        double out_b = (double)model[N_PARAMS-1];
        // out = out_w * hid_out + out_b
        probas[i] = (dot_df_avx(hid_out, out_ws, N_HIDDEN) + out_b) / tmp;
        // Sigmoid activation function
        probas[i] = 1.0/(1.0 + std::exp(-(probas[i])));
    }
    delete[] scaled;
}

bool model::rbf_train(double *params, int size, double *p_scores, double *d_scores) {
    svm_problem prob;
    double mins[DIM_S], maxs[DIM_S];
    for (int j = 0; j < DIM_S; ++j) {
        mins[j] =  std::numeric_limits<double>::infinity();
        maxs[j] = -std::numeric_limits<double>::infinity();
    }

    // split data into training set
    int n_pos = 0, n_neg = 0;
    std::vector<int> train_indices;
    prob.y = new double[size];
    prob.l = 0;
    for (int i = 0; i < size; ++i) {
        double s = p_scores[i];
        if (s > UP_PROBA || (s < 0.01 && s > FLOAT_ZERO)) {
            double* row = params + i*DIM_A;
            train_indices.push_back(i);

            for (int j = 0; j < DIM_S; ++j) {
                double v = row[j];
                if (v < mins[j]) mins[j] = v;
                if (v > maxs[j]) maxs[j] = v;
            }
            if (s > UP_PROBA) { prob.y[prob.l] = 1.0; n_pos ++; } 
            else { prob.y[prob.l] = -1.0; n_neg ++; }
            ++prob.l;
        }
    }

    if (n_pos < MIN_SET_SIZE || n_neg < MIN_SET_SIZE) return false; 

    // preprocess data (scale data to [0, 1])
    for (int j = 0; j < DIM_S; ++j) {
        double intv = maxs[j] - mins[j];
        if (intv == 0.0) {
            for (int i = 0; i < size; ++i)
                params[i*DIM_A + j] = 0.0;
        } else {
            double inv_intv = 1.0 / intv;
            double minv = mins[j];

            for (int i = 0; i < size; ++i) {
                double* row = params + i*DIM_A;
                row[j] = (row[j] - minv) * inv_intv;
            }
        }
    }

    // calculate gamma
    double sum = 0.0;
    double sqsum = 0.0;
    const int total = prob.l * DIM_S;
    for (int i = 0; i < prob.l; ++i) {
        double* row = params + train_indices[i]*DIM_A;
        for (int j = 0; j < DIM_S; ++j) {
            double v = row[j];
            sum += v;
            sqsum += v * v;
        }
    }
    double mean = sum / total;
    double mean_sq = mean * mean;
    double var = (sqsum / total) - mean_sq;
    if (var <= 0) var = FLOAT_ZERO;
    param.gamma = 0.5 / (DIM_S * var);
    // train svm model
    prob.x = new double *[prob.l];
    for (int i = 0; i < prob.l; i ++)
        prob.x[i] = params + train_indices[i]*DIM_A;
    svm_model *model = svm_train(&prob, &param);
    if (!model) return false;
    // predict scores
#ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
#endif
    for (int i = 0; i < size; i ++)
        d_scores[i] = svm_predict_score(model, params + i*DIM_A);
    svm_free_and_destroy_model(&model);

    delete[] prob.x;
    delete[] prob.y;
    return true;
}