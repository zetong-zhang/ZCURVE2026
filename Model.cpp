#include "Model.hpp"
extern const unsigned char _binary_meta_bin_start[];
extern const unsigned char _binary_meta_bin_end[];
#define FLOAT_ZERO 1E-25
// Minimum set size for SVM training.
static const int MIN_SET_SIZE = 5;
/* SVM parameters. */
static svm_parameter param = {
    C_SVC,   /* svm_type     */
    RBF,  /* kernel_type  */
    3,       /* degree       */
    1/DIM_A,   /* gamma        */
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
/* Start codon frequency (ATG:GTG:TTG:CTG) */
double START_FREQ[] = { 0.78, 0.14, 0.07, 0.01 };
/**
 * @brief           Calculate dot product of two vectors using AVX instructions.
 * @param x         The first vector.
 * @param y         The second vector.
 * @param dim       The dimension of the vectors.
 * @return          The dot product of the two vectors.
 */
inline double dot_product_avx(const double* x, const float* y, const int dim) {
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

static double* fisher_train(
    const double* data,
    int size,
    int dim,
    const int* labels
) {
    double* mu0 = new double[dim]();
    double* mu1 = new double[dim]();
    int n0 = 0, n1 = 0;

    for (int i = 0; i < size; ++i) {
        const double* x = data + i * dim;
        if (labels[i] == 0) {
            for (int d = 0; d < dim; ++d)
                mu0[d] += x[d];
            n0++;
        } else {
            for (int d = 0; d < dim; ++d)
                mu1[d] += x[d];
            n1++;
        }
    }

    if (n0 < MIN_SET_SIZE || n1 < MIN_SET_SIZE)
        return nullptr;

    for (int d = 0; d < dim; ++d) {
        mu0[d] /= n0;
        mu1[d] /= n1;
    }

    double* Sw = new double[dim * dim]();
    for (int i = 0; i < size; ++i) {
        const double* x = data + i * dim;
        const double* mu = (labels[i] == 0 ? mu0 : mu1);

        for (int r = 0; r < dim; ++r) {
            for (int c = 0; c < dim; ++c) {
                Sw[r * dim + c] += (x[r] - mu[r]) * (x[c] - mu[c]);
            }
        }
    }

    int D = dim;
    double* A = new double[D * 2 * D]();

    for (int i = 0; i < D; ++i) {
        for (int j = 0; j < D; ++j)
            A[i * 2 * D + j] = Sw[i * D + j];
        A[i * 2 * D + (D + i)] = 1.0;
    }

    for (int i = 0; i < D; ++i) {
        double pivot = A[i * 2 * D + i];
        if (std::fabs(pivot) < FLOAT_ZERO)
            return nullptr;

        for (int j = 0; j < 2 * D; ++j)
            A[i * 2 * D + j] /= pivot;

        for (int k = 0; k < D; ++k) {
            if (k == i) continue;
            double factor = A[k * 2 * D + i];
            for (int j = 0; j < 2 * D; ++j)
                A[k * 2 * D + j] -= factor * A[i * 2 * D + j];
        }
    }

    double* Sw_inv = new double[D * D];
    for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
            Sw_inv[i * D + j] = A[i * 2 * D + (D + j)];

    double* w = new double[dim]();
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j)
            w[i] += Sw_inv[i * dim + j] * (mu1[j] - mu0[j]);
    }

    delete[] mu0;
    delete[] mu1;
    delete[] Sw;
    delete[] A;
    delete[] Sw_inv;

    return w;
}

double fisher_score(
    const double* data,
    int dim,
    const double* w
) {
    double score = 0.0;
    for (int i = 0; i < dim; ++i)
        score += w[i] * data[i];
    return score;
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

static double *fisher_model(bio::orf_array &seeds, std::mt19937_64 &ran_eng) {
    const int n_seeds = (int) seeds.size();
    double *params_seeds = new double[DIM_S*n_seeds*2];
    int *labels = new int[n_seeds*2]();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < n_seeds; i ++) {
        const int seqlen = seeds[i].len;
        char *orgn_seq = seeds[i].pstr;
        char *shuf_seq = new char[seqlen];
        std::memcpy(shuf_seq, orgn_seq, seqlen*sizeof(char));
        std::shuffle(shuf_seq, shuf_seq+seqlen, ran_eng);
        encoding::encode(orgn_seq, seqlen, params_seeds+i*DIM_S, 3);
        encoding::encode(shuf_seq, seqlen, params_seeds+(n_seeds+i)*DIM_S, 3);
        labels[i] = 1;
        delete[] shuf_seq;
    }
    double *fisher_coef = fisher_train(params_seeds, n_seeds*2, DIM_S, labels);
    delete[] params_seeds;
    delete[] labels;
    return fisher_coef;
}

bool model::GS_Finder(bio::orf_array &seeds, int flanking, std::mt19937_64 &ran_eng, bool QUIET) {
    int n_alter = 0, i;
    for (i = 0; i < seeds.size(); i ++) n_alter += seeds[i].starts.size();
    std::cerr << "Total Alter Starts: " << std::setw(33) << n_alter << "\n";

    /* plot z-curves for the upstream flanking regions */
    int plot_range = flanking*2;
    double *curves = new double[plot_range*3]();
    for (i = 0; i < seeds.size(); i ++) {
        int up_start = seeds[i].t_start - flanking;
        int up_end = seeds[i].t_start + flanking;
        if (up_start < 0 || up_end > seeds[i].host_len) continue;
        char *pstr = seeds[i].pstr - flanking;
        encoding::z_curve(pstr, plot_range, curves);
    }
    return false;
}

void model::mlp_predict(int index, double *data, int size, double *probas) {
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
            hid_out[j] = dot_product_avx(x, hid_w, DIM_A);
            // ReLU activation function
            hid_out[j] = std::max(0.0, hid_out[j]+(double)hid_bs[j]);
        }
        // Output Layer (200, ) -> (1, )
        double out_b = (double)model[N_PARAMS-1];
        // out = out_w * hid_out + out_b
        probas[i] = dot_product_avx(hid_out, out_ws, N_HIDDEN) + out_b;
        // Sigmoid activation function
        probas[i] = 1.0/(1.0 + std::exp(-(probas[i])));
    }
    delete[] scaled;
}

bool model::train_predict(double *params, int size, double *init_score, double *score) {
    svm_problem prob;
    double mins[DIM_A], maxs[DIM_A];
    for (int j = 0; j < DIM_A; ++j) {
        mins[j] =  std::numeric_limits<double>::infinity();
        maxs[j] = -std::numeric_limits<double>::infinity();
    }

    // split data into training set
    int n_pos = 0, n_neg = 0;
    std::vector<int> train_indices;
    prob.y = new double[size];
    prob.l = 0;
    for (int i = 0; i < size; ++i) {
        double s = init_score[i];
        if (s > UP_PROBA || (s < 1e-4 && s > FLOAT_ZERO)) {
            double* row = params + i*DIM_A;
            train_indices.push_back(i);

            for (int j = 0; j < DIM_A; ++j) {
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
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int j = 0; j < DIM_A; ++j) {
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
    const int total = prob.l * DIM_A;
    for (int i = 0; i < prob.l; ++i) {
        double* row = params + train_indices[i]*DIM_A;
        for (int j = 0; j < DIM_A; ++j) {
            double v = row[j];
            sum += v;
            sqsum += v * v;
        }
    }
    double mean = sum / total;
    double mean_sq = mean * mean;
    double var = (sqsum / total) - mean_sq;
    if (var <= 0) var = FLOAT_ZERO;
    param.gamma = 1.0 / (DIM_A * var);
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
        score[i] = svm_predict_score(model, params + i*DIM_A);
    svm_free_and_destroy_model(&model);

    delete[] prob.x;
    delete[] prob.y;
    return true;
}