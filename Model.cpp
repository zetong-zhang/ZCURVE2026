#include "Model.hpp"

static const int MIN_SET_SIZE = 30;

static svm_parameter param = {
    C_SVC,   /* svm_type     */
    RBF,     /* kernel_type  */
    3,       /* degree       */
    1/DIM,   /* gamma        */
    0.0,     /* coef0        */
    200,     /* cache_size   */
    1e-4,    /* eps          */
    1.0,     /* C            */
    0,       /* nr_weight    */
    nullptr, /* weight_label */
    nullptr, /* weight       */
    0.0,     /* nu           */
    0.0,     /* p            */
    true,    /* shrinking    */
    false    /* probability  */
};

float MODELS[N_MODELS][N_PARAMS];

inline double dot_product_avx(const double* x, const float* y, const int dim) {
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
}

bool model::init_models(const fs::path model_path) {
    fs::path models_path = model_path / "meta.bin";
    std::ifstream file(models_path, std::ios::binary);
	if (!file.is_open()) {
		std::cerr << "\nError: unable to open " << models_path << '\n';
		return false;
	}
    int chunk_size = N_PARAMS*sizeof(float);
    for (int i = 0; i < N_MODELS; i ++) {
    	file.read((char*)MODELS[i], chunk_size);
		if (!file || file.gcount() != chunk_size) {
			std::cerr << "\nError: data corruption in " << models_path << "\n";
			return false;
		}
	}
    return true;
}

void model::mlp_predict(int index, double *data, int size, double *probas) {
    if (!size) return;
    /* scale data */
    float *means = MODELS[index], *stds = means + DIM;
    double *scaled = encoding::std_trans(data, size, DIM, means, stds);
    /* mlp process */
    float *model = stds + DIM; // first hidden w
    float *hid_bs = model + N_HIDDEN*DIM;  // first hidden b
    float *out_ws = hid_bs + N_HIDDEN;  // output w
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < size; i ++) {
        double *x = scaled + i*DIM;
        // Hidden Layer (765, ) -> (200, )
        double hid_out[N_HIDDEN] = {0.0};
        for (int j = 0; j < N_HIDDEN; j ++) {
            float *hid_w = model + j*DIM;
            // hid_out = hid_w * x + hid_b
            hid_out[j] = dot_product_avx(x, hid_w, DIM);
            // ReLU activation function
            hid_out[j] = std::max(0.0, hid_out[j]+(double)hid_bs[j]);
        }
        // Output Layer (200, ) -> (1, )
        double out_b = (double)model[N_PARAMS-1];
        // out = out_w * hid_out + out_b
        probas[i] = dot_product_avx(hid_out, out_ws, N_HIDDEN);
        // Sigmoid activation function
        probas[i] = 1.0/(1.0 + std::exp(-(probas[i])));
    }
    delete[] scaled;
}

bool model::train_predict(double *params, int size, double *init_score, double *score) {
    svm_problem prob;

    double mins[DIM], maxs[DIM];
    for (int j = 0; j < DIM; ++j) {
        mins[j] =  std::numeric_limits<double>::infinity();
        maxs[j] = -std::numeric_limits<double>::infinity();
    }

    int n_pos = 0, n_neg = 0;
    std::vector<int> train_indices;
    prob.y = new double[size];
    prob.l = 0;
    for (int i = 0; i < size; ++i) {
        double s = init_score[i];
        if (s > 0.9 || s < 1e-4) {
            double* row = params + i*DIM;
            train_indices.push_back(i);

            for (int j = 0; j < DIM; ++j) {
                double v = row[j];
                if (v < mins[j]) mins[j] = v;
                if (v > maxs[j]) maxs[j] = v;
            }
            if (s > 0.9) { prob.y[prob.l] = 1.0; n_pos ++; } 
            else { prob.y[prob.l] = -1.0; n_neg ++; }
            ++prob.l;
        }
    }

    if (n_pos < MIN_SET_SIZE || n_neg < MIN_SET_SIZE) return false; 
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int j = 0; j < DIM; ++j) {
        double intv = maxs[j] - mins[j];
        if (intv == 0.0) {
            for (int i = 0; i < size; ++i)
                params[i*DIM + j] = 0.0;
        } else {
            double inv_intv = 1.0 / intv;
            double minv = mins[j];

            for (int i = 0; i < size; ++i) {
                double* row = params + i*DIM;
                row[j] = (row[j] - minv) * inv_intv;
            }
        }
    }

    double sum = 0.0;
    double sqsum = 0.0;
    const int total = prob.l * DIM;
    for (int i = 0; i < prob.l; ++i) {
        double* row = params + train_indices[i]*DIM;
        for (int j = 0; j < DIM; ++j) {
            double v = row[j];
            sum += v;
            sqsum += v * v;
        }
    }

    double mean = sum / total;
    double mean_sq = mean * mean;
    double var = (sqsum / total) - mean_sq;
    if (var <= 0) var = 1e-12;

    param.gamma = 1.0 / (DIM * var);

    prob.x = new double *[prob.l];
    for (int i = 0; i < prob.l; i ++) {
        prob.x[i] = params + train_indices[i]*DIM;
    }

    svm_model *model = svm_train(&prob, &param);
    if (!model) return false;

#ifdef _OPENMP
    #pragma omp parallel for schedule(guided)
#endif
    for (int i = 0; i < size; i ++) {
        score[i] = svm_predict_score(model, params + i*DIM);
    }
    svm_free_and_destroy_model(&model);

    delete[] prob.x;
    delete[] prob.y;
    return true;
}