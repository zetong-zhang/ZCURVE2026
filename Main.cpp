/*** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                 *
 *      PROTEIN-CODING GENE RECOGNITION SYSTEM OF ZCURVE 2026      *
 *                                                                 *    
 *      @copyright:   (C) 2003-2026 TUBIC, Tianjin University      *
 *      @author:      Zetong Zhang, Yan Lin, Feng Gao              *
 *      @version:     0.0.6-SNAPSHOT                               *
 *      @date:        2025-11-30                                   *
 *      @license:     GNU GPLv3                                    *
 *      @contact:     ylin@tju.edu.cn | fgao@tju.edu.cn            *
 *                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
    #include <limits.h>
#endif

#include <omp.h>
#include <chrono>

#include "cxxopts.hpp"

#include "BioIO.hpp"
#include "BioUtil.hpp"
#include "Encoding.hpp"
#include "Model.hpp"

/* run the program in quiet mode */
static bool QUIET = false;
/* set the random seed */
static int RAN_SEED = 32522;
/* get the home path */
static fs::path get_home_path();

int main(int argc, char *argv[]) {
    auto start_t = std::chrono::system_clock::now();  // program start time
    auto exe_name = fs::path(argv[0]).filename().string();  // executable file name
    fs::path home_path = get_home_path();  // home path
    std::mt19937_64 ran_eng(RAN_SEED);  // random engine
    
    /* set and parse parameters */
    cxxopts::Options options(exe_name);

    /* general parameters */
    options.add_options("General")
        ("h,help",     "Print help menu and exit.")

        ("q,quiet",    "Run quietly with no stderr output.")

        ("T,threads",  "Number of threads to use. (default: all)",
         cxxopts::value<uint32_t>());
    
    /* input/output parameters */
    options.add_options("Input/Output")
        ("i,input",    "Specify FASTA/Genbank input file. (default: stdin)",
         cxxopts::value<std::string>())
        
        ("o,output",   "Specify output file. (default: stdout)",
         cxxopts::value<std::string>())
        
        ("f,format",   "Select output format (gff, gbk).",
         cxxopts::value<std::string>()->default_value("gff"))
        
        ("a,faa",      "Write protein translations to the selected file.",
         cxxopts::value<std::string>())
 
        ("d,fna",      "Write nucleotide sequences of genes to the selected file.",
         cxxopts::value<std::string>());
    
    /* orf-finder parameters */
    options.add_options("ORFfinder")
        ("g,table",     "Specify a translation table to use.",
         cxxopts::value<uint32_t>()->default_value("11"))

        ("l,minlen",    "Specify the mininum length of ORFs.",
         cxxopts::value<uint32_t>()->default_value("90"))

        ("c,circ",      "Treat topology as circular.");
    
    /* prediction parameters */
    options.add_options("Prediction")
        ("b,bypass",   "Bypass semi-supervised SVM training.")

        ("s,thres",    "Specify putative gene score threshold.",
         cxxopts::value<double>()->default_value("0"));
    
    cxxopts::ParseResult args;
    try {
        args = options.parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\nUse '-h' or '--help' to show help info\n";
        return 1;
    }

    /* show help information and exit with code 0 */
    if (argc <= 1 || args.count("help")) {
        std::cerr << "- - - - - - - - - - - - - - - - - - - - - - - - - - -\n"
                  << "PROTEIN-CODING GENE RECOGNITION SYSTEM OF ZCURVE 2026\n\n"
                  << "Copyright:  (C) 2003-2026 TUBIC, Tianjin University  \n"
                  << "Authors:    Zetong Zhang, Yan Lin, Feng Gao          \n"
                  << "Date:       November 30, 2025                        \n"
                  << "Contact:    ylin@tju.edu.cn | fgao@tju.edu.cn        \n"
                  << "- - - - - - - - - - - - - - - - - - - - - - - - - - -"
                  << options.help() << "\nExample: " << exe_name << ' '
                  << "-i example.fa -o example.gff -c -f gff               \n";
        
        return 0;
    }

    /* set using quiet mode or not */
    if (args.count("quiet")) QUIET = true;
    if (!QUIET) std::cerr << "\nPROTEIN-CODING GENE RECOGNITION SYSTEM OF ZCURVE 2026\n\n";
    
    /* convert time to timestamp for the calculation of used seconds */
    auto start_time_t = std::chrono::system_clock::to_time_t(start_t);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&start_time_t), "%Y-%m-%d %H:%M:%S");
    if (!QUIET) std::cerr << "Program Start Time: " << std::setw(33) << oss.str() << '\n';

    /* set the number of cpus for multi-thread process (openmp) */
    auto num_threads = args.count("threads") ? 
    args["threads"].as<uint32_t>() : omp_get_max_threads();
    omp_set_num_threads(num_threads);
    if (!QUIET) std::cerr << "Threads Used:" << std::setw(40) << num_threads << '\n';
    
    /* initialize translation table and start/stop codon list */
    uint32_t table_code = args["table"].as<uint32_t>();
    str_array starts{ "ATG", "GTG", "TTG" };
    str_array stops{ "TAA", "TAG", "TGA" };

    if (table_code == 1) {
        starts.assign({ "ATG "});
    } else if (table_code == 4) {
        stops.assign({ "TAA", "TAG" });
        trans_tbl[14] = 'W';
    } else if (table_code != 11) {
        std::cerr << "\nError: unsupported translation table " << table_code 
                  << " (options: 1 (standard), 4 (mycoplasma, spiroplasma), 11 (bacteria, archaea))\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Translation Table:" << std::setw(35) << table_code << '\n';

    /* set mininum gene length */
    auto minlen = args["minlen"].as<uint32_t>();
    if (minlen < 30) {
        std::cerr << "\nError: mininum gene length shorter than 30\n";
        return 1;
    }
    if (!QUIET) std::cerr << "Mininum Gene Length:" << std::setw(30) << minlen << " nt\n";

    /* check output format */
    auto format = args["format"].as<std::string>();
    if (format != "gff" && format != "gbk") {
        std::cerr << "\nError: unknown output format (options: gff (GFF3), gbk (GenBank))\n";
        return 1;
    }

    /* handle input */
    std::string input = "-";
    if (args.count("input")) input = args["input"].as<std::string>();
    bio::record_array scaffolds(0);
    if (!bio_io::read_source(input, scaffolds)) {
        if (input == "example.fa") {
            auto example = (home_path / "example.fa").string();
            if(bio_io::read_source(example, scaffolds)) goto CHECK;
        } return 1;
    }
    CHECK: if (!scaffolds.size()) {
        std::cerr << "Error: no valid sequence was read\n";
        return 0;
    }
    if (!QUIET) std::cerr << "Number of Scaffolds: " << std::setw(32) << scaffolds.size() << '\n';

    size_t total_len = 0;
    double gc_cont = 0.0;

    /* extract all the orfs */
    bool circ = (bool) args.count("circ");
    bio::orf_array orfs;
    for (int i = 0; i < scaffolds.size(); i ++) {
        const auto sublen = scaffolds[i].sequence.length();
        total_len += sublen;
        gc_cont += bio_util::gc_count(scaffolds[i].sequence.c_str(), sublen);

        if (sublen >= minlen) {
            bio_util::get_orfs(scaffolds[i], starts, stops, minlen, circ, orfs);
        } else if (!QUIET) {
            std::cerr << "Warning: skip " << scaffolds[i].name << " (too short)\n";
        }
    }
    const uint32_t n_orfs = (uint32_t) orfs.size();

    /* calculate basic genome stats */
    gc_cont /= total_len;
    if (!QUIET) {
        std::cerr << "Genome Size: " << std::setw(37) << total_len << " bp\n";
        std::cerr << "G+C Content: " << std::setw(38) << std::fixed << std::setprecision(2) 
                  << gc_cont * 100 << " % \n";
        std::cerr << "Number of ORFs: " << std::setw(37) << n_orfs << '\n';
    }

    /* sort orfs by G+C content */
    std::sort(orfs.begin(), orfs.end(), [](bio::orf& a, bio::orf& b) { return a.gc_frac < b.gc_frac; });
    int gc_intv_count[N_MODELS+1] = {0};
    {
        double max_gc = 0.20;
        for (int i = 0, j = 0; i < n_orfs; i ++) {
            if (orfs.at(i).gc_frac > max_gc) {
                max_gc += 0.01;
                if ((++j) >= (N_MODELS+1)) break;
            }
            gc_intv_count[j] ++;
        }
    }
    
    /* convert orfs to zcurve params */
    double *params = new double[DIM_A*n_orfs];
    encoding::encode_orfs(orfs, params);

    /* load pre-trained models and calculate scores */
    if (!model::init_models(home_path)) return 1;
    double *probas = new double[n_orfs]();
    int off = 0, num_seeds = 0;
    if (!QUIET) std::cerr << "Initialization:" << std::setw(38) << "0 %";
    for (int i = 0; i < N_MODELS; i ++) {
        int size = gc_intv_count[i];
        model::mlp_predict(i, params+off*DIM_A, size, probas+off);
        if (!QUIET) std::cerr << "\rInitialization:" << std::setw(36) << (int)(i*1.67) << " %";
        off += gc_intv_count[i];
    }
    int_array seeds;
    for (int i=0;i<n_orfs;i++) if (probas[i]>UP_PROBA) { seeds.push_back(i); num_seeds ++; }
    if (!QUIET) std::cerr << "\nNumber of Seed ORFs: " << std::setw(32) << seeds.size() << "\n";
    
    /* train rbs-svm models */
    if (!QUIET) std::cerr << "Training Model ...";
    double *scores = new double[n_orfs]();
    bool training = !((bool) args.count("bypass"));
    if (training) {
        training = model::train_predict(params, n_orfs, probas, scores);
        if (!QUIET) std::cerr << std::setw(36) << (training ? "Done\n" : "Skipped\n");
    } else if (!QUIET) std::cerr << std::setw(36) << "Bypassed\n";
    
    /* classifying orfs */
    auto thres = args["thres"].as<double>();
    int num_putative = 0;
    for (int i = 0; i < n_orfs; i ++) {
        if (!training || (scores[i] < thres && probas[i] >= 0.5)) orfs[i].score = probas[i] - 0.5;
        else orfs[i].score = scores[i];
        if (orfs[i].score > thres) num_putative ++;
    }
    bio::orf_array putative(num_putative);
    std::copy_if(orfs.begin(), orfs.end(), putative.begin(), [thres](const bio::orf& orf) { 
        return orf.score > thres; 
    });
    if (!QUIET) std::cerr << "Number of Putative Genes:" << std::setw(28) << num_putative << "\n";

    /* write results to output */
    std::string output = "-";
    if (args.count("output")) output = args["output"].as<std::string>();
    std::sort(putative.begin(), putative.end(), [](bio::orf& a, bio::orf& b) {
        int host_cmp = std::strcmp(a.host, b.host);
        if (host_cmp < 0) return true;
        else if (host_cmp > 0) return false;
        else {
            int a_end = (a.strand=='+')?a.end:a.host_len-a.t_start;
            int b_end = (b.strand=='+')?b.end:b.host_len-b.t_start;
            return a_end < b_end;
        }
    });
    bio_io::write_result(putative, output, format);

    /* write protein sequences */
    if (args.count("faa")) {
        auto faa = args["faa"].as<std::string>();
        if(!bio_io::write_faa(putative, faa)) return 1;
    }

    /* write nucleotide sequences */
    if (args.count("fna")) {
        auto fna = args["fna"].as<std::string>();
        if(!bio_io::write_fna(putative, fna)) return 1;
    }

    if (!QUIET) {
        auto end_t = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t);
        auto seconds = duration.count() / 1000.0;
        std::cerr << "\nFinished in " << std::fixed << std::setprecision(3) << seconds << " s\n";
    }
}

static fs::path get_home_path() {
#ifdef _WIN32
    wchar_t buffer[MAX_PATH];
    GetModuleFileNameW(NULL, buffer, MAX_PATH);
    return fs::path(buffer).parent_path();
#else
    char buffer[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer)-1);
    if (len != -1) {
        buffer[len] = '\0';
        return fs::path(buffer).parent_path();
    }
    return fs::current_path();
#endif
}