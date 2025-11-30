/** 
 * @brief       Bioinformatics data structures. 
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.0.3-SNAPSHOT
 * @date        2025-11-30
 * @license     GNU GPLv3
 * @contact     fgao@tju.edu.cn
 */
#ifndef BIO_STRUCT
#define BIO_STRUCT

#include <string>
#include <vector>
#include <filesystem>

/* utility array types. */
typedef std::vector<int>          int_array;
typedef std::vector<std::string>  str_array;
typedef std::vector<char *>       pch_array;

/* bioinformatics data structures */
namespace bio {
    /*  scaffold record type */
    typedef struct record {
        std::string name;            // scaffold name
        std::string sequence;        // scaffold sequence
        char *complement = nullptr;  // complement sequence
        record(): name("anonymous") {}
        record(std::string &&name, std::string &&sequence):
        name(std::move(name)), sequence(std::move(sequence)) {}
        ~record() { delete[] complement; }
    } record;
    /*  ORF type */
    typedef struct orf {
        char *     host;        // host scaffold name
        int_array  starts;      // alter start positions
        int_array  types;       // alter start types 
        int        t_start;     // true start position
        int        end;         // true end position
        int        len;         // length
        char       strand;      // strand direction
        char *     pstr;        // nucleotide sequence
        float      gc_frac;     // GC fraction
        double     score=0;     // zcurve score
        bool       edge=false;  // at edge or not
        orf(): host((char*)"anonymous") {}
        orf(char *host, int_array &&starts, int_array &&types, int end, 
            int t_start, int len, char strand, char *pstr, float gc_frac):
        host(host),starts(std::move(starts)),types(std::move(types)), len(len), 
        t_start(t_start),end(end),strand(strand),pstr(pstr),gc_frac(gc_frac){}
    } orf;

    /*  scaffold record array type */
    typedef std::vector<record> record_array;
    /*  ORF array type */
    typedef std::vector<orf>    orf_array;
}

namespace fs = std::filesystem;

#endif