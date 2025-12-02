/** 
 * @brief       Bioinformatics I/O operations. 
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.0.4-SNAPSHOT
 * @date        2025-11-30
 * @license     GNU GPLv3
 * @contact     fgao@tju.edu.cn
 */

#ifndef BIO_IO
#define BIO_IO

#include <iostream>
#include <fstream>
#include <sstream>
#ifdef ZLIB
#include <zlib.h>
#endif
#include <assert.h>
#include <iomanip>

#include "BioStruct.hpp"
#include "BioUtil.hpp"
/* The buffer size for file reading. */
#define BUFF_SIZE 65536
/* The version info of the software. */
#define VERSION  "ZCURVE_2026"
/* The namespace for bioinformatics I/O operations. */
namespace bio_io {
    /** 
     * @brief   Read a FASTA or GenBank file from a stream.
     * 
     * @param filename  The path to the file.
     * @param scaffolds The record array to store the scaffolds.
     * @return true If the file is successfully read.
     *         false If the file is not successfully read.
     */
    bool read_source(const std::string &filename, bio::record_array &scaffolds);
    /** 
     * @brief   Write the ORFs to a result file.
     * 
     * @param orfs      The ORF array to be written.
     * @param filename  The path to the output file.
     * @param format    The format of the output file.
     * @param thres     The threshold score for ORF selection.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_result(bio::orf_array &orfs, const std::string &filename, 
                      const std::string &format, const double thres);
    /** 
     * @brief   Write the ORFs to a FASTA file.
     * 
     * @param orfs      The ORF array to be written.
     * @param filename  The path to the output file.
     * @param thres     The threshold score for ORF selection.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_faa(bio::orf_array &orfs, const std::string &filename, const double thres);
    /** 
     * @brief   Write the ORFs to a FASTA file.
     * 
     * @param orfs      The ORF array to be written.
     * @param filename  The path to the output file.
     * @param thres     The threshold score for ORF selection.
     * @return true If the file is successfully written.
     *         false If the file is not successfully written.
     */
    bool write_fna(bio::orf_array &orfs, const std::string &filename, const double thres);
}

#endif