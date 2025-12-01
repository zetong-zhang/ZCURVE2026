/**
 * @brief   Utility functions for bioinformatics.
 * 
 * @author      Zetong Zhang, Yan Lin, Feng Gao
 * @version     0.0.3-SNAPSHOT
 * @date        2025-11-30
 * @license     GNU GPLv3
 * @contact     fgao@tju.edu.cn
 */
#ifndef BIO_UTIL
#define BIO_UTIL

#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "BioStruct.hpp"

/*  codon to amino acid map */
extern char trans_tbl[64];

namespace bio_util {
    /**
     * @brief   count GC number in a nucleotide sequence.
     * 
     * @param   seq nucleotide sequence.
     * @param   len length of the sequence.
     * @return  float GC number.
     */
    float gc_count(const char *seq, const size_t len);
    /**
     * @brief   calculate GC fraction in a nucleotide sequence.
     * 
     * @param   seq nucleotide sequence.
     * @param   len length of the sequence.
     * @return  float GC fraction.
     */
    float gc_fraction(const char *seq, const size_t len);
    /**
     * @brief   get the complement sequence of a nucleotide sequence.
     * 
     * @param   seq nucleotide sequence.
     * @return  char* complement sequence.
     */
    char *get_complement(const std::string &seq);
    /**
     * @brief   get ORFs in a nucleotide sequence.
     * 
     * @param   record  scaffold record.
     * @param   starts  start codon array.
     * @param   stops   stop codon array.
     * @param   circ    treat as circular.
     * @param   minlen  minimum length of the ORF.
     * @param   orfs    ORF array to be written.
     */
    void  get_orfs(bio::record &record, const str_array &starts, const str_array &stops, 
        const int minlen, const bool circ, bio::orf_array &orfs);
    /**
     * @brief   translate a nucleotide sequence to a protein sequence.
     * 
     * @param   orf    ORF to be translated.
     * @param   prolen length of the protein sequence.
     * @return  char* protein sequence.
     */
    char *gene2protein(bio::orf &orf, int prolen);
}

#endif