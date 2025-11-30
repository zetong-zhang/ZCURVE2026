#include "BioUtil.hpp"
/* complement base table */
static char COMP[] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 
    'N',
  /* A    B    C    D    E    F    G    H  */
    'T', 'V', 'G', 'H', 'N', 'N', 'C', 'D', 
  /* I    J    K    L    M    N    O    P  */
    'N', 'N', 'M', 'N', 'K', 'N', 'N', 'N',
  /* Q    R    S    T    U    V    W    X  */
    'N', 'Y', 'W', 'A', 'A', 'B', 'S', 'N', 
  /* Y    Z */
    'R', 'T', 'N', 'N', 'N', 'N', 'N', 'N',
  /* a    b    c    d    e    f    g    h  */
    't', 'v', 'g', 'h', 'n', 'n', 'c', 'd', 
  /* i    j    k    l    m    n    o    p  */
    'n', 'n', 'm', 'n', 'k', 'n', 'n', 'n',
  /* q    r    s    t    u    v    w    x  */
    'n', 'y', 'w', 'a', 'a', 'b', 's', 'n', 
  /* y    z */
    'r', 't'
};
/*  base to offset map */
static const std::map<char, int> base2offset = {
    { 'T', 0 }, { 'C', 1 }, { 'A', 2 }, { 'G', 3 },
    { 't', 0 }, { 'c', 1 }, { 'a', 2 }, { 'g', 3 }
};
/*  codon to amino acid map */
char trans_tbl[64] = {
    'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S',
    'Y', 'Y', '*', '*', 'C', 'C', '*', 'W',
    'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P',
    'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 
    'I', 'I', 'I', 'M', 'T', 'T', 'T', 'T',
    'N', 'N', 'K', 'K', 'S', 'S', 'R', 'R', 
    'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A',
    'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'
};

float bio_util::gc_count(
    const char *pstr,
    const size_t len
) {
    float count = 0.0F;
    for (size_t i = 0; i < len; i ++) {
        char c = std::tolower(pstr[i]);
        count += (int)(c == 'g' || c == 'c');
    }
    return count;
}

float bio_util::gc_fraction(
    const char *pstr,
    const size_t len
) { 
    return gc_count(pstr, len) / len; 
}

char *bio_util::get_complement(
    const std::string &origin
) {
    size_t len = origin.size();
    char *comp = new char[len];
    if (!comp) return nullptr;
    for (size_t i = 0; i < len; i ++)
        comp[i] = COMP[origin.at(len - i - 1)];
    return comp;
}

static int match_codon(
    const char *origin,
    const str_array &codon_lst
) {
    int num_codons = (int) codon_lst.size();
    for (int index = 0; index < num_codons; index ++) {
        const char *pstr = codon_lst.at(index).c_str();
        if (!std::strncmp(pstr, origin, 3))
            return (int) index;
    }
    return -1;
}

char *bio_util::gene2protein(bio::orf &gene, int prolen) {
    char *protein = new char[prolen]; 
    try {
        for (int i = 0; i < prolen; i ++) {
            if (!(gene.edge || i)) {
                protein[0] = 'M';
                continue;
            }
            int f = base2offset.at(gene.pstr[i*3+0]);
            int s = base2offset.at(gene.pstr[i*3+1]);
            int t = base2offset.at(gene.pstr[i*3+2]);
            protein[i] = trans_tbl[f*16+s*4+t];
        }
    } catch (const std::out_of_range& e) {
        delete[] protein;
        return nullptr;
    }
    return protein;
}

void bio_util::get_orfs(
    bio::record &scaffold,
    const str_array &starts,
    const str_array &stops,
    const int minlen,
    bio::orf_array &orfs
) {
    char *genome = (char *) scaffold.sequence.c_str();
    int length = (int) scaffold.sequence.size();
    for (int neg_strand = 0; neg_strand < 2; neg_strand ++) {
        if (neg_strand) {
            genome = get_complement(scaffold.sequence);
            scaffold.complement = genome;
        }
        for (int phase = 0; phase < 3; phase ++) {
            int_array slocs(0), types(0);
            int ps;

            for (ps = phase; ps < length; ps += 3) {
                int match = match_codon(genome+ps, starts);
                if (match > -1 && slocs.size() < 6) {
                    slocs.push_back(ps);
                    types.push_back(match);
                } else if (slocs.size() && match_codon(genome+ps, stops) > -1) {
                    int end = ps + 3;
                    int t_start = slocs.at(0);
                    for (int k = 0; k < types.size(); k ++) {
                        int alter = slocs.at(k);
                        if (types[k] == 0 && (alter-t_start)<=30) {
                            t_start = alter;
                            break;
                        }
                    }
                    int seqlen = end - t_start;
                    if (seqlen >= minlen) {
                        char *pstr = genome + t_start;
                        float gc_frac = gc_fraction(pstr, seqlen);
                        char strand = neg_strand ? '-' : '+';
                        if (neg_strand) {
                            end = length - end;
                            t_start = length - t_start;
                            for (int i = 0; i < slocs.size(); i ++)
                                slocs.at(i) = length - slocs.at(i);
                        }
                        orfs.emplace_back(
                            (char*)scaffold.name.c_str(), std::move(slocs), std::move(types), 
                            end, t_start, seqlen, strand, pstr, gc_frac
                        );
                    }
                    slocs.clear(), types.clear();
                }
            }

            if (slocs.size()) {
                int end = ps; if (ps > length) end -= 3;
                int t_start = slocs.at(0);
                for (int k = 0; k < types.size(); k ++) {
                    int alter = slocs.at(k);
                    if (types[k] == 0 && (alter-t_start)<=30) {
                        t_start = alter;
                        break;
                    }
                }
                int seqlen = end - t_start;
                if (seqlen >= minlen) {
                    char *pstr = genome + t_start;
                    float gc_frac = gc_fraction(pstr, seqlen);
                    char strand = neg_strand ? '-' : '+';
                    if (neg_strand) {
                        end = length - end;
                        t_start = length - t_start;
                        for (int i = 0; i < slocs.size(); i ++)
                            slocs.at(i) = length - slocs.at(i);
                    }
                    orfs.emplace_back(
                        (char*)scaffold.name.c_str(), std::move(slocs), std::move(types), 
                        end, t_start, seqlen, strand, pstr, gc_frac
                    );
                    orfs[orfs.size()-1].edge = true;
                }
            }
        }
    }
}

void  bio_util::get_edge_orfs(
    bio::record &scaffold,
    const str_array &starts,
    const str_array &stops,
    const int minlen,
    bio::orf_array &orfs
) {
    char *genome = (char *) scaffold.sequence.c_str();
    int length = (int) scaffold.sequence.size();
    for (int neg_strand = 0; neg_strand < 2; neg_strand ++) {
        if (neg_strand) {
            genome = scaffold.complement;
            if (!genome) return;
        }
        for (int phase = 0; phase < 3; phase ++) {
            bool has_start = false;
            for (int ps = phase; ps < length; ps += 3) {
                if (!has_start && match_codon(genome+ps, starts) > -1) {
                    has_start = true;
                    continue;
                }
                if (match_codon(genome+ps, stops) > -1) {
                    int end = ps + 3, t_start = phase;
                    int seqlen = end - t_start;
                    if (!has_start && seqlen >= minlen) {
                        char *pstr = genome + t_start;
                        float gc_frac = gc_fraction(pstr, seqlen);
                        char strand = neg_strand ? '-' : '+';
                        if (neg_strand) {
                            end = length - end;
                            t_start = length - t_start;
                        }
                        orfs.emplace_back(
                            (char*)scaffold.name.c_str(), int_array(0), int_array(0), 
                            end, t_start, seqlen, strand, pstr, gc_frac
                        );
                        orfs[orfs.size()-1].edge = true;
                    }
                    break;
                }
            }
        }
    }
}