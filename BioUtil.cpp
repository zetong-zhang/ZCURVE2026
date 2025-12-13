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

double bio_util::gc_count(
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

double bio_util::gc_fraction(
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
            if (!(gene.partial5 || i)) {
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

static int refine_start(int_array &start_locs, int_array &start_types) {
    int t_start = start_locs.at(0);
    for (int k = 0; k < start_locs.size(); k ++) {
        int alter = start_locs.at(k);
        if (start_types[k] == 0 && (alter-t_start)<=60) {
            t_start = alter;
            break;
        }
    }
    return t_start;
}

static op_type check_overprint(const bio::orf &a, const bio::orf &b) {
    static const int MIN_OLEN = 60;
    if (!strcmp(a.host, b.host) && a.strand == b.strand) {
        int olen = 0;
        if (a.strand == '+' && a.end > b.t_start && b.end > a.t_start)
            olen = std::min(a.end, b.end) - std::max(a.t_start, b.t_start);
        if (a.strand == '-' && a.end < b.t_start && b.end < a.t_start)
            olen = std::min(a.t_start, b.t_start) - std::max(a.end, b.end);
        if (olen == a.len || olen == b.len) return op_type::INCLUDE;
        else if (olen >= MIN_OLEN) return op_type::INTERSECT;
    }
    return op_type::DISJOINT;
}

void bio_util::get_orfs(
    bio::record &scaffold,
    const str_array &starts,
    const str_array &stops,
    const int minlen,
    const bool circ,
    bio::orf_array &orfs
) {
    char *genome = (char *) scaffold.sequence.c_str();
    int length = (int) scaffold.sequence.size();
    char *hostname = (char*)scaffold.name.c_str();
    for (int neg_strand = 0; neg_strand < 2; neg_strand ++) {
        if (neg_strand) {
            genome = get_complement(scaffold.sequence);
            scaffold.complement = genome;
        }
        std::vector<int> phase_first_orfs(3);
        std::vector<int> first_orfs_to_del(0); 
        for (int phase = 0; phase < 3; phase ++) {
            int_array slocs(0), types(0);
            int ps;
            phase_first_orfs.at(phase) = orfs.size();

            // search for standard ORFs
            for (ps = phase; ps < length; ps += 3) {
                int match = match_codon(genome+ps, starts);
                if (match > -1 && slocs.size() < 6) {
                    slocs.push_back(ps);
                    types.push_back(match);
                } else if (slocs.size() && match_codon(genome+ps, stops) > -1) {
                    int end = ps + 3;
                    int t_start = refine_start(slocs, types);
                    int seqlen = end - t_start;
                    if (seqlen >= minlen) {
                        char *pstr = genome + t_start;
                        double gc_frac = gc_fraction(pstr, seqlen);
                        char strand = neg_strand ? '-' : '+';
                        orfs.emplace_back(
                            hostname, length, std::move(slocs), std::move(types), 
                            end, t_start, seqlen, strand, pstr, gc_frac
                        );
                    }
                    slocs.clear(), types.clear();
                }
            }

            bool partial3 = (bool) slocs.size();
            // search for partial 3'-end ORFs
            if (partial3) {
                int end = ps; if (ps > length) end = length;
                int t_start = refine_start(slocs, types);
                int seqlen = end - t_start;
                if (seqlen >= minlen || circ) {
                    char *pstr = genome + t_start;
                    double gc_frac = gc_fraction(pstr, seqlen);
                    char strand = neg_strand ? '-' : '+';
                    orfs.emplace_back(
                        hostname, length, std::move(slocs), std::move(types), 
                        end, t_start, seqlen, strand, pstr, gc_frac
                    );
                    orfs.at(orfs.size()-1).partial3 = true;
                }
            }

            // search for partial 5'-end ORFs
            bool has_start = false;
            int ps0 = (phase + ((3 - length % 3) % 3)) % 3;
            for (int ps = ps0; ps < length; ps += 3) {
                if (!has_start && match_codon(genome+ps, starts) > -1) {
                    has_start = true;
                    continue;
                }
                if (match_codon(genome+ps, stops) > -1) {
                    int end = ps + 3;
                    if (circ && partial3) {
                        bio::orf &orf3 = orfs.at(orfs.size()-1);
                        int seqlen = orf3.len + end;
                        if (seqlen >= minlen) {
                            orf3.seq = new char[seqlen];
                            std::strncpy(orf3.seq, orf3.pstr, orf3.len);
                            std::strncpy(orf3.seq+orf3.len, genome, end);
                            orf3.end = end;
                            orf3.len = seqlen;
                            orf3.pstr = orf3.seq;
                            orf3.gc_frac = gc_fraction(orf3.seq, seqlen);
                            orf3.partial3 = false;
                        } else orfs.pop_back();
                        if (has_start) first_orfs_to_del.push_back(ps0);
                    } else {
                        int t_start = ps0;
                        int seqlen = end - t_start;
                        if (!has_start && seqlen >= minlen) {
                            char *pstr = genome + t_start;
                            double gc_frac = gc_fraction(pstr, seqlen);
                            char strand = neg_strand ? '-' : '+';
                            orfs.emplace_back(
                                hostname, length, int_array(0), int_array(0), 
                                end, t_start, seqlen, strand, pstr, gc_frac
                            );
                            orfs.at(orfs.size()-1).partial5 = true;
                        }
                    }
                    break;
                }
            }
        }

        for (auto index : first_orfs_to_del) {
            orfs.erase(orfs.begin()+phase_first_orfs.at(index));
        }
    }
}