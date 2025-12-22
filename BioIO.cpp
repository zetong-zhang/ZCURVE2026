#include "BioIO.hpp"

// static function declarations

/**
 * @brief   Check if a file is gzipped.
 * 
 * @param filename  The path to the file.
 * @return true If the file is gzipped.
 *         false If the file is not gzipped.
 */
static bool is_gzip(const std::string &filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) return false;
    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);
    // 1st and 2nd magic of gzip are 0x1F and 0x8B
    return (magic[0] == 0x1F && magic[1] == 0x8B);
}

/**
 * @brief   Read a FASTA or GenBank file from a stream.
 * 
 * @param handle    The input stream.
 * @param scaffolds The record array to store the scaffolds.
 * @return true If the file is successfully read.
 *         false If the file is not successfully read.
 */
static bool read_stream(
    std::istream &handle,
    bio::record_array &scaffolds
) {
    std::string header;
    // skip empty lines
    while(std::getline(handle, header))
        if (!header.empty()) break;
    if (header.empty()) {
        std::cerr << "\nError: empty content \n";
        return false;
    }
    std::string line, buffer;
    // FASTA files should start with '>'
    if (header.at(0) == '>') {
        while (std::getline(handle, line)) {
            if (line.empty()) continue;
            else if (line.at(0) == '>') {
                size_t sep = header.find(' ');
                std::string entry = sep>1 ? header.substr(1, sep-1):header.substr(1);
                if (!entry.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    header = std::move(line);
                    buffer.clear();
                } else {
                    std::cerr << "\nError: invalid FASTA header \n";
                    return false;
                }
            } else {
                if (!std::isalpha(line.back())) line.pop_back();
                buffer.append(std::move(line));
            }
            line.clear();
        }
        if (!buffer.empty()) {
            size_t sep = header.find(' ');
            std::string entry = sep>1 ? header.substr(1, sep-1):header.substr(1);
            scaffolds.emplace_back(std::move(entry), std::move(buffer));
        }
        return true;
    // GenBank files should start with 'LOCUS ... '
    } else if (!header.find("LOCUS")) {
        bool in_origin = false;
        std::string entry;
        std::istringstream header_stream(header);
        header_stream >> entry;
        if (!(header_stream >> entry)) {
            std::cerr << "\nError: invalid LOCUS line \n";
            return false;
        }
        while (std::getline(handle, line)) {
            if (line.empty()) continue;
            else if (!line.find("ORIGIN")) in_origin = true;
            else if (!line.find("//")) {
                if (!buffer.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    entry.clear(), buffer.clear();
                }
                in_origin = false;
            } else if (!line.find("LOCUS")) {
                if (!buffer.empty()) {
                    scaffolds.emplace_back(std::move(entry), std::move(buffer));
                    entry.clear(), buffer.clear();
                }
                in_origin = false; 

                std::istringstream header_stream(line);
                header_stream >> entry;
                if (!(header_stream >> entry)) {
                    std::cerr << "\nError: invalid LOCUS line \n";
                    return false;
                }
            } else if (in_origin) {
                for (char c : line) if (std::isalpha(c)) buffer += std::toupper(c);
            }
        }
        if (!buffer.empty()) {
            scaffolds.emplace_back(std::move(entry), std::move(buffer));
        }
        return true;
    } else {
        std::cerr << "\nError: unsupported file type (options: GenBank, FASTA)\n";
        return false;
    }
}
/**
 * @brief   Write the ORFs to a result file.
 * 
 * @param handle    The output stream.
 * @param orfs      The ORF array to be written.
 * @param format    The format of the output file.
 * @return true If the file is successfully written.
 *         false If the file is not successfully written.
 */
static bool write_stream(
    std::ostream &handle,
    bio::orf_array &orfs,
    const std::string &format
) {
    int count = (int) orfs.size();
    if (format == "gff") {
        handle << "##gff-version 3 \n";
        for (int i = 0, j = 0; i < count; i ++) {
            // seqid + source + type
            handle << orfs[i].host << '\t' << VERSION << "\tCDS\t";
            int rstart, rend, host_len = orfs[i].host_len;
            if (orfs[i].strand == '+') {
                rstart = orfs[i].t_start + 1;
                rend = orfs[i].end;
            } else {
                rstart = host_len - orfs[i].end + 1;
                rend = host_len - orfs[i].t_start;
            }
            // start + end
            handle << rstart << '\t' << rend << '\t';
            // score
            handle << std::fixed << std::setprecision(3) 
                   << orfs[i].score << '\t';
            // strand + phase
            handle << orfs[i].strand << "\t0\tID=orf";
            // attributes
            handle << std::setw(6) << std::setfill('0') << (++j);
            if (orfs[i].partial3) handle << ";partial3=true\n";
            else if (orfs[i].partial5) handle << ";partial5=true\n";
            else handle << "\n";
        }
        return true;
    } else if (format == "gbk") {
        char *last_scaffold = nullptr;
        for (int i = 0, j = 0; i < count; i ++) {
            if (!last_scaffold || std::strcmp(orfs[i].host, last_scaffold)) {
                if (last_scaffold) handle << "ORIGIN\n//\n";
                handle << "LOCUS       " << orfs[i].host << "\n";
                handle << "DEFINITION  " << orfs[i].host << "\n";
                last_scaffold = orfs[i].host;
                handle << "FEATURES             Location/Qualifiers\n";
            }
            bool neg_strand = (bool)(orfs[i].strand == '-');
            handle << "     CDS             ";
            int rstart, rend, host_len = orfs[i].host_len;
            if (!neg_strand) {
                rstart = orfs[i].t_start + 1;
                rend = orfs[i].end;
            } else {
                rstart = host_len - orfs[i].end + 1;
                rend = host_len - orfs[i].t_start;
                handle << "complement(";
            }
            if (rstart > rend) {
                handle << "join(" << rstart << ".." 
                << orfs[i].host_len << ",1.." << rend << ')';
            } else {
                if (orfs[i].partial5 && !neg_strand || 
                    orfs[i].partial3 && neg_strand) handle << '<';
                handle << rstart << "..";
                if (orfs[i].partial3 && !neg_strand || 
                    orfs[i].partial5 && neg_strand) handle << '>';
                handle << rend;
            }
            if (neg_strand) handle << ")";
            handle << "\n                     /note=\"version=" << VERSION 
                   << ";ID=orf" << std::setw(6) << std::setfill('0') << (++j)
                   << ";score=" << std::fixed << std::setprecision(3) 
                   << orfs[i].score << "\"\n";
        }
        handle << "ORIGIN\n//\n";
        return true;
    }
    std::cerr << "Error: unsupported output format '" 
              << format << "' (options: gff, gbk)\n";
    return false;
}
// API functions 

bool bio_io::read_source(
    const std::string &filename,
    bio::record_array &scaffolds
) {
    if (filename == "-") return read_stream(std::cin, scaffolds);
#ifdef ZLIB
    else if (is_gzip(filename)) {
        gzFile file = gzopen(filename.c_str(), "rb");
        assert(file != nullptr);
        char buff[BUFF_SIZE];
        int bytes_read, err_code;
        std::string content;
        while ((bytes_read = gzread(file, buff, BUFF_SIZE)))
            content.append(buff, bytes_read);
        const char *err_msg = gzerror(file, &err_code);
        if (err_code != Z_OK && err_code != Z_STREAM_END) {
            gzclose(file);
            std::cerr << "\nError: failed to decompress " << filename << "\n";
            return false;
        }
        gzclose(file);
        std::istringstream handle(std::move(content));
        return read_stream(handle, scaffolds);
    }  
#endif
    else {
        std::ifstream handle(filename);
        if (!handle.is_open()) {
            if (filename != "example.fa") 
                std::cerr << "\nError: cannot open " << filename << "\n";
            return false;
        }
        return read_stream(handle, scaffolds);
    }
}

bool bio_io::write_result(
    bio::orf_array &orfs, 
    const std::string &filename,
    const std::string &format
) {
    if (filename == "-") return write_stream(std::cout, orfs, format);
    else {
        std::ofstream handle(filename);
        if (!handle.is_open()) {
            std::cerr << "\nError: failed to open " << filename << std::endl;
            return false;
        }
        return write_stream(handle, orfs, format);
    }
}

bool bio_io::write_faa(
    bio::orf_array &orfs, 
    const std::string &filename
) {
    std::ofstream handle(filename);
    if (!handle.is_open()) {
        std::cerr << "\nError: failed to open " << filename << std::endl;
        return false;
    }
    int count = (int) orfs.size();
    for (int i = 0; i < count; i ++) {
        int rstart, rend, host_len = orfs[i].host_len;
        if (orfs[i].strand == '+') {
            rstart = orfs[i].t_start + 1;
            rend = orfs[i].end;
        } else {
            rstart = host_len - orfs[i].end + 1;
            rend = host_len - orfs[i].t_start;
        }
        handle << '>' << orfs[i].host << ':' << rstart << ".." << rend << '(' 
               << orfs[i].strand << ") score=" << orfs[i].score << '\n';
        int prolen = orfs[i].len / 3;
        char *protein = bio_util::gene2protein(orfs[i], prolen);
        if (!protein) {
            std::cerr << "Warning: invalid characters encountered in "
                      << orfs[i].host << ':' << rstart 
                      << ".." << rend << '(' << orfs[i].strand << ")\n";
            continue;
        }
        for (int j = 0; j < prolen; j += 80) {
            int line_length = std::min(80, prolen - j);
            handle.write(protein + j, line_length);
            handle << '\n';
        }
        delete[] protein;
    }
    handle.close();
    return true;
}

bool bio_io::write_fna(
    bio::orf_array &orfs, 
    const std::string &filename
) {
    std::ofstream handle(filename);
    if (!handle.is_open()) {
        std::cerr << "\nError: failed to open " << filename << std::endl;
        return false;
    }
    int count = (int) orfs.size();
    for (int i = 0; i < count; i ++) {
        int rstart, rend, host_len = orfs[i].host_len;;
        if (orfs[i].strand == '+') {
            rstart = host_len - orfs[i].t_start + 1;
            rend = host_len - orfs[i].end;
        } else {
            rstart = orfs[i].end + 1;
            rend = orfs[i].t_start;
        }
        handle << '>' << orfs[i].host << ':' << rstart << ".." << rend << '(' 
               << orfs[i].strand << ") score=" << orfs[i].score << '\n';
        for (int j = 0; j < orfs[i].len; j += 80) {
            int line_length = std::min(80, orfs[i].len - j);
            handle.write(orfs[i].pstr + j, line_length);
            handle << '\n';
        }
    }
    handle.close();
    return true;
}