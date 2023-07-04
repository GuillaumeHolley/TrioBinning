#include <iostream>
#include <sstream>
#include <set>

#include "ColoredCDBG.hpp"
#include "PairID.hpp"

using namespace std;

struct TrioBinning_Opt {

	vector<string> in_fn_read_parentA;
    vector<string> in_fn_read_parentB;
    vector<string> in_fn_read_proband;

    vector<string> in_fn_asm_parentA;
    vector<string> in_fn_asm_parentB;
    vector<string> in_fn_asm_proband;

    string out_fn;

    string out_graph_fn;
    string in_graph_fn;

    int k;

    int min_cov;
    int max_cov;

    int max_len_uniq_run;
    int min_len_shared_run;

    double min_qual_score;
    double threshold_partioning;

    size_t nb_threads;

    bool no_seq_out;
    bool force_partioning;
    bool verbose;

    bool phase;
    bool hom;

	TrioBinning_Opt() : k(31), nb_threads(1), min_cov(2), max_cov(-1), max_len_uniq_run(97), min_len_shared_run(1000),
                        min_qual_score(0.0), threshold_partioning(0.8), no_seq_out(false), force_partioning(false), verbose(false),
                        phase(false), hom(false) {}
};

struct km_query_t {

    double tot;
    double found;
    double shared;
    double uniq_a;
    double uniq_b;

    size_t k;

    bool unphased;

    km_query_t() {

        clear();
    };

    void clear() {

        tot = 0.0;
        found = 0.0;
        shared = 0.0;
        uniq_a = 0.0;
        uniq_b = 0.0;

        k = 0;

        unphased = false;
    }

    void add(const km_query_t& o){

        tot += o.tot;
        found += o.found;
        shared += o.shared;
        uniq_a += o.uniq_a;
        uniq_b += o.uniq_b;

        return;
    }
};

class UnitigData : public CDBG_Data_t<UnitigData> {

    public:

        inline UnitigData() {

            clear();
        }

        inline void clear() {

            km_cov = 0;

            pid_a.clear();
            pid_b.clear();
        }

        inline void clear(const UnitigMap<UnitigData>& um){

            clear();
        }

        void concat(const UnitigMap<UnitigData>& um_dest, const UnitigMap<UnitigData>& um_src){

            const size_t k = um_dest.getGraph()->getK();

            const UnitigData* data_dest = um_dest.getData();
            const UnitigData* data_src = um_src.getData();

            const size_t len_dest = (um_dest.size - k + 1);
            const size_t len_src = (um_src.size - k + 1);

            km_cov = data_dest->getKmerCoverage() + data_src->getKmerCoverage();

            if (um_dest.strand) {

                pid_a = data_dest->getPosA();
                pid_b = data_dest->getPosB();
            }
            else {

                for (const size_t pos_a : data_dest->getPosA()) pid_a.add(len_dest - pos_a - 1);
                for (const size_t pos_b : data_dest->getPosB()) pid_b.add(len_dest - pos_b - 1);
            }

            if (um_src.strand) {

                for (const size_t pos_a : data_src->getPosA()) pid_a.add(len_dest + pos_a);
                for (const size_t pos_b : data_src->getPosB()) pid_b.add(len_dest + pos_b);
            }
            else {

                for (const size_t pos_a : data_src->getPosA()) pid_a.add(len_dest + len_src - pos_a - 1);
                for (const size_t pos_b : data_src->getPosB()) pid_b.add(len_dest + len_src - pos_b - 1);
            }
        }

        void merge(const UnitigMap<UnitigData>& um_dest, const const_UnitigMap<UnitigData>& um_src){

            const size_t k = um_dest.getGraph()->getK();
            const UnitigData* data_src = um_src.getData();

            km_cov = data_src->getKmerCoverage();

            if (um_src.strand) {

                for (const size_t pos_a : data_src->getPosA()){

                    if ((pos_a >= um_src.dist) && (pos_a < (um_src.dist + um_src.len))) pid_a.add(pos_a - um_src.dist + um_dest.dist);
                }

                for (const size_t pos_b : data_src->getPosB()){

                    if ((pos_b >= um_src.dist) && (pos_b < (um_src.dist + um_src.len))) pid_b.add(pos_b - um_src.dist + um_dest.dist);
                }
            }
            else {

                for (const size_t pos_a : data_src->getPosA()){

                    if ((pos_a >= um_src.dist) && (pos_a < (um_src.dist + um_src.len))) pid_a.add(um_src.len - (pos_a - um_src.dist - 1) + um_dest.dist);
                }

                for (const size_t pos_b : data_src->getPosB()){

                    if ((pos_b >= um_src.dist) && (pos_b < (um_src.dist + um_src.len))) pid_b.add(um_src.len - (pos_b - um_src.dist - 1) + um_dest.dist);
                }
            }
        }

        inline string serialize(const const_UnitigMap<UnitigData>& um) const {

            return string(  "KC:Z:" + std::to_string(getKmerCoverage()) + '\t' +
                            "UC:Z:" + std::to_string(static_cast<size_t>(getUnitigCoverage(um))));
        }

        inline void increaseKmerCoverage(const size_t count) {

            if (count > (0xffffffffffffffffULL - km_cov)) km_cov = 0xffffffffffffffffULL;
            else km_cov += count;
        }

        inline void decreaseKmerCoverage(const size_t count) {
            
            if (count > km_cov) km_cov = 0;
            else km_cov -= count;
        }

        inline size_t getKmerCoverage() const {

            return km_cov;
        }

        inline double getUnitigCoverage(const const_UnitigMap<UnitigData>& um) const {

            return (static_cast<double>(km_cov) / static_cast<double>(um.size - um.getGraph()->getK() + 1));
        }

        inline PairID& getPosA() {

            return pid_a;
        }

        inline const PairID& getPosA() const {

            return pid_a;
        }

        inline PairID& getPosB() {

            return pid_b;
        }

        inline const PairID& getPosB() const {

            return pid_b;
        }

    private:

        size_t km_cov;

        PairID pid_a;
        PairID pid_b;
};

void PrintUsage(const TrioBinning_Opt& opt) {

    cout << endl << "TrioBinning " << endl << endl;
    cout << "Usage: TrioBinning [PARAMETERS]" << endl << endl;
    cout << "[COMMAND]:" << endl << endl;
    cout << "   phase               Phase reads or contigs" << endl;
    cout << "   hom                 Detect homozygous regions" << endl << endl;

    cout << "[PARAMETERS]: phase" << endl << endl;
    cout << "   > Mandatory with required argument:" << endl << endl <<
    "   -a, --input-read-parentA    Input read file/s in fasta/fastq(.gz) from parent A" << endl <<
    "   -A, --input-asm-parentA     Input assembly file/s fasta(.gz) from parent A" << endl <<
    "   -b, --input-read-parentB    Input read file/s in fasta/fastq(.gz) from parent B" << endl <<
    "   -B, --input-asm-parentB     Input assembly file/s fasta(.gz) from parent B" << endl <<
    "   -c, --input-read-proband    Input read file/s in fasta/fastq(.gz) from proband." << endl <<
    "   -C, --input-asm-proband     Input assembly file/s in fasta(.gz) from proband." << endl <<
    "   -o, --output-file           Prefix of the output file(s)" << endl <<
    endl << "   > Optional with required argument:" << endl << endl <<
    "   -k, --kmer                  Length of k-mers (default: " << opt.k << ")" << endl <<
    "   -t, --threads               Number of threads (default: " << opt.nb_threads << ")" << endl <<
    "   -m, --min-cov               Minimum coverage (default: " << ((opt.min_cov == -1) ? string("No minimum") : std::to_string(opt.min_cov)) << ")" << endl <<
    "   -M, --max-cov               Maximum coverage (default: " << ((opt.max_cov == -1) ? string("No maximum") : std::to_string(opt.max_cov)) << ")" << endl <<
    "   -U, --max-len-uniq-run      Maximum length of unique k-mer runs (default: " << opt.max_len_uniq_run << ")" << endl <<
    "   -p, --partioning            Unique k-mer similarity threshold for partioning a read (default: " << opt.threshold_partioning << ")" << endl <<
    "   -w, --output-graph-file     Output graph file (default: no output graph)" << endl <<
    "   -r, --input-graph-files     Input prefix for graph and color file (default: no reading, graph is built from input files)" << endl <<
    endl << "   > Optional with no argument:" << endl << endl <<
    "   -S, --no-seq-out            Only output the summary file, not the sequences" << endl <<
    "   -f, --force-partioning      Reads with only homozygous kmers are binned in both haplotypes (default: " << (opt.force_partioning ? "true" : "false") << ")" << endl <<
    "   -v, --verbose               Print information messages during execution" << endl << endl;

    cout << "[PARAMETERS]: hom" << endl << endl;
    cout << "   > Mandatory with required argument:" << endl << endl <<
    "   -a, --input-read-parentA    Input read file/s in fasta/fastq(.gz) from parent A" << endl <<
    "   -A, --input-asm-parentA     Input assembly file/s fasta(.gz) from parent A" << endl <<
    "   -b, --input-read-parentB    Input read file/s in fasta/fastq(.gz) from parent B" << endl <<
    "   -B, --input-asm-parentB     Input assembly file/s fasta(.gz) from parent B" << endl <<
    "   -c, --input-read-proband    Input read file/s in fasta/fastq(.gz) from proband." << endl <<
    "   -C, --input-asm-proband     Input assembly file/s in fasta(.gz) from proband." << endl <<
    "   -o, --output-file           Prefix of the output file(s)" << endl <<
    endl << "   > Optional with required argument:" << endl << endl <<
    "   -k, --kmer                  Length of k-mers (default: " << opt.k << ")" << endl <<
    "   -t, --threads               Number of threads (default: " << opt.nb_threads << ")" << endl <<
    "   -m, --min-cov               Minimum coverage (default: " << ((opt.min_cov == -1) ? string("No minimum") : std::to_string(opt.min_cov)) << ")" << endl <<
    "   -M, --max-cov               Maximum coverage (default: " << ((opt.max_cov == -1) ? string("No maximum") : std::to_string(opt.max_cov)) << ")" << endl <<
    "   -s, --min-len-shared-run    Minimum length of shared k-mer runs (default: " << opt.min_len_shared_run << ")" << endl <<
    "   -w, --output-graph-file     Output graph file (default: no output graph)" << endl <<
    "   -r, --input-graph-files     Input prefix for graph and color file (default: no reading, graph is built from input files)" << endl <<
    endl << "   > Optional with no argument:" << endl << endl <<
    "   -v, --verbose               Print information messages during execution" << endl << endl;
}

bool parse_ProgramOptions(int argc, char **argv, TrioBinning_Opt& opt) {

    int option_index = 0, c;

    const char* opt_string = "a:A:b:B:c:C:o:k:t:m:M:U:p:s:w:r:Sfv";

    static struct option long_options[] = {

        {"input-read-parentA",  required_argument,  0, 'a'},
        {"input-asm-parentA",   required_argument,  0, 'A'},
        {"input-read-parentB",  required_argument,  0, 'b'},
        {"input-asm-parentB",   required_argument,  0, 'B'},
        {"input-read-proband",  required_argument,  0, 'c'},
        {"input-asm-proband",   required_argument,  0, 'C'},
        {"output-file",         required_argument,  0, 'o'},
        {"kmer",                required_argument,  0, 'k'},
        {"threads",             required_argument,  0, 't'},
        {"min-cov",             required_argument,  0, 'm'},
        {"max-cov",             required_argument,  0, 'M'},
        {"max-len-uniq-run",    required_argument,  0, 'U'},
        {"partioning",          required_argument,  0, 'p'},
        {"min-len-shared-run",  required_argument,  0, 's'},
        {"output-graph-file",   required_argument,  0, 'w'},
        {"input-graph-files",   required_argument,  0, 'r'},
        {"no-seq-out",          no_argument,        0, 'S'},
        {"force-partioning",    no_argument,        0, 'f'},
        {"verbose",             no_argument,        0, 'v'},
        {0,                     0,                  0,  0 }
    };

    if (strcmp(argv[1], "phase") == 0) opt.phase = true;
    else if (strcmp(argv[1], "hom") == 0) opt.hom = true;

    if (opt.phase || opt.hom) {

        while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

            switch (c) {

                case 'a':
                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) opt.in_fn_read_parentA.push_back(argv[optind]);
                    break;
                case 'b':
                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) opt.in_fn_read_parentB.push_back(argv[optind]);
                    break;
                case 'c':
                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) opt.in_fn_read_proband.push_back(argv[optind]);
                    break;
                case 'A':
                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) opt.in_fn_asm_parentA.push_back(argv[optind]);
                    break;
                case 'B':
                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) opt.in_fn_asm_parentB.push_back(argv[optind]);
                    break;
                case 'C':
                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind) opt.in_fn_asm_proband.push_back(argv[optind]);
                    break;
                case 'o':
                    opt.out_fn = optarg;
                    break;
                case 'k':
                    opt.k = atoi(optarg);
                    break;
                case 't':
                    opt.nb_threads = atoi(optarg);
                    break;
                case 'm':
                    opt.min_cov = atoi(optarg);
                    break;
                case 'M':
                    opt.max_cov = atoi(optarg);
                    break;
                case 'U':
                    opt.max_len_uniq_run = atoi(optarg);
                    break;
                case 'p':
                    opt.threshold_partioning = atof(optarg);
                    break;
                case 's':
                    opt.min_len_shared_run = atof(optarg);
                    break;
                case 'w':
                    opt.out_graph_fn = optarg;
                    break;
                case 'r':
                    opt.in_graph_fn = optarg;
                    break;
                case 'S':
                    opt.no_seq_out = true;
                    break;
                case 'f':
                    opt.force_partioning = true;
                    break;
                case 'v':
                    opt.verbose = true;
                    break;
                default: return false;
            }
        }

        return true;
    }

    return false;
}

bool check_ProgramOptions(TrioBinning_Opt& opt) {

    bool ret = true;

    const bool input_graph_prefix = (opt.in_graph_fn.length() != 0);

    const size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [](vector<string>& v_files) -> bool {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        bool ret = true;

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "TrioBinning::TrioBinning(): File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const string s_ext = file.substr(file.find_last_of(".") + 1);

                if ((s_ext == "txt")){

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        while (i_file_txt.getline(buffer, 4096)){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "TrioBinning::TrioBinning(): Could not open file " << buffer << " for reading." << endl;
                                ret = false;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "TrioBinning::TrioBinning(): Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
                else files_tmp.push_back(file);
            }
        }

        delete[] buffer;

        v_files = move(files_tmp);

        return ret;
    };

    if (opt.nb_threads <= 0){

        cerr << "TrioBinning::TrioBinning(): Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "TrioBinning::TrioBinning(): Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "TrioBinning::TrioBinning(): Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << "." << endl;
        ret = false;
    }

    if (opt.min_cov < 1){

        cerr << "TrioBinning::TrioBinning(): Minimum coverage (" << opt.min_cov << ") must be at least 1." << endl;
        ret = false;
    }

    if ((opt.max_cov != -1) && (opt.min_cov > opt.max_cov)){

        cerr << "TrioBinning::TrioBinning(): Minimum coverage (" << opt.min_cov << ") cannot be greater than maximum coverage (" << opt.max_cov << ") ." << endl;
        ret = false;
    }

    if ((opt.threshold_partioning <= 0.0) || (opt.threshold_partioning > 1.0)) {

        cerr << "TrioBinning::TrioBinning(): Unique k-mer similarity threshold to partition a read (" << opt.min_cov << ") must be >= 0.0 and < 1.0." << endl;
        ret = false;
    }

    if (opt.max_len_uniq_run < 1) {

        cerr << "TrioBinning::TrioBinning(): Maximum unique k-mer run length (" << opt.max_len_uniq_run << ") must be at least 1." << endl;
        ret = false;
    }

    if (opt.min_len_shared_run < 1) {

        cerr << "TrioBinning::TrioBinning(): Minimum shared k-mer run length (" << opt.min_len_shared_run << ") must be at least 1." << endl;
        ret = false;
    }

    if (!input_graph_prefix && opt.in_fn_read_parentA.empty() && opt.in_fn_asm_parentA.empty()) {

        cerr << "TrioBinning::TrioBinning(): Missing input file/s for parent A." << endl;
        ret = false;
    }
    else if (!input_graph_prefix && opt.in_fn_read_parentB.empty() && opt.in_fn_asm_parentB.empty()) {

        cerr << "TrioBinning::TrioBinning(): Missing input file/s for parent B." << endl;
        ret = false;
    }
    else if (!input_graph_prefix && !opt.in_fn_read_parentA.empty() && !opt.in_fn_asm_parentA.empty()) {

        cerr << "TrioBinning::TrioBinning(): Cannot use reads and assembly for parent A." << endl;
        ret = false;
    }
    else if (!input_graph_prefix && !opt.in_fn_read_parentB.empty() && !opt.in_fn_asm_parentB.empty()) {

        cerr << "TrioBinning::TrioBinning(): Cannot use reads and assembly for parent B." << endl;
        ret = false;
    }
    else if (opt.in_fn_read_proband.empty() && opt.in_fn_asm_proband.empty()) {

        cerr << "TrioBinning::TrioBinning(): Missing input file/s for proband." << endl;
        ret = false;
    }
    else if (!opt.in_fn_read_proband.empty() && !opt.in_fn_asm_proband.empty()) {

        cerr << "TrioBinning::TrioBinning(): Cannot use reads and assembly for proband." << endl;
        ret = false;
    }
    else {

        if (!input_graph_prefix) {

            ret = ret && check_files(opt.in_fn_asm_parentA.empty() ? opt.in_fn_read_parentA : opt.in_fn_asm_parentA);
            ret = ret && check_files(opt.in_fn_asm_parentB.empty() ? opt.in_fn_read_parentB : opt.in_fn_asm_parentB);
        }

        ret = ret && check_files(opt.in_fn_asm_proband.empty() ? opt.in_fn_read_proband : opt.in_fn_asm_proband);
    }

    if (opt.out_fn.length() == 0){

        cerr << "TrioBinning::TrioBinning(): Missing output filename." << endl;
        ret = false;
    }
    else {

	    FILE* fp = fopen(opt.out_fn.c_str(), "w");

	    if (fp == NULL) {

	        cerr << "TrioBinning::TrioBinning(): Could not open output filename for writing: " << opt.out_fn << "." << endl;
	        ret = false;
	    }
	    else {

	        fclose(fp);

	        if (remove(opt.out_fn.c_str()) != 0) cerr << "TrioBinning::TrioBinning(): Could not remove temporary file " << opt.out_fn << "." << endl;
	    }
	}

    if (opt.out_graph_fn.length() != 0) {

        FILE* fp = fopen(opt.out_graph_fn.c_str(), "w");

        if (fp == NULL) {

            cerr << "TrioBinning::TrioBinning(): Could not open output graph filename for writing: " << opt.out_graph_fn << "." << endl;
            ret = false;
        }
        else {

            fclose(fp);

            if (remove(opt.out_graph_fn.c_str()) != 0) cerr << "TrioBinning::TrioBinning(): Could not remove temporary file " << opt.out_graph_fn << "." << endl;
        }
    }

    if (input_graph_prefix) {

        if (!check_file_exists(string(opt.in_graph_fn + "_parentA.gfa")) && !check_file_exists(string(opt.in_graph_fn + "_parentA.gfa.gz"))) {

            cerr << "TrioBinning::TrioBinning(): Could not open input proband graph filename for reading: " << string(opt.in_graph_fn + "._parentA.gfa(.gz)") << "." << endl;
            ret = false;
        }

        if (!check_file_exists(string(opt.in_graph_fn + "_parentB.gfa")) && !check_file_exists(string(opt.in_graph_fn + "_parentB.gfa.gz"))) {

            cerr << "TrioBinning::TrioBinning(): Could not open input proband graph filename for reading: " << string(opt.in_graph_fn + "._parentB.gfa(.gz)") << "." << endl;
            ret = false;
        }
    }

    return ret;
}

inline void toUpperCase(char* s, const size_t len) {

    for (char* s_tmp = s; s_tmp < (s + len); ++s_tmp) *s_tmp &= 0xDF;
}

void addCoverage(CompactedDBG<UnitigData>& cdbg, const CDBG_Build_opt& opt, const vector<string>& v_in_fn) {

    if (cdbg.isInvalid() || (cdbg.length() == 0)) return;

    const size_t k = cdbg.getK();
    const size_t thread_seq_buf_sz = 1048576; // 1Mb
    const size_t thread_name_buf_sz = (thread_seq_buf_sz / (k + 1)) + 1;

    string s;

    size_t len_read;
    size_t pos_read;

    FileParser fp(v_in_fn);

    LockGraph lck_g(opt.nb_threads * 1024);

    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz) {

        const char* str_end = seq_buf + seq_buf_sz;

        while (seq_buf < str_end) { // for each input

            const int len = strlen(seq_buf);

            toUpperCase(seq_buf, len);

            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;

            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

                const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
                const UnitigMap<UnitigData> um = cdbg.findUnitig(seq_buf, p_.second, len);

                if (!um.isEmpty) { // Read maps to a Unitig

                    const uint64_t h = um.getUnitigHead().hash();
                    
                    UnitigData* ud = um.getData();

                    lck_g.lock_unitig(h);

                    ud->increaseKmerCoverage(um.len);

                    lck_g.unlock_unitig(h);

                    it_kmer_h += um.len - 1;
                }
            }

            seq_buf += len + 1;
        }
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_seq_buf = thread_seq_buf_sz - k;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_seq_buf){

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k){

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';

                        pos_read += sz_seq_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    char** buffer_seq = new char*[opt.nb_threads];
    size_t* buffer_seq_sz = new size_t[opt.nb_threads];

    len_read = 0;
    pos_read = 0;

    if (opt.nb_threads == 1){

        buffer_seq[0] = new char[thread_seq_buf_sz];

        while (!reading_function(buffer_seq[0], buffer_seq_sz[0])) worker_function(buffer_seq[0], buffer_seq_sz[0]);
    }
    else {

        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            buffer_seq[t] = new char[thread_seq_buf_sz];

            workers.emplace_back(

                [&, t]{

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) return;

                            stop = reading_function(buffer_seq[t], buffer_seq_sz[t]);
                        }

                        worker_function(buffer_seq[t], buffer_seq_sz[t]);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    for (size_t t = 0; t < opt.nb_threads; ++t) delete[] buffer_seq[t];

    delete[] buffer_seq;
    delete[] buffer_seq_sz;

    fp.close();
}

void trimCoverage(CompactedDBG<UnitigData>& cdbg, const TrioBinning_Opt& opt) {

    if (opt.verbose){

        cout << "TrioBinning::trimCoverage(): Trimming graph using minimum coverage " << opt.min_cov << "." << endl;

        if (opt.max_cov != -1) cout << "TrioBinning::trimCoverage(): Trimming graph using maximum coverage " << opt.max_cov << "." << endl;
    }

    vector<Kmer> v_km_delete;

    for (const auto& um : cdbg) {

        const size_t cov = static_cast<size_t>(um.getData()->getUnitigCoverage(um));

        if ((cov < opt.min_cov) || ((opt.max_cov != -1) && (cov > opt.max_cov))) v_km_delete.push_back(um.getUnitigHead());
    }

    for (const auto& km : v_km_delete) {

        const UnitigMap<UnitigData> um = cdbg.find(km);

        if (!um.isEmpty) {

            const size_t cov = static_cast<size_t>(um.getData()->getUnitigCoverage(um));

            if ((cov < opt.min_cov) || ((opt.max_cov != -1) && (cov > opt.max_cov))) cdbg.remove(um);
        }
    }

    if (opt.verbose) cout << "TrioBinning::trimCoverage(): " << cdbg.nbKmers() << " left in the graph." << endl;
}

void trimKmers(CompactedDBG<UnitigData>& cdbg, const TrioBinning_Opt& opt, const size_t nb_km_min) {

    double min_cov = 2.0;

    size_t nb_km = cdbg.nbKmers();

    vector<Kmer> v_km_delete;

    while (nb_km > nb_km_min) {

        nb_km = 0;

        for (const auto& um : cdbg) {

            if (um.getData()->getUnitigCoverage(um) >= min_cov) nb_km += um.len;
        }

        if (nb_km > nb_km_min) min_cov += 1.0;
    }

    if (nb_km < nb_km_min) min_cov = max(min_cov - 1.0, 2.0);

    if (opt.verbose) cout << "TrioBinning::trimKmers(): Trimming graph using coverage " << min_cov << "." << endl;

    for (const auto& um : cdbg) {

        if (um.getData()->getUnitigCoverage(um) < min_cov) v_km_delete.push_back(um.getUnitigHead());
    }

    for (const auto& km : v_km_delete) {

        const UnitigMap<UnitigData> um = cdbg.find(km);

        if (!um.isEmpty && (um.getData()->getUnitigCoverage(um) < min_cov)) cdbg.remove(um);
    }

    if (opt.verbose) cout << "TrioBinning::trimKmers(): " << cdbg.nbKmers() << " left in the graph." << endl;
}

void colorParents(CompactedDBG<UnitigData>& cdbg, const CDBG_Build_opt& opt, const string& prefix_fn_graph_parentA, const string& prefix_fn_graph_parentB) {

    if (cdbg.isInvalid() || (cdbg.length() == 0)) return;

    const size_t k = cdbg.getK();
    const size_t thread_seq_buf_sz = 1048576; // 1Mb
    const size_t thread_name_buf_sz = (thread_seq_buf_sz / (k + 1)) + 1;

    string s;

    size_t len_read = 0;
    size_t pos_read = 0;

    LockGraph lck_g(opt.nb_threads * 1024);

    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz, const bool isParentA) {

        const char* str_end = seq_buf + seq_buf_sz;

        while (seq_buf < str_end) { // for each input

            const int len = strlen(seq_buf);

            toUpperCase(seq_buf, len);

            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;

            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

                const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
                const UnitigMap<UnitigData> um = cdbg.findUnitig(seq_buf, p_.second, len);

                if (!um.isEmpty) { // Read maps to a Unitig

                    const uint64_t h = um.getUnitigHead().hash();
                    
                    PairID& pid = isParentA ? um.getData()->getPosA() : um.getData()->getPosB();

                    lck_g.lock_unitig(h);

                    for (size_t i = um.dist; i < um.dist + um.len; ++i) pid.add(i);

                    lck_g.unlock_unitig(h);

                    it_kmer_h += um.len - 1;
                }
            }

            seq_buf += len + 1;
        }
    };

    auto reading_function = [&](FileParser& fp, char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_seq_buf = thread_seq_buf_sz - k;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_seq_buf){

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k){

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';

                        pos_read += sz_seq_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    char** buffer_seq = new char*[opt.nb_threads];
    size_t* buffer_seq_sz = new size_t[opt.nb_threads];

    if (opt.nb_threads == 1){

        buffer_seq[0] = new char[thread_seq_buf_sz];

        {
            const vector<string> v_fn_graph_parentA(1, string(prefix_fn_graph_parentA + ".gfa.gz"));

            FileParser fp(v_fn_graph_parentA);

            while (!reading_function(fp, buffer_seq[0], buffer_seq_sz[0])) worker_function(buffer_seq[0], buffer_seq_sz[0], true);
        }

        len_read = 0;
        pos_read = 0;

        s.clear();

        {
            const vector<string> v_fn_graph_parentB(1, string(prefix_fn_graph_parentB + ".gfa.gz"));

            FileParser fp(v_fn_graph_parentB);

            while (!reading_function(fp, buffer_seq[0], buffer_seq_sz[0])) worker_function(buffer_seq[0], buffer_seq_sz[0], false);
        }

    }
    else {

        {
            const vector<string> v_fn_graph_parentA(1, string(prefix_fn_graph_parentA + ".gfa.gz"));

            FileParser fp(v_fn_graph_parentA);

            bool stop = false;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file;

            for (size_t t = 0; t < opt.nb_threads; ++t){

                buffer_seq[t] = new char[thread_seq_buf_sz];

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_file);

                                if (stop) {

                                    delete[] buffer_seq[t];
                                    return;
                                }

                                stop = reading_function(fp, buffer_seq[t], buffer_seq_sz[t]);
                            }

                            worker_function(buffer_seq[t], buffer_seq_sz[t], true);
                        }

                        delete[] buffer_seq[t];
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        len_read = 0;
        pos_read = 0;

        s.clear();

        {
            const vector<string> v_fn_graph_parentB(1, string(prefix_fn_graph_parentB + ".gfa.gz"));

            FileParser fp(v_fn_graph_parentB);

            bool stop = false;

            vector<thread> workers; // need to keep track of threads so we can join them

            mutex mutex_file;

            for (size_t t = 0; t < opt.nb_threads; ++t){

                buffer_seq[t] = new char[thread_seq_buf_sz];

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            {
                                unique_lock<mutex> lock(mutex_file);

                                if (stop) {

                                    delete[] buffer_seq[t];
                                    return;
                                }

                                stop = reading_function(fp, buffer_seq[t], buffer_seq_sz[t]);
                            }

                            worker_function(buffer_seq[t], buffer_seq_sz[t], false);
                        }

                        delete[] buffer_seq[t];
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }

    delete[] buffer_seq;
    delete[] buffer_seq_sz;
}

inline char getQual(const double score, const size_t qv_min = 0) {

    const char phred_base_std = static_cast<char>(33);
    const char phred_scale_std = static_cast<char>(40);

    const double qv_score = score * static_cast<double>(phred_scale_std - qv_min);

    return static_cast<char>(qv_score + phred_base_std + qv_min);
}

inline double getScore(const char c, const size_t qv_min = 0) {

    const char phred_base_std =  static_cast<char>(33);
    const char phred_scale_std =  static_cast<char>(40);

    const double qv_score = static_cast<double>(c - phred_base_std - qv_min);

    return (qv_score / static_cast<double>(phred_scale_std - qv_min));
}

inline size_t getMeanScore(const char* s, const size_t l){

    double sum = 0;

    for (const char* s_end = s + l; s < s_end; ++s) sum += getScore(*s);

    return static_cast<size_t>(sum * 40 / l);
}

void phaseSeq(CompactedDBG<UnitigData>& cdbg, const TrioBinning_Opt& opt_tb) {

    const size_t k = cdbg.getK(); // K-mer size
    const size_t buffer_sz = 1048576; // 1 MB reading buffer per thread

    const char min_qual = getQual(opt_tb.min_qual_score);
    const char min_qual_0 = getQual(0.0); // Lowest possible quality

    const string out_fn_ext = (opt_tb.in_fn_asm_proband.empty() ? "fq" : "fa");

    ofstream outfile_tsv, outfile_h1, outfile_h2, outfile_unphased;
    ostream out_tsv(0), out_h1(0), out_h2(0), out_unphased(0);

    FileParser fp(opt_tb.in_fn_asm_proband.empty() ? opt_tb.in_fn_read_proband : opt_tb.in_fn_asm_proband);

    size_t file_id;

    auto randGenerator = std::bind(std::uniform_int_distribution<>(0,1), std::default_random_engine());

    outfile_tsv.open(string(opt_tb.out_fn + ".tsv").c_str());
    out_tsv.rdbuf(outfile_tsv.rdbuf());

    if (!opt_tb.no_seq_out) {

        outfile_h1.open(string(opt_tb.out_fn + "_h1." + out_fn_ext).c_str());
        outfile_h2.open(string(opt_tb.out_fn + "_h2." + out_fn_ext).c_str());
        outfile_unphased.open(string(opt_tb.out_fn + "_unphased." + out_fn_ext).c_str());

        out_h1.rdbuf(outfile_h1.rdbuf());
        out_h2.rdbuf(outfile_h2.rdbuf());
        out_unphased.rdbuf(outfile_unphased.rdbuf());
    }

    auto queryRead = [&](const string& q_read) {

        const char* s_str = q_read.c_str();
        const size_t s_len = q_read.length();

        KmerHashIterator<RepHash> it_kmer_h(s_str, s_len, k), it_kmer_h_end;

        PairID pid_a, pid_b;

        for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

            const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
            const UnitigMap<UnitigData> um = cdbg.findUnitig(s_str, p_.second, s_len);

            if (!um.isEmpty) { // Read maps to a Unitig

                const PairID& l_pid_a = um.getData()->getPosA();
                const PairID& l_pid_b = um.getData()->getPosB();

                for (size_t i = um.dist; i < um.dist + um.len; ++i) {

                    if (l_pid_a.contains(i)) pid_a.add(p_.second + i - um.dist);
                    if (l_pid_b.contains(i)) pid_b.add(p_.second + i - um.dist);
                }

                it_kmer_h += um.len - 1;
            }
        }

        return pair<PairID, PairID>(pid_a, pid_b);
    };

    auto countRuns = [&](const pair<PairID, PairID>& pid, const size_t min_len_run, const size_t len_km_query){

        km_query_t km_query;

        km_query.tot = len_km_query;

        for (size_t i = 0; i < len_km_query; ++i){

            if (pid.first.contains(i) && !pid.second.contains(i)){

                size_t len_run = 1;

                for (size_t j = i+1; j < len_km_query; ++j){

                    if (pid.first.contains(j) && !pid.second.contains(j)) ++len_run;
                    else break;
                }

                if (len_run >= min_len_run) km_query.uniq_a += static_cast<double>(len_run - min_len_run + 1);

                i += len_run - 1;
            }
        }

        for (size_t i = 0; i < len_km_query; ++i){

            if (!pid.first.contains(i) && pid.second.contains(i)){

                size_t len_run = 1;

                for (size_t j = i+1; j < len_km_query; ++j){

                    if (!pid.first.contains(j) && pid.second.contains(j)) ++len_run;
                    else break;
                }

                if (len_run >= min_len_run) km_query.uniq_b += static_cast<double>(len_run - min_len_run + 1);

                i += len_run - 1;
            }
        }

        for (size_t i = 0; i < len_km_query; ++i){

            if (pid.first.contains(i) && pid.second.contains(i)){

                size_t len_run = 1;

                for (size_t j = i+1; j < len_km_query; ++j){

                    if (pid.first.contains(j) && pid.second.contains(j)) ++len_run;
                    else break;
                }

                if (len_run >= min_len_run) km_query.shared += static_cast<double>(len_run - min_len_run + 1);

                i += len_run - 1;
            }
        }

        return km_query;
    };

    if (opt_tb.nb_threads == 1){

        string in_read, in_name, in_qual;

        while (fp.read(in_read, file_id)) {

            in_name = string(fp.getNameString());
            in_qual = string();

            if (fp.getQualityScoreString() != nullptr) in_qual = string(fp.getQualityScoreString());

            std::transform(in_read.begin(), in_read.end(), in_read.begin(), ::toupper);

            km_query_t km_query, km_query_min_run;

            if (in_read.length() < k) km_query.unphased = true;
            else {

                const pair<PairID, PairID> p_pid = queryRead(in_read);

                size_t min_len_run = 1;

                bool phasing_unclear = false;

                for (size_t len_run = min_len_run; len_run <= opt_tb.max_len_uniq_run; ++len_run){

                    km_query = countRuns(p_pid, len_run, in_read.length() - k + 1);
                    km_query.k = len_run + k - 1;

                    const double uniq_km = km_query.uniq_a + km_query.uniq_b;

                    if (uniq_km != 0.0) {

                        if ((km_query.uniq_a >= opt_tb.threshold_partioning * uniq_km) || (km_query.uniq_b >= opt_tb.threshold_partioning * uniq_km)) break;
                        else {

                            if (!phasing_unclear) km_query_min_run = km_query;

                            phasing_unclear = true;
                        }
                    }
                }

                if ((km_query.k == (opt_tb.max_len_uniq_run + k - 1)) && phasing_unclear) {

                    if (opt_tb.force_partioning) km_query = km_query_min_run;

                    km_query.unphased = true;
                }
            }

            out_tsv << in_name << "\t" << km_query.shared << "\t" << km_query.uniq_a << "\t" << km_query.uniq_b << "\t" << km_query.k << "\t" << endl;

            if (!opt_tb.no_seq_out) {

                if (in_qual.empty()) {

                    if (km_query.unphased || (km_query.uniq_a + km_query.uniq_b == 0.0)) {

                        if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                            if (km_query.uniq_a > km_query.uniq_b) out_h1 << '>' << in_name << '\n' << in_read << endl;
                            else out_h2 << '>' << in_name << '\n' << in_read << endl;
                        }
                        else out_unphased << '>' << in_name << '\n' << in_read << endl;
                    }
                    // More than 80% of unique k-mers are specific to parent A -> H1
                    else if (km_query.uniq_a >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) out_h1 << '>' << in_name << '\n' << in_read << endl;
                    // More than 80% of unique k-mers are specific to parent B -> H2
                    else if (km_query.uniq_b >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) out_h2 << '>' << in_name << '\n' << in_read << endl;
                    // Otherwise -> Unphased
                    else if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                        if (km_query.uniq_a > km_query.uniq_b) out_h1 << '>' << in_name << '\n' << in_read << endl;
                        else out_h2 << '>' << in_name << '\n' << in_read << endl;
                    }
                    else out_unphased << '>' << in_name << '\n' << in_read << endl;
                }
                else {

                    if (km_query.unphased || (km_query.uniq_a + km_query.uniq_b == 0.0)) {

                        if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                            if (km_query.uniq_a > km_query.uniq_b) out_h1 << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                            else out_h2 << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                        }
                        else out_unphased << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                    }
                    else if (km_query.uniq_a >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) { // More than 80% of unique k-mers are specific to parent A -> H1

                        out_h1 << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                    }
                    else if (km_query.uniq_b >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) { // More than 80% of unique k-mers are specific to parent B -> H2

                        out_h2 << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                    }
                    else if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                        if (km_query.uniq_a > km_query.uniq_b) out_h1 << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                        else out_h2 << '@' << in_name << '\n' << in_qual << '\n' << '+' << '\n' << in_read << endl;
                    }
                    else out_unphased << '@' << in_name << '\n' << in_read << '\n' << '+' << '\n' << in_qual << endl;
                }
            }
        }
    }
    else {

        bool stop = false;

        size_t nb_reads_proc = 0;

        vector<thread> workers; // need to keep track of threads so we can join them

        SpinLock mutex_file_in;
        SpinLock mutex_file_out;

        for (size_t t = 0; t < opt_tb.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    vector<string> v_in_read, v_in_name, v_in_qual;
                    vector<km_query_t> v_out;

                    string in_read;

                    size_t in_read_len = 0;

                    while (true) {

                        mutex_file_in.acquire();

                        if (stop){

                            mutex_file_in.release();
                            return;
                        }
                        else {

                            while (in_read_len < buffer_sz) {

                                stop = !fp.read(in_read, file_id);

                                if (!stop){

                                    in_read_len += in_read.length();

                                    const char* qual = fp.getQualityScoreString();

                                    v_in_read.push_back(move(in_read));
                                    v_in_name.push_back(string(fp.getNameString()));
                                    v_in_qual.push_back((qual != nullptr) ? string(qual) : string());
                                }
                                else break;
                            }

                            mutex_file_in.release();

                            if (!v_in_read.empty()){

                                for (size_t i = 0; i < v_in_read.size(); ++i) {

                                    std::transform(v_in_read[i].begin(), v_in_read[i].end(), v_in_read[i].begin(), ::toupper); // Make sure input string is in upper char

                                    km_query_t km_query, km_query_min_run;

                                    if (v_in_read[i].length() < k) km_query.unphased = true;
                                    else {

                                        const pair<PairID, PairID> p_pid = queryRead(v_in_read[i]);

                                        size_t min_len_run = 1;

                                        bool phasing_unclear = false;

                                        for (size_t len_run = min_len_run; len_run <= opt_tb.max_len_uniq_run; ++len_run){

                                            km_query = countRuns(p_pid, len_run, v_in_read[i].length() - k + 1);
                                            km_query.k = len_run + k - 1;

                                            const double uniq_km = km_query.uniq_a + km_query.uniq_b;

                                            if (uniq_km != 0.0) {

                                                if ((km_query.uniq_a >= opt_tb.threshold_partioning * uniq_km) || (km_query.uniq_b >= opt_tb.threshold_partioning * uniq_km)) break;
                                                else {

                                                    if (!phasing_unclear) km_query_min_run = km_query;

                                                    phasing_unclear = true;
                                                }
                                            }
                                        }

                                        if ((km_query.k == (opt_tb.max_len_uniq_run + k - 1)) && phasing_unclear) {

                                            if (opt_tb.force_partioning) km_query = km_query_min_run;

                                            km_query.unphased = true;
                                        }
                                    }

                                    v_out.push_back(km_query);
                                }

                                mutex_file_out.acquire();

                                for (size_t i = 0; i < v_in_read.size(); ++i){

                                    const km_query_t& km_query = v_out[i];

                                    out_tsv << v_in_name[i] << "\t" << km_query.shared << "\t" << km_query.uniq_a << "\t" << km_query.uniq_b << "\t" << km_query.k << "\t" << endl;

                                    if (!opt_tb.no_seq_out) {

                                        if (v_in_qual[i].empty()) {

                                            if (km_query.unphased || (km_query.uniq_a + km_query.uniq_b == 0.0)) {

                                                if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                                                    if (km_query.uniq_a > km_query.uniq_b) out_h1 << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                                    else out_h2 << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                                }
                                                else out_unphased << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                            }
                                            // More than 80% of unique k-mers are specific to parent A -> H1
                                            else if (km_query.uniq_a >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) out_h1 << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                            // More than 80% of unique k-mers are specific to parent B -> H2
                                            else if (km_query.uniq_b >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) out_h2 << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                            // Otherwise -> Unphased
                                            else if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                                                if (km_query.uniq_a > km_query.uniq_b) out_h1 << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                                else out_h2 << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                            }
                                            else out_unphased << '>' << v_in_name[i] << '\n' << v_in_read[i] << endl;
                                        }
                                        else {

                                            if (km_query.unphased || (km_query.uniq_a + km_query.uniq_b == 0.0)) {

                                                if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                                                    if (km_query.uniq_a > km_query.uniq_b) out_h1 << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                                    else out_h2 << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                                }
                                                else out_unphased << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                            }
                                            else if (km_query.uniq_a >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) { // More than 80% of unique k-mers are specific to parent A -> H1

                                                out_h1 << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                            }
                                            else if (km_query.uniq_b >= opt_tb.threshold_partioning * (km_query.uniq_a + km_query.uniq_b)) { // More than 80% of unique k-mers are specific to parent B -> H2

                                                out_h2 << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                            }
                                            else if (opt_tb.force_partioning && (km_query.uniq_a != km_query.uniq_b)) {

                                                if (km_query.uniq_a > km_query.uniq_b) out_h1 << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                                else out_h2 << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                            }
                                            else out_unphased << '@' << v_in_name[i] << '\n' << v_in_read[i] << '\n' << '+' << '\n' << v_in_qual[i] << endl;
                                        }
                                    }
                                }

                                mutex_file_out.release();
                            }

                            in_read_len = 0;

                            v_in_read.clear();
                            v_in_name.clear();
                            v_in_qual.clear();
                            v_out.clear();

                            in_read.clear();
                        }
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    outfile_tsv.close();

    if (!opt_tb.no_seq_out) {

        outfile_h1.close();
        outfile_h2.close();
        outfile_unphased.close();
    }

    fp.close();
}

void detectHomSeq(CompactedDBG<UnitigData>& cdbg, const TrioBinning_Opt& opt_tb) {

    const size_t k = cdbg.getK(); // K-mer size
    const size_t buffer_sz = 1048576; // 1 MB reading buffer per thread

    const char min_qual = getQual(opt_tb.min_qual_score);
    const char min_qual_0 = getQual(0.0); // Lowest possible quality

    const string out_fn_ext = (opt_tb.in_fn_asm_proband.empty() ? "fq" : "fa");

    ofstream outfile_bed;
    ostream out_bed(0);

    FileParser fp(opt_tb.in_fn_asm_proband.empty() ? opt_tb.in_fn_read_proband : opt_tb.in_fn_asm_proband);

    size_t file_id;

    auto randGenerator = std::bind(std::uniform_int_distribution<>(0,1),std::default_random_engine());

    outfile_bed.open(string(opt_tb.out_fn + ".bed").c_str());
    out_bed.rdbuf(outfile_bed.rdbuf());

    auto querySeq = [&](const string& q_read) {

        const char* s_str = q_read.c_str();
        const size_t s_len = q_read.length();

        KmerHashIterator<RepHash> it_kmer_h(s_str, s_len, k), it_kmer_h_end;

        PairID pid_a, pid_b;

        for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

            const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
            const UnitigMap<UnitigData> um = cdbg.findUnitig(s_str, p_.second, s_len);

            if (!um.isEmpty) { // Read maps to a Unitig

                const PairID& l_pid_a = um.getData()->getPosA();
                const PairID& l_pid_b = um.getData()->getPosB();

                for (size_t i = um.dist; i < um.dist + um.len; ++i) {

                    if (l_pid_a.contains(i)) pid_a.add(p_.second + i - um.dist);
                    if (l_pid_b.contains(i)) pid_b.add(p_.second + i - um.dist);
                }

                it_kmer_h += um.len - 1;
            }
        }

        return pair<PairID, PairID>(pid_a, pid_b);
    };

    auto countRuns = [&](const pair<PairID, PairID>& pid, const size_t min_len_run, const size_t len_km_query){

        vector<pair<size_t, size_t>> v_out;

        for (size_t i = 0; i < len_km_query; ++i) {

            bool parentA = pid.first.contains(i);
            bool parentB = pid.second.contains(i);

            if (parentA && parentB){

                size_t len_run = 1;

                for (size_t j = i+1; j < len_km_query; ++j){

                    parentA = pid.first.contains(j);
                    parentB = pid.second.contains(j);

                    if ((parentA && parentB) /*|| (!parentA && !parentB)*/) ++len_run;
                    else break;
                }

                if (len_run >= min_len_run) v_out.push_back({i, i+len_run});

                i += len_run - 1;
            }
        }

        return v_out;
    };

    if (opt_tb.nb_threads == 1){

        string in_read, in_name;

        while (fp.read(in_read, file_id)) {

            if (in_read.length() >= k) {

                in_name = string(fp.getNameString());

                std::transform(in_read.begin(), in_read.end(), in_read.begin(), ::toupper);

                const pair<PairID, PairID> p_pid = querySeq(in_read);
                const vector<pair<size_t, size_t>> v_pos = countRuns(p_pid, opt_tb.min_len_shared_run, in_read.length() - k + 1);

                for (const auto& p : v_pos) out_bed << in_name << "\t" << p.first << "\t" << p.second << endl;
            }
        }
    }
    else {

        bool stop = false;

        size_t nb_reads_proc = 0;

        vector<thread> workers; // need to keep track of threads so we can join them

        SpinLock mutex_file_in;
        SpinLock mutex_file_out;

        for (size_t t = 0; t < opt_tb.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    vector<string> v_in_read, v_in_name;
                    vector<pair<string, pair<size_t, size_t>>> v_out;

                    string in_read;

                    size_t in_read_len = 0;

                    while (true) {

                        mutex_file_in.acquire();

                        if (stop){

                            mutex_file_in.release();
                            return;
                        }
                        else {

                            while (in_read_len < buffer_sz) {

                                stop = !fp.read(in_read, file_id);

                                if (!stop){

                                    if (in_read.length() >= k) {

                                        in_read_len += in_read.length();

                                        v_in_read.push_back(move(in_read));
                                        v_in_name.push_back(string(fp.getNameString()));
                                    }
                                }
                                else break;
                            }

                            mutex_file_in.release();

                            if (!v_in_read.empty()){

                                for (size_t i = 0; i < v_in_read.size(); ++i) {

                                    std::transform(v_in_read[i].begin(), v_in_read[i].end(), v_in_read[i].begin(), ::toupper); // Make sure input string is in upper char

                                    const pair<PairID, PairID> p_pid = querySeq(v_in_read[i]);
                                    const vector<pair<size_t, size_t>> v_pos = countRuns(p_pid, opt_tb.min_len_shared_run, v_in_read[i].length() - k + 1);

                                    for (const auto& p : v_pos) v_out.push_back({v_in_name[i], p});
                                }

                                if (!v_out.empty()) {

                                    mutex_file_out.acquire();

                                    for (const auto& p : v_out) out_bed << p.first << "\t" << p.second.first << "\t" << p.second.second << endl;

                                    mutex_file_out.release();
                                }
                            }

                            in_read_len = 0;

                            v_in_read.clear();
                            v_in_name.clear();
                            v_out.clear();

                            in_read.clear();
                        }
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    outfile_bed.close();
    fp.close();
}

CDBG_Build_opt getCDBG_Build_opt(const TrioBinning_Opt& opt_tb, const vector<string>& v_in_fn, const bool in_reads = true) {

    CDBG_Build_opt opt_out;

    opt_out.verbose = opt_tb.verbose;
    opt_out.nb_threads = opt_tb.nb_threads;
    opt_out.k = opt_tb.k;
    opt_out.build = true;

    if ((opt_tb.min_cov == 1) || !in_reads) opt_out.filename_ref_in = v_in_fn;
    else opt_out.filename_seq_in = v_in_fn;

    return opt_out;
}

int main(int argc, char *argv[]) {

    const TrioBinning_Opt opt_tb_default;

    if (argc < 2) PrintUsage(opt_tb_default);
    else {

        TrioBinning_Opt opt_tb;

		if (!parse_ProgramOptions(argc, argv, opt_tb)) {

            PrintUsage(opt_tb_default);

            return 1;
        }
        else if (!check_ProgramOptions(opt_tb)) return 1;
		else {

            size_t nb_km_proband = 0;

            // No graph loading or proband graph is missing
            if (opt_tb.in_graph_fn.empty() || (!check_file_exists(string(opt_tb.in_graph_fn + "_proband.gfa")) && !check_file_exists(string(opt_tb.in_graph_fn + "_proband.gfa.gz")))) {

                const bool in_asm = !opt_tb.in_fn_asm_proband.empty();

                CDBG_Build_opt opt_cdbg;

                if (!in_asm) opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_read_proband, true);
                else opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_asm_proband, false);

                CompactedDBG<UnitigData> cdbg(opt_tb.k);

                cdbg.build(opt_cdbg);
                cdbg.write(string((opt_tb.in_graph_fn.empty() ? opt_tb.out_fn : opt_tb.in_graph_fn) + "_proband"), opt_cdbg.nb_threads, true, false, false, true, true, opt_cdbg.verbose);

                nb_km_proband = cdbg.nbKmers();
            }

            // Build graph of parent A
            if (opt_tb.in_graph_fn.empty()) {

                const bool in_asm = !opt_tb.in_fn_asm_parentA.empty();

                CDBG_Build_opt opt_cdbg;

                if (!in_asm) opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_read_parentA, true);
                else opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_asm_parentA, false);

                CompactedDBG<UnitigData> cdbg(opt_tb.k);

                cdbg.build(opt_cdbg);
                cdbg.write(string(opt_tb.out_fn + "_parentA"), opt_cdbg.nb_threads, true, false, false, true, true, opt_cdbg.verbose);
            }

            // Build graph of parent B
            if (opt_tb.in_graph_fn.empty()) {

                const bool in_asm = !opt_tb.in_fn_asm_parentB.empty();

                CDBG_Build_opt opt_cdbg;

                if (!in_asm) opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_read_parentB, true);
                else opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_asm_parentB, false);

                CompactedDBG<UnitigData> cdbg(opt_tb.k);

                cdbg.build(opt_cdbg);
                cdbg.write(string(opt_tb.out_fn + "_parentB"), opt_cdbg.nb_threads, true, false, false, true, true, opt_cdbg.verbose);
            }

            // Bin reads from proband
            {
                const bool in_asm = !opt_tb.in_fn_asm_proband.empty();

                CDBG_Build_opt opt_cdbg;

                if (!in_asm) opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_read_proband, true);
                else opt_cdbg = getCDBG_Build_opt(opt_tb, opt_tb.in_fn_asm_proband, false);

                CompactedDBG<UnitigData> cdbg(opt_tb.k);

                if (opt_tb.in_graph_fn.empty()) {

                    cdbg.read(string(opt_tb.out_fn + "_proband.gfa.gz"), opt_cdbg.nb_threads, opt_cdbg.verbose);

                    colorParents(cdbg, opt_cdbg, string(opt_tb.out_fn + "_parentA"), string(opt_tb.out_fn + "_parentB"));
                }
                else {

                    cdbg.read(string(opt_tb.in_graph_fn + "_proband.gfa.gz"), opt_cdbg.nb_threads, opt_cdbg.verbose);

                    colorParents(cdbg, opt_cdbg, string(opt_tb.in_graph_fn + "_parentA"), string(opt_tb.in_graph_fn + "_parentB"));
                }

                if (opt_tb.phase) phaseSeq(cdbg, opt_tb);
                else if (opt_tb.hom) detectHomSeq(cdbg, opt_tb);
            }
        }
	}
}