#include <iostream>
#include <chrono>
#include <cstdlib>

#include "br_index.hpp"
#include "br_index_nplcp.hpp"
#include "utils.hpp"
#include "nucleotide.h"

using namespace bri;
using namespace std;

string check = string();
long allowed = 0;
bool nplcp = false;

struct ReadRecord {
    string id;
    string read;

    string qual;

    /**
     * Deep copies the strings
     */
    ReadRecord(string id, string read, string qual)
        : id(id), read(read), qual(qual) {
    }
};

string getFileExt(const string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != string::npos) {
        return (s.substr(i + 1, s.length() - i));
    }

    return ("");
}

vector<ReadRecord> getReads(const string& file) {
    vector<ReadRecord> reads;
    reads.reserve(200000);

    const auto& extension = getFileExt(file);

    bool fasta =
        (extension == "FASTA") || (extension == "fasta") || (extension == "fa");
    bool fastq = (extension == "fq") || (extension == "fastq");

    ifstream ifile(file.c_str());
    if (!ifile) {
        throw runtime_error("Cannot open file " + file);
    }
    if (!fasta && !fastq) {
        // this is a not readable

        throw runtime_error("extension " + extension +
                            " is not a valid extension for the readsfile");
    } else if (fasta) {
        // fasta file
        string read = "";
        string id = "";
        string qual = ""; // empty quality string for fasta
        string line;

        while (getline(ifile, line)) {
            if (line.empty()) {
                continue; // Skip empty lines
            }

            if (line[0] == '>' || line[0] == '@') {
                // This is an ID line
                if (!id.empty()) {
                    // If we already have data, process it and clear
                    reads.emplace_back(id, read, qual);
                    reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
                    id.clear();
                    read.clear();
                }
                id = line.substr(1); // Extract ID (skip '>')
            } else {
                // This is a sequence line
                read += line;
            }
        }

        // Process the last entry if it exists
        if (!id.empty()) {
            reads.emplace_back(id, read, qual);
            reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
        }
    } else {
        // fastQ
        string read = "";
        string id = "";
        string qual = "";
        string plusLine = ""; // Skip the '+' line
        string line;

        while (getline(ifile, id) && getline(ifile, read) &&
               getline(ifile, plusLine) && // Skip the '+' line
               getline(ifile, qual)) {
            if (!id.empty() && id[0] != '@') {
                throw runtime_error("File " + file +
                                    "doesn't appear to be in FastQ format");
            }

            if (id.back() == '\n') {
                id.pop_back();
            }
            if (!read.empty() && read.back() == '\n') {
                read.pop_back();
            }

            assert(id.size() > 1);
            id = (id.substr(1));
            reads.emplace_back(id, read, qual);
            reverse(qual.begin(), qual.end());
            reads.emplace_back(id, Nucleotide::getRevCompl(read), qual);
            id.clear(), read.clear(), qual.clear();
        }
    }

    return reads;
}

void help()
{
	cout << "bri-locate: locate all occurrences of the input patterns" << endl;
    cout << "             allowing some mismatched characters."        << endl << endl;

	cout << "Usage: bri-locate [options] <index> <patterns>" << endl;
    cout << "   -nplcp       use the version without PLCP." << endl;
    cout << "   -m <number>  max number of mismatched characters allowed (0 by default)" << endl;
	cout << "   -c <text>    check correctness of each pattern occurrence on this text file (must be the same indexed)" << endl;
	cout << "   <index>      index file (with extension .bri)" << endl;
	cout << "   <patterns>   file in pizza&chili format containing the patterns." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if (s.compare("-c") == 0) 
    {

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -c option." << endl;
			help();
		}

		check = string(argv[ptr]);
		ptr++;
    
    }
    else if (s.compare("-m") == 0)
    {

        if(ptr>=argc-1){
            cout << "Error: missing parameter after -m option." << endl;
            help();
        }

        char* e;
        allowed = strtol(argv[ptr],&e,10);

        if(*e != '\0' || allowed < 0){
            cout << "Error: invalid negative value after -m option." << endl;
            help();
        }

        ptr++;

	}
    else if (s.compare("-nplcp") == 0)
    {

        nplcp = true;

    }
    else
    {

		cout << "Error: unknown option " << s << endl;
		help();

	}

}

template<class T>
void locate_all(ifstream& in, string patterns)
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    string text;
    bool c = false;

    if (check.compare(string()) != 0)
    {
        c = true;

        ifstream ifs1(check);
        stringstream ss;
        ss << ifs1.rdbuf();
        text = ss.str();
    }

    auto t1 = high_resolution_clock::now();

    T idx;

    idx.load(in);

    auto t2 = high_resolution_clock::now();

    cout << "searching patterns with mismatches at most " << allowed << " ... " << endl;

    cout << "Reading in reads from " << patterns << endl;
    vector<ReadRecord> reads;
    try {
        reads = getReads(patterns);
    } catch (const exception& e) {
        string er = e.what();
        er += " Did you provide a valid reads file?";
        throw runtime_error(er);
    }
    // ifstream ifs(patterns);

    //read header of the pizza&chilli input file
	//header example:
	//# number=7 length=10 file=genome.fasta forbidden=\n\t
    // string header;
    // getline(ifs,header);

    ulint n = reads.size();
    // ulint m = get_patterns_length(header);

    ulint last_perc = 0;

    ulint occ_tot = 0;

    auto t3 = high_resolution_clock::now();
    auto t4 = high_resolution_clock::now();
    auto t5 = high_resolution_clock::now();
    ulint count_time = 0;
    ulint locate_time = 0;
    ulint tot_time = 0;

    // extract patterns from file and search them in the index
    for (unsigned int i = 0; i < n; i += 2) {
        ulint perc = 100 * i / n;
        if (perc > last_perc)
        {
            cout << perc << "% done ..." << endl;
            last_perc = perc;
        }

        string p = reads[i].read;

        t3 = high_resolution_clock::now();
        auto samples = idx.search_with_mismatch(p,allowed);
        t4 = high_resolution_clock::now();
        auto occs = idx.locate_samples(samples);
        t5 = high_resolution_clock::now();

        count_time += duration_cast<microseconds>(t4-t3).count();
        locate_time += duration_cast<microseconds>(t5-t4).count();
        occ_tot += occs.size();
        tot_time += duration_cast<microseconds>(t4-t3).count() + duration_cast<microseconds>(t5-t4).count();

        // Now also match the reverse complement
        p = reads[i+1].read;

        t3 = high_resolution_clock::now();
        samples = idx.search_with_mismatch(p,allowed);
        t4 = high_resolution_clock::now();
        occs = idx.locate_samples(samples);
        t5 = high_resolution_clock::now();

        count_time += duration_cast<microseconds>(t4-t3).count();
        locate_time += duration_cast<microseconds>(t5-t4).count();
        occ_tot += occs.size();
        tot_time += duration_cast<microseconds>(t4-t3).count() + duration_cast<microseconds>(t5-t4).count();

        if (c) // check occurrences
        {
            cout << "number of occs with at most " << allowed << " mismatch   : " << occs.size() << endl;
            for (auto o : occs)
            {
                int mismatches = 0;
                for (size_t i = 0; i < p.size(); ++i)
                {
                    if (text[o+i] != p[i]) mismatches++;
                }
                if (mismatches > allowed) 
                {
                    cout << "Error: wrong occurrence:  " << o << endl;
                    cout << "       original pattern:  " << p << endl;
                    cout << "       wrong    pattern:  " << text.substr(o,p.size()) << endl;
                }
            }
        }

        
    }

    double occ_avg = (double)occ_tot / n*2;
    
    cout << endl << occ_avg << " average occurrences per pattern" << endl;

    // ifs.close();

    ulint load = duration_cast<milliseconds>(t2-t1).count();
    cout << "Load time  : " << load << " milliseconds" << endl;

    cout << "Number of patterns             n = " << n/2 << endl;
	// cout << "Pattern length                 m = " << m << endl;
	cout << "Total number of occurrences  occ = " << occ_tot << endl << endl;

    cout << "LF-mapping time: " << count_time << " microseconds" << endl;
    cout << "Phi        time: " << locate_time << " microseconds" << endl;
    cout << "Total time     : " << tot_time << " microseconds" << endl;
	cout << "Search time    : " << (double)tot_time/n*2 << " microseconds/pattern (total: " << n/2 << " patterns)" << endl;
	cout << "Search time    : " << (double)tot_time/occ_tot << " microseconds/occurrence (total: " << occ_tot << " occurrences)" << endl;
}



int main(int argc, char** argv)
{
    if (argc < 3) help();

    int ptr = 1;

    while (ptr < argc - 2) parse_args(argv, argc, ptr);

    string idx_file(argv[ptr]);
    string patt_file(argv[ptr+1]);

    ifstream in(idx_file);

    cout << "Loading br-index" << endl;

    if (nplcp)
        locate_all<br_index_nplcp<> >(in, patt_file);
    else 
        locate_all<br_index<> >(in, patt_file);

    in.close();

}