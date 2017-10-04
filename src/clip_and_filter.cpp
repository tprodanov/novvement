#include <fstream>
#include <iostream>
#include <vector>
#include <string>

size_t hamming_distance(std::string const& seq1, std::string const& seq2) {
    size_t res = std::max(seq1.length(), seq2.length()) - std::min(seq1.length(), seq2.length());
        
    for (size_t i = 0; i < std::min(seq1.length(), seq2.length()); ++i) {
        if (seq1[i] != seq2[i]) {
            res += 1;
        }
    }
    return res;
}

bool compare_with_target(std::string const& seq, std::vector<std::string> const& target, size_t min_dist) {
    for (std::string const& seq2 : target) {
        if (hamming_distance(seq, seq2) < min_dist) {
            return false;
        }
    }
    return true;
}


std::vector<std::string> load_segments(std::istream& inp, size_t start, size_t length) {
    std::vector<std::string> res;
    std::string seq;
    std::string current;

    while (inp) {
        getline(inp, current);
        if (current[0] == '>') {
            if (!seq.empty() && seq.length() > start) {
                res.push_back(seq.substr(start, length));
            }
            seq = "";
        } else {
            seq += current;
        }
    }

    if (!seq.empty() && seq.length() > start) {
        res.push_back(seq.substr(start, length));
    }
    return res;
}


void filter_sequences(std::istream& sequences,
                     std::vector<std::string> const& segments,
                     size_t start, size_t length, size_t min_dist,
                     std::ostream& outp) {
    outp << "sequence\tclipped_sequence\n";
    std::string seq;

    while (sequences) {
        getline(sequences, seq);
        if (seq.length() <= start) {
            continue;
        }

        std::string clipped_seq = seq.substr(start, length);
    
        if (compare_with_target(clipped_seq, segments, min_dist)) {
            outp << seq << '\t' << clipped_seq << '\n';
        }
    }
}


int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "argv should contain exactly 4 arguments: segments, min_dist, left, right" << std::endl;
        exit(1);
    }

    std::ifstream segments_f(argv[1]);
    
    size_t min_dist;
    sscanf(argv[2], "%zu", &min_dist);

    size_t left;
    sscanf(argv[3], "%zu", &left);

    size_t right;
    sscanf(argv[4], "%zu", &right);

    --left;
    size_t length = right - left;

    std::vector<std::string> segments = load_segments(segments_f, left, length);
    filter_sequences(std::cin, segments, left, length, min_dist, std::cout);
}
