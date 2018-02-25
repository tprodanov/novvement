#include<iostream>
#include<thread>
#include<fstream>
#include<vector>
#include<string>
#include<algorithm>
#include<thread>
#include<cassert>


template<typename T, typename U=T>
using pair_vec = std::vector<std::pair<T, U>>;


void string_toupper(std::string& str) {
    std::transform(str.begin(), str.end(), str.begin(), toupper);
}


pair_vec<std::string> read_fasta(std::istream& fin) {
    pair_vec<std::string> result;
    std::string name;
    std::string seq;

    while (fin) {
        std::string tmp;
        getline(fin, tmp);
        if (tmp[0] == '>') {
            if (!name.empty()) {
                string_toupper(seq);
                result.push_back(std::make_pair(name, seq));
            }
            name = tmp.substr(1);
            seq = "";
        } else {
            seq += tmp;
        }
    }
    
    string_toupper(seq);
    result.push_back(std::make_pair(name, seq));
    return result;
}


template<typename T>
struct Matrix {

    Matrix(size_t height, size_t width)
        : data_(new T[width * height])
        , height_(height)
        , width_(width)
    {
        //EMPTY
    }

    ~Matrix() {
        delete[] data_;
    }

    T& at(size_t i, size_t j) {
        assert(i < height_ && j < width_);
        return data_[i * width_ + j];
    }

    T const& at(size_t i, size_t j) const {
        assert(i < height_ && j < width_);
        return data_[i * width_ + j];
    }

    void print() const {
        for (size_t i = 0; i < height_; ++i) {
            printf("%-3d ", i);
            for (size_t j = 0; j < width_; ++j) {
                printf("%3.0f ", at(i, j));
            }
            std::cout << std::endl;
        }
    }

private:

    T* data_;
    size_t width_;
    size_t height_;

};


using alignment_res = std::tuple<float, size_t, size_t>;


alignment_res fitting_align(std::string const& long_seq, std::string const& short_seq,
                            float mismatch_penalty, float indel_penalty) {
    size_t n = long_seq.size();
    size_t m = short_seq.size();
    Matrix<float> score(n + 1, m + 1);
    Matrix<size_t> indels(n + 1, m + 1);
    Matrix<size_t> mismatches(n + 1, m + 1);

    for (size_t i = 0; i <= n; ++i) {
        score.at(i, 0) = 0;
        indels.at(i, 0) = 0;
        mismatches.at(i, 0) = 0;
    }
    for (size_t j = 1; j <= m; ++j) {
        score.at(0, j) = j * indel_penalty;
        indels.at(0, j) = j;
        mismatches.at(0, j) = 0;
    }

    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = 1; j <= m; ++j) {
            if (long_seq[i - 1] == short_seq[j - 1]) {
                score.at(i, j) = score.at(i - 1, j - 1) + 1;
                mismatches.at(i, j) = mismatches.at(i - 1, j - 1);
            } else {
                score.at(i, j) = score.at(i - 1, j - 1) + mismatch_penalty;
                mismatches.at(i, j) = mismatches.at(i - 1, j - 1) + 1;
            }
            indels.at(i, j) = indels.at(i - 1, j - 1);

            if (score.at(i, j) < score.at(i - 1, j) + indel_penalty) {
                score.at(i, j) = score.at(i - 1, j) + indel_penalty;
                mismatches.at(i, j) = mismatches.at(i - 1, j);
                indels.at(i, j) = indels.at(i - 1, j) + 1; 
            }
            if (score.at(i, j) < score.at(i, j - 1) + indel_penalty) {
                score.at(i, j) = score.at(i, j - 1) + indel_penalty;
                mismatches.at(i, j) = mismatches.at(i, j - 1);
                indels.at(i, j) = indels.at(i, j - 1) + 1;
            }
        }
    }

    float best_score = -10000;
    size_t best_i = 0;
    for (size_t i = 0; i <= n; ++i) {
        if (score.at(i, m) > best_score) {
            best_score = score.at(i, m);
            best_i = i;
        }
    }
    return std::make_tuple(best_score, mismatches.at(best_i, m), indels.at(best_i, m));
}


std::pair<std::string, alignment_res>
            align_to_genes(std::string const& short_seq, pair_vec<std::string> const& genes,
                           float mismatch_penalty, float indel_penalty) {
    alignment_res best;
    std::string best_gene;
    for (auto const& entry : genes) {
        std::string const& gene_name = entry.first;
        std::string const& gene_seq = entry.second;
        
        auto current = fitting_align(gene_seq, short_seq, mismatch_penalty, indel_penalty);
        if (std::get<0>(current) > std::get<0>(best)) {
            best = current;
            best_gene = gene_name;
        }
    }

    return std::make_pair(best_gene, best);    
}


void align_all(pair_vec<std::string> const& novel,
               pair_vec<std::string> const& genes,
               std::ostream& fout,
               float mismatch_penalty,
               float indel_penalty) {
    for (auto const& entry : novel) {
        auto res = align_to_genes(entry.second, genes, mismatch_penalty, indel_penalty);
        std::string const& gene = res.first;
        float score;
        size_t mismatches;
        size_t indels;
        std::tie(score, mismatches, indels) = res.second;

        fout << entry.first << '\t' << gene << '\t' << score << '\t' << mismatches << '\t' << indels << '\n';
    }
}


int main(int argc, char *argv[]) {
    // Arguments:
    // stdin -- cropped sequences fasta, stdout -- output
    // 1 -- germline genes
    // 2 -- mismatch
    // 3 -- indel
    
    pair_vec<std::string> novel = read_fasta(std::cin);
    std::ifstream genes_f(argv[1]);
    pair_vec<std::string> genes = read_fasta(genes_f);

    float mismatch = atof(argv[2]);
    float indel = atof(argv[3]);    

    align_all(novel, genes, std::cout, mismatch, indel);
    return 0;
}

