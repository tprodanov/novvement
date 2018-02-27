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


using consensus_t = pair_vec<size_t, char>;


struct alignment_args {
    float mismatch;
    float open_gap;
    float extend_gap;
};


struct alignment_res {

    alignment_res()
        : mismatches(0)
        , indels(0)
    {
        // EMPTY
    }

    size_t mismatches;
    size_t indels;
    consensus_t consensus;
};


template<typename T>
T const& max(T const& a, T const& b, T const& c) {
    return std::max(a, std::max(b, c));
}


enum Plane {
    main_plane = 0,
    deletions,
    inserts
};


void construct_consensus(std::string const& long_seq, std::string const& short_seq,
                         Matrix<float> const& score,
                         Matrix<float> const& deletions,
                         Matrix<float> const& inserts,
                         size_t i,
                         alignment_args const& penalties,
                         alignment_res& result) {
    size_t j = short_seq.size();

    Plane plane = Plane::main_plane;
    if (deletions.at(i, j) > std::max(score.at(i, j), inserts.at(i, j))) {
        plane = Plane::deletions;
    } else if (inserts.at(i, j) > std::max(score.at(i, j), deletions.at(i, j))) {
        plane = Plane::inserts;
    }

    while (j > 0) {
        if (plane == Plane::main_plane) {
            assert(i);
            float match = long_seq[i - 1] == short_seq[j - 1] ? 1 : penalties.mismatch;
            if (score.at(i, j) == deletions.at(i - 1, j - 1) + match) {
                plane = Plane::deletions;
            } else if(score.at(i, j) == inserts.at(i - 1, j - 1) + match) {
                plane = Plane::inserts;
            } else {
                assert(score.at(i, j) == score.at(i - 1, j - 1) + match);
            }
            if (long_seq[i - 1] != short_seq[j - 1]) {
                result.mismatches += 1;
                result.consensus.push_back(std::make_pair(i, short_seq[j - 1]));
            }
            --i;
            --j;
        }

        else if (plane == Plane::deletions) {
            if (deletions.at(i, j) == score.at(i - 1, j) + penalties.open_gap + penalties.extend_gap) {
                plane = Plane::main_plane;
            } else {
                assert(deletions.at(i, j) == deletions.at(i - 1, j) + penalties.extend_gap);
            }
            result.indels += 1;
            result.consensus.push_back(std::make_pair(i, '-'));
            --i;
        }

        else {
            assert(plane == Plane::inserts);
            if (inserts.at(i, j) == score.at(i, j - 1) + penalties.open_gap + penalties.extend_gap) {
                plane = Plane::main_plane;
            } else {
                assert(inserts.at(i, j) == inserts.at(i, j - 1) + penalties.extend_gap);
            }
            result.indels +=1;
            result.consensus.push_back(std::make_pair(i, tolower(short_seq[j - 1])));
            --j;
        }
    }
}


float fitting_alignment(std::string const& long_seq, std::string const& short_seq,
                        alignment_args const& penalties,
                        alignment_res* result = nullptr) {
    size_t n = long_seq.size();
    size_t m = short_seq.size();
    Matrix<float> score(n + 1, m + 1);
    Matrix<float> deletions(n + 1, m + 1);
    Matrix<float> inserts(n + 1, m + 1);

    for (size_t i = 0; i <= n; ++i) {
        score.at(i, 0) = 0;
        deletions.at(i, 0) = penalties.open_gap + penalties.extend_gap;
        inserts.at(i, 0) = penalties.open_gap + penalties.extend_gap;
    }

    for (size_t j = 1; j <= m; ++j) {
        score.at(0, j) = penalties.open_gap + j * penalties.extend_gap;
        deletions.at(0, j) = score.at(0, j);
        inserts.at(0, j) = score.at(0, j);
    }

    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = 1; j <= m; ++j) {
            float match = long_seq[i - 1] == short_seq[j - 1] ? 1 : penalties.mismatch;
            score.at(i, j) = max(score.at(i - 1, j - 1),
                                 deletions.at(i - 1, j - 1),
                                 inserts.at(i - 1, j - 1))
                            + match;
            deletions.at(i, j) = std::max(score.at(i - 1, j) + penalties.open_gap,
                                          deletions.at(i - 1, j))
                            + penalties.extend_gap;
            inserts.at(i, j) = std::max(score.at(i, j - 1) + penalties.open_gap,
                                        inserts.at(i, j - 1))
                            + penalties.extend_gap;
        }
    }

    size_t best_i = 0;
    float best_score = -10000;
    for (size_t i = 0; i <= n; ++i) {
        if (best_score < max(score.at(i, m),
                             deletions.at(i, m),
                             inserts.at(i, m))) {
            best_score = max(score.at(i, m),
                             deletions.at(i, m),
                             inserts.at(i, m));
            best_i = i;
        }
    }

    if (result != nullptr) {
        construct_consensus(long_seq, short_seq,
                            score, deletions,  inserts, best_i,
                            penalties, *result);
    }
    return best_score;
}


std::string align_to_genes(std::string const& short_seq, pair_vec<std::string> const& genes,
                           alignment_args const& penalties, alignment_res& result,
                           float& score) {
    std::pair<std::string, std::string> const* best_gene = nullptr;
    score = -1000;

    for (auto const& entry : genes) {
        std::string const& gene_name = entry.first;
        std::string const& gene_seq = entry.second;
        
        float current = fitting_alignment(gene_seq, short_seq, penalties, nullptr);
        if (current > score) {
            score = current;
            best_gene = &entry;
        }
    }

    fitting_alignment(best_gene->second, short_seq, penalties, &result);
    return best_gene->first;
}


void align_all(pair_vec<std::string> const& novel,
               pair_vec<std::string> const& genes,
               std::ostream& fout,
               alignment_args const& penalties) {
    for (auto const& entry : novel) {
        alignment_res result;
        float score;
        std::string gene = align_to_genes(entry.second, genes, penalties, result, score);
        
        fout << entry.first << '\t' << gene << '\t' << score
             << '\t' << result.mismatches << '\t' << result.indels << '\t';
        if (result.consensus.empty()) {
            fout << '*';
        } else {
            for (size_t i = 0; i < result.consensus.size(); ++i) {
                if (i) {
                    fout << ',';
                }
                fout << result.consensus[i].first << ':' << result.consensus[i].second;
            }
        }
        fout << '\n';
        fout << std::flush;
    }
}


int main(int argc, char *argv[]) {
    // Arguments:
    // stdin -- cropped sequences fasta, stdout -- output
    // 1 -- germline genes
    // 2 -- mismatch
    // 3 -- open_gap
    // 4 -- extend_gap
    
    pair_vec<std::string> novel = read_fasta(std::cin);
    std::ifstream genes_f(argv[1]);
    pair_vec<std::string> genes = read_fasta(genes_f);

    alignment_args penalties;
    penalties.mismatch = atof(argv[2]);
    penalties.open_gap = atof(argv[3]);
    penalties.extend_gap = atof(argv[4]);

    align_all(novel, genes, std::cout, penalties);
    return 0;
}

