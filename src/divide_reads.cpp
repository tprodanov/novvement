#include<iostream>
#include<fstream>
#include<vector>
#include<tuple>
#include<string>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<unordered_map>
#include<unordered_set>
#include<stdio.h>


template<size_t J, typename ... Ts, typename U>
std::tuple<Ts...> const* tuple_array_search(std::vector<std::tuple<Ts...>> const& arr, U const& x) {
    for (auto const& v : arr) {
        if (std::get<J>(v) == x) {
            return &v;
        }
    }
    return nullptr;
}


std::vector<std::string> split(std::string str, char delimiter) {
    std::vector<std::string> res;
    size_t start = 0;
    size_t pos = 0;
    while ((pos = str.find(delimiter, start)) != std::string::npos) {
        res.push_back(str.substr(start, pos - start));
        start = pos + 1;
    }
    res.push_back(str.substr(start));
    return res;
}


using error_tup = std::tuple<size_t, std::string>;
using errors_vec = std::vector<error_tup>;


struct error_hash : public std::unary_function<error_tup, size_t> {
    size_t operator()(error_tup const& key) const {
        return std::get<0>(key) * 719495873 + hasher_(std::get<1>(key)) * 1134029473;
    }

private:
    std::hash<std::string> hasher_;
};

template<typename T>
using errors_map = std::unordered_map<error_tup, T, error_hash>;

using errors_pair = std::pair<error_tup, error_tup>;
struct errors_pair_hash : public std::unary_function<errors_pair, size_t> {
    size_t operator()(errors_pair const& key) const {
        return 710503121 * hasher_(key.first) + 75219103 * hasher_(key.second);
    }

private:
    error_hash hasher_;
};

template<typename T>
using errors_pair_map = std::unordered_map<errors_pair, T, errors_pair_hash>;


template<typename T, typename U>
struct compare_by_second {
    bool operator()(std::pair<T, U> const& a, std::pair<T, U> const& b) const {
        return a.second < b.second;
    }
};


template<typename ... Ts>
bool contains(std::vector<std::tuple<Ts...>> const& arr, std::tuple<Ts...> const& el) {
    for (auto const& x : arr) {
        if (x == el) {
            return true;
        }
    }
    return false;
}


void write_errors_vec(errors_vec const& errors, std::ostream& outp) {
    for (size_t i = 0; i < errors.size(); ++i) {
        outp << std::get<0>(errors[i]) << ':' << std::get<1>(errors[i]);
        if (i < errors.size() - 1) {
            outp << ',';
        }
    }
    if (!errors.size()) {
        outp << '*';
    }
}


struct read {
    read(std::string const& name, std::string const& errors, std::string const& label)
        : name_(name)
        , label_(label) {
        if (errors != "*") {
            for (auto const& error : split(errors, ',')) {
                auto error_split = split(error, ':');
                original_errors_.push_back(std::make_tuple(stoi(error_split[0]),
                                                           error_split[1]));
            }
        }
        update(errors_vec());
    }

    void update(errors_vec const& consensus) {
        current_errors_.clear();
        for (auto const& entry : consensus) {
            if (tuple_array_search<0>(original_errors_, std::get<0>(entry)) == nullptr) {
                current_errors_.push_back(std::make_tuple(std::get<0>(entry), "*"));
            }
        }

        for (auto const& error : original_errors_) {
            if (!contains(consensus, error)) {
                current_errors_.push_back(error);
            }
        }
        std::sort(current_errors_.begin(), current_errors_.end());
    }

    void write_to(std::ostream& outp, std::string const& segment, std::string const& subset) const {
        outp << name_ << '\t' << segment << '\t';
        write_errors_vec(original_errors_, outp);
        outp << '\t' << label_ << '\t' << subset << '\n';
    }

    errors_vec const& get_original_errors() const {
        return original_errors_;
    }

    errors_vec const& get_current_errors() const {
        return current_errors_;
    }

    std::string const& get_label() const {
        return label_;
    }

private:

    std::string name_;
    std::string label_;
    errors_vec original_errors_;
    errors_vec current_errors_;

};


std::vector<double> small_log_factorials() {
    unsigned long long f = 1;
    std::vector<double> res;
    res.push_back(0);
    for (unsigned long long i = 1; i <= 20; ++i) {
        f *= i;
        res.push_back(log(f));
    }
    return res;
}


double log_factorial(size_t n) {
    static std::vector<double> values = small_log_factorials();
    if (n >= values.size()) {
        for (size_t i = values.size(); i <= n; ++i) {
            values.push_back(values[values.size() - 1] + log(i));
        }
    }
    return values[n];
}


errors_vec consensus(std::vector<read> const& reads) {
    errors_map<size_t> counter;
    for (read const& r : reads) {
        for (auto const& err : r.get_original_errors()) {
            counter[err] += 1;
        }
    }
    size_t s = reads.size() / 2;


    errors_vec res;
    for (auto const& entry : counter) {
        if (entry.second > s) {
            res.push_back(entry.first);
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}


double significance(size_t n, size_t n1, size_t n2, size_t n12) {
    double base = log_factorial(n1) + log_factorial(n - n1) + log_factorial(n2) + log_factorial(n - n2) - log_factorial(n);

    double s = 0;
    for (size_t i = n12; i <= std::min(n1, n2); ++i) {
        double current = log_factorial(i) + log_factorial(n1 - i) + log_factorial(n2 - i) + log_factorial(n - n1 - n2 + i);
        s += std::exp(base - current);
    }
    return s;
}


struct reads_subset {
    reads_subset(std::string const& segment, std::vector<read>&& reads)
        : segment_(segment)
        , reads_(std::move(reads))
        , consensus_(consensus(reads_)) {
        for (read& r : reads_) {
            r.update(consensus_);
        }
    }

    std::vector<std::pair<errors_pair, size_t>> top_pairs(size_t min_subset_size, size_t n_pairs) const {
        errors_pair_map<size_t> counter;
        for (auto const& r : reads_) {
            std::vector<error_tup> const& r_errors = r.get_current_errors();
            for (size_t i = 0; i < r_errors.size(); ++i) {
                for (size_t j = i + 1; j < r_errors.size(); ++j) {
                    counter[std::make_pair(r_errors[i], r_errors[j])] += 1;
                }
            }
        }

        std::vector<std::pair<errors_pair, size_t>> elems(counter.begin(), counter.end());
        std::sort(elems.rbegin(), elems.rend(), compare_by_second<errors_pair, size_t>());

        std::vector<std::pair<errors_pair, size_t>> res;
        for (auto const& entry : elems) {
            if (entry.second < min_subset_size || reads_.size() - entry.second < min_subset_size) {
                continue;
            }
            res.push_back(entry);
            if (res.size() >= n_pairs) {
                return res;
            }
        }
        return res;
    }

    bool top_peak(float peak_rate, float noise_rate,
                  error_tup& peak, size_t& error_count) const {
        errors_map<size_t> counter;
        for (auto const& r : reads_) {
            for (auto const& err : r.get_current_errors()) {
                counter[err] += 1;
            }
        }

        float total = reads_.size();

        bool found = false;
        for (auto const& entry : counter) {
            if (entry.second / total >= peak_rate && entry.second / total <= 1 - peak_rate) {
                if (found) {
                    return false;
                }
                found = true;
                peak = entry.first;
                error_count = entry.second;
            } else if (entry.second / total >= noise_rate && entry.second / total <= 1 - noise_rate) {
                return false;
            }
        }
        return found;
    }

    size_t count_reads_with(error_tup const& err) const {
        size_t count = 0;
        for (read const& r : reads_) {
            if (contains(r.get_current_errors(), err)) {
                count += 1;
            }
        }
        return count;
    }

    std::vector<reads_subset> divide_if_possible(double max_significance, size_t min_subset_size, size_t n_pairs) {
        auto pairs = top_pairs(min_subset_size, n_pairs);
        std::vector<reads_subset> res;

        for (auto const& entry : pairs) {
            error_tup const& err1 = entry.first.first;
            error_tup const& err2 = entry.first.second;
            size_t n12 = entry.second;

            size_t n1 = count_reads_with(err1);
            size_t n2 = count_reads_with(err2);

            double s = significance(reads_.size(), n1, n2, n12);
            if (s > max_significance / pairs.size()) {
                continue;
            }

            std::vector<read> reads1;
            std::vector<read> reads2;
            divide_by(err1, err2, reads1, reads2);

            res.emplace_back(segment_, std::move(reads1));
            res.emplace_back(segment_, std::move(reads2));
            return res;
        }
        return res;
    }

    std::vector<reads_subset> divide_if_possible_1_peak(float peak_rate, float noise_rate, size_t min_subset_size) {
        error_tup peak;
        size_t peak_count;
        std::vector<reads_subset> res;

        if (!top_peak(peak_rate, noise_rate, peak, peak_count)) {
            return res;
        }

        if (peak_count < min_subset_size || reads_.size() - peak_count < min_subset_size) {
            return res;
        }

        std::vector<read> reads1;
        std::vector<read> reads2;
        for (size_t i = 0; i < reads_.size(); ++i) {
            if (contains(reads_[i].get_current_errors(), peak)) {
                reads1.push_back(std::move(reads_[i]));
            } else {
                reads2.push_back(std::move(reads_[i]));
            }
        }

        res.emplace_back(segment_, std::move(reads1));
        res.emplace_back(segment_, std::move(reads2));
        return res;
    }

    size_t different_labels() const {
        std::unordered_set<size_t> lengths;
        for (auto const& r : reads_) {
            lengths.insert(r.get_label().size());
        }
        return lengths.size();
    }

    void write(size_t index, size_t total_reads, std::ostream& reads_outp, std::ostream& summary_outp) const {
        char name[50];
        sprintf(name, "%s:%zd", segment_.c_str(), index);

        size_t labels = different_labels();
        char long_name[100];
        sprintf(long_name, "\"%s:%2zd. Coverage: %zd. Labels: %zd\"", segment_.c_str(), index, reads_.size(), labels);

        for (read const& r : reads_) {
            r.write_to(reads_outp, segment_, long_name);
        }

        summary_outp << name << '\t' << reads_.size() << '\t' << static_cast<float>(reads_.size()) / total_reads
                     << '\t' << labels << '\t';
        write_errors_vec(consensus_, summary_outp);
        summary_outp << '\n';
    }

    struct subset_compare {
        bool operator()(reads_subset const& a, reads_subset const& b) {
            return a.reads_.size() < b.reads_.size();
        }
    };

private:

    void divide_by(error_tup const& err1, error_tup const& err2, std::vector<read>& reads1, std::vector<read>& reads2) {
        for (size_t i = 0; i < reads_.size(); ++i) {
            if (contains(reads_[i].get_current_errors(), err1) && contains(reads_[i].get_current_errors(), err2)) {
                reads1.push_back(std::move(reads_[i]));
            } else {
                reads2.push_back(std::move(reads_[i]));
            }
        }
    }

    std::string segment_;
    std::vector<read> reads_;
    errors_vec consensus_;

};


std::vector<reads_subset> divide_while_possible(std::string const& segment, std::vector<read>&& reads,
                                    double max_significance, size_t min_subset_size, size_t n_pairs,
                                    float peak_rate, float noise_rate) {
    std::vector<reads_subset> subsets { reads_subset(segment, std::move(reads)) };
    std::vector<reads_subset> res;
    while (subsets.size()) {
        reads_subset& subset = subsets[subsets.size() - 1];
        auto division = subset.divide_if_possible(max_significance, min_subset_size, n_pairs);
        if (division.size()) {
            subsets.pop_back();
            subsets.push_back(std::move(division[0]));
            subsets.push_back(std::move(division[1]));

        } else {
            division = subset.divide_if_possible_1_peak(peak_rate, noise_rate, min_subset_size);
            if (division.size()) {
                res.push_back(std::move(division[0]));
                res.push_back(std::move(division[1]));
            } else {
                res.push_back(std::move(subset));
            }
            subsets.pop_back();
        }
    }
    return res;
}


void write_subsets(std::vector<reads_subset> const& subsets, size_t total_reads,
                   std::ostream& reads_outp, std::ostream& summary_outp) {
    for (size_t i = 0; i < subsets.size(); ++i) {
        subsets[i].write(i + 1, total_reads, reads_outp, summary_outp);
    }
}


std::unordered_map<std::string, std::vector<read>> load_errors(std::istream& inp) {
    std::string current;
    do {
        std::getline(inp, current);
    } while (current[0] == '#');

    std::unordered_map<std::string, std::vector<read>> reads;
    while (inp) {
        std::getline(inp, current);
        if (!current.size()) {
            break;
        }
        auto line_split = split(current, '\t');

        std::string name = line_split[0];
        std::string segment = line_split[1];
        std::string errors = line_split[2];
        std::string label = line_split[3];

        reads[segment].emplace_back(name, errors, label);
    }
    return reads;
}


size_t count_reads(std::unordered_map<std::string, std::vector<read>> const& reads) {
    size_t count = 0;
    for (auto const& entry : reads) {
        count += entry.second.size();
    }
    return count;
}


int main(int argc, char *argv[]) {
    // Arguments:
    // stdin -- input reads
    // stdout -- output reads
    // stderr -- output summary
    // 1 -- significance
    // 2 -- pairs
    // 3 -- coverage
    // 4 -- noise_rate
    // 5 -- peak_rate

    double max_significance = std::stof(argv[1]);
    size_t n_pairs = std::stoi(argv[2]);
    size_t min_subset_size = std::stoi(argv[3]);
    float noise_rate = std::stof(argv[4]);
    float peak_rate = std::stof(argv[5]);

    auto reads = load_errors(std::cin);
    size_t total_reads = count_reads(reads);

    for (auto entry : reads) {
        auto subsets = divide_while_possible(entry.first, std::move(entry.second), max_significance, min_subset_size, n_pairs,
                                             peak_rate, noise_rate);
        std::sort(subsets.rbegin(), subsets.rend(), reads_subset::subset_compare());
        write_subsets(subsets, total_reads, std::cout, std::cerr);
    }
}
