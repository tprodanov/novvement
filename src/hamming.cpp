#include <fstream>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>


struct Vertex {

    static bool diff_lengths_accepted;
    
    Vertex(std::string const& seq) {
        this->seq = seq;
    }

    void add_read(std::string const& read_name) {
        reads.push_back(read_name);       
    }

    size_t coverage() const {
        return reads.size();
    }

    size_t distance(std::string const& other) const {
        size_t res = 0;
        if (seq.length() != other.length()) {
            if (!Vertex::diff_lengths_accepted) {
                return 1000;
            } else {
                res += std::max(seq.length(), other.length())
                     - std::min(seq.length(), other.length());
            }
        }
        
        for (size_t i = 0; i < std::min(seq.length(), other.length()); ++i) {
            if (seq[i] != other[i]) {
                res += 1;
            }
        }
        return res;
    }

    void add_edge(Vertex& other) {
        edges.insert(&other);
        other.edges.insert(this);
    }

    std::string const& get_seq() const {
        return seq;
    }

    std::unordered_set<Vertex const*> const& get_edges() const {
        return edges;
    }

    std::vector<std::string> const& get_reads() const {
        return reads;
    }

private:

    std::string seq;
    std::unordered_set<Vertex const*> edges;
    std::vector<std::string> reads;

};

bool Vertex::diff_lengths_accepted = false;


bool compare_vertices(Vertex const* first, Vertex const* second) {
    return first->coverage() > second->coverage();
}


struct HammingGraph {

    HammingGraph(size_t tau)
            : tau(tau) {
    }

    ~HammingGraph() {
        for (auto& item: vertices) {
            delete item.second;
        }
    }

    void add_read(std::string const& read_name, std::string const& seq) {
        auto it = vertices.find(seq);
        if (it != vertices.end()) {
            it->second->add_read(read_name);
            return;
        }

        Vertex* current = new Vertex(seq);
        current->add_read(read_name);

        for (auto& item : vertices) {
            if (item.second->distance(current->get_seq()) <= tau) {
                item.second->add_edge(*current);
            }
        }

        vertices[seq] = current;
    }

    std::vector<Vertex const*> bfs(Vertex const* start, std::unordered_set<Vertex const*>& visited) const {
        std::queue<Vertex const*> to_traverse;
        to_traverse.push(start);
        visited.insert(start);

        std::vector<Vertex const*> result;

        while (to_traverse.size()) {
            Vertex const* current = to_traverse.front();
            to_traverse.pop();
            result.push_back(current);

            for (Vertex const* neigh : current->get_edges()) {
                if (visited.count(neigh) == 0) {
                    visited.insert(neigh);
                    to_traverse.push(neigh);
                }
            }
        }
        return result;
    }

    std::vector<std::vector<Vertex const*>> components() const {
        std::vector<std::vector<Vertex const*>> result;
        std::unordered_set<Vertex const*> visited;
        
        for (auto const& item : vertices) {
            Vertex const* vertex = item.second;
            if (visited.count(vertex) == 0) {
                std::vector<Vertex const*> component = bfs(vertex, visited);
                std::sort(component.begin(), component.end(), compare_vertices);
                result.push_back(std::move(component));
            }
        }
        return result;
    }

private:

    std::unordered_map<std::string, Vertex*> vertices;
    size_t tau;

};


void create_graph(std::istream& inp, HammingGraph& graph) {
    std::string read_name;
    std::string read_seq;

    // header
    getline(inp, read_name);

    while (inp >> read_name >> read_seq) {
        graph.add_read(read_name, read_seq);
    }
}


void print_labels(HammingGraph const& graph, std::ostream& outp) {
    outp << "read\tlabel\n";

    for (auto const& component : graph.components()) {
        std::string label = component[0]->get_seq();
        for (auto const& vertex : component) {
            for (auto const& read : vertex->get_reads()) {
                outp << read << '\t' << label << '\n';
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "argv should contain exactly 2 arguments: tau, accept_diff_len" << std::endl;
        exit(1);
    }
    
    size_t tau;
    sscanf(argv[1], "%zu", &tau);

    Vertex::diff_lengths_accepted = std::string(argv[2]) == "True";

    HammingGraph graph(tau);
    create_graph(std::cin, graph);
    print_labels(graph, std::cout);
}

