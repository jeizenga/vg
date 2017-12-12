#ifndef VG_UTILITY_HPP_INCLUDED
#define VG_UTILITY_HPP_INCLUDED

#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <signal.h>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <unordered_set>
#include <random>
#include <unistd.h>
#include "vg.pb.h"
#include "sha1.hpp"
#include "Variant.h"

namespace vg {

using namespace std;

char reverse_complement(const char& c);
string reverse_complement(const string& seq);
void reverse_complement_in_place(string& seq);
/// Return True if the given string is entirely Ns of either case, and false
/// otherwise.
bool is_all_n(const string& seq);
int get_thread_count(void);
string wrap_text(const string& str, size_t width);
bool is_number(const string& s);

// split a string on any character found in the string of delimiters (delims)
std::vector<std::string>& split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems);
std::vector<std::string> split_delims(const std::string &s, const std::string& delims);

const std::string sha1sum(const std::string& data);
const std::string sha1head(const std::string& data, size_t head);

bool allATGC(const string& s);
string nonATGCNtoN(const string& s);
// Convert ASCII-encoded DNA to upper case
string toUppercase(const string& s);
double median(std::vector<int> &v);
double stdev(const std::vector<double>& v);

template<typename T>
double stdev(const T& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return std::sqrt(sq_sum / v.size());
}

// Φ is the normal cumulative distribution function
// https://en.wikipedia.org/wiki/Cumulative_distribution_function
double phi(double x1, double x2);
    
/// Inverse CDF of a standard normal distribution. Must have 0 < quantile < 1.
double normal_inverse_cdf(double quantile);

/*
 * Return the log of the sum of two log-transformed values without taking them
 * out of log space.
 */
inline double add_log(double log_x, double log_y) {
    return log_x > log_y ? log_x + log(1.0 + exp(log_y - log_x)) : log_y + log(1.0 + exp(log_x - log_y));
}
    
/*
 * Return the log of the difference of two log-transformed values without taking
 * them out of log space.
 */
inline double subtract_log(double log_x, double log_y) {
    return log_x + log(1.0 - exp(log_y - log_x));
}
 
/**
 * Convert a number ln to the same number log 10.
 */   
inline double ln_to_log10(double ln) {
    return ln / log(10);
}

/**
 * Convert a number log 10 to the same number ln.
 */   
inline double log10_to_ln(double l10) {
    return l10 * log(10);
}
    
// Convert a probability to a natural log probability.
inline double prob_to_logprob(double prob) {
    return log(prob);
}
// Convert natural log probability to a probability
inline double logprob_to_prob(double logprob) {
    return exp(logprob);
}
// Add two probabilities (expressed as logprobs) together and return the result
// as a logprob.
inline double logprob_add(double logprob1, double logprob2) {
    // Pull out the larger one to avoid underflows
    double pulled_out = max(logprob1, logprob2);
    return pulled_out + prob_to_logprob(logprob_to_prob(logprob1 - pulled_out) + logprob_to_prob(logprob2 - pulled_out));
}
// Invert a logprob, and get the probability of its opposite.
inline double logprob_invert(double logprob) {
    return prob_to_logprob(1.0 - logprob_to_prob(logprob));
}

// Convert integer Phred quality score to probability of wrongness.
inline double phred_to_prob(int phred) {
    return pow(10, -((double)phred) / 10);
}

// Convert probability of wrongness to integer Phred quality score.
inline double prob_to_phred(double prob) {
    return -10.0 * log10(prob);
}

// Convert a Phred quality score directly to a natural log probability of wrongness.
inline double phred_to_logprob(int phred) {
    return (-((double)phred) / 10) / log10(exp(1.0));
}

// Convert a natural log probability of wrongness directly to a Phred quality score.
inline double logprob_to_phred(double logprob ) {
    return -10.0 * logprob * log10(exp(1.0));
}

// Take the geometric mean of two logprobs
inline double logprob_geometric_mean(double lnprob1, double lnprob2) {
    return log(sqrt(exp(lnprob1 + lnprob2)));
}

// Same thing in phred
inline double phred_geometric_mean(double phred1, double phred2) {
    return prob_to_phred(sqrt(phred_to_prob(phred1 + phred2)));
}

// normal pdf, from http://stackoverflow.com/a/10848293/238609
template <typename T>
T normal_pdf(T x, T m, T s)
{
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-T(0.5) * a * a);
}

// Emit a stack trace when something bad happens.
void emit_stacktrace(int signalNumber, siginfo_t *signalInfo, void *signalContext);

// This is an internal function used by the above.
string demangle_frame(string mangled);

template<typename T, typename V>
set<T> map_keys_to_set(const map<T, V>& m) {
    set<T> r;
    for (auto p : m) r.insert(p.first);
    return r;
}

// pairwise maximum
template<typename T>
vector<T> pmax(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> c;
    assert(a.size() == b.size());
    c.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(c),
                   [](T a, T b) { return std::max<T>(a, b); });
    return c;
}

// maximum of all vectors
template<typename T>
vector<T> vpmax(const std::vector<std::vector<T>>& vv) {
    std::vector<T> c;
    if (vv.empty()) return c;
    c = vv.front();
    typename std::vector<std::vector<T> >::const_iterator v = vv.begin();
    ++v; // skip the first element
    for ( ; v != vv.end(); ++v) {
        c = pmax(c, *v);
    }
    return c;
}

/**
 * Compute the sum of the values in a collection. Values must be default-
 * constructable (like numbers are).
 */
template<typename Collection>
typename Collection::value_type sum(const Collection& collection) {

    // Set up an alias
    using Item = typename Collection::value_type;

    // Make a new zero-valued item to hold the sum
    auto total = Item();
    for(auto& to_sum : collection) {
        total += to_sum;
    }

    return total;

}

/**
 * Compute the sum of the values in a collection, where the values are log
 * probabilities and the result is the log of the total probability. Items must
 * be convertible to/from doubles for math.
 */
template<typename Collection>
typename Collection::value_type logprob_sum(const Collection& collection) {

    // Set up an alias
    using Item = typename Collection::value_type;

    // Pull out the minimum value
    auto min_iterator = min_element(begin(collection), end(collection));

    if(min_iterator == end(collection)) {
        // Nothing there, p = 0
        return Item(prob_to_logprob(0));
    }

    auto check_iterator = begin(collection);
    ++check_iterator;
    if(check_iterator == end(collection)) {
        // We only have a single element anyway. We don't want to subtract it
        // out because we'll get 0s.
        return *min_iterator;
    }

    // Pull this much out of every logprob.
    Item pulled_out = *min_iterator;

    if(logprob_to_prob(pulled_out) == 0) {
        // Can't divide by 0!
        // TODO: fix this in selection
        pulled_out = prob_to_logprob(1);
    }

    Item total(0);
    for(auto& to_add : collection) {
        // Sum up all the scaled probabilities.
        total += logprob_to_prob(to_add - pulled_out);
    }

    // Re-log and re-scale
    return pulled_out + prob_to_logprob(total);
}

/// Find the system temp directory using defaults and environment variables
string find_temp_dir();

/// Create a temporary file starting with the given base name
string tmpfilename(const string& base);

/// Create a temporary file in the appropriate system temporary directory
string tmpfilename();

// Code to detect if a variant lacks an ID and give it a unique but repeatable
// one.
string get_or_make_variant_id(const vcflib::Variant& variant);
string make_variant_id(const vcflib::Variant& variant);

// TODO: move these to genotypekit on a VCF emitter?

/**
 * Create the reference allele for an empty vcflib Variant, since apaprently
 * there's no method for that already. Must be called before any alt alleles are
 * added.
 */
void create_ref_allele(vcflib::Variant& variant, const std::string& allele);

/**
 * Add a new alt allele to a vcflib Variant, since apaprently there's no method
 * for that already.
 *
 * If that allele already exists in the variant, does not add it again.
 *
 * Retuerns the allele number (0, 1, 2, etc.) corresponding to the given allele
 * string in the given variant. 
 */
int add_alt_allele(vcflib::Variant& variant, const std::string& allele);

/**
 * We have a transforming map function that we can chain.
 */ 
template <template <class T, class A = std::allocator<T>> class Container, typename Input, typename Output>
Container<Output> map_over(const Container<Input>& in, const std::function<Output(const Input&)>& lambda) {
    Container<Output> to_return;
    for (const Input& item : in) {
        to_return.push_back(lambda(item));
    }
    return to_return;
}

/**
 * We have a wrapper of that to turn a container reference into a container of pointers.
 */
template <template <class T, class A = std::allocator<T>> class Container, typename Item>
Container<const Item*> pointerfy(const Container<Item>& in) {
    return map_over<Container, Item, const Item*>(in, [](const Item& item) -> const Item* {
        return &item;
    });
}

// Simple little tree
template<typename T>
struct TreeNode {
    T v;
    vector<TreeNode<T>*> children;
    TreeNode<T>* parent;
    TreeNode() : parent(0) {}
    ~TreeNode() { for (auto c : children) { delete c; } }
    void for_each_preorder(function<void(TreeNode<T>*)> lambda) {
        lambda(this);
        for (auto c : children) {
            c->for_each_preorder(lambda);
        }
    }
    void for_each_postorder(function<void(TreeNode<T>*)> lambda) {
        for (auto c : children) {
            c->for_each_postorder(lambda);
        }
        lambda(this);
    }
};
    
/**
 * A custom Union-Find data structure that supports merging a set of indices in
 * disjoint sets in amortized nearly linear time. This implementation also supports
 * querying the size of the group containing an index in constant time and querying
 * the members of the group containing an index in linear time in the size of the group.
 */
class UnionFind {
public:
    /// Construct UnionFind for this many indices
    UnionFind(size_t size);
    
    /// Destructor
    ~UnionFind();
    
    /// Returns the number of indices in the UnionFind
    size_t size();
    
    /// Returns the group ID that index i belongs to (can change after calling union)
    size_t find_group(size_t i);
    
    /// Merges the group containing index i with the group containing index j
    void union_groups(size_t i, size_t j);
    
    /// Returns the size of the group containing index i
    size_t group_size(size_t i);
    
    /// Returns a vector of the indices in the same group as index i
    vector<size_t> group(size_t i);
    
    /// Returns all of the groups, each in a separate vector
    vector<vector<size_t>> all_groups();
    
    /// A string representation of the current state for debugging
    string current_state();
    
private:
    
    struct UFNode;
    vector<UFNode> uf_nodes;
};

struct UnionFind::UFNode {
    UFNode(size_t index) : head(index), rank(0), size(1) {}
    ~UFNode() {}
    
    size_t rank;
    size_t size;
    size_t head;
    unordered_set<size_t> children;
};

    
template<typename T>
struct Tree {
    typedef TreeNode<T> Node;
    Node* root;
    Tree(Node* r = 0) : root(r) { }
    ~Tree() { delete root; }
    void for_each_preorder(function<void(Node*)> lambda) {
        if (root) root->for_each_preorder(lambda);
    }
    void for_each_postorder(function<void(Node*)> lambda) {
       if (root) root->for_each_postorder(lambda);
    }

};
    
vector<size_t> range_vector(size_t begin, size_t end);

struct IncrementIter {
public:
    IncrementIter(size_t number) : current(number) {
        
    }
    
    inline IncrementIter& operator=(const IncrementIter& other) {
        current = other.current;
        return *this;
    }
    
    inline bool operator==(const IncrementIter& other) const {
        return current == other.current;
    }
    
    inline bool operator!=(const IncrementIter& other) const {
        return current != other.current;
    }
    
    inline IncrementIter operator++() {
        current++;
        return *this;
    }
    
    inline IncrementIter operator++( int ) {
        IncrementIter temp = *this;
        current++;
        return temp;
    }
    
    inline size_t operator*(){
        return current;
    }
    
private:
    size_t current;
};
    
size_t integer_power(size_t x, size_t power);

// Get a callback with an istream& to an open file if a file name argument is
// present after the parsed options, or print an error message and exit if one
// is not. Handles "-" as a filename as indicating standard input. The reference
// passed is guaranteed to be valid only until the callback returns. Bumps up
// optind to the next argument if a filename is found.
void get_input_file(int& optind, int argc, char** argv, function<void(istream&)> callback);

// Parse out the name of an input file (i.e. the next positional argument), or
// throw an error. File name must be nonempty, but may be "-" or may not exist.
string get_input_file_name(int& optind, int argc, char** argv);

// Parse out the name of an output file (i.e. the next positional argument), or
// throw an error. File name must be nonempty.
string get_output_file_name(int& optind, int argc, char** argv);

// Get a callback with an istream& to an open file. Handles "-" as a filename as
// indicating standard input. The reference passed is guaranteed to be valid
// only until the callback returns.
void get_input_file(const string& file_name, function<void(istream&)> callback);

double slope(const std::vector<double>& x, const std::vector<double>& y);
double fit_zipf(const vector<double>& y);

/// Computes base^exponent in log(exponent) time
size_t integer_power(uint64_t base, uint64_t exponent);
/// Computes base^exponent mod modulus in log(exponent) time without requiring more
/// than 64 bits to represent exponentiated number
size_t modular_exponent(uint64_t base, uint64_t exponent, uint64_t modulus);

/// Returns a uniformly random DNA sequence of the given length
string random_sequence(size_t length);
    
}

#endif
