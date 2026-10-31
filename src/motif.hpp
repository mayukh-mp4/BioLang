#pragma once
#include <array>
#include <bitset>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

// How find_all scans: Auto picks per motif; the others force one strategy (for tests/benchmarks).
enum class ScanMode { Auto, Restart, Linear };
extern ScanMode scanMode;

// Thrown when a motif pattern is malformed. `pos` is the index into the pattern.
struct MotifError {
    std::string msg;
    size_t pos;
};

// A compiled motif: a minimised DFA over compressed alphabet classes.
struct Dfa {
    std::string pattern;
    bool shortest = false;                      // lazy quantifier present -> shortest match
    std::array<int, 256> classMap{};            // byte -> alphabet class
    int numClasses = 0;
    std::vector<std::bitset<256>> classChars;   // members of each class
    std::vector<int> trans;                     // numStates * numClasses, -1 = dead
    std::vector<char> accepting;
    int start = 0;
    std::array<char, 256> startByte{};          // bytes that can begin a match
    size_t nfaStates = 0;
    size_t dfaStatesRaw = 0;                    // before minimisation

    // Reverse unanchored DFA: after reading s[i] (right to left) it accepts iff a match starts at i.
    // Built on first use, because it can be much larger than the forward DFA.
    mutable std::vector<int> rtrans;
    mutable std::vector<char> racc;
    mutable int reverseState = 0;               // 0 = not built, 1 = built, -1 = too large
    mutable size_t revStatesRaw = 0;

    int numStates() const { return static_cast<int>(accepting.size()); }

    // End of the match anchored at `from` (longest, or shortest if lazy), or -1.
    long long matchEnd(const std::string& seq, size_t from, size_t* steps = nullptr) const;
    bool ensureReverse() const;
    std::vector<long long> scanLinear(const std::string& seq) const;
    std::vector<long long> scanRestart(const std::string& seq, size_t budget, bool& exceeded) const;

    // Start positions of non-overlapping, leftmost matches in `seq`.
    std::vector<long long> findAll(const std::string& seq) const;

    // Human-readable description of the automaton (used by --dfa).
    void dump(std::ostream& os) const;
};

// pattern -> parse -> Thompson NFA -> subset construction -> minimised DFA.
// Throws MotifError on a malformed pattern.
std::shared_ptr<Dfa> compileMotif(const std::string& pattern);
