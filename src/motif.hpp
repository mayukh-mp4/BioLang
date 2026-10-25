#pragma once
#include <array>
#include <bitset>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

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

    int numStates() const { return static_cast<int>(accepting.size()); }

    // Start positions of non-overlapping, leftmost matches in `seq`.
    std::vector<long long> findAll(const std::string& seq) const;

    // Human-readable description of the automaton (used by --dfa).
    void dump(std::ostream& os) const;
};

// pattern -> parse -> Thompson NFA -> subset construction -> minimised DFA.
// Throws MotifError on a malformed pattern.
std::shared_ptr<Dfa> compileMotif(const std::string& pattern);
