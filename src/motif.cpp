#include "motif.hpp"
#include <algorithm>
#include <cstdio>
#include <map>
#include <ostream>
#include <unordered_map>

namespace {

// ============================================================
// 1. Motif pattern -> syntax tree
// ============================================================
struct MNode {
    enum Kind { EMPTY, CHARS, CAT, ALT, STAR, PLUS, OPT, REPEAT } kind;
    std::bitset<256> set;           // CHARS
    std::unique_ptr<MNode> a, b;    // children
    int min = 0, max = -1;          // REPEAT (max = -1 means unbounded)
    explicit MNode(Kind k) : kind(k) {}
};
using MPtr = std::unique_ptr<MNode>;

const int MAX_REPEAT = 1000;

// Codon classes usable inside a motif: <codon>, <sense>, <stop> (DNA alphabet).
// <sense> is every codon except the three stop codons TAA, TAG and TGA.
const std::map<std::string, std::string>& codonClasses() {
    static const std::map<std::string, std::string> table = {
        {"codon", "[ACGT][ACGT][ACGT]"},
        {"stop",  "(TAA|TAG|TGA)"},
        {"sense", "([ACG][ACGT][ACGT]|T([CT][ACGT]|A[CT]|G[CGT]))"},
    };
    return table;
}

class MotifParser {
public:
    explicit MotifParser(const std::string& s) : s_(s) {}
    bool lazy = false;  // set if any lazy quantifier (`*?`, `+?`, `{n,m}?` ...) is used

    MPtr parse() {
        MPtr n = parseAlt();
        if (more()) fail("unmatched ')'");
        return n;
    }

private:
    const std::string& s_;
    size_t pos_ = 0;

    bool more() const { return pos_ < s_.size(); }
    char cur() const { return s_[pos_]; }
    [[noreturn]] void fail(const std::string& m) const { throw MotifError{m, pos_}; }

    static MPtr make(MNode::Kind k) { return std::make_unique<MNode>(k); }
    static MPtr binary(MNode::Kind k, MPtr a, MPtr b) {
        MPtr n = make(k);
        n->a = std::move(a);
        n->b = std::move(b);
        return n;
    }

    // alt = cat { '|' cat }
    MPtr parseAlt() {
        MPtr left = parseCat();
        while (more() && cur() == '|') {
            pos_++;
            left = binary(MNode::ALT, std::move(left), parseCat());
        }
        return left;
    }

    // cat = { repeat }
    MPtr parseCat() {
        MPtr left;
        while (more() && cur() != '|' && cur() != ')') {
            MPtr r = parseRepeat();
            left = left ? binary(MNode::CAT, std::move(left), std::move(r)) : std::move(r);
        }
        return left ? std::move(left) : make(MNode::EMPTY);
    }

    // repeat = atom [ ('*' | '+' | '?' | '{n}' | '{n,}' | '{n,m}') ['?'] ]
    MPtr parseRepeat() {
        MPtr atom = parseAtom();
        if (!more()) return atom;

        MPtr n;
        char c = cur();
        if (c == '*')      { pos_++; n = make(MNode::STAR); }
        else if (c == '+') { pos_++; n = make(MNode::PLUS); }
        else if (c == '?') { pos_++; n = make(MNode::OPT); }
        else if (c == '{') { n = parseBraces(); }
        else return atom;

        n->a = std::move(atom);
        if (more() && cur() == '?') { pos_++; lazy = true; }
        if (more() && (cur() == '*' || cur() == '+' || cur() == '?' || cur() == '{'))
            fail("nothing to repeat");
        return n;
    }

    int readNumber() {
        if (!more() || cur() < '0' || cur() > '9') fail("expected a number in '{...}'");
        long v = 0;
        while (more() && cur() >= '0' && cur() <= '9') {
            v = v * 10 + (cur() - '0');
            pos_++;
            if (v > MAX_REPEAT) fail("repeat count too large (max 1000)");
        }
        return static_cast<int>(v);
    }

    MPtr parseBraces() {
        pos_++;  // '{'
        int mn = readNumber(), mx = mn;
        if (more() && cur() == ',') {
            pos_++;
            mx = (more() && cur() >= '0' && cur() <= '9') ? readNumber() : -1;
        }
        if (!more() || cur() != '}') fail("missing '}'");
        pos_++;
        if (mx != -1 && mx < mn) fail("invalid repeat range");
        MPtr n = make(MNode::REPEAT);
        n->min = mn;
        n->max = mx;
        return n;
    }

    static char unescape(char e) {
        switch (e) {
            case 'n': return '\n';
            case 't': return '\t';
            default:  return e;
        }
    }

    MPtr charNode(unsigned char c) {
        MPtr n = make(MNode::CHARS);
        n->set.set(c);
        return n;
    }

    MPtr parseAtom() {
        char c = cur();
        if (c == '(') {
            pos_++;
            MPtr n = parseAlt();
            if (!more() || cur() != ')') fail("missing ')'");
            pos_++;
            return n;
        }
        if (c == '[') return parseClass();
        if (c == '<') {
            size_t close = s_.find('>', pos_);
            if (close == std::string::npos) fail("missing '>' after '<'");
            std::string name = s_.substr(pos_ + 1, close - pos_ - 1);
            auto it = codonClasses().find(name);
            if (it == codonClasses().end())
                fail("unknown codon class <" + name + "> (use <codon>, <sense> or <stop>)");
            pos_ = close + 1;
            MotifParser sub(it->second);
            return sub.parse();
        }
        if (c == '.') {
            pos_++;
            MPtr n = make(MNode::CHARS);
            n->set.set();
            return n;
        }
        if (c == '\\') {
            pos_++;
            if (!more()) fail("pattern ends with '\\'");
            char e = cur();
            pos_++;
            return charNode(static_cast<unsigned char>(unescape(e)));
        }
        if (c == '*' || c == '+' || c == '?' || c == '{') fail("nothing to repeat");
        pos_++;
        return charNode(static_cast<unsigned char>(c));
    }

    // '[' ['^'] item { item } ']'     item = char | char '-' char
    MPtr parseClass() {
        size_t open = pos_;
        pos_++;  // '['
        bool negate = false;
        if (more() && cur() == '^') { negate = true; pos_++; }

        std::bitset<256> set;
        bool any = false;
        while (more() && cur() != ']') {
            unsigned char lo = readClassChar();
            unsigned char hi = lo;
            if (more() && cur() == '-' && pos_ + 1 < s_.size() && s_[pos_ + 1] != ']') {
                pos_++;
                hi = readClassChar();
                if (hi < lo) fail("invalid range in character class");
            }
            for (int ch = lo; ch <= hi; ch++) set.set(ch);
            any = true;
        }
        if (!more()) { pos_ = open; fail("missing ']'"); }
        if (!any) fail("empty character class");
        pos_++;  // ']'
        if (negate) set.flip();
        MPtr n = make(MNode::CHARS);
        n->set = set;
        return n;
    }

    unsigned char readClassChar() {
        char c = cur();
        pos_++;
        if (c == '\\') {
            if (!more()) fail("pattern ends with '\\'");
            c = unescape(cur());
            pos_++;
        }
        return static_cast<unsigned char>(c);
    }
};

// ============================================================
// 2. Thompson's construction: syntax tree -> NFA
// ============================================================
struct NState {
    bool isChar = false;
    std::bitset<256> chars;   // consumed characters (if isChar)
    int to = -1;
    std::vector<int> eps;     // epsilon transitions
};
struct Frag { int start, end; };

class NfaBuilder {
public:
    std::vector<NState> st;

    Frag build(const MNode& n) {
        switch (n.kind) {
            case MNode::EMPTY: {
                int s = newState(), e = newState();
                st[s].eps.push_back(e);
                return {s, e};
            }
            case MNode::CHARS: {
                int s = newState(), e = newState();
                st[s].isChar = true;
                st[s].chars = n.set;
                st[s].to = e;
                return {s, e};
            }
            case MNode::CAT: {
                Frag a = build(*n.a), b = build(*n.b);
                st[a.end].eps.push_back(b.start);
                return {a.start, b.end};
            }
            case MNode::ALT: {
                Frag a = build(*n.a), b = build(*n.b);
                int s = newState(), e = newState();
                st[s].eps = {a.start, b.start};
                st[a.end].eps.push_back(e);
                st[b.end].eps.push_back(e);
                return {s, e};
            }
            case MNode::STAR: return star(build(*n.a));
            case MNode::OPT:  return opt(build(*n.a));
            case MNode::PLUS: {
                Frag a = build(*n.a);
                int e = newState();
                st[a.end].eps = {a.start, e};
                return {a.start, e};
            }
            case MNode::REPEAT: return repeat(n);
        }
        return {0, 0};
    }

private:
    int newState() { st.emplace_back(); return static_cast<int>(st.size()) - 1; }

    Frag star(Frag a) {
        int s = newState(), e = newState();
        st[s].eps = {a.start, e};
        st[a.end].eps = {a.start, e};
        return {s, e};
    }
    Frag opt(Frag a) {
        int s = newState(), e = newState();
        st[s].eps = {a.start, e};
        st[a.end].eps.push_back(e);
        return {s, e};
    }

    // x{n,m}: n mandatory copies followed by (m-n) optional copies, or a star.
    Frag repeat(const MNode& n) {
        Frag result{-1, -1};
        auto append = [&](Frag f) {
            if (result.start < 0) result = f;
            else { st[result.end].eps.push_back(f.start); result.end = f.end; }
        };
        for (int i = 0; i < n.min; i++) append(build(*n.a));
        if (n.max == -1) append(star(build(*n.a)));
        else for (int i = n.min; i < n.max; i++) append(opt(build(*n.a)));
        if (result.start < 0) {
            int s = newState(), e = newState();
            st[s].eps.push_back(e);
            return {s, e};
        }
        return result;
    }
};

// ============================================================
// 3. Subset construction + minimisation
// ============================================================
std::vector<int> closure(const std::vector<NState>& st, const std::vector<int>& seeds) {
    std::vector<char> seen(st.size(), 0);
    std::vector<int> stack, out;
    for (int s : seeds) if (!seen[s]) { seen[s] = 1; stack.push_back(s); }
    while (!stack.empty()) {
        int s = stack.back();
        stack.pop_back();
        out.push_back(s);
        for (int t : st[s].eps) if (!seen[t]) { seen[t] = 1; stack.push_back(t); }
    }
    std::sort(out.begin(), out.end());
    return out;
}

const size_t MAX_DFA_STATES = 50000;

}  // namespace

std::shared_ptr<Dfa> compileMotif(const std::string& pattern) {
    // --- parse and build the NFA ---
    MotifParser parser(pattern);
    MPtr ast = parser.parse();
    NfaBuilder nb;
    Frag frag = nb.build(*ast);
    const std::vector<NState>& st = nb.st;

    auto dfa = std::make_shared<Dfa>();
    dfa->pattern = pattern;
    dfa->shortest = parser.lazy;
    dfa->nfaStates = st.size();

    // --- alphabet compression: bytes that no NFA transition can tell apart share a class ---
    std::unordered_map<std::string, int> seenSets;
    std::vector<const std::bitset<256>*> sets;
    for (const NState& s : st)
        if (s.isChar && seenSets.emplace(s.chars.to_string(), 0).second) sets.push_back(&s.chars);

    std::unordered_map<std::string, int> sigToClass;
    std::vector<int> reps;  // one representative byte per class
    for (int b = 0; b < 256; b++) {
        std::string sig(sets.size(), '0');
        for (size_t k = 0; k < sets.size(); k++) if ((*sets[k])[b]) sig[k] = '1';
        auto it = sigToClass.emplace(sig, static_cast<int>(sigToClass.size()));
        if (it.second) { reps.push_back(b); dfa->classChars.emplace_back(); }
        dfa->classMap[b] = it.first->second;
        dfa->classChars[it.first->second].set(b);
    }
    const int K = static_cast<int>(reps.size());
    dfa->numClasses = K;

    // --- subset construction ---
    std::map<std::vector<int>, int> ids;
    std::vector<std::vector<int>> subsets;
    auto addState = [&](const std::vector<int>& s) {
        auto it = ids.find(s);
        if (it != ids.end()) return it->second;
        int id = static_cast<int>(subsets.size());
        ids.emplace(s, id);
        subsets.push_back(s);
        return id;
    };
    addState(closure(st, {frag.start}));

    std::vector<int> trans;
    for (size_t i = 0; i < subsets.size(); i++) {
        if (subsets.size() > MAX_DFA_STATES)
            throw MotifError{"pattern too complex: DFA exceeds 50000 states", 0};
        std::vector<int> cur = subsets[i];
        for (int c = 0; c < K; c++) {
            std::vector<int> moved;
            for (int s : cur)
                if (st[s].isChar && st[s].chars[reps[c]]) moved.push_back(st[s].to);
            trans.push_back(moved.empty() ? -1 : addState(closure(st, moved)));
        }
    }
    const int n = static_cast<int>(subsets.size());
    dfa->dfaStatesRaw = n;

    std::vector<char> acc(n, 0);
    for (int i = 0; i < n; i++)
        acc[i] = std::binary_search(subsets[i].begin(), subsets[i].end(), frag.end);

    // --- minimisation (Moore partition refinement) ---
    std::vector<int> cls(n);
    for (int i = 0; i < n; i++) cls[i] = acc[i] ? 1 : 0;
    size_t count = 0;
    while (true) {
        std::map<std::vector<int>, int> sigIds;
        std::vector<int> next(n);
        for (int i = 0; i < n; i++) {
            std::vector<int> sig(K + 1);
            sig[0] = cls[i];
            for (int c = 0; c < K; c++) {
                int t = trans[i * K + c];
                sig[c + 1] = t < 0 ? -1 : cls[t];
            }
            next[i] = sigIds.emplace(sig, static_cast<int>(sigIds.size())).first->second;
        }
        cls = next;
        if (sigIds.size() == count) break;
        count = sigIds.size();
    }

    // Renumber classes in BFS order from the start state, so the start state is S0.
    std::vector<int> repState(count, -1);
    for (int i = 0; i < n; i++) if (repState[cls[i]] < 0) repState[cls[i]] = i;
    std::vector<int> newId(count, -1);
    std::vector<int> order;
    newId[cls[0]] = 0;
    order.push_back(cls[0]);
    for (size_t q = 0; q < order.size(); q++) {
        int rep = repState[order[q]];
        for (int c = 0; c < K; c++) {
            int t = trans[rep * K + c];
            if (t < 0) continue;
            int tc = cls[t];
            if (newId[tc] < 0) { newId[tc] = static_cast<int>(order.size()); order.push_back(tc); }
        }
    }
    for (int cl : order) {
        int rep = repState[cl];
        dfa->accepting.push_back(acc[rep]);
        for (int c = 0; c < K; c++) {
            int t = trans[rep * K + c];
            dfa->trans.push_back(t < 0 ? -1 : newId[cls[t]]);
        }
    }
    dfa->start = 0;
    for (int b = 0; b < 256; b++)
        dfa->startByte[b] = dfa->trans[dfa->start * K + dfa->classMap[b]] >= 0;
    return dfa;
}

// ============================================================
// 4. Scanning
// ============================================================
// At each start position the DFA is run until it dies; every character costs two
// table lookups and no backtracking. Matches are non-overlapping and leftmost; the
// longest match is taken, or the shortest if the pattern uses a lazy quantifier.
std::vector<long long> Dfa::findAll(const std::string& seq) const {
    std::vector<long long> out;
    const size_t n = seq.size();
    const int K = numClasses;
    size_t i = 0;
    while (i < n) {
        // Fast skip over bytes that cannot start a match.
        while (i < n && !startByte[static_cast<unsigned char>(seq[i])]) i++;
        if (i >= n) break;
        int state = start;
        size_t j = i;
        long long lastEnd = -1;
        while (j < n) {
            state = trans[state * K + classMap[static_cast<unsigned char>(seq[j])]];
            if (state < 0) break;
            j++;
            if (accepting[state]) {
                lastEnd = static_cast<long long>(j);
                if (shortest) break;
            }
        }
        if (lastEnd > static_cast<long long>(i)) {
            out.push_back(static_cast<long long>(i));
            i = static_cast<size_t>(lastEnd);
        } else {
            i++;
        }
    }
    return out;
}

// ============================================================
// 5. Pretty printer for --dfa
// ============================================================
static std::string charName(int c) {
    if (c == '\n') return "\\n";
    if (c == '\t') return "\\t";
    if (c == ' ')  return "space";
    if (c > 32 && c < 127) return std::string(1, static_cast<char>(c));
    char buf[8];
    std::snprintf(buf, sizeof buf, "\\x%02X", c);
    return buf;
}

static std::string rangeString(const std::bitset<256>& s) {
    std::string out;
    int i = 0;
    while (i < 256) {
        if (!s[i]) { i++; continue; }
        int j = i;
        while (j + 1 < 256 && s[j + 1]) j++;
        if (!out.empty()) out += ",";
        if (j == i)          out += charName(i);
        else if (j == i + 1) out += charName(i) + "," + charName(j);
        else                 out += charName(i) + "-" + charName(j);
        i = j + 1;
    }
    return out;
}

void Dfa::dump(std::ostream& os) const {
    size_t largest = 0;
    int largestId = -1;
    for (int c = 0; c < numClasses; c++)
        if (classChars[c].count() > largest) { largest = classChars[c].count(); largestId = c; }

    std::vector<std::string> label(numClasses);
    for (int c = 0; c < numClasses; c++) {
        size_t cnt = classChars[c].count();
        if (cnt == 1) label[c] = rangeString(classChars[c]);
        else if (c == largestId && cnt > 128) label[c] = "other";
        else {
            std::string r = rangeString(classChars[c]);
            if (r.size() > 28) r = r.substr(0, 25) + "...";
            label[c] = "{" + r + "}";
        }
    }

    os << "Motif /" << pattern << "/\n";
    os << "  match rule       : "
       << (shortest ? "shortest match (lazy quantifier present)" : "longest match") << "\n";
    os << "  NFA states       : " << nfaStates << "\n";
    os << "  DFA states       : " << dfaStatesRaw << " -> " << numStates() << " after minimisation\n";
    os << "  alphabet classes : " << numClasses << "\n";
    os << "  transitions:\n";

    const int limit = 60;
    for (int s = 0; s < numStates() && s < limit; s++) {
        os << "    S" << s << (s == start ? " (start)" : "") << (accepting[s] ? " [accept]" : "") << ":";
        // Group alphabet classes that lead to the same target state.
        std::map<int, std::vector<int>> byTarget;
        for (int c = 0; c < numClasses; c++) {
            int t = trans[s * numClasses + c];
            if (t >= 0) byTarget[t].push_back(c);
        }
        bool any = !byTarget.empty();
        for (const auto& kv : byTarget) {
            std::string lab;
            if (static_cast<int>(kv.second.size()) == numClasses) lab = "any";
            else for (int c : kv.second) lab += (lab.empty() ? "" : "|") + label[c];
            os << "  " << lab << "->S" << kv.first;
        }
        if (!any) os << "  (none)";
        os << "\n";
    }
    if (numStates() > limit) os << "    ... " << (numStates() - limit) << " more states\n";
}
