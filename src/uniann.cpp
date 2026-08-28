#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <array>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <stdexcept>
#include <utility>

using namespace std;

//------------------------------------------------------------
// States:
// 0: N
// 1: E0  2: E1  3: E2
// 4: I0  5: I1  6: I2
//------------------------------------------------------------

static const int NUM_STATES = 7;

static const array<string, NUM_STATES> state_name = {
    "N", "E0", "E1", "E2", "I0", "I1", "I2"
};

static const double NEG_INF = -1e9;
static const int MIN_INTRON = 40;
static const int MIN_EXON   = 3;
static const int MIN_INTER  = 30;
static const int MIN_SINGLE = 100;
static const double SCORE_TOLERANCE = 1e-8;

enum class AnchorType { Donor, Acceptor };
enum class CandidateStatus {
    Complete,
    NoValidUpstreamBoundary,
    NoValidDownstreamBoundary,
    MaximumDistanceReached,
    ScoreValidationFailure,
    AnnotationConversionFailure
};

struct PathMetadata {
    int intron_len = 0;
    int exon_len = 0;
    int inter_len = 0;
    int exon_from = 0;
    // The predecessor of the cell carrying this metadata.  The legacy start
    // rule checks this value to ensure the intergenic run itself came from N.
    int predecessor_state = -1;
};

struct GenomicInterval {
    int start_0based = 0;
    int end_0based_exclusive = 0;
    bool operator==(const GenomicInterval &other) const {
        return start_0based == other.start_0based &&
               end_0based_exclusive == other.end_0based_exclusive;
    }
    bool operator<(const GenomicInterval &other) const {
        return tie(start_0based, end_0based_exclusive) <
               tie(other.start_0based, other.end_0based_exclusive);
    }
};

struct SpliceJunction {
    // Both coordinates are the first intronic bases of their motifs.  The
    // interval is [donor_0based, acceptor_0based + 2).
    int donor_0based = -1;
    int acceptor_0based = -1;
    char strand = '+';
    bool operator==(const SpliceJunction &other) const {
        return donor_0based == other.donor_0based &&
               acceptor_0based == other.acceptor_0based &&
               strand == other.strand;
    }
    bool operator<(const SpliceJunction &other) const {
        return tie(donor_0based, acceptor_0based, strand) <
               tie(other.donor_0based, other.acceptor_0based, other.strand);
    }
};

struct SpliceAnchor {
    string sequence_id;
    char strand = '+';
    int path_index = -1;
    int dp_position = -1;
    // First base of the GT or AG dinucleotide.
    int genomic_position_0based = -1;
    int genomic_position_1based = -1;
    AnchorType type = AnchorType::Donor;
    int left_state = -1;
    int right_state = -1;
    int frame = -1;
    double splice_score = NEG_INF;
    int reference_transcript_index = -1;
};

struct TranscriptAnnotation {
    string sequence_id;
    char strand = '+';
    int genomic_start_0based = -1;
    int genomic_end_0based_exclusive = -1;
    int path_start = -1;
    int path_end = -1;
    vector<GenomicInterval> exons;
    vector<GenomicInterval> cds_intervals;
    vector<SpliceJunction> junctions;
    double score = NEG_INF;
};

struct AlternativeCandidate {
    SpliceAnchor anchor;
    CandidateStatus status = CandidateStatus::NoValidUpstreamBoundary;
    bool upstream_complete = false;
    bool downstream_complete = false;
    bool valid = false;
    bool identical_to_global = false;
    bool kept_after_deduplication = false;
    int local_interval_start = -1;
    int local_interval_end = -1;
    double upstream_score = NEG_INF;
    double anchor_score = NEG_INF;
    double downstream_score = NEG_INF;
    double total_score = NEG_INF;
    double reference_interval_score = NEG_INF;
    double score_delta = NEG_INF;
    string classification = "complex";
    string failure_reason;
    vector<int> path_positions;
    vector<int> path_states;
    TranscriptAnnotation transcript;
};

struct LocalKBestOptions {
    bool enabled = false;
    int k = 5;
    string output_filename;
    string report_filename;
    double start_score_drop = std::numeric_limits<double>::infinity();
    bool no_dp_dump = false;
};

struct AlternativeOptions {
    bool enabled = false;
    string output_filename;
    string report_filename;
    int max_distance = 1000000;
    double min_score_delta = -numeric_limits<double>::infinity();
    bool include_incomplete = false;
    bool debug = false;
};

//------------------------------------------------------------
// Helpers
//------------------------------------------------------------

inline bool is_exon(int s) {
    return (s == 1 || s == 2 || s == 3);
}

inline bool is_intron(int s) {
    return (s == 4 || s == 5 || s == 6);
}

inline bool is_label_exon(string s) {
    return (s == "E0" || s == "E1" || s == "E2");
}

inline bool is_label_intron(string s) {
    return (s == "I0" || s == "I1" || s == "I2");
}


//------------------------------------------------------------
// Read FASTA (single sequence)
//------------------------------------------------------------
string read_fasta(const string &file) {
    ifstream in(file);
    if (!in) {
        cerr << "Cannot open FASTA " << file << "\n";
        exit(1);
    }

    string seq, line;
    while (getline(in, line)) {
        if (!line.empty() && line[0] == '>') continue;
        for (char c : line) {
            if (!isspace(c)) seq.push_back(c);
        }
    }
    return seq;
}

//------------------------------------------------------------
// Load emissions: pos \t 7 values
//------------------------------------------------------------
vector<array<float, NUM_STATES>> load_emissions(const string &file, int L) {
    vector<array<float, NUM_STATES>> emit(L);
    for (int i = 0; i < L; i++)
        for (int s = 0; s < NUM_STATES; s++)
            emit[i][s] = NEG_INF;

    ifstream in(file);
    if (!in) {
        cerr << "Cannot open emissions " << file << "\n";
        exit(1);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        int pos;
        if (!(ss >> pos) || pos < 0 || pos >= L) {
            cerr << "Invalid emission position in " << file << ": " << line << "\n";
            exit(1);
        }
        // upstream: file now carries a single intron column; states 5 and 6
        // share the intron emission
        for (int s = 0; s < NUM_STATES-2; s++) {
            if (!(ss >> emit[pos][s]) || isnan(emit[pos][s])) {
                cerr << "Invalid emission row in " << file << ": " << line << "\n";
                exit(1);
            }
        }
        emit[pos][5]=emit[pos][4];
        emit[pos][6]=emit[pos][4];
    }
    return emit;
}

//------------------------------------------------------------
// Load sparse ATG/GT/AG scores
//------------------------------------------------------------
vector<double> load_sparse_scores(const string &file, int L) {
    vector<double> scores(L, NEG_INF);

    ifstream in(file);
    if (!in) {
        cerr << "Cannot open scores " << file << "\n";
        exit(1);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        int pos;
        double score;
        if (!(ss >> pos >> score) || pos < 0 || pos >= L || isnan(score)) {
            cerr << "Invalid sparse score row in " << file << ": " << line << "\n";
            exit(1);
        }
        scores[pos] = score;
    }
    return scores;
}

//------------------------------------------------------------
// Extract FASTA header
//------------------------------------------------------------
string get_fasta_header(const string &file) {
    ifstream in(file);
    if (!in) return "sequence";

    string line;
    while (getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            stringstream ss(line.substr(1));
            string id;
            ss >> id;
            return id;
        }
    }
    return "sequence";
}

//------------------------------------------------------------
// Safety check: verify GT/AG coordinates match the sequence
//------------------------------------------------------------
void safety_check_gt_ag_atg(
    const vector<char> &seq,
    const vector<double> &gt_score,
    const vector<double> &ag_score,
    const vector<double> &atg_score,
    const vector<double> &stop_score
) {
    int L = seq.size();
    for (int pos = 0; pos < L - 1; pos++) {

        // Check GT
        if (gt_score[pos] > NEG_INF) {
            string dinuc;
            dinuc.push_back(toupper(seq[pos]));
            dinuc.push_back(toupper(seq[pos + 1]));
            if (dinuc != "GT") {
                cerr << "WARNING: GT score "<< gt_score[pos] << " at position "
                     << pos << " but sequence has " << dinuc << "\n";
            }
        }

        // Check AG
        if (ag_score[pos] > NEG_INF) {
            string dinuc;
            dinuc.push_back(toupper(seq[pos]));
            dinuc.push_back(toupper(seq[pos + 1]));
            if (dinuc != "AG") {
                cerr << "WARNING: AG score "<< ag_score[pos] << " at position "
                     << pos << " but sequence has " << dinuc << "\n";
            }
        }

        // Check ATG
        if (atg_score[pos] > NEG_INF && pos + 2 < L) {
            string trinuc;
            trinuc.push_back(toupper(seq[pos]));
            trinuc.push_back(toupper(seq[pos + 1]));
            trinuc.push_back(toupper(seq[pos + 2]));
            if (trinuc != "ATG") {
                cerr << "WARNING: ATG score "<< atg_score[pos] << " at position "
                     << pos << " but sequence has " << trinuc << "\n";
            }
        }

        // Check STOP
        if (stop_score[pos] > NEG_INF && pos + 2 < L) {
            string trinuc;
            trinuc.push_back(toupper(seq[pos]));
            trinuc.push_back(toupper(seq[pos + 1]));
            trinuc.push_back(toupper(seq[pos + 2]));
            if (trinuc != "TAG"  && trinuc != "TGA" && trinuc != "TAA") {
                cerr << "WARNING: STOP score "<< stop_score[pos] << " at position "
                     << pos << " but sequence has " << trinuc << "\n";
            }
        }
    }
}

//------------------------------------------------------------
// Baseline transition matrix (log-space), excluding GT/AG overrides
//------------------------------------------------------------
vector<vector<double>> init_transitions() {
    vector<vector<double>> trans(NUM_STATES, vector<double>(NUM_STATES, NEG_INF));

    double self_prob = 0.0;

    // Noncoding self
    trans[0][0] = self_prob;

    // Exon frame cycling: none (self only)
    trans[1][1] = self_prob;
    trans[2][2] = self_prob;
    trans[3][3] = self_prob;

    // Intron states: self only
    trans[4][4] = self_prob;
    trans[5][5] = self_prob;
    trans[6][6] = self_prob;

    return trans;
}

//------------------------------------------------------------
// DP structures
//------------------------------------------------------------
struct DPCell {
    double dp;
    short int bt;
    int intron_len;
    int exon_len;
    int inter_len;
    short int exon_from;
};

vector<vector<DPCell>> init_dp(int L,
                               const vector<array<float, NUM_STATES>> &emit)
{
    vector<vector<DPCell>> dp(L, vector<DPCell>(NUM_STATES));

    // Initialization at position 0 — force start in N
    for (int s = 0; s < NUM_STATES; s++) {
        double start_prob = (s == 0 ? 0.0 : NEG_INF);
        double e = emit[0][s];

        dp[0][s].dp = start_prob + e;
        dp[0][s].bt = -1;
        dp[0][s].intron_len = is_intron(s) ? 1 : 0;
        dp[0][s].exon_len   = is_exon(s)   ? 1 : 0;
        dp[0][s].inter_len  = (s == 0)   ? 1 : 0;
        dp[0][s].exon_from   = is_exon(s)   ? 0 : 0;
    }

    return dp;
}

//------------------------------------------------------------
// Authoritative transition evaluation
//------------------------------------------------------------
struct ModelInputs {
    const vector<array<float, NUM_STATES>> &emit;
    const vector<double> &gt_score;
    const vector<double> &ag_score;
    const vector<double> &atg_score;
    const vector<double> &stop_score;
    const vector<char> &seq;
    const vector<vector<double>> &trans;
};

struct TransitionResult {
    bool allowed = false;
    double transition_score = NEG_INF;
    double emission_score = NEG_INF;
    PathMetadata next_metadata;
    string failure_reason;
};

static PathMetadata metadata_from_cell(const DPCell &cell) {
    return {cell.intron_len, cell.exon_len, cell.inter_len,
            cell.exon_from, cell.bt};
}

static int splice_motif_start_0based(int dp_position) {
    return dp_position - 1;
}

static int codon_start_0based(int dp_position) {
    return std::max(0, dp_position - 2);
}

static bool is_splice_transition(int from, int to) {
    return (is_exon(from) && is_intron(to)) ||
           (is_intron(from) && is_exon(to));
}

static bool sequence_has_dinucleotide(const vector<char> &seq, int start,
                                      char first, char second) {
    return start >= 0 && start + 1 < static_cast<int>(seq.size()) &&
           toupper(seq[start]) == first && toupper(seq[start + 1]) == second;
}

static string sequence_codon_ending_at(const vector<char> &seq, int position) {
    if (position < 2 || position >= static_cast<int>(seq.size())) return "";
    string codon;
    codon.push_back(toupper(seq[position - 2]));
    codon.push_back(toupper(seq[position - 1]));
    codon.push_back(toupper(seq[position]));
    return codon;
}

static bool is_stop_codon(const string &codon) {
    return codon == "TAA" || codon == "TAG" || codon == "TGA";
}

static double destination_emission(int position, int to,
                                   const ModelInputs &inputs) {
    double emission = inputs.emit[position][to];
    // Preserve the legacy TAG/AG overlap convention exactly.
    if (emission <= -1e6 && position > 0 && position + 1 <
            static_cast<int>(inputs.seq.size()) &&
        toupper(inputs.seq[position - 1]) == 'A' &&
        toupper(inputs.seq[position]) == 'G') {
        emission = inputs.emit[position + 1][to];
    }
    return emission;
}

static TransitionResult evaluate_transition(
    int position,
    int from,
    int to,
    const PathMetadata &previous,
    const ModelInputs &inputs
) {
    TransitionResult result;
    if (position <= 0 || position >= static_cast<int>(inputs.seq.size()) ||
        from < 0 || from >= NUM_STATES || to < 0 || to >= NUM_STATES) {
        result.failure_reason = "invalid_score_array_index";
        return result;
    }

    double log_t = inputs.trans[from][to];
    const string codon = sequence_codon_ending_at(inputs.seq, position);
    const bool stop = is_stop_codon(codon);

    if (is_exon(to) && is_exon(from) && from == to && stop &&
        codon_start_0based(position) % 3 == to - 1) {
        log_t = NEG_INF;
    }

    if (is_exon(from) && is_intron(to) &&
        sequence_has_dinucleotide(inputs.seq,
                                  splice_motif_start_0based(position), 'G', 'T')) {
        const int length = previous.exon_len - 2;
        log_t = NEG_INF;
        if ((to - 4) == (from - 1) &&
            (length >= MIN_EXON || previous.exon_from == 0)) {
            log_t = inputs.gt_score[splice_motif_start_0based(position)];
        } else {
            result.failure_reason = "minimum_exon_or_frame_violation";
        }
    }

    if (is_intron(from) && is_exon(to) &&
        sequence_has_dinucleotide(inputs.seq,
                                  splice_motif_start_0based(position), 'A', 'G')) {
        const int length = previous.intron_len + 2;
        log_t = NEG_INF;
        if (length >= MIN_INTRON) {
            const int intron_frame = from - 4;
            const int exon_frame = to - 1;
            if ((exon_frame - intron_frame + 3) % 3 == length % 3) {
                log_t = inputs.ag_score[splice_motif_start_0based(position)];
            } else {
                result.failure_reason = "frame_violation";
            }
        } else {
            result.failure_reason = "minimum_intron_violation";
        }
    }

    if (is_exon(from) && to == 0 && stop) {
        const int length = previous.exon_len - 2;
        const int frame = from - 1;
        if (codon_start_0based(position) % 3 == frame) {
            if (is_intron(previous.exon_from) ||
                (previous.exon_from == 0 && length > MIN_SINGLE)) {
                log_t = inputs.stop_score[codon_start_0based(position)];
            } else {
                // This finite legacy penalty is a model rule, not an invalid
                // transition sentinel, and is retained for compatibility.
                log_t = -1e3;
            }
        } else {
            result.failure_reason = "stop_codon_frame_violation";
        }
    }

    if (is_exon(to) && from == 0 && position >= 2) {
        const int length = previous.inter_len + 2;
        if (codon == "ATG" &&
            (length >= MIN_INTER || position < MIN_INTER) &&
            previous.predecessor_state == 0) {
            const int frame = to - 1;
            if (codon_start_0based(position) % 3 == frame) {
                log_t = position < 25
                    ? 1.0
                    : inputs.atg_score[codon_start_0based(position)];
            } else {
                result.failure_reason = "start_codon_frame_violation";
            }
        } else if (codon == "ATG") {
            result.failure_reason = "minimum_intergenic_violation";
        }
    }

    result.transition_score = log_t;
    result.emission_score = destination_emission(position, to, inputs);
    result.allowed = log_t > NEG_INF && result.emission_score > NEG_INF;
    if (!result.allowed && result.failure_reason.empty())
        result.failure_reason = "transition_or_emission_forbidden";

    PathMetadata next;
    next.predecessor_state = from;
    if (is_intron(to))
        next.intron_len = is_intron(from) ? previous.intron_len + 1 : 1;
    if (is_exon(to)) {
        if (is_exon(from)) {
            next.exon_len = previous.exon_len + 1;
            next.exon_from = previous.exon_from;
        } else {
            next.exon_len = from == 0 ? 3 : 1;
            next.exon_from = from;
        }
    }
    if (to == 0)
        next.inter_len = from == 0 ? previous.inter_len + 1 : 1;
    result.next_metadata = next;
    return result;
}

static void assign_metadata(DPCell &cell, const PathMetadata &metadata) {
    cell.intron_len = metadata.intron_len;
    cell.exon_len = metadata.exon_len;
    cell.inter_len = metadata.inter_len;
    cell.exon_from = metadata.exon_from;
}

//------------------------------------------------------------
// Full Viterbi DP
//------------------------------------------------------------
void run_viterbi(
    vector<vector<DPCell>> &dp,
    const vector<array<float, NUM_STATES>> &emit,
    const vector<double> &gt_score,
    const vector<double> &ag_score,
    const vector<double> &atg_score,
    const vector<double> &stop_score,
    const vector<char> &seq,
    const vector<vector<double>> &trans
) {
    const int L = seq.size();
    const ModelInputs inputs{emit, gt_score, ag_score, atg_score,
                             stop_score, seq, trans};

    for (int i = 1; i < L; i++) {
        for (int to = 0; to < NUM_STATES; to++) {
            double best = -1e18;
            int best_from = -1;
            PathMetadata best_metadata;

            for (int from = 0; from < NUM_STATES; from++) {
                const TransitionResult evaluated = evaluate_transition(
                    i, from, to, metadata_from_cell(dp[i - 1][from]), inputs);
                // The legacy DP deliberately propagates its finite NEG_INF
                // sentinel into unreachable cells; retaining that arithmetic is
                // required for byte-identical stderr and tie breaking.
                const double candidate = dp[i - 1][from].dp +
                    evaluated.transition_score + evaluated.emission_score;
                if (candidate > best) {
                    best = candidate;
                    best_from = from;
                    best_metadata = evaluated.next_metadata;
                }
            }

            dp[i][to].dp = best;
            dp[i][to].bt = best_from;
            assign_metadata(dp[i][to], best_metadata);
        }
    }
}

//------------------------------------------------------------
// Termination: find best final state
//------------------------------------------------------------
pair<double, int> viterbi_termination(
    const vector<vector<DPCell>> &dp,
    int L
) {
    double best_final = -1e18;
    int best_state = -1;

    for (int s = 0; s < NUM_STATES; s++) {
        if (dp[L - 1][s].dp > best_final) {
            best_final = dp[L - 1][s].dp;
            best_state = s;
        }
    }
    return {best_final, best_state};
}

//------------------------------------------------------------
// Backtrace: reconstruct optimal path
//------------------------------------------------------------
vector<int> viterbi_backtrace(
    const vector<vector<DPCell>> &dp,
    int L,
    int best_state
) {
    vector<int> path_states(L);
    int cur = best_state;

    for (int i = L - 1; i >= 0; i--) {
        path_states[i] = cur;
        cur = dp[i][cur].bt;
        if (cur < 0) break;
    }
    return path_states;
}

//------------------------------------------------------------
// Backtrace: get scores for optimal path
//------------------------------------------------------------
vector<double> viterbi_backtrace_scores(
    const vector<vector<DPCell>> &dp,
    int L,
    int best_state
) {
    vector<double> path_scores(L);
    int cur = best_state;

    for (int i = L - 1; i >= 0; i--) {
        path_scores[i] = dp[i][cur].dp;
        cur = dp[i][cur].bt;
        if (cur < 0) break;
    }
    return path_scores;
}

//------------------------------------------------------------
// Convert state numbers to labels
//------------------------------------------------------------
vector<string> states_to_labels(const vector<int> &path_states) {
    vector<string> labels(path_states.size());
    for (size_t i = 0; i < path_states.size(); i++) {
        labels[i] = state_name[path_states[i]];
    }
    return labels;
}

//------------------------------------------------------------
// Global index counter (mirrors Perl's $index)
//------------------------------------------------------------
int gff_index = 0;

//------------------------------------------------------------
// Write a single GFF3 feature
//------------------------------------------------------------
void write_gff_feature(
    const string &seqid,
    const string &state,
    int start0,
    int end0,
    const string &f_fasta,
    double best_final
) {
    // Extract filename from path
    size_t slash = f_fasta.find_last_of("/\\");
    string tname = (slash == string::npos)
        ? f_fasta
        : f_fasta.substr(slash + 1);

    int start = start0;   // 0-based
    int end   = end0;

    string type;

    if (state[0] == 'E') {
        type = "CDS";
        start = start0;
    }
    else if (state[0] == 'I') {
        type = "intron";
        end = end0;
    }
    else {
        type = "region";   // noncoding
        end = end0;
        gff_index++;
    }

    // Perl rounds best_final to 2 decimals
    double bf = floor(best_final * 100.0) / 100.0;

    if(type != "region")
      cout << seqid << "\t"
         << "UniAnn" << "\t"
         << type << "\t"
         << start << "\t"
         << end << "\t"
         << bf << "\t"
         << "+" << "\t"
         << "." << "\t"
         << "Parent=" << tname << "." << gff_index
         << ";state=" << state
         << "\n";
}

//------------------------------------------------------------
// Emit all GFF3 features from the state path
//------------------------------------------------------------
void write_gff_from_path(
    const vector<string> &labels,
    const vector<double> &scores,
    const string &seqid,
    const string &f_fasta,
    double best_final
) {
  cout << "##gff-version 3\n";
  string current_state = labels[0];
  string previous_state = labels[0];
  double score_offset = scores[0];

  int start = 0;
  int end = 0;

  for (size_t i = 1; i < labels.size(); i++) {
    if (labels[i] != current_state) {
      end = i - 1;
      if (current_state == "N")  { // transitioned to exon
        end = i - 2;
      }
      if (labels[i] == "N") { //transitioned to non-coding
        end = i + 1;
      }
      if(end > start) {
        if(is_label_exon(current_state) && is_label_intron(previous_state)) {
          write_gff_feature(seqid, current_state, start+2, end, f_fasta, (scores[i]-score_offset)/(end-start+1));
        }else if(is_label_intron(current_state) && is_label_exon(previous_state)) {
          write_gff_feature(seqid, current_state, start, end+2, f_fasta, (scores[i]-score_offset)/(end-start+1));
        }else{
          write_gff_feature(seqid, current_state, start, end, f_fasta, (scores[i]-score_offset)/(end-start+1));
        }
      }
      if (current_state == "N"){
        start = i - 1;
        score_offset = scores[i];
      } else if (labels[i] == "N" ){
        start = i + 2;
      } else {
        start = i;
      }
      previous_state = current_state;
      current_state = labels[i];
    }
  }

  // Final segment
  write_gff_feature(seqid, current_state, start, labels.size() - 1, f_fasta, best_final);
}


//------------------------------------------------------------
// Transcript conversion and splice-anchor extraction
//------------------------------------------------------------
static string anchor_type_name(AnchorType type) {
    return type == AnchorType::Donor ? "donor" : "acceptor";
}

static string candidate_status_name(CandidateStatus status) {
    switch (status) {
        case CandidateStatus::Complete: return "complete";
        case CandidateStatus::NoValidUpstreamBoundary:
            return "no_valid_upstream_boundary";
        case CandidateStatus::NoValidDownstreamBoundary:
            return "no_valid_downstream_boundary";
        case CandidateStatus::MaximumDistanceReached:
            return "maximum_distance_reached";
        case CandidateStatus::ScoreValidationFailure:
            return "score_validation_failure";
        case CandidateStatus::AnnotationConversionFailure:
            return "annotation_conversion_failure";
    }
    return "unknown";
}

static vector<TranscriptAnnotation> build_transcript_annotations(
    const vector<int> &positions,
    const vector<int> &states,
    const string &seqid
) {
    vector<TranscriptAnnotation> transcripts;
    if (positions.size() != states.size() || states.size() < 2) return transcripts;

    bool inside = false;
    int exon_start = -1;
    int pending_donor = -1;
    TranscriptAnnotation current;
    current.sequence_id = seqid;

    for (size_t k = 1; k < states.size(); ++k) {
        const int from = states[k - 1];
        const int to = states[k];
        const int dp_position = positions[k];

        if (!inside && from == 0 && is_exon(to)) {
            inside = true;
            current = TranscriptAnnotation{};
            current.sequence_id = seqid;
            current.strand = '+';
            current.path_start = dp_position;
            current.genomic_start_0based = codon_start_0based(dp_position);
            exon_start = current.genomic_start_0based;
            pending_donor = -1;
        }
        if (!inside) continue;

        if (is_exon(from) && is_intron(to)) {
            const int donor = splice_motif_start_0based(dp_position);
            if (exon_start < 0 || donor <= exon_start) {
                inside = false;
                continue;
            }
            current.exons.push_back({exon_start, donor});
            current.cds_intervals.push_back({exon_start, donor});
            pending_donor = donor;
        } else if (is_intron(from) && is_exon(to)) {
            const int acceptor = splice_motif_start_0based(dp_position);
            if (pending_donor < 0 || acceptor < pending_donor) {
                inside = false;
                continue;
            }
            current.junctions.push_back({pending_donor, acceptor, '+'});
            exon_start = acceptor + 2;
            pending_donor = -1;
        } else if (is_exon(from) && to == 0) {
            const int transcript_end = dp_position + 1;
            if (exon_start < 0 || transcript_end <= exon_start ||
                pending_donor >= 0) {
                inside = false;
                continue;
            }
            current.exons.push_back({exon_start, transcript_end});
            current.cds_intervals.push_back({exon_start, transcript_end});
            current.path_end = dp_position;
            current.genomic_end_0based_exclusive = transcript_end;
            transcripts.push_back(current);
            inside = false;
            exon_start = -1;
        }
    }
    return transcripts;
}

static vector<TranscriptAnnotation> build_transcript_annotations(
    const vector<int> &states,
    const string &seqid
) {
    vector<int> positions(states.size());
    for (size_t i = 0; i < positions.size(); ++i) positions[i] = i;
    return build_transcript_annotations(positions, states, seqid);
}

static int transcript_for_position(const vector<TranscriptAnnotation> &transcripts,
                                   int dp_position) {
    for (size_t i = 0; i < transcripts.size(); ++i) {
        if (transcripts[i].path_start <= dp_position &&
            dp_position <= transcripts[i].path_end)
            return static_cast<int>(i);
    }
    return -1;
}

static vector<SpliceAnchor> extract_splice_anchors(
    const vector<int> &path_states,
    const vector<TranscriptAnnotation> &transcripts,
    const string &seqid,
    const ModelInputs &inputs,
    bool debug
) {
    vector<SpliceAnchor> anchors;
    for (int i = 1; i < static_cast<int>(path_states.size()); ++i) {
        const int from = path_states[i - 1];
        const int to = path_states[i];
        if (!is_splice_transition(from, to)) continue;
        const int transcript_index = transcript_for_position(transcripts, i);
        if (transcript_index < 0) continue;

        SpliceAnchor anchor;
        anchor.sequence_id = seqid;
        anchor.path_index = i;
        anchor.dp_position = i;
        anchor.genomic_position_0based = splice_motif_start_0based(i);
        anchor.genomic_position_1based = anchor.genomic_position_0based + 1;
        anchor.type = is_exon(from) ? AnchorType::Donor : AnchorType::Acceptor;
        anchor.left_state = from;
        anchor.right_state = to;
        anchor.frame = anchor.type == AnchorType::Donor ? from - 1 : to - 1;
        anchor.splice_score = anchor.type == AnchorType::Donor
            ? inputs.gt_score[anchor.genomic_position_0based]
            : inputs.ag_score[anchor.genomic_position_0based];
        anchor.reference_transcript_index = transcript_index;
        anchors.push_back(anchor);

        if (debug) {
            cerr << "ALT anchor " << anchor_type_name(anchor.type)
                 << " seq=" << seqid
                 << " motif0=" << anchor.genomic_position_0based
                 << " dp=" << anchor.dp_position
                 << " transition=" << state_name[from] << "->"
                 << state_name[to] << " frame=" << anchor.frame
                 << " score=" << anchor.splice_score << "\n";
        }
    }
    return anchors;
}

static string intron_chain_key(const TranscriptAnnotation &transcript) {
    ostringstream out;
    out << transcript.sequence_id << '|' << transcript.strand << '|';
    for (const auto &junction : transcript.junctions)
        out << junction.donor_0based << '-' << junction.acceptor_0based << ',';
    return out.str();
}

static string full_transcript_key(const TranscriptAnnotation &transcript) {
    ostringstream out;
    out << intron_chain_key(transcript) << '|' << transcript.genomic_start_0based
        << '-' << transcript.genomic_end_0based_exclusive << '|';
    for (const auto &exon : transcript.exons)
        out << exon.start_0based << '-' << exon.end_0based_exclusive << ',';
    out << '|';
    for (const auto &cds : transcript.cds_intervals)
        out << cds.start_0based << '-' << cds.end_0based_exclusive << ',';
    return out.str();
}

//------------------------------------------------------------
// Explicit path validation and scoring
//------------------------------------------------------------
struct ExplicitScoreResult {
    bool valid = false;
    double score = NEG_INF;
    int anchor_count = 0;
    string failure_reason;
};

static ExplicitScoreResult score_path_explicitly(
    const vector<int> &positions,
    const vector<int> &states,
    const PathMetadata &initial_metadata,
    const ModelInputs &inputs,
    const SpliceAnchor *required_anchor,
    bool require_complete
) {
    ExplicitScoreResult result;
    if (positions.size() != states.size() || positions.size() < 2) {
        result.failure_reason = "non_contiguous_path";
        return result;
    }
    PathMetadata metadata = initial_metadata;
    double score = 0.0;
    bool transcript_started = false;
    bool transcript_ended = false;

    for (size_t k = 1; k < states.size(); ++k) {
        if (positions[k] != positions[k - 1] + 1) {
            result.failure_reason = "non_contiguous_path";
            return result;
        }
        const int from = states[k - 1];
        const int to = states[k];
        const TransitionResult transition = evaluate_transition(
            positions[k], from, to, metadata, inputs);

        if (required_anchor && positions[k] == required_anchor->dp_position &&
            from == required_anchor->left_state &&
            to == required_anchor->right_state) {
            ++result.anchor_count;
        }

        if (!transcript_started && from == 0 && to == 0) {
            metadata = transition.next_metadata;
            continue;  // virtual intergenic source is outside the scored model
        }
        if (!transition.allowed) {
            result.failure_reason = transition.failure_reason;
            return result;
        }
        if (transcript_ended && from == 0 && to == 0) {
            score += transition.transition_score + transition.emission_score;
            metadata = transition.next_metadata;
            continue;
        }
        if (!transcript_started) {
            if (!(from == 0 && is_exon(to))) {
                result.failure_reason = "no_valid_upstream_boundary";
                return result;
            }
            transcript_started = true;
        }
        if (transcript_ended) {
            result.failure_reason = "multiple_transcripts_in_candidate";
            return result;
        }
        score += transition.transition_score + transition.emission_score;
        metadata = transition.next_metadata;
        if (is_exon(from) && to == 0) transcript_ended = true;
    }

    if (required_anchor && result.anchor_count != 1) {
        result.failure_reason = "anchor_not_present_exactly_once";
        return result;
    }
    if (require_complete && (!transcript_started || !transcript_ended)) {
        result.failure_reason = "no_valid_downstream_boundary";
        return result;
    }
    result.valid = true;
    result.score = score;
    return result;
}

static ExplicitScoreResult score_interval_explicitly(
    const vector<int> &positions,
    const vector<int> &states,
    const PathMetadata &initial_metadata,
    const ModelInputs &inputs
) {
    ExplicitScoreResult result;
    if (positions.size() != states.size() || positions.size() < 2) {
        result.failure_reason = "non_contiguous_path";
        return result;
    }
    PathMetadata metadata = initial_metadata;
    double score = 0.0;
    for (size_t k = 1; k < states.size(); ++k) {
        if (positions[k] != positions[k - 1] + 1) {
            result.failure_reason = "non_contiguous_path";
            return result;
        }
        const TransitionResult transition = evaluate_transition(
            positions[k], states[k - 1], states[k], metadata, inputs);
        if (!transition.allowed) {
            result.failure_reason = transition.failure_reason;
            return result;
        }
        score += transition.transition_score + transition.emission_score;
        metadata = transition.next_metadata;
    }
    result.valid = true;
    result.score = score;
    return result;
}

static bool valid_transcript_geometry(const TranscriptAnnotation &transcript,
                                      string &reason) {
    if (transcript.exons.empty() ||
        transcript.exons.size() != transcript.cds_intervals.size() ||
        transcript.junctions.size() + 1 != transcript.exons.size()) {
        reason = "annotation_conversion_failure";
        return false;
    }
    if (transcript.junctions.empty() &&
        transcript.exons.front().end_0based_exclusive -
            transcript.exons.front().start_0based - 3 <= MIN_SINGLE) {
        reason = "single_exon_minimum_length_violation";
        return false;
    }
    for (const auto &junction : transcript.junctions) {
        if (junction.acceptor_0based + 2 - junction.donor_0based < MIN_INTRON) {
            reason = "minimum_intron_violation";
            return false;
        }
    }
    int previous_end = -1;
    for (const auto &exon : transcript.exons) {
        if (exon.start_0based < 0 || exon.end_0based_exclusive <= exon.start_0based ||
            exon.start_0based < previous_end) {
            reason = "impossible_overlapping_segments";
            return false;
        }
        previous_end = exon.end_0based_exclusive;
    }
    for (size_t i = 0; i < transcript.junctions.size(); ++i) {
        const auto &junction = transcript.junctions[i];
        if (junction.donor_0based != transcript.exons[i].end_0based_exclusive ||
            junction.acceptor_0based + 2 != transcript.exons[i + 1].start_0based) {
            reason = "ambiguous_splice_coordinates";
            return false;
        }
    }
    return true;
}

//------------------------------------------------------------
// Bounded anchored constrained-forward search
//------------------------------------------------------------
struct LocalCell {
    bool valid = false;
    double score = NEG_INF;
    PathMetadata metadata;
    int bt_state = -1;
    int bt_flag = -1;
    double score_at_anchor = NEG_INF;
    double anchor_contribution = NEG_INF;
};

static bool differs_from_reference_splice(
    int position, int from, int to, const vector<int> &reference_states
) {
    const bool candidate_splice = is_splice_transition(from, to);
    const bool reference_splice = position > 0 &&
        is_splice_transition(reference_states[position - 1],
                             reference_states[position]);
    if (!candidate_splice && !reference_splice) return false;
    return from != reference_states[position - 1] ||
           to != reference_states[position];
}

static AlternativeCandidate search_anchor(
    const SpliceAnchor &anchor,
    const vector<int> &reference_states,
    const vector<vector<DPCell>> &global_dp,
    const vector<TranscriptAnnotation> &references,
    const ModelInputs &inputs,
    const AlternativeOptions &options
) {
    AlternativeCandidate candidate;
    candidate.anchor = anchor;
    const TranscriptAnnotation &reference =
        references[anchor.reference_transcript_index];
    const int L = inputs.seq.size();
    int locus_start = 0;
    int locus_end = L - 1;
    if (anchor.reference_transcript_index > 0)
        locus_start = references[anchor.reference_transcript_index - 1]
                          .genomic_end_0based_exclusive;
    if (anchor.reference_transcript_index + 1 < static_cast<int>(references.size()))
        locus_end = references[anchor.reference_transcript_index + 1]
                        .genomic_start_0based - 1;

    const int distance_start = max(0, anchor.dp_position - options.max_distance);
    const int distance_end = min(L - 1, anchor.dp_position + options.max_distance);
    const int start = max(locus_start, distance_start);
    const int end = min(locus_end, distance_end);
    candidate.local_interval_start = start;
    candidate.local_interval_end = end;

    if (start >= anchor.dp_position || reference_states[start] != 0) {
        candidate.status = distance_start > locus_start
            ? CandidateStatus::MaximumDistanceReached
            : CandidateStatus::NoValidUpstreamBoundary;
        candidate.failure_reason = candidate_status_name(candidate.status);
        return candidate;
    }

    using FlagCells = array<LocalCell, 4>; // bit 1: anchor; bit 0: splice-diverged
    vector<array<FlagCells, NUM_STATES>> local(end - start + 1);
    LocalCell &source = local[0][0][0];
    source.valid = true;
    source.score = 0.0;
    source.metadata = metadata_from_cell(global_dp[start][0]);
    if (start > 0 && source.metadata.predecessor_state < 0)
        source.metadata.predecessor_state = 0;

    double best_endpoint_score[2] = {NEG_INF, NEG_INF};
    int best_endpoint_offset[2] = {-1, -1};
    int best_endpoint_flag[2] = {-1, -1};

    for (int position = start + 1; position <= end; ++position) {
        const int offset = position - start;
        for (int to = 0; to < NUM_STATES; ++to) {
            for (int from = 0; from < NUM_STATES; ++from) {
                for (int old_flag = 0; old_flag < 4; ++old_flag) {
                    const LocalCell &previous = local[offset - 1][from][old_flag];
                    if (!previous.valid) continue;
                    const bool seen = (old_flag & 2) != 0;
                    if (seen && from == 0) continue; // transcript already terminated

                    const TransitionResult transition = evaluate_transition(
                        position, from, to, previous.metadata, inputs);
                    if (!transition.allowed) continue;

                    const bool exact_anchor = position == anchor.dp_position &&
                        from == anchor.left_state && to == anchor.right_state;
                    const bool new_seen = seen || exact_anchor;
                    if (is_exon(from) && to == 0 && !new_seen) continue;

                    const bool new_diverged = (old_flag & 1) != 0 ||
                        differs_from_reference_splice(position, from, to,
                                                      reference_states);
                    const int new_flag = (new_seen ? 2 : 0) |
                                         (new_diverged ? 1 : 0);
                    double score;
                    if (!seen && from == 0 && to == 0) {
                        // All legal intergenic starting boundaries are represented
                        // without rewarding or penalizing unmodeled flanking bases.
                        score = 0.0;
                    } else {
                        score = previous.score + transition.transition_score +
                                transition.emission_score;
                    }

                    LocalCell &destination = local[offset][to][new_flag];
                    if (!destination.valid || score > destination.score) {
                        destination.valid = true;
                        destination.score = score;
                        destination.metadata = transition.next_metadata;
                        destination.bt_state = from;
                        destination.bt_flag = old_flag;
                        destination.score_at_anchor = exact_anchor
                            ? score : previous.score_at_anchor;
                        destination.anchor_contribution = exact_anchor
                            ? transition.transition_score + transition.emission_score
                            : previous.anchor_contribution;
                    }
                }
            }
        }

        for (int divergent = 0; divergent <= 1; ++divergent) {
            const int flag = 2 | divergent;
            const LocalCell &endpoint = local[offset][0][flag];
            if (endpoint.valid && endpoint.score > best_endpoint_score[divergent]) {
                best_endpoint_score[divergent] = endpoint.score;
                best_endpoint_offset[divergent] = offset;
                best_endpoint_flag[divergent] = flag;
            }
        }
    }

    // Prefer the explicitly splice-divergent hypothesis.  Retain the reference
    // reconstruction in diagnostics when no divergent complete path exists.
    int chosen_divergence = best_endpoint_offset[1] >= 0 ? 1 : 0;
    if (best_endpoint_offset[chosen_divergence] < 0) {
        candidate.status = distance_end < locus_end
            ? CandidateStatus::MaximumDistanceReached
            : CandidateStatus::NoValidDownstreamBoundary;
        candidate.upstream_complete = true;
        candidate.failure_reason = candidate_status_name(candidate.status);
        return candidate;
    }

    const int endpoint_offset = best_endpoint_offset[chosen_divergence];
    int state = 0;
    int flag = best_endpoint_flag[chosen_divergence];
    vector<int> reversed_positions;
    vector<int> reversed_states;
    for (int offset = endpoint_offset; offset >= 0; --offset) {
        reversed_positions.push_back(start + offset);
        reversed_states.push_back(state);
        if (offset == 0) break;
        const LocalCell &cell = local[offset][state][flag];
        state = cell.bt_state;
        flag = cell.bt_flag;
        if (state < 0 || flag < 0) {
            candidate.failure_reason = "impossible_backpointer";
            candidate.status = CandidateStatus::ScoreValidationFailure;
            return candidate;
        }
    }
    reverse(reversed_positions.begin(), reversed_positions.end());
    reverse(reversed_states.begin(), reversed_states.end());
    candidate.path_positions = move(reversed_positions);
    candidate.path_states = move(reversed_states);

    const LocalCell &endpoint = local[endpoint_offset][0]
                                    [best_endpoint_flag[chosen_divergence]];
    candidate.total_score = endpoint.score;
    candidate.upstream_score = endpoint.score_at_anchor;
    candidate.anchor_score = endpoint.anchor_contribution;
    candidate.downstream_score = endpoint.score - endpoint.score_at_anchor;
    candidate.upstream_complete = true;
    candidate.downstream_complete = true;

    const ExplicitScoreResult explicit_score = score_path_explicitly(
        candidate.path_positions, candidate.path_states,
        metadata_from_cell(global_dp[start][0]), inputs, &anchor, true);
    const double scale = max({1.0, fabs(candidate.total_score),
                              fabs(explicit_score.score)});
    if (!explicit_score.valid ||
        fabs(candidate.total_score - explicit_score.score) >
            SCORE_TOLERANCE * scale) {
        candidate.status = CandidateStatus::ScoreValidationFailure;
        candidate.failure_reason = explicit_score.valid
            ? "dp_explicit_score_disagreement" : explicit_score.failure_reason;
        return candidate;
    }

    auto annotations = build_transcript_annotations(
        candidate.path_positions, candidate.path_states, anchor.sequence_id);
    if (annotations.size() != 1) {
        candidate.status = CandidateStatus::AnnotationConversionFailure;
        candidate.failure_reason = "annotation_conversion_failure";
        return candidate;
    }
    candidate.transcript = annotations.front();
    candidate.transcript.score = candidate.total_score;
    if (!valid_transcript_geometry(candidate.transcript,
                                   candidate.failure_reason)) {
        candidate.status = CandidateStatus::AnnotationConversionFailure;
        return candidate;
    }

    const int comparison_start = candidate.transcript.path_start - 1;
    const int comparison_end = candidate.transcript.path_end;
    vector<int> reference_positions;
    vector<int> reference_slice;
    for (int position = comparison_start; position <= comparison_end; ++position) {
        reference_positions.push_back(position);
        reference_slice.push_back(reference_states[position]);
    }
    const int initial_reference_state = reference_states[comparison_start];
    const ExplicitScoreResult reference_score = score_interval_explicitly(
        reference_positions, reference_slice,
        metadata_from_cell(global_dp[comparison_start][initial_reference_state]),
        inputs);
    if (!reference_score.valid) {
        candidate.status = CandidateStatus::ScoreValidationFailure;
        candidate.failure_reason = "reference_interval_" +
                                   reference_score.failure_reason;
        return candidate;
    }
    candidate.reference_interval_score = reference_score.score;
    candidate.score_delta = candidate.total_score - reference_score.score;

    candidate.identical_to_global =
        intron_chain_key(candidate.transcript) == intron_chain_key(reference);
    candidate.valid = true;
    candidate.status = CandidateStatus::Complete;
    return candidate;
}

//------------------------------------------------------------
// Classification and deterministic deduplication
//------------------------------------------------------------
static bool interval_contains(const GenomicInterval &outer,
                              const GenomicInterval &inner) {
    return outer.start_0based <= inner.start_0based &&
           inner.end_0based_exclusive <= outer.end_0based_exclusive;
}

static string classify_splice_difference(
    const TranscriptAnnotation &reference,
    const TranscriptAnnotation &candidate
) {
    if (reference.junctions == candidate.junctions) return "identical";
    if (reference.junctions.size() == candidate.junctions.size()) {
        int donor_changes = 0;
        int acceptor_changes = 0;
        for (size_t i = 0; i < reference.junctions.size(); ++i) {
            donor_changes += reference.junctions[i].donor_0based !=
                             candidate.junctions[i].donor_0based;
            acceptor_changes += reference.junctions[i].acceptor_0based !=
                                candidate.junctions[i].acceptor_0based;
        }
        if (donor_changes == 1 && acceptor_changes == 0)
            return "alternative_donor";
        if (donor_changes == 0 && acceptor_changes == 1)
            return "alternative_acceptor";
    }

    for (size_t i = 1; i + 1 < reference.exons.size(); ++i) {
        bool covered = false;
        for (const auto &exon : candidate.exons)
            covered = covered || interval_contains(exon, reference.exons[i]);
        if (!covered && candidate.exons.size() < reference.exons.size())
            return "exon_skipping";
    }
    for (size_t i = 1; i + 1 < candidate.exons.size(); ++i) {
        bool covered = false;
        for (const auto &exon : reference.exons)
            covered = covered || interval_contains(exon, candidate.exons[i]);
        if (!covered && candidate.exons.size() > reference.exons.size())
            return "alternative_exon";
    }
    for (const auto &junction : reference.junctions) {
        GenomicInterval intron{junction.donor_0based,
                               junction.acceptor_0based + 2};
        for (const auto &exon : candidate.exons)
            if (interval_contains(exon, intron)) return "intron_retention";
    }
    return "complex";
}

static vector<int> deduplicate_candidates(
    vector<AlternativeCandidate> &candidates,
    const vector<TranscriptAnnotation> &references,
    double min_score_delta
) {
    vector<int> order;
    for (size_t i = 0; i < candidates.size(); ++i) {
        auto &candidate = candidates[i];
        if (!candidate.valid) continue;
        const auto &reference =
            references[candidate.anchor.reference_transcript_index];
        candidate.classification = classify_splice_difference(
            reference, candidate.transcript);
        candidate.identical_to_global = candidate.classification == "identical";
        if (!candidate.identical_to_global &&
            candidate.score_delta >= min_score_delta)
            order.push_back(i);
    }
    sort(order.begin(), order.end(), [&](int left, int right) {
        const auto &a = candidates[left];
        const auto &b = candidates[right];
        if (a.total_score != b.total_score) return a.total_score > b.total_score;
        const string ak = intron_chain_key(a.transcript);
        const string bk = intron_chain_key(b.transcript);
        if (ak != bk) return ak < bk;
        if (a.transcript.genomic_start_0based != b.transcript.genomic_start_0based)
            return a.transcript.genomic_start_0based <
                   b.transcript.genomic_start_0based;
        if (a.transcript.genomic_end_0based_exclusive !=
            b.transcript.genomic_end_0based_exclusive)
            return a.transcript.genomic_end_0based_exclusive <
                   b.transcript.genomic_end_0based_exclusive;
        return a.anchor.genomic_position_0based < b.anchor.genomic_position_0based;
    });

    set<pair<int, string>> seen;
    vector<int> kept;
    for (int index : order) {
        const auto key = make_pair(
            candidates[index].anchor.reference_transcript_index,
            intron_chain_key(candidates[index].transcript));
        if (seen.insert(key).second) {
            candidates[index].kept_after_deduplication = true;
            kept.push_back(index);
        }
    }
    return kept;
}

//------------------------------------------------------------
// Alternative GFF3 and diagnostic report
//------------------------------------------------------------
static string format_score(double score) {
    if (!isfinite(score) || score <= NEG_INF) return ".";
    ostringstream out;
    out << setprecision(12) << score;
    return out.str();
}

static void write_hierarchical_transcript(
    ostream &out,
    const TranscriptAnnotation &transcript,
    const string &gene_id,
    const string &transcript_id
) {
    out << transcript.sequence_id << "\tUniAnn\ttranscript\t"
        << transcript.genomic_start_0based + 1 << '\t'
        << transcript.genomic_end_0based_exclusive << '\t'
        << format_score(transcript.score) << "\t+\t.\tID=" << transcript_id
        << ";Parent=" << gene_id << "\n";

    int cumulative_cds = 0;
    for (size_t i = 0; i < transcript.exons.size(); ++i) {
        const auto &exon = transcript.exons[i];
        out << transcript.sequence_id << "\tUniAnn\texon\t"
            << exon.start_0based + 1 << '\t' << exon.end_0based_exclusive
            << "\t.\t+\t.\tID=" << transcript_id << ".exon" << i + 1
            << ";Parent=" << transcript_id << "\n";
        const auto &cds = transcript.cds_intervals[i];
        const int phase = i == 0 ? 0 : (3 - cumulative_cds % 3) % 3;
        out << transcript.sequence_id << "\tUniAnn\tCDS\t"
            << cds.start_0based + 1 << '\t' << cds.end_0based_exclusive
            << "\t.\t+\t" << phase << "\tID=" << transcript_id
            << ".cds" << i + 1 << ";Parent=" << transcript_id << "\n";
        cumulative_cds += cds.end_0based_exclusive - cds.start_0based;
        if (i < transcript.junctions.size()) {
            const auto &junction = transcript.junctions[i];
            out << transcript.sequence_id << "\tUniAnn\tintron\t"
                << junction.donor_0based + 1 << '\t'
                << junction.acceptor_0based + 2
                << "\t.\t+\t.\tID=" << transcript_id << ".intron"
                << i + 1 << ";Parent=" << transcript_id << "\n";
        }
    }
}

static bool write_alternative_gff3(
    const string &filename,
    vector<TranscriptAnnotation> references,
    const vector<AlternativeCandidate> &candidates,
    const vector<int> &kept
) {
    ofstream out(filename);
    if (!out) return false;
    out << "##gff-version 3\n";
    for (size_t locus = 0; locus < references.size(); ++locus) {
        vector<const AlternativeCandidate *> alternatives;
        int gene_start = references[locus].genomic_start_0based;
        int gene_end = references[locus].genomic_end_0based_exclusive;
        for (int index : kept) {
            if (candidates[index].anchor.reference_transcript_index !=
                static_cast<int>(locus)) continue;
            alternatives.push_back(&candidates[index]);
            gene_start = min(gene_start,
                candidates[index].transcript.genomic_start_0based);
            gene_end = max(gene_end,
                candidates[index].transcript.genomic_end_0based_exclusive);
        }
        ostringstream gene_builder;
        gene_builder << "UniAnnGene" << setw(6) << setfill('0') << locus + 1;
        const string gene_id = gene_builder.str();
        out << references[locus].sequence_id << "\tUniAnn\tgene\t"
            << gene_start + 1 << '\t' << gene_end
            << "\t.\t+\t.\tID=" << gene_id << "\n";
        references[locus].score = references[locus].score <= NEG_INF
            ? 0.0 : references[locus].score;
        write_hierarchical_transcript(out, references[locus], gene_id,
                                      gene_id + ".t1");
        for (size_t i = 0; i < alternatives.size(); ++i)
            write_hierarchical_transcript(out, alternatives[i]->transcript,
                gene_id, gene_id + ".t" + to_string(i + 2));
    }
    return static_cast<bool>(out);
}

static string tsv_escape_chain(const TranscriptAnnotation &transcript) {
    ostringstream out;
    for (size_t i = 0; i < transcript.junctions.size(); ++i) {
        if (i) out << ',';
        out << transcript.junctions[i].donor_0based << '-'
            << transcript.junctions[i].acceptor_0based;
    }
    return out.str();
}

static bool write_alt_report(const string &filename,
                             const vector<AlternativeCandidate> &candidates,
                             bool include_incomplete) {
    ofstream out(filename);
    if (!out) return false;
    out << "sequence_id\toriginal_gene_id\toriginal_transcript_id\t"
        << "candidate_id\tanchor_position_0based\tanchor_position_1based\t"
        << "anchor_type\tanchor_left_state\tanchor_right_state\t"
        << "candidate_start_0based\tcandidate_end_0based_exclusive\t"
        << "candidate_start_1based\tcandidate_end_1based_inclusive\tstrand\t"
        << "upstream_complete\tdownstream_complete\tstatus\tfailure_reason\t"
        << "candidate_total_score\treference_interval_score\tscore_delta\t"
        << "num_exons\tnum_introns\tintron_chain\tclassification\t"
        << "identical_to_global\tdeduplication_key\tkept_after_deduplication\n";
    int candidate_id = 0;
    for (const auto &candidate : candidates) {
        if (!candidate.valid && !include_incomplete) continue;
        ++candidate_id;
        ostringstream gene;
        gene << "UniAnnGene" << setw(6) << setfill('0')
             << candidate.anchor.reference_transcript_index + 1;
        const auto &t = candidate.transcript;
        out << candidate.anchor.sequence_id << '\t' << gene.str() << '\t'
            << gene.str() << ".t1\tCandidate" << setw(6) << setfill('0')
            << candidate_id << setfill(' ') << '\t'
            << candidate.anchor.genomic_position_0based << '\t'
            << candidate.anchor.genomic_position_1based << '\t'
            << anchor_type_name(candidate.anchor.type) << '\t'
            << state_name[candidate.anchor.left_state] << '\t'
            << state_name[candidate.anchor.right_state] << '\t'
            << t.genomic_start_0based << '\t'
            << t.genomic_end_0based_exclusive << '\t'
            << (t.genomic_start_0based >= 0 ? t.genomic_start_0based + 1 : -1)
            << '\t' << t.genomic_end_0based_exclusive << "\t+\t"
            << candidate.upstream_complete << '\t'
            << candidate.downstream_complete << '\t'
            << candidate_status_name(candidate.status) << '\t'
            << candidate.failure_reason << '\t'
            << format_score(candidate.total_score) << '\t'
            << format_score(candidate.reference_interval_score) << '\t'
            << format_score(candidate.score_delta) << '\t'
            << t.exons.size() << '\t' << t.junctions.size() << '\t'
            << tsv_escape_chain(t) << '\t' << candidate.classification << '\t'
            << candidate.identical_to_global << '\t'
            << (candidate.valid ? intron_chain_key(t) : "") << '\t'
            << candidate.kept_after_deduplication << '\n';
    }
    return static_cast<bool>(out);
}

struct AlternativeSummary {
    int donor_anchors = 0;
    int acceptor_anchors = 0;
    int complete = 0;
    int incomplete = 0;
    int invalid = 0;
    int identical = 0;
    map<string, int> classifications;
    set<string> full_keys;
    set<string> chain_keys;
};

static void print_alternative_summary(
    size_t original_count,
    const vector<SpliceAnchor> &anchors,
    const vector<AlternativeCandidate> &candidates,
    const vector<int> &kept
) {
    AlternativeSummary summary;
    for (const auto &anchor : anchors) {
        if (anchor.type == AnchorType::Donor) ++summary.donor_anchors;
        else ++summary.acceptor_anchors;
    }
    for (const auto &candidate : candidates) {
        if (candidate.valid) {
            ++summary.complete;
            summary.full_keys.insert(full_transcript_key(candidate.transcript));
            summary.chain_keys.insert(intron_chain_key(candidate.transcript));
            if (candidate.identical_to_global) ++summary.identical;
        } else if (candidate.upstream_complete ||
                   candidate.status == CandidateStatus::MaximumDistanceReached) {
            ++summary.incomplete;
        } else {
            ++summary.invalid;
        }
    }
    for (int index : kept) ++summary.classifications[candidates[index].classification];
    cerr << "Number of original transcripts: " << original_count << '\n'
         << "Number of splice-site anchors: " << anchors.size() << '\n'
         << "Number of donor anchors: " << summary.donor_anchors << '\n'
         << "Number of acceptor anchors: " << summary.acceptor_anchors << '\n'
         << "Number of complete anchored candidates: " << summary.complete << '\n'
         << "Number of incomplete anchored candidates: " << summary.incomplete << '\n'
         << "Number of invalid candidates: " << summary.invalid << '\n'
         << "Number of candidates before deduplication: " << candidates.size() << '\n'
         << "Number of unique full transcript structures: "
         << summary.full_keys.size() << '\n'
         << "Number of unique intron chains: " << summary.chain_keys.size() << '\n'
         << "Number identical to original transcripts: " << summary.identical << '\n'
         << "Number of alternative-splicing candidates: " << kept.size() << '\n'
         << "Number of alternative donors: "
         << summary.classifications["alternative_donor"] << '\n'
         << "Number of alternative acceptors: "
         << summary.classifications["alternative_acceptor"] << '\n'
         << "Number of exon-skipping candidates: "
         << summary.classifications["exon_skipping"] << '\n'
         << "Number of intron-retention candidates: "
         << summary.classifications["intron_retention"] << '\n'
         << "Number of complex candidates: "
         << summary.classifications["complex"] << '\n';
}

static bool parse_local_k_best_options(int argc, char **argv,
                                      LocalKBestOptions &options,
                                      string &error) {
    bool has_drop = false;
    for (int i = 7; i < argc; ++i) {
        const string argument = argv[i];
        if (argument == "--local-k-best") {
            options.enabled = true;
            if (i + 1 < argc && argv[i+1][0] != '-') {
                try {
                    options.k = stoi(argv[++i]);
                    if (options.k < 1) throw invalid_argument("k < 1");
                } catch (...) {
                    error = "Invalid value for " + argument + ": " + argv[i];
                    return false;
                }
            }
        } else if (argument == "--local-k-output" || argument == "--local-k-report" || argument == "--local-k-start-score-drop") {
            if (i + 1 >= argc) {
                error = "Missing value for " + argument;
                return false;
            }
            const string value = argv[++i];
            if (argument == "--local-k-output") options.output_filename = value;
            else if (argument == "--local-k-report") options.report_filename = value;
            else {
                try {
                    options.start_score_drop = stod(value);
                    if (!isfinite(options.start_score_drop) ||
                        options.start_score_drop < 0.0) throw invalid_argument("range");
                    has_drop = true;
                } catch (...) {
                    error = "Invalid value for " + argument + ": " + value;
                    return false;
                }
            }
        } else if (argument == "--no-dp-dump") {
            options.no_dp_dump = true;
        } else if (argument == "--alternative-splicing" || argument == "--alt-include-incomplete" || argument == "--alt-debug") {
            // handled elsewhere
        } else if (argument == "--alt-output" || argument == "--alt-report" ||
                   argument == "--alt-max-distance" || argument == "--alt-min-score-delta") {
            ++i; // skip value
        } else {
            error = "Unknown option " + argument;
            return false;
        }
    }
    if (!options.enabled && (!options.output_filename.empty() || !options.report_filename.empty() || has_drop)) {
        error = "local-k options require --local-k-best";
        return false;
    }
    return true;
}

static bool parse_alternative_options(int argc, char **argv,
                                      AlternativeOptions &options,
                                      string &error) {
    for (int i = 7; i < argc; ++i) {
        const string argument = argv[i];
        if (argument == "--alternative-splicing") {
            options.enabled = true;
        } else if (argument == "--alt-output" || argument == "--alt-report" ||
                   argument == "--alt-max-distance" ||
                   argument == "--alt-min-score-delta") {
            if (i + 1 >= argc) {
                error = "Missing value for " + argument;
                return false;
            }
            const string value = argv[++i];
            try {
                if (argument == "--alt-output") options.output_filename = value;
                else if (argument == "--alt-report") options.report_filename = value;
                else if (argument == "--alt-max-distance") {
                    options.max_distance = stoi(value);
                    if (options.max_distance <= 0) throw invalid_argument("range");
                } else options.min_score_delta = stod(value);
            } catch (const exception &) {
                error = "Invalid value for " + argument + ": " + value;
                return false;
            }
        } else if (argument == "--alt-include-incomplete") {
            options.include_incomplete = true;
        } else if (argument == "--alt-debug") {
            options.debug = true;
        } else if (argument == "--local-k-best") {
            if (i + 1 < argc && argv[i+1][0] != '-') {
                ++i;
            }
        } else if (argument == "--no-dp-dump") {
            // handled elsewhere
        } else if (argument == "--local-k-output" || argument == "--local-k-report" || argument == "--local-k-start-score-drop") {
            if (i + 1 < argc) ++i; // skip value if present
        } else {
            error = "Unknown option " + argument;
            return false;
        }
    }
    if (!options.enabled && (!options.output_filename.empty() ||
        !options.report_filename.empty() || options.include_incomplete ||
        options.debug)) {
        error = "Alternative-splicing options require --alternative-splicing";
        return false;
    }
    return true;
}

static string default_alternative_output(const string &fasta) {
    return fasta + ".alternative_splicing.gff";
}

//------------------------------------------------------------
// MAIN
//------------------------------------------------------------
struct LocalKBestCell {
    double dp = NEG_INF;
    int bt_state = -1;
    int bt_k = -1;
    PathMetadata metadata;
};

struct PathCandidate {
    double score;
    int from_state;
    int from_k;
    TransitionResult evaluated;

    bool operator<(const PathCandidate &other) const {
        if (score != other.score) return score > other.score;
        if (from_state != other.from_state) return from_state < other.from_state;
        return from_k < other.from_k;
    }
};

struct LocalPathResult {
    vector<int> path_positions;
    vector<int> path_states;
    double score;
    double actual_start_score = NEG_INF;
    double inherited_start_score = NEG_INF;
};


struct LocalKBestCandidate {
    int reference_gene_index;
    int rank;
    double score;
    double score_delta;
    double actual_start_score;
    double inherited_start_score;
    string classification;
    bool identical_to_reference;
    TranscriptAnnotation transcript;
    bool kept = false;
};

static bool run_local_k_best(
    const string &seqid,
    const vector<TranscriptAnnotation> &references,
    const vector<int> &global_path_states,
    const ModelInputs &inputs,
    const LocalKBestOptions &options
) {
    int L = inputs.seq.size();
    int K = options.k;

    ofstream gff_out;
    if (!options.output_filename.empty()) {
        gff_out.open(options.output_filename);
        if (!gff_out) {
            cerr << "Error opening " << options.output_filename << " for writing\n";
            return false;
        }
        gff_out << "##gff-version 3\n";
    }

    ofstream tsv_out;
    if (!options.report_filename.empty()) {
        tsv_out.open(options.report_filename);
        if (!tsv_out) {
            cerr << "Error opening " << options.report_filename << " for writing\n";
            return false;
        }
        tsv_out << "reference_gene_index\trank\treference_start_1based\t"
                << "reference_end_1based\tcandidate_start_1based\t"
                << "candidate_end_1based\tfixed_interval_score\t"
                << "reference_fixed_interval_score\tscore_delta\t"
                << "actual_start_score\tinherited_start_score\tintron_chain\t"
                << "reference_intron_chain\t"
                << "identical_to_reference\tclassification\n";
    }

    struct ActiveCell {
        double dp = NEG_INF;
        PathMetadata metadata;
    };

    struct CompactBP {
        int16_t bt_state;
        int16_t bt_k;
    };

    for (size_t ref_idx = 0; ref_idx < references.size(); ++ref_idx) {
        const auto &ref = references[ref_idx];

        int L_bound = (ref_idx == 0) ? 0 : references[ref_idx - 1].path_end;
        int R_bound = (ref_idx == references.size() - 1) ? L - 1 : references[ref_idx + 1].path_start - 1;
        int local_len = R_bound - L_bound + 1;

        if (local_len <= 0) continue;

        double ref_atg_score = (ref.path_start < 25) ? 1.0 : inputs.atg_score[codon_start_0based(ref.path_start)];

        double ref_fixed_score = inputs.emit[L_bound][0];
        PathMetadata ref_meta;
        ref_meta.predecessor_state = 0;
        ref_meta.inter_len = 1;
        int ref_cur_state = 0;
        for (int p = L_bound + 1; p <= R_bound; p++) {
            int next_state = global_path_states[p];
            TransitionResult eval = evaluate_transition(p, ref_cur_state, next_state, ref_meta, inputs);
            ref_fixed_score += eval.transition_score + eval.emission_score;
            ref_meta = eval.next_metadata;
            ref_cur_state = next_state;
        }

        vector<int> candidate_starts;
        for (int p = max(L_bound, 2); p <= ref.path_start; p++) {
            if (p == ref.path_start) {
                candidate_starts.push_back(p);
            } else {
                double score = inputs.atg_score[codon_start_0based(p)];
                if (score > NEG_INF && score >= ref_atg_score - options.start_score_drop) {
                    candidate_starts.push_back(p);
                }
            }
        }

        vector<vector<ActiveCell>> prev_dp(8, vector<ActiveCell>(K));
        vector<vector<ActiveCell>> curr_dp(8, vector<ActiveCell>(K));
        vector<CompactBP> backpointers(local_len * 8 * K, CompactBP{-1, -1});

        prev_dp[0][0].dp = 0.0 + inputs.emit[L_bound][0];
        prev_dp[0][0].metadata.inter_len = 1;
        prev_dp[0][0].metadata.predecessor_state = 0;

        for (int i = 1; i < local_len; i++) {
            int p = L_bound + i;
            bool is_cand_start = (find(candidate_starts.begin(), candidate_starts.end(), p) != candidate_starts.end());

            for (int to = 0; to < 8; to++) {
                vector<PathCandidate> candidates;

                for (int from = 0; from < 8; from++) {
                    if (to == 0 && from != 0) continue;
                    if (to == 7 && from != 7 && !is_exon(from)) continue;
                    if (to == 7 && is_exon(from) && p < ref.path_start)
                        continue;
                    if (is_exon(to) && from == 7) continue;
                    if (is_intron(to) && from == 7) continue;
                    if (from == 0 && is_exon(to) && !is_cand_start) continue;

                    int to_orig = (to == 7) ? 0 : to;
                    int from_orig = (from == 7) ? 0 : from;

                    for (int k = 0; k < K; k++) {
                        if (prev_dp[from][k].dp <= NEG_INF) continue;

                        TransitionResult evaluated = evaluate_transition(p, from_orig, to_orig, prev_dp[from][k].metadata, inputs);

                        if (from == 0 && is_exon(to) && is_cand_start && evaluated.allowed) {
                            evaluated.transition_score = ref_atg_score;
                        }

                        if (evaluated.allowed) {
                            double score = prev_dp[from][k].dp + evaluated.transition_score + evaluated.emission_score;
                            candidates.push_back({score, from, k, evaluated});
                        }
                    }
                }

                sort(candidates.begin(), candidates.end());
                for (int k = 0; k < min(K, (int)candidates.size()); k++) {
                    curr_dp[to][k].dp = candidates[k].score;
                    curr_dp[to][k].metadata = candidates[k].evaluated.next_metadata;
                    backpointers[i * 8 * K + to * K + k] = CompactBP{static_cast<int16_t>(candidates[k].from_state), static_cast<int16_t>(candidates[k].from_k)};
                }
            }
            prev_dp = curr_dp;
            for (int to = 0; to < 8; to++) {
                for (int k = 0; k < K; k++) {
                    curr_dp[to][k].dp = NEG_INF;
                }
            }
        }

        vector<LocalPathResult> paths;
        for (int k = 0; k < K; k++) {
            if (prev_dp[7][k].dp <= NEG_INF) continue;

            LocalPathResult res;
            res.score = prev_dp[7][k].dp;

            int cur_state = 7;
            int cur_k = k;
            bool valid = true;

            for (int i = local_len - 1; i >= 0; i--) {
                int p = L_bound + i;
                res.path_positions.push_back(p);
                res.path_states.push_back(cur_state == 7 ? 0 : cur_state);

                if (i > 0) {
                    CompactBP bp = backpointers[i * 8 * K + cur_state * K + cur_k];
                    int next_state = bp.bt_state;
                    int next_k = bp.bt_k;
                    if (next_state < 0) { valid = false; break; }

                    if (next_state == 0 && is_exon(cur_state)) {
                        res.actual_start_score = inputs.atg_score[codon_start_0based(p)];
                        res.inherited_start_score = ref_atg_score;
                    }

                    cur_state = next_state;
                    cur_k = next_k;
                }
            }

            if (valid) {
                reverse(res.path_positions.begin(), res.path_positions.end());
                reverse(res.path_states.begin(), res.path_states.end());
                paths.push_back(res);
            }
        }

        vector<LocalKBestCandidate> unique_candidates;
        set<string> seen;

        for (const auto &path : paths) {
            auto anns = build_transcript_annotations(path.path_positions, path.path_states, seqid);
            if (anns.empty()) continue;
            TranscriptAnnotation ann = anns[0];
            ann.score = path.score;

            string key = full_transcript_key(ann);
            if (seen.insert(key).second) {
                LocalKBestCandidate cand;
                cand.transcript = ann;
                cand.score = path.score;
                cand.score_delta = path.score - ref_fixed_score;
                cand.actual_start_score = path.actual_start_score;
                cand.inherited_start_score = path.inherited_start_score;
                cand.classification = classify_splice_difference(ref, ann);
                cand.identical_to_reference = (key == full_transcript_key(ref));
                unique_candidates.push_back(cand);
            }
            if (unique_candidates.size() >= static_cast<size_t>(K)) break;
        }


        if (unique_candidates.size() > 0 && gff_out) {
            int min_start = unique_candidates[0].transcript.genomic_start_0based;
            int max_end = unique_candidates[0].transcript.genomic_end_0based_exclusive;
            for (const auto &cand : unique_candidates) {
                min_start = min(min_start, cand.transcript.genomic_start_0based);
                max_end = max(max_end, cand.transcript.genomic_end_0based_exclusive);
            }
            string g_id = "gene" + to_string(ref_idx + 1);
            gff_out << ref.sequence_id << "\tUniAnn\tgene\t"
                << min_start + 1 << '\t'
                << max_end << "\t.\t+\t.\tID="
                << g_id << "\n";
        }

        for (size_t rank = 0; rank < unique_candidates.size(); ++rank) {

            const auto &cand = unique_candidates[rank];
            string t_id = "gene" + to_string(ref_idx + 1) + ".local_k" + to_string(rank + 1);
            string g_id = "gene" + to_string(ref_idx + 1);

            if (gff_out) {
                write_hierarchical_transcript(gff_out, cand.transcript, g_id, t_id);
            }

            if (tsv_out) {
                tsv_out << (ref_idx + 1) << '\t'
                        << (rank + 1) << '\t'
                        << (ref.genomic_start_0based + 1) << '\t'
                        << ref.genomic_end_0based_exclusive << '\t'
                        << (cand.transcript.genomic_start_0based + 1) << '\t'
                        << cand.transcript.genomic_end_0based_exclusive << '\t'
                        << cand.score << '\t'
                        << ref_fixed_score << '\t'
                        << cand.score_delta << '\t'
                        << cand.actual_start_score << '\t'
                        << cand.inherited_start_score << '\t'
                        << tsv_escape_chain(cand.transcript) << '\t'
                        << tsv_escape_chain(ref) << '\t'
                        << (cand.identical_to_reference ? "1" : "0") << '\t'
                        << cand.classification << '\n';
            }
        }
    }
    if (!options.output_filename.empty() && !gff_out.good()) return false;
    if (!options.report_filename.empty() && !tsv_out.good()) return false;
    return true;
}


int main(int argc, char** argv) {

    if (argc < 7) {
        cerr << "Usage: " << argv[0]
             << " seq.fasta emissions.txt gt.txt ag.txt atg.txt stop.txt"
             << " [--alternative-splicing] [--alt-output FILE]"
             << " [--alt-report FILE] [--alt-max-distance BP]"
             << " [--alt-min-score-delta SCORE] [--alt-include-incomplete]"
             << " [--alt-debug]\n"
             << "             [--local-k-best K] [--local-k-output FILE]\n"
             << "             [--local-k-report FILE] [--local-k-start-score-drop SCORE]\n"
             << "             [--no-dp-dump]\n";
        return 1;
    }


    LocalKBestOptions local_k_options;
    string local_k_error;
    if (!parse_local_k_best_options(argc, argv, local_k_options, local_k_error)) {
        cerr << local_k_error << "\n";
        return 1;
    }
    AlternativeOptions alternative_options;
    string option_error;
    if (!parse_alternative_options(argc, argv, alternative_options, option_error)) {
        cerr << option_error << "\n";
        return 1;
    }

    string f_fasta = argv[1];
    string f_emit  = argv[2];
    string f_gt    = argv[3];
    string f_ag    = argv[4];
    string f_atg   = argv[5];
    string f_stop   = argv[6];

    //--------------------------------------------------------
    // Load FASTA
    //--------------------------------------------------------
    string seq_str = read_fasta(f_fasta);
    int L = seq_str.size();
    if (L == 0) {
        cerr << "FASTA contains no sequence bases " << f_fasta << "\n";
        return 1;
    }

    vector<char> seq(L);
    for (int i = 0; i < L; i++)
        seq[i] = seq_str[i];

    //--------------------------------------------------------
    // Load emissions and splice scores
    //--------------------------------------------------------
    auto emit     =  load_emissions(f_emit, L);
    auto gt_score =  load_sparse_scores(f_gt, L);
    auto ag_score =  load_sparse_scores(f_ag, L);
    auto atg_score = load_sparse_scores(f_atg, L);
    auto stop_score = load_sparse_scores(f_stop, L);

    //--------------------------------------------------------
    // Safety check: GT/AG coordinates match sequence
    //--------------------------------------------------------
    safety_check_gt_ag_atg(seq, gt_score, ag_score, atg_score, stop_score);

    //--------------------------------------------------------
    // Initialize transitions
    //--------------------------------------------------------
    auto trans = init_transitions();

    //--------------------------------------------------------
    // Initialize DP
    //--------------------------------------------------------
    auto dp = init_dp(L, emit);

    //--------------------------------------------------------
    // Run full Viterbi
    //--------------------------------------------------------
    run_viterbi(dp, emit, gt_score, ag_score, atg_score, stop_score, seq, trans);

    if (!local_k_options.no_dp_dump) {
        //print DP and BT matrices
        for (int i = 0; i < L; i++) {
          fprintf(stderr,"%d\tdp\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",i,int(dp[i][0].dp),int(dp[i][1].dp),int(dp[i][2].dp),int(dp[i][3].dp),int(dp[i][4].dp),int(dp[i][5].dp),int(dp[i][6].dp));
          fprintf(stderr,"%d\tbt\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",i,int(dp[i][0].bt),int(dp[i][1].bt),int(dp[i][2].bt),int(dp[i][3].bt),int(dp[i][4].bt),int(dp[i][5].bt),int(dp[i][6].bt));
        }
    }


    //--------------------------------------------------------
    // Termination
    //--------------------------------------------------------
    auto [best_final, best_state] = viterbi_termination(dp, L);

    //--------------------------------------------------------
    // Backtrace
    //--------------------------------------------------------
    auto path_states = viterbi_backtrace(dp, L, best_state);
    auto path_scores = viterbi_backtrace_scores(dp, L, best_state);
    auto path_labels = states_to_labels(path_states);

    //--------------------------------------------------------
    // Write GFF3 output
    //--------------------------------------------------------
    string seqid = get_fasta_header(f_fasta);
    write_gff_from_path(path_labels, path_scores, seqid, f_fasta, best_final);

    if (local_k_options.enabled) {
        const ModelInputs inputs{emit, gt_score, ag_score, atg_score, stop_score, seq, trans};
        auto references = build_transcript_annotations(path_states, seqid);
        for (auto &reference : references) {
            if (reference.path_start > 0 && reference.path_end >= reference.path_start)
                reference.score = dp[reference.path_end][0].dp -
                                  dp[reference.path_start - 1][0].dp;
        }
        if (!run_local_k_best(seqid, references, path_states, inputs,
                              local_k_options)) {
            return 1;
        }
    }


    if (alternative_options.enabled) {
        const auto started = chrono::steady_clock::now();
        const ModelInputs inputs{emit, gt_score, ag_score, atg_score,
                                 stop_score, seq, trans};
        auto references = build_transcript_annotations(path_states, seqid);
        for (auto &reference : references) {
            if (reference.path_start > 0 && reference.path_end >= reference.path_start)
                reference.score = dp[reference.path_end][0].dp -
                                  dp[reference.path_start - 1][0].dp;
        }
        auto anchors = extract_splice_anchors(path_states, references, seqid,
                                              inputs, alternative_options.debug);
        vector<AlternativeCandidate> candidates;
        candidates.reserve(anchors.size());
        for (const auto &anchor : anchors) {
            auto candidate = search_anchor(anchor, path_states, dp, references,
                                           inputs, alternative_options);
            if (alternative_options.debug) {
                cerr << "ALT candidate anchor="
                     << anchor.genomic_position_0based
                     << " status=" << candidate_status_name(candidate.status)
                     << " score=" << candidate.total_score
                     << " reference=" << candidate.reference_interval_score
                     << " delta=" << candidate.score_delta
                     << " reason=" << candidate.failure_reason << "\n";
            }
            candidates.push_back(move(candidate));
        }
        const vector<int> kept = deduplicate_candidates(
            candidates, references, alternative_options.min_score_delta);
        const string output_filename = alternative_options.output_filename.empty()
            ? default_alternative_output(f_fasta)
            : alternative_options.output_filename;
        if (!write_alternative_gff3(output_filename, references,
                                    candidates, kept)) {
            cerr << "Cannot write alternative annotation "
                 << output_filename << "\n";
            return 1;
        }
        if (!alternative_options.report_filename.empty() &&
            !write_alt_report(alternative_options.report_filename, candidates,
                              alternative_options.include_incomplete)) {
            cerr << "Cannot write alternative report "
                 << alternative_options.report_filename << "\n";
            return 1;
        }
        print_alternative_summary(references.size(), anchors, candidates, kept);
        if (alternative_options.debug) {
            const double seconds = chrono::duration<double>(
                chrono::steady_clock::now() - started).count();
            cerr << "ALT elapsed_seconds=" << seconds << "\n";
        }
    }

    return 0;
}

