// Independent reference retains every raw metadata key: no length caps or IDs.
#include "../src/metadata_viterbi.h"
#include <cassert>
#include <iostream>
#include <map>
#include <random>
#include <tuple>
using namespace std;
using namespace uniann;

vector<vector<double>> init_transitions() {
    vector<vector<double>> transitions(NUM_STATES, vector<double>(NUM_STATES, NEG_INF));
    for (int state = 0; state < NUM_STATES; ++state) transitions[state][state] = 0;
    return transitions;
}

// Fixtures still describe scores per motif file; fold them into the single
// per-position array the decoder reads, keeping only sequence-backed entries.
vector<double> merged_site_scores(const vector<char> &seq, const vector<double> &gt,
                                  const vector<double> &ag, const vector<double> &atg,
                                  const vector<double> &stop) {
    vector<double> site(seq.size(), NEG_INF);
    for (int p = 0; p < static_cast<int>(seq.size()); ++p) {
        switch (motif_starting_at(seq, p)) {
            case SiteMotif::Donor: site[p] = gt[p]; break;
            case SiteMotif::Acceptor: site[p] = ag[p]; break;
            case SiteMotif::Start: site[p] = atg[p]; break;
            case SiteMotif::Stop: site[p] = stop[p]; break;
            default: break;
        }
    }
    return site;
}

using RawKey = tuple<int, int, int, int, int, int>;

double raw_metadata_optimum(const ModelInputs &inputs) {
    map<RawKey, double> frontier;
    frontier[{0, 0, 0, 1, 0, -1}] = inputs.emit[0][0];
    for (int p = 1; p < static_cast<int>(inputs.seq.size()); ++p) {
        map<RawKey, double> next;
        for (const auto &entry : frontier) {
            const auto &[from, il, el, nl, ef, pred] = entry.first;
            PathMetadata metadata{il, el, nl, ef, pred};
            for (int to = 0; to < NUM_STATES; ++to) {
                const auto e = evaluate_transition(p, from, to, metadata, inputs);
                if (!e.allowed) continue;
                const auto &m = e.next_metadata;
                const RawKey key{to, m.intron_len, m.exon_len, m.inter_len,
                                 m.exon_from, m.predecessor_state};
                const double score = entry.second + e.transition_score + e.emission_score;
                auto it = next.find(key);
                if (it == next.end()) next.emplace(key, score);
                else if (score > it->second) it->second = score;
            }
        }
        frontier.swap(next);
    }
    double best = -numeric_limits<double>::infinity();
    for (const auto &entry : frontier) best = max(best, entry.second);
    return best;
}

void check_instance(int length, unsigned seed, bool negative = false) {
    mt19937 rng(seed);
    vector<char> seq(length, 'A');
    const vector<string> motifs{"ATG", "GT", "AG", "TAA", "TAG", "TGA"};
    if (length < 1000) {
        for (int p = 0; p + 3 < length; p += 7) {
            const auto &m = motifs[rng() % motifs.size()];
            copy(m.begin(), m.end(), seq.begin() + p);
        }
    } else {
        // A splice crosses each boundary of the 4096-base future traceback flush intervals.
        for (int boundary = 4096; boundary < length + 50; boundary += 4096) {
            for (const auto &entry : vector<pair<int,string>>{
                     {boundary - 120, "ATG"}, {boundary - 15, "GT"},
                     {boundary + 37, "AG"}, {boundary + 102, "TAA"}}) {
                if (entry.first >= 0 && entry.first + 3 < length)
                    copy(entry.second.begin(), entry.second.end(), seq.begin() + entry.first);
            }
        }
    }
    vector<EmissionRow> emit(length);
    vector<double> gt(length, NEG_INF), ag = gt, atg = gt, stop = gt;
    for (int p = 0; p < length; ++p) {
        for (int s = 0; s < EMISSION_COLUMNS; ++s)
            emit[p][s] = negative ? -200000000.0f : static_cast<int>(rng()%17) - 8;
        if (p + 1 < length && seq[p] == 'G' && seq[p+1] == 'T') gt[p] = rng()%30;
        if (p + 1 < length && seq[p] == 'A' && seq[p+1] == 'G') ag[p] = rng()%30;
        if (p + 2 < length) {
            string codon(seq.begin() + p, seq.begin() + p + 3);
            if (codon == "ATG") atg[p] = rng()%30;
            if (is_stop_codon(codon)) stop[p] = rng()%30;
        }
    }
    const auto trans = init_transitions();
    const auto site = merged_site_scores(seq, gt, ag, atg, stop);
    const ModelInputs inputs{emit, site, seq, trans};
    const double expected = raw_metadata_optimum(inputs);
    const auto states = make_metadata_states();
    MetadataScores previous, next;
    previous.fill(-numeric_limits<double>::infinity());
    previous[initial_metadata_state()] = emit[0][0];
    // An independent, deliberately dense traceback is cheap for these small
    // fixtures and detects mistakes in sharing, pruning, flushing and RLE.
    vector<MetadataParents> parents(length);
    for (int p = 1; p < length; ++p) {
        advance_metadata_scores(p, previous, next, states, inputs, &parents[p]);
        previous.swap(next);
    }
    int current = max_element(previous.begin(), previous.end()) - previous.begin();
    assert(previous[current] == expected);
    vector<int> dense_path(length);
    for (int p = length - 1; p > 0; --p) {
        dense_path[p] = states[current].state;
        current = parents[p][current];
        assert(current != NO_METADATA_PARENT);
    }
    assert(current == initial_metadata_state());
    vector<double> dense_scores(length);
    dense_scores[0] = emit[0][0];
    PathMetadata replay_metadata{0, 0, 1, 0, -1};
    for (int p = 1; p < length; ++p) {
        const auto edge = evaluate_transition(p, dense_path[p - 1], dense_path[p], replay_metadata, inputs);
        assert(edge.allowed);
        dense_scores[p] = dense_scores[p - 1] + edge.transition_score + edge.emission_score;
        replay_metadata = edge.next_metadata;
    }
    for (bool compress : {false, true}) {
        const auto path = decode_metadata_viterbi(inputs, compress);
        assert(path.score == expected);
        int covered = 0;
        for (const auto &run : path.runs) {
            assert(run.start == covered);
            assert(run.score_at_start == dense_scores[run.start]);
            for (int p = run.start; p <= run.end; ++p) assert(run.state == dense_path[p]);
            covered = run.end + 1;
        }
        assert(covered == length);
    }
}

void check_state_ids() {
    const auto states = make_metadata_states();
    assert(METADATA_STATE_COUNT == 511);
    for (int id = 0; id < METADATA_STATE_COUNT; ++id) {
        assert(metadata_state_id(states[id].state, states[id].metadata) == id);
        assert(states[id].self_destination >= 0);
        assert(states[id].self_destination < METADATA_STATE_COUNT);
    }
}


void check_older_intron_survives() {
    const int length = 200;
    vector<char> seq(length, 'A');
    for (const auto &entry : vector<pair<int, string>>{
             {0, "ATG"}, {30, "GT"}, {60, "GT"}, {80, "AG"}, {151, "TAA"}})
        copy(entry.second.begin(), entry.second.end(), seq.begin() + entry.first);
    vector<EmissionRow> emit(length);
    for (int p = 0; p < length; ++p) {
        emit[p].fill(-1000);
        emit[p][0] = 0;
        emit[p][4] = 0;
        if (p >= 2 && p <= 60) emit[p][1] = 1;
        if (p >= 81 && p <= 152) emit[p][2] = 10;
    }
    vector<double> gt(length, NEG_INF), ag = gt, atg = gt, stop = gt;
    gt[30] = 0; gt[60] = 100; ag[80] = 100; atg[0] = 100; stop[151] = 100;
    const auto trans = init_transitions();
    const auto site = merged_site_scores(seq, gt, ag, atg, stop);
    const ModelInputs inputs{emit, site, seq, trans};
    const auto states = make_metadata_states();
    MetadataScores previous, next;
    previous.fill(-numeric_limits<double>::infinity());
    previous[initial_metadata_state()] = emit[0][0];
    for (int p = 1; p < length; ++p) {
        advance_metadata_scores(p, previous, next, states, inputs);
        previous.swap(next);
        if (p == 80) {
            // The stronger donor wins in a single I0 cell, but its intron is
            // too short for this acceptor. Both histories must still exist.
            const int older = metadata_state_id(4, canonicalize_metadata(4, {50, 0, 0, 0, -1}));
            const int younger = metadata_state_id(4, canonicalize_metadata(4, {20, 0, 0, 0, -1}));
            assert(isfinite(previous[older]));
            assert(previous[older] < previous[younger]);
        }
    }
    assert(*max_element(previous.begin(), previous.end()) == 950);
    assert(raw_metadata_optimum(inputs) == 950);
}

int main() {
    check_older_intron_survives();
    check_state_ids();
    for (unsigned seed = 0; seed < 80; ++seed) check_instance(160, seed);
    for (int length : {1, 2, 4095, 4096, 4097, 8192, 8193, 8400})
        check_instance(length, 91);
    check_instance(40, 92, true);
    cout << "PASS raw metadata oracle: 80 randomized, 8 boundary, 1 below-sentinel cases\n";
}
