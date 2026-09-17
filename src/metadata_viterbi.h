#pragma once

#include "metadata_states.h"
#include "shared_trace.h"
#include <stdexcept>

namespace uniann {

struct ScoredRun {
    int state;
    int start; // Inclusive, zero-based genomic positions.
    int end;
    double score_at_start = 0;
};

struct MetadataPath {
    double score = 0;
    std::vector<ScoredRun> runs;
    size_t max_active_states = 1;
    size_t peak_trace_nodes = 0;
    size_t trace_capacity_bytes = 0;
};

// Recover only the boundary scores needed by the GFF writer. Raw metadata
// replay also checks that state canonicalization preserved every chosen edge.
inline void score_metadata_path(MetadataPath &path, const ModelInputs &inputs) {
    double score = inputs.emit[0][0];
    PathMetadata metadata{0, 0, 1, 0, -1};
    int previous_state = 0;
    for (auto &run : path.runs) {
        for (int position = run.start; position <= run.end; ++position) {
            if (position > 0) {
                const auto edge = evaluate_transition(position, previous_state, run.state, metadata, inputs);
                if (!edge.allowed) throw std::logic_error("Illegal transition in metadata traceback");
                score = score + edge.transition_score + edge.emission_score;
                metadata = edge.next_metadata;
            }
            if (position == run.start) run.score_at_start = score;
            previous_state = run.state;
        }
    }
    if (score != path.score) throw std::logic_error("Metadata traceback score differs from optimum");
}

// Keep two score frontiers and the ancestors of their surviving hypotheses.
// Completed common prefixes become runs in the final path; no chromosome-wide
// DP, backpointer matrix, per-base score path or string-label path is retained.
inline MetadataPath decode_metadata_viterbi(const ModelInputs &inputs, bool compress_runs = true) {
    if (inputs.seq.empty() || !std::isfinite(inputs.emit[0][0]) || inputs.emit[0][0] <= NEG_INF)
        throw std::invalid_argument("Metadata Viterbi requires a valid initial N emission");
    const int length = inputs.seq.size();
    const auto states = make_metadata_states();
    const int initial = initial_metadata_state();
    // Two frontiers swapped by pointer; only active entries are ever written
    // or reset, so a step costs O(active states) rather than O(slots).
    MetadataScores frontier[2];
    frontier[0].fill(-std::numeric_limits<double>::infinity());
    frontier[1].fill(-std::numeric_limits<double>::infinity());
    MetadataScores *previous = &frontier[0], *next = &frontier[1];
    (*previous)[initial] = inputs.emit[0][0];
    ActiveMetadataStates active_from{static_cast<uint16_t>(initial)}, active_to;
    active_from.reserve(METADATA_STATE_COUNT);
    active_to.reserve(METADATA_STATE_COUNT);

    MetadataPath path;
    uint64_t emitted_positions = 0;
    SharedTrace trace(METADATA_STATE_COUNT, initial,
        [&](uint8_t state, uint32_t start, uint32_t end) {
            if (start != emitted_positions || end < start)
                throw std::logic_error("Noncontiguous metadata traceback");
            emitted_positions = uint64_t(end) + 1;
            if (!path.runs.empty() && path.runs.back().state == state) {
                path.runs.back().end = end;
            } else {
                path.runs.push_back({state, static_cast<int>(start), static_cast<int>(end), 0});
            }
        });

    std::vector<uint8_t> coarse_states(METADATA_STATE_COUNT);
    std::vector<int> parent_slots(METADATA_STATE_COUNT, -1);
    std::vector<int> trace_slots;
    trace_slots.reserve(METADATA_STATE_COUNT);
    for (int id = 0; id < METADATA_STATE_COUNT; ++id) coarse_states[id] = states[id].state;
    MetadataParents parents;
    constexpr int flush_interval = 4096;
    for (int position = 1; position < length; ++position) {
        advance_metadata_scores(position, *previous, active_from, *next, active_to, states, inputs, parents);
        if (active_to.empty()) throw std::runtime_error("No valid metadata path reaches this position");
        trace_slots.clear();
        for (const uint16_t id : active_to) {
            parent_slots[id] = parents[id];
            trace_slots.push_back(id);
        }
        path.max_active_states = std::max(path.max_active_states, active_to.size());
        trace.advance(trace_slots, parent_slots, coarse_states, position, compress_runs);
        for (const uint16_t id : active_from) (*previous)[id] = -std::numeric_limits<double>::infinity();
        std::swap(previous, next);
        active_from.swap(active_to);
        if (position % flush_interval == 0) trace.flush_common();
    }

    int winner = active_from.front();
    for (const uint16_t id : active_from)
        if ((*previous)[id] > (*previous)[winner]) winner = id;
    path.score = (*previous)[winner];
    path.peak_trace_nodes = trace.peak_live_nodes();
    path.trace_capacity_bytes = trace.pool_capacity() * trace.bytes_per_node();
    trace.finish(winner);
    if (trace.live_nodes() != 0 || emitted_positions != inputs.seq.size())
        throw std::logic_error("Incomplete metadata traceback");
    score_metadata_path(path, inputs);
    return path;
}

} // namespace uniann
