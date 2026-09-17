#pragma once

#include "viterbi_model.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace uniann {

// Keep only metadata that can affect a future transition. Lengths beyond a
// minimum collapse, except intron length must still distinguish reading frames.
inline PathMetadata canonicalize_metadata(int state, const PathMetadata &metadata) {
    PathMetadata key;
    if (state == 0) {
        key.inter_len = std::min(metadata.inter_len, MIN_INTER - 2);
        key.predecessor_state = metadata.predecessor_state == 0 ? 0 : 1;
    } else if (is_exon(state)) {
        key.exon_from = metadata.exon_from == 0 ? 0 : 4;
        const int cap = key.exon_from == 0 ? MIN_SINGLE + 3 : MIN_EXON + 2;
        key.exon_len = std::min(metadata.exon_len, cap);
    } else {
        const int threshold = MIN_INTRON - 2;
        key.intron_len = metadata.intron_len < threshold
            ? metadata.intron_len : threshold + (metadata.intron_len - threshold) % 3;
    }
    return key;
}

// Reserve the unused zero-length entries to keep the ID mapping arithmetic.
// These widths describe storage slots, not the number of active hypotheses.
constexpr int METADATA_N_WIDTH = 2 * (MIN_INTER - 1);
constexpr int METADATA_FIRST_EXON_WIDTH = MIN_SINGLE + 4;
constexpr int METADATA_EXON_WIDTH = METADATA_FIRST_EXON_WIDTH + MIN_EXON + 3;
constexpr int METADATA_INTRON_WIDTH = MIN_INTRON + 1;
constexpr int METADATA_INTRON_OFFSET = METADATA_N_WIDTH + 3 * METADATA_EXON_WIDTH;
constexpr int METADATA_STATE_COUNT = METADATA_INTRON_OFFSET + 3 * METADATA_INTRON_WIDTH;
constexpr uint16_t NO_METADATA_PARENT = std::numeric_limits<uint16_t>::max();
static_assert(METADATA_STATE_COUNT < NO_METADATA_PARENT, "metadata IDs must fit uint16_t");
using MetadataScores = std::array<double, METADATA_STATE_COUNT>;
using MetadataParents = std::array<uint16_t, METADATA_STATE_COUNT>;

struct MetadataState {
    int state;
    PathMetadata metadata;
    int self_destination;
};
using MetadataStates = std::array<MetadataState, METADATA_STATE_COUNT>;

inline int metadata_state_id(int state, const PathMetadata &key) {
    if (state == 0)
        return key.inter_len * 2 + (key.predecessor_state == 0 ? 0 : 1);
    if (is_exon(state))
        return METADATA_N_WIDTH + (state - 1) * METADATA_EXON_WIDTH +
            (key.exon_from == 0 ? 0 : METADATA_FIRST_EXON_WIDTH) + key.exon_len;
    return METADATA_INTRON_OFFSET + (state - 4) * METADATA_INTRON_WIDTH + key.intron_len;
}

inline MetadataStates make_metadata_states() {
    MetadataStates states;
    for (int id = 0; id < METADATA_STATE_COUNT; ++id) {
        int state;
        PathMetadata metadata;
        if (id < METADATA_N_WIDTH) {
            state = 0;
            metadata.inter_len = id / 2;
            metadata.predecessor_state = id % 2;
        } else if (id < METADATA_INTRON_OFFSET) {
            const int offset = id - METADATA_N_WIDTH;
            state = 1 + offset / METADATA_EXON_WIDTH;
            const int length = offset % METADATA_EXON_WIDTH;
            metadata.exon_from = length < METADATA_FIRST_EXON_WIDTH ? 0 : 4;
            metadata.exon_len = length < METADATA_FIRST_EXON_WIDTH
                ? length : length - METADATA_FIRST_EXON_WIDTH;
        } else {
            const int offset = id - METADATA_INTRON_OFFSET;
            state = 4 + offset / METADATA_INTRON_WIDTH;
            metadata.intron_len = offset % METADATA_INTRON_WIDTH;
        }
        PathMetadata next = metadata;
        if (state == 0) {
            ++next.inter_len;
            next.predecessor_state = 0;
        } else if (is_exon(state)) {
            ++next.exon_len;
        } else {
            ++next.intron_len;
        }
        states[id] = {state, metadata, metadata_state_id(state, canonicalize_metadata(state, next))};
    }
    return states;
}

inline int initial_metadata_state() {
    return metadata_state_id(0, canonicalize_metadata(0, PathMetadata{0, 0, 1, 0, -1}));
}

// Ascending IDs with finite scores in a frontier. Iterating this list instead
// of all METADATA_STATE_COUNT slots keeps each step proportional to the number
// of surviving hypotheses, which is usually far below the slot count.
using ActiveMetadataStates = std::vector<uint16_t>;

// One exact Viterbi step. A state ID includes metadata, so a higher-scoring
// short intron cannot discard an older intron that can use the next acceptor.
//
// `active_from` must list every finite entry of `previous` in ascending order
// and `next` must hold -infinity everywhere on entry. On return `active_to`
// lists the finite entries of `next` in ascending order, and `parents` is
// valid for exactly those entries. Visiting sources in ascending order with
// strict improvement keeps ties repeatable.
inline void advance_metadata_scores(
    int position, const MetadataScores &previous, const ActiveMetadataStates &active_from,
    MetadataScores &next, ActiveMetadataStates &active_to,
    const MetadataStates &states, const ModelInputs &inputs, MetadataParents &parents
) {
    active_to.clear();
    const SiteContext site = site_context(inputs.seq, position);

    // Self-edge scores do not depend on metadata. IDs are grouped by coarse
    // state, so ascending `active_from` lets one evaluation per coarse state
    // serve all of its hypotheses; the precomputed destination still advances
    // each hypothesis' metadata correctly.
    TransitionResult self;
    int self_state = -1;

    for (const uint16_t from_id : active_from) {
        const MetadataState &source = states[from_id];
        const int from = source.state;
        if (from != self_state) {
            self = evaluate_transition(position, from, from, PathMetadata{}, inputs, site);
            self_state = from;
        }
        const auto retain = [&](int to_id, double score) {
            if (score > next[to_id]) {
                if (next[to_id] == -std::numeric_limits<double>::infinity())
                    active_to.push_back(static_cast<uint16_t>(to_id));
                next[to_id] = score;
                parents[to_id] = from_id;
            }
        };
        if (self.allowed) {
            retain(source.self_destination,
                   previous[from_id] + self.transition_score + self.emission_score);
        }
        const auto cross_edge = [&](int to) {
            const auto edge = evaluate_transition(position, from, to, source.metadata, inputs, site);
            if (!edge.allowed) return;
            const auto key = canonicalize_metadata(to, edge.next_metadata);
            retain(metadata_state_id(to, key),
                   previous[from_id] + edge.transition_score + edge.emission_score);
        };
        if (from == 0 && site.start) {
            for (int to = 1; to <= 3; ++to) cross_edge(to);
        } else if (is_exon(from)) {
            if (site.donor) cross_edge(from + 3);
            if (site.stop) cross_edge(0);
        } else if (is_intron(from) && site.acceptor) {
            for (int to = 1; to <= 3; ++to) cross_edge(to);
        }
    }
    std::sort(active_to.begin(), active_to.end());
}

// Dense convenience form: scans `previous` for its active entries, resets
// `next` and marks every unreached parent. The decoder uses the sparse form.
inline void advance_metadata_scores(
    int position, const MetadataScores &previous, MetadataScores &next,
    const MetadataStates &states, const ModelInputs &inputs,
    MetadataParents *parents = nullptr
) {
    ActiveMetadataStates active_from, active_to;
    for (int id = 0; id < METADATA_STATE_COUNT; ++id)
        if (std::isfinite(previous[id])) active_from.push_back(static_cast<uint16_t>(id));
    next.fill(-std::numeric_limits<double>::infinity());
    MetadataParents local_parents;
    MetadataParents &target = parents ? *parents : local_parents;
    target.fill(NO_METADATA_PARENT);
    advance_metadata_scores(position, previous, active_from, next, active_to, states, inputs, target);
}

} // namespace uniann
