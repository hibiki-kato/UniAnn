#include "../src/shared_trace.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <random>

using namespace uniann;

struct Run {
    uint8_t state;
    uint32_t start_pos;
    uint32_t end_pos;
    bool operator==(const Run& other) const {
        return state == other.state && start_pos == other.start_pos && end_pos == other.end_pos;
    }
};

std::vector<Run> merge_runs(const std::vector<Run>& runs) {
    std::vector<Run> merged;
    for (const auto& r : runs) {
        if (!merged.empty() && merged.back().state == r.state && merged.back().end_pos + 1 == r.start_pos) {
            merged.back().end_pos = r.end_pos;
        } else {
            merged.push_back(r);
        }
    }
    return merged;
}

void test_basic() {
    std::vector<Run> emitted;
    auto emit_run = [&](uint8_t s, uint32_t start, uint32_t end) {
        emitted.push_back({s, start, end});
    };

    SharedTrace tree(3, 0, emit_run);

    tree.advance({0, -1, -1}, {1, 0, 0}, 1);
    tree.advance({0, 0, -1}, {1, 2, 0}, 2);
    tree.advance({0, 1, -1}, {1, 2, 0}, 3);

    tree.flush_common();

    auto merged = merge_runs(emitted);
    assert(merged.size() == 1);
    assert(merged[0] == (Run{0, 0, 0}));

    tree.finish(1);

    merged = merge_runs(emitted);
    assert(merged.size() == 3);
    assert(merged[0] == (Run{0, 0, 0}));
    assert(merged[1] == (Run{1, 1, 1}));
    assert(merged[2] == (Run{2, 2, 3}));

    assert(tree.live_nodes() == 0);

    std::cout << "test_basic passed\n";
}

void test_rle_alias_safety() {
    std::vector<Run> emitted;
    auto emit_run = [&](uint8_t s, uint32_t start, uint32_t end) {
        emitted.push_back({s, start, end});
    };

    SharedTrace tree(3, 0, emit_run);
    tree.advance({0, 0, -1}, {1, 1, 0}, 1);
    tree.advance({0, 1, -1}, {1, 1, 0}, 2);

    tree.finish(0);

    auto merged = merge_runs(emitted);
    assert(merged.size() == 2);
    assert(merged[0] == (Run{0, 0, 0}));
    assert(merged[1] == (Run{1, 1, 2}));
    assert(tree.live_nodes() == 0);
    std::cout << "test_rle_alias_safety passed\n";
}

void test_invalid_arguments() {
    std::vector<Run> emitted;
    auto emit_run = [&](uint8_t s, uint32_t start, uint32_t end) {
        emitted.push_back({s, start, end});
    };
    SharedTrace tree(3, 0, emit_run);

    // Test invalid negative argument (not -1)
    bool threw = false;
    try {
        tree.advance({0, -2, -1}, {1, 0, 0}, 1);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    // Test invalid positive argument
    threw = false;
    try {
        tree.advance({0, 3, -1}, {1, 0, 0}, 1);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    // Test inactive positive parent
    tree.advance({0, -1, -1}, {1, 0, 0}, 1);
    threw = false;
    try {
        tree.advance({0, 1, -1}, {1, 0, 0}, 2); // 1 is inactive
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    // Dimensions mismatch
    threw = false;
    try {
        tree.advance({0, -1}, {1, 0}, 2);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    // Invalid winning slot
    threw = false;
    try {
        tree.finish(2);
    } catch (const std::invalid_argument&) { threw = true; }
    assert(threw);

    std::cout << "test_invalid_arguments passed\n";
}

void test_randomized(bool compress) {
    size_t slots = 10;
    size_t steps = 1000;
    std::mt19937 rng(42 + compress);

    std::vector<Run> emitted;
    auto emit_run = [&](uint8_t s, uint32_t start, uint32_t end) {
        emitted.push_back({s, start, end});
    };

    SharedTrace tree(slots, 0, emit_run);

    std::vector<std::vector<uint8_t>> oracle(slots);
    oracle[0].push_back(0);

    for (uint32_t pos = 1; pos <= steps; ++pos) {
        std::vector<int> parent_slots(slots, -1);
        std::vector<uint8_t> states(slots, 0);

        std::vector<int> active;
        for (size_t i = 0; i < slots; ++i) {
            if (!oracle[i].empty()) active.push_back(i);
        }
        if (active.empty()) break;

        std::vector<std::vector<uint8_t>> next_oracle(slots);

        for (size_t to = 0; to < slots; ++to) {
            if (rng() % 100 < 20 && to != 0) {
                continue;
            }
            int p = active[rng() % active.size()];
            uint8_t state = rng() % 4;

            parent_slots[to] = p;
            states[to] = state;

            next_oracle[to] = oracle[p];
            next_oracle[to].push_back(state);
        }

        tree.advance(parent_slots, states, pos, compress);
        oracle = std::move(next_oracle);

        if (pos % 10 == 0) tree.flush_common();
    }

    std::vector<int> active;
    for (size_t i = 0; i < slots; ++i) {
        if (!oracle[i].empty()) active.push_back(i);
    }
    int winner = active[rng() % active.size()];

    tree.finish(winner);

    auto merged = merge_runs(emitted);

    std::vector<Run> expected;
    for (size_t i = 0; i < oracle[winner].size(); ++i) {
        uint8_t s = oracle[winner][i];
        if (!expected.empty() && expected.back().state == s && expected.back().end_pos + 1 == i) {
            expected.back().end_pos = i;
        } else {
            expected.push_back({s, (uint32_t)i, (uint32_t)i});
        }
    }

    if (merged.size() != expected.size()) {
        std::cerr << "Mismatch sizes: merged " << merged.size() << " expected " << expected.size() << "\n";
    }
    assert(merged.size() == expected.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        if (!(merged[i] == expected[i])) {
            std::cerr << "Mismatch at " << i << ": merged " << (int)merged[i].state << " [" << merged[i].start_pos << "," << merged[i].end_pos << "] expected " << (int)expected[i].state << " [" << expected[i].start_pos << "," << expected[i].end_pos << "]\n";
        }
        assert(merged[i] == expected[i]);
    }

    assert(tree.live_nodes() == 0);

    std::cout << "test_randomized(compress=" << compress << ") passed\n";
    std::cout << "  Stats: peak_live=" << tree.peak_live_nodes()
              << ", created=" << tree.created_nodes()
              << ", reused=" << tree.reused_extensions() << "\n";
}


void test_long_same_state_run() {
    std::vector<Run> emitted;
    SharedTrace tree(1, 0, [&](uint8_t state, uint32_t start, uint32_t end) {
        emitted.push_back({state, start, end});
    });
    for (uint32_t position = 1; position <= 100000; ++position) {
        tree.advance({0}, {0}, position);
        if (position % 4096 == 0) tree.flush_common();
    }
    assert(tree.peak_live_nodes() == 1);
    tree.finish(0);
    assert(emitted.size() == 1);
    assert(emitted[0] == (Run{0, 0, 100000}));
    assert(tree.live_nodes() == 0);
}

int main() {
    test_long_same_state_run();
    test_basic();
    test_rle_alias_safety();
    test_invalid_arguments();
    test_randomized(false);
    test_randomized(true);
    std::cout << "All tests passed successfully.\n";
    return 0;
}
