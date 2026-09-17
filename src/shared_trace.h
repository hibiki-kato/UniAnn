#pragma once

#include <vector>
#include <functional>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <utility>

namespace uniann {

constexpr uint32_t NULL_NODE = 0xFFFFFFFF;

struct TraceNode {
    uint32_t parent = NULL_NODE;
    uint32_t first_child = NULL_NODE;
    uint32_t next_sibling = NULL_NODE;
    uint32_t prev_sibling = NULL_NODE;

    uint32_t head_refs = 0;
    uint32_t child_count = 0;

    uint32_t start_pos = 0;
    uint32_t end_pos = 0;
    uint8_t state = 0;
};

class SharedTrace {
private:
    std::vector<TraceNode> pool;
    uint32_t free_head = NULL_NODE;

    std::vector<uint32_t> slot_heads;
    std::vector<uint32_t> parent_usages_buffer;
    std::vector<uint32_t> next_heads_buffer;
    // Slots whose head is set, ascending. Every step touches only these and
    // the incoming active slots, never the whole slot range.
    std::vector<int> active_slots;
    std::vector<int> next_active_slots;
    uint32_t root = NULL_NODE;

    uint32_t last_position = 0;
    size_t live_node_count = 0;
    size_t peak_live_node_count = 0;
    size_t created_node_count = 0;
    size_t reused_extension_count = 0;

    std::function<void(uint8_t, uint32_t, uint32_t)> emit_run;

    uint32_t allocate_node() {
        if (live_node_count >= 0xFFFFFFFE) {
            throw std::runtime_error("Node pool exhausted (active nodes)");
        }
        uint32_t id;
        if (free_head != NULL_NODE) {
            id = free_head;
            free_head = pool[id].next_sibling;
        } else {
            if (pool.size() >= 0xFFFFFFFE) {
                throw std::runtime_error("Node pool capacity exhausted");
            }
            id = static_cast<uint32_t>(pool.size());
            pool.push_back(TraceNode{});
        }
        pool[id] = TraceNode{};
        live_node_count++;
        if (live_node_count > peak_live_node_count) {
            peak_live_node_count = live_node_count;
        }
        created_node_count++;
        return id;
    }

    void free_node(uint32_t id) {
        pool[id].next_sibling = free_head;
        free_head = id;
        live_node_count--;
    }

    bool can_extend(uint32_t parent, uint32_t uses, uint8_t state) const {
        const auto &node = pool[parent];
        return uses == 1 && node.head_refs == 1 && node.child_count == 0 && node.state == state;
    }

    // Iteratively prune dead leaves up the tree
    void prune_dead_branches(uint32_t id) {
        // A node can be freed if it's no longer an active head and has no children.
        // The root node is never pruned this way; it is only freed when flushed or finished.
        while (id != NULL_NODE && id != root && pool[id].head_refs == 0 && pool[id].child_count == 0) {
            uint32_t parent_id = pool[id].parent;

            uint32_t prev = pool[id].prev_sibling;
            uint32_t next = pool[id].next_sibling;

            // Unlink from siblings/parent
            if (prev != NULL_NODE) {
                pool[prev].next_sibling = next;
            } else if (parent_id != NULL_NODE) {
                pool[parent_id].first_child = next;
            }

            if (next != NULL_NODE) {
                pool[next].prev_sibling = prev;
            }

            if (parent_id != NULL_NODE) {
                pool[parent_id].child_count--;
            }

            free_node(id);
            id = parent_id;
        }
    }

public:
    SharedTrace(size_t slots, size_t initial_slot, std::function<void(uint8_t, uint32_t, uint32_t)> emit_run_)
        : slot_heads(slots, NULL_NODE), emit_run(std::move(emit_run_))
    {
        if (initial_slot >= slots) {
            throw std::invalid_argument("Invalid initial slot");
        }
        parent_usages_buffer.assign(slots, 0);
        next_heads_buffer.assign(slots, NULL_NODE);
        active_slots.reserve(slots);
        next_active_slots.reserve(slots);
        root = allocate_node();
        pool[root].state = 0;
        pool[root].start_pos = 0;
        pool[root].end_pos = 0;
        pool[root].head_refs = 1;
        slot_heads[initial_slot] = root;
        active_slots.push_back(static_cast<int>(initial_slot));
    }

    // Dense form: `parent_slots[to]` is -1 for an inactive destination.
    void advance(const std::vector<int>& parent_slots, const std::vector<uint8_t>& states, uint32_t position, bool compress_runs = true) {
        if (parent_slots.size() != slot_heads.size() || states.size() != slot_heads.size()) {
            throw std::invalid_argument("Frontier dimensions differ from slot count");
        }
        std::vector<int> active_to;
        for (size_t to = 0; to < parent_slots.size(); ++to) {
            if (parent_slots[to] >= 0) {
                active_to.push_back(static_cast<int>(to));
            } else if (parent_slots[to] != -1) {
                throw std::invalid_argument("Invalid negative parent slot referenced (must be -1 for inactive)");
            }
        }
        advance(active_to, parent_slots, states, position, compress_runs);
    }

    // Sparse form: only `parent_slots[to]` and `states[to]` for `to` in
    // `active_to` are read. `active_to` must be ascending and duplicate-free.
    void advance(const std::vector<int>& active_to, const std::vector<int>& parent_slots,
                 const std::vector<uint8_t>& states, uint32_t position, bool compress_runs = true) {
        if (root == NULL_NODE || position != last_position + 1) {
            throw std::invalid_argument("Trace positions must advance by one");
        }
        if (active_to.empty()) throw std::invalid_argument("Empty next frontier");
        // Validate all parents before changing the tree or its usage counts.
        for (int to : active_to) {
            if (to < 0 || to >= static_cast<int>(slot_heads.size()) ||
                to >= static_cast<int>(parent_slots.size()) || to >= static_cast<int>(states.size())) {
                throw std::invalid_argument("Active slot outside the frontier");
            }
            int p = parent_slots[to];
            if (p < 0 || p >= static_cast<int>(slot_heads.size()) || slot_heads[p] == NULL_NODE) {
                throw std::invalid_argument("Invalid positive or inactive parent slot referenced");
            }
        }
        for (int to : active_to) parent_usages_buffer[parent_slots[to]]++;

        // Create or extend nodes for the new frontier
        for (int to : active_to) {
            int p = parent_slots[to];
            uint32_t p_id = slot_heads[p];
            uint8_t s = states[to];

            // Extend only an unshared leaf. Otherwise another surviving path
            // could observe an end position belonging to this candidate alone.
            if (compress_runs && can_extend(p_id, parent_usages_buffer[p], s)) {
                pool[p_id].end_pos = position;
                pool[p_id].head_refs++; // Will be decremented below when releasing old heads
                next_heads_buffer[to] = p_id;
                reused_extension_count++;
            } else {
                uint32_t new_id = allocate_node();
                pool[new_id].state = s;
                pool[new_id].start_pos = position;
                pool[new_id].end_pos = position;
                pool[new_id].parent = p_id;
                pool[new_id].head_refs = 1;

                // Link to parent's children list
                pool[new_id].next_sibling = pool[p_id].first_child;
                if (pool[p_id].first_child != NULL_NODE) {
                    pool[pool[p_id].first_child].prev_sibling = new_id;
                }
                pool[p_id].first_child = new_id;
                pool[p_id].child_count++;

                next_heads_buffer[to] = new_id;
            }
        }

        // Acquire every next head before releasing old ones: a selected
        // parent must not disappear while another branch is being pruned.
        for (int old_slot : active_slots) {
            uint32_t old_id = slot_heads[old_slot];
            pool[old_id].head_refs--;
            prune_dead_branches(old_id);
            parent_usages_buffer[old_slot] = 0;
        }

        slot_heads.swap(next_heads_buffer);
        // The buffer now holds the released heads; clear them so it is empty
        // again for the next step.
        for (int old_slot : active_slots) next_heads_buffer[old_slot] = NULL_NODE;
        active_slots.assign(active_to.begin(), active_to.end());
        last_position = position;
    }

    // Flush unique prefix of the trace
    void flush_common() {
        while (root != NULL_NODE && pool[root].head_refs == 0 && pool[root].child_count == 1) {
            emit_run(pool[root].state, pool[root].start_pos, pool[root].end_pos);

            uint32_t child = pool[root].first_child;
            pool[child].parent = NULL_NODE;
            pool[child].prev_sibling = NULL_NODE;
            pool[child].next_sibling = NULL_NODE;

            free_node(root);
            root = child;
        }
    }

    // Finish trace reconstruction for a chosen winning slot
    void finish(size_t winning_slot) {
        if (winning_slot >= slot_heads.size() || slot_heads[winning_slot] == NULL_NODE) {
            throw std::invalid_argument("Invalid or inactive winning slot");
        }

        // Prune all other branches
        for (int slot : active_slots) {
            if (static_cast<size_t>(slot) != winning_slot) {
                uint32_t old_id = slot_heads[slot];
                pool[old_id].head_refs--;
                prune_dead_branches(old_id);
                slot_heads[slot] = NULL_NODE;
            }
        }
        active_slots.clear();

        // Now there should be a single unbranched path from root to the winning leaf
        flush_common();

        // Flush the winning leaf
        if (root != NULL_NODE) {
            emit_run(pool[root].state, pool[root].start_pos, pool[root].end_pos);
            pool[root].head_refs = 0; // Release root's remaining reference
            free_node(root);
            root = NULL_NODE;
        }

        slot_heads[winning_slot] = NULL_NODE;
    }

    size_t live_nodes() const { return live_node_count; }
    size_t peak_live_nodes() const { return peak_live_node_count; }
    size_t pool_slots() const { return pool.size(); }
    size_t pool_capacity() const { return pool.capacity(); }
    size_t bytes_per_node() const { return sizeof(TraceNode); }
    size_t created_nodes() const { return created_node_count; }
    size_t reused_extensions() const { return reused_extension_count; }
};

} // namespace uniann
