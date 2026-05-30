#ifndef union_find_h
#define union_find_h

#include <cassert>
#include "../utils/types.hpp"

// see https://github.com/kartikkukreja/blog-codes/blob/master/src/Union%20Find%20%28Disjoint%20Set%29%20Data%20Structure.cpp
// template<std::size_t S>
class UnionFind {
private:
    uint32_t set_count;
    std::vector<uint32_t> set_size, parent;

public:
// Create an empty union find data structure with initial_count isolated sets.
UnionFind(uint32_t initial_count) {
    set_count = initial_count;

    set_size.resize(initial_count, 1);

    parent.resize(initial_count);
    for (uint32_t i = 0; i < initial_count; i++){
        parent[i] = i;
    }
}
~UnionFind() {
    set_size.clear();
    parent.clear();
}

// Return the parent of component corresponding to object p.
uint32_t find_root(uint32_t x) {
    uint32_t root = x;
    while (root != parent[root]){
        root = parent[root];
    }
    // set new root/shortcuts
    while (x != root) { 
        uint32_t old_parent = parent[x]; 
        parent[x] = root; 
        x = old_parent; 
    }
    return root;
}

// Replace sets containing x and y with their union.
void merge(uint32_t x, uint32_t y) {
    uint32_t root_x = find_root(x);
    uint32_t root_y = find_root(y);
    if(root_x == root_y) return;
    if(set_size[root_x] > set_size[root_y]) std::swap(root_x, root_y);
    parent[root_x] = root_y;
    set_size[root_y] += set_size[root_x];
    set_count--;
}

// Are objects x and y in the same set?
bool is_connected(uint32_t x, uint32_t y) { 
    return find_root(x) == find_root(y); 
}

// Return the number of disjoint sets.
uint32_t get_set_count() { return set_count; }

uint32_t get_count_of_set(uint32_t x) { return set_size[x]; }
uint32_t operator[](uint32_t x) { return set_size[x]; }

void set_count_of_set(uint32_t x, uint32_t new_set_size) { set_size[x] = new_set_size; }
};

#endif