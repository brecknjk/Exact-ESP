#ifndef EXTRACTOR_SPECIAL_GRAPH_CONSTRUCTOR_HPP
#define EXTRACTOR_SPECIAL_GRAPH_CONSTRUCTOR_HPP

#include "../utils/types.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

// Construct a simple "special" graph:
//
// Params:
// - coords: output vector of node coordinates (cleared and filled)
// - edges:  output vector of edges (cleared and filled)
// - n:      number of nodes to create (must be >= 1)
inline void constructWorstCasePolyAnyaGraph(std::vector<Coordinates>& coords, std::vector<Edge>& edges, uint32_t n){
    if (n == 0) throw std::invalid_argument("n must be >= 1");
    assert(coords.size() == 0);
    assert(edges.size() == 0);

    long width = 2 + n + n + n + n + 3;
    long height = 3 + n + n + 1 + 1 + 2;
    uint32_t connectFrom = 0;
    uint32_t numberComponents = 0;

    auto _insertConnections = [&]() { 
        for(long i = connectFrom; i < coords.size() - 1; i++){
            edges.emplace_back(i, i + 1);
        }
        edges.emplace_back(coords.size() - 1, connectFrom);
        connectFrom = coords.size();
        numberComponents++;
    };

    // Box
    coords.emplace_back(width, 0);
    coords.emplace_back(width, height);
    coords.emplace_back(0, height);
    coords.emplace_back(0, 0);
    coords.emplace_back(2, 0);
    coords.emplace_back(2, height - 2);
    coords.emplace_back(width - 2, height - 2);
    coords.emplace_back(width - 2, 2);
    coords.emplace_back(width - 3, 2);
    coords.emplace_back(width - 3, height - 3);
    for(long i = width - n - 3; i >= width - 3 - n - n; i--){
        coords.emplace_back(i, height - 3);
    }
    coords.emplace_back(3, height - 3);
    coords.emplace_back(3, 0);
    for(long i = width - n - 3 - n; i <= width - n - 3; i++){
        coords.emplace_back(i, 0);
    }
    _insertConnections();

    //Small Box in corner to confuse polyanyas skipahead
    coords.emplace_back(5, height - 0.5);
    coords.emplace_back(5, height - 1.5);
    coords.emplace_back(6, height - 1.5);
    coords.emplace_back(6, height - 0.5);
    _insertConnections();

    //Line
    coords.emplace_back(5, 2);
    coords.emplace_back(5, 2 + 1);
    coords.emplace_back(n, 2 + 1);
    coords.emplace_back(n, 2);
    _insertConnections();

    //Obsticals left
    for(long i = 4; i < height - 4; i += 2){
        coords.emplace_back(5 - 0.5, i);
        coords.emplace_back(5 + 0.5, i + 0.5);
        coords.emplace_back(5 - 0.5, i + 1);
        _insertConnections();
    }

    //Obsticals right
    for(long i = 4; i < height - 4; i += 2){
        coords.emplace_back(n + 0.5, i);
        coords.emplace_back(n - 0.5, i + 0.5);
        coords.emplace_back(n + 0.5, i + 1);
        _insertConnections();
    }

    if(settings.VERBOSE_LEVEL > 0){
        printf("size domain = %li x %li\n", width, height);
        printf("generated: %i obsticals, %li nodes, %li edges\n", numberComponents, coords.size(), edges.size());
    }
    assert(edges.size() == coords.size());

    if(settings.DEBUG_LEVEL > 0){
        for(auto e : edges){
            assert(e.source >= 0 && e.source < coords.size());
            assert(e.target >= 0 && e.target < coords.size());
        }
    }
}

#endif // EXTRACTOR_SPECIAL_GRAPH_CONSTRUCTOR_HPP