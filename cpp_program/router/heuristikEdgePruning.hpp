#ifndef heuristik_edgePruning_h
#define heuristik_edgePruning_h

#include "../utils/types.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

#include "graph.hpp"

constexpr double PRUNING_EPSILON_HEURISTIK = 1.00000001;

class HeuristikEdgePruning{
private:
    Graph* g = nullptr;
    uint32_t edgeOffset = 0;
    std::vector<std::pair<NT_Distance, uint32_t>> bestWavelets_viaEdge;

    uint32_t const inline map_edge(uint32_t e) const { return e - (edgeOffset << 1); };
    uint32_t const inline rev_map_edge(uint32_t e) const { return e + (edgeOffset << 1); };
public:
    HeuristikEdgePruning(){}
    HeuristikEdgePruning(Graph* g, uint32_t edgeCount, uint32_t edgeOffset) : g(g), edgeOffset(edgeOffset){
        assert(settings.polyanya.heuristikEdgePruningProof);
        assert(!settings.polyanya.merge_cells); // assume triangles
        bestWavelets_viaEdge.resize(((edgeCount - edgeOffset) << 1) + 12, {MAX_DOUBLE, UINT32_MAX});
        if(settings.VERBOSE_LEVEL > 1) printf("HeuristikEdgePruning initialized with %i edges and offset %i\n", edgeCount - edgeOffset, edgeOffset);
    }
    bool isSet(){return g != nullptr && bestWavelets_viaEdge.size() > 0;};

    void resetAfterSearch(std::list<uint32_t> & visited){
        if(!isSet()) return;
        for(auto & e : visited) {
            if(g->isEdgeConstrained(e << 1)) continue;
            bestWavelets_viaEdge[map_edge(e << 1)] = {MAX_DOUBLE, UINT32_MAX};
            bestWavelets_viaEdge[map_edge(e << 1) ^ 1] = {MAX_DOUBLE, UINT32_MAX};
        }
        // test proper reset
        if(settings.DEBUG_LEVEL > 2){
            for(uint32_t i = 0; i < bestWavelets_viaEdge.size(); i++){
                assert(bestWavelets_viaEdge[i].first == MAX_DOUBLE);
            }
        }
    }

    void insertCandidate(uint32_t edge_id, NT_Distance & dist, uint32_t generator){
        if(!isSet()) return;
        if(generator == 7){
            auto temp = 0;
        }
        if(bestWavelets_viaEdge[map_edge(edge_id)].first > dist){
            bestWavelets_viaEdge[map_edge(edge_id)] = {dist, generator};
        }
    }

    bool isDominated(const Wavelet & w, NT_Distance & costs_secondNode, NT_Distance & costs_firstNode){
        uint32_t n0 = g->firstNode(w.item), n1 = g->secondNode(w.item);
        uint32_t edge_n0 = g->nextEdge(g->otherEdge(w.item));
        uint32_t edge_n1 = g->nextEdge(edge_n0);
        assert(g->hasNode(edge_n0, n0) && g->hasNode(edge_n1, n1));
        assert(!g->isEdgeConstrained(w.item));

        if(!g->isEdgeConstrained(edge_n0)){
            auto & entry_n0 = bestWavelets_viaEdge[map_edge(edge_n0)];
            if(entry_n0.second != UINT32_MAX && entry_n0.second != w.generator && entry_n0.first * PRUNING_EPSILON_HEURISTIK < costs_secondNode){
                if(g->rightTurn(w.generator, entry_n0.second, n1) && (w.right_b == n1 || g->rightTurn(w.generator, entry_n0.second, n0))){
                    return true;
                }
            }
        }
        if(!g->isEdgeConstrained(edge_n1)){
            auto & entry_n1 = bestWavelets_viaEdge[map_edge(edge_n1)];
            if(entry_n1.second != UINT32_MAX && entry_n1.second != w.generator && entry_n1.first * PRUNING_EPSILON_HEURISTIK < costs_firstNode){
                if(g->leftTurn(w.generator, entry_n1.second, n0) && (w.left_b == n0 || g->leftTurn(w.generator, entry_n1.second, n1))){
                    return true;
                }
            }
        }
        return false;
    }
};

#endif