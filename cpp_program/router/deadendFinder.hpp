#ifndef DeadendFinder_h
#define DeadendFinder_h

#include "../utils/types.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"
#include "../utils/unionFind.hpp"

#include "graph.hpp"

class DeadendFinder{
    std::vector<uint32_t> edge_to_exit;
    std::vector<bool> edge_is_deadend;
    uint32_t edgeOffset;
    uint32_t last_lcp, last_first_start, last_first_goal;
    bool isSet = false, isInit = false;

    inline bool edgeIsSaved(uint32_t e) const {
        return (e >> 1) - edgeOffset >= 0 && (e >> 1) - edgeOffset < edge_to_exit.size();
    };
    inline uint32_t map_edge(uint32_t e) const {
        assert(edgeIsSaved(e));
        return (e >> 1) - edgeOffset;
    };
    inline bool edgeIsDeadend(uint32_t e) const {
        return (e >> 1) < edgeOffset ? true : (e >> 1) - edgeOffset >= edge_to_exit.size() ? false : edge_is_deadend[map_edge(e)];
    };
    uint32_t rev_map_edge(uint32_t e) const {return (e + edgeOffset) << 1;};

    uint32_t getFirstDeadEndEdge(uint32_t e) const {
        return edge_is_deadend[map_edge(e)] ? map_edge(e) : edge_to_exit[map_edge(e)];
    }

    // starting from start and goal use edge_to_exit to set edge_is_deadend to false until sea
    uint32_t giveLCP_noReset(uint32_t startEdge, uint32_t goalEdge, bool mark_deadEnd = false){
        while(true){
            if(startEdge == goalEdge) return startEdge;
            if(startEdge == edge_to_exit[startEdge] && goalEdge == edge_to_exit[goalEdge]) return UINT32_MAX;
            if(startEdge != edge_to_exit[startEdge]) {
                if(edge_is_deadend[startEdge] == mark_deadEnd) return startEdge;
                edge_is_deadend[startEdge] = mark_deadEnd;
                startEdge = edge_to_exit[startEdge];
            }
            if(goalEdge != edge_to_exit[goalEdge]) {
                if(edge_is_deadend[goalEdge] == mark_deadEnd) return goalEdge;
                edge_is_deadend[goalEdge] = mark_deadEnd;
                goalEdge = edge_to_exit[goalEdge];
            }
        }
    }

    // starting from start and goal use edge_to_exit to set edge_is_deadend to false until sea
    uint32_t lcp(uint32_t startEdge, uint32_t goalEdge){
        auto lcp = giveLCP_noReset(startEdge, goalEdge, false);
        giveLCP_noReset(startEdge, goalEdge, true); // revert
        return lcp;
    }

    // starting from start and goal use edge_to_exit to set edge_is_deadend to false until sea
    void markPath(uint32_t start_index, uint32_t lcp, bool mark_deadEnd = false){
        uint32_t current_index = start_index;
        while(true){
            assert(current_index != UINT32_MAX);
            assert(edge_is_deadend[current_index] != mark_deadEnd || current_index == lcp);
            edge_is_deadend[current_index] = mark_deadEnd;
            if(current_index == lcp) break;
            if(current_index == edge_to_exit[current_index]) break;
            current_index = edge_to_exit[current_index];
        }
    }

public:
    // inline const bool isInititalized() const {return isInit;}
    inline const bool isDeadend(uint32_t e) const { return isInit ? edgeIsDeadend(e) : false;}

    DeadendFinder(){};
	DeadendFinder(Graph* g, uint32_t edge_in_sea){
        isInit = true;
        edgeOffset = g->constrainedEdgeDivider;
        assert(edge_in_sea != UINT32_MAX);
        assert(!g->edges[edge_in_sea >> 1].isConstrained());
        assert(edge_in_sea >> 1 >= edgeOffset);
        edge_to_exit.resize(g->edges.size() - edgeOffset, UINT32_MAX);

        {
            auto uf = UnionFind(g->nodeCount());
            for(uint32_t i = 0; i < edgeOffset; i++){
                assert(g->edges[i].isConstrained());
                assert(g->edges[i].N[0] != UINT32_MAX && g->edges[i].N[1] != UINT32_MAX);
                uf.merge(g->edges[i].N[0], g->edges[i].N[1]);
            }
            edge_is_deadend.resize(g->edges.size() - edgeOffset);
            for(uint32_t e = edgeOffset; e < g->edges.size(); e++){
                assert(!g->edges[e].isConstrained());
                assert(g->edges[e].N[0] != UINT32_MAX && g->edges[e].N[1] != UINT32_MAX);
                edge_is_deadend[map_edge(e << 1)] = uf.is_connected(g->edges[e].N[0], g->edges[e].N[1]);
            }
        }
        
        edge_to_exit[map_edge(edge_in_sea)] = map_edge(edge_in_sea);
        std::list<std::pair<uint32_t, uint32_t>> queue;
        queue.emplace_back(edge_in_sea, map_edge(edge_in_sea));
        queue.emplace_back(edge_in_sea ^ 1, map_edge(edge_in_sea));
        uint32_t e = edgeOffset << 1; //used as iterator to find new components when queue is empty
        uint32_t debug_components = 1; //number of components, if all is connected == 1
        assert(edgeIsSaved(edge_in_sea));
        
        while(!queue.empty()){
            auto current_edge = queue.front().first; 
            auto lastExit = queue.front().second; queue.pop_front();
            assert(lastExit >> 1 < g->edges.size());
            assert(current_edge != UINT32_MAX);
            assert(g->nextEdge(current_edge) != UINT32_MAX);
            if(settings.VERBOSE_LEVEL > 2){
                printf("DeadendFinder: processing edge %u (%u %u), exit edge %u (%u %u)\n", 
                    current_edge, g->firstNode(current_edge), g->secondNode(current_edge),
                    lastExit, g->firstNode(rev_map_edge(lastExit)), g->secondNode(rev_map_edge(lastExit)));
            }

            for(auto f_iter = FaceCirculator(g, current_edge); f_iter.hasNext(); f_iter++){
                if(g->edges[f_iter.current_e >> 1].isConstrained()) continue; // skip constrained edges
                uint32_t current_index = map_edge(f_iter.current_e);
                if(edge_to_exit[current_index] != UINT32_MAX) continue; // deduplication
                //assert(edge_is_deadend[current_index]); 
                assert(lastExit != current_index); // don't go back
                
                edge_to_exit[current_index] = lastExit;
                queue.emplace_back(f_iter.current_e ^ 1, edge_is_deadend[current_index] ? current_index : lastExit);

                if(settings.VERBOSE_LEVEL > 2){
                    printf("  setting edge %u (%u %u) to exit via edge %u (%u %u), deadend=%s\n", 
                        rev_map_edge(current_index), g->firstNode(rev_map_edge(current_index)), g->secondNode(rev_map_edge(current_index)),
                        rev_map_edge(edge_to_exit[current_index]),
                        g->firstNode(rev_map_edge(edge_to_exit[current_index])),
                        g->secondNode(rev_map_edge(edge_to_exit[current_index])),
                        edge_is_deadend[current_index] ? "true" : "false");
                }
            }
            while(e < g->edges.size() * 2 && queue.empty()){
                if(edge_to_exit[map_edge(e)] == UINT32_MAX){
                    edge_to_exit[map_edge(e)] = map_edge(e);
                    queue.emplace_back(e, map_edge(e));
                    queue.emplace_back(e ^ 1, map_edge(e));
                    debug_components++;
                }
                e += 2;
            }
        }
        if(settings.DEBUG_LEVEL > 0){
            // counts the count and the maximal size of each freespace region
            auto countEdgesInFreespace = std::vector<uint32_t>(edge_to_exit.size(), 0);
            uint32_t countDeadendEdges = 0, countFreeSpaceEdges = 0;
            for(uint32_t i = 0; i < edge_to_exit.size(); i++){
                assert(edge_to_exit[i] != UINT32_MAX);
                if(edge_is_deadend[i]){
                    countDeadendEdges++;
                    continue;
                }
                countFreeSpaceEdges++;
                countEdgesInFreespace[edge_to_exit[i]]++;
            }
            uint32_t max = *std::max_element(countEdgesInFreespace.begin(), countEdgesInFreespace.end());
            uint32_t freespacesCount = 0;
            for(auto item : countEdgesInFreespace){
                if(item > 0) freespacesCount++;
            }
            printf("DeadendFinder: found %u components\n", debug_components);
            printf("DeadendFinder: #main freespace: %u\n", max);
            printf("DeadendFinder: #main freespace / #total edges: %f\n", max / (float)g->edges.size());
            printf("DeadendFinder: count off free spaces: %u\n", freespacesCount);
            printf("#TotalEdges: %u\n", (uint32_t) g->edges.size());
            printf("#ObsticalEdges: %u\n", edgeOffset);
            printf("#DeadendEdges: %u\n", countDeadendEdges);
            printf("#FreeSpaceEdges: %u\n", countFreeSpaceEdges);
        }else if(settings.VERBOSE_LEVEL > 1) {
            printf("DeadendFinder: found %u components\n", debug_components);
        }
    };

    bool markBeforeSearch(uint32_t nearestEdgeStart, uint32_t nearestEdgeGoal){
        if(nearestEdgeStart == UINT32_MAX || nearestEdgeGoal == UINT32_MAX) return false;
        if(nearestEdgeStart == nearestEdgeGoal) return true;
        assert(!isSet);

        last_first_start = getFirstDeadEndEdge(nearestEdgeStart);
        last_first_goal = getFirstDeadEndEdge(nearestEdgeGoal);
        last_lcp = lcp(last_first_start, last_first_goal);
        if(last_lcp == UINT32_MAX) return false;
        markPath(last_first_start, last_lcp, false);
        markPath(last_first_goal, last_lcp, false);
        isSet = true;
        return true;
    }

    void resetAfterSearch(){
        if(!isSet) return; // no need for reset
        markPath(last_first_start, last_lcp, true);
        markPath(last_first_goal, last_lcp, true);
        isSet = false;
    }
    
};
#endif