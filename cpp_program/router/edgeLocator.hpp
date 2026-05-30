#ifndef locator_h
#define locator_h

#include <stack>

#include "../utils/types.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

#include "graph.hpp"

class Locator{
    struct BBox{
        Coordinates min, max;
        uint32_t item = UINT32_MAX;

        static uint8_t getDim() { return 4; };
        bool isVaid() const { return item != UINT32_MAX; }
        
        void merge(const BBox & b0, const BBox & b1){
            min.lat = std::min(b0.min.lat, b1.min.lat);
            min.lon = std::min(b0.min.lon, b1.min.lon);
            max.lat = std::max(b0.max.lat, b1.max.lat);
            max.lon = std::max(b0.max.lon, b1.max.lon);
            item = std::min(b0.item, b1.item);
        };
        BBox() {};
        BBox(const Coordinates & min, const Coordinates & max, uint32_t item) : min(min), max(max), item(item){};
        
        inline bool is_smaller_in_dim(const BBox & b, uint32_t dim) const {
            assert(dim < 4);
            switch (dim){
            case 0:
                return min.lat < b.min.lat;
            case 1:
                return min.lon < b.min.lon;
            case 2:
                return max.lat < b.max.lat;
            case 3:
                return max.lon < b.max.lon;
            default:
                throw std::invalid_argument("Invalid dimension");
            }
        };
    };

    template<class T>
    class KDTreeInplace {
    public:
        uint32_t index_firstLeaf = UINT32_MAX;
        uint8_t maxDim = 0;
        std::vector<T> nodes;
        
        bool inline isLeaf(uint32_t index){return index >= index_firstLeaf;};

        const inline uint32_t size(){return nodes.size();}
        T& operator[](uint32_t index) {return nodes[index];};


        void buildTree(const std::vector<T> & items_to_store){
            index_firstLeaf = items_to_store.size() - 1;
            nodes.resize(items_to_store.size() + index_firstLeaf, T());
            maxDim = T::getDim();

            std::vector<std::vector<uint32_t>> sorted_ids_per_dim(maxDim, std::vector<uint32_t>(items_to_store.size()));
            #pragma omp parallel for
            for(uint8_t d = 0; d < maxDim; d++){
                auto comparator = [&items_to_store, d](uint32_t a, uint32_t b) {
                    return items_to_store[a].is_smaller_in_dim(items_to_store[b], d);
                };
                
                std::iota(sorted_ids_per_dim[d].begin(), sorted_ids_per_dim[d].end(), 0);
                std::sort(sorted_ids_per_dim[d].begin(), sorted_ids_per_dim[d].end(), comparator);
                if(settings.VERBOSE_LEVEL >= 1) printf("sorted in dim %i\n", d);
            }
            
            //id, depth -> calc size, dim
            std::stack<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>> build_stack;
            std::vector<bool> lookUp_order = std::vector<bool>(items_to_store.size());
            build_stack.push({0, {0, items_to_store.size() - 1}});
            while(!build_stack.empty()){
                uint32_t id = build_stack.top().first;
                uint8_t dim = int_log2(id + 1) % maxDim;
                auto range = build_stack.top().second;
                uint32_t itemCount = range.second - range.first + 1;
                
                //leaf case
                if(itemCount == 1){
                    nodes[id] = items_to_store[sorted_ids_per_dim[dim][range.first]];
                    build_stack.pop();
                    continue;
                }
                if(nodes[2*id + 1].isVaid() && nodes[2*id + 2].isVaid()){
                    nodes[id].merge(nodes[2*id + 1], nodes[2*id + 2]);
                    build_stack.pop();
                    continue;
                }
                //sort other dimensions
                auto buffer = std::vector<uint32_t>(itemCount);
                //create lookup for sorting
                
                uint32_t pivot;
                if(itemCount >> (int_log2(itemCount) - 1) == 3){
                    pivot = range.first + (1 << int_log2(itemCount));
                }else{
                    pivot = range.first + itemCount - (1 << (int_log2(itemCount) - 1));
                }
                for(uint32_t i = range.first; i <= range.second; i++){
                    lookUp_order[sorted_ids_per_dim[dim][i]] = i < pivot;
                }
                for(uint8_t d = 0; d < maxDim; d++){
                    if(d == dim) continue;
                    // copy into buffer
                    for(uint32_t i = 0; i < itemCount; i++){
                        buffer[i] = sorted_ids_per_dim[d][i + range.first];
                    }
                    //sort
                    uint32_t start_l = range.first, start_r = pivot;
                    for(uint32_t i = 0; i < itemCount; i++){
                        assert(start_r <= sorted_ids_per_dim[d].size());
                        sorted_ids_per_dim[d][(lookUp_order[buffer[i]]) ? start_l++ : start_r++] = buffer[i];
                    }
                }
                
                build_stack.push({2*id + 1, {range.first, pivot - 1}});
                build_stack.push({2*id + 2, {pivot, range.second}});
            }
        }

        const inline bool hit_in_dim(const uint32_t id, const T & min_search, const T & max_search, const uint8_t dim) const{
            return !nodes[id].is_smaller_in_dim(min_search, dim) && nodes[id].is_smaller_in_dim(max_search, dim);
        }
        const inline bool hit_in_all_dim(const uint32_t id, const T & min_search, const T & max_search) const{
            for(uint8_t dim = 0; dim < maxDim; dim++){
                if(!hit_in_dim(id, min_search, max_search, dim)) return false;
            }
            return true;
        }

        std::vector<uint32_t> findItems(T & min_search, T & max_search){
            std::vector<uint32_t> found;

            std::stack<uint32_t> search_stack;
            search_stack.push(0);
            while(!search_stack.empty()){
                uint32_t id = search_stack.top();
                uint8_t dim = int_log2(id + 1) % maxDim;
                search_stack.pop();

                //leaf case
                if(isLeaf(id)){
                    if(hit_in_all_dim(id, min_search, max_search)){
                        found.push_back(nodes[id].item);
                        if(settings.VERBOSE_LEVEL > 3) printf("found candidate: min(%lf,%lf) max(%lf,%lf)\n", nodes[id].min.lat, nodes[id].min.lon, nodes[id].max.lat, nodes[id].max.lon);
                    } 
                }else{
                    if(hit_in_dim(2*id + 1, min_search, max_search, dim)) search_stack.push(2*id + 1);
                    if(hit_in_dim(2*id + 2, min_search, max_search, dim)) search_stack.push(2*id + 2);
                }
            }
            return found;
        }
    };
	KDTreeInplace<BBox> cells;
    Graph* g;
    bool isSet = false;

	std::vector<uint32_t> getCellContent(const Coordinates & c){
		auto min_search = BBox(Coordinates(-MAX_FLOAT, -MAX_FLOAT), c, UINT32_MAX);
		auto max_search = BBox(c, Coordinates(MAX_FLOAT, MAX_FLOAT), UINT32_MAX);
		return cells.findItems(min_search, max_search);
	};

public:
    Locator(){};
    Locator(Graph* g) : g(g){
        isSet = true;
        assert(g->nodeCount() > 0);
        std::vector<bool> alreadyProcessed(g->edges.size() * 2, false);
        std::vector<BBox> bboxes;
        for(uint32_t i = 0; i < g->edges.size() * 2; i++){
            if(alreadyProcessed[i]) continue;
            assert(g->edges[i >> 1].valid());
            if(g->edges[i >> 1].E[i & 1] == UINT32_MAX) continue; // skip one direction of constraind Edges

            Coordinates min(MAX_FLOAT, MAX_FLOAT);
            Coordinates max(-MAX_FLOAT, -MAX_FLOAT);

            for(FaceCirculator f_iter(g, i); f_iter.hasNext(); f_iter++){
                auto node = f_iter.secondNode();
                if(g->coords[node].x() < min.lat) min.lat = to_double(g->coords[node].x());
                if(g->coords[node].x() > max.lat) max.lat = to_double(g->coords[node].x());
                if(g->coords[node].y() < min.lon) min.lon = to_double(g->coords[node].y());
                if(g->coords[node].y() > max.lon) max.lon = to_double(g->coords[node].y());
                alreadyProcessed[f_iter.current_e] = true;
            }
            bboxes.emplace_back(min, max, i);
        }
        alreadyProcessed.clear();
        cells.buildTree(bboxes);
    };
    bool isLocatorSet(){return isSet;};

	const uint32_t getNearestEdge(const Point & p){
        uint32_t debug_countOpsCGAL = 0;
        BBox min_search(Coordinates(-MAX_FLOAT, -MAX_FLOAT), Coordinates(to_double(p.x()) - EPS_DOUBLE, to_double(p.y()) - EPS_DOUBLE), UINT32_MAX);
        BBox max_search(Coordinates(to_double(p.x()) + EPS_DOUBLE, to_double(p.y()) + EPS_DOUBLE), Coordinates(MAX_FLOAT, MAX_FLOAT), UINT32_MAX);

        for(uint32_t e_id : cells.findItems(min_search, max_search)){
            if(g->isInCell(p, e_id, debug_countOpsCGAL)){
                if(settings.VERBOSE_LEVEL > 2) printf("Cell Found: needed %d steps\n", debug_countOpsCGAL);
                return e_id;
            }
        }
        if(settings.VERBOSE_LEVEL > 2) printf("No Cell Found: needed %d steps\n", debug_countOpsCGAL);
        return UINT32_MAX; // no edge found
    };

    void permute(std::vector<uint32_t> & permutation){
        for(uint32_t i = 0; i < cells.size(); i++){
            cells[i].item = permutation[cells[i].item];
        }
    }
};
#endif