#ifndef tests_graph_h
#define tests_graph_h

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/timer.hpp"
#include "../utils/section.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"
#include "../utils/unionFind.hpp"

class Graph{
public:
    std::vector<Point> coords;
    std::vector<Edge_t> edges;
    uint32_t constrainedEdgeDivider = UINT32_MAX, constrainedNodeDivider = 0;
// public:
    inline const uint32_t nodeCount() const { return coords.size(); }
    inline const uint32_t edgeCount() const { return edges.size(); }
    inline const bool isNodeConstrained(const uint32_t n_id) const { return n_id < constrainedNodeDivider; }
    inline const bool isEdgeConstrained(const uint32_t e_id) const { return edges[e_id >> 1].isConstrained(); }
    inline const bool leftTurn(const uint32_t n0, const uint32_t n1, const uint32_t n2) const { 
        return CGAL::left_turn(coords[n0], coords[n1], coords[n2]);
    }
    inline const bool rightTurn(const uint32_t n0, const uint32_t n1, const uint32_t n2) const { 
        return CGAL::right_turn(coords[n0], coords[n1], coords[n2]); 
    }
    //inline const bool hasNext(const uint32_t e_id){ return nextEdge(e_id) != UINT32_MAX; }
    inline const uint32_t nextEdge(const uint32_t e_id) const {
        assert(e_id != UINT32_MAX);
        return edges[e_id >> 1].E[e_id & 1];
    }
    inline const uint32_t otherEdge(const uint32_t e_id) const { return e_id ^ 1; }
    inline const bool hasNode(const uint32_t e_id, const uint32_t n_id) const { return firstNode(e_id) == n_id || secondNode(e_id) == n_id; }
    inline const uint32_t firstNode(const uint32_t e_id) const { return edges[e_id >> 1].N[1 ^ e_id & 1]; }
    inline const uint32_t secondNode(const uint32_t e_id) const { return edges[e_id >> 1].N[e_id & 1]; }

    // go over all edges and find cell containin Point p, runs in O(n) time
    inline const uint32_t getNearestEdge_naiv(const Point & p);
    inline const bool isInCell(const Point & p, const uint32_t e_id, uint32_t & debug_countOpsCGAL);

    //test functions
    bool allTriangles();
    bool test_allCells_finite();
    uint32_t get_predecessorEdge(uint32_t e);
    void convexifier(const std::vector<bool> & constrained_edges);
    std::vector<uint32_t> deleteUnusedEdges();
    uint32_t test_all_cellsConvex(bool ccw = true);
    uint32_t test_no_cellMergePossible(const std::vector<bool> & constrained_edges);
    uint32_t calc_countTriangles(uint32_t constrainedEdgeDivider);

    bool test_connectedGraph(){
        std::vector<bool> visited_edges(edges.size(), false);
        std::list<uint32_t> queue;
        queue.push_back(0);
        visited_edges[0] = true;
        while(!queue.empty()){
            uint32_t current_edge = queue.front() >> 1;
            queue.pop_front();
            if(settings.VERBOSE_LEVEL > 3) printf("visiting edge %u (%u %u)\n", current_edge * 2, firstNode(current_edge * 2), secondNode(current_edge * 2));
            if(edges[current_edge].E[0] != UINT32_MAX && !visited_edges[edges[current_edge].E[0] >> 1]){
                queue.push_back(edges[current_edge].E[0]);
                visited_edges[edges[current_edge].E[0] >> 1] = true;
            }
            if(edges[current_edge].E[1] != UINT32_MAX && !visited_edges[edges[current_edge].E[1] >> 1]){
                queue.push_back(edges[current_edge].E[1]);
                visited_edges[edges[current_edge].E[1] >> 1] = true;
            }
        }
        for(uint32_t n = 0; n < visited_edges.size(); n++){
            if(!visited_edges[n]){
                std::cout << "edge " << n << " not reachable\n";
                return false;
            }
        }
        return true;
    }

    void do_indexShift(){
        bool needIndexShift = false;
        for(uint32_t e = 0; e < edges.size() * 2; e++){
            auto next_e = nextEdge(e);
            if(next_e == UINT32_MAX) continue;
            if(!edges[next_e >> 1].hasNode(firstNode(e)) && !edges[next_e >> 1].hasNode(secondNode(e))){
                needIndexShift = true;
                break;
            }
        }
        if(needIndexShift){
            auto section = Section("let edges point to one side of other edge, via id = id * 2 + direction");
            for(uint32_t e = 0; e < edges.size() * 2; e++){
                if(nextEdge(e) != UINT32_MAX && !edges[nextEdge(e)].hasNode(secondNode(e))) std::swap(edges[e >> 1].E[0], edges[e >> 1].E[1]);
                assert(nextEdge(e) == UINT32_MAX || edges[nextEdge(e)].hasNode(secondNode(e)));
            }
            for(uint32_t e = 0; e < edges.size() * 2; e++){
                if(nextEdge(e) != UINT32_MAX) edges[e >> 1].E[e & 1] = nextEdge(e) * 2 + (1 ^ edges[nextEdge(e)].getIndex_n(secondNode(e)));
            }
            section.end_section();
        }
    };

    void set_constrainedEdgeDivider(){
        bool needPermutation = false;
        constrainedEdgeDivider = 0;
        for(uint32_t e = 0; e < edges.size(); e++){
            if(edges[e].isConstrained()) {
                if(constrainedEdgeDivider != e) needPermutation  = true; 
                constrainedEdgeDivider++;
            }
        }
        if(needPermutation){
            for(uint32_t e = 0; e < edges.size() * 2; e++){
                assert(nextEdge(e) == UINT32_MAX || edges[nextEdge(e) >> 1].hasNode(secondNode(e)));
            }
            auto section = Section("reorder Nodes so that constrained/obstacle Nodes get lower id", 2);
            std::vector<uint32_t> permutation(edges.size());
            for(uint32_t i = 0; i < permutation.size(); i++) permutation[i] = i;
            for(uint32_t l = 0, r = constrainedEdgeDivider; l + r < edges.size() + constrainedEdgeDivider;){
                if(edges[l].isConstrained()){
                    assert(l < constrainedEdgeDivider);
                    l++;
                }else if(!edges[r].isConstrained()){
                    assert(r < permutation.size());
                    r++;
                }else{
                    std::swap(edges[l], edges[r]);
                    std::swap(permutation[l], permutation[r]);
                    assert(!edges[r].isConstrained());
                    l++; r++;
                }
            }
            for(uint32_t e = 0; e < edges.size() * 2; e++){
                if(nextEdge(e) == UINT32_MAX) continue;
                edges[e >> 1].E[e & 1] = (permutation[nextEdge(e) >> 1] << 1) + (nextEdge(e) & 1);
                assert(edges[nextEdge(e) >> 1].hasNode(secondNode(e)) || edges[nextEdge(e) >> 1].hasNode(firstNode(e)));
                assert(edges[nextEdge(e) >> 1].hasNode(secondNode(e)));
            }
            section.end_section();
        }
    }

    void set_constrainedNodeDivider(){
        bool needPermutation = false;
        constrainedNodeDivider = 0;
        std::vector<bool> node_is_constrained(coords.size(), false);
        for(uint32_t e = 0; e < edges.size() * 2; e++){
            if(edges[e >> 1].isConstrained()){
                node_is_constrained[edges[e >> 1].N[e & 1]] = true;
            }
        }
        for(uint32_t n = 0; n < coords.size(); n++){
            if(node_is_constrained[n]) {
                if(constrainedNodeDivider != n) needPermutation = true; 
                constrainedNodeDivider++;
            }
        }
        if(needPermutation){
            auto section = Section("reorder Nodes so that constrained/obstacle Nodes get lower id", 2);
            std::vector<uint32_t> permutation(coords.size());
            for(uint32_t i = 0; i < permutation.size(); i++) permutation[i] = i;
            for(uint32_t l = 0, r = constrainedNodeDivider; l + r < permutation.size() + constrainedNodeDivider;){
                if(node_is_constrained[l]){
                    l++;
                }else if(!node_is_constrained[r]) {
                    r++;
                }else{
                    std::swap(coords[l], coords[r]);
                    std::swap(permutation[l], permutation[r]);
                    l++; r++;
                }
            }
            for(uint32_t e = 0; e < edges.size() * 2; e++){
                if(secondNode(e) != UINT32_MAX) edges[e >> 1].N[e & 1] = permutation[secondNode(e)];
            }
            section.end_section();
        }
    }

    Graph(const std::string & path);
    void write_polyanyaStyle(const std::string & outputFile);
    void write_CGAL_off_style(const std::string & outputFile);

    const double inline getDist_inexact_index(const uint32_t p0, const uint32_t p1) { return getDist_inexact(coords[p0], coords[p1]); };
    const NT_Distance inline getDist_exact_index(const uint32_t p0, const uint32_t p1) { return getDist_exact(coords[p0], coords[p1]); };
    const bool inline collinear_index(const uint32_t p0, const uint32_t p1, const uint32_t p2) { 
        if(p0 == p1 || p0 == p2 || p1 == p2) return true;
        return CGAL::collinear(coords[p0], coords[p1], coords[p2]);
    };

    void flipEdgesForMoreDeadends();

    void saveGraphData(std::string fileName){
        try {
            if(boost::algorithm::ends_with(fileName, ".bin")){
                FILE* fp = fopen(fileName.c_str(), "w+b");
                save_to_bin(fp, coords, edges);
                fclose(fp);
            }else if(boost::algorithm::ends_with(fileName, ".mesh")){
                write_polyanyaStyle(fileName);
            }else if(boost::algorithm::ends_with(fileName, ".off")){
                write_CGAL_off_style(fileName);
            }else{
                std::ofstream out_stream (fileName);
                save_to_text(out_stream, coords, edges);
                out_stream.close();
            };
        }catch (...) {
            std::cout << "coudn't read " << fileName << std::endl;
        }
    };

    void loadGraphData(std::string fileName){
        try {
            if(boost::algorithm::ends_with(fileName, ".bin")){
                FILE* fp = fopen(fileName.c_str(), "r+b");
                load_from_bin(fp, coords, edges);
                fclose(fp);
            }else if(boost::algorithm::ends_with(fileName, ".mesh")){
                std::ifstream in_stream (fileName);
                load_from_mesh(in_stream, coords, edges);
                in_stream.close();
            }else{
                std::ifstream in_stream (fileName);
                load_from_text(in_stream, coords, edges);
                in_stream.close();
            };
        }catch (...) {
            std::cout << "coudn't read " << fileName << std::endl;
            throw std::runtime_error("Failed to load data from " + fileName);
        }
    };
};

struct FaceCirculator{
    uint32_t first_e, current_e;
    Graph* g;
    bool hasNext_b = true;

    FaceCirculator(Graph* g, const uint32_t e_id) : g(g), first_e(e_id), current_e(e_id) {};
    inline const uint32_t firstNode(){ return g->firstNode(current_e); }
    inline const uint32_t secondNode(){ return g->secondNode(current_e); }
    void next(){
        current_e = g->nextEdge(current_e);
        hasNext_b = current_e != first_e && current_e != UINT32_MAX;
    };
    bool hasNext(){ return hasNext_b; }
    FaceCirculator& operator++(){
        next();
        return *this;
    }
    FaceCirculator operator++(int){
		FaceCirculator tmp(*this);
		next();
		return tmp;
	}
};

struct NodeCirculator{
    uint32_t first_e, current_e;
    Graph* g;
    bool hasNext_b = true;

    NodeCirculator(Graph* g, const uint32_t e_id) : g(g), first_e(e_id), current_e(e_id) {};
    inline const uint32_t firstNode(){ return g->firstNode(current_e); }
    inline const uint32_t secondNode(){ return g->secondNode(current_e); }
    void next(){
        current_e = g->otherEdge(g->nextEdge(current_e));
        hasNext_b = current_e != first_e && (current_e | 1) != UINT32_MAX;
    };
    bool hasNext(){ return hasNext_b; }
    NodeCirculator& operator++(){
        next();
        return *this;
    }
    NodeCirculator operator++(int){
		NodeCirculator tmp(*this);
		next();
		return tmp;
	}
};

Graph::Graph(const std::string & path){
    auto section = Section("load Nodes and Edges", 1);
    loadGraphData(path);
    section.end_section();
    
    if(settings.use_merkartor){
        for(const auto & c : coords){
            if(c.y() < -180.5 || c.y() > 180.5 || c.x() < -90.5 || c.x() > 90.5){
                settings.use_merkartor = false;
                if(settings.VERBOSE_LEVEL > 0){
                    std::cout << "warning: disabling merkartor coordinates, as coordinates are out of lat/lon bounds\n";
                    printf("         coordinate (%.6f, %.6f) is out of bounds\n", to_double(c.x()), to_double(c.y()));
                }
                break;
            }
        }
    }

    settings.setBounds(coords);
    if(settings.use_merkartor){
        for(uint32_t i = 0; i < coords.size(); i++){
            coords[i] = bounded_mercator(coords[i]);
        }
    }

    do_indexShift();
    //assert(test_connectedGraph());
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        //if(nextEdge(e) != UINT32_MAX && !edges[nextEdge(e)].hasNode(secondNode(e))) std::swap(edges[e >> 1].N[0], edges[e >> 1].N[1]);
        assert(nextEdge(e) == UINT32_MAX || edges[nextEdge(e) >> 1].hasNode(secondNode(e)));
    }
    set_constrainedEdgeDivider();
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        //if(nextEdge(e) != UINT32_MAX && !edges[nextEdge(e)].hasNode(secondNode(e))) std::swap(edges[e >> 1].N[0], edges[e >> 1].N[1]);
        assert(nextEdge(e) == UINT32_MAX || edges[nextEdge(e) >> 1].hasNode(secondNode(e)));
    }
    set_constrainedNodeDivider();

    if(settings.polyanya.flipEdgesForMoreDeadends) flipEdgesForMoreDeadends();

    assert(coords.size() > 0);
    if(settings.DEBUG_LEVEL > 1) {
        for(uint32_t e = 0; e < edges.size() * 2; e++){
            //if(nextEdge(e) != UINT32_MAX && !edges[nextEdge(e)].hasNode(secondNode(e))) std::swap(edges[e >> 1].N[0], edges[e >> 1].N[1]);
            assert(nextEdge(e) == UINT32_MAX || edges[nextEdge(e) >> 1].hasNode(secondNode(e)));
        }

        for(uint32_t i = 0; i < edges.size() * 2; i++){
            if(edges[i >> 1].E[i & 1] == UINT32_MAX) continue;
            uint32_t n0 = edges[i >> 1].N[i & 1];
            uint32_t n1 = edges[i >> 1].other_n(n0);
            uint32_t otherEdge = edges[i >> 1].E[i & 1] >> 1;
            assert(edges[otherEdge].hasNode(n0) || edges[otherEdge].hasNode(n1));
            uint32_t n2 = edges[otherEdge].other_n(n0);
            assert(!CGAL::right_turn(coords[n1], coords[n0], coords[n2]));
        }
    }
    
    uint32_t errors;
    if(settings.DEBUG_LEVEL > 1) assert(test_allCells_finite());
    if(settings.DEBUG_LEVEL > 1) assert(allTriangles());
    if(settings.DEBUG_LEVEL > 1) errors = test_all_cellsConvex(true);
    if(settings.VERBOSE_LEVEL >= 0) printf("loaded %li Nodes and %li Edges\n", coords.size(), edges.size());
};

void Graph::flipEdgesForMoreDeadends(){
    if(!allTriangles()) {
        if(settings.VERBOSE_LEVEL > 0) printf("skipping edge flipping, graph not triangulated\n");
        return;
    }
    auto uf = UnionFind(nodeCount());
    for(uint32_t i = 0; i < constrainedEdgeDivider; i++){
        assert(edges[i].isConstrained());
        uf.merge(edges[i].N[0], edges[i].N[1]);
    }
    std::list<uint32_t> queue;
    std::vector<bool> edge_is_deadend(edges.size(), false);
    for(uint32_t e = 0; e < edges.size(); e++){
        edge_is_deadend[e] = uf.is_connected(edges[e].N[0], edges[e].N[1]);
        if(!edge_is_deadend[e]) queue.push_back(e << 1);
    }
    uint32_t flipped_edges = 0;
    while(!queue.empty()){
        uint32_t current_edge = queue.front(); queue.pop_front();

        uint32_t n0 = firstNode(current_edge), n1 = secondNode(current_edge);
        if(edge_is_deadend[current_edge >> 1]) continue;

        assert(firstNode(nextEdge(current_edge)) == n1);
        assert(firstNode(nextEdge(otherEdge(current_edge))) == n0);
        uint32_t next_n = secondNode(nextEdge(current_edge));
        uint32_t next_n_other = secondNode(nextEdge(otherEdge(current_edge)));

        if(!uf.is_connected(next_n, next_n_other)) continue;
        if(!CGAL::left_turn(coords[next_n_other], coords[next_n], coords[n0])) continue;
        if(!CGAL::right_turn(coords[next_n_other], coords[next_n], coords[n1])) continue;

        //do edgeflip
        uint32_t f0_0 = current_edge, f1_0 = otherEdge(current_edge);
        uint32_t f0_1 = nextEdge(f0_0), f1_1 = nextEdge(f1_0);
        uint32_t f0_2 = nextEdge(f0_1), f1_2 = nextEdge(f1_1);

        // f0_1 -> f1_0
        edges[f0_1 >> 1].setNeighbour(f1_0, next_n);
        // f1_1 -> f0_0
        edges[f1_1 >> 1].setNeighbour(f0_0, next_n_other);
        // f0_0 -> f0_2
        edges[f0_0 >> 1].setNeighbour(f0_2, n1);
        // f1_0 -> f1_2
        edges[f1_0 >> 1].setNeighbour(f1_2, n0);
        // f0_2 -> f1_1
        edges[f0_2 >> 1].setNeighbour(f1_1, n0);
        // f1_2 -> f0_1
        edges[f1_2 >> 1].setNeighbour(f0_1, n1);
        // f0_0 (n0,n1) -> (next_n_other, next_n)
        edges[f0_0 >> 1].N[f0_0 & 1] = next_n;
        edges[f1_0 >> 1].N[f1_0 & 1] = next_n_other;

        edge_is_deadend[current_edge >> 1] = true;
        flipped_edges++;

        if(!edge_is_deadend[f0_1 >> 1]) queue.push_back(f0_1);
        if(!edge_is_deadend[f1_1 >> 1]) queue.push_back(f1_1);
        if(!edge_is_deadend[f0_2 >> 1]) queue.push_back(f0_2);
        if(!edge_is_deadend[f1_2 >> 1]) queue.push_back(f1_2);
    }
    if(settings.VERBOSE_LEVEL > 1) printf("flipped %i edges to increase deadends\n", flipped_edges);
}

uint32_t Graph::get_predecessorEdge(uint32_t e){
    for(FaceCirculator iter(this, e); iter.hasNext(); iter.next()){
        if(nextEdge(iter.current_e) == e) return iter.current_e;
    }
    return UINT32_MAX;
}
inline const uint32_t Graph::getNearestEdge_naiv(const Point & p){
    for(uint32_t e_id = 0; e_id < edges.size() * 2; e_id++){
        if(isEdgeConstrained(e_id)) continue; // skip constrained edges
        uint32_t cgalOps = 0;
        if(isInCell(p, e_id, cgalOps)) return e_id;
    }
    return UINT32_MAX; // no edge found
};

inline const bool Graph::isInCell(const Point & p, const uint32_t e_id, uint32_t & debug_countOpsCGAL){
    for(FaceCirculator f_iter(this, e_id); f_iter.hasNext(); f_iter++){
        assert(f_iter.current_e != UINT32_MAX);
        if(settings.VERBOSE_LEVEL > 3) 
            printf("node %s and %s\n",
                to_string(coords[f_iter.firstNode()]).c_str(),
                to_string(coords[f_iter.secondNode()]).c_str());
        
        assert(nextEdge(f_iter.current_e) != UINT32_MAX);
        assert(firstNode(nextEdge(f_iter.current_e)) == f_iter.secondNode());
        assert(firstNode(f_iter.current_e) == f_iter.firstNode());
        assert(secondNode(f_iter.current_e) == f_iter.secondNode());
        
        debug_countOpsCGAL++;
        if(CGAL::right_turn(coords[f_iter.firstNode()], coords[f_iter.secondNode()], p)){
            return false;
        }
    }
    return true;
};

bool Graph::allTriangles(){
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        if(edges[e >> 1].E[e & 1] == UINT32_MAX) continue;

        uint32_t edges_per_cell = 0;
        for(FaceCirculator iter(this, e); iter.hasNext(); iter++){
            edges_per_cell++;
        }
        if(edges_per_cell != 3){
            printf("ERROR: Cell has %i edges, should be 3\n", edges_per_cell);
            return false;
        }
    }
    return true;
}

bool Graph::test_allCells_finite(){
    uint32_t maxCellSize = 0;
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        if(nextEdge(e) == UINT32_MAX) continue;
        assert(firstNode(e) < coords.size());
        assert(secondNode(e) < coords.size());
        uint32_t currentCellSize = 0;
        
        uint32_t prevNode = firstNode(e);
        for(FaceCirculator f_iter(this, e); f_iter.hasNext(); f_iter.next()){
            assert(prevNode == firstNode(f_iter.current_e));
            prevNode = secondNode(f_iter.current_e);
            if(currentCellSize > edges.size()){
                printf("error at edge %d", e);
                return false;
            } 
            currentCellSize++;
        }
        maxCellSize = std::max(maxCellSize, currentCellSize);
    }
    if(settings.VERBOSE_LEVEL > 2) printf("max cell size is %i\n", maxCellSize);
    return true;
}

// tries to reduce number of cells by joining them, ensures the convexity
// don't join cells if edge ist constrained
// removed edges get invalid values but are still in memory
bool compare_convexifier(const std::pair<double, uint32_t> & a, const std::pair<double, uint32_t> & b){ return a.first > b.first; };
void Graph::convexifier(const std::vector<bool> & constrained_edges){
    if(settings.DEBUG_LEVEL > 1) printf("found %i non-convex cells\n", test_all_cellsConvex());
    if(settings.DEBUG_LEVEL > 1) printf("found %i cell that can be merged\n", test_no_cellMergePossible(constrained_edges));

    uint32_t removed_edges = 0;
    typedef std::priority_queue<std::pair<double, uint32_t>, std::vector<std::pair<double, uint32_t>>, bool(*)(const std::pair<double, uint32_t>&, const std::pair<double, uint32_t>&)> PQ;
    PQ queue(compare_convexifier);
    for(uint32_t e = 0; e < edges.size(); e++){
        if(edges[e].isConstrained()) continue;
        queue.emplace(to_double(CGAL::squared_distance(coords[edges[e].N[0]], coords[edges[e].N[1]])), e << 1);
    } 
    while(!queue.empty()){
        uint32_t e = queue.top().second;
        queue.pop();
        assert(e != UINT32_MAX && !edges[e >> 1].isConstrained());
        if(constrained_edges[e >> 1]){
            bool whole_cell_constrained = true;
            for(FaceCirculator iter(this, e); iter.hasNext(); iter.next()){
                if(!constrained_edges[iter.current_e >> 1]){
                    whole_cell_constrained = false;
                    break;
                }
            }
            for(FaceCirculator iter(this, otherEdge(e)); iter.hasNext(); iter.next()){
                if(!constrained_edges[iter.current_e >> 1]){
                    whole_cell_constrained = false;
                    break;
                }
            }
            if(!whole_cell_constrained) continue;
        }

        uint32_t n0 = edges[e >> 1].N[0];
        uint32_t n1 = edges[e >> 1].N[1];
        auto pre1 = get_predecessorEdge(e);
        assert(FaceCirculator(this, pre1).secondNode() == n1);
        auto pre2 = get_predecessorEdge(e ^ 1);
        assert(FaceCirculator(this, pre2).secondNode() == n0);

        if(!CGAL::left_turn(coords[firstNode(pre1)], coords[n1], coords[secondNode(nextEdge(otherEdge(e)))])) continue;
        if(!CGAL::left_turn(coords[firstNode(pre2)], coords[n0], coords[secondNode(nextEdge(e))])) continue;

        edges[pre1 >> 1].E[pre1 & 1] = nextEdge(otherEdge(e));
        edges[pre2 >> 1].E[pre2 & 1] = nextEdge(e);

        edges[e >> 1] = Edge_t(UINT32_MAX, UINT32_MAX, UINT32_MAX, UINT32_MAX);
        removed_edges++;
    }
    if(settings.VERBOSE_LEVEL > 0) printf("removed %i out of %li edges\n\t%f %% removed\n", removed_edges, edges.size(), double(removed_edges) / double(edges.size()) * 100.0);
    if(settings.DEBUG_LEVEL > 1) printf("found %i non-convex cells\n", test_all_cellsConvex());
    if(settings.DEBUG_LEVEL > 1) printf("found %i cell that can be merged\n", test_no_cellMergePossible(constrained_edges));
};

// deletes edges with n0 == n1 == UINT32_MAX or e0 == e1 == UINT32_MAX
std::vector<uint32_t> Graph::deleteUnusedEdges(){
    std::vector<uint32_t> permutation(edges.size(), UINT32_MAX);
    uint32_t newEdgeCount = 0;
    for(uint32_t e = 0; e < edges.size(); e++){
        if((edges[e].N[0] == UINT32_MAX && edges[e].N[1] == UINT32_MAX) || (edges[e].E[0] == UINT32_MAX && edges[e].E[1] == UINT32_MAX)) continue;
        permutation[e] = newEdgeCount++;
    }
    for(uint32_t e = 0; e < edges.size(); e++){
        if((edges[e].N[0] == UINT32_MAX && edges[e].N[1] == UINT32_MAX) || (edges[e].E[0] == UINT32_MAX && edges[e].E[1] == UINT32_MAX)) continue;
        if(edges[e].E[0] != UINT32_MAX){
            assert(permutation[edges[e].E[0] >> 1] != UINT32_MAX);
        }
        if(edges[e].E[1] != UINT32_MAX && permutation[edges[e].E[1] >> 1] == UINT32_MAX){
            auto temp = edges[edges[e].E[1] >> 1];
            auto temp2 = edges[e].E[1] >> 1;
            assert(false);
        } 
    }
    for(uint32_t e = 0; e < edges.size(); e++){
        if((edges[e].N[0] == UINT32_MAX && edges[e].N[1] == UINT32_MAX) || (edges[e].E[0] == UINT32_MAX && edges[e].E[1] == UINT32_MAX)) continue;
        assert(permutation[e] != UINT32_MAX);
        if(edges[e].E[0] != UINT32_MAX){
            assert(permutation[edges[e].E[0] >> 1] != UINT32_MAX);
            edges[e].E[0] = (permutation[edges[e].E[0] >> 1] << 1) ^ (edges[e].E[0] & 1);
        }
        if(edges[e].E[1] != UINT32_MAX){
            assert(permutation[edges[e].E[1] >> 1] != UINT32_MAX);
            edges[e].E[1] = (permutation[edges[e].E[1] >> 1] << 1) ^ (edges[e].E[1] & 1);
        } 
        assert(permutation[e] <= e);
        edges[permutation[e]] = edges[e];
    }
    if(settings.VERBOSE_LEVEL > 0) printf("removed %ld unused edges -- new edgeCount: %d", edges.size() - newEdgeCount, newEdgeCount);
    edges.resize(newEdgeCount);
    set_constrainedEdgeDivider();
    return permutation;
};

uint32_t Graph::test_all_cellsConvex(bool ccw){
    uint32_t errors = 0;
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        if(edges[e >> 1].E[e & 1] == UINT32_MAX) continue;

        for(FaceCirculator iter(this, e); iter.hasNext(); iter.next()){
            uint32_t n0 = iter.firstNode();
            uint32_t n1 = iter.secondNode();
            uint32_t n2 = secondNode(nextEdge(iter.current_e));

            assert(n0 != n1 && n1 != n2 && n0 != n2);
            assert(firstNode(nextEdge(iter.current_e)) == n1);
            if(ccw && CGAL::right_turn(coords[n0], coords[n1], coords[n2]) || !ccw && CGAL::left_turn(coords[n0], coords[n1], coords[n2])){
                errors++;
                if(settings.VERBOSE_LEVEL > 1) 
                    printf("ERROR: Cell is not convex %s, %s, %s\n", 
                        to_string(coords[n0]).c_str(), 
                        to_string(coords[n1]).c_str(),
                        to_string(coords[n2]).c_str());
            }
        }
    }
    return errors;
};

uint32_t Graph::test_no_cellMergePossible(const std::vector<bool> & constrained_edges){
    uint32_t errors = 0;
    std::vector<bool> possible_toMerge(edges.size() * 2, false);
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        if(edges[e >> 1].E[e & 1] == UINT32_MAX) continue;
        assert(nextEdge(e) != UINT32_MAX);

        FaceCirculator iter(this, e);
        uint32_t n0 = iter.firstNode();
        uint32_t n1 = iter.secondNode();
        iter.next();
        assert(iter.current_e != UINT32_MAX);
        if(constrained_edges[nextEdge(iter.current_e) >> 1]) continue;
        if(nextEdge(otherEdge(iter.current_e)) == UINT32_MAX) continue;
        
        uint32_t n2 = secondNode(nextEdge(otherEdge(iter.current_e)));

        if(CGAL::left_turn(coords[n0], coords[n1], coords[n2])){
            assert(!possible_toMerge[iter.current_e]);
            if(possible_toMerge[iter.current_e ^ 1] == true){
                //printf("ERROR: Cell can be merged %f %f, %f %f, %f %f\n", coords[n0].x(), coords[n0].y(), coords[n1].x(), coords[n1].y(), coords[n2].x(), coords[n2].y());
                errors++;
            };
            possible_toMerge[iter.current_e] = true;
        }
    }
    return errors;
};

uint32_t Graph::calc_countTriangles(uint32_t constrainedEdgeDivider){
    assert(allTriangles());
    uint32_t edgeSides = (edges.size() * 2) - constrainedEdgeDivider;
    uint32_t debugCounter = 0;
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        assert(edges[e >> 1].E[e & 1] != edges[e >> 1].E[(e & 1) ^ 1]);
        if(edges[e >> 1].E[e & 1] == UINT32_MAX) continue;
        if(e < constrainedEdgeDivider * 2) assert(edges[e >> 1].E[(e & 1) ^ 1] == UINT32_MAX);
        debugCounter++;
    }
    assert(!edges[constrainedEdgeDivider].isConstrained());
    assert(edges[constrainedEdgeDivider - 1].isConstrained());
    assert(debugCounter == edgeSides);
    assert(edgeSides % 3 == 0);
    return edgeSides / 3;
};


// mesh
// 2
// numberNodes numberPolygons
// x y neighbourCount p0 p1 ... // points
// cornerCount n0 n1 ... p0 p1 ...                   //polygons n0 n1 -> is edge to p0
// 2 2 2 -1 4
// 1 3 2 15 -1
// 2 3 6 1 15 -1 4 0 14
// 3 2 5 2 0 4 -1 5
// 3 1 2 5 -1
// 15 1 3 -1 2 5
// 15 3 6 6 3 14 0 2 -1
// ...
// 3 60 23 91 96 115 93
// 3 61 62 63 119 116 -1
// 3 61 63 76 82 118 84
void Graph::write_polyanyaStyle(const std::string & outputFile){
    // preprocessing
    std::vector<uint32_t> node_to_first_edge(nodeCount(), UINT32_MAX);
    std::vector<uint32_t> nodeDegree(nodeCount(), 0);
    for(uint32_t e = 0; e < edgeCount() * 2; e++){
        uint32_t n = edges[e >> 1].N[(e & 1)], otherNextEdge = edges[e >> 1].E[(e & 1) ^ 1];
        nodeDegree[n]++;
        if(nextEdge(e) == UINT32_MAX) continue;
        if(node_to_first_edge[n] == UINT32_MAX || otherNextEdge == UINT32_MAX){
            node_to_first_edge[n] = e;
        }
    }
    std::vector<uint32_t> edge_to_face(edgeCount() * 2, UINT32_MAX);
    std::vector<uint32_t> face_to_edge, faceSize;
    for(uint32_t e = 0; e < edgeCount() * 2; e++){
        if(nextEdge(e) == UINT32_MAX) continue;
        if(edge_to_face[e] != UINT32_MAX) continue;
        faceSize.emplace_back(0);
        for(FaceCirculator f_iter(this, e); f_iter.hasNext(); f_iter++){
            edge_to_face[f_iter.current_e] = face_to_edge.size();
            faceSize.back()++;
        }
        face_to_edge.emplace_back(e);
    }
    // output
    try{
        std::ofstream outStream (outputFile, std::ios::trunc);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << "mesh\n2\n" << coords.size() << " " << face_to_edge.size() << '\n';

        for(uint32_t n = 0; n < coords.size(); n++){
            nodeDegree[n] = 0;
            for(NodeCirculator n_iter(this, node_to_first_edge[n]); n_iter.hasNext(); n_iter++){nodeDegree[n]++;}
            assert(nodeDegree[n] != 1);
            outStream << coords[n].x() << " " << coords[n].y() << " " << nodeDegree[n];
            for(NodeCirculator n_iter(this, node_to_first_edge[n]); n_iter.hasNext(); n_iter++){
                if(edge_to_face[n_iter.current_e] == UINT32_MAX){
                    outStream << " -1";
                }else{
                    outStream << " " << edge_to_face[n_iter.current_e];
                }
            }
            outStream << '\n';
        }
        for(uint32_t f = 0; f < face_to_edge.size(); f++){
            assert(faceSize[f] > 2);
            outStream << faceSize[f];
            for(FaceCirculator f_iter(this, face_to_edge[f]); f_iter.hasNext(); f_iter++){
                outStream << " " << f_iter.secondNode(); // ab: a -> abc
            }
            for(FaceCirculator f_iter(this, face_to_edge[f]); f_iter.hasNext(); f_iter++){
                if(edge_to_face[otherEdge(f_iter.current_e)] == UINT32_MAX){
                    outStream << " -1";
                }else{
                    outStream << " " << edge_to_face[otherEdge(f_iter.current_e)];
                }
            }
            outStream << '\n';
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}


//OFF
//nodeCount faceCount 0?
//x y z
//...
//faceSize n0 n1 n2 n3 ...
void Graph::write_CGAL_off_style(const std::string & outputFile){
    // preprocessing
    uint32_t faceCount = 0;
    std::vector<bool> edgeCovered(edges.size() * 2, false);
    for(uint32_t e = 0; e < edgeCount() * 2; e++){
        if(edgeCovered[e]) continue;
        if(nextEdge(e) == UINT32_MAX) continue;
        faceCount++;
        for(FaceCirculator f_iter(this, e); f_iter.hasNext(); f_iter++){
            assert(edgeCovered[f_iter.current_e] == false);
            edgeCovered[f_iter.current_e] = true;
        }
    }
    // output
    try{
        std::ofstream outStream (outputFile, std::ios::trunc);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << "OFF\n" << coords.size() + 1 << " " << faceCount << " 0\n";

        for(uint32_t n = 0; n < coords.size(); n++){
            outStream << coords[n].x() << " " << coords[n].y() << " 0\n";
        }
        outStream << "0 0 100000000\n";

        edgeCovered = std::vector<bool>(edges.size() * 2, false);
        for(uint32_t e = 0; e < edgeCount() * 2; e++){
            if(edgeCovered[e]) continue;
            if(nextEdge(e) == UINT32_MAX){
                continue;
                // outStream << "3 " << firstNode(e) << " " << secondNode(e) << " " << coords.size() << '\n';
            }else{
                uint32_t faceSize = 0;
                for(FaceCirculator f_iter(this, e); f_iter.hasNext(); f_iter++){
                    edgeCovered[f_iter.current_e] = true;
                    faceSize++;
                }
                outStream << faceSize;
                for(FaceCirculator f_iter(this, e); f_iter.hasNext(); f_iter++){
                    outStream << " " << f_iter.secondNode(); // ab: a -> abc
                }
                outStream << '\n';
            }
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

#endif