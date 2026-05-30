#ifndef polyanya_h
#define polyanya_h

#include <queue>
#include <list>
#include <optional>
#include <variant>
#include <sstream>

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

#include "graph.hpp"
#include "deadendFinder.hpp"
#include "edgeLocator.hpp"
#include "landmarks.hpp"
#include "edgePruning.hpp"
#include "heuristikEdgePruning.hpp"

#define LEFT true
#define RIGHT false

struct PolyanyaStats{
    public:
    uint32_t nodes_generated = 0, wavelets_pushed = 0, pq_pops = 0, wavelets_processed = 0, serach_id = UINT32_MAX, prunedOnNotGenAgain = 0;
    double time_ms = 0.0;
    double distance = -1.0;
    Point start, goal;
    uint32_t hops = 0, numberEdges = 0, max = UINT32_MAX;
    double mean = 0.0;
    uint32_t maxPQSize = 0;
    void reset(){
        nodes_generated = 0; wavelets_pushed = 0; pq_pops = 0; wavelets_processed = 0; prunedOnNotGenAgain = 0;
        hops = 0, numberEdges = 0, mean = 0.0, max = UINT32_MAX;
        time_ms = 0.0; 
        distance = -1.0;
        maxPQSize = 0;
    }
    const void setStats(uint32_t hops, uint32_t numberEdges, double mean, uint32_t max){
        this->hops = hops;
        this->numberEdges = numberEdges;
        this->mean = mean;
        this->max = max;
    }
    const std::string getStats() const {
        return (std::stringstream() << std::setprecision(6)
                                    << serach_id
                            << "\t" << time_ms
                            << "\t" << std::setprecision(std::numeric_limits<double>::max_digits10) << distance << std::setprecision(6)
                            << "\t" << hops
                            << "\t" << numberEdges
                            << "\t" << mean
                            << "\t" << max
                            << "\t" << "(" << start.x() << "," << start.y() << ")"
                            << "\t" << "(" << goal.x() << "," << goal.y() << ")"
                            << "\t" << wavelets_pushed
                            << "\t" << maxPQSize
                            << "\t" << pq_pops
                            << "\t" << wavelets_processed
                            << "\t" << nodes_generated << '\n').str();
    };
    static const std::string getStatsHeadder(){ 
        return (std::stringstream() << "serach_id  "
                           << "\t" << "time_ms    "
                           << "\t" << "distance   "
                           << "\t" << "hops       "
                           << "\t" << "numberEdges"
                           << "\t" << "mean       "
                           << "\t" << "max        "
                           << "\t" << "start      "
                           << "\t" << "goal       "
                           << "\t" << "pq_pushs   "
                           << "\t" << "maxPQSize  "
                           << "\t" << "pq_pops    "
                           << "\t" << "processed  "
                           << "\t" << "generated  " << '\n').str();
    };
};

struct TotalStats{
    public:
    uint32_t count = 0;
    uint64_t nodes_generated = 0, wavelets_pushed = 0, pq_pops = 0, wavelets_processed = 0, numberEdges = 0, hops = 0;
    uint32_t max = 0;
    double time_ms = 0.0;
    uint32_t maxPQSize = 0;

    void addStats(PolyanyaStats& stats){
        count++;
        maxPQSize = std::max(maxPQSize, stats.maxPQSize);
        nodes_generated += stats.nodes_generated;
        wavelets_pushed += stats.wavelets_pushed;
        pq_pops += stats.pq_pops;
        wavelets_processed += stats.wavelets_processed;
        numberEdges += stats.numberEdges;
        hops += stats.hops;
        max = std::max(max, stats.max);
        time_ms += stats.time_ms;
    }

    void reset(){
        count = 0;
        maxPQSize = 0;
        nodes_generated = 0;
        wavelets_pushed = 0;
        pq_pops = 0;
        wavelets_processed = 0;
        numberEdges = 0;
        hops = 0;
        max = 0;
        time_ms = 0.0;
    }

    const std::string getStats(){
        return (std::stringstream() << count
                           << "\t" << time_ms / count
                            << "\t" << get_mem_usage()
                           << "\t" << (double) nodes_generated / count
                           << "\t" << (double) wavelets_pushed / count
                           << "\t" << maxPQSize
                           << "\t" << (double) pq_pops / count
                           << "\t" << (double) wavelets_processed / count
                           << "\t" << (double) numberEdges / count
                           << "\t" << (double) hops / count
                           << "\t" << max
                           << "\t" << (double) wavelets_processed / numberEdges).str();
    };
    const std::string getStatsHeadder(){ 
        return (std::stringstream()<< "count       "
                           << "\t" << "time_ms     "
                           << "\t" << "mem_usage(KB)"
                           << "\t" << "generated   "
                           << "\t" << "pq_pushs    "
                           << "\t" << "maxPQSize   "
                           << "\t" << "pq_pops     "
                           << "\t" << "processed   "
                           << "\t" << "numberEdges "
                           << "\t" << "hops        "
                           << "\t" << "max         "
                           << "\t" << "mean        " << '\n').str();
    };
};

class PolyAnya : public Graph{
protected:
    DeadendFinder deadendFinder;
    Locator locator;
    Landmarks landmarks;
    HeuristikEdgePruning heuristikEdgePruning;
    EdgePruning edgePruning;
public:
    // stats
    PolyanyaStats stats;
    TotalStats accumulated_stats;
    uint32_t waveletIds = 0;
    uint32_t startID = UINT32_MAX, goalID = UINT32_MAX;
    bool startPointInserted = false, goalPointInserted = false;

    //for faster searching
    static bool compare(const std::pair<double, Wavelet> & a, const std::pair<double, Wavelet> & b) { return a.first > b.first; }
    typedef std::priority_queue<std::pair<double, Wavelet>, std::vector<std::pair<double, Wavelet>>, bool(*)(const std::pair<double, Wavelet>&, const std::pair<double, Wavelet>&)> PQ;
    PQ pq = PQ(compare);
    std::vector<uint32_t> previous_node;
	std::vector<NT_Distance> costs;
    // for reset
    std::list<uint32_t> visited;
    std::vector<uint32_t> count_visited;

    // for prunin aka dont generate node with same costs twice
    std::vector<bool> has_generated;
    

    virtual uint32_t insertPoint(const Point & p, uint32_t first_edge_id, bool is_goal);
    virtual void deleteLastPoint();

    virtual void resetSearch();
    const double getDistance_no_heuristik(const Wavelet & w) const;
    const double getDistance(const Wavelet & w) const;
    inline const double getDistance_simpler(const Wavelet & w) const;
    inline const double getDistance_noRefelection(const Wavelet & w) const;

    void testDeadend(){
        std::vector<uint32_t> edges_on_opt_path;
        uint32_t last_n = UINT32_MAX;
        for(auto n : getNodeIdOnRoute()){
            if(last_n != UINT32_MAX) {
                auto seg1 = Segment(coords[n], coords[last_n]);
                for (uint32_t e : visited){
                    auto seg2 = Segment(coords[firstNode(e << 1)], coords[secondNode(e << 1)]);
                    if(CGAL::do_intersect(seg1, seg2) && firstNode(e << 1) != n && firstNode(e << 1) != last_n && secondNode(e << 1) != n && secondNode(e << 1) != last_n){
                        edges_on_opt_path.push_back(e);
                    }
                }
            }
            last_n = n;
        }

        assert(!settings.polyanya.use_deadend_pruning);
        settings.polyanya.use_deadend_pruning = true;
        prepareSearch(stats.start, stats.goal);

        for(auto e : edges_on_opt_path){
            assert(!deadendFinder.isDeadend(e << 1));
        }
    }

    void continuous_dijkstra();

    void initLandmarks(){
        assert(settings.polyanya.use_landmarks);
        landmarks = Landmarks(this, settings.polyanya.landmark_count);
        uint32_t landmark = landmarks.getInitalLandmark();
        for(uint16_t l = 0; l < settings.polyanya.landmark_count; l++){
            auto dists = routeToAll(landmark);
            landmark = landmarks.addLandmark(landmark, dists);
        }
        resetSearch();
    }
    
    const bool isDeadend(const uint32_t edge_id);
    inline const bool isTurn(bool dir, uint32_t n0, uint32_t n1, Point & test) const;
    const inline bool needTurn(uint32_t n, NT_Distance & dist, uint32_t bound, uint32_t generator){
        if(n == generator) return false;
        if(bound != n || costs[n] * onePlusEps < dist) return false;
        bool updated = update_dist_to_node(n, generator, UINT32_MAX);
        //old TODO:
        // if(!has_generated[n] && costs[n] * onePlusEps > dist && previous_node[n] != generator) previous_node[n] = generator;
        // if(previous_node[n] != generator) return false;
        // if(!isNodeConstrained(n)) return false;
        //new TODO:
        if(!isNodeConstrained(n)) return false;
        if(!has_generated[n] && previous_node[n] != generator) assert(costs[n] < dist);
        if(!has_generated[n] && previous_node[n] != generator && updated) previous_node[n] = generator;
        if(previous_node[n] != generator) return false;

        if(has_generated[n]) stats.prunedOnNotGenAgain++;
        return settings.polyanya.dont_gen_twice ? !has_generated[n] : true;
    }
    bool update_dist_to_node(uint32_t n, uint32_t generator, uint32_t edge_id);
    void visitEdge(uint32_t e_id){
        if(e_id == UINT32_MAX) return;
        assert(e_id < count_visited.size());
        if(count_visited[e_id] != UINT32_MAX) return;
        visited.push_back(e_id);
        count_visited[e_id] = 0;
    }

public:
    PolyAnya(const std::string & path);

    const std::list<Coordinates> getRoute();
    const std::list<Coordinates> getRoute_withoutMercator();
    const NT_Distance evalutatePath(std::list<Coordinates>);
    const Exact_Distance evalutatePath_exact(std::list<Coordinates>);
    inline const std::list<uint32_t> getNodeIdOnRoute();
    const std::list<std::pair<Coordinates, Coordinates>> getSearchSpace();
    double getRouteCosts() {return stats.distance;};
    inline void processEdge_new(Wavelet & w, bool addRight, bool addLeft, std::vector<Wavelet>& waveletsToCreate);
    inline void processEdge(Wavelet & w, bool addRight, bool addLeft, std::vector<Wavelet>& waveletsToCreate);
    void routePolyanya(Point start_point, Point goal_point);
    void routePolyanya(uint32_t startID, uint32_t goalID);
    std::vector<NT_Distance> routeToAll(uint32_t start_point);
    inline void addStartPoint(uint32_t start, uint32_t start_edge);
    inline const uint32_t getNearestEdge(const Point & p){
        if(!locator.isLocatorSet()){
            if(settings.VERBOSE_LEVEL > 0) printf("warining: running an O(n) edgeLocation, because locator isn't set!");
            return getNearestEdge_naiv(p);
        } else {
            return locator.getNearestEdge(p);
        }
    };
    inline const uint32_t getAdjacentEdge(uint32_t n){
        assert(n < coords.size());
        for(FaceCirculator f_iter(this, getNearestEdge(coords[n])); f_iter.hasNext(); ++f_iter){
            if(f_iter.secondNode() == n) return  f_iter.current_e;
        }
        throw std::runtime_error("node isn't part of triangulation");
    };
    bool inDomain(const Point & p) { 
        return ((settings.use_merkartor ? getNearestEdge(bounded_mercator(p)) : getNearestEdge(p)) | 1) != UINT32_MAX; 
    }

    std::vector<uint32_t> getCount_of_visited(){ 
        std::vector<uint32_t> count;
        count.reserve(visited.size());
        for(auto & e : visited){
            if(count_visited[e] == UINT32_MAX) continue;
            if(count_visited[e] == 0) continue;
            assert(e < count_visited.size());
            count.push_back(count_visited[e]);
        }
        return count;
    };

    void printHeadder(std::ostream& out = std::cout){
        if(settings.VERBOSE_LEVEL < 0) return;
        out << stats.getStatsHeadder();
    };

    void printStats(std::ostream& out = std::cout){
        auto count_of_visited = getCount_of_visited();
        uint32_t hops = getNodeIdOnRoute().size();
        uint32_t countEdgesVisited = count_of_visited.size();
        double mean = count_of_visited.size() == 0 ? 0.0 : calc_mean(count_of_visited);
        uint32_t max = count_of_visited.size() == 0 ? 0 : *std::max_element(count_of_visited.begin(), count_of_visited.end());
        stats.setStats(hops, countEdgesVisited, mean, max);
        accumulated_stats.addStats(stats);
        if(settings.VERBOSE_LEVEL < 1) return;
        out << stats.getStats();
    };

    void printTotalHeadder(std::ostream& out = std::cout){
        if(settings.VERBOSE_LEVEL < 0) return;
        out << accumulated_stats.getStatsHeadder() << std::endl;
    };
    void printTotalStats(std::ostream& out = std::cout){
        if(settings.VERBOSE_LEVEL < -1) return;
        resetSearch();
        printTotalHeadder(out);
        out << accumulated_stats.getStats() << std::endl;
    };

    std::pair<uint32_t, uint32_t> prepareSearch(const Point & start_point, const Point & goal_point);
};

PolyAnya::PolyAnya(const std::string & path) : Graph(path) {
    Section section;
    //finding point in sea
    Point seaPoint(settings.seaPoint.lat, settings.seaPoint.lon);
    uint32_t s_edge = getNearestEdge_naiv(settings.use_merkartor ? bounded_mercator(seaPoint) : seaPoint);
    if(s_edge == UINT32_MAX){
        if(settings.VERBOSE_LEVEL > 0) printf("Warning: sea point %s is not in triangulation, choosing random point\n", to_string(seaPoint).c_str());
        s_edge = edges.size() * 2 - 1;
        settings.seaPoint = Coordinates(coords[firstNode(s_edge)]);
        seaPoint = to_point(settings.seaPoint);
        seaPoint = settings.use_merkartor ? rev_bounded_mercator(seaPoint) : seaPoint;

        assert(!isEdgeConstrained(s_edge));
    }
    
    if(settings.DEBUG_LEVEL > 2) assert(test_allCells_finite());


    if(settings.polyanya.use_deadend_pruning && !settings.polyanya.merge_cells_before_deadend_pruning){
        section = Section("init deadend finder assuming sea at " + to_string(seaPoint), 1); // (30, -35)
        deadendFinder = DeadendFinder(this, s_edge);
        section.end_section();
    }

    if(settings.polyanya.merge_cells){
        section = Section("merging cells", 1);
        std::vector<bool> constrained_edges(edges.size());
        for(uint32_t e = 0; e < edges.size(); e++){constrained_edges[e] = isDeadend(e << 1);}
        convexifier(constrained_edges);
        section.end_section();
        auto permutation = deleteUnusedEdges();
    }

    if(settings.polyanya.use_deadend_pruning && settings.polyanya.merge_cells){
        s_edge = getNearestEdge_naiv(settings.use_merkartor ? bounded_mercator(seaPoint) : seaPoint);
        section = Section("init deadend finder assuming sea at " + to_string(seaPoint), 1); // (30, -35)
        deadendFinder = DeadendFinder(this, s_edge);
        section.end_section();
    }

    section = Section("init edge locator", 1);
    locator = Locator(this);
    section.end_section();

    // init search structs
    // make space for soure and target node
    previous_node.reserve(coords.size() + 8);
    costs.reserve(coords.size() + 8);
    has_generated.reserve(coords.size() + 8);

    previous_node.resize(coords.size(), UINT32_MAX);
    costs.resize(coords.size(), MAX_DOUBLE);
    has_generated.resize(coords.size(), false);

    count_visited.resize(edges.size(), UINT32_MAX);

    if(settings.polyanya.use_edgePruning){
        Section section("init EdgePruning");
        edgePruning = EdgePruning(this, edges.size(), constrainedEdgeDivider);
        section.end_section();
    }

    if(settings.polyanya.heuristikEdgePruningProof){
        heuristikEdgePruning = HeuristikEdgePruning(this, edges.size(), constrainedEdgeDivider);
    }

    if(settings.polyanya.use_landmarks){
        Section section("init landmarks");
        initLandmarks();
        section.end_section();
    }

    if(settings.DEBUG_LEVEL > 2) assert(test_allCells_finite());
    if(settings.DEBUG_LEVEL > 1 && !settings.polyanya.merge_cells) assert(allTriangles());
    if(settings.DEBUG_LEVEL > 1) assert(test_all_cellsConvex(true) == 0);
};

inline const std::list<Coordinates> PolyAnya::getRoute(){
    std::list<Coordinates> route;
    for(auto n : getNodeIdOnRoute()){
        route.emplace_front(settings.use_merkartor ? rev_bounded_mercator(coords[n]) : coords[n]);
    }
    return route;
}

inline const std::list<Coordinates> PolyAnya::getRoute_withoutMercator(){
    std::list<Coordinates> route;
    for(auto n : getNodeIdOnRoute()){
        route.emplace_front(coords[n]);
    }
    return route;
}

inline const NT_Distance PolyAnya::evalutatePath(std::list<Coordinates> route){
    NT_Distance distance = 0.0;
    Coordinates last_point = route.front();
    for(auto &p : route){
        distance += getDist_exact(to_point(last_point), to_point(p));
        last_point = p;
    }
    return distance;
}

inline const Exact_Distance PolyAnya::evalutatePath_exact(std::list<Coordinates> route){
    Exact_Distance distance = 0.0;
    Coordinates last_point = route.front();
    for(auto &p : route){
        distance += getDist_exact_exact(Exact_Point(last_point.lat, last_point.lon), Exact_Point(p.lat, p.lon));
        last_point = p;
    }
    return distance;
}

inline const std::list<uint32_t> PolyAnya::getNodeIdOnRoute(){
    std::list<uint32_t> route;
    uint32_t currentNode = goalID;
    if(currentNode == UINT32_MAX || stats.distance == -1) return route;
    if(visited.size() == 0) {
        route.emplace_front(startID);
        route.emplace_back(goalID);
        return route;
    }
    while (previous_node[currentNode] != UINT32_MAX) {
        route.emplace_front(currentNode);
        assert(currentNode != previous_node[currentNode]);
        currentNode = previous_node[currentNode];
    }
    route.emplace_front(startID);
    return route;
}

inline const std::list<std::pair<Coordinates, Coordinates>> PolyAnya::getSearchSpace() {
    std::list<std::pair<Coordinates, Coordinates>> searchSpace;
    for (uint32_t e : visited){
        searchSpace.emplace_back(Coordinates(coords[firstNode(e << 1)]), Coordinates(coords[secondNode(e << 1)]));
        if(settings.use_merkartor) searchSpace.back().first = rev_bounded_mercator(searchSpace.back().first);
        if(settings.use_merkartor) searchSpace.back().second = rev_bounded_mercator(searchSpace.back().second);
    }
    return searchSpace;
}

inline bool PolyAnya::update_dist_to_node(uint32_t n, uint32_t generator, uint32_t edge_id){
    assert(edge_id < edges.size() || edge_id == UINT32_MAX);

    auto newDist = costs[generator] + getDist_exact_index(n, generator);
    uint32_t prev_n = previous_node[n];

    // #ifdef CGAL_USE_EPECK && !CGAL_USE_EPECK
    // if(previous_node[n] == generator && !has_generated[n]){
    //     costs[n] = newDist;
    //     return true;
    // }
    // // edge case if points are collinear check this instead of costly exact number s
    // if(prev_n != UINT32_MAX && costs[n] * 1.01 >= newDist && newDist * 1.01 >= costs[n] && collinear_index(generator, prev_n, n)){
    //     // I prev(g) == prev(n) || coll prev(g) prev(n), n
    //     uint32_t prev_prev_n = previous_node[prev_n], prev_g = previous_node[generator];
    //     if(prev_prev_n != UINT32_MAX && CGAL::collinear_are_ordered_along_line(coords[generator], coords[prev_n], coords[n])){
    //         if(prev_prev_n == generator || collinear_index(prev_prev_n, generator, n)){
    //             assert(abs(to_double(costs[n]) - to_double(newDist)) < 0.000001);
    //             return true;
    //         }
    //     }
    //     // II prev(prev(n)) == g || coll prev(prev(n)) g, n
    //     if(prev_g != UINT32_MAX && CGAL::collinear_are_ordered_along_line(coords[prev_n], coords[generator], coords[n])){
    //         if(prev_g == previous_node[n] || collinear_index(prev_g, prev_n, n)){
    //             assert(abs(to_double(costs[n]) - to_double(newDist)) < 0.000001);
    //             return true;
    //         }
    //     }
    // }
    // assert(to_double(costs[n]) * 1.0001 >= to_double(newDist) || previous_node[n] != generator);
    // if(previous_node[n] != generator && costs[n] < newDist) return false;
    // #else
    //new TODO:
    if(costs[n] < newDist) return false;
    if(costs[n] == newDist) return true;

    //old
    // if(costs[n] <= newDist * onePlusEps) return false;
    // #endif

    if(n != generator && previous_node[n] < previous_node.size()){
        if(settings.DEBUG_LEVEL > 1 && settings.VERBOSE_LEVEL > 2) 
            printf("changed gen from %u: from %u to %u\n", n, previous_node[n], generator);
    }
    previous_node[n] = generator;
    costs[n] = newDist;
    has_generated[n] = false;
    if(edge_id != UINT32_MAX) visitEdge(edge_id);
    return true;
}

inline void PolyAnya::processEdge(Wavelet & w, bool addRight, bool addLeft, std::vector<Wavelet>& waveletsToCreate){
    FaceCirculator iter(this, w.item); 
    uint32_t n0 = iter.firstNode(), n1 = iter.secondNode();
    iter.next();
    
    uint32_t newRightBound = n1, newLeftBound;
    // propagate right until no right turn
    while(iter.hasNext()){
        newLeftBound = iter.secondNode();
        
        if(w.generator == w.right_b || isTurn(LEFT, w.generator, w.right_b, coords[newLeftBound])){
            if(w.generator != w.right_b && CGAL::collinear(coords[w.generator], coords[w.right_b != UINT32_MAX ? w.right_b : previous_node[w.generator]], coords[newRightBound]) && !CGAL::collinear(coords[w.generator], coords[w.right_b != UINT32_MAX ? w.right_b : previous_node[w.generator]], coords[newLeftBound])){
                break;
            }
            newLeftBound = w.right_b;
        }
        if(addRight && newRightBound != newLeftBound){
            assert(previous_node[n1] == w.generator && w.right_b == n1);
            if(isDeadend(otherEdge(iter.current_e))){
                if(newLeftBound != w.right_b) update_dist_to_node(newLeftBound, n1, iter.current_e >> 1);
            }else{
                waveletsToCreate.emplace_back(waveletIds++, n1, otherEdge(iter.current_e), newLeftBound == n1 ? UINT32_MAX : newLeftBound, newRightBound, costs[n1], w.generator);
            }
        }
        newRightBound = newLeftBound;
        if(newLeftBound == w.right_b) break;
        iter.next();
    }
    
    // propagate mid until no right turn
    while(iter.hasNext()){
        newLeftBound = iter.secondNode();
        if(w.generator != w.left_b && !isTurn(RIGHT, w.generator, w.left_b, coords[newLeftBound])){
            newLeftBound = w.left_b;
        }
        if(newLeftBound != w.left_b) update_dist_to_node(newLeftBound, w.generator, iter.current_e >> 1);
        if(!isDeadend(otherEdge(iter.current_e))){
            assert(newRightBound != newLeftBound);
            auto left = newLeftBound;
            if(w.generator != w.left_b && CGAL::collinear(coords[w.left_b != UINT32_MAX ? w.left_b : previous_node[w.generator]], coords[iter.secondNode()], coords[w.generator])){
                left = iter.secondNode();
            }
            waveletsToCreate.emplace_back(waveletIds++, w.generator, otherEdge(iter.current_e), left, newRightBound, costs[w.generator], previous_node[w.generator]);
        }
        newRightBound = newLeftBound;
        if(newLeftBound == w.left_b) break;
        iter.next();
    }

    if(w.generator != w.left_b && CGAL::collinear(coords[w.left_b != UINT32_MAX ? w.left_b : previous_node[w.generator]], coords[iter.secondNode()], coords[w.generator])){
        iter.next();
        newRightBound = iter.firstNode();
    }
    //propagate left no breaking condition
    while(iter.hasNext() && addLeft){
        newLeftBound = iter.secondNode();
        if(addLeft && newRightBound != newLeftBound){
            assert(previous_node[n0] == w.generator);
            if(isDeadend(otherEdge(iter.current_e))){
                if(newLeftBound != n0) update_dist_to_node(newLeftBound, n0, iter.current_e >> 1);
            } else {
                waveletsToCreate.emplace_back(waveletIds++, n0, otherEdge(iter.current_e), newLeftBound, newRightBound == n0 ? UINT32_MAX : newRightBound, costs[n0], w.generator);
            }
        }
        newRightBound = newLeftBound;
        iter.next();
    }
}

inline void PolyAnya::processEdge_new(Wavelet & w, bool addRight, bool addLeft, std::vector<Wavelet>& waveletsToCreate){
    FaceCirculator iter(this, w.item); 
    uint32_t n0 = iter.firstNode(), n1 = iter.secondNode();
    iter.next();
    
    uint32_t newRightBound = n1, newLeftBound;
    bool foundRightTurn = false, foundLeftTurn = false;
    // propagate right until no right turn
    while(iter.hasNext()){
        if(!foundRightTurn){
            newLeftBound = iter.secondNode();
            if(w.generator == w.right_b || !isTurn(RIGHT, w.generator, w.right_b, coords[newLeftBound])){
                foundRightTurn = true;
                newLeftBound = w.right_b;
            }
            if(addRight && newRightBound != newLeftBound){
                assert(previous_node[n1] == w.generator && w.right_b == n1);
                if(isDeadend(otherEdge(iter.current_e))){
                    if(newLeftBound != w.right_b) update_dist_to_node(newLeftBound, n1, iter.current_e >> 1);
                }else{
                    waveletsToCreate.emplace_back(waveletIds++, n1, otherEdge(iter.current_e), newLeftBound == n1 ? UINT32_MAX : newLeftBound, newRightBound, costs[n1], w.generator);
                }
            }
            newRightBound = newLeftBound;
        }
        if(foundRightTurn && !foundLeftTurn){
            newLeftBound = iter.secondNode();
            if(w.generator != w.left_b && isTurn(LEFT, w.generator, w.left_b, coords[newLeftBound])){
                foundLeftTurn = true;
                newLeftBound = w.left_b;
            }
            //assert(newRightBound != newLeftBound);
            if(newRightBound != newLeftBound){
                if(isDeadend(otherEdge(iter.current_e))){
                    if(newLeftBound != w.left_b) update_dist_to_node(newLeftBound, w.generator, iter.current_e >> 1);
                }else{
                    waveletsToCreate.emplace_back(waveletIds++, w.generator, otherEdge(iter.current_e), newLeftBound, newRightBound, costs[w.generator], previous_node[w.generator]);
                }
            }
            newRightBound = newLeftBound;
        }
        if(foundLeftTurn){
            newLeftBound = iter.secondNode();
            if(addLeft && newRightBound != newLeftBound){
                if(isDeadend(otherEdge(iter.current_e))){
                    if(newLeftBound != n0) update_dist_to_node(newLeftBound, n0, iter.current_e >> 1);
                }else{
                    waveletsToCreate.emplace_back(waveletIds++, n0, otherEdge(iter.current_e), newLeftBound, newRightBound == n0 ? UINT32_MAX : newRightBound, costs[n0], w.generator);
                }
            }
            newRightBound = newLeftBound;
        }
        iter.next();
    }
}

inline void PolyAnya::addStartPoint(uint32_t start, uint32_t start_edge){
    assert(visited.size() == 0);
    assert(start_edge != UINT32_MAX);
    costs[start] = 0;
    has_generated[start] = true; stats.nodes_generated++;
    previous_node[start] == UINT32_MAX;
    // fill starting edges
    if(!hasNode(start_edge, start)){
        for(FaceCirculator f_iter(this, start_edge); f_iter.hasNext(); f_iter++){
            if(hasNode(f_iter.current_e, start)) start_edge = f_iter.current_e;
        }
    }
    assert(hasNode(start_edge, start));
    if(firstNode(start_edge) == start) start_edge = otherEdge(start_edge);
    std::vector<uint32_t> startingEdges;
    if(isNodeConstrained(start)) {
        while(nextEdge(otherEdge(start_edge)) != UINT32_MAX){
            assert(hasNode(start_edge, start));
            start_edge = get_predecessorEdge(otherEdge(start_edge));
        }
    }
    for(NodeCirculator n_iter(this, start_edge); n_iter.hasNext(); n_iter++){
        for(FaceCirculator f_iter(this, otherEdge(n_iter.current_e)); f_iter.hasNext(); f_iter++){
            if(hasNode(f_iter.current_e, start)) continue;
            startingEdges.push_back(otherEdge(f_iter.current_e));
        }
    }
    for(uint32_t edge : startingEdges){
        assert(edge != UINT32_MAX);
        if(costs[firstNode(edge)] < getDist_exact_index(firstNode(edge), start)) continue;
        previous_node[firstNode(edge)] = start;
        costs[firstNode(edge)] = getDist_exact_index(firstNode(edge), start);
        visitEdge(edge >> 1);

        if(isDeadend(edge)) continue; // skip deadends
        
        Wavelet wavelet(waveletIds++, start, edge, firstNode(edge), secondNode(edge), costs[start]);
        if(!settings.polyanya.use_edgePruning || edgePruning.insertWave(wavelet)){
            assert(count_visited[edge >> 1] != UINT32_MAX);
            // visitEdge(wavelet.item >> 1);
            pq.emplace(getDistance(wavelet), wavelet); stats.wavelets_pushed++;
            if(settings.DEBUG_LEVEL > 1 && settings.VERBOSE_LEVEL > 2){
                printf("\tDEBUG: push edge %u (%u %u); bounds: (%u,%u); gen=%u; g=%f f=%f\n", 
                    wavelet.item, firstNode(wavelet.item), secondNode(wavelet.item), 
                    wavelet.left_b, wavelet.right_b, wavelet.generator,
                    to_double(wavelet.cost_generator), to_double(getDistance(wavelet)));
            }
        }
    }
}

inline void PolyAnya::continuous_dijkstra(){
    std::vector<Wavelet> waveletsToCreate;
    double current_costs;
    Wavelet wavelet;

    // main loop
    while (!pq.empty() || waveletsToCreate.size() > 0) {
        if(settings.polyanya.skip_ahead && waveletsToCreate.size() == 1){
            assert(!settings.polyanya.use_edgePruning);
            wavelet = waveletsToCreate[0];
            current_costs = getDistance(wavelet);
            visitEdge(wavelet.item >> 1);
        }else{
            for(auto & w : waveletsToCreate){
                if(firstNode(w.item) != w.generator && secondNode(w.item) != w.generator) assert(!CGAL::left_turn(coords[firstNode(w.item)], coords[secondNode(w.item)], coords[w.generator]));
                if(!settings.polyanya.use_edgePruning || edgePruning.insertWave(w)){
                visitEdge(w.item >> 1);
                if(!settings.polyanya.use_edgePruning){
                    pq.emplace(getDistance(w), w); stats.wavelets_pushed++;
                }else{
                    pq.emplace(edgePruning.getDisanceHeuristik(w, goalID), w); stats.wavelets_pushed++;
                }
                if(settings.DEBUG_LEVEL > 1 && settings.VERBOSE_LEVEL > 2){
                printf("\tDEBUG: push edge %u (%u %u); bounds: (%u,%u); gen=%u; g=%f f=%f\n", 
                    w.item, firstNode(w.item), secondNode(w.item), 
                    w.left_b, w.right_b, w.generator,
                        to_double(w.cost_generator), to_double(getDistance(w)));
                    }
                }
            }

            if(pq.empty()) break;
            stats.maxPQSize = std::max(stats.maxPQSize, (uint32_t)pq.size());

            current_costs = pq.top().first;
            wavelet = pq.top().second;
            pq.pop(); stats.pq_pops++;

            if(settings.VERBOSE_LEVEL > 1 && stats.pq_pops % 100000 == 0){
                std::cout<< "DEBUG iteration: " << stats.pq_pops << ", current costs: " << current_costs << ", pq size: " << pq.size() << "\n";
            }
        }
        waveletsToCreate.clear();
        uint32_t generator = wavelet.generator;
        
        assert(wavelet.item != UINT32_MAX);
        assert(edges[wavelet.item >> 1].E[wavelet.item & 1] != UINT32_MAX);
        assert(!isDeadend(wavelet.item));
        assert(to_double(wavelet.cost_generator) >= to_double(costs[generator]));
        assert(wavelet.left_b != UINT32_MAX || wavelet.right_b != UINT32_MAX);
        assert(count_visited[wavelet.item >> 1] != UINT32_MAX);

        if(settings.DEBUG_LEVEL > 1 && settings.VERBOSE_LEVEL > 2){
            printf("DEBUG: pop edge %u (%u %u); bounds: (%u,%u); gen=%u; g=%f f=%f\n", 
                wavelet.item, firstNode(wavelet.item), secondNode(wavelet.item), 
                wavelet.left_b, wavelet.right_b, generator,
                to_double(wavelet.cost_generator), to_double(current_costs));
        }
        
        //duplicate detection aka outdated
        if(to_double(wavelet.cost_generator) > to_double(costs[generator] * onePlusEps)){
            continue;
        }
        if(wavelet.prev_generator != previous_node[generator]){
            if(settings.VERBOSE_LEVEL > 1) printf("pruned because gen %u changed from %u to %u\n", generator, wavelet.prev_generator, previous_node[generator]);
            continue;
        }

        if(settings.polyanya.use_edgePruning){
            if(!edgePruning.isWavePresent(wavelet)){
                if(settings.VERBOSE_LEVEL > 2) printf("pruned because edge is dominated\n");
                continue;
            }
            auto edgePruningDist = edgePruning.getDisanceHeuristik(wavelet, goalID);
            assert(edgePruningDist != MAX_DOUBLE);
            if(edgePruningDist > current_costs * 1.00001){
                pq.emplace(edgePruningDist, wavelet);
                continue;
            }
        }

        uint32_t n0 = firstNode(wavelet.item);
        uint32_t n1 = secondNode(wavelet.item);
        NT_Distance costs_secondNode = getDist_exact_index(n1, generator) + costs[generator];
        NT_Distance costs_firstNode = getDist_exact_index(n0, generator) + costs[generator];
        // simple pruning
        if(settings.polyanya.use_simple_pruning_heuristik){
            double edgeLength = getDist_inexact(coords[n1], coords[n0]);
            edgeLength *= 1.000001; // add small epsilon as compensation for rounding errors
            if(costs[n1] < costs_firstNode - edgeLength) continue;
            if(costs[n0] < costs_secondNode - edgeLength) continue;
        }

        if(settings.polyanya.heuristikEdgePruningProof){
            if(heuristikEdgePruning.isDominated(wavelet, costs_secondNode, costs_firstNode)){
                continue;
            }
        }

        if(settings.polyanya.use_advanced_pruning){
            if(previous_node[n1] != UINT32_MAX && previous_node[n1] != generator){
                bool pruneCondition = to_double(costs[n1]) * PRUNING_EPSILON_HEURISTIK < to_double(costs_secondNode);
                bool isRightN1N0 = CGAL::right_turn(coords[n1], coords[previous_node[n1]], coords[n0]);
                bool isRightGN0 = CGAL::left_turn(coords[previous_node[n1]], coords[generator], coords[n0]);
                bool isRightGN1 = CGAL::left_turn(coords[previous_node[n1]], coords[generator], coords[n1]);
                if(pruneCondition && isRightN1N0 && isRightGN0 && isRightGN1){
                    continue;
                }
            }
            if(previous_node[n0] != UINT32_MAX && previous_node[n0] != generator){
                bool pruneCondition = to_double(costs[n0]) * PRUNING_EPSILON_HEURISTIK < to_double(costs_firstNode);
                bool isRightN1N0 = CGAL::right_turn(coords[n1], coords[previous_node[n0]], coords[n0]);
                bool isRightN0G = CGAL::right_turn(coords[n0], coords[previous_node[n0]], coords[generator]);
                bool isRightN1G = CGAL::right_turn(coords[n1], coords[previous_node[n0]], coords[generator]);
                if(pruneCondition && isRightN1N0 && isRightN0G && isRightN1G){
                    continue;
                }   
            }
        }

        // goal test
        if (settings.polyanya.breakCondition && goalID != UINT32_MAX){
            if(!settings.polyanya.skip_ahead && current_costs >= costs[goalID] * 1.0001) break; //TODO maybe remove
            if(settings.polyanya.skip_ahead && (pq.empty() || pq.top().first >= costs[goalID]) && current_costs >= costs[goalID]) break;
        }

        // Actual Processing: scanning cell ccw and propagate right, straigt or left.
        //      just propagate left and right if possible see bool: addLeft & addRight
        assert(count_visited[wavelet.item >> 1] != UINT32_MAX);
        count_visited[wavelet.item >> 1]++;
        stats.wavelets_processed++;
        bool addRight = needTurn(n1, costs_secondNode, wavelet.right_b, generator);
        bool addLeft = needTurn(n0, costs_firstNode, wavelet.left_b, generator);

        if(settings.polyanya.ignore_collinear){
            processEdge_new(wavelet, addRight, addLeft, waveletsToCreate);
        } else {
            processEdge(wavelet, addRight, addLeft, waveletsToCreate);
        }

        for(auto & wtc : waveletsToCreate){
            if(settings.polyanya.heuristikEdgePruningProof){
                visitEdge(wtc.item >> 1);
                if(wtc.left_b != n0 && firstNode(wtc.item) == wtc.left_b){
                    auto dist = wtc.cost_generator + getDist_exact_index(wtc.generator, firstNode(wtc.item));
                    heuristikEdgePruning.insertCandidate(wavelet.item, dist, wtc.generator);
                }
                if(wtc.right_b != n1 && secondNode(wtc.item) == wtc.right_b){
                    auto dist = wtc.cost_generator + getDist_exact_index(wtc.generator, secondNode(wtc.item));
                    heuristikEdgePruning.insertCandidate(wavelet.item, dist, wtc.generator);
                }
            }
            
            if(!has_generated[wtc.generator]){
                stats.nodes_generated++;
                has_generated[wtc.generator] = true;
            }
        }
    }
}

std::pair<uint32_t, uint32_t> PolyAnya::prepareSearch(const Point & start_point, const Point & goal_point){
    resetSearch();
    assert(waveletIds == 0);
    stats.serach_id++;
    stats.start = settings.use_merkartor ? rev_bounded_mercator(start_point) : start_point;;
    stats.goal = settings.use_merkartor ? rev_bounded_mercator(goal_point) : goal_point;
    if(settings.DEBUG_LEVEL > 3) assert(test_all_cellsConvex(true) == 0);

    uint32_t nearestEdgeStart = getNearestEdge(start_point);
    uint32_t nearestEdgeGoal = getNearestEdge(goal_point);
    
    if(settings.DEBUG_LEVEL > 1 && nearestEdgeStart != UINT32_MAX){
        uint32_t dump_CGAL_ops = 0;
        assert(isInCell(start_point, nearestEdgeStart, dump_CGAL_ops));
        assert(isInCell(goal_point, nearestEdgeGoal, dump_CGAL_ops));
    }

    if(settings.VERBOSE_LEVEL > 2 && nearestEdgeStart != UINT32_MAX){
        for(FaceCirculator f_iter(this, nearestEdgeStart); f_iter.hasNext(); f_iter++){
            auto p = settings.use_merkartor ? rev_bounded_mercator(coords[f_iter.firstNode()]) : coords[f_iter.firstNode()];
            std::cout << "(" << p.x() << "," << p.y() << ") ";
        }
        auto p = settings.use_merkartor ? rev_bounded_mercator(start_point) : start_point;
        std::cout << ": point inside: (" << p.x() << "," << p.y() << ") \n";
    }
    if(settings.VERBOSE_LEVEL > 2 && nearestEdgeGoal != UINT32_MAX){
        for(FaceCirculator f_iter(this, nearestEdgeGoal); f_iter.hasNext(); f_iter++){
            auto p = settings.use_merkartor ? rev_bounded_mercator(coords[f_iter.firstNode()]) : coords[f_iter.firstNode()];
            std::cout << "(" << p.x() << "," << p.y() << ") ";
        }
        auto p = settings.use_merkartor ? rev_bounded_mercator(goal_point) : goal_point;
        std::cout << ": point inside: (" << p.x() << "," << p.y() << ") \n";
    }

    if(nearestEdgeGoal != UINT32_MAX && nearestEdgeStart != UINT32_MAX && nearestEdgeGoal != nearestEdgeStart){

        for(FaceCirculator f_iter(this, nearestEdgeStart); f_iter.hasNext(); f_iter++){ //TODO fix if no insert needed
            if(!edges[f_iter.current_e >> 1].isConstrained() && !edges[f_iter.current_e >> 1].hasNode(startID)){
                nearestEdgeStart = f_iter.current_e;
                break;
            }
        }
        for(FaceCirculator f_iter(this, nearestEdgeGoal); f_iter.hasNext(); f_iter++){
            if(!edges[f_iter.current_e >> 1].isConstrained() && !edges[f_iter.current_e >> 1].hasNode(goalID)){
                nearestEdgeGoal = f_iter.current_e;
                break;
            }
        }
        goalID = insertPoint(goal_point, nearestEdgeGoal, false);
        goalPointInserted = true;
        startID = insertPoint(start_point, nearestEdgeStart, true);
        startPointInserted = true;

        if(settings.DEBUG_LEVEL > 3) assert(test_allCells_finite());

        if(settings.polyanya.use_deadend_pruning) {
            assert(!edges[nearestEdgeStart >> 1].isConstrained());
            assert(!edges[nearestEdgeGoal >> 1].isConstrained());
            bool success = deadendFinder.markBeforeSearch(nearestEdgeStart, nearestEdgeGoal);
            if(!success){
                if(settings.VERBOSE_LEVEL > 0) 
                    printf("No route possible due to deadend pruning between %s and %s\n", to_string(start_point).c_str(), to_string(goal_point).c_str());
                nearestEdgeStart = UINT32_MAX;
                nearestEdgeGoal = UINT32_MAX;
            }
        }

        if(settings.polyanya.use_landmarks){
            landmarks.select_landmarks(startID, goalID, settings.polyanya.landmarks_to_use);
        }
    };

    return std::make_pair(nearestEdgeStart, nearestEdgeGoal);
}

inline void PolyAnya::routePolyanya(Point start_point, Point goal_point){
    assert(settings.polyanya.use_edgePruning == false || edgePruning.isSet());
    if(settings.use_merkartor) start_point = bounded_mercator(start_point);
    if(settings.use_merkartor) goal_point = bounded_mercator(goal_point);

    auto start_goal_edges = prepareSearch(start_point, goal_point);

    if(start_goal_edges.first == UINT32_MAX || start_goal_edges.second == UINT32_MAX){ // one point is not in bound -> no route
        if(settings.VERBOSE_LEVEL > 1) printf("One point is not in bound: %s or %s\n", to_string(start_point).c_str(), to_string(goal_point).c_str());
        stats.distance = -1.0;
    } else if(start_goal_edges.first == start_goal_edges.second){ // both points are in the same cell -> directly visible
        if(settings.VERBOSE_LEVEL > 1) printf("Both points are in the same cell: %s or %s\n", to_string(start_point).c_str(), to_string(goal_point).c_str());
        stats.distance = getDist_inexact(start_point, goal_point);
    } else {
        //search
        Section section("search route from " + to_string(start_point) + " to " + to_string(goal_point), 2);
        try{
            auto timer = Timer();
            uint32_t s_e = (edges.size() - 1) * 2;
            addStartPoint(startID, s_e);
            continuous_dijkstra();
            timer.stop_timer();
            stats.time_ms = timer.getTimeInMiliSec();
            if(costs[goalID] >= MAX_DOUBLE){
                if(settings.VERBOSE_LEVEL >= 0) printf("warning no route found! distance is -1\n");
                stats.distance = -1;
            }else{
                stats.distance = to_double(costs[goalID]);
                // Exact_Distance a = 0; //TODO remove
                // // a.
                // a
                if(settings.polyanya.print_out_exact_pathCosts) stats.distance = CGAL::to_double(evalutatePath_exact(getRoute_withoutMercator()));
            }
        }catch(const std::exception & e){
            std::cerr << "Error during continuous_dijkstra: " << e.what() << std::endl;
            return;
        }
        section.end_section();
    }
    if(settings.VERBOSE_LEVEL > 1) std::cout << "Memory usage: " << get_mem_usage() << " KB" << std::endl;
    printStats();
}

inline void PolyAnya::routePolyanya(uint32_t startID, uint32_t goalID){
    resetSearch();
    assert(startID < coords.size() && goalID < coords.size());
    assert(waveletIds == 0);
    stats.serach_id++;
    stats.start = settings.use_merkartor ? rev_bounded_mercator(coords[startID]) : coords[startID];;
    stats.goal = settings.use_merkartor ? rev_bounded_mercator(coords[goalID]) : coords[goalID];

    //search
    Section section("search route from " + to_string(coords[startID]) + " to " + to_string(coords[goalID]), 2);
    try{
        auto timer = Timer();
        addStartPoint(startID, getNearestEdge(coords[startID]));
        continuous_dijkstra();
        timer.stop_timer();
        stats.time_ms = timer.getTimeInMiliSec();
        if(costs[goalID] == MAX_DOUBLE){
            if(settings.VERBOSE_LEVEL >= 0) printf("warning no route found! distance is -1\n");
            stats.distance = -1;
        }else{
            stats.distance = to_double(costs[goalID]);
        }
    }catch(const std::exception & e){
        std::cerr << "Error during continuous_dijkstra: " << e.what() << std::endl;
        return;
    }
    section.end_section();

    printStats();
    if(settings.VERBOSE_LEVEL > 1) printf("Route has distance %lf\n", stats.distance);
}

inline std::vector<NT_Distance> PolyAnya::routeToAll(uint32_t start_point) {
    resetSearch();
    stats.serach_id++;

    bool old_breakCondition = settings.polyanya.breakCondition;
    settings.polyanya.breakCondition = false;
    bool old_use_deadend_pruning = settings.polyanya.use_deadend_pruning;
    settings.polyanya.use_deadend_pruning = false;
    bool old_edge_pruning = settings.polyanya.use_edgePruning;
    settings.polyanya.use_edgePruning = false;
    goalID = UINT32_MAX;

    costs[start_point] = 0;
    has_generated[start_point] = true;
    stats.nodes_generated++;

    //TODO: clean up code
    for(uint32_t e = 0; e < edges.size() * 2; e++){
        uint32_t start = firstNode(e);
        if(!has_generated[start]) continue;
        for(FaceCirculator f_iter(this, e); f_iter.hasNext(); ++f_iter){
            if((f_iter.current_e | 1) == UINT32_MAX || f_iter.secondNode() == start) break;
            if(costs[f_iter.secondNode()] > getDist_exact_index(f_iter.secondNode(), start)){
                previous_node[f_iter.secondNode()] = start;
                costs[f_iter.secondNode()] = getDist_exact_index(f_iter.secondNode(), start);
            }
            visitEdge(f_iter.current_e >> 1);
            if(f_iter.firstNode() == start) continue;
            if(isDeadend(f_iter.current_e)) continue;

            auto edge = otherEdge(f_iter.current_e);
            Wavelet wavelet(waveletIds++, start, edge, firstNode(edge), secondNode(edge), costs[start]);
            if(!settings.polyanya.use_edgePruning || edgePruning.insertWave(wavelet)){
                visitEdge(wavelet.item >> 1);
                pq.emplace(getDistance_no_heuristik(wavelet), wavelet); stats.wavelets_pushed++;
            }
        }
        
    }
    
    continuous_dijkstra();
    auto dists = costs;

    //reset
    for(uint32_t n = 0; n < coords.size(); n++){
        assert(costs[n] != MAX_DOUBLE || previous_node[n] == UINT32_MAX);
        //assert(costs[n] != MAX_DOUBLE);
        costs[n] = MAX_DOUBLE;
        previous_node[n] = UINT32_MAX;
        has_generated[n] = false;
    }
    resetSearch();
    //restore old settings
    settings.polyanya.breakCondition = old_breakCondition;
    settings.polyanya.use_deadend_pruning = old_use_deadend_pruning;
    settings.polyanya.use_edgePruning = old_edge_pruning;
    return dists;
};

void PolyAnya::resetSearch(){
    //reset stats
    stats.reset();
    waveletIds = 0;

    for(auto e : visited){
        count_visited[e] = UINT32_MAX;

        auto n0 = firstNode(e << 1);
        auto n1 = secondNode(e << 1);
        costs[n0] = MAX_DOUBLE;
        has_generated[n0] = false;
        previous_node[n0] = UINT32_MAX;
        costs[n1] = MAX_DOUBLE;
        has_generated[n1] = false;
        previous_node[n1] = UINT32_MAX;
    }
    edgePruning.resetAfterSearch(visited);
    heuristikEdgePruning.resetAfterSearch(visited);
    visited.clear();
    pq = PQ(compare);

    if(settings.polyanya.use_deadend_pruning) deadendFinder.resetAfterSearch();
    if(settings.polyanya.use_landmarks) landmarks.reset();

    if(startID != UINT32_MAX){
        costs[startID] = MAX_DOUBLE;
        has_generated[startID] = false;
        previous_node[startID] = UINT32_MAX;
        startID = UINT32_MAX;
    }

    if(startPointInserted) deleteLastPoint();
    startPointInserted = false;
    if(goalPointInserted) deleteLastPoint();
    goalPointInserted = false;

    // test proper reset
    if(settings.DEBUG_LEVEL > 1){
        for(auto c : costs) assert(c == MAX_DOUBLE);
        for(auto p : previous_node) assert(p == UINT32_MAX);
        assert(visited.size() == 0);
        for(auto g : has_generated) assert(g == false);
        for(auto c : count_visited) assert(c == UINT32_MAX);
    }
    if(settings.DEBUG_LEVEL > 3) assert(test_allCells_finite());
};

inline const bool PolyAnya::isDeadend(const uint32_t edge_id){
    assert((edge_id | 1) != UINT32_MAX);
    if((edge_id >> 1) < constrainedEdgeDivider) return true;
    return settings.polyanya.use_deadend_pruning ? deadendFinder.isDeadend(edge_id) : false;
};

inline const double PolyAnya::getDistance_no_heuristik(const Wavelet & w) const {
    const Point &n0 = coords[firstNode(w.item)];
    const Point &n1 = coords[secondNode(w.item)];
    return to_double(w.cost_generator) + CGAL::sqrt(to_double(CGAL::squared_distance(coords[w.generator], K::Segment_2 (n0, n1))));
}

// like in polyanya paper:
//      if left Turn: g -> left -> goal
//      if right Turn: g -> right -> goal
// else g -> goal
const double PolyAnya::getDistance(const Wavelet & w) const {
    assert(goalID == UINT32_MAX || getDistance_no_heuristik(w) <= getDistance_simpler(w) + 1e-9);
    if(settings.polyanya.no_heuristik || goalID == UINT32_MAX) return getDistance_no_heuristik(w);
    if(settings.polyanya.use_landmarks){
        assert(goalID == UINT32_MAX || getDistance_no_heuristik(w) <= landmarks.getHeuristik(w) + 1e-9);
        return landmarks.getHeuristik(w);
    } 
    // if(settings.polyanya.use_landmarks) return std::max(landmarks.getHeuristik(w), getDistance_simpler(w));
    if(settings.polyanya.use_simple_distance_heuristik) return getDistance_simpler(w);
    
    const Point &n0 = coords[firstNode(w.item)];
    const Point &n1 = coords[secondNode(w.item)];
    const Point &g = coords[w.generator];
    Point goal = coords[goalID];
    Point p = g;

    if(CGAL::right_turn(n0, n1, goal)){
        goal = Transformation(CGAL::REFLECTION, Line_2(n0, n1))(goal);
    }

    double distance = to_double(costs[w.generator]);      
    if(isTurn(LEFT, w.generator, w.left_b, goal)){
        if(w.left_b == firstNode(w.item)){
            p = n0;
        }else{
            p = calc_intersection(g, coords[w.left_b != UINT32_MAX ? w.left_b : previous_node[w.generator]], n0, n1);
        }
    }else if(isTurn(RIGHT, w.generator, w.right_b, goal)){
        if(w.right_b == secondNode(w.item) ){
            p = n1;
        }else{
            p = calc_intersection(g, coords[w.right_b != UINT32_MAX ? w.right_b : previous_node[w.generator]], n0, n1);
        }
    }
    if(settings.DEBUG_LEVEL > 1){
        double dist_no_heuristik = getDistance_no_heuristik(w);
        double dist_simpler = getDistance_simpler(w);
        double dist_current = distance + getDist_inexact(goal, p) + getDist_inexact(g, p);
        assert(dist_no_heuristik <= dist_current * (1 + 1e-9));
        // assert(dist_simpler <= dist_current * (1 + 1e-9));
    }
    return distance + getDist_inexact(goal, p) + getDist_inexact(g, p);
}

// // gets distance via trianglei inequality to save CGAL ops
// inline const double PolyAnya::getDistance_simpler(const Wavelet & w){
//     const Point &n0 = coords[firstNode(w.item)], &n1 = coords[secondNode(w.item)];
//     const Point &g = coords[w.generator], &goal = coords[goalID];
//     double edgeLength = getDist_inexact(n0, n1);
//     double distOvern0 = getDist_inexact(goal, n0) + (w.left_b == UINT32_MAX ? getDist_inexact(n0, g) : getDist_inexact(n0, coords[w.left_b]) + getDist_inexact(coords[w.left_b], g));
//     double distOvern1 = getDist_inexact(goal, n1) + (w.right_b == UINT32_MAX ? getDist_inexact(n1, g) : getDist_inexact(n1, coords[w.right_b]) + getDist_inexact(coords[w.right_b], g));
//     double hValue = std::max(std::min(distOvern0, distOvern1) - (edgeLength), std::max(distOvern0, distOvern1) - 2 * edgeLength);
//     hValue = std::max(hValue, 0.0);
//     return to_double(costs[w.generator]) + hValue;
// };

inline const double PolyAnya::getDistance_simpler(const Wavelet & w) const {
    const Point &n0 = coords[firstNode(w.item)];
    const Point &n1 = coords[secondNode(w.item)];
    const Point &g = coords[w.generator];
    Point goal = coords[goalID];
    Point p = g;

    if(CGAL::right_turn(n0, n1, goal)){
        goal = Transformation(CGAL::REFLECTION, Line_2(n0, n1))(goal); 
    }

    double distance = to_double(costs[w.generator]);     
    if(isTurn(LEFT, w.generator, w.left_b, goal)){
        if(w.left_b != UINT32_MAX) p = coords[w.left_b];
        if(w.left_b != UINT32_MAX && !isTurn(RIGHT, w.left_b, firstNode(w.item), goal)) {
            return distance + getDist_inexact(goal, n0) + getDist_inexact(n0, p) + getDist_inexact(p, g);
        }
        if(w.left_b == UINT32_MAX && !isTurn(RIGHT, w.generator, firstNode(w.item), goal)){
            return distance + getDist_inexact(goal, n0) + getDist_inexact(n0, g);
        }
    }else if(isTurn(RIGHT, w.generator, w.right_b, goal)){
        if(w.right_b != UINT32_MAX) p = coords[w.right_b];
        if(w.right_b != UINT32_MAX && !isTurn(LEFT, w.right_b, secondNode(w.item), goal)) {
            return distance + getDist_inexact(goal, n1) + getDist_inexact(n1, p) + getDist_inexact(p, g);
        }
        if(w.right_b == UINT32_MAX && !isTurn(LEFT, w.generator, secondNode(w.item), goal)){
            return distance + getDist_inexact(goal, n1) + getDist_inexact(n1, g);
        }
    }
    return distance + getDist_inexact(goal, p) + getDist_inexact(p, g);
};

inline uint32_t PolyAnya::insertPoint(const Point & p, uint32_t first_edge_id, bool is_start_point = false){
    coords.push_back(p);
    costs.push_back(MAX_DOUBLE);
    has_generated.push_back(false);
    previous_node.push_back(UINT32_MAX);
    uint32_t newNodeindex = coords.size() - 1;

    if(settings.polyanya.use_landmarks) landmarks.insertPoint(newNodeindex, first_edge_id, !is_start_point);

    if(first_edge_id == UINT32_MAX){
        if(settings.VERBOSE_LEVEL > 1) std::cerr << "Point is not in bound: " << p.x() << ", " << p.y() << std::endl;
    } else{
        uint32_t lastNewEdge = UINT32_MAX;
        for(FaceCirculator f_iter(this, first_edge_id); f_iter.hasNext(); ++f_iter){
            edges.emplace_back(newNodeindex, f_iter.secondNode(), lastNewEdge, nextEdge(f_iter.current_e));
            count_visited.emplace_back(UINT32_MAX);
            lastNewEdge = (2 * (edges.size() - 1)) ^ 1;
            edges[f_iter.current_e >> 1].set_next_e(f_iter.secondNode(), lastNewEdge ^ 1);
            f_iter.current_e = lastNewEdge;
        }

        // link first new edge with last 
        edges[nextEdge(first_edge_id) >> 1].set_next_e(newNodeindex, lastNewEdge);
    }
    return newNodeindex;
};

inline void PolyAnya::deleteLastPoint(){
    uint32_t id = coords.size() - 1, first_e = (edges.size() - 1) * 2;
    coords.pop_back();
    costs.pop_back();
    has_generated.pop_back();
    previous_node.pop_back();

    assert(coords.size() == id);
    assert(costs.size() == id);
    assert(has_generated.size() == id);
    assert(previous_node.size() == id);

    // nothing to do here
    if(!edges.back().hasNode(id)) return;
    // restore links
    for(NodeCirculator n_iter(this, first_e); n_iter.hasNext(); n_iter++){
        assert(n_iter.firstNode() == id || n_iter.secondNode() == id);
        uint32_t cw_edge = nextEdge(nextEdge(n_iter.current_e));
        assert(edges[cw_edge >> 1].E[cw_edge & 1] >> 1 == n_iter.current_e >> 1);
        assert(edges[cw_edge >> 1].E[cw_edge & 1] == n_iter.current_e);
        assert(firstNode(nextEdge(otherEdge(n_iter.current_e))) != id);
        assert(secondNode(nextEdge(otherEdge(n_iter.current_e))) != id);
        edges[cw_edge >> 1].E[cw_edge & 1] = nextEdge(otherEdge(n_iter.current_e));
    }
    // delete edges to deleted Node
    while(edges.back().hasNode(id)){
        edges.pop_back();
        count_visited.pop_back();
    }

    if(settings.polyanya.use_landmarks) landmarks.deleteLastPoint();
    if(settings.DEBUG_LEVEL > 0 && settings.polyanya.use_landmarks) assert(landmarks.getNodeCount() == coords.size());
};

inline const bool PolyAnya::isTurn(bool dir, uint32_t n0, uint32_t n1, Point & test) const {
    if(n0 == n1) return false;
    if(n1 == UINT32_MAX){
        n1 = previous_node[n0];
        std::swap(n0, n1);
    }
    if(n0 == UINT32_MAX){
        n0 = previous_node[n1];
        std::swap(n0, n1);
    }
    return dir == LEFT ? CGAL::left_turn(coords[n0], coords[n1], test) : CGAL::right_turn(coords[n0], coords[n1], test);
};

#endif