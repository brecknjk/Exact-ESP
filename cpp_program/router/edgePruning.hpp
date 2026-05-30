#ifndef edgePruning_h
#define edgePruning_h

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

#include "graph.hpp"

constexpr double PRUNING_EPSILON = 0.00000001;
constexpr double onePlusEps = 1 + PRUNING_EPSILON;
constexpr double oneMinusEps = 1 - PRUNING_EPSILON;

inline void filterAndSortIntersections(std::list<Coordinates> & intersections, const Coordinates & left_b, const Coordinates & right_b){
    assert(intersections.size() <= 2);
    if(intersections.size() == 0) return;
    //filter results after intervall
    double length_sqrd = left_b.squared_distance(right_b) * 1.0001;
    auto r_it = intersections.begin();
    if(intersections.size() == 1){
        left_b.squared_distance(*r_it) > length_sqrd || right_b.squared_distance(*r_it) > length_sqrd ? r_it = intersections.erase(r_it) : ++r_it;
    }else{
        left_b.squared_distance(*r_it) > length_sqrd || right_b.squared_distance(*r_it) > length_sqrd ? r_it = intersections.erase(r_it) : ++r_it;
        left_b.squared_distance(*r_it) > length_sqrd || right_b.squared_distance(*r_it) > length_sqrd ? r_it = intersections.erase(r_it) : ++r_it;
        //sort insections along edge
        if(intersections.size() > 1 && left_b.squared_distance(intersections.front()) > left_b.squared_distance(intersections.back())){
            std::swap(intersections.front(), intersections.back());
        }
    }
}

inline std::list<Coordinates> getHyperbelIntersection(const Coordinates & g0, const Coordinates & g1, double difference, const Coordinates & n0, const Coordinates & n1, bool filter = false){
    std::list<Coordinates> result;

    double e = g0.distance(g1) / 2.0;
    double a_sqrd = difference / 2.0 * difference / 2.0;
    double b_sqrd = e*e - a_sqrd;

    Coordinates translate((g1.lat + g0.lat) / 2, (g1.lon + g0.lon) / 2);

    double cos_alpha = (g1.lat - g0.lat) / (2*e);
    double sin_alpha = (g1.lon - g0.lon) / (2*e);
    Matrix2x2 rotate(cos_alpha, -sin_alpha, sin_alpha, cos_alpha);
    Matrix2x2 rev_rotate(cos_alpha, sin_alpha, -sin_alpha, cos_alpha);

    auto n0_trans = rotate * (n0 - translate);
    auto n2_trans_DEBUG = rotate * (n1 - translate);
    auto diff_trans = rotate * (n1 - translate) - n0_trans;

    if(a_sqrd == 0){
        // edge case hyperbel is line
        if(diff_trans.lat == 0 && n0_trans.lat == 0) return result;
        assert(diff_trans.lat != 0);
        auto x = - n0_trans.lat / diff_trans.lat;
        result.emplace_back(rev_rotate * Coordinates(diff_trans.lat * x + n0_trans.lat, diff_trans.lon * x + n0_trans.lon) + translate);
    }else{
        //mitternachts formel
        auto a_temp = diff_trans.lat*diff_trans.lat/a_sqrd - diff_trans.lon*diff_trans.lon/b_sqrd;
        auto b_temp = 2 * (diff_trans.lat*n0_trans.lat/a_sqrd - diff_trans.lon*n0_trans.lon/b_sqrd);
        auto c_temp = n0_trans.lat*n0_trans.lat/a_sqrd - n0_trans.lon*n0_trans.lon/b_sqrd - 1;
        auto discriminante = b_temp*b_temp - 4*a_temp*c_temp;
        if(discriminante >= 0){
            auto x = (-b_temp + sqrt(discriminante)) / (2*a_temp);
            result.emplace_back(rev_rotate * Coordinates(diff_trans.lat * x + n0_trans.lat, diff_trans.lon * x + n0_trans.lon) + translate);
        }
        if(discriminante > 0){
            auto x = (-b_temp - sqrt(discriminante)) / (2*a_temp);
            result.emplace_back(rev_rotate * Coordinates(diff_trans.lat * x + n0_trans.lat, diff_trans.lon * x + n0_trans.lon) + translate);
        }
    }
    if(filter){
        for(auto it = result.begin(); it != result.end();){
            it = abs(g0.distance(*it) + difference - g1.distance(*it)) < 0.00001 ? ++it : result.erase(it);
        }
    }
    return result;
}

// test function to test if getHyperbelIntersection && filterAndSortIntersections works as intended
bool test_HyperbelIntersection(){
    std::list<Coordinates> result;

    // no intersection
    result = getHyperbelIntersection(Coordinates(-2,0), Coordinates(3,0), 3, Coordinates(0, -5), Coordinates(0,5));
    assert(result.size() == 0);

    //horizontal line
    result = getHyperbelIntersection(Coordinates(-2,0), Coordinates(3,0), 3, Coordinates(-2, 2), Coordinates(3,2));
    assert(result.size() == 2);
    filterAndSortIntersections(result, Coordinates(-2, 2), Coordinates(3,2));
    assert(result.size() == 2); // not filtered
    assert(result.front().distance(Coordinates(-1.63, 2)) < 0.01);
    assert(result.back().distance(Coordinates(2.62, 2)) < 0.01);
    std::swap(result.front(), result.back());
    filterAndSortIntersections(result, Coordinates(-2, 2), Coordinates(3,2));
    assert(result.size() == 2); // not filtered
    assert(result.front().distance(Coordinates(-1.63, 2)) < 0.01);
    assert(result.back().distance(Coordinates(2.62, 2)) < 0.01);
    filterAndSortIntersections(result, Coordinates(-1, 2), Coordinates(2,2));
    assert(result.size() == 0); // all filtered

    //vertical line
    result = getHyperbelIntersection(Coordinates(-2,0), Coordinates(3,0), 3, Coordinates(3, 2), Coordinates(3,-5));
    assert(result.size() == 2);
    filterAndSortIntersections(result, Coordinates(3, 2), Coordinates(3,-5));
    assert(result.size() == 1); // one filtered
    assert(result.front().distance(Coordinates(3, -2.67)) < 0.01);

    //collinear line
    result = getHyperbelIntersection(Coordinates(-2,0), Coordinates(3,0), 3, Coordinates(-3,0), Coordinates(4,0));
    assert(result.size() == 2);
    filterAndSortIntersections(result, Coordinates(-3,0), Coordinates(4,0));
    assert(result.size() == 2);
    assert(result.front().distance(Coordinates(-1, 0)) < 0.01);
    assert(result.back().distance(Coordinates(2, 0)) < 0.01);

    result = getHyperbelIntersection(Coordinates(-2,1), Coordinates(2,0), 3.0, Coordinates(-14, 18), Coordinates(10,6));
    assert(result.size() == 2);
    filterAndSortIntersections(result, Coordinates(-14, 18), Coordinates(8,7));
    assert(result.size() == 1);
    assert(result.front().distance(Coordinates(-10, 16)) < 0.01);
    return true;
}

#ifdef CGAL_USE_EPECK
class EdgePruning{
public:
    EdgePruning(){}
    EdgePruning(Graph* g, uint32_t edgeCount, uint32_t edgeOffset){
        assert(!settings.polyanya.use_edgePruning);
        printf("EdgePruning is disabled!\n");
    }
    bool isSet(){return false;};
    bool insertWave(const Wavelet & w){return true;}
    bool isWavePresent(const Wavelet & w){return true;}
    void resetAfterSearch(std::list<uint32_t> & visited){};
    double getDisanceHeuristik(const Wavelet & w, uint32_t goalID){return 0;};
};
#else
class EdgePruning{
protected:
    struct Wave{
        uint32_t id, generator;
        Coordinates left_bound, right_bound;
        double cost_generator = 0.0;
        Wave(uint32_t id, uint32_t generator, Coordinates left_bound, Coordinates right_bound, double cost_generator) : id(id), generator(generator), left_bound(left_bound), right_bound(right_bound), cost_generator(cost_generator){};
        void printWave(Graph* g){
            if(generator == UINT32_MAX){
                auto l = settings.use_merkartor ? rev_bounded_mercator(left_bound) : left_bound;
                auto r = settings.use_merkartor ? rev_bounded_mercator(right_bound) : right_bound;
                printf("wavefront: gen: UINT32_MAX, dist gen: %lf, left(%lf,%lf), right(%lf,%lf)\n", cost_generator, l.lat, l.lon, r.lat, r.lon);
            }else{
                auto gen = settings.use_merkartor ? rev_bounded_mercator(Coordinates(g->coords[generator].x(), g->coords[generator].y())) : Coordinates(g->coords[generator].x(), g->coords[generator].y());
                auto l = settings.use_merkartor ? rev_bounded_mercator(left_bound) : left_bound;
                auto r = settings.use_merkartor ? rev_bounded_mercator(right_bound) : right_bound;
                printf("wavefront: gen: (%lf,%lf), dist gen: %lf, left(%lf,%lf), right(%lf,%lf)\n", gen.lat, gen.lon, cost_generator, l.lat, l.lon, r.lat, r.lon);
            }
            
        }
    };
    Graph* g = nullptr;
    std::vector<std::list<Wave>> wavesOnEdges;
    uint32_t edgeOffset = 0, debugMaxWavefrontsOnEdge = 0, debug_MaxSameWavefrontOnEdge = 0;

public:
    EdgePruning(){}
    EdgePruning(Graph* g, uint32_t edgeCount, uint32_t edgeOffset) : g(g), edgeOffset(edgeOffset){
        wavesOnEdges.resize(edgeCount - edgeOffset);
        if(settings.VERBOSE_LEVEL > 1) printf("EdgePruning initialized with %i edges and offset %i\n", edgeCount - edgeOffset, edgeOffset);
        if(settings.DEBUG_LEVEL > 0) assert(test_HyperbelIntersection());
    }
    bool isBetterAtPoint(Wave &old_w, Wave& new_w, Coordinates & test, double error = 0.00001){
        return getDiffrenceAtPoint(old_w, new_w, test) > error;
    }
    double getDiffrenceAtPoint(Wave &old_w, Wave& new_w, Coordinates & test){
        if(old_w.generator == UINT32_MAX) return MAX_DOUBLE;
        if(new_w.generator == UINT32_MAX) return -MAX_DOUBLE;
        auto dist_old_w = getDist_inexact(g->coords[old_w.generator], Point(test.lat, test.lon));
        auto dist_new_w = getDist_inexact(g->coords[new_w.generator], Point(test.lat, test.lon));
        return (dist_old_w + old_w.cost_generator) - (dist_new_w + new_w.cost_generator);
    }
    bool isSet(){return g != nullptr && wavesOnEdges.size() > 0;};
    bool edgeIsSaved(uint32_t e){
        return (e >> 1) - edgeOffset >= 0 && (e >> 1) - edgeOffset < wavesOnEdges.size();
    };
    uint32_t map_edge(uint32_t e){
        assert(edgeIsSaved(e));
        return (e >> 1) - edgeOffset;
    };
    uint32_t rev_map_edge(uint32_t e){return (e + edgeOffset) << 1;};
    void reset(std::vector<uint32_t> & visited){
        for(auto e : visited){
            wavesOnEdges[map_edge(e)].clear();
        }
    }
    void resetAfterSearch(std::list<uint32_t> & visited){
        if(!isSet()) return;
        for(auto & e : visited) {
            if(edgeIsSaved(e << 1)) wavesOnEdges[map_edge(e << 1)].clear();
        }
        if(settings.DEBUG_LEVEL > 1) assert(debug_testProperReset());
    }
    bool isWavePresent(const Wavelet & w){
        if(!edgeIsSaved(w.item)) return true;
        for(auto & wave : wavesOnEdges[map_edge(w.item)]){
            if(wave.id == w.id) return true;
        }
        return false;
    }
    double getDisanceHeuristik(const Wavelet & w, uint32_t goalID){
        if(!edgeIsSaved(w.item)) return true;
        double best_distance = MAX_DOUBLE;
        for(auto & wave : wavesOnEdges[map_edge(w.item)]){
            if(wave.id != w.id) continue;
            Point n0 = to_point(wave.left_bound);
            Point n1 = to_point(wave.right_bound);
            const Point &generator = g->coords[w.generator];
            if(CGAL::left_turn(n0, n1, generator)){
                std::swap(n0, n1);
            }
            Point goal = g->coords[goalID];
            Point p = generator;

            if(CGAL::right_turn(n0, n1, goal)){
                goal = Transformation(CGAL::REFLECTION, Line_2(n0, n1))(goal);
            }
            if(CGAL::left_turn(generator, n0, goal)){
                best_distance = std::min(best_distance, wave.cost_generator + getDist_inexact(generator, n0) + getDist_inexact(goal, n0));
            } else if(CGAL::right_turn(generator, n1, goal)){
                best_distance = std::min(best_distance, wave.cost_generator + getDist_inexact(generator, n1) + getDist_inexact(goal, n1));
            } else{
                best_distance = std::min(best_distance, wave.cost_generator + getDist_inexact(generator, goal));
            }
        }
        return best_distance;
    }
    //return if wavefront could be inserted aka if not dominated
    bool insertWave(const Wavelet & w){
        if(!edgeIsSaved(w.item)) return true;
        auto e = map_edge(w.item);

        Point n0 = g->coords[g->firstNode(rev_map_edge(e))];
        Point n1 = g->coords[g->secondNode(rev_map_edge(e))];
        Coordinates n0_c(n0);
        Coordinates n1_c(n1);

        Coordinates newLeft_b = Coordinates(calc_intersection(g->coords[w.generator], g->coords[w.left_b != UINT32_MAX ? w.left_b : w.prev_generator], n0, n1));
        Coordinates newRight_b = Coordinates(calc_intersection(g->coords[w.generator], g->coords[w.right_b != UINT32_MAX ? w.right_b : w.prev_generator], n0, n1));
        if(w.left_b == g->firstNode(w.item)) newLeft_b = Coordinates(g->coords[w.left_b]);
        if(w.right_b == g->secondNode(w.item)) newRight_b = Coordinates(g->coords[w.right_b]);

        assert(w.left_b != g->secondNode(w.item));
        assert(w.right_b != g->firstNode(w.item));
        if(CGAL::collinear(n0, n1, g->coords[w.generator])){ // edge case if generator is collinear -> take closest bound
            if(w.generator != g->firstNode(rev_map_edge(e)) && w.generator != g->secondNode(rev_map_edge(e))){
                return false;
            }
            assert(w.left_b != UINT32_MAX && w.right_b != UINT32_MAX);
            newLeft_b = Coordinates(g->coords[w.left_b]);
            newRight_b = Coordinates(g->coords[w.right_b]);
        }
        if((w.item & 1) == 1) std::swap(newLeft_b, newRight_b);
        
        if(wavesOnEdges[e].size() == 0) wavesOnEdges[e].emplace_back(UINT32_MAX, UINT32_MAX, Coordinates(n0), Coordinates(n1), MAX_DOUBLE);

        //for left -> right: 
        // if new_I lies completely right -> continue
        // if new_I lies completely left:
            //if foundLeftbound: add intervall && set to false else continue
        // -> intersects!
        // if it lies completly in new_I and is better -> dominated -> remove it
        // if it starts left of new_I -> set right bound to new_I.left
        // if it ends right of new_I -> set left bound to new_I.right

        bool foundLeftbound = false;
        for (auto it = wavesOnEdges[e].begin(); it != wavesOnEdges[e].end(); ++it){
            if(it->id == 571026){
                auto temp = wavesOnEdges[e];
                // settings.VERBOSE_LEVEL = 4;
            }
            if(it->id == w.id) continue;
            if(!foundLeftbound && n1_c.squared_distance(it->right_bound) * onePlusEps >= n1_c.squared_distance(newLeft_b)){
                //new_I lies completely right
                continue;
            }

            if(n1_c.squared_distance(it->left_bound) <= n1_c.squared_distance(newRight_b) * onePlusEps){
                //new_I lies completely left
                if(foundLeftbound){
                    it = wavesOnEdges[e].insert(it, Wave(w.id, w.generator, newLeft_b, newRight_b, w.cost_generator));
                    it++;
                    foundLeftbound = false;
                }
                break;
            }
            Wave itCopy = *it;
            //assert(n1_c.squared_distance(it->left_bound) <= n1_c.squared_distance(newLeft_b) * onePlusEps == isOrdered_onLine(n1_c, it->left_bound, newLeft_b, onePlusEps));
            bool new_is_left_bigger = n1_c.squared_distance(it->left_bound) <= n1_c.squared_distance(newLeft_b) * onePlusEps;
            bool new_is_right_bigger = n0_c.squared_distance(it->right_bound) <= n0_c.squared_distance(newRight_b) * onePlusEps;
            Coordinates left = new_is_left_bigger ? it->left_bound : newLeft_b;
            Coordinates right = new_is_right_bigger ? it->right_bound : newRight_b;
            //not inbound?
            if(n0_c.squared_distance(right) - n0_c.squared_distance(left) < -PRUNING_EPSILON * 0.001) continue;
            assert(n0_c.squared_distance(right) - n0_c.squared_distance(left) >= -PRUNING_EPSILON * 0.001);
            if(abs(n0_c.squared_distance(right) - n0_c.squared_distance(left)) < PRUNING_EPSILON * 0.001){
                if(foundLeftbound && !new_is_right_bigger){
                    it = wavesOnEdges[e].insert(it, Wave(w.id, w.generator, newLeft_b, newRight_b, w.cost_generator));
                    it++;
                    foundLeftbound = false;
                    break;
                }
                if(!foundLeftbound) continue;
            }
            assert(n0_c.squared_distance(right) - n0_c.squared_distance(left) >= 0);
            if(n0_c.squared_distance(right) - n0_c.squared_distance(left) < 0) continue;

            double difference = it->cost_generator - w.cost_generator;
            bool isBetterLeft = it->generator == UINT32_MAX ? true : left.distance(g->coords[w.generator]) < left.distance(g->coords[it->generator]) + difference;
            
            //bool isBetterRight = getDist_inexact(g->coords[w.generator], to_point(right)) + difference < getDist_inexact(to_point(g->coords[it->generator]), to_point(right));
            
            //->intersects?
            std::list<Coordinates> result = (it->generator == UINT32_MAX) ? std::list<Coordinates>() : getHyperbelIntersection(g->coords[it->generator], g->coords[w.generator], difference, n0, n1, true);
            filterAndSortIntersections(result, left, right);

            if(isBetterLeft && !new_is_left_bigger){
                it = wavesOnEdges[e].insert(it, Wave(it->id, it->generator, it->left_bound, left, it->cost_generator));
                it++;
                it->left_bound = left;
                // continue;
                // assert(false); //can happen?
            }

            if(foundLeftbound && !new_is_left_bigger){
                it = wavesOnEdges[e].insert(it, Wave(w.id, w.generator, newLeft_b, left, w.cost_generator));
                it++;
                foundLeftbound = false;
                break;
            }

            if(isBetterLeft && !foundLeftbound){
                newLeft_b = left;
                foundLeftbound = true;
            }
            if(!isBetterLeft && foundLeftbound){
                assert(new_is_left_bigger);
                it = wavesOnEdges[e].insert(it, Wave(w.id, w.generator, newLeft_b, left, w.cost_generator));++it;
                newLeft_b = left;
                it->left_bound = left;
                foundLeftbound = false;
            }
            if(isBetterLeft || foundLeftbound){
                it->left_bound = result.size() == 0 ? newRight_b : result.front();
                if(result.size() != 0){
                    //insert wave with leftbound = newLeft_b, rightbound = result.front()
                    it = wavesOnEdges[e].insert(it, Wave(w.id, w.generator, newLeft_b, result.front(), w.cost_generator));++it;
                    if(result.size() == 2){
                        it->right_bound = result.back();
                        newLeft_b = result.back();
                        if(!new_is_right_bigger){
                            it = wavesOnEdges[e].insert(++it, Wave(w.id, w.generator, newLeft_b, right, w.cost_generator));
                            it = wavesOnEdges[e].insert(++it, Wave(itCopy.id, itCopy.generator, right, itCopy.right_bound, itCopy.cost_generator));
                            foundLeftbound = false;
                            break;
                        }
                    }else{
                        foundLeftbound = false;
                    }
                    
                }
            }else{
                newLeft_b = result.size() == 0 ? right : result.front();
                assert(!foundLeftbound);
                if(result.size() != 0){
                    it->right_bound = result.front();
                    if(result.size() == 2){
                        it = wavesOnEdges[e].insert(++it, Wave(w.id, w.generator, result.front(), result.back(), w.cost_generator)); ++it;
                        it = wavesOnEdges[e].insert(it, Wave(itCopy.id, itCopy.generator, result.back(), itCopy.right_bound, itCopy.cost_generator));
                        newLeft_b = right;
                    }else{
                        foundLeftbound = true;
                    }
                }
            }

            if(foundLeftbound && isBetterLeft && !new_is_right_bigger){
                //add intervall with leftbound = right, rightbound = old_it->right_bound
                it = ++wavesOnEdges[e].insert(it, Wave(w.id, w.generator, newLeft_b, right, w.cost_generator));
                foundLeftbound = false;
                break;
            }

            if(foundLeftbound && !isBetterLeft && !new_is_right_bigger){
                //add intervall with leftbound = right, rightbound = old_it->right_bound
                it = wavesOnEdges[e].insert(++it, Wave(w.id, w.generator, newLeft_b, right, w.cost_generator));
                it = wavesOnEdges[e].insert(++it, Wave(itCopy.id, itCopy.generator, right, itCopy.right_bound, itCopy.cost_generator));
                foundLeftbound = false;
                break;
            }
        }

        if(foundLeftbound) wavesOnEdges[e].emplace_back(w.id, w.generator, newLeft_b, newRight_b, w.cost_generator);

        auto temp = wavesOnEdges[e];

        
        //afterward iterate over list and delete all empty intervalls
        for(auto it = wavesOnEdges[e].begin(); it != wavesOnEdges[e].end(); ){
            if(it->left_bound.distance(n0) > it->right_bound.distance(n0)) it->left_bound = it->right_bound;
            auto length = it->left_bound.distance(it->right_bound);
            if(it->left_bound.distance(n0) * onePlusEps >= it->right_bound.distance(n0)){
                if(it != wavesOnEdges[e].begin() && getDiffrenceAtPoint(*it, *std::prev(it), it->left_bound) > length){
                    auto tempPrevIter = std::prev(it);
                    assert(it->left_bound.distance(n0) * onePlusEps >= it->right_bound.distance(n0));
                    it = wavesOnEdges[e].erase(it);
                }else if(std::next(it) != wavesOnEdges[e].end() && getDiffrenceAtPoint(*it, *std::next(it), it->right_bound) > length){
                    auto tempNextIter = std::next(it);
                    assert(it->left_bound.distance(n0) * onePlusEps >= it->right_bound.distance(n0));
                    it = wavesOnEdges[e].erase(it);
                }else if(length == 0){
                    it = wavesOnEdges[e].erase(it);
                }else{
                    it++;
                }
            }else{
                it++;
            }
            //it->left_bound.distance(n0) * onePlusEps >= it->right_bound.distance(n0) ? it = wavesOnEdges[e].erase(it) : ++it;
        }

        auto lastID = UINT32_MAX - 1;
        for(auto it = wavesOnEdges[e].begin(); it != wavesOnEdges[e].end(); it++){
            if(lastID == it->id && it->left_bound.distance(n0) < it->right_bound.distance(n0)){
                std::prev(it)->right_bound = it->right_bound;
                it = std::prev(wavesOnEdges[e].erase(it));
            }
            lastID = it->id;
        }

        auto temp3232 = wavesOnEdges[e];

        //if(!inserted && settings.VERBOSE_LEVEL > 2) printf("pruned edge %i from gen %i, leftb: %i, rightb: %i\n", w.item, w.generator, w.left_b, w.right_b);
        if(settings.VERBOSE_LEVEL > 2 && !isWavePresent(w)) {
            printf("\nPruned edge %i from gen %i, leftb: %i, rightb: %i\n", w.item, w.generator, w.left_b, w.right_b);
        }
        if(settings.DEBUG_LEVEL > 1) assert(debug_testWholeEdgeCovered(w.item));
        return isWavePresent(w);
    }

    bool debug_testWholeEdgeCovered(uint32_t e_id){
        if(wavesOnEdges[map_edge(e_id)].size() == 0) return true;
        auto e = map_edge(e_id);
        Point n0 = g->coords[g->firstNode(rev_map_edge(e))];
        Point n1 = g->coords[g->secondNode(rev_map_edge(e))];

        auto temp = wavesOnEdges[e];

        Coordinates left_b(n0);
        uint32_t last_id = UINT32_MAX - 1;
        uint32_t debugIter = 0;
        for(auto & wave : wavesOnEdges[e]){
            last_id = wave.id;
            if(left_b.distance(wave.left_bound) < 0.01){ //relaxed check
                left_b = wave.right_bound;
            }else{
                assert(false); return false;
            }
            debugIter++;
        }
        if(left_b.distance(n1) > 0.01){ //relaxed check
            assert(false); return false;
        }
        return true;
    }

    bool debug_testProperReset(){
        for(uint32_t e = 0; e < wavesOnEdges.size(); e++){
            if(!wavesOnEdges[e].empty()){
                printf("EdgePruning not properly reset\n");
                printf("Edge %i (%i, %i) not reset\n", rev_map_edge(e), g->firstNode(rev_map_edge(e)), g->secondNode(rev_map_edge(e)));
                assert(false); return false;
            }
        }
        return true;
    }
};
#endif

#endif