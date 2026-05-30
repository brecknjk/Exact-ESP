#ifndef landmarks_h
#define landmarks_h

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

#include "graph.hpp"


//to_lowerFloat()

// LandMarks?
//     each node: successsor + dist

//     heuristik: dist(w,t) + cost[w] >= dist(w,l) + cost[w] - dist(t,l)
//         let t be in cell abc: dist(t,l) <= dist(a,l) + dist(t,a)
//         dist(w,l) + cost[w] >= cost[gen[w]] + dist(gen[w],suc_l[gen[w]])over(w) + dist(suc_l[gen[w]],l)
//     alt: dist(w,l) + cost[w] >= max(a,b): cost[gen[w]] + dist(a,l) + dist(a,gen[w]) - 2*dist(a,b)
//     alt: dist(w,l) + cost[w] >= min(a,b): cost[gen[w]] + dist(a,l) + dist(a,gen[w]) - dist(a,b)

class Landmarks{
protected:
    Graph* g;
    uint16_t LANDMARK_COUNT = 8;
    std::vector<std::vector<double>> landmark_distances;
    
    // for search
    std::vector<double> upperbounds_goals_distances;
    uint32_t goal_id;
    std::vector<uint16_t> selectedLandmarks;

    // debug value
    std::vector<uint32_t> landmarks;
public:
    inline uint32_t getNodeCount(){return landmark_distances[0].size();}
    //inline void routeLandmarks(Point start_point, Point goal_point);
    inline void printLandmarks(){
        if(settings.VERBOSE_LEVEL < 0) return;
        printf("%i landmarks: \n", LANDMARK_COUNT);
        for(uint32_t l : landmarks) printf("\t(%lf, %lf) \n", to_double(g->coords[l].x()), to_double(g->coords[l].y()));
    }
    Landmarks(){}
    Landmarks(Graph* g, uint16_t LANDMARK_COUNT = 8) : g(g), LANDMARK_COUNT(LANDMARK_COUNT){
        landmark_distances.reserve(LANDMARK_COUNT);
        landmarks.reserve(LANDMARK_COUNT);
        upperbounds_goals_distances.reserve(LANDMARK_COUNT);
    }
    //selection: 1 min y-value
    uint32_t getInitalLandmark(){
        uint32_t landmark = 0;
        // find first landmark
        for(uint32_t i = 1; i < g->coords.size(); i++){
            if(g->coords[i].y() < g->coords[landmark].y()) landmark = i;
        }
        return landmark;
    }
    //choose l: max: dist(w,l) - dist(t,l)
    uint32_t addLandmark(uint32_t landmark, std::vector<NT_Distance>& dists){
        // set dist of landmark
        std::vector<double> dists_double;
        for(auto & dist : dists){
            dists_double.push_back(to_double(dist));
        }
        landmark_distances.emplace_back(dists_double);
        landmark_distances.back().reserve(landmark_distances.back().size() + 2);
        landmarks.emplace_back(landmark);
        // find next Landmark
        uint32_t nextLandmark = UINT32_MAX;
        double nextLandmark_dist = -MAX_DOUBLE;
        for(uint32_t n = 0; n < dists_double.size(); n++){
            double current_min = MAX_DOUBLE;
            for(uint32_t l = 0; l < landmark_distances.size(); l++){
                current_min = std::min(current_min, landmark_distances[l][n]);
            }
            if(current_min > nextLandmark_dist){
                nextLandmark_dist = current_min;
                nextLandmark = n;
            }
        }
        return nextLandmark;
    }
    // choose l: max: dist(s,l) - dist(t,l)
    void select_landmarks(uint32_t start, uint32_t goal, uint16_t landmarks_to_use = 1){
        assert(selectedLandmarks.size() == 0);
        assert(landmarks_to_use <= LANDMARK_COUNT);
        assert(goal == goal_id);

        std::vector<bool> is_selected(LANDMARK_COUNT, false);
        for(uint16_t c = 0; c < landmarks_to_use; c++){
            double best_v = -MAX_DOUBLE;
            uint16_t best_l = UINT8_MAX;
            for(uint16_t l = 0; l < LANDMARK_COUNT; l++){
                if(is_selected[l]) continue;
                if(best_v < landmark_distances[l][start] - landmark_distances[l][goal]) {
                    best_v = landmark_distances[l][start] - landmark_distances[l][goal];
                    best_l = l;
                }
            }
            if(settings.VERBOSE_LEVEL > 1) printf("selected landmark %i\n", best_l);
            assert(best_l != UINT8_MAX);
            is_selected[best_l] = true;
            selectedLandmarks.emplace_back(best_l);
        }
        assert(selectedLandmarks.size() == landmarks_to_use);
    }

//TODO select best heuristic calc
//     heuristik: dist(w,t) + cost[w] >= dist(w,l) + cost[w] - dist(t,l)
//         let t be in cell abc: dist(t,l) <= dist(a,l) + dist(t,a)
//         dist(w,l) + cost[w] >= cost[gen[w]] + dist(gen[w],suc_l[gen[w]])over(w) + dist(suc_l[gen[w]],l)
//     alt: dist(w,l) + cost[w] >= max(a,b): cost[gen[w]] + dist(a,l) + dist(a,gen[w]) - 2*dist(a,b)
//     alt: dist(w,l) + cost[w] >= min(a,b): cost[gen[w]] + dist(a,l) + dist(a,gen[w]) - dist(a,b)
    double getHeuristik(const Wavelet & w) const {
        uint32_t n0 = g->firstNode(w.item), n1 = g->secondNode(w.item);

        auto temp2 = std::min(g->getDist_exact_index(w.generator, n0), g->getDist_exact_index(w.generator, n1));
        auto temp3 = CGAL::sqrt(to_double(CGAL::squared_distance(g->coords[w.generator], K::Segment_2 (g->coords[n0], g->coords[n1]))));
        assert(temp3 <= temp2 + EPS_FLOAT);

        double distToEdge = CGAL::sqrt(to_double(CGAL::squared_distance(g->coords[w.generator], K::Segment_2 (g->coords[n0], g->coords[n1]))));  
        double max_heuristik = -MAX_DOUBLE;
        for(auto l : selectedLandmarks){
            // std::cout << to_string(rev_bounded_mercator(g->coords[n0])) << std::endl; //DEBUG
            auto c = landmark_distances[l][n0];
            auto d = landmark_distances[l][n1];
            auto a = landmark_distances[l][n0] + g->getDist_inexact_index(n0, w.generator);
            auto b = landmark_distances[l][n1] + g->getDist_inexact_index(n1, w.generator);
            auto temp0 = (std::max(a, b) - 2 * g->getDist_inexact_index(n0, n1)) - upperbounds_goals_distances[l];
            auto temp5 = (std::min(a, b) - g->getDist_inexact_index(n0, n1)) - upperbounds_goals_distances[l];
            auto temp6 = (std::max(c, d) - g->getDist_inexact_index(n0, n1)) + distToEdge - upperbounds_goals_distances[l];
            auto temp7 = (std::min(c, d) - 0.5*g->getDist_inexact_index(n0, n1)) + distToEdge - upperbounds_goals_distances[l];
            auto temp4 = std::max(temp5, temp7); 
            auto temp8 = std::max(temp0, temp6);
            auto temp1 = landmark_distances[l][goal_id] - (std::min(a, b));
            auto current_heuristik = std::max(temp1, std::max(temp4, temp8));

            // auto e = landmark_distances[l][n0] - g->getDist_exact_index(n0, n1) + distToEdge;
            // auto f = landmark_distances[l][n1] - g->getDist_exact_index(n0, n1) + distToEdge;
            // double temp0 = abs(std::max(e, f) - upperbounds_goals_distances[l]);
            // double temp1 = abs(std::min(a, b) - landmark_distances[l][goal_id]);
            
            // if(temp4 == max_heuristik && temp4 + w.cost_generator > distToEdge) assert(temp5 <= temp7);
            // if(temp8 == max_heuristik && temp8 + w.cost_generator > distToEdge) assert(temp0 <= temp6);
            max_heuristik = std::max(max_heuristik, current_heuristik);
        }
        return to_double(w.cost_generator) + std::max(distToEdge, max_heuristik);
    }

    void reset(){
        selectedLandmarks.clear();
        upperbounds_goals_distances.clear();
        goal_id = UINT32_MAX;
    };
    void insertPoint(uint32_t p, uint32_t edgeOfFace, bool is_goal = false){
        if(is_goal && goal_id != UINT32_MAX && settings.VERBOSE_LEVEL > 0) printf("warning adding goal, point without reset!\n");
        if(is_goal) goal_id = p;
        for(uint32_t l = 0; l < LANDMARK_COUNT; l++){
            double lowerbound = MAX_DOUBLE, upperbound = 0;
            for(FaceCirculator f_iter(g, edgeOfFace); f_iter.hasNext(); ++f_iter){
                auto n = f_iter.firstNode();
                lowerbound = std::min(lowerbound, landmark_distances[l][n] - g->getDist_inexact_index(p, n));
                upperbound = std::max(upperbound, landmark_distances[l][n] + g->getDist_inexact_index(p, n));
            }
            assert(lowerbound <= upperbound);
            landmark_distances[l].push_back(lowerbound);
            if(is_goal) upperbounds_goals_distances.push_back(upperbound);
        }
    };
    void deleteLastPoint(){
        for(uint32_t l = 0; l < LANDMARK_COUNT; l++){
            landmark_distances[l].pop_back();
        }
    };
};

#endif