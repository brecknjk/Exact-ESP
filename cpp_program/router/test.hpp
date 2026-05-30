#ifndef testsPolayanja_h
#define testsPolayanja_h

#include <unordered_map>
#include <stdexcept>

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

#include "edgeLocator.hpp"
#include "deadendFinder.hpp"
#include "polyanya.hpp"

Coordinates randomPoint(){
    double lat = (rand() % (uint64_t(settings.boundDomain.max_lat - settings.boundDomain.min_lat) * 1000)) / 1000.0 + settings.boundDomain.min_lat;
    assert(lat >= settings.boundDomain.min_lat && lat <= settings.boundDomain.max_lat);
    double lon = (rand() % (uint64_t(settings.boundDomain.max_lon - settings.boundDomain.min_lon) * 1000)) / 1000.0 + settings.boundDomain.min_lon;
    assert(lon >= settings.boundDomain.min_lon && lon <= settings.boundDomain.max_lon);
    return Coordinates(lat, lon);
}

std::vector<std::pair<Coordinates, Coordinates>> givePointsInDomain(PolyAnya* pa, uint32_t count = 1000, long seed = -1){
    std::vector<std::pair<Coordinates, Coordinates>> points;
    uint32_t debug_iter = 0;
    (seed != -1) ? srand(seed) : srand(time(NULL) * getpid());
    while(points.size() < count){
        if(debug_iter++ > 100 * count) throw(std::runtime_error("to many points are out of bound"));
        auto s = randomPoint();
        auto t = randomPoint();
        if(pa->inDomain(Point(s.lat, s.lon)) && pa->inDomain(Point(t.lat, t.lon))) {
            points.push_back(std::make_pair(s, t));
        }
    }
    return points;
}

bool compareDistVectors(const std::vector<PolyanyaStats>& a, const std::vector<PolyanyaStats>& b, double error = 0.0000001){
    assert(a.size() == b.size());
    for(uint32_t i = 0; i < a.size(); i++){
        if(fabs(a[i].distance - b[i].distance) > error) {
            std::cout << "difference at index " << i << ": " << a[i].distance << " vs. " << b[i].distance << "\n";
            std::cout << "stats a: " << a[i].getStats() << "\n";
            std::cout << "stats b: " << b[i].getStats() << "\n";
            return false;
        }
    }
    return true;
}

std::vector<PolyanyaStats> runPolyanya(PolyAnya* pa, 
                                std::vector<std::pair<Coordinates, Coordinates>> & points, 
                                std::string sectionName = "",
                                std::ostream &out = std::cout, 
                                int verbose_level = 1,
                                std::vector<std::list<Coordinates>>* routes = nullptr
                            ){
    std::vector<PolyanyaStats> statsList;
    pa->resetSearch();
    pa->accumulated_stats.reset();
    if(settings.VERBOSE_LEVEL > 0) settings.polyanya.printPolyanyaSettings();
    if(verbose_level < settings.VERBOSE_LEVEL) pa->printHeadder(out);
    auto section = Section(sectionName, verbose_level);
    for(uint32_t i = 0; i < points.size(); i++){
        auto s = points[i].first;
        auto t = points[i].second;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        if(routes) routes->push_back(pa->getRoute_withoutMercator());
        if(verbose_level < settings.VERBOSE_LEVEL) out << pa->stats.getStats();
        statsList.push_back(pa->stats);
    }
    if(verbose_level <= settings.VERBOSE_LEVEL) pa->printTotalHeadder(out);
    if(verbose_level <= settings.VERBOSE_LEVEL) pa->printTotalStats(out);
    section.end_section();
    assert(statsList.size() == points.size());
    return statsList;
}

bool test_001_simpleDistHeuristic_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    srand(time(NULL) * getpid());
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.no_heuristik = false;
    settings.polyanya.use_simple_pruning_heuristik = true;
    settings.polyanya.use_deadend_pruning = true;
    settings.polyanya.use_landmarks = false;
    settings.polyanya.use_edgePruning = false;
    
    auto section = Section("run test_001_simpleDistHeuristic_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.use_simple_distance_heuristik = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));

        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;

        settings.polyanya.use_simple_distance_heuristik = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));

        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;

        if(dist1 != dist2 ){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_001_simpleDistHeuristic_correctnes\n");
        };
    }
    section.end_section();
    printf("time no_simple_h: %lf ms vs. time simple_h: %lf ms\n", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops no_simple_h: %lf vs. pq-Pops simple_h: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    return true;
}

bool test_002_noHeuristic_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    srand(time(NULL) * getpid());
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.use_simple_distance_heuristik = true;
    settings.polyanya.use_simple_pruning_heuristik = true;
    settings.polyanya.use_deadend_pruning = true;
    settings.polyanya.use_landmarks = false;
    settings.polyanya.use_edgePruning = false;

    auto section = Section("run test_002_noHeuristic_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.no_heuristik = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));

        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;

        settings.polyanya.no_heuristik = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));

        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;

        if(dist1 != dist2){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_002_noHeuristic_correctnes");
        };
    }
    section.end_section();
    printf("time heuristik: %lf ms vs. time no_heuristik: %lf ms", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops heuristik: %lf vs. pq-Pops no_heuristik: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    settings.polyanya.no_heuristik = false;
    return true;
}

bool test_003_pruning_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    srand(time(NULL) * getpid());
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.no_heuristik = false;
    settings.polyanya.use_simple_distance_heuristik = true;
    settings.polyanya.use_landmarks = false;
    settings.polyanya.use_deadend_pruning = true;
    settings.polyanya.use_edgePruning = false;
    settings.polyanya.use_advanced_pruning = false;

    settings.polyanya.printPolyanyaSettings();
    auto section = Section("run test_003_pruning_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.use_simple_pruning_heuristik = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;

        settings.polyanya.use_simple_pruning_heuristik = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;

        if(dist1 != dist2){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_003_pruning_correctnes");
        };
    }
    section.end_section();
    printf("time noPrune: %lf ms vs. time prune: %lf ms", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops noPrune: %lf vs. pq-Pops prune: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    return true;
}

bool test_003b_advanced_pruning_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    srand(time(NULL) * getpid());
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.no_heuristik = false;
    settings.polyanya.use_simple_distance_heuristik = true;
    settings.polyanya.use_landmarks = false;
    settings.polyanya.use_deadend_pruning = true;
    settings.polyanya.use_edgePruning = false;
    settings.polyanya.use_simple_pruning_heuristik = false;

    auto section = Section("run test_003b_advanced_pruning_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.use_advanced_pruning = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;

        settings.polyanya.use_advanced_pruning = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;

        if(dist1 != dist2){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_003_pruning_correctnes");
        };
    }
    section.end_section();
    printf("time noPrune: %lf ms vs. time prune: %lf ms", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops noPrune: %lf vs. pq-Pops prune: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    return true;
}

bool test_004_landmarks_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    srand(time(NULL) * getpid());
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.no_heuristik = false;
    settings.polyanya.use_simple_distance_heuristik = true;
    settings.polyanya.use_simple_pruning_heuristik = true;
    settings.polyanya.use_deadend_pruning = false;

    auto section = Section("run test_004_landmarks_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.use_landmarks = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;
        pa->resetSearch(); // needed for switching use_landmarks

        settings.polyanya.use_landmarks = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;
        pa->resetSearch(); // needed for switching use_landmarks

        if(dist1 != dist2){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_004_landmarks_correctnes");
        };
    }
    section.end_section();
    printf("time noLandmark: %lf ms vs. time Landmark: %lf ms\n", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops noLandmark: %lf vs. pq-Pops Landmark: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    return true;
}

bool test_004_addingRemovingPoints(PolyAnya* pa, uint32_t rounds = 100){
    srand(time(NULL) * getpid());
    try{
        for(uint32_t i = 0; i < rounds; i++){
            auto c = randomPoint();
            auto p = Point(c.lat, c.lon);
            if(settings.use_merkartor) p = bounded_mercator(p);
            auto nearest_edge = pa->getNearestEdge(p);
            if(nearest_edge == UINT32_MAX) continue;
            pa->insertPoint(p, nearest_edge, true);
            pa->deleteLastPoint();
        }
    }catch(...){
        printf("error at test_004_addingRemovingPoints");
        return false;
    }
    return true;
}

bool test_005_addingRemovingPoints_repeated(PolyAnya* pa, uint32_t rounds = 100){
    srand(time(NULL) * getpid());
    std::unordered_map<uint32_t, Point> points_to_insert;
    for(uint32_t i = 0; i < rounds; i++){
        auto c = randomPoint();
        auto p = Point(c.lat, c.lon);
        if(settings.use_merkartor) p = bounded_mercator(p);
        auto nearest_edge = pa->getNearestEdge(p);
        if(nearest_edge == UINT32_MAX) continue;
        points_to_insert[nearest_edge] = p;
    }
    try{
        for(auto [k,v] : points_to_insert){
            pa->insertPoint(v, k);
        }
        for(uint32_t i = 0; i < points_to_insert.size(); i++){
            pa->deleteLastPoint();
        }
    }catch(...){
        printf("error at test_005_addingRemovingPoints_repeated");
        return false;
    }
    return true;
}

bool test_006_perfect_pruning_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    auto seed = time(NULL) * getpid();
    srand(seed);
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.no_heuristik = false;
    settings.polyanya.use_simple_distance_heuristik = true;
    settings.polyanya.use_landmarks = false;
    settings.polyanya.breakCondition = true;

    auto section = Section("run test_006_perfect_pruning_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.use_simple_pruning_heuristik = false;
        settings.polyanya.use_edgePruning = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;

        settings.polyanya.use_simple_pruning_heuristik = false;
        settings.polyanya.use_edgePruning = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;

        if(abs(dist1 - dist2) > NT(0.0000001)){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_006_perfect_pruning_correctnes");
        };
    }
    section.end_section();
    printf("time noPrune: %lf ms vs. time prune: %lf ms\n", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops noPrune: %lf vs. pq-Pops prune: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    return true;
}

bool test_007_deadend_pruning_correctnes(PolyAnya* pa, uint32_t rounds = 1000){
    double total_time_a = 0, total_time_b = 0;
    uint64_t total_pq_pops_a = 0, total_pq_pops_b = 0;
    srand(time(NULL) * getpid());
    auto p = givePointsInDomain(pa, rounds);

    settings.polyanya.no_heuristik = false;
    settings.polyanya.use_simple_distance_heuristik = true;
    settings.polyanya.use_landmarks = false;
    settings.polyanya.breakCondition = true;
    settings.polyanya.use_edgePruning = false;
    settings.polyanya.use_simple_pruning_heuristik = true;

    auto section = Section("run test_007_deadend_pruning_correctnes");
    if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStatsHeadder();
    for(uint32_t i = 0; i < rounds; i++){
        Coordinates s = p[i].first;
        Coordinates t = p[i].second;

        settings.polyanya.use_deadend_pruning = false;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist1 = pa->stats.distance;
        total_time_a += pa->stats.time_ms;
        total_pq_pops_a += pa->stats.pq_pops;
        pa->resetSearch(); // needed for deadend pruning

        settings.polyanya.use_deadend_pruning = true;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
        auto dist2 = pa->stats.distance;
        total_time_b += pa->stats.time_ms;
        total_pq_pops_b += pa->stats.pq_pops;
        pa->resetSearch(); // needed for deadend pruning

        if(abs(dist1 - dist2) > NT(0.00001)){
            printf("error at iter %d: from %s to %s dist: %lf vs. %lf\n", 
                i, 
                to_string(s).c_str(),
                to_string(t).c_str(),
                to_double(dist1),
                to_double(dist2));
            throw std::runtime_error("Distances do not match in test_007_deadend_pruning_correctnes");
        };
    }
    section.end_section();
    printf("time noDeadend: %lf ms vs. time deadend: %lf ms\n", total_time_a / (double)rounds, total_time_b / (double)rounds);
    printf("pq-Pops noDeadend: %lf vs. pq-Pops deadend: %lf\n", total_pq_pops_a / (double)rounds, total_pq_pops_b / (double)rounds);
    return true;
}

#endif