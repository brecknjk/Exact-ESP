#include <limits>
#include <sstream>

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/settings.hpp"

#include "polyanya.hpp"
#include "test.hpp"

void help_router() {
	std::cout << "USAGE: <inFile> <mode> <optinoalFlags>\n";

	std::cout << "<mode>: \t availible modes: benchmark, server, local, test, batch, worstCasePolyanya <secound input file>\n";
	std::cout << "<inFile>: \t path to graphFile\n";

	std::cout << "<optinoalFlags>:\n";
	std::cout << "\t\t-h or --help: prints help\n";
	std::cout << "\t\t-s <path_to_settings>: change default settings path\n";
	std::cout << "\t\t-v <VerboseLevel>: select verbose level, default: 1\n";
	std::cout << "\t\t-d <DebugLevel>: select debug level, default: 0\n";
	std::cout << "\t\t-tr <rounds_per_test>: number of rounds for test mode\n";
	std::cout << "\t\tfor more options modify settings.json\n";
}

void inline print_graphs(PolyAnya& polyanya) {
	// prints the graph in a format that can be used for visualization in OpenGL
	if (settings.pringGL_graph) {
		Section section = Section("write polygon to gl file: " + settings.pringGL_graph_path);
		std::vector<Coordinates> output_coords;
		std::vector<std::pair<Edge_t, uint8_t>> output_edges;
		for (auto p : polyanya.coords) {
			output_coords.emplace_back(settings.use_merkartor ? rev_bounded_mercator(Coordinates(p)) : Coordinates(p));
		}
		for (uint32_t i = 0; i < polyanya.edges.size(); i++) {
			auto& e = polyanya.edges[i];
			if(e.isConstrained()) {
				output_edges.emplace_back(std::make_pair(e, 1));
			} else if(polyanya.isDeadend(i*2)){
				if(settings.print_deadend_edges) output_edges.emplace_back(std::make_pair(e, 2));
			} else {
				if(settings.print_freeSpace_edges) output_edges.emplace_back(std::make_pair(e, 3));
			}
		}
		write_glFile(settings.pringGL_graph_path, output_coords, output_edges);
		section.end_section();
	}

	// prints the graph in a format that can be used by polyanya (and also this program)
	if (settings.pringMesh_graph) {
		Section section = Section("write mesh to file: " + settings.pringMesh_graph_path);
		polyanya.write_polyanyaStyle(settings.pringMesh_graph_path);
		section.end_section();
	}

	// prints the graph in the .off format that can be used by CGAL for shortest paths on polyhedral surfaces
	// here the obstacles are represented as holes in the mesh, so the shortest path on the mesh corresponds to the shortest path in the polygon with holes
	if(settings.pringOFF_graph) {
		Section section = Section("write graph to OFF file: " + settings.pringOFF_graph_path);
		polyanya.write_CGAL_off_style(settings.pringOFF_graph_path);
		section.end_section();
	}
}

void printSearchGraphAndRoute(PolyAnya& polyanya) {
	if (settings.pringGL_route) {
		auto route = polyanya.getRoute();
		auto section = Section("write route to gl file");
		write_route_to_glFile(settings.pringGL_route_path, route);
		section.end_section();
	}

	if (settings.pringGL_searchspace) {
		auto visited = polyanya.getSearchSpace();
		auto section = Section("write searchSpace to gl file");
		write_Searchspace_to_glFile(settings.pringGL_searchspace_path, visited);
		section.end_section();
	}

	if (settings.pringGL_searchspace && settings.pringGL_route && polyanya.stats.distance >= 0) {
		auto route = polyanya.getRoute();
		auto visited = polyanya.getSearchSpace();
		auto section = Section("write searchSpace and route to gl file");
		write_SearchspaceAndRoute_to_glFile(settings.pringGL_searchspaceAndRoute_path, visited, route);
		section.end_section();
	}
}

// For polynjaworstCase
// auto n = 1000;
// auto start = Point((n+5)/2, 1);
// auto goal = Point(1, 1);
std::pair<Coordinates, Coordinates> getQueryForWorstCasePolyanya(uint32_t n = 1000) {
	return { Coordinates(n / 2 + 2.5, 1), Coordinates(1, 1) };
}

void runPolyanyaBatchMode(PolyAnya* pa, std::vector<std::pair<Coordinates, Coordinates>> & points){
    auto section = Section("Running Batch Mode");
	pa->accumulated_stats.reset();
    if(settings.VERBOSE_LEVEL > 0) settings.polyanya.printPolyanyaSettings();
    std::ofstream os(settings.batch.outputFile);
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    
    if(settings.VERBOSE_LEVEL > 0) pa->printHeadder();
    pa->printHeadder(os);
    for(uint32_t i = 0; i < points.size(); i++){
        auto s = points[i].first;
        auto t = points[i].second;
        pa->routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
		if(settings.VERBOSE_LEVEL > 0) std::cout << pa->stats.getStats();
		os << pa->stats.getStats();
    }
    if(settings.VERBOSE_LEVEL > 0) pa->printTotalStats();
    pa->printTotalStats(os);
    section.end_section();
    
    os.close();
    if (settings.VERBOSE_LEVEL > 0) printf("written batch results to %s\n", settings.batch.outputFile.c_str());
    section.end_section();
}

void runBoxPlot(PolyAnya* pa, std::vector<std::pair<Coordinates, Coordinates>> & points){
	//fastes baseline
	std::cout << "Running Baseline" << std::endl;
	pa->accumulated_stats.reset();
	settings.polyanya.fastesSettings();
	auto fastes_stats = runPolyanya(pa, points, "Fastes");
	printSearchGraphAndRoute(*pa);

	// default polyanya
	std::cout << "Running defaultPolyanya" << std::endl;
	settings.polyanya.defaultPolyanyaSettings();
	auto defaultPA_stats = runPolyanya(pa, points, "PA");
	printSearchGraphAndRoute(*pa);

	// no deadend pruning
	std::cout << "Running noDeadend" << std::endl;
	settings.polyanya.fastesSettings();
	settings.polyanya.use_deadend_pruning = false;
	auto noDeadend_stats = runPolyanya(pa, points, "FastesNoDeadend");
	printSearchGraphAndRoute(*pa);

	std::ofstream outStream (settings.path_to_parentFolder + "/results/boxplot_data.csv", std::ios::trunc);
	outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
	outStream << "type\t" << PolyanyaStats::getStatsHeadder();
	for(const auto& stat : fastes_stats) outStream << "Fastes\t" << stat.getStats();
	for(const auto& stat : defaultPA_stats) outStream << "PA\t" << stat.getStats();
	for(const auto& stat : noDeadend_stats) outStream << "FastesNoDeadend\t" << stat.getStats();
	outStream.close();
}

int main(int argc, char** argv) {
	Section section;

	std::string inputFile, path_to_settings = "settings.json";
	std::string mode;
	try {
		inputFile = std::string(argv[1]);
		mode = std::string(argv[2]);

		for (int i = 3; i < argc; ++i) {
			std::string token = argv[i];
			if (token == "-h" || token == "--help") {
				help_router();
				return 0;
			} else if (token == "-s" && i + 1 < argc) {
				path_to_settings = std::string(argv[++i]);
			} else if (token == "-v" && i + 1 < argc) {
				settings.VERBOSE_LEVEL = std::stoi(argv[++i]);
			} else if (token == "-d" && i + 1 < argc) {
				settings.DEBUG_LEVEL = std::stoi(argv[++i]);
			} else if (token == "-m") {
				settings.use_merkartor = true;
			} else if (token == "-fp" && i + 1 < argc) {
				settings.path_to_parentFolder = std::string(argv[++i]);
			} else if (token == "-tr" && i + 1 < argc) {
				settings.tests.rounds_per_test = std::stoi(argv[++i]);
			} else {
				throw std::runtime_error("Unknown argument: " + token);
			}
		}
		settings.readSettings(path_to_settings);
	} catch (std::exception const& e) {
		std::cerr << "Trouble parsing arguments: " << e.what() << std::endl;
		help_router();
		return -1;
	}

	section = Section("init polyanya, path is: " + inputFile);
	PolyAnya polyanya(inputFile); //PolyAnyaLandmarks
	section.end_section();
	if (settings.pringMesh_graph || settings.pringGL_graph) print_graphs(polyanya);

	//get points
	std::vector<std::pair<Coordinates, Coordinates>> points;
	if (mode == "worstCasePolyanya") {
		std::cout << "Running worstCasePolyanya Mode" << std::endl;
		polyanya.printHeadder();
		uint32_t n = (polyanya.nodeCount() - 22) / 8;
		assert(n * 8 + 22 == polyanya.nodeCount());
		points.emplace_back(getQueryForWorstCasePolyanya(n));
	} else if(settings.batch.read_batchFile) {
		points = read_batchFile(settings.batch.inputFile);
		if (settings.VERBOSE_LEVEL > 0) printf("read %s\nrun %li single querys\n", settings.batch.inputFile.c_str(), points.size());
	} else if(settings.tests.longQuerys){
		auto candidatePoints = givePointsInDomain(&polyanya, settings.tests.rounds_per_test * settings.tests.multiplier_for_longQuerys);
		settings.polyanya.fastesSettings();
		auto stats_allRoutes = runPolyanya(&polyanya, candidatePoints, "Finding Benchmarks", std::cout, 2);
		std::sort(stats_allRoutes.begin(), stats_allRoutes.end(), [] (const PolyanyaStats& a, const PolyanyaStats& b) { return a.wavelets_processed < b.wavelets_processed; });
		for(uint32_t i = 0; i < settings.tests.rounds_per_test; i++) {
			points.emplace_back(stats_allRoutes[i].start, stats_allRoutes[i].goal);
		}
	} else {
		points = givePointsInDomain(&polyanya, settings.tests.rounds_per_test);
	}
	if (settings.printScenarioFile && !settings.batch.read_batchFile) createScenarioFile(settings.printScenarioFile_path, points);
	
	if (mode == "debug") {
		std::cout << "Running debug Mode" << std::endl;
		settings.polyanya.no_heuristik = false;
		settings.polyanya.use_simple_distance_heuristik = true;
		settings.polyanya.use_landmarks = false;
		settings.polyanya.use_deadend_pruning = true;
		settings.polyanya.use_edgePruning = false;
		settings.polyanya.use_advanced_pruning = false;

		settings.polyanya.use_simple_pruning_heuristik = true;
		polyanya.routePolyanya(Point(1505.193000,763.558433), Point(352.260000,385.054433));
		printSearchGraphAndRoute(polyanya);

		settings.polyanya.use_simple_pruning_heuristik = false;
		polyanya.routePolyanya(Point(1505.193000,763.558433), Point(352.260000,385.054433));
		printSearchGraphAndRoute(polyanya);

	} else if (mode == "boxplot") {
		std::cout << "Running boxplot Mode" << std::endl;
		runBoxPlot(&polyanya, points);

	} else if (mode == "server") {
		std::cout << "Running server Mode" << std::endl;
		std::cout << "TODO" << std::endl;
		//startServer(&graph, &polyanya);

	} else if (mode == "batch") {
		std::cout << "Running Batch Mode" << std::endl;
		runPolyanyaBatchMode(&polyanya, points);

	} else if (mode == "local") {
		std::cout << "Local Mode" << std::endl;
		Coordinates s, t;
		while (true) {
			try {
				std::cout << "start point x/lat?" << std::endl;
				std::cin >> s.lat;
				if (s.lat < settings.boundDomain.min_lat || s.lat > settings.boundDomain.max_lat) throw std::invalid_argument("start point x/lat out of bounds");
				std::cout << "start point y/lon?" << std::endl;
				std::cin >> s.lon;
				if (s.lon < settings.boundDomain.min_lon || s.lon > settings.boundDomain.max_lon) throw std::invalid_argument("start point y/lon out of bounds");
				std::cout << "goal point x/lat?" << std::endl;
				std::cin >> t.lat;
				if (t.lat < settings.boundDomain.min_lat || t.lat > settings.boundDomain.max_lat) throw std::invalid_argument("goal point x/lat out of bounds");
				std::cout << "goal point y/lon?" << std::endl;
				std::cin >> t.lon;
				if (t.lon < settings.boundDomain.min_lon || t.lon > settings.boundDomain.max_lon) throw std::invalid_argument("goal point y/lon out of bounds");

				section = Section("run polyanya");
				polyanya.routePolyanya(Point(s.lat, s.lon), Point(t.lat, t.lon));
				section.end_section();

				printSearchGraphAndRoute(polyanya);
			} catch (std::invalid_argument& e) {
				std::cerr << e.what() << std::endl;
			} catch (...) {
				std::cerr << "some error occured in Local Mode" << std::endl;
			}
		}

	} else if (mode == "test") {
		std::cout << "Test Mode" << std::endl;
		try {
			uint32_t rounds = settings.tests.rounds_per_test;

			assert(settings.polyanya.use_deadend_pruning);
			bool use_edgePruning = settings.polyanya.use_edgePruning;
			bool use_landmarks = settings.polyanya.use_landmarks;
			settings.polyanya.use_edgePruning = false;
			settings.polyanya.use_landmarks = false;

			if (!test_003_pruning_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_003_pruning_correctnes failed!"));
			if (!test_003b_advanced_pruning_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_003b_advanced_pruning_correctnes failed!"));
			if (!test_007_deadend_pruning_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_007_deadend_pruning_correctnes failed!"));
			if (!test_001_simpleDistHeuristic_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_001_simpleDistHeuristic_correctnes failed!"));
			if (!test_002_noHeuristic_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_002_noHeuristic_correctnes failed!"));
			if (!test_004_addingRemovingPoints(&polyanya, rounds)) throw(std::runtime_error("test_004_addingRemovingPoints failed!"));
			if (!test_005_addingRemovingPoints_repeated(&polyanya, rounds)) throw(std::runtime_error("test_005_addingRemovingPoints_repeated failed!"));
			if (use_edgePruning) {
				if (!test_006_perfect_pruning_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_006_perfect_pruning_correctnes failed!"));
			}
			if (use_landmarks) {
				if (!test_004_landmarks_correctnes(&polyanya, rounds)) throw(std::runtime_error("test_004_landmarks_correctnes failed!"));
			}

		} catch (std::runtime_error& e) {
			std::cerr << e.what() << std::endl;
		} catch (...) {
			std::cerr << "some error occured in Test Mode" << std::endl;
		}

	} else if (mode == "benchmark") {
		std::cout << "Benchmark Mode" << std::endl;

		assert(settings.polyanya.use_deadend_pruning);
		bool use_edgePruning = settings.polyanya.use_edgePruning;
		bool use_landmarks = settings.polyanya.use_landmarks;
		settings.polyanya.use_edgePruning = false;
		settings.polyanya.use_landmarks = false;

		try {
			//fastes baseline
			std::cout << "Running Baseline" << std::endl;
			settings.polyanya.fastesSettings();
			auto dist_baseline = runPolyanya(&polyanya, points, "Baseline: Fastes");

			// justSimplePruning
			std::cout << "Running justSimplePruning" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.use_advanced_pruning = false;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "JustSimplePruning"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_justSimplePruning"));
			}

			// dist_advancedPruning_noSimple
			std::cout << "Running dist_advancedPruning_noSimple" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.use_simple_pruning_heuristik = false;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "JustAdvancedPruning"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_advancedPruning_noSimple"));
			}

			// fastes with skipahead
			std::cout << "Running fastes with skipahead" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.skip_ahead = true;
			runPolyanya(&polyanya, points, "FastSkipAhead");
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "FastSkipAhead"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at fastes with skipahead"));
			}

			// no deadend pruning
			std::cout << "Running noDeadend" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.use_deadend_pruning = false;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "FastesNoDeadend"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_noDeadend"));
			}

			// fastes no break condition
			std::cout << "Running fastes no break condition" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.breakCondition = false;
			settings.polyanya.no_heuristik = true;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "FastNoBreak"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at no break condition"));
			}

			// fastes no distance heuristik
			std::cout << "Running fastes no distance heuristik" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.no_heuristik = true;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "FastesNoA*"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_fastesNoDistanceHeuristic"));
			}

			// fastes no break condition
			std::cout << "Running fastes no break no deadend condition" << std::endl;
			settings.polyanya.fastesSettings();
			settings.polyanya.breakCondition = false;
			settings.polyanya.use_deadend_pruning = false;
			settings.polyanya.no_heuristik = true;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "dist_noBreakNoDeadendCondition"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_noBreakNoDeadendCondition"));
			}

			// default polyanya
			std::cout << "Running defaultPolyanya" << std::endl;
			settings.polyanya.defaultPolyanyaSettings();
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "PA"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_defaultPolyanya"));
			}

			// no distance heuristik
			std::cout << "Running no distance heuristik" << std::endl;
			settings.polyanya.defaultPolyanyaSettings();
			settings.polyanya.no_heuristik = true;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "PA_NoA*"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_noDistanceHeuristic"));
			}

			// no distance heuristik
			std::cout << "Running PA no break heuristik" << std::endl;
			settings.polyanya.defaultPolyanyaSettings();
			settings.polyanya.no_heuristik = true;
			settings.polyanya.breakCondition = false;
			if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "PA_NoA*NoBreak"))) {
				if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_PAnoBreakHeuristic"));
			}

			if (use_landmarks) {
				//landmarks
				std::cout << "Running landmarks" << std::endl;
				settings.polyanya.fastesSettings();
				settings.polyanya.use_landmarks = true;
				if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "Landmarks"))) {
					if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_landmarks"));
				}
			}

			if (use_edgePruning) {
				//edgePruning
				std::cout << "Running edgePruning" << std::endl;
				settings.polyanya.edgePruningSettings();
				if(!compareDistVectors(dist_baseline, runPolyanya(&polyanya, points, "EdgePruning"))) {
					if(settings.DEBUG_LEVEL > 0) throw(std::runtime_error("distances error at dist_edgePruning"));
				}
			}

		} catch (std::runtime_error& e) {
			std::cerr << e.what() << std::endl;
		} catch (...) {
			std::cerr << "some error occured in Benchmark Mode" << std::endl;
		}
	} else {
		std::cerr << "Unknown Mode" << std::endl;
		help_router();
		return -1;
	}
	std::cout << "Final Memory usage: " << get_mem_usage() << " KB" << std::endl;
	return 0;
};