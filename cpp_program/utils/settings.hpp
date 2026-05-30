#ifndef settings_hpp
#define settings_hpp

#include <boost/property_tree/ptree.hpp>                                        
#include <boost/property_tree/json_parser.hpp> 
#include <fstream>

#include "types.hpp"
#include "io.hpp"
#include "mathLib.hpp"

struct Settings {
    struct BoundDomain{
        double min_lat = -89.5;
        double min_lon = -179.9;
        double max_lat = 89.5;
        double max_lon = 179.9;
    };
    struct PolyanyaSettings{
        bool print_out_exact_pathCosts = false;
        bool use_simple_distance_heuristik = true;
        bool no_heuristik = false;
        bool use_deadend_pruning = true;
        bool use_simple_pruning_heuristik = true;
        bool use_advanced_pruning = true;
        bool heuristikEdgePruningProof = false;
        bool use_edgePruning = true;
        bool merge_cells = true;
        bool flipEdgesForMoreDeadends = true;
        bool merge_cells_before_deadend_pruning = false;
        bool dont_gen_twice = true;
        bool breakCondition = true; // break if current_costs >= costs[goalID]
        bool skip_ahead = false;
        bool use_landmarks = true;
        bool ignore_collinear = true;
        uint16_t landmark_count = 8;
        uint16_t landmarks_to_use = 4;

        void defaultPolyanyaSettings(){
            use_simple_distance_heuristik = false;
            no_heuristik = false;
            use_deadend_pruning = false;
            use_simple_pruning_heuristik = false;
            use_advanced_pruning = false;
            heuristikEdgePruningProof = false;
            use_edgePruning = false;
            dont_gen_twice = true;
            breakCondition = true;
            skip_ahead = true;
            use_landmarks = false;
            ignore_collinear = false;
        };

        void fastesSettings(){
            use_simple_distance_heuristik = true;
            no_heuristik = false;
            use_deadend_pruning = true;
            use_simple_pruning_heuristik = true;
            use_advanced_pruning = true;
            heuristikEdgePruningProof = false;
            use_edgePruning = false;
            dont_gen_twice = true;
            breakCondition = true;
            skip_ahead = false;
            use_landmarks = false;
            ignore_collinear = true;
        };

        void edgePruningSettings(){
            use_simple_distance_heuristik = false;
            no_heuristik = false;
            use_deadend_pruning = true;
            use_simple_pruning_heuristik = false;
            use_advanced_pruning = false;
            heuristikEdgePruningProof = false;
            use_edgePruning = true;
            dont_gen_twice = true;
            breakCondition = true;
            skip_ahead = false;
            use_landmarks = false;
            ignore_collinear = false;
        };

        void printPolyanyaSettings(){
            std::cout << "Polyanya Settings:" << std::endl;
            if(no_heuristik){
                std::cout << "  Heuristik: none" << std::endl;
            }else{
                std::cout << "  use_simple_distance_heuristik: " << use_simple_distance_heuristik << std::endl;
            }
            std::cout << "  use_deadend_pruning: " << use_deadend_pruning << std::endl;
            std::cout << "  use_simple_pruning_heuristik: " << use_simple_pruning_heuristik << std::endl;
            std::cout << "  use_advanced_pruning: " << use_advanced_pruning << std::endl;
            std::cout << "  heuristikEdgePruningProof: " << heuristikEdgePruningProof << std::endl;
            std::cout << "  use_edgePruning: " << use_edgePruning << std::endl;
            std::cout << "  dont_gen_twice: " << dont_gen_twice << std::endl;
            std::cout << "  breakCondition: " << breakCondition << std::endl;
            std::cout << "  skip_ahead: " << skip_ahead << std::endl;
            std::cout << "  ignore_collinear: " << ignore_collinear << std::endl;
            if(merge_cells){
                std::cout << "  merge_cells: " << merge_cells << std::endl;
                std::cout << "  merge_cells_before_deadend_pruning: " << merge_cells_before_deadend_pruning << std::endl;
            }
            if(use_landmarks){
                std::cout << "  use_landmarks: " << use_landmarks << std::endl;
                std::cout << "  landmark_count: " << landmark_count << std::endl;
                std::cout << "  landmarks_to_use: " << landmarks_to_use << std::endl;
            }
        };
    };
    struct ExtractorSettings{
        bool pruneTriangulation = true;
        float pruneRatio_longestToShortestEdge2 = 20.0;
        float pruneRatio_maxVertexCount = 1.5;
    };
    struct BatchSettings{
        std::string inputFile = "";
        std::string outputFile = "";
        bool read_batchFile = false;
        bool write_batchFile = false;
    };
    struct ServerSettings{
        std::string host = "localhost";
        int port = 8080;
    };
    struct TestSettings{
        uint32_t rounds_per_test = 0;
        bool longQuerys = false;
        uint32_t multiplier_for_longQuerys = 2;
    };
    struct LocalSettings{};

    int DEBUG_LEVEL = -99;
    int VERBOSE_LEVEL = -99;

    std::string path_to_parentFolder = "";

    Coordinates seaPoint = Coordinates(INFINITY, INFINITY);
    BoundDomain boundDomain;
    bool use_merkartor = false;

    PolyanyaSettings polyanya;

    bool pringGL_searchspace = false;
    std::string pringGL_searchspace_path = "";
    bool pringGL_route = false;
    std::string pringGL_route_path = "";
    bool pringGL_graph = false;
    std::string pringGL_searchspaceAndRoute_path = "";
    std::string pringGL_graph_path = "";
    bool print_deadend_edges = true;
    bool print_freeSpace_edges = true;

    bool pringMesh_graph = false;
    std::string pringMesh_graph_path = "";
    bool printScenarioFile = false;
    std::string printScenarioFile_path = "";

    bool pringOFF_graph = false;
    std::string pringOFF_graph_path = "";

    ExtractorSettings extractor;

    BatchSettings batch;
    ServerSettings server;
    TestSettings tests;
    LocalSettings local;

    Settings() {}
    void readSettings(const std::string & settings_path) {
        namespace pt = boost::property_tree;
        pt::ptree settings;
        // Load settings from a file or initialize with default values
        try {
            pt::read_json(settings_path, settings);
        } catch (std::exception const& e) {
            std::cerr << "Trouble reading settings.json: " << e.what() << std::endl;
            throw std::runtime_error("Failed to read settings file.");
        }

        try {
            if(DEBUG_LEVEL == -99) DEBUG_LEVEL = settings.get_child("DEBUG_LEVEL").get_value<int>();
            if(VERBOSE_LEVEL == -99) VERBOSE_LEVEL = settings.get_child("VERBOSE_LEVEL").get_value<int>();

            if(path_to_parentFolder == "") path_to_parentFolder = settings.get_child("path_to_parentFolder").get_value<std::string>();

            if(seaPoint.lat == INFINITY) seaPoint.lat = settings.get_child("seaPoint").get_child("lat").get_value<double>();
            if(seaPoint.lon == INFINITY) seaPoint.lon = settings.get_child("seaPoint").get_child("lon").get_value<double>();
            boundDomain.min_lat = settings.get_child("bound_domain").get_child("min_lat").get_value<double>();
            boundDomain.min_lon = settings.get_child("bound_domain").get_child("min_lon").get_value<double>();
            boundDomain.max_lat = settings.get_child("bound_domain").get_child("max_lat").get_value<double>();
            boundDomain.max_lon = settings.get_child("bound_domain").get_child("max_lon").get_value<double>();
            if(!use_merkartor) use_merkartor = settings.get_child("use_merkartor").get_value<bool>();

            auto polyanya_settings = settings.get_child("polyanya_settings");
            polyanya.use_simple_distance_heuristik = polyanya_settings.get_child("use_simple_distance_heuristik").get_value<bool>();
            polyanya.no_heuristik = polyanya_settings.get_child("no_heuristik").get_value<bool>();
            polyanya.use_deadend_pruning = polyanya_settings.get_child("use_deadend_pruning").get_value<bool>();
            polyanya.use_simple_pruning_heuristik = polyanya_settings.get_child("use_simple_pruning_heuristik").get_value<bool>();
            polyanya.use_advanced_pruning = polyanya_settings.get_child("use_advanced_pruning").get_value<bool>();
            polyanya.use_edgePruning = polyanya_settings.get_child("use_edgePruning").get_value<bool>();
            polyanya.merge_cells = polyanya_settings.get_child("merge_cells").get_value<bool>();
            polyanya.merge_cells_before_deadend_pruning = polyanya_settings.get_child("merge_cells_before_deadend_pruning").get_value<bool>();
            polyanya.dont_gen_twice = polyanya_settings.get_child("dont_gen_twice").get_value<bool>();
            polyanya.breakCondition = polyanya_settings.get_child("breakCondition").get_value<bool>();
            polyanya.skip_ahead = polyanya_settings.get_child("skip_ahead").get_value<bool>();
            polyanya.use_landmarks = polyanya_settings.get_child("use_landmarks").get_value<bool>();
            polyanya.landmark_count = polyanya_settings.get_child("landmark_count").get_value<uint16_t>();
            polyanya.landmarks_to_use = polyanya_settings.get_child("landmarks_to_use").get_value<uint16_t>();
            polyanya.ignore_collinear = polyanya_settings.get_child("ignore_collinear").get_value<bool>();

            pringGL_searchspace = settings.get_child("pringGL_searchspace").get_value<bool>();
            pringGL_searchspace_path = path_to_parentFolder + settings.get_child("pringGL_searchspace_path").get_value<std::string>();
            pringGL_route = settings.get_child("pringGL_route").get_value<bool>();
            pringGL_route_path = path_to_parentFolder + settings.get_child("pringGL_route_path").get_value<std::string>();
            pringGL_searchspaceAndRoute_path = path_to_parentFolder + settings.get_child("pringGL_searchspaceAndRoute_path").get_value<std::string>();
            pringGL_graph = settings.get_child("pringGL_graph").get_value<bool>();
            pringGL_graph_path = path_to_parentFolder + settings.get_child("pringGL_graph_path").get_value<std::string>();
            print_deadend_edges = settings.get_child("print_deadend_edges").get_value<bool>();
            print_freeSpace_edges = settings.get_child("print_freeSpace_edges").get_value<bool>();

            pringMesh_graph = settings.get_child("pringMesh_graph").get_value<bool>();
            pringMesh_graph_path = path_to_parentFolder + settings.get_child("pringMesh_graph_path").get_value<std::string>();
            printScenarioFile = settings.get_child("printScenarioFile").get_value<bool>();
            printScenarioFile_path = path_to_parentFolder + settings.get_child("printScenarioFile_path").get_value<std::string>();
            pringOFF_graph = settings.get_child("pringOFF_graph").get_value<bool>();
            pringOFF_graph_path = path_to_parentFolder + settings.get_child("pringOFF_graph_path").get_value<std::string>();

            auto extractor_settings = settings.get_child("extractor");
            extractor.pruneTriangulation = extractor_settings.get_child("pruneTriangulation").get_value<bool>();
            extractor.pruneRatio_longestToShortestEdge2 = extractor_settings.get_child("pruneRatio_longestToShortestEdge2").get_value<float>();
            extractor.pruneRatio_maxVertexCount = extractor_settings.get_child("pruneRatio_maxVertexCount").get_value<float>();

            auto batch_settings = settings.get_child("batch");
            batch.inputFile = path_to_parentFolder + batch_settings.get_child("inputFile").get_value<std::string>();
            batch.outputFile = path_to_parentFolder + batch_settings.get_child("outputFile").get_value<std::string>();
            batch.read_batchFile = batch_settings.get_child("read_batchFile").get_value<bool>();
            batch.write_batchFile = batch_settings.get_child("write_batchFile").get_value<bool>();

            auto server_settings = settings.get_child("server");
            server.host = server_settings.get_child("host").get_value<std::string>();
            server.port = server_settings.get_child("port").get_value<int>();

            auto test_settings = settings.get_child("tests");
            if(tests.rounds_per_test == 0) tests.rounds_per_test = test_settings.get_child("rounds_per_test").get_value<uint32_t>();

            auto local_settings = settings.get_child("local");
            // No specific local settings in the JSON, but can be extended if needed
        } catch (std::exception const& e) {
            std::cerr << "Trouble parsing settings.json: " << e.what() << std::endl;
            throw std::runtime_error("Failed to parse settings file.");
        }

        if(VERBOSE_LEVEL > 0) {
            std::cout << "Settings loaded from " << settings_path << std::endl;
            std::cout << "DEBUG_LEVEL: " << DEBUG_LEVEL << std::endl;
            std::cout << "VERBOSE_LEVEL: " << VERBOSE_LEVEL << std::endl;
            std::cout << "Sea Point: (" << seaPoint.lat << ", " << seaPoint.lon << ")" << std::endl;
            std::cout << "Use Merkartor: " << (use_merkartor ? "true" : "false") << std::endl;
        }
    }

    void setBounds(std::vector<Coordinates> & coords){
        // test and correct bounds
		double min_lat = MAX_DOUBLE, min_lon = MAX_DOUBLE, max_lat = MIN_DOUBLE, max_lon = MIN_DOUBLE;
		for(auto & c : coords){
			min_lat = std::min(min_lat, c.lat);
			min_lon = std::min(min_lon, c.lon);
			max_lat = std::max(max_lat, c.lat);
			max_lon = std::max(max_lon, c.lon);
		}
		if(min_lat < boundDomain.min_lat || min_lon < boundDomain.min_lon || max_lat > boundDomain.max_lat || max_lon > boundDomain.max_lon){
			if(VERBOSE_LEVEL > 0) printf("correcting bound\n");
            if(VERBOSE_LEVEL > 1) printf("old bounds: %lf %lf %lf %lf\n", boundDomain.min_lat, boundDomain.min_lon, boundDomain.max_lat, boundDomain.max_lon);
			boundDomain.min_lat = min_lat - 0.00001 * (max_lat - min_lat);
			boundDomain.min_lon = min_lon - 0.00001 * (max_lon - min_lon);
			boundDomain.max_lat = max_lat + 0.00001 * (max_lat - min_lat);
			boundDomain.max_lon = max_lon + 0.00001 * (max_lon - min_lon);
            if(VERBOSE_LEVEL > 1) printf("new bounds: %lf %lf %lf %lf\n", boundDomain.min_lat, boundDomain.min_lon, boundDomain.max_lat, boundDomain.max_lon);
		}
    };

    void setBounds(std::vector<Point> & points){
        std::vector<Coordinates> coords;
        for(auto &p : points){coords.emplace_back(p);}
        setBounds(coords);
    };
};
Settings settings; // global settings object

#endif