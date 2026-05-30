#include "osmExtractor.hpp"
#include "triangulator.hpp"
#include "specialGraphConstructor.hpp"

#include "../utils/types.hpp"
#include "../utils/io.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/settings.hpp"

std::vector<std::vector<uint32_t>> connectEdges(std::vector<Edge> & edges, uint32_t nodeCopunt){
	auto lookup = std::vector<std::pair<uint32_t, uint32_t>>(nodeCopunt, {UINT32_MAX, UINT32_MAX});
	for(auto & e : edges){
		assert(lookup[e.source].second == UINT32_MAX);
		assert(lookup[e.target].second == UINT32_MAX);
		lookup[e.source] = {e.target, lookup[e.source].first};
		lookup[e.target] = {e.source, lookup[e.target].first};
	}
	edges.clear();

	std::vector<std::vector<uint32_t>> segments;
	for(uint32_t n = 0; n < nodeCopunt; n++){
		if(lookup[n].first == UINT32_MAX && lookup[n].second == UINT32_MAX) continue;
		assert((lookup[n].first == UINT32_MAX ) == (lookup[n].second == UINT32_MAX));


		uint32_t current_node = n, lastNode = UINT32_MAX;
		std::vector<uint32_t> segment;
		// not closed!
		// segment.push_back(current_node);
		do{
			assert(current_node != UINT32_MAX);
			assert((lookup[current_node].first != UINT32_MAX ) || (lookup[current_node].second != UINT32_MAX));
			// edge case if not closed -> close it
			if((lookup[current_node].first == UINT32_MAX ) || (lookup[current_node].second == UINT32_MAX)){
				segment.push_back(n);
				lookup[current_node] = {UINT32_MAX, UINT32_MAX};
				break;
			}
			uint32_t next_node = (lookup[current_node].first == lastNode) ? lookup[current_node].second : lookup[current_node].first;

			lookup[current_node] = {UINT32_MAX, UINT32_MAX};
			lastNode = current_node;
			current_node = next_node;

			segment.push_back(current_node);
		} while(current_node != n);
		segments.push_back(segment);
	}
	return segments;
}

void insertBoundingBoxEdges(std::vector<Coordinates> & coords, std::vector<Edge> & simple_edges){
	settings.setBounds(coords);
	// Insert the corners of the bounding box
	coords.push_back(Coordinates(settings.boundDomain.min_lat, settings.boundDomain.min_lon));
	coords.push_back(Coordinates(settings.boundDomain.min_lat, settings.boundDomain.max_lon));
	coords.push_back(Coordinates(settings.boundDomain.max_lat, settings.boundDomain.max_lon));
	coords.push_back(Coordinates(settings.boundDomain.max_lat, settings.boundDomain.min_lon));

	// Insert the vertices of the bounding box into the triangulation
	simple_edges.push_back(Edge(coords.size() - 4, coords.size() - 3));
	simple_edges.push_back(Edge(coords.size() - 3, coords.size() - 2));
	simple_edges.push_back(Edge(coords.size() - 2, coords.size() - 1));
	simple_edges.push_back(Edge(coords.size() - 1, coords.size() - 4));
}

void help_extractor() {
	std::cout << "USAGE: <infile> <outfile> <OPTIONAL_FLAGS>" << std::endl;
	std::cout << "if generating worstCasePolyAnya: use -wcPa instead of <infile>" << std::endl;
	std::cout << "OPTIONAL:" << std::endl;
	std::cout << "-swcPa <number>: size of instance, if genrating worst case for polyanya" << std::endl;
	std::cout << "--connect: if coastlines should be connected" << std::endl;
	std::cout << "-scl <second_outfile>" << std::endl;
	std::cout << "set boundingbox via settings.json" << std::endl;
	std::cout << "set seapoint via settings.json" << std::endl;
}

int main(int argc, char ** argv) {
	if (argc < 3) {
		help_extractor();
		return -1;
	}
	std::string inputFile, outputFile, secondOutputFile = "", path_to_settings = "settings.json";
	bool connect = false, worstCasePolyAnya = false;
	uint32_t sizeOfWorstCasePolyAnya = 100;

	try {
		inputFile = argv[1];
		outputFile = argv[2];

		for (int i = 3; i < argc; ++i) {
			std::string token(argv[i]);
			if (token == "-h" || token == "--help") {
				help_extractor();
				return 0;
			} else if (token == "-swcPa" && i + 1 < argc) {
				sizeOfWorstCasePolyAnya = std::stoi(argv[++i]);
			} else if (token == "-s" && i + 1 < argc) {
				path_to_settings = std::string(argv[++i]);
			} else if (token == "-v" && i + 1 < argc){
				settings.VERBOSE_LEVEL = std::stoi(argv[++i]);
			} else if (token == "-d" && i + 1 < argc){
				settings.DEBUG_LEVEL = std::stoi(argv[++i]);
			} else if (token == "-m"){
				settings.use_merkartor  = true;
			} else if (token == "-fp" && i + 1 < argc) {
				settings.path_to_parentFolder = std::string(argv[++i]);
			} else if (token == "--connect") {
				connect = true;
			} else if (token == "-scl" && i + 1 < argc) {
				secondOutputFile = argv[++i];
				assert(!boost::algorithm::ends_with(secondOutputFile, ".pbf"));
			} else {
				throw std::runtime_error("Unknown argument: " + token);
			}
		}
		settings.readSettings(path_to_settings);

		if(inputFile == "-wcPa"){
			worstCasePolyAnya = true;
			assert(connect == false);
			assert(settings.use_merkartor == false);
		} else{
			assert(access( inputFile.c_str(), F_OK ) != -1 );
		}
	} catch (std::exception const& e) {
		std::cerr << e.what() << std::endl;
		help_extractor();
		return -1;
	}

	Section section;
	std::vector<Edge> simple_edges;
	std::vector<Coordinates> coords;

	if (boost::algorithm::ends_with(inputFile, ".pbf")){
		std::vector<std::string> inputFileNames = {argv[1]}; // need list for some strange reason

		section = Section("Gathering Nodes and Edges");
		gatherEdgesAndNodes(inputFileNames, coords, simple_edges);
		section.end_section();

		if(secondOutputFile != ""){
			section = Section("save Nodes and Edges to second file");
			saveData(secondOutputFile, coords, simple_edges);
			section.end_section();
		}
	}else if(worstCasePolyAnya) {
		section = Section("constructing worst case scenario for Polyanya [O(n³) vs O(n*log(n))]");
		constructWorstCasePolyAnyaGraph(coords, simple_edges, sizeOfWorstCasePolyAnya);
		section.end_section();
	}else if(boost::algorithm::ends_with(inputFile, "polygons.fmi")){
		section = Section("read polygons in fmi style to: " + inputFile);
		read_coastLine(inputFile, coords, simple_edges);
		section.end_section();
	} else {
		section = Section("read nodes and edges from file to: " + inputFile);
		loadData(inputFile, coords, simple_edges);
		section.end_section();
	}

	if(settings.VERBOSE_LEVEL > 0) printf("loaded %li coords and %li edges\n", coords.size(), simple_edges.size());

	if (connect) {
		section = Section("connect Edges to Segments");
		auto segments = connectEdges(simple_edges, coords.size());
		section.end_section();

		if (boost::algorithm::ends_with(outputFile, ".fmi")){
			section = Section("print Segments in fmi style");
			write_coastLine(outputFile, coords, segments);
			section.end_section();
		}
	} else {
		if(!worstCasePolyAnya) insertBoundingBoxEdges(coords, simple_edges);
		if(worstCasePolyAnya) settings.seaPoint = Coordinates(1, 1);

		section = Section("calc constrained DT");
		auto edges = calc_constrained_DT(coords, simple_edges, settings.seaPoint);
		section.end_section();
		simple_edges.clear();

		section = Section("removing nodes with deg = 0");
		std::vector<uint32_t> node_permutation(coords.size(), 0);
		for(uint32_t e = 0; e < edges.size() * 2; e++){
			node_permutation[edges[e >> 1].N[(e & 1)]]++;
		}
		uint32_t newNodeCount = 0;
		for(uint32_t n = 0; n < node_permutation.size(); n++){
			node_permutation[n] = node_permutation[n] == 0 ? UINT32_MAX : newNodeCount++;
		}
		for(uint32_t n = 0; n < coords.size(); n++){
			if(node_permutation[n] != UINT32_MAX) coords[node_permutation[n]] = coords[n];
		}
		coords.resize(newNodeCount);
		for(uint32_t e = 0; e < edges.size(); e++){
			assert(node_permutation[edges[e].N[0]] != UINT32_MAX);
			edges[e].N[0] = node_permutation[edges[e].N[0]];
			assert(node_permutation[edges[e].N[1]] != UINT32_MAX);
			edges[e].N[1] = node_permutation[edges[e].N[1]];
		}
		if(settings.VERBOSE_LEVEL > 0) printf("removed %li nodes with degree = 0\n", node_permutation.size() - newNodeCount);
		node_permutation.clear();
		section.end_section();

        section = Section("let edges point to one side of other edge, via id = id * 2 + direction");
        for(uint32_t i = 0; i < edges.size(); i++){
            if(edges[i].E[0] != UINT32_MAX) edges[i].E[0] = edges[i].E[0] * 2 + (1 ^ edges[edges[i].E[0]].getIndex_n(edges[i].N[0]));
            if(edges[i].E[1] != UINT32_MAX) edges[i].E[1] = edges[i].E[1] * 2 + (1 ^ edges[edges[i].E[1]].getIndex_n(edges[i].N[1]));
        }
        section.end_section();

		for(uint32_t e_id = 0; e_id < edges.size() * 2; e_id++){
       		if(edges[e_id >> 1].E[e_id & 1] == UINT32_MAX) continue;

			uint32_t n0 = edges[e_id >> 1].N[1 ^ e_id & 1];
			uint32_t n1 = edges[e_id >> 1].N[e_id & 1];
			auto next_e = edges[e_id >> 1].E[e_id & 1];
			uint32_t n2 = edges[next_e >> 1].N[next_e & 1];

			assert(n0 != n1 && n1 != n2 && n0 != n2);
			assert(edges[next_e >> 1].N[1 ^ next_e & 1] == n1);
			assert(!CGAL::right_turn(Point(coords[n0].lat, coords[n0].lon), Point(coords[n1].lat, coords[n1].lon), Point(coords[n2].lat, coords[n2].lon)));
		}

		if(boost::algorithm::ends_with(outputFile, ".gl")) {
			section = Section("print Edges in gl style");
			write_glFile(outputFile, coords, edges);
			section.end_section();
		} else {
			section = Section("print nodes and linked edges");
			std::vector<Point> points;
			points.reserve(coords.size());
			for(const Coordinates & c : coords){
				points.emplace_back(c.lat, c.lon);
			}
			saveData(outputFile, points, edges);
			points.clear();
			section.end_section();
		}
	}
	return 0;
}