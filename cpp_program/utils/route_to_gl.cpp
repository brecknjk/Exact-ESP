#include "types.hpp"
#include "io.hpp"
#include "mathLib.hpp"
#include "settings.hpp"

void help_route_to_gl() {
	std::cout << "USAGE: <infile> <outfile> <OPTIONAL_FLAGS>" << std::endl;
	std::cout << "OPTIONAL:" << std::endl;
	std::cout << "set settings via settings.json" << std::endl;
}

int main(int argc, char ** argv) {
	if (argc < 3) {
		help_route_to_gl();
		return -1;
	}
	std::string inputFile, outputFile, path_to_settings = "settings.json";

	try {
		inputFile = argv[1];
		assert(access( inputFile.c_str(), F_OK ) != -1 );
		outputFile = argv[2];

		for (int i = 3; i < argc; ++i) {
			std::string token(argv[i]);
			if (token == "-h" || token == "--help") {
				help_route_to_gl();
				return 0;
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
			} else {
				throw std::runtime_error("Unknown argument: " + token);
			}
		}

		settings.readSettings(path_to_settings);
	} catch (std::exception const& e) {
		std::cerr << e.what() << std::endl;
		help_route_to_gl();
		return -1;
	}

	Section section;
    std::list<Coordinates> coords;

    section = Section("read route file");
    read_routeFile(inputFile, coords);
    section.end_section();

    for(auto & c: coords){
        if(settings.use_merkartor) c = rev_bounded_mercator(c);
    }

    section = Section("write route file");
    write_route_to_glFile(outputFile, coords);
    section.end_section();

	return 0;
}