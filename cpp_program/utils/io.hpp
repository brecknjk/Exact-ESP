#ifndef io_h
#define io_h

#include "types.hpp"
#include "mathLib.hpp"
#include "timer.hpp"
#include "section.hpp"

#include <boost/asio.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <sys/resource.h>
long get_mem_usage(){
    struct rusage usage;
    int ret = getrusage(RUSAGE_SELF, &usage);
    if (ret != 0) {
        std::cerr << "Error getting memory usage" << std::endl;
        return -1;
    }
    return usage.ru_maxrss; // in KB
}

void print_memUsage_KB(){
    std::cout << "Memory usage: " << get_mem_usage() << " KB" << std::endl;
}
void print_memUsage_MB(){
    std::cout << "Memory usage: " << get_mem_usage() / 1024 << " MB" << std::endl;
}
void print_memUsage_GB(){
    std::cout << "Memory usage: " << get_mem_usage() / (1024 * 1024) << " GB" << std::endl;
}

void save_to_text(std::ofstream & out_stream) {};
template<typename T, typename... Args>
void save_to_text(std::ofstream & out_stream, std::vector<T>& data, Args & ... args) {
    std::string line;	
    out_stream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);

    uint64_t dataSize = data.size();
    out_stream << dataSize << std::endl;

    for(uint64_t i = 0; i < dataSize; i++) {
        out_stream << data[i] << std::endl;
    }
    save_to_text(out_stream, args...);
};


void load_from_text(std::ifstream & in_stream){};
template<typename T, typename... Args>
void load_from_text(std::ifstream & in_stream, std::vector<T>& data, Args & ... args) {
    std::string line;	
    if (!in_stream.is_open()) throw -1;

    //read header	
	while(in_stream.peek() == '#' || in_stream.peek() == '\n'){
		std::getline(in_stream, line);
	};
    std::getline(in_stream, line);
    uint64_t dataSize = atol(line.c_str());
    data.resize(dataSize);
    
    for(uint64_t i = 0; i < dataSize; i++) {
        T value;
        in_stream >> value;
        data[i] = value;
    }
    load_from_text(in_stream, args...);
};

void load_from_mesh(std::ifstream & in_stream, std::vector<Point>& nodes, std::vector<Edge_t>& edges){
    std::string line;	
    if (!in_stream.is_open()) throw -1;

    //read header	
    while(in_stream.peek() == '#' || in_stream.peek() == '\n'){
        std::getline(in_stream, line);
    };
    std::getline(in_stream, line); //mesh
    std::getline(in_stream, line); //2
    uint32_t nodeCount = 0, faceCount = 0;
    in_stream >> nodeCount >> faceCount;

    nodes.resize(nodeCount);
    for(uint32_t n = 0; n < nodeCount; n++){
        double x, y;
        in_stream >> x >> y;
        nodes[n] = Point(x, y);
        std::getline(in_stream, line); //skip rest of the line
    }

    edges.clear();
    std::vector<uint32_t> faceToEdge(faceCount, UINT32_MAX);
    for(uint32_t f = 0; f < faceCount; f++){
        uint32_t edgePerFace = 0;
        in_stream >> edgePerFace;
        faceToEdge[f] = edges.size();
        std::vector<uint32_t> faceNodes(edgePerFace);
        //read nodes
        for(uint32_t e = 0; e < edgePerFace; e++){
            uint32_t n;
            in_stream >> n;
            faceNodes[e] = n;
        }
        std::vector<uint32_t> edgeIDs_of_face(edgePerFace);
        //read faces
        for(uint32_t e = 0; e < edgePerFace; e++){
            uint32_t f_next;
            in_stream >> f_next;
            uint32_t n0 = faceNodes[e];
            uint32_t n1 = faceNodes[(e + edgePerFace - 1) % edgePerFace];

            if(f_next == UINT32_MAX || faceToEdge[f_next] == UINT32_MAX){
                edgeIDs_of_face[e] = edges.size();
                edges.emplace_back(Edge_t(n0, n1, UINT32_MAX, UINT32_MAX));
            }else{
                for(uint32_t e_nextf = faceToEdge[f_next]; e_nextf < faceToEdge[f_next + 1]; e_nextf++){
                    uint32_t n0_e_iter = edges[e_nextf].N[0];
                    uint32_t n1_e_iter = edges[e_nextf].N[1];
                    if(n0_e_iter == n1 && n1_e_iter == n0){
                        edges[e_nextf].flip();
                        edgeIDs_of_face[e] = e_nextf;
                        break;
                    }else if(n0_e_iter == n0 && n1_e_iter == n1){
                        assert(false); // face is not ccw
                    };
                }
            }
        }
        // link edges to each other
        for(uint32_t e = 0; e < edgePerFace; e++){
            uint32_t next_e = edgeIDs_of_face[(e + 1) % edgePerFace];
            assert(edges[edgeIDs_of_face[e]].E[0] == UINT32_MAX);
            edges[edgeIDs_of_face[e]].E[0] = next_e;
        }
    }
};

void save_to_mesh(std::ofstream & out_stream, std::vector<Point>& nodes, std::vector<Edge_t>& edges) {
    printf("save_to_mesh not implemented yet, use graph function write_polyanyaStyle\n");
    assert(false); // not implemented yet
};

void save_to_bin(FILE* fp) {};
template<typename T, typename... Args>
void save_to_bin(FILE* fp, std::vector<T>& data, Args & ... args) {
    uint64_t dataSize = data.size();
    fwrite(&dataSize, sizeof(uint64_t), 1, fp);
    fwrite(data.data(), sizeof(T), dataSize, fp);
    save_to_bin(fp, args...);
};

void load_from_bin(FILE* fp) {};
template<typename T, typename... Args>
void load_from_bin(FILE* fp, std::vector<T>& data, Args & ... args) {
    uint64_t dataSize;
    std::ignore = fread(&dataSize, sizeof(uint64_t), 1, fp);
    data.resize(dataSize);
    std::ignore = fread(data.data(), sizeof(T), dataSize, fp);
    load_from_bin(fp, args...);
};

template<typename... Args>
void loadData(std::string fileName, Args & ... args){
    try {
        if(boost::algorithm::ends_with(fileName, ".bin")){
            FILE* fp = fopen(fileName.c_str(), "r+b");
            load_from_bin(fp, args...);
            fclose(fp);
        }else {
            std::ifstream in_stream (fileName);
            load_from_text(in_stream, args...);
            in_stream.close();
        };
    }catch (...) {
        std::cout << "coudn't read " << fileName << std::endl;
        throw std::runtime_error("Failed to load data from " + fileName);
    }
};

template<typename... Args>
void saveData(std::string fileName, Args & ... args){
    try {
        if(boost::algorithm::ends_with(fileName, ".bin")){
            FILE* fp = fopen(fileName.c_str(), "w+b");
            save_to_bin(fp, args...);
            fclose(fp);
        }else if(boost::algorithm::ends_with(fileName, ".mesh")){
            std::ofstream out_stream (fileName);
            //save_to_mesh(out_stream, args...);
            out_stream.close();
        }else{
            std::ofstream out_stream (fileName);
            save_to_text(out_stream, args...);
            out_stream.close();
        };
    }catch (...) {
        std::cout << "coudn't read " << fileName << std::endl;
    }
};

void inline write_coastLine(const std::string & outputFile, const std::vector<Coordinates> & coords, const std::vector<std::vector<uint32_t>> & segments, bool printClosed = false){
    try{
        std::ofstream outCoastlines (outputFile);
        outCoastlines << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outCoastlines << segments.size() << '\n';
        for(const std::vector<uint32_t> & s : segments){
            outCoastlines << s.size() << '\n';
            for(uint32_t n : s){
                outCoastlines << coords[n] << '\n';
            }
            if(printClosed) outCoastlines << coords[s[0]] << '\n';
        }
        outCoastlines.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

void inline read_coastLine(const std::string & inputFile, std::vector<Coordinates> & coords, std::vector<std::vector<uint32_t>> & segments){
    try{
        std::ifstream inCoastlines (inputFile);
        if (!inCoastlines.is_open()) throw std::runtime_error("can't open file: " + inputFile);
        uint32_t segmentCount, segmentSize;
        inCoastlines >> segmentCount;
        segments.resize(segmentCount);
        coords.clear();
        for(uint32_t i = 0; i < segmentCount; i++){
            inCoastlines >> segmentSize;
            for(uint32_t j = 0; j < segmentSize; j++){
                Coordinates c;
                inCoastlines >> c;
                coords.push_back(c);
            }
        }
        inCoastlines.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << inputFile << std::endl;
        assert(false); // abort if no release
    }
}

void inline read_coastLine(const std::string & inputFile, std::vector<Coordinates> & coords, std::vector<Edge> & simpleEdges){
    try{
        std::ifstream inCoastlines (inputFile);
        if (!inCoastlines.is_open()) throw std::runtime_error("can't open file: " + inputFile);
        uint32_t segmentCount, segmentSize;
        Coordinates c;
        inCoastlines >> segmentCount;
        coords.clear();
        simpleEdges.clear();
        for(uint32_t i = 0; i < segmentCount; i++){
            assert(inCoastlines.eof() == false);
            inCoastlines >> segmentSize;
            assert(segmentSize > 1);

            for(uint32_t j = 0; j < segmentSize; j++){
                assert(inCoastlines.eof() == false);
                inCoastlines >> c;
                coords.push_back(c);
                simpleEdges.push_back(Edge(coords.size() - 1, coords.size() + (segmentSize <= j + 1 ? - segmentSize : 0)));
            }
        }
        inCoastlines.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << inputFile << std::endl;
        assert(false); // abort if no release
    }
}

void inline read_routeFile(const std::string & inputFile, std::list<Coordinates> & coords){
    try{
        std::ifstream inSt (inputFile);
        if (!inSt.is_open()) throw std::runtime_error("can't open file: " + inputFile);
        while(inSt.eof() == false){
            Coordinates c;
            inSt >> c.lat >> c.lon;
            coords.emplace_back(c);
        }
        inSt.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading route file: " << inputFile << std::endl;
        assert(false); // abort if no release
    }
}

std::vector<std::pair<Coordinates, Coordinates>> readScenarioFile(const std::string & inputFile, std::vector<double> & dist){
    std::vector<std::pair<Coordinates, Coordinates>> startGoalPairs;
    try{
        std::string line;
        std::ifstream inScen (inputFile);
        if (!inScen.is_open()) throw std::runtime_error("can't open file: " + inputFile);
        getline(inScen, line); // skip version line

        while(!inScen.eof()){
        	uint32_t id;
            std::string mapString;
            Coordinates start, goal;
            double distance, mapWidth, maxHeight;
            inScen >> id >> mapString >> mapWidth >> maxHeight; // dump
            if(inScen.eof()) break;
            inScen >> start.lat >> start.lon >> goal.lat >> goal.lon >> distance;
            startGoalPairs.emplace_back(start, goal);
            dist.push_back(distance);
        }
        
        inScen.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << inputFile << std::endl;
        assert(false); // abort if no release
    }
    return startGoalPairs;
}

std::vector<std::pair<Coordinates, Coordinates>> inline read_batchFile(const std::string & inputFile){
    if(boost::algorithm::ends_with(inputFile, ".scen")){
        std::vector<double> dist;
        return readScenarioFile(inputFile, dist);
    }
    std::vector<std::pair<Coordinates, Coordinates>> startGoalPairs;
    try{
        std::ifstream inSt (inputFile);
        if (!inSt.is_open()) throw std::runtime_error("can't open file: " + inputFile);
        while(inSt.eof() == false){
            Coordinates start, goal;
            inSt >> start.lat >> start.lon >> goal.lat >> goal.lon;
            std::swap(start.lat, start.lon);
            std::swap(goal.lat, goal.lon);  
            if(settings.use_merkartor) start = rev_bounded_mercator(start);
            if(settings.use_merkartor) goal = rev_bounded_mercator(goal);
            startGoalPairs.emplace_back(start, goal);
        }
        inSt.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << inputFile << std::endl;
        assert(false); // abort if no release
    }
    return startGoalPairs;
}

void inline write_glFile(const std::string & outputFile, std::vector<Coordinates> & c, const std::vector<std::pair<Edge_t, uint8_t>> & e){
    double minLat = MAX_DOUBLE, minLon = MAX_DOUBLE, maxLat = -MAX_DOUBLE, maxLon = -MAX_DOUBLE;
    for(auto & coord : c){
        minLat = std::min(minLat, coord.lat);
        minLon = std::min(minLon, coord.lon);
        maxLat = std::max(maxLat, coord.lat);
        maxLon = std::max(maxLon, coord.lon);
    }
    if(minLat < -181 || minLon < -181 || maxLat > 181 || maxLon > 181){
        auto size = std::max(maxLon - minLon, maxLat - minLat);
        for(uint32_t i = 0; i < c.size(); i++){
            c[i].lat = (c[i].lat - minLat) * 10 / size;
            c[i].lon = (c[i].lon - minLon) * 10 / size;
        }
    }
    try{
        std::ofstream outStream (outputFile);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << c.size() << '\n';
        outStream << e.size() << '\n';
        for(const Coordinates & coord : c){
            outStream << coord.lat << " " << coord.lon << '\n';
        }
        for(auto & edge_pair : e){
            if(edge_pair.first.N[0] == UINT32_MAX || edge_pair.first.N[1] == UINT32_MAX) continue; // skip invalid edges
            assert(edge_pair.second > 0); // color should not be 0
            outStream << edge_pair.first.N[0] << " " << edge_pair.first.N[1] << " " << 1 << " " << (int) edge_pair.second << '\n';
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble writing to " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

void inline write_glFile(const std::string & outputFile, std::vector<Coordinates> & c, const std::vector<Edge_t> & e, uint8_t color = 1){
    std::vector<std::pair<Edge_t, uint8_t>> e_with_color(e.size());
    for(uint32_t i = 0; i < e.size(); i++){
        e_with_color[i] = std::make_pair(e[i], color);
    }
    write_glFile(outputFile, c, e_with_color);
}

void inline write_route_to_glFile(const std::string & outputFile, const std::list<Coordinates> & c, uint8_t color = 5){
    if(c.size() == 0) return;
    if(settings.VERBOSE_LEVEL > 0){
        std::cout << "route with " << c.size() << " hops: ";
        for(Coordinates coord : c){
            if(settings.use_merkartor) coord = bounded_mercator(coord);
            std::cout << "(" << coord.lat << "," << coord.lon << ") ";
        }
        std::cout << std::endl;
    }
    try{
        std::ofstream outStream (outputFile);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << c.size() << '\n';
        outStream << c.size() - 1 << '\n';
        for(const Coordinates & coord : c){
            outStream << coord.lat << " " << coord.lon << '\n';
        }
        for(uint32_t i = 0; i < c.size() - 1; i++){
            outStream << i << " " << i + 1 << " " << 3 << " " << color << '\n';
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble writing to " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

void inline write_Searchspace_to_glFile(const std::string & outputFile, const std::list<std::pair<Coordinates, Coordinates>> & sp, int color = 4){
    try{
        std::ofstream outStream (outputFile);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << sp.size() * 2 << '\n';
        outStream << sp.size() << '\n';
        for(auto [c0, c1] : sp){
            outStream << c0.lat << " " << c0.lon << '\n' << c1.lat << " " << c1.lon << '\n';
        }
        for(uint32_t e = 0; e < sp.size(); e++){
            outStream << e*2 << " " << e*2+1 << " " << 2 << " " << color << '\n';
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble writing to " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

void inline write_SearchspaceAndRoute_to_glFile(const std::string & outputFile, const std::list<std::pair<Coordinates, Coordinates>> & sp, const std::list<Coordinates> & route, uint8_t color = 4){
    try{
        std::ofstream outStream (outputFile);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << sp.size() * 2 + route.size() << '\n';
        outStream << sp.size() + route.size() - 1 << '\n';
        for(auto [c0, c1] : sp){
            outStream << c0.lat << " " << c0.lon << '\n' << c1.lat << " " << c1.lon << '\n';
        }
        for(const Coordinates & coord : route){
            outStream << coord.lat << " " << coord.lon << '\n';
        }
        for(uint32_t e = 0; e < sp.size(); e++){
            outStream << e*2 << " " << e*2+1 << " " << 2 << " " << color << '\n';
        }
        for(uint32_t i = 0; i < route.size() - 1; i++){
            outStream << i + sp.size() * 2 << " " << i + sp.size() * 2 + 1 << " " << 3 << " " << color + 1 << '\n';
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble writing to " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

// s.bucket >> map >> s.xsize >> s.ysize >> s.start.x >> s.start.y >> s.goal.x >> s.goal.y >> s.gridcost
void createScenarioFile(const std::string & outputFile, const std::vector<std::pair<Coordinates, Coordinates>> & points){
    Coordinates mercator_min_min = Coordinates(-89.0, -179.0);
    Coordinates mercator_min_max = Coordinates(-89.0, 179.0);
    Coordinates mercator_max_min = Coordinates(89.0, -179.0);
    Coordinates mercator_max_max = Coordinates(89.0, 179.0);
    if(settings.use_merkartor){
        mercator_min_min = bounded_mercator(mercator_min_min);
        mercator_min_max = bounded_mercator(mercator_min_max);
        mercator_max_min = bounded_mercator(mercator_max_min);
        mercator_max_max = bounded_mercator(mercator_max_max);
    }
    double min_x = std::min(mercator_min_min.lat, mercator_max_min.lat);
    double max_x = std::max(mercator_min_max.lat, mercator_max_max.lat);
    double min_y = std::min(mercator_min_min.lon, mercator_min_max.lon);
    double max_y = std::max(mercator_max_min.lon, mercator_max_max.lon);
    if(max_y < min_y) std::swap(max_y, min_y);
    double xsize = max_x - min_x;
    double ysize = max_y - min_y;
    try{
        std::ofstream outStream (outputFile, std::ios::trunc);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << "version 1\n";
        for(uint32_t i = 0; i < points.size(); i++){
            auto s_p = settings.use_merkartor ? bounded_mercator(points[i].first) : points[i].first;
            auto t_p = settings.use_merkartor ? bounded_mercator(points[i].second) : points[i].second;
            outStream << "0\texample.mesh\t" << xsize << '\t' << ysize << '\t' << s_p.lat << '\t' << s_p.lon << '\t' << t_p.lat << '\t' << t_p.lon << "\t1.0\n";
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble writing to " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

void writePathsToFile(const std::string & outputFile, const std::vector<std::list<Coordinates>> & paths){
    if(settings.VERBOSE_LEVEL > 0 && settings.use_merkartor) printf("Warning: write path without rev_mercartor");
    try{
        std::ofstream outStream (outputFile, std::ios::trunc);
        outStream << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
        outStream << paths.size() << std::endl;
        for(auto & path : paths){
            outStream << path.size() << " ";
            for(auto & p : path){
                outStream << p.lat << " " << p.lon << " ";
            }
            outStream << std::endl;
        }
        outStream.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble writing to " << outputFile << std::endl;
        assert(false); // abort if no release
    }
}

std::vector<std::list<Coordinates>> readPathsFromFile(const std::string & inputFile){
    std::vector<std::list<Coordinates>> paths;
    try{
        std::string line;
        std::ifstream inScen (inputFile);
        if (!inScen.is_open()) throw std::runtime_error("can't open file: " + inputFile);
        uint32_t count = 0;
        inScen >> count;
        paths.reserve(count);
        for(uint32_t i = 0; i < count; i++){
            paths.push_back(std::list<Coordinates>());
            uint32_t pathLength = 0;
            inScen >> pathLength;
            for(uint32_t j = 0; j < count; i++){
                Coordinates p;
                inScen >> p.lat >> p.lon;
            }
        }
        inScen.close();
    }catch(std::exception& e){
        std::cerr << e.what() << std::endl;
        std::cerr << "error: trouble reading batch file: " << inputFile << std::endl;
        assert(false); // abort if no release
    }
    return paths;
}

#endif