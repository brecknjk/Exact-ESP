#ifndef triangulator_h
#define triangulator_h

#define CGAL_NO_CDT_2_WARNING
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_conformer_2.h>

// #include <CGAL/draw_triangulation_2.h>

#include "../utils/types.hpp"
#include "../utils/mathLib.hpp"
#include "../utils/io.hpp"

struct Face_info {
    bool in_domain = false; // true if face is in domain
    uint32_t edgeID = UINT32_MAX;
    Face_info() {}
    Face_info(bool in_domain) : in_domain(in_domain) {}
    void set_in_domain(bool in_domain) {this->in_domain = in_domain;}
    bool is_in_domain() const {return in_domain;}
    void set_edgeID(uint32_t edgeID) {
        this->edgeID = edgeID;
    }
    uint32_t get_edgeID() const {
        return edgeID;
    }
};
typedef CGAL::Exact_intersections_tag                                   Itag;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info, K>         Fb2;
typedef CGAL::Triangulation_vertex_base_with_info_2<uint32_t, K>        Vb2;
typedef CGAL::Constrained_triangulation_face_base_2<K, Fb2>             CFb2;
typedef CGAL::Triangulation_data_structure_2<Vb2,CFb2>                  Tds2;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds2, Itag>       CDT;

typedef CDT::Point                                                      Point_CDT;
typedef CDT::Face_handle                                                Face_handle;

uint32_t markInDomain(CDT & cdt, const Point & pointInSea) {
    std::list<Face_handle> queue;
    uint32_t totalEdgesInDomain = 0;
    // mark all faces as not in domain
    for(CDT::Face_handle f : cdt.all_face_handles()){
        f->info().set_in_domain(false);
    }
    // put first face in queue
    queue.push_back(cdt.locate(pointInSea));
    while(!queue.empty()){
        Face_handle fh = queue.front();
        queue.pop_front();
        if(!fh->info().is_in_domain()){
            fh->info().set_in_domain(true);

            for(uint32_t i = 0; i < 3; i++){
                Face_handle fi = fh->neighbor(i);
                if(!fi->info().is_in_domain()){
                    assert(fi->info().is_in_domain() == false);
                    totalEdgesInDomain++;
                    if(!cdt.is_constrained(CDT::Edge(fh,i))) queue.push_back(fi);
                }
            }
        }
    }
    return totalEdgesInDomain;
}

void refineTriangulation(const CDT &cdt, std::vector<Point> & toInsert, std::vector<double> & radius2_toInsert, double longestEdge2) {
    for (CDT::Face_handle fh : cdt.all_face_handles()) {
	    if (!fh->info().is_in_domain() || (cdt.is_infinite(fh))) continue;
        //if ((cdt.is_infinite(fh))) continue;

        // idea: do not refine triangles where the longest edge is a boundary
        // as this might lead to even skinnier triangles
        double minLen2 = MAX_DOUBLE;
        double maxLen2 = 0.0;
        bool doRefine = true;
        for(uint32_t i = 0; i < 3; i++){
            double edgeLen2 = CGAL::to_double(K::Segment_2(fh->vertex((i+1) % 3)->point(),fh->vertex((i+2) % 3)->point()).squared_length());
            if(maxLen2 < edgeLen2) doRefine = !cdt.is_constrained(CDT::Edge(fh,i));
            minLen2 = std::min(minLen2, edgeLen2);
            maxLen2 = std::max(maxLen2, edgeLen2);
        }
        if(!doRefine) continue;

        K::Circle_2 myCircle(fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point());	
		Point center = myCircle.center();
		double radius2 = CGAL::to_double(myCircle.squared_radius());

        //if((radius2 < 2.0 * minLen2) && (maxLen2 <= longestEdge2)) continue;
        if(maxLen2 <= longestEdge2) continue;

        //Point center_rev_mercator = (settings.use_merkartor) ? Point(rev_bounded_mercator(center)) : center;
        if(settings.boundDomain.min_lat > center.x() || settings.boundDomain.max_lat < center.x()) continue;
        if(settings.boundDomain.min_lon > center.y() || settings.boundDomain.max_lon < center.y()) continue;
        
        toInsert.emplace_back(center);
        radius2_toInsert.push_back(radius2);
    }
};

void pruneSteiner(std::vector<Point> &toPrune, std::vector<double> & radius2_toPrune, int GRIDSIZE = 1000){
	std::vector<Point> survivors;
    Coordinates c_min(settings.boundDomain.min_lat, settings.boundDomain.min_lon);
    Coordinates c_max(settings.boundDomain.max_lat, settings.boundDomain.max_lon);
    c_min = (settings.use_merkartor) ? bounded_mercator(c_min) : c_min;
    c_max = (settings.use_merkartor) ? bounded_mercator(c_max) : c_max;

	double xDelta = (c_max.lat - c_min.lat) / GRIDSIZE;
	double yDelta = (c_max.lon - c_min.lon) / GRIDSIZE;
	// we have a one-dimensional vector of vectors index (simulating 2-dim vector)
	std::vector<std::vector<uint32_t>> grid(GRIDSIZE * GRIDSIZE, std::vector<uint32_t>()); 
	// sort toPrune according to squared radius?
    std::vector<uint32_t> indices(toPrune.size());
    std::iota (std::begin(indices), std::end(indices), 0); // Fill with 0, 1, ...
    std::cout << "Before pruning: " << toPrune.size() << std::endl;
	std::sort(indices.begin(), indices.end(), [&radius2_toPrune](uint32_t &left, uint32_t &right) { return radius2_toPrune[left] < radius2_toPrune[right]; });
	std::cout << "Sorted Steiner points according to radius2" << std::endl;

    // check if there is already some other point at distance sqrt(sqRadius)
    auto hasCollision = [&] (uint32_t i, int xGrid, int yGrid) {
        // go over all relevant grid cells and check all points in there
        int gridRadius = uint32_t(sqrt(radius2_toPrune[i]) / std::min(xDelta, yDelta)) + 1;
		for(uint32_t x = std::max(0, xGrid-gridRadius); x <= std::min(GRIDSIZE-1, xGrid+gridRadius); x++){
			for(uint32_t y = std::max(0, yGrid-gridRadius); y <= std::min(GRIDSIZE-1, yGrid+gridRadius); y++){
				uint32_t cellNr = GRIDSIZE*y + x;
				for(uint32_t j : grid[cellNr]){
					double sqDist = CGAL::to_double(K::Segment_2(toPrune[i], toPrune[j]).squared_length());
					if ((sqDist < radius2_toPrune[i]) && (sqDist < radius2_toPrune[j])) return true;
				}
			}
        }
        return false;
    };

	for(uint32_t i : indices){
		int xGrid = std::min(std::max(int((CGAL::to_double(toPrune[i].x()) - c_min.lat) / xDelta), 0), GRIDSIZE-1);
		int yGrid = std::min(std::max(int((CGAL::to_double(toPrune[i].y()) - c_min.lon) / yDelta), 0), GRIDSIZE-1);
        bool foundTooClose = hasCollision(i, xGrid, yGrid);
		
		if (!foundTooClose){
			grid[GRIDSIZE*yGrid + xGrid].push_back(i);
			survivors.push_back(toPrune[i]);
		}
	}
	toPrune = survivors;
	std::cout << "After pruning: " << toPrune.size() << std::endl;
};

// assumes bounding box is already inserted
std::vector<Edge_t> calc_constrained_DT(std::vector<Coordinates> & coords, std::vector<Edge> & constrains, const Coordinates & seaCoords) {
    Point seaPoint(seaCoords.lat, seaCoords.lon);
    if(settings.use_merkartor) seaPoint = Point(bounded_mercator(seaPoint));
    
    CDT cdt;

    std::vector<std::pair<uint32_t, uint32_t>> constrainsIndices;
    constrainsIndices.reserve(constrains.size());
    for(const Edge & e : constrains){
        constrainsIndices.emplace_back(e.source, e.target);
    }
    std::vector<Point> points;
    points.reserve(coords.size());
    for(const Coordinates & c : coords){
        points.emplace_back(c.lat, c.lon);
        if(settings.use_merkartor) points.back() = bounded_mercator(points.back());
    }
    settings.setBounds(points);

    cdt.insert_constraints(points.begin(), points.end(), constrainsIndices.begin(), constrainsIndices.end());
    assert(cdt.is_valid());
    std::cout << "inserted " << constrains.size() << " constraints" << std::endl;

    constrainsIndices.clear();
    points.clear();

    // if(settings.DEBUG_LEVEL > 2) CGAL::draw(cdt);

    // std::cout << "Number of vertices before: " << cdt.number_of_vertices() << std::endl;
    // // make it conforming Delaunay
    // CGAL::make_conforming_Delaunay_2(cdt);
    // std::cout << "Number of vertices after make_conforming_Delaunay_2: " << cdt.number_of_vertices() << std::endl;

    // CGAL::draw(cdt);

    // //https://doc.cgal.org/latest/Mesh_2/index.html
    // // // then make it conforming Gabriel
    // CGAL::make_conforming_Gabriel_2(cdt);
    // std::cout << "Number of vertices after make_conforming_Gabriel_2: "<< cdt.number_of_vertices() << std::endl;

    // CGAL::draw(cdt);

    if(settings.extractor.pruneTriangulation){
        auto section = Section("refine triangulation");
        uint32_t initaleVertexCount = cdt.number_of_vertices();
        uint32_t oldVertexCount = initaleVertexCount;
        uint32_t MAX_ROUNDS = 100;
        double size_domain = (settings.boundDomain.max_lat - settings.boundDomain.min_lat) * (settings.boundDomain.max_lon - settings.boundDomain.min_lon);

        double shortestEdge2 = MAX_DOUBLE;
        for(CDT::Edge e : cdt.finite_edges()){
            double edgeLen2 = CGAL::to_double(cdt.segment(e).squared_length());
            if(shortestEdge2 > edgeLen2 && edgeLen2 > 0) shortestEdge2 = edgeLen2;
        }

        for(uint32_t i = 0; i < MAX_ROUNDS; i++){
            if(float(oldVertexCount) / float(initaleVertexCount) >= float(settings.extractor.pruneRatio_maxVertexCount)) break;
            markInDomain(cdt, seaPoint);

            std::vector<Point> toInsert;
            std::vector<double> radius2_toInsert;
            std::cout << "Longest sqrd edge " << size_domain / initaleVertexCount * settings.extractor.pruneRatio_longestToShortestEdge2 << std::endl;
            refineTriangulation(cdt, toInsert, radius2_toInsert, size_domain / initaleVertexCount * settings.extractor.pruneRatio_longestToShortestEdge2);
            std::cout << "Found " << toInsert.size() << " potential refinement Points " << std::endl;
            pruneSteiner(toInsert, radius2_toInsert, ceil(sqrt(initaleVertexCount / 100)));

            cdt.insert(toInsert.begin(), toInsert.end());
            assert(cdt.is_valid());
            std::cout << "After insertion of " << toInsert.size() << " Steiner points we have " << cdt.number_of_vertices() << std::endl;
                
            // if(settings.DEBUG_LEVEL > 2) CGAL::draw(cdt);

            if(cdt.number_of_vertices() == oldVertexCount) break;
            oldVertexCount = cdt.number_of_vertices();
        }
        section.end_section();

        // if(settings.DEBUG_LEVEL > 2) CGAL::draw(cdt);
    }

	uint32_t totalEdgesInDomain = markInDomain(cdt, seaPoint);
    std::cout << "marked " << totalEdgesInDomain << " edges in domain" << std::endl;
    
    coords.clear();
    coords.reserve(cdt.number_of_vertices());
    for(CDT::Vertex_handle vh : cdt.finite_vertex_handles()){
        vh->info() = coords.size();
        if(settings.use_merkartor){
            coords.emplace_back(rev_bounded_mercator(Coordinates(vh->point())));
        }else{
            coords.emplace_back(vh->point());
        }
    }

    std::vector<Edge_t> edgesToDraw(totalEdgesInDomain);
    uint32_t constrainedEdgeCount = 0, otherEdgeCount = 0;
    for (CDT::Face_handle fh : cdt.all_face_handles()){
        if (!fh->info().is_in_domain() || cdt.is_infinite(fh)) continue;
        assert(fh->info().get_edgeID() == UINT32_MAX);
        uint32_t index[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
        
        for(uint32_t i = 0; i < 3; i++){
            if(fh->neighbor(i)->info().edgeID == UINT32_MAX){
                index[i] = cdt.is_constrained(CDT::Edge(fh,i)) ? constrainedEdgeCount++ : edgesToDraw.size() - ++otherEdgeCount;
                edgesToDraw[index[i]] = Edge_t(fh->vertex((i+1)%3)->info(), fh->vertex((i+2)%3)->info());
            }else{
                auto neighbor = fh->neighbor(i);
                index[i] = neighbor->info().get_edgeID();
                for(uint32_t j = 0; j < 3; j++){
                    if(!fh->has_vertex(neighbor->vertex(j))) break;
                    index[i] = edgesToDraw[index[i]].next_e(neighbor->vertex((j + 2) % 3)->info() );
                }
            }
        }
        if(settings.DEBUG_LEVEL > 1) {
            auto n0 = fh->vertex(0)->info();
            auto n1 = fh->vertex(1)->info();
            auto n2 = fh->vertex(2)->info();
            auto n0m = bounded_mercator(coords[n0]);
            auto n1m = bounded_mercator(coords[n1]);
            auto n2m = bounded_mercator(coords[n2]);
            assert(!CGAL::collinear(Point(n0m.lat, n0m.lon), Point(n1m.lat, n1m.lon), Point(n2m.lat, n2m.lon)));
            if(!CGAL::left_turn(Point(n0m.lat, n0m.lon), Point(n1m.lat, n1m.lon), Point(n2m.lat, n2m.lon))){
                printf("ERROR: Cell is not convex %f %f, %f %f, %f %f\n", n0m.lat, n0m.lon, n1m.lat, n1m.lon, n2m.lat, n2m.lon);
            }
        }

        assert(index[0] != UINT32_MAX && index[1] != UINT32_MAX && index[2] != UINT32_MAX);
        edgesToDraw[index[0]].setNeighbour(index[1], fh->vertex(2)->info());
        edgesToDraw[index[1]].setNeighbour(index[2], fh->vertex(0)->info());
        edgesToDraw[index[2]].setNeighbour(index[0], fh->vertex(1)->info());
        fh->info().edgeID = index[0];
    }
    cdt.clear();
    cdt = CDT(); // free memory
    std::cout << "found " << constrainedEdgeCount << " constrained edges and " << otherEdgeCount << " other edges" << std::endl;
    assert(constrainedEdgeCount + otherEdgeCount == totalEdgesInDomain);
    return edgesToDraw;
}

#endif