//
// Created by kochda on 13.12.23.
//
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <queue>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Real_timer.h>

#include <sys/time.h>

// HOLGER 22Jan06 : changed all suseconds_t to off_t 

//
//! A SIMPLE CLASS FOR TIME MEASUREMENTS.
//
class Timer
{
  private:

    //! The timer value (initially zero)
    off_t _usecs;

    //! The timer value at the last mark set (initially zero)
    off_t _usecs_at_mark;

    //! Used by the gettimeofday command.
    struct timeval _tstart;
    
    //! Used by the gettimeofday command.
    struct timeval _tend;
    
    //! Used by the gettimeofday command.
    struct timezone _tz;

    //! Indicates whether a measurement is running.
    bool _running;

  public:

    //! The default constructor.
    Timer() { reset(); }
    
    //! Resets the timer value to zero and stops the measurement.
    void reset() { _usecs = _usecs_at_mark = 0; _running = false; }

    //! Mark the current point in time (to be considered by next usecs_since_mark)
    void mark() { stop(); _usecs_at_mark = _usecs; cont(); }

    //! Resets the timer value to zero and starts the measurement.
    void start()
    {
      _usecs = _usecs_at_mark = 0;
      gettimeofday(&_tstart, &_tz);
      _running = true;
    }

    //! Continues the measurement without resetting the timer value (no effect it running)
    void cont()
    { 
      if (_running == false)
      {
	gettimeofday(&_tstart, &_tz);
	_running = true;
      }
    }
    
    //! Stops the measurement (does *not* return the timer value anymore)
    void stop()
    {
      gettimeofday(&_tend, &_tz);
      _running = false;
      _usecs += (off_t)(1000000) * (off_t)(_tend.tv_sec - _tstart.tv_sec) + (off_t)(_tend.tv_usec - _tstart.tv_usec);
        //return _usecs;
    }

  bool isRunning() const {
      return _running;
    }

    //! Time at last stop (initially zero)
    off_t value() const { return _usecs; } /* in microseconds */
    off_t usecs() const { return _usecs; } /* in microseconds */
    off_t msecs() const { return _usecs/1000; } /* in milliseconds */
    float msecs_float() const { return _usecs/1000.0; } 
    float secs() const { return _usecs/1000000.0; } /* in seconds */

    //! Time from last mark to last stop (initally zero)
    off_t usecs_since_mark() const { return _usecs - _usecs_at_mark; }


};

using namespace std;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_intersections_tag Itag;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CDT::Point Point;
typedef CDT::Edge Edge;
typedef K::Ray_2 Ray;
typedef CDT::Face_handle Face_handle;

using namespace std;

void WGS84toGoogleBing(double lat, double lon, double &x, double &y)
{
    x = lon * 20037508.34 / 180;
    y = log(tan((90 + lat) * M_PI / 360)) / (M_PI / 180);
    y = y * 20037508.34 / 180;
}
void GoogleBingtoWGS84Mercator(double x, double y, double &lat, double &lon)
{
    lon = (x / 20037508.34) * 180;
    lat = (y / 20037508.34) * 180;

    lat = 180 / M_PI * (2 * atan(exp(lat * M_PI / 180)) - M_PI / 2);
}

class LinkPair;

class GlobePoint
{
public:
    double lat, lon;
    int i;
    int poly_id;
    bool forceKeep = false;
    LinkPair *currPair;
    GlobePoint(double _lat, double _lon, int _i, int _poly_id = 0)
    {
        lat = _lat;
        lon = _lon;
        i = _i;
        poly_id = _poly_id;
    }
};

struct Coordinates
{
    double lat, lon;
    Coordinates() {}
    Coordinates(double lat_, double lon_) : lat(lat_), lon(lon_) {}
};

double dist(GlobePoint *p1, GlobePoint *p2, const bool &isMercator)
{
    double srcX1, srcY1, srcX2, srcY2;
    if (isMercator)
    {
        srcX1 = p1->lat;
        srcX2 = p2->lat;
        srcY1 = p1->lon;
        srcY2 = p2->lon;
    }
    else
    {
        WGS84toGoogleBing(p1->lat, p1->lon, srcX1, srcY1);
        WGS84toGoogleBing(p2->lat, p2->lon, srcX2, srcY2);
    }
    return sqrt((srcX1 - srcX2) * (srcX1 - srcX2) + (srcY1 - srcY2) * (srcY1 - srcY2));
}

class LinkedPoint
{
public:
    GlobePoint *point;
    LinkedPoint *prev, *next;
    LinkedPoint(GlobePoint *_point)
    {
        point = _point;
    }
};

class LinkPair
{
public:
    LinkedPoint *lpoint;
    double prio;
    LinkPair(LinkedPoint *_lpoint, const bool &isMercator)
    {
        lpoint = _lpoint;
        lpoint->point->currPair = this;
        GlobePoint *p1 = lpoint->prev->point;
        GlobePoint *gp = lpoint->point;
        GlobePoint *p3 = lpoint->next->point;
        double a = dist(p1, gp, isMercator);
        double b = dist(gp, p3, isMercator);
        double c = dist(p3, p1, isMercator);
        double s = 0.5 * (a + b + c);
        prio = s * (s - a) * (s - b) * (s - c); // sqrt not necessary
    }
};

struct compare
{
    bool operator()(const LinkPair *a, const LinkPair *b)
    {
        return a->prio > b->prio;
    }
};

void writeEntireGL(const vector<vector<GlobePoint>> &myPolysShrunk, const string &output)
{
    ofstream writer(output);
    writer << setprecision(20);
    u_int amountPoints = 0;
    for (auto polygon : myPolysShrunk)
    {
        amountPoints += polygon.size();
    }
    writer << amountPoints << "\n"
           << amountPoints << "\n";
    cout << "Writing points" << endl;
    for (auto polygon : myPolysShrunk)
    {
        for (auto point : polygon)
        {
            double lat, lon;
            GoogleBingtoWGS84Mercator(point.lat, point.lon, lat, lon);
            writer << lat << " " << lon << "\n";
        }
    }
    cout << "Writing edges" << endl;
    u_int ptsSoFar = 0;
    for (auto polygon : myPolysShrunk)
    {
        for (int i = 0; i < polygon.size(); ++i)
        {
            int add = (i + 1) % (polygon.size());
            int writeAmount = ptsSoFar + add;
            // cout << add << " " << writeAmount << endl;
            writer << ptsSoFar + i << " " << writeAmount << " 4 2 " << "\n";
        }
        ptsSoFar += polygon.size();
    }
    cout << "Done writing gl" << endl;
    writer.close();
}

int findLen(LinkedPoint *lp)
{
    int counter = 1;
    LinkedPoint *curr = lp->next;
    while (curr != lp)
    {
        curr = curr->next;
        counter++;
    }
    return counter;
}

void writeGL(LinkedPoint *lp, string path, int color)
{
    int len = findLen(lp);
    ofstream outFile(path);
    outFile << setprecision(20);
    outFile << len << endl;
    outFile << len << endl;

    // points
    LinkedPoint *curr = lp;
    for (int i = 0; i < len; i++)
    {
        double lat, lon;
        GoogleBingtoWGS84Mercator(curr->point->lat, curr->point->lon, lat, lon);
        outFile << setprecision(20) << lat << " " << lon << "\n";
        curr = curr->next;
    }
    // edges
    outFile << len - 1 << " 0 4 " << color << endl;
    for (int i = 1; i < len; i++)
    {
        outFile << i - 1 << " " << i << " 4 " << color << endl;
    }
}

void printTriangle(const CDT::Vertex_handle &v1, const CDT::Vertex_handle &v2, const CDT::Vertex_handle &v3, int offset)
{
    double x1 = v1->point().x();
    double y1 = v1->point().y();
    double x2 = v2->point().x();
    double y2 = v2->point().y();
    double x3 = v3->point().x();
    double y3 = v3->point().y();
    double lat1, lon1, lat2, lon2, lat3, lon3;
    GoogleBingtoWGS84Mercator(x1, y1, lat1, lon1);
    GoogleBingtoWGS84Mercator(x2, y2, lat2, lon2);
    GoogleBingtoWGS84Mercator(x3, y3, lat3, lon3);
    cout << setprecision(20) << lat1 << " " << lon1 << endl;
    cout << setprecision(20) << lat2 << " " << lon2 << endl;
    cout << setprecision(20) << lat3 << " " << lon3 << endl;
    cout << offset << " " << offset + 1 << " 4 2" << endl;
    cout << offset + 1 << " " << offset + 2 << " 4 2" << endl;
    cout << offset + 2 << " " << offset << " 4 2" << endl;
}

bool pointsSeeEachOther(const CDT::Vertex_handle &startVertex, const CDT::Vertex_handle &targetVertex, const CDT::Vertex_handle &skipVertex, const CDT &cdt)
{
    //	cout << startVertex->info() << " to " << targetVertex->info() << " via " << skipVertex->info() << endl;
    //	printTriangle(startVertex, targetVertex, skipVertex, 0);

    //	cout << "start/skip/target" << endl;
    //	cout << startVertex->point() << endl;
    //	cout << skipVertex->point() << endl;
    //	cout << targetVertex->point() << endl;

    if (CGAL::orientation(startVertex->point(), skipVertex->point(), targetVertex->point()) == CGAL::RIGHT_TURN)
    {
        //		cout << "concave" << endl;
        return false; // concavity
    }
    if (CGAL::orientation(startVertex->point(), skipVertex->point(), targetVertex->point()) == CGAL::COLLINEAR)
    {
        //		cout << "collinear" << endl;
        return true;
    }

    auto edgeCirculator = startVertex->incident_edges();

    // if vertex is adjacent to infinite vertex, start from there so we don't see it again
    if (cdt.is_edge(startVertex, cdt.infinite_vertex()))
    {
        //		cout << "sees inf" << endl;
        while (edgeCirculator->first->vertex((edgeCirculator->second + 1) % 3) != cdt.infinite_vertex())
        {
            edgeCirculator++;
        }
        edgeCirculator++;
        //		cout << "found inf" << endl;
    }

    // rotate until we find the two edges between which skipVertex lies
    while (
        !(
            (CGAL::orientation(startVertex->point(), skipVertex->point(), edgeCirculator->first->vertex(edgeCirculator->second)->point()) != CGAL::LEFT_TURN) &&
            (CGAL::orientation(startVertex->point(), skipVertex->point(), edgeCirculator->first->vertex((edgeCirculator->second + 1) % 3)->point()) == CGAL::LEFT_TURN)))
    {
        edgeCirculator++;
    }
    // found triangle at startVertex through which the ray from startVertex to skipVertex goes
    // (if skipVertex is a vertex of this triangle, it is the right vertex as seen from startVertex)
    //	printTriangle(edgeCirculator->first->vertex(0), edgeCirculator->first->vertex(1), edgeCirculator->first->vertex(2), 0);
    //	cout << (CGAL::orientation(startVertex->point(), skipVertex->point(), edgeCirculator->first->vertex(edgeCirculator->second)->point()) != CGAL::LEFT_TURN) << endl;
    //	cout << (CGAL::orientation(startVertex->point(), skipVertex->point(), edgeCirculator->first->vertex( (edgeCirculator->second+1)%3 )->point()) == CGAL::LEFT_TURN) << endl;
    //	printTriangle(startVertex, skipVertex, edgeCirculator->first->vertex(edgeCirculator->second), 6);

    // find first edge ccw of targetVertex (or edge to targetVertex) (but don't go across constraint edges)
    while (CGAL::orientation(startVertex->point(), targetVertex->point(), edgeCirculator->first->vertex((edgeCirculator->second + 1) % 3)->point()) == CGAL::RIGHT_TURN)
    {
        //		cout << "rotate" << endl;
        if (cdt.is_constrained(*edgeCirculator))
        { // we would walk across a constraint edge now
            return false;
        }
        edgeCirculator++;
    }
    // found triangle at startVertex through which the ray from startVertex to targetVertex goes

    // walk along ray until either finding a constraint edge, or finding targetVertex
    //	printTriangle(startVertex, skipVertex, targetVertex, 0);
    //	cout << "start walk" << endl;
    CDT::Face_handle currFace = edgeCirculator->first;
    CDT::Vertex_handle currVertex = edgeCirculator->first->vertex((edgeCirculator->second + 1) % 3);
    while (currVertex != targetVertex)
    {
        //		printTriangle(currFace->vertex(0), currFace->vertex(1), currFace->vertex(2), 3);
        //		cout << "walking" << endl;
        // go left or right of currVertex?
        int nextFaceIndex = -1;
        if (CGAL::orientation(startVertex->point(), targetVertex->point(), currVertex->point()) == CGAL::RIGHT_TURN)
        {
            nextFaceIndex = (currFace->index(currVertex) + 2) % 3; // go left
                                                                   //			cout << "go left" << endl;
        }
        else
        {
            nextFaceIndex = (currFace->index(currVertex) + 1) % 3; // go right
                                                                   //			cout << "go right" << endl;
        }
        if (currFace->is_constrained(nextFaceIndex))
        { // constraint edge blocking the way
            //			cout << "constraint edge in the way" << endl;
            return false;
        }
        CDT::Face_handle oldFace = currFace;
        currFace = currFace->neighbor(nextFaceIndex);
        currVertex = currFace->vertex(currFace->index(oldFace));
    }
    //	cout << "visible" << endl;
    return true;
}

// returns one of the remaining points
LinkedPoint *shrinkPolygon(vector<GlobePoint> &poly, const std::vector<CDT::Vertex_handle> &myVHs, const CDT &cdt, const bool &isMercator)
{
    LinkedPoint *first = new LinkedPoint(&poly.at(0));
    LinkedPoint *curr = first;
    for (int j = 1; j < poly.size(); j++)
    {
        LinkedPoint *newPoint = new LinkedPoint(&poly.at(j));
        curr->next = newPoint;
        newPoint->prev = curr;
        curr = newPoint;
    }
    curr->next = first;
    first->prev = curr;

    //	writeGL(first, std::string("coastline-big.gl"), 4);

    int removeTarget = 0;
    if (poly.size() > 2000)
    {
        removeTarget = poly.size() * 0.995;
    }
    else
    {
        removeTarget = poly.size() * 0.95;
    }
    if (poly.size() - removeTarget < 3)
    {
        removeTarget = poly.size() - 3;
    }
    priority_queue<LinkPair *, vector<LinkPair *>, compare> pq;
    curr = first;
    for (int j = 0; j < poly.size(); j++)
    {
        if (pointsSeeEachOther(myVHs[curr->prev->point->i], myVHs[curr->next->point->i], myVHs[curr->point->i], cdt))
        {
            LinkPair *pair = new LinkPair(curr, isMercator);
            pq.push(pair);
        }
        curr = curr->next;
    }

    int removed = 0;
    LinkedPoint *somePoint = curr; // somePoint is always part of the remaining polygon
    while (pq.size() > 0 && removed < removeTarget)
    {
        LinkPair *pair = pq.top();
        pq.pop();
        GlobePoint *point = pair->lpoint->point;
        if (!point->forceKeep && pair == point->currPair)
        { // ignore duplicate insertions
            // remove from coastline and update neighbors
            LinkedPoint *prev = pair->lpoint->prev;
            LinkedPoint *next = pair->lpoint->next;
            prev->next = next;
            next->prev = prev;
            somePoint = next;
            removed++;

            // create new pairs so old ones are invalidated
            LinkPair *pair1 = new LinkPair(prev, isMercator);
            LinkPair *pair2 = new LinkPair(next, isMercator);
            // pair1->prio += pair->prio;
            // pair2->prio += pair->prio;
            // reinsert neighbors into pq if necessary
            if (pointsSeeEachOther(myVHs[prev->prev->point->i], myVHs[prev->next->point->i], myVHs[prev->point->i], cdt))
            {
                pq.push(pair1);
            }
            if (pointsSeeEachOther(myVHs[next->prev->point->i], myVHs[next->next->point->i], myVHs[next->point->i], cdt))
            {
                pq.push(pair2);
            }
        }
    }
    if (removed < removeTarget)
    {
        cout << "pq ran empty before removal target, removed only " << removed << " instead of " << removeTarget << endl;
        //		writeGL(somePoint, std::string("coastline-small.gl"), 2);
        //		exit(0);
    }

    return somePoint;
}

// returns one of the remaining points
vector<LinkedPoint *> shrinkAllPolygonsAtOnce(long reduced_size, vector<vector<GlobePoint>> &polys, const vector<std::vector<CDT::Vertex_handle>> &vhs_per_poly, const vector<CDT> &cdts_per_poly, const bool &isMercator)
{
    // DEBUG
    cout << "Start shrinking" << endl;
    Timer timer;
    timer.start();
    // DBEUG END
    vector<LinkedPoint *> first_points_of_polys(polys.size());
    long current_size = 0;
    vector<int> size_per_poly(polys.size(), 0);
    for (int i = 0; i < polys.size(); i++)
    {
        LinkedPoint *first = new LinkedPoint(&polys.at(i).at(0));
        LinkedPoint *curr = first;
        for (int j = 1; j < polys.at(i).size(); j++)
        {
            LinkedPoint *newPoint = new LinkedPoint(&polys.at(i).at(j));
            curr->next = newPoint;
            newPoint->prev = curr;
            curr = newPoint;
        }
        curr->next = first;
        first->prev = curr;
        first_points_of_polys.at(i) = first;
        if (polys.at(i).size() > 2)
        {
            current_size += polys.at(i).size();
            size_per_poly.at(i) = polys.at(i).size();
        }
    }
    timer.stop();
    cout << "Setting up linked points took " << timer.secs() << " seconds" << endl;
    timer.reset();
    timer.start();
    priority_queue<LinkPair *, vector<LinkPair *>, compare> pq;
    for (int i = 0; i < polys.size(); i++)
    {
        LinkedPoint *curr = first_points_of_polys.at(i);
        for (int j = 0; j < polys.at(i).size() && size_per_poly.at(i) >= 3; j++)
        {
            if (pointsSeeEachOther(vhs_per_poly.at(curr->point->poly_id).at(curr->prev->point->i), vhs_per_poly.at(curr->point->poly_id).at(curr->next->point->i), vhs_per_poly.at(curr->point->poly_id).at(curr->point->i), cdts_per_poly.at(curr->point->poly_id)))
            {
                LinkPair *pair = new LinkPair(curr, isMercator);
                pq.push(pair);
            }
            curr = curr->next;
        }
    }
    timer.stop();
    cout << "Setting up the queue took " << timer.secs() << " seconds" << endl;
    while (pq.size() > 0 && current_size > reduced_size)
    {
        LinkPair *pair = pq.top();
        pq.pop();
        GlobePoint *point = pair->lpoint->point;
        if (size_per_poly.at(point->poly_id) < 3)
        {
            continue;
        }
        if (!point->forceKeep && pair == point->currPair)
        { // ignore duplicate insertions
            // remove from coastline and update neighbors
            LinkedPoint *prev = pair->lpoint->prev;
            LinkedPoint *next = pair->lpoint->next;
            prev->next = next;
            next->prev = prev;
            first_points_of_polys.at(next->point->poly_id) = next;
            current_size--;
            size_per_poly.at(next->point->poly_id)--;
            if (size_per_poly.at(next->point->poly_id) < 3)
            {
                current_size -= size_per_poly.at(next->point->poly_id);
                continue;
            }
            // create new pairs so old ones are invalidated
            LinkPair *pair1 = new LinkPair(prev, isMercator);
            LinkPair *pair2 = new LinkPair(next, isMercator);
            // pair1->prio += pair->prio;
            // pair2->prio += pair->prio;
            // reinsert neighbors into pq if necessary
            if (pointsSeeEachOther(vhs_per_poly.at(prev->point->poly_id)[prev->prev->point->i], vhs_per_poly.at(prev->point->poly_id)[prev->next->point->i], vhs_per_poly.at(prev->point->poly_id)[prev->point->i], cdts_per_poly.at(prev->point->poly_id)))
            {
                pq.push(pair1);
            }
            if (pointsSeeEachOther(vhs_per_poly.at(prev->point->poly_id)[next->prev->point->i], vhs_per_poly.at(prev->point->poly_id)[next->next->point->i], vhs_per_poly.at(prev->point->poly_id)[next->point->i], cdts_per_poly.at(prev->point->poly_id)))
            {
                pq.push(pair2);
            }
        }
    }
    if (current_size > reduced_size)
    {
        cout << "pq ran empty before removal target, reduced size is " << current_size << " instead of " << reduced_size << endl;
    }

    cout << "Finished shrinking" << endl;

    // DEBUG
    {
        long zero_counter = 0;
        long max_shrinks = 0;
        long max_shrinks_poly = 0;
        for (int i = 0; i < polys.size(); i++)
        {
            long num_shrinks = polys.at(i).size() - size_per_poly.at(i);
            if (num_shrinks == 0)
            {
                zero_counter++;
            }
            if (num_shrinks > max_shrinks)
            {
                max_shrinks = num_shrinks;
                max_shrinks_poly = i;
            }
        }
        cout << "Number of polys with no removed nodes: " << zero_counter << endl;
        cout << "Max removed nodes per poly: " << max_shrinks << " out of " << polys.at(max_shrinks_poly).size() << endl;
    }
    // DEBUG END

    return first_points_of_polys;
}

// if isMercator, file is expected to be mercator and the output will be mercator as well
void landshrinker(vector<vector<GlobePoint>> &myPolys, long target_number_of_nodes, const string &output, const bool &isMercator)
{
    // for each polygon: compute CGAL CDT, shrink polygon, and write to output file
    Timer localTimer;
    localTimer.start();
    int outputCount = 0;
    int smallCount = 0;
    int remainingPolygonCount = 0;
    vector<vector<GlobePoint>> myPolysShrunk;
    long num_corners_all_polys = 0;
    vector<CDT> cdts_per_poly(myPolys.size());
    vector<std::vector<CDT::Vertex_handle>> vhs_per_poly(myPolys.size());
    for (int i = 0; i < myPolys.size(); i++)
    {
        // cout << "coastline " << i << ": " << myPolys.at(i).size() << endl;

        // CGAL triangulation

        std::vector<Point> myPoints;

        int counter = 0;
        for (int j = 0; j < myPolys.at(i).size(); j++)
        {
            double srcX, srcY;
            if (isMercator)
            {
                srcX = myPolys.at(i).at(j).lat;
                srcY = myPolys.at(i).at(j).lon;
            }
            else
            {
                double curLat = myPolys.at(i).at(j).lat;
                double curLon = myPolys.at(i).at(j).lon;
                WGS84toGoogleBing(myPolys.at(i).at(j).lat, myPolys.at(i).at(j).lon, srcX, srcY);
            }
            //			cout << "[" << setprecision(20) << curLon << ", " << curLat << "]," << endl;
            //			cout << setprecision(20) << srcX << " " << srcY << endl;
            myPoints.push_back(Point(srcX, srcY));
            vhs_per_poly.at(i).push_back(cdts_per_poly.at(i).insert(myPoints.back()));
            vhs_per_poly.at(i).back()->info() = counter++;
        }
        for (int j = 0; j < myPolys.at(i).size(); j++)
        {
            cdts_per_poly.at(i).insert_constraint(vhs_per_poly.at(i).at(j), vhs_per_poly.at(i).at((j + 1) % myPolys.at(i).size()));
        }
        assert(cdts_per_poly.at(i).is_valid());
        num_corners_all_polys += myPolys.at(i).size();
    }
    localTimer.stop();
    cout << "File: " << output << " . Setting up triangulation in " << localTimer.secs() << " s." << "\n";
    localTimer.cont();
    //  shrink polygon
    vector<LinkedPoint *> remaining_point_per_poly = shrinkAllPolygonsAtOnce(target_number_of_nodes, myPolys, vhs_per_poly, cdts_per_poly, isMercator);
    for (int p = 0; p < remaining_point_per_poly.size(); p++)
    {
        LinkedPoint *remainingPoint = remaining_point_per_poly.at(p);
        int remainingSize = findLen(remainingPoint);
        if (remainingSize > 2)
        {
            vector<GlobePoint> curPolyShrunk;
            LinkedPoint *curr = remainingPoint;
            for (int j = 0; j < remainingSize; j++)
            {
                curPolyShrunk.push_back(GlobePoint(curr->point->lat, curr->point->lon, 0));
                curr = curr->next;
            }
            myPolysShrunk.push_back(curPolyShrunk);
            outputCount += remainingSize;
            if (myPolys.at(remainingPoint->point->poly_id).size() < 1000)
            {
                smallCount += remainingSize;
            }
            remainingPolygonCount++;
        }
    }
    localTimer.stop();
    cout << "File: " << output << " . Shrunk in " << localTimer.secs() << " s." << endl;
    //writeEntireGL(myPolysShrunk, "shrunk.gl");

    ofstream outFile(output);
    outFile << setprecision(20);
    outFile << remainingPolygonCount << endl;
    for (int i = 0; i < remainingPolygonCount; i++)
    {
        outFile << myPolysShrunk.at(i).size() << endl;
        for (int j = 0; j < myPolysShrunk.at(i).size(); j++)
        {
            outFile << myPolysShrunk.at(i).at(j).lat << " " << myPolysShrunk.at(i).at(j).lon << endl;
        }
    }
    outFile.close();
    cout << "shrunk to " << outputCount << " points in " << remainingPolygonCount << " coastlines and wrote to coastlinesShrunk.txt" << endl;
    cout << "s: " << smallCount << endl;
}

vector<vector<GlobePoint>> readInPolygons(string input, const bool keepCornerPolygons)
{
    vector<vector<GlobePoint>> myPolys;
    ifstream inFile(input);
    int nofPolys;
    inFile >> nofPolys;
    cout << "We have " << nofPolys << " polygons to read" << endl;

    CGAL::Real_timer myTimer;
    myTimer.start();

    int inputSize = 0;
    for (int i = 0; i < nofPolys; i++)
    {
        if (i % 100000 == 0)
            cout << i << " " << flush;
        int polySize;
        inFile >> polySize;
        if (polySize == 0)
            cout << i << ": " << polySize << " ";
        vector<GlobePoint> curPoly;
        for (int j = 0; j < polySize; j++)
        {
            double lat, lon;
            // store x y as lat lon in case we get mercator input
            inFile >> lat >> lon;
            if (j > 0 && lat == curPoly.back().lat && lon == curPoly.back().lon)
            {
                cout << "identical points: " << i << " " << j << " " << inputSize << endl;
            }
            curPoly.push_back(GlobePoint(lat, lon, j, i));
        }
        /*
        if (i == 624383 || i == 626210) { // hack fix for polygons with wrong orientation in coastlines3.txt
            reverse(curPoly.begin(), curPoly.end());
            for (int j = 0; j < curPoly.size(); j++) {
                curPoly.at(j).i = j;
            }
        }
        */
        inputSize += curPoly.size();
        myPolys.push_back(curPoly);
    }

    inFile.close();
    // The last four polygons are the corner polygons and are important to keep. So we forbid to delete these
    if (keepCornerPolygons) {
        for(int i = myPolys.size() - 1; i >= myPolys.size() - 4; i--) {
            for(int j = 0; j < myPolys[i].size(); j++) {
                myPolys[i][j].forceKeep = true;
            }
        }
    }

    myTimer.stop();
    cout << "Reading data in " << myTimer.time() << "s\n";
    cout << "read " << inputSize << " corners and " << myPolys.size() << " polygons" << endl;
    return myPolys;
}

vector<vector<Coordinates>> readInPaths(string path_to_path_file)
{
    ifstream reader(path_to_path_file);
    vector<vector<Coordinates>> paths;
    string tmp;
    reader >> tmp;
    if (tmp.find("Query") != string::npos)
    {
        while(reader >> tmp)
        {
            while(tmp.find("Distance") == string::npos)
            {
                reader >> tmp;
            }
            int path_size;
            reader >> tmp >> path_size;
            vector<Coordinates> path;
            for(int i = 0; i < path_size; i++)
            {
                Coordinates coords;
                reader >> coords.lat >> coords.lon;
                path.push_back(coords);
            }
            paths.push_back(path);
        }
        cout << "Finished reading " << paths.size() << " paths" << endl;
    }
    else
    {
        int num_paths = stoi(tmp);
        paths.assign(num_paths, vector<Coordinates>());
        cout << "Start reading " << num_paths << " paths" << endl;
        for (int i = 0; i < num_paths; i++)
        {
            reader >> tmp >> tmp >> tmp;
            int path_size;
            reader >> path_size;
            paths.at(i).reserve(path_size);
            for (int j = 0; j < path_size; j++)
            {
                double lat, lon;
                reader >> tmp >> lat >> lon;
                paths.at(i).push_back(Coordinates(lat, lon));
            }
        }
    }

    reader.close();
    return paths;
}

double getDistanceBetweenPointAndCoordinates(const GlobePoint &point, const Coordinates &coord)
{
    double diff_lat = point.lat - coord.lat;
    double diff_lon = point.lon - coord.lon;
    return diff_lat * diff_lat + diff_lon * diff_lon;
}

double getDistanceBetweenPoints(const GlobePoint &point1, const GlobePoint &point2)
{
    double diff_lat = point1.lat - point2.lat;
    double diff_lon = point1.lon - point2.lon;
    return diff_lat * diff_lat + diff_lon * diff_lon;
}

vector<vector<int>> getClosestPoints(const vector<vector<Coordinates>> &coords_per_path, const vector<GlobePoint> &points)
{
    const int NO_ENTRY = -1;
    int max_cluster_size = 100;
    vector<vector<int>> closest_point_ids(coords_per_path.size());
    struct TreeNode
    {
        int point_id;
        int child1, child2;
        double max_distance;
        vector<int> cluster;
        TreeNode()
        {
            max_distance = 0;
            child1 = NO_ENTRY, child2 = NO_ENTRY;
        }
        TreeNode(int point_id_) : point_id(point_id_)
        {
            max_distance = 0;
            child1 = NO_ENTRY, child2 = NO_ENTRY;
        }
    };
    TreeNode start_node;
    {
        for (int i = 0; i < points.size(); i++)
        {
            start_node.cluster.push_back(i);
        }
    }
    vector<TreeNode> tree_nodes = {start_node};
    int index = 0;
    long debug_dist_counter = 0;
    long debug_node_counter = 0;
    while (index < tree_nodes.size())
    {
        if (tree_nodes.at(index).cluster.size() < max_cluster_size)
        {
            index++;
            continue;
        }
        int pointid_child1 = tree_nodes.at(index).cluster.at(rand() % tree_nodes.at(index).cluster.size());
        int pointid_child2 = tree_nodes.at(index).cluster.at(rand() % tree_nodes.at(index).cluster.size());
        TreeNode child1 = TreeNode(pointid_child1);
        TreeNode child2 = TreeNode(pointid_child2);
        for (int node : tree_nodes.at(index).cluster)
        {
            double dist1 = getDistanceBetweenPoints(points.at(node), points.at(pointid_child1));
            double dist2 = getDistanceBetweenPoints(points.at(node), points.at(pointid_child2));
            debug_dist_counter += 2;
            if (dist1 < dist2)
            {
                child1.cluster.push_back(node);
                if (child1.max_distance < dist1)
                {
                    child1.max_distance = dist1;
                }
            }
            else
            {
                child2.cluster.push_back(node);
                if (child2.max_distance < dist2)
                {
                    child2.max_distance = dist2;
                }
            }
        }
        tree_nodes.at(index).cluster.clear();
        tree_nodes.at(index).child1 = tree_nodes.size();
        tree_nodes.push_back(child1);
        tree_nodes.at(index).child2 = tree_nodes.size();
        tree_nodes.push_back(child2);
        index++;
    }
    // DEBUG
    // cout << "Tree construciton finished" << endl;
    // DEBUG END
    for (int i = 0; i < coords_per_path.size(); i++)
    {
        // DEBUG
        // cout << "Path " << i << endl;
        // DEBUG END
        debug_node_counter += coords_per_path.at(i).size();
        for (const auto &coord : coords_per_path.at(i))
        {
            typedef pair<double, int> PQ_Pair;
            priority_queue<PQ_Pair, vector<PQ_Pair>, greater<PQ_Pair>> pq;
            pq.push(make_pair(0, 0));
            double min_dist = INT_MAX;
            int closest_node = NO_ENTRY;
            while (!pq.empty())
            {
                const TreeNode &next_node = tree_nodes.at(pq.top().second);
                if (closest_node != NO_ENTRY && pq.top().first - next_node.max_distance > min_dist)
                {
                    pq.pop();
                    continue;
                }
                pq.pop();
                for (int node : next_node.cluster)
                {
                    double dist = getDistanceBetweenPointAndCoordinates(points.at(node), coord);
                    if (closest_node == NO_ENTRY || dist < min_dist)
                    {
                        min_dist = dist;
                        closest_node = node;
                    }
                }
                debug_dist_counter += next_node.cluster.size();
                if (next_node.child1 != NO_ENTRY)
                {
                    double dist1 = getDistanceBetweenPointAndCoordinates(points.at(tree_nodes.at(next_node.child1).point_id), coord);
                    if (closest_node == NO_ENTRY || dist1 < min_dist)
                    {
                        min_dist = dist1;
                        closest_node = tree_nodes.at(next_node.child1).point_id;
                    }
                    if (closest_node == NO_ENTRY || dist1 - tree_nodes.at(next_node.child1).max_distance < min_dist)
                    {
                        pq.push(make_pair(dist1, next_node.child1));
                    }
                    double dist2 = getDistanceBetweenPointAndCoordinates(points.at(tree_nodes.at(next_node.child2).point_id), coord);
                    if (closest_node == NO_ENTRY || dist2 < min_dist)
                    {
                        min_dist = dist2;
                        closest_node = tree_nodes.at(next_node.child2).point_id;
                    }

                    if (closest_node == NO_ENTRY || dist2 - tree_nodes.at(next_node.child2).max_distance < min_dist)
                    {
                        pq.push(make_pair(dist2, next_node.child2));
                    }
                    debug_dist_counter += 2;
                }
            }
            closest_point_ids.at(i).push_back(closest_node);
        }
    }
    // DEBUG
    if (debug_node_counter > 0)
    {
        cout << "Conducted " << debug_dist_counter / debug_node_counter << " distance computations on average" << endl;
    }
    // DEBUG END
    return closest_point_ids;
}

int main(int argc, char *argv[])
{
    string coastlines, output, path_to_path_file;
    long num_remaining_nodes = -1;
    double fraction_of_remaining_nodes = -1;
    double fraction_of_poi_nodes = 0.5; //If path file is provided, important nodes are determined. The value of this variable determines the fraction of important nodes among the remaining nodes.
    vector<vector<Coordinates>> paths;
    bool keepCornerPolygons = false;
    while (1)
    {
        int c = getopt(argc, argv, "c:o:n:f:p:i:");
        if (c == -1)
            break; /* end of list */
        switch (c)
        {
        case '?':
            break;
        case ':':
            fprintf(stderr, "missing argument.\n");
            break;
        case 'c':
            coastlines = optarg;
            break;
        case 'o':
            output = optarg;
            break;
        case 'f':
            fraction_of_remaining_nodes = stod(optarg);
            break;
        case 'n':
            num_remaining_nodes = stoi(optarg);
            break;
        case 'p':
            path_to_path_file = optarg;
            paths = readInPaths(path_to_path_file);
            break;
        case 'i':
            fraction_of_poi_nodes = stod(optarg);
            break;
        case 'k':
            keepCornerPolygons = true;
            break;
        }
    }

    Timer timer;
    timer.start();
    vector<vector<GlobePoint>> polygons = readInPolygons(coastlines, keepCornerPolygons);

    long target_num_remaining_nodes = 0;
    if (num_remaining_nodes >= 0)
    {
        target_num_remaining_nodes = num_remaining_nodes;
    }
    else
    {
        long num_nodes = 0;
        for (int i = 0; i < polygons.size(); i++)
        {
            if (polygons.at(i).size() >= 3)
            {
                num_nodes += polygons.at(i).size();
            }
        }
        if (fraction_of_remaining_nodes >= 0)
        {
            target_num_remaining_nodes = fraction_of_remaining_nodes * num_nodes;
        }
        else
        {
            cout << "Warning: target size unspecified. Use -n to set an absolute value or -f to set a relative value." << endl;
            target_num_remaining_nodes = 0.01 * num_nodes;
        }
    }
    if (paths.size() > 0)
    {
        int num_important_nodes = target_num_remaining_nodes * fraction_of_poi_nodes;
        cout << "Use path file to keep " << num_important_nodes << " important nodes" << endl;
        vector<GlobePoint> all_points;
        vector<vector<int>> local_to_global_index(polygons.size());
        vector<int> sorted_indices;
        int index = 0;
        for (int i = 0; i < polygons.size(); i++)
        {
            for (int j = 0; j < polygons.at(i).size(); j++)
            {
                all_points.push_back(polygons.at(i).at(j));
                sorted_indices.push_back(index);
                local_to_global_index.at(i).push_back(index++);
            }
        }
        vector<vector<int>> closest_points_per_path = getClosestPoints(paths, all_points);
        vector<int> counter_per_point(all_points.size(), 0);
        long num_nodes_in_paths = 0;
        for (int i = 0; i < closest_points_per_path.size(); i++)
        {
            num_nodes_in_paths += closest_points_per_path.at(i).size();
            for (int j = 0; j < closest_points_per_path.at(i).size(); j++)
            {
                counter_per_point.at(closest_points_per_path.at(i).at(j))++;
            }
        }
        sort(sorted_indices.begin(), sorted_indices.end(), [&](int i, int j)
             { return counter_per_point.at(i) > counter_per_point.at(j); });
        vector<bool> is_important_node(all_points.size(), false);
        for (int i = 0; i < num_important_nodes; i++)
        {
            is_important_node.at(sorted_indices.at(i)) = true;
        }
        for (int i = 0; i < polygons.size(); i++)
        {
            for (int j = 0; j < polygons.at(i).size(); j++)
            {
                if (is_important_node.at(local_to_global_index.at(i).at(j)))
                {
                    polygons.at(i).at(j).forceKeep = true;
                }
            }
        }
    }

    cout << "Try to shrink polygons to " << target_num_remaining_nodes << " nodes" << endl;
    landshrinker(polygons, target_num_remaining_nodes, output, true);
    timer.stop();
    cout << "TIME to shrink: " << coastlines << " is "<< timer.secs() << "\n"
         << endl;
}