#ifndef MY_TYPES_H
#define MY_TYPES_H

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <limits>
#include <queue>
#include <ranges>
#include <ostream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <type_traits>
#include <functional>
#include <stdint.h>
#include <cstdlib>
#include <cmath>

#include "cgalTypes.hpp"

struct Edge_t;
struct Face_t;
struct HalfEdge_t;
class Triangulation;

struct Pair { uint32_t a, b; };
struct Triple { uint32_t a, b, c; };
enum Direction {LEFT, RIGHT, MID};

char ccw(char index){return (index + 1) % 3;}
char cw(char index){return (index + 2) % 3;}

struct Wavelet{
    public:
    uint32_t left_b = UINT32_MAX, right_b = UINT32_MAX, generator = UINT32_MAX;
    NT_Distance cost_generator = 0.0;
    uint32_t prev_generator = UINT32_MAX;
    uint32_t item, id;
    Wavelet(){};
    Wavelet(uint32_t generator, uint32_t node) : generator(generator), item(node) {};
    Wavelet(uint32_t id, uint32_t generator, uint32_t edge, uint32_t left_b, uint32_t right_b, NT_Distance & cost_generator, uint32_t prev_generator = UINT32_MAX) : 
            id(id),
            generator(generator), 
            item(edge), 
            left_b(left_b), 
            right_b(right_b), 
            cost_generator(cost_generator),
            prev_generator(prev_generator){};
    bool isEdge(){
        return left_b == UINT32_MAX && right_b == UINT32_MAX;
    };
};

struct Edge {
    Edge() {}
    Edge(uint32_t source, uint32_t target) : source(source), target(target) {}
    uint32_t source, target;
};
std::ostream& operator <<(std::ostream& stream, const Edge& e) {
    return stream << e.source << " " << e.target;
};
std::istream& operator >>(std::istream& stream, Edge& e) {
    return stream >> e.source >> e.target;
};

uint8_t next(uint8_t i){return (i+1)%2;}
struct Edge_t {
    uint32_t E[2] = {UINT32_MAX, UINT32_MAX};
    uint32_t N[2] = {UINT32_MAX, UINT32_MAX};
    Edge_t() {}
    Edge_t(uint32_t n0, uint32_t n1) : N{ n0, n1} {}
    Edge_t(uint32_t n0, uint32_t n1, uint32_t e0, uint32_t e1) : E{ e0, e1}, N{ n0, n1} {}
    
    void flip(){ std::swap(N[0], N[1]); std::swap(E[0], E[1]); };
    bool valid() const {return !hasNode(UINT32_MAX);};
    void setNeighbour(uint32_t e, uint32_t n) {
        E[getIndex_n(n)] = e;
    }
    bool hasNode(const uint32_t n) const { return N[0] == n || N[1] == n; }
    bool isinfinte() const { return hasNode(UINT32_MAX); }
    char getIndex_n(uint32_t n) const {
        assert(hasNode(n));
        return N[0] == n ? 0 : 1;
    }
    uint32_t other_n(uint32_t n) const { return N[next(getIndex_n(n))]; }
    bool hasNeigbour(const uint32_t e) const { return E[0] == e || E[1] == e; }
    bool isConstrained() const { return hasNeigbour(UINT32_MAX); }
    char getIndex_e(uint32_t e) const {
        assert(hasNeigbour(e));
        return E[0] == e ? 0 : 1;
    }
    uint32_t next_e(uint32_t n) const { return E[getIndex_n(n)]; }
    void set_next_e(uint32_t n, uint32_t e) { E[getIndex_n(n)] = e; }
};
std::ostream& operator <<(std::ostream& stream, const Edge_t& e) {
    return stream << e.N[0] << ' ' << e.N[1] << ' ' << e.E[0] << ' ' << e.E[1];
};
std::istream& operator >>(std::istream& stream, Edge_t& e) {
    return stream >> e.N[0] >> e.N[1] >> e.E[0] >> e.E[1];
};

struct Face_t {
    uint32_t F[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};
    uint32_t N[3] = {UINT32_MAX, UINT32_MAX, UINT32_MAX};

    Face_t(uint32_t n0, uint32_t n1, uint32_t n2) : N{ n0, n1 , n2 } {}
    Face_t(uint32_t n0, uint32_t n1, uint32_t n2, uint32_t f0, uint32_t f1, uint32_t f2) : F{ f0, f1, f2 }, N{ n0, n1 , n2 } {}

    bool hasNode(const uint32_t n) { return N[0] == n || N[1] == n || N[2] == n; }
    bool isinfinte() { return hasNode(UINT32_MAX); }
    char getIndex_n(uint32_t n) {
        assert(hasNode(n));
        return N[0] == n ? 0 : (N[1] == n ? 1 : 2);
    }
    uint32_t next_cw_n(uint32_t n) { return N[cw(getIndex_n(n))]; }
    uint32_t next_ccw_n(uint32_t n) { return N[ccw(getIndex_n(n))]; }

    bool hasNeigbour(const uint32_t f) { return F[0] == f || F[1] == f || F[2] == f; }
    bool isWellDefined() { return hasNeigbour(UINT32_MAX); }
    char getIndex_f(uint32_t f) {
        assert(hasNeigbour(f));
        return F[0] == f ? 0 : (F[1] == f ? 1 : 2);
    }
    uint32_t next_cw_f(uint32_t n) { return F[cw(getIndex_n(n))]; }
    uint32_t next_ccw_f(uint32_t n) { return F[ccw(getIndex_n(n))]; }

    char find_index_uncommonNode(Face_t* f) {
        return f->hasNode(N[0]) ? f->hasNode(N[1]) ? 2 : 1 : 0;
    }
};
std::ostream& operator <<(std::ostream& stream, const Face_t& f) {
    return stream << f.N[0] << ' ' << f.N[1] << ' ' << f.N[2] << ' ' << f.F[0] << ' ' << f.F[1] << ' ' << f.F[2];
};
std::istream& operator >>(std::istream& stream, Face_t& f) {
    return stream >> f.N[0] >> f.N[1] >> f.N[2] >> f.F[0] >> f.F[1] >> f.F[2];
};

struct KruskalLog {
    KruskalLog(uint32_t x, uint32_t y, uint32_t index) : x(x), y(y), index(index){}
    uint32_t x, y, index;
};

struct Coordinates {
    double lat = 0, lon = 0;
    Coordinates() {}
    Coordinates(const Point & p) : lat(to_double(p.x())), lon(to_double(p.y())) {}
    Coordinates(double lat, double lon) : lat(lat), lon(lon) {}

    const Coordinates operator +(const Coordinates & c) const {return {lat + c.lat, lon + c.lon};};
    const Coordinates operator -(const Coordinates & c) const {return {lat - c.lat, lon - c.lon};};
    const Coordinates operator *(double c) const {return {lat*c, lon*c};};
    const double inline getLength() const {return sqrt((double) lat * lat + lon * lon);}
    void normalize(){
        auto length = getLength();
        lat /= length;
        lon /= length;
        assert(getLength() == 1);
    }
    const uint8_t get_quadrand(const Coordinates c){
        if(c.lat < lat){
            return c.lon < lon ? 2 : 1;
        }else{
            return c.lon < lon ? 3 : 0;
        }
    };
    const inline double distance(const Coordinates & c) const {return sqrt(squared_distance(c));};
    const inline double distance(const Point & p) const {return distance(Coordinates(p));};
    const inline double squared_distance(const Coordinates & c) const {return (lat - c.lat)*(lat - c.lat) + (lon - c.lon)*(lon - c.lon);};
    const inline double squared_distance(const Point & p) const {return squared_distance(Coordinates(p));};
};
// assumes poins are collinear, return true if b is between a and c
bool isOrdered_onLine(const Coordinates & a, const Coordinates & b, const Coordinates & c, double error = 1.0001){
    return (a.lat <= b.lat * error && b.lat <= c.lat * error && a.lat < c.lat) 
        || (a.lat * error >= b.lat && b.lat * error >= c.lat && a.lat > c.lat) 
        || (a.lon <= b.lon * error && b.lon <= c.lon * error && a.lon < c.lon) 
        || (a.lon * error >= b.lon && b.lon * error >= c.lon && a.lon > c.lon);
}
std::ostream& operator <<(std::ostream& stream, const Coordinates& c) {
    return stream << c.lat << " " << c.lon;
};
std::istream& operator >>(std::istream& stream, Coordinates& c) {
    return stream >> c.lat >> c.lon;
};
std::string to_string(const Coordinates & c){
    return "(" + std::to_string(c.lat) + "," + std::to_string(c.lon) + ")";
}
Point to_point(const Coordinates & c){ return Point(c.lat, c.lon); };

struct Matrix2x2 {
    double x0 = 1, y0 = 0, x1 = 0, y1 = 1;
    Matrix2x2() {};
    Matrix2x2(double x0, double y0, double x1, double y1) : x0(x0), y0(y0), x1(x1), y1(y1){};
    Matrix2x2(const Coordinates & c0, const Coordinates & c1) : x0(c0.lat), y0(c0.lon), x1(c1.lat), y1(c1.lon) {};
    Matrix2x2 operator *(Matrix2x2 m) {
        return Matrix2x2(
            x0 * m.x0 + x1 * m.y0, 
            x0 * m.x1 + x1 * m.y1, 
            y0 * m.x0 + y1 * m.y0, 
            y0 * m.x1 + y1 * m.y1
        );
    };
    Coordinates operator *(Coordinates c) {
        return Coordinates(
            x0 * c.lat + x1 * c.lon, 
            y0 * c.lat + y1 * c.lon
        );
    };
};

struct Vector {
    Vector() {}
    Vector(float x, float y, float z) : x(x), y(y), z(z) {}
    float x,y,z;
};

#endif
