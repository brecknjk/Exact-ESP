#ifndef MATHLIB_H
#define MATHLIB_H

#include <memory>
#include <string>
#include <stdint.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>

#include "types.hpp"

constexpr uint32_t RADIUS_EARTH = 6371;
constexpr uint32_t EDGECOST_MULT = 100;
constexpr uint32_t MAX_EDGE_RANGE = 30 * EDGECOST_MULT;
constexpr float INF_FLOAT = std::numeric_limits<float>::infinity();
constexpr float MAX_FLOAT = std::numeric_limits<float>::max();
constexpr float MIN_FLOAT = std::numeric_limits<float>::min();
constexpr float EPS_FLOAT = std::numeric_limits<float>::epsilon();
constexpr double INF_DOUBLE = std::numeric_limits<double>::infinity();
constexpr double MAX_DOUBLE = std::numeric_limits<double>::max();
constexpr double MIN_DOUBLE = std::numeric_limits<double>::min();
constexpr double EPS_DOUBLE = std::numeric_limits<float>::epsilon();

Vector o = Vector(0,0,0);
std::vector<std::pair<double, double>> const directions = {std::pair(1.0,1.0),std::pair(1.0,-1.0),std::pair(-1.0,1.0),std::pair(-1.0,-1.0)};

uint8_t int_log2(uint32_t x){
    assert(x != 0);
    return 31 - (__builtin_clz(x));
};

const float inline to_lowerFloat(double d){
    assert(d >= 0);
    float f = d;
    return ((double)f > d) ? nextafterf(f, 0) : f;
}

double inline sqrd_r_to_d(double sqrd_radius) { return sqrt(sqrd_radius * 2); }

template<typename T>
double calc_mean(const std::vector<T> &vec) {
    if (vec.size() == 0) return 0.0;
    return (double) std::accumulate(vec.begin(), vec.end(), 0.0) / (double) vec.size();
}

//https://stackoverflow.com/questions/33268513/calculating-standard-deviation-variance-in-c
template<typename T>
double variance(const std::vector<T> &vec) {
    const size_t sz = vec.size();
    if (sz <= 1) return 0.0;

    // Calculate the mean
    const double m = calc_mean<T>(vec);

    // Now calculate the variance
    auto variance_func = [&m, &sz](T accumulator, const T& val) {
        return accumulator + ((val - m)*(val - m) / (sz - 1));
    };

    return std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

template<typename T>
T stDeviation(const std::vector<T> &vec) {
    return sqrt(variance(vec));
}

// source: https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
Vector inline generateRandomPoint() {
    double z = (((double) std::rand() / (RAND_MAX)) - 0.5) * 2 * RADIUS_EARTH;
    double phi = ((double) std::rand() / (RAND_MAX)) * 2 * M_PI;
    double factor = sqrt(RADIUS_EARTH * RADIUS_EARTH - z * z);
    return Vector(factor * cos(phi), factor * sin(phi), z);
}

// source: http://www.movable-type.co.uk/scripts/latlong-vectors.html
Vector inline toVector(const Coordinates & c) {
    double lat_rad = c.lat * M_PI / 180.0;
    double lon_rad = c.lon * M_PI / 180.0;
    return Vector(RADIUS_EARTH * cos(lat_rad)*cos(lon_rad), RADIUS_EARTH * cos(lat_rad)*sin(lon_rad), RADIUS_EARTH * sin(lat_rad));
}

// source: http://www.movable-type.co.uk/scripts/latlong-vectors.html
Coordinates inline toCoords(const Vector & v){
    return Coordinates(180 * M_1_PI * (atan2(v.z, sqrt(v.x * v.x + v.y * v.y))), 180 * M_1_PI * (atan2(v.y, v.x)));
}
    
Vector inline getCross(const Vector v1, const Vector v2) {
    return Vector(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}
    
double inline getProduct(const Vector & v1, const Vector & v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
    
double inline getLength(const Vector & v) {
    return sqrt((double) v.x * v.x + v.y * v.y + v.z * v.z);
}

double inline getLength(const Coordinates & c) {
    return sqrt((double) c.lat * c.lat + c.lon * c.lon);
}
Coordinates inline normalize(Coordinates & c){
    return c * (1 / getLength(c));
}
    
bool inline inBound(const Coordinates & c_min, const Coordinates & c_max, const Coordinates & c_test) {
    return c_min.lat <= c_test.lat < c_max.lat and c_min.lon <= c_test.lon < c_max.lon;
}

// source: http://www.movable-type.co.uk/scripts/latlong-vectors.html
uint32_t inline getDist(const Vector & v1, const Vector & v2) {
    return EDGECOST_MULT * RADIUS_EARTH * atan2(getLength(getCross(v1, v2)), getProduct(v1, v2));
}

// source: http://www.movable-type.co.uk/scripts/latlong-vectors.html
double inline getDist_double(const Vector & v1, const Vector & v2) {
    return atan2(getLength(getCross(v1, v2)), getProduct(v1, v2));
}
    
//check if c und d are on oppostite sides of (OAxAB)*AC < or > 0
bool inline checkIntersect(const Vector & a, const Vector & b, const Vector & c, const Vector & d) {
    auto x1 = getCross(Vector(a.x - o.x, a.y - o.y, a.z - o.z), Vector(b.x - a.x, b.y - a.y, b.z - a.z));
    auto w_c = getProduct(x1, Vector(c.x - a.x, c.y - a.y, c.z - a.z));
    auto w_d = getProduct(x1, Vector(d.x - a.x, d.y - a.y, d.z - a.z));
    //if(w_c != 0 or w_d != 0) hitCorner = True;
    if ((w_c < 0) == (w_d < 0)){
        return false;
    }
        
    auto x2 = getCross(Vector(c.x - o.x, c.y - o.y, c.z - o.z), Vector(d.x - c.x, d.y - c.y, d.z - c.z));
    auto w_a = getProduct(x2, Vector(a.x - c.x, a.y - c.y, a.z - c.z));
    auto w_b = getProduct(x2, Vector(b.x - c.x, b.y - c.y, b.z - c.z));
    //if(w_a != 0 and w_b != 0) hitCorner = True;
    if ((w_a < 0) == (w_b < 0)){
        return false;
    }
    return true;
};

bool inline checkIntersect_opt(const Vector & a, const Vector & b, const Vector & c, const Vector & d, const Vector & abX, const Vector & cdX) {
    auto cMinusA = Vector(c.x - a.x, c.y - a.y, c.z - a.z);
    auto w_c = getProduct(abX, cMinusA);
    auto w_d = getProduct(abX, Vector(d.x - a.x, d.y - a.y, d.z - a.z));
    //if(w_c != 0 or w_d != 0) hitCorner = True;
    if ((w_c < 0) == (w_d < 0)){
        return false;
    }
        
    auto w_a = -getProduct(cdX, cMinusA);
    auto w_b = getProduct(cdX, Vector(b.x - c.x, b.y - c.y, b.z - c.z));
    //if(w_a != 0 and w_b != 0) hitCorner = True;
    if ((w_a < 0) == (w_b < 0)){
        return false;
    }
    return true;
};

// Funktion zur Berechnung der Mercator-Projektion
Coordinates bounded_mercator(const Coordinates & geo, uint32_t mapWidth=2000, uint32_t mapHeight=1000) {
    Coordinates cart;
    assert(geo.lon >= -180.1 && geo.lon <= 180.1);
    assert(geo.lat >= -89 && geo.lat <= 89);
    cart.lat = (geo.lon+180)*(mapWidth/360);
    double latRad = geo.lat*M_PI/180;
    double mercN = log(tan((M_PI/4)+(latRad/2)));
    cart.lon = (mapHeight/2)-(mapHeight*mercN/(2*M_PI));
    assert(!isnan(cart.lat) && !isnan(cart.lon));
    return cart;
}
Point bounded_mercator(const Point & geo, uint32_t mapWidth=2000, uint32_t mapHeight=1000) {
    Coordinates c = bounded_mercator(Coordinates(to_double(geo.x()), to_double(geo.y())), mapWidth, mapHeight);
    return Point(c.lat, c.lon);
}

Coordinates rev_bounded_mercator(const Coordinates & cart, uint32_t mapWidth=2000, uint32_t mapHeight=1000) {
    Coordinates geo;
    double mercN = ((mapHeight/2) - cart.lon)* (2*M_PI) / mapHeight;
    double latRad = (atan(exp(mercN)) - (M_PI/4)) * 2;
    geo.lat = latRad * 180 / M_PI;
    geo.lon = (cart.lat / (mapWidth / 360)) - 180;
    return geo;
}
Point rev_bounded_mercator(const Point & cart, uint32_t mapWidth=2000, uint32_t mapHeight=1000) {
    Coordinates c = rev_bounded_mercator(Coordinates(to_double(cart.x()), to_double(cart.y())), mapWidth, mapHeight);
    return Point(c.lat, c.lon);
}

Vector inline getCross_intermediate(const Vector & a, const Vector & b) {
    return getCross(Vector(a.x - o.x, a.y - o.y, a.z - o.z), Vector(b.x - a.x, b.y - a.y, b.z - a.z));
};

#endif