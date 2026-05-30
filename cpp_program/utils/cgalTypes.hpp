#ifndef CGAL_My_TYPES_HPP
#define CGAL_My_TYPES_HPP

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Segment_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     K_Exact;
typedef K_Exact::FT                                                     Exact_Distance;
typedef K_Exact::Point_2                                                Exact_Point;
const Exact_Distance inline getDist_exact_exact(const Exact_Point & p0, const Exact_Point & p1) { return CGAL::sqrt(CGAL::squared_distance(p0, p1)); };

#ifdef CGAL_USE_EPECK
    typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt K;
    // typedef CGAL::Lazy_exact_nt<Exact_Distance>
    // template <typename Lazy_FT>
    // double cast_to_double(const Lazy_FT& x) {
    // x.exact();
    // return to_double(x);
    // }
#else
    typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
#endif

typedef K::FT                                                           NT;
typedef K::Line_2                                                       Line_2;
typedef K::Segment_2                                                    Segment;
typedef K::Point_2                                                      Point;
typedef CGAL::Aff_transformation_2<K>                                   Transformation;

#ifdef CGAL_USE_INTERVALL
    #include <CGAL/Interval_nt.h>
    typedef CGAL::Interval_nt<true>             NT_Distance; //default true, false is faster but need Protector
    const NT_Distance inline getSqrdDist_exact(const Point & p0, const Point & p1) { return (NT_Intervall(p0.x()) * NT_Intervall(p1.x())) + (NT_Intervall(p0.y()) * NT_Intervall(p1.y()));  };
#else
    typedef NT                                                          NT_Distance;
    const NT_Distance inline getSqrdDist_exact(const Point & p0, const Point & p1) { return CGAL::squared_distance(p0, p1); };
#endif

#ifdef CGAL_USE_EPECK
const inline double to_double(const NT_Distance& nt){
    auto cgal_to_double = CGAL::to_double(nt);
    // auto call_doubleValue = nt.doubleValue();
    // auto cgal_to_double_approx_first = CGAL::to_double(nt.approx().doubleValue());
    // assert(cgal_to_double == cgal_to_double_approx_first && cgal_to_double_approx_first == call_doubleValue);
    return cgal_to_double;
};
#else
const inline double to_double(const NT_Distance& nt){return nt;};
#endif


const double inline getSqrdDist_inexact(const Point & p0, const Point & p1) {return to_double(CGAL::squared_distance(p0, p1));};
const NT_Distance inline getDist_exact(const Point & p0, const Point & p1) { return CGAL::sqrt(getSqrdDist_exact(p0, p1)); };
const double inline getDist_inexact(const Point & p0, const Point & p1) {return CGAL::sqrt(getSqrdDist_inexact(p0, p1));};

inline Point calc_intersection(const Point & g, const Point & bound, const Point & n0, const Point & n1){
    try{
        using boost::get; using std::get;
        return get<Point>(CGAL::intersection(Line_2(n0, n1), Line_2(g, bound)).value());
    }catch(const std::exception & e){
        return getDist_inexact(g, n0) < getDist_inexact(g, n1) ? n0 : n1;
    }
};

std::string to_string(const Point & p){return "(" + std::to_string(to_double(p.x())) + "," + std::to_string(to_double(p.y())) + ")";}

#endif