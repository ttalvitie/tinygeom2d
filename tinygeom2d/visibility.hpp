#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "geometry.hpp"
#include "domain.hpp"

namespace tinygeom2d {

// Implementation details
namespace visibility_detail {

typedef std::pair<Point, Point> Edge;

struct Event {
    Edge edge;
    bool add;
    
    Point vertex() const {
        return add ? edge.first : edge.second;
    }
};

// Ordering of events based on angle around center
struct EventAngleCmp {
    Point center;
    
    bool operator()(Event a, Event b) {
        Point pointA = a.add ? a.edge.first : a.edge.second;
        Point pointB = b.add ? b.edge.first : b.edge.second;
        if(pointA == pointB) {
            // Removals first
            return (int)a.add < (int)b.add;
        } else {
            return angleLT(center, pointA, center, pointB);
        }
    }
};

// Ordering of events based on angle around center with reference point giving
// zero angle direction
struct EventAngleCmpWithReference {
    Point center;
    Point reference;
    
    bool operator()(Event a, Event b) {
        Point pointA = a.vertex();
        Point pointB = b.vertex();
        if(pointA == pointB) {
            // Removals first
            return (int)a.add < (int)b.add;
        } else {
            int sideA = (int)angleLT(center, pointA, center, reference);
            int sideB = (int)angleLT(center, pointB, center, reference);
            if(sideA == sideB) {
                return angleLT(center, pointA, center, pointB);
            } else {
                return sideA < sideB;
            }
        }
    }
};

// Ordering of edges visible to a ray from center by distance
struct DistLT {
    Point center;
    
    bool operator()(Edge a, Edge b) const {
        if(a == b) {
            return false;
        }
        if(a.first == b.first) {
            return isCCW(a.second, a.first, b.second);
        }
        if(isCCW(center, a.first, b.first)) {
            return isCCW(a.first, b.first, a.second);
        } else {
            return isCCW(b.first, b.second, a.first);
        }
    }
};

// Is b between a and c in CCW ordered angles around center
inline bool anglesOrdered(Point center, Point a, Point b, Point c) {
    if(angleLT(center, a, center, c)) {
        return angleLT(center, a, center, b) && angleLT(center, b, center, c);
    } else {
        return angleLT(center, a, center, b) || angleLT(center, b, center, c);
    }
}

}

// The visibile vertices and edges of a domain from a center point inside the
// domain.
class PointVisibility {
public:
    // The center point, i.e. the point for which the visible vertices and edges
    // are given in this structure.
    Point center() const {
        return center_;
    }
    
    // The vertices visible from the center, ordered by angle (in the sense of
    // angleLT).
    const std::vector<Point>& verts() const {
        return verts_;
    }
    
    // The parts of edges visible from the center. The edges are ordered such
    // that the visible part of edges()[i] is between vertices verts()[i] and
    // verts()[i + 1] (or verts()[0] if i + 1 == size()) when looking from
    // center(). Note that the endpoints of a visible edge are not always
    // visible - the edge may be behind the limiting vertices verts()[i]
    // and verts()[i + 1]. Each edge (a, b) is ordered such that (center, a, b)
    // is CCW oriented, and a appears before b in the boundary of the domain.
    const std::vector<std::pair<Point, Point>>& edges() const {
        return edges_;
    }
    
    // The number of both vertices and parts of edges visible from the center.
    std::size_t size() const {
        return verts_.size();
    }
    
private:
    PointVisibility() { }
    
    Point center_;
    std::vector<Point> verts_;
    std::vector<std::pair<Point, Point>> edges_;
    
    friend PointVisibility computePointVisibility(const Domain& domain, Point center);
};

// Compute the visible vertices and edges in the domain from given center point
// in the interior of the domain, i.e. domain.isInteriorPoint(center) returns
// true.  Throws std::domain_error if the center point is not in the interior.
inline PointVisibility computePointVisibility(const Domain& domain, Point center) {
    using namespace visibility_detail;
    
    const char* notInteriorMsg =
        "tinygeom2d::computePointVisibility: Center point is not an interior "
        "point of the domain";
    
    PointVisibility ret;
    ret.center_ = center;
    
    // We do a rotational ray sweep around the center, maintaining the set of
    // edges intersecting the ray ordered by distance
    
    // First, enumerate all edge addition and removal events and edges that
    // are initially intersecting the sweepray
    std::vector<Edge> initialEdges;
    std::vector<Event> events;
    for(const std::vector<Point>& poly : domain.boundary()) {
        Point a = poly.back();
        if(a == center) {
            throw std::domain_error(notInteriorMsg);
        }
        for(Point b : poly) {
            if(b == center) {
                throw std::domain_error(notInteriorMsg);
            }
            
            Point x = a;
            Point y = b;
            
            if(!isCCW(center, x, y)) {
                std::swap(x, y);
            }
            
            if(yCoordLT(x, center) && yCoordLT(center, y)) {
                initialEdges.emplace_back(x, y);
            }
            events.push_back({{x, y}, true});
            events.push_back({{x, y}, false});
            
            a = b;
        }
    }
    
    // Check interiority by the parity of the number of edges to the right of
    // center
    if(initialEdges.size() % 2 == 0) {
        throw std::domain_error(notInteriorMsg);
    }
    
    // Sort events by angle
    std::sort(events.begin(), events.end(), EventAngleCmp{center});
    
    typedef std::set<Edge, DistLT> Sweepray;
    DistLT distLT{center};
    Sweepray sweepray(distLT);
    
    // Add initial edges to the sweep ray
    for(Edge edge : initialEdges) {
        sweepray.insert(edge);
    }
    
    // Add new vertex and edge to the output only when we move to a new
    // frontmost vertex, in order to not see through vertices between their
    // removal and re-addition
    bool hasPrevFrontVertex = false;
    Point prevFrontVertex;
    auto addVisibility = [&]() {
        ret.verts_.push_back(prevFrontVertex);
        ret.edges_.push_back(*sweepray.begin());
    };
    
    // Rotate the sweepray, handling the edge add and del events
    for(Event event : events) {
        // See if the vertex of this event is frontmost
        bool isFrontVertex;
        if(event.add) {
            isFrontVertex =
                sweepray.empty() || distLT(event.edge, *sweepray.begin());
        } else {
            isFrontVertex = *sweepray.begin() == event.edge;
        }
        
        if(isFrontVertex) {
            // If we move to a new front vertex, add visibility vertex and edge
            if(hasPrevFrontVertex && event.vertex() != prevFrontVertex) {
                addVisibility();
            }
            hasPrevFrontVertex = true;
            prevFrontVertex = event.vertex();
        }
        
        // Change sweepray according to the event
        if(event.add) {
            sweepray.insert(event.edge);
        } else {
            sweepray.erase(event.edge);
        }
    }
    
    // Add the last vertex and edge directly to the right from center
    addVisibility();
    
    return ret;
}

// The visibile vertices and edges of a domain from a center vertex in the
// boundary of the domain.
class VertexVisibility {
public:
    // The center vertex, i.e. the vertex for which the visible vertices and
    // edges are given in this structure.
    Point center() const {
        return center_;
    }
    
    // The vertices visible from the center, in CCW order around the center.
    // The first and last vertices are always the next and previous vertices
    // from the center in the boundary polygon, respectively. The center vertex
    // is not listed as visible.
    const std::vector<Point>& verts() const {
        return verts_;
    }
    
    // The parts of edges visible from the center. The edges are ordered such
    // that the visible part of edges()[i] is between vertices verts()[i] and
    // verts()[i + 1] when looking from center(). Note that the endpoints of a
    // visible edge are not always visible - the edge may be behind the limiting
    // vertices verts()[i] and verts()[i + 1]. Each edge (a, b) is ordered such
    // that (center, a, b) is CCW oriented, and a appears before b in the
    // boundary of the domain. There is always one less elements in edges()
    // than in verts().
    const std::vector<std::pair<Point, Point>>& edges() const {
        return edges_;
    }
    
private:
    VertexVisibility() { }
    
    Point center_;
    std::vector<Point> verts_;
    std::vector<std::pair<Point, Point>> edges_;
    
    friend VertexVisibility computeVertexVisibility(const Domain& domain, Point center);
};

// Compute the vertices and edges in the domain visible from a given vertex of
// the domain. Throws std::domain_error if the center point is not a vertex of
// the boundary.
inline VertexVisibility computeVertexVisibility(const Domain& domain, Point center) {
    using namespace visibility_detail;
    
    VertexVisibility ret;
    ret.center_ = center;
    
    std::pair<std::size_t, std::size_t> id = domain.vertexID(center);
    
    const std::vector<Point> centerPoly = domain.boundary()[id.first];
    Point next = centerPoly[(id.second + 1) % centerPoly.size()];
    Point prev = centerPoly[(id.second + centerPoly.size() - 1) % centerPoly.size()];
    
    // We do a rotational ray sweep around the center, maintaining the set of
    // edges intersecting the ray ordered by distance. The ray starts from the
    // the direction of next and ends in prev, because elsewhere visibility is
    // blocked by the obstacle polygon
    
    // First, enumerate all edge addition and removal events and edges that
    // are initially intersecting the sweepray
    std::vector<Edge> initialEdges;
    std::vector<Event> events;
    for(const std::vector<Point>& poly : domain.boundary()) {
        Point a = poly.back();
        for(Point b : poly) {
            if(a != center && b != center) {
                Point x = a;
                Point y = b;
                
                if(!isCCW(center, x, y)) {
                    std::swap(x, y);
                }
                
                if(x == next || (isCCW(center, x, next) && isCCW(center, next, y))) {
                    initialEdges.emplace_back(x, y);
                }
                
                // Filter out events that happen outside angle range between
                // next and prev
                if(x != next && x != prev && anglesOrdered(center, next, x, prev)) {
                    events.push_back({{x, y}, true});
                }
                if(y != next && y != prev && anglesOrdered(center, next, y, prev)) {
                    events.push_back({{x, y}, false});
                }
            }
            
            a = b;
        }
    }
    
    // Sort events by angle
    std::sort(events.begin(), events.end(), EventAngleCmpWithReference{center, next});
    
    typedef std::set<Edge, DistLT> Sweepray;
    DistLT distLT{center};
    Sweepray sweepray(distLT);
    
    // Add initial edges to the sweep ray
    for(Edge edge : initialEdges) {
        sweepray.insert(edge);
    }
    
    // First vertex is always next
    ret.verts_.push_back(next);
    if(sweepray.empty()) throw 0;
    ret.edges_.push_back(*sweepray.begin());
    
    // Add new vertex and edge to the output only when we move to a new
    // frontmost vertex, in order to not see through vertices between their
    // removal and re-addition
    bool hasPrevFrontVertex = false;
    Point prevFrontVertex;
    auto addVisibility = [&]() {
        ret.verts_.push_back(prevFrontVertex);
        ret.edges_.push_back(*sweepray.begin());
    };
    
    // Rotate the sweepray, handling the edge add and del events
    for(Event event : events) {
        // See if the vertex of this event is frontmost
        bool isFrontVertex;
        if(event.add) {
            isFrontVertex =
                sweepray.empty() || distLT(event.edge, *sweepray.begin());
        } else {
            isFrontVertex = *sweepray.begin() == event.edge;
        }
        
        if(isFrontVertex) {
            // If we move to a new front vertex, add visibility vertex and edge
            if(hasPrevFrontVertex && event.vertex() != prevFrontVertex) {
                addVisibility();
            }
            hasPrevFrontVertex = true;
            prevFrontVertex = event.vertex();
        }
        
        // Change sweepray according to the event
        if(event.add) {
            if(!sweepray.insert(event.edge).second) throw 0;
        } else {
            if(!sweepray.erase(event.edge)) throw 0;
        }
    }
    
    // Add the possible second-to-last vertex
    if(hasPrevFrontVertex) {
        addVisibility();
    }
    
    // Last vertex is always prev
    ret.verts_.push_back(prev);
    
    return ret;
}

// Returns the vector of results of computeVertexVisibility for every vertex of
// the domain.
inline std::vector<VertexVisibility> computeAllVertexVisibilities(const Domain& domain) {
    std::vector<VertexVisibility> ret;
    // TODO: dummy implementation, write fast simultaneous sweepray algorithm
    for(const std::vector<Point>& poly : domain.boundary()) {
        for(Point center : poly) {
            ret.push_back(computeVertexVisibility(domain, center));
        }
    }
    return ret;
}

}
