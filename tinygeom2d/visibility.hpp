#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "geometry.hpp"
#include "domain.hpp"

namespace tinygeom2d {

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
    // center(). Please note that the endpoints of a visible edge are not always
    // visible - the edge may be behind the limiting vertices verts()[i]
    // and verts()[i + 1]. Each edge (a, b) such that (center, a, b) is CCW
    // oriented, and a appears before b in the boundary of the domain.
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
    const char* notInteriorMsg =
        "tinygeom2d::computePointVisibility: Center point is not an interior "
        "point of the domain";
    
    PointVisibility ret;
    ret.center_ = center;
    
    typedef std::pair<Point, Point> Edge;
    
    // We do a rotational ray sweep around the center, maintaining the set of
    // edges intersecting the ray ordered by distance
    
    // First, enumerate all edge addition and removal events and edges that
    // are initially intersecting the sweepray
    
    // In events, (edge, true) means addition and (edge, false) means removal
    typedef std::pair<Edge, bool> Event;
    
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
    std::sort(events.begin(), events.end(), [&](Event a, Event b) {
        Point pointA = a.second ? a.first.first : a.first.second;
        Point pointB = b.second ? b.first.first : b.first.second;
        if(pointA == pointB) {
            // Removals first
            return (int)a.second < (int)b.second;
        } else {
            return angleLT(center, pointA, center, pointB);
        }
    });
    
    // Comparison of edges visible to the ray by distance
    auto distLT = [&](Edge a, Edge b) {
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
    };
    
    typedef std::set<Edge, decltype(distLT)> Sweepray;
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
        bool isAddEvent = event.second;
        Edge eventEdge = event.first;
        Point eventVertex = isAddEvent ? eventEdge.first : eventEdge.second;
        
        // See if the vertex of this event is frontmost
        bool isFrontVertex;
        if(isAddEvent) {
            isFrontVertex =
                sweepray.empty() || distLT(eventEdge, *sweepray.begin());
        } else {
            isFrontVertex = *sweepray.begin() == eventEdge;
        }
        
        if(isFrontVertex) {
            // If we move to a new front vertex, add visibility vertex and edge
            if(hasPrevFrontVertex && eventVertex != prevFrontVertex) {
                addVisibility();
            }
            hasPrevFrontVertex = true;
            prevFrontVertex = eventVertex;
        }
        
        // Change sweepray according to the event
        if(isAddEvent) {
            sweepray.insert(eventEdge);
        } else {
            sweepray.erase(eventEdge);
        }
    }
    
    // Add the last vertex and edge directly to the right from center
    addVisibility();
    
    return ret;
}

}
