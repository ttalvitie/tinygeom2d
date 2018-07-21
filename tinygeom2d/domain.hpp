#pragma once

#include <cstddef>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "geometry.hpp"
#include "intersection.hpp"

namespace tinygeom2d {

// The identifier of a vertex of the boundary of a polygonal domain. Contains
// indexes to both the polygon and the vertex of the polygon.
struct VertexID {
    std::size_t polygon;
    std::size_t vertex;
};

// A bounded domain in the plane with polygonal boundaries. The domain may be
// disconnected and contain holes.
class Domain {
public:
    // Creates an empty domain.
    Domain() { }
    
    // Given boundary polygons, create a bounded domain, i.e. points inside an
    // odd number of polygons are inside the domain. The orientations of the
    // boundary polygons do not matter. No vertex may appear twice in the
    // boundary, and boundary edges may not intersect each other. Each boundary
    // polygon must contain at least three vertices.
    // Throws std::invalid_argument if the boundary is invalid.
    Domain(std::vector<std::vector<Point>> boundary)
        : boundary_(std::move(boundary))
    {
        // Check for sizes of boundary polygons.
        for(const std::vector<Point>& polygon : boundary_) {
            if(polygon.size() < 3) {
                throw std::invalid_argument("tinygeom2d::Domain::create*: boundary polygon contains less than 3 vertices");
            }
        }
        
        // Check for vertices that appear twice
        std::vector<Point> vertices;
        for(const std::vector<Point>& polygon : boundary_) {
            for(Point vertex : polygon) {
                vertices.push_back(vertex);
            }
        }
        std::sort(vertices.begin(), vertices.end(), yCoordLT);
        for(std::size_t i = 1; i < vertices.size(); ++i) {
            if(vertices[i - 1] == vertices[i]) {
                throw std::invalid_argument("tinygeom2d::Domain::create*: point appears twice on the boundary");
            }
        }
        
        // Check for intersecting boundary edges
        std::vector<std::pair<Point, Point>> edges;
        for(const std::vector<Point>& polygon : boundary_) {
            Point a = polygon.back();
            for(Point b : polygon) {
                edges.emplace_back(a, b);
                a = b;
            }
        }
        if(intersects(std::move(edges))) {
            throw std::invalid_argument("tinygeom2d::Domain::create*: intersecting boundary edges found");
        }
        
        orientBoundary_();
    }
    
    // Returns the boundary polygons of the domain. The polygons are oriented
    // such that outer boundaries are counterclockwise oriented and inner
    // boundaries (holes in the domain) are clockwise oriented. The ordering
    // of the boundary polygons is guaranteed to be the same as in the argument
    // given to the constructor, but the ordering of the vertices in each
    // polygon may have been reversed.
    const std::vector<std::vector<Point>>& boundary() const {
        return boundary_;
    }
    
    // Returns the position of the vertex with given identifier. The indices
    // are not checked for overflows.
    Point vertex(VertexID id) const {
        return boundary_[id.polygon][id.vertex];
    }
    
    // Returns the position of the predecessor of a vertex with given
    // identifier. The indices are not checked for overflows.
    Point prevVertex(VertexID id) const {
        const std::vector<Point>& poly = boundary_[id.polygon];
        if(id.vertex == 0) {
            return poly.back();
        } else {
            return poly[id.vertex - 1];
        }
    }
    
    // Returns the position of the successor of a vertex with given identifier.
    // The indices are not checked for overflows.
    Point nextVertex(VertexID id) const {
        const std::vector<Point>& poly = boundary_[id.polygon];
        if(id.vertex == poly.size() - 1) {
            return poly.front();
        } else {
            return poly[id.vertex + 1];
        }
    }
    
    // Returns the identifier of the predecessor of a vertex with given
    // identifier. The indices are not checked for overflows.
    VertexID prevVertexID(VertexID id) const {
        const std::vector<Point>& poly = boundary_[id.polygon];
        if(id.vertex == 0) {
            return {id.polygon, poly.size() - 1};
        } else {
            return {id.polygon, id.vertex - 1};
        }
    }
    
    // Returns the identifier of the successor of a vertex with given
    // identifier. The indices are not checked for overflows.
    VertexID nextVertexID(VertexID id) const {
        const std::vector<Point>& poly = boundary_[id.polygon];
        if(id.vertex == poly.size() - 1) {
            return {id.polygon, 0};
        } else {
            return {id.polygon, id.vertex + 1};
        }
    }
    
private:
    void orientBoundary_() {
        std::vector<bool> seen(boundary_.size());
        std::vector<bool> outer(boundary_.size());
        
        // Run a sweepline algorithm to determine for each polygon whether there
        // is an even or odd number of points to the left of an arbitrary
        // vertex, which can be used to determine whether the polygon is an
        // inner or outer boundary.
        typedef std::pair<Point, Point> Edge;
        enum EventType {
            // Two edges appear simultaneously on the sweepline (edge1 is to the
            // left of edge2)
            Add,
            // A edge1 is replaced by edge2
            Replace,
            // Two edges disappear simultaneously from the sweepline (edge1 is
            // to the left of edge2)
            Remove
        };
        struct Event {
            Point pos;
            EventType type;
            std::size_t polyIdx;
            Edge edge1;
            Edge edge2;
            
            bool operator<(const Event& other) const {
                return yCoordLT(pos, other.pos);
            }
        };
        
        // Construct all events
        std::vector<Event> events;
        for(std::size_t polyIdx = 0; polyIdx < boundary_.size(); ++polyIdx) {
            const std::vector<Point>& poly = boundary_[polyIdx];
            Point a = poly[poly.size() - 2];
            Point b = poly.back();
            for(Point c : poly) {
                if(yCoordLT(a, b)) {
                    if(yCoordLT(b, c)) {
                        // c
                        // b
                        // a
                        events.push_back({b, Replace, polyIdx, {a, b}, {b, c}});
                    } else {
                        if(isCCW(a, b, c)) {
                            //  b
                            // c a
                            events.push_back({b, Remove, polyIdx, {c, b}, {a, b}});
                        } else {
                            //  b
                            // a c
                            events.push_back({b, Remove, polyIdx, {a, b}, {c, b}});
                        }
                    }
                } else {
                    if(yCoordLT(b, c)) {
                        if(isCCW(a, b, c)) {
                            // a c
                            //  b
                            events.push_back({b, Add, polyIdx, {b, a}, {b, c}});
                        } else {
                            // c a
                            //  b
                            events.push_back({b, Add, polyIdx, {b, c}, {b, a}});
                        }
                    } else {
                        // a
                        // b
                        // c
                        events.push_back({b, Replace, polyIdx, {c, b}, {b, a}});
                    }
                }
                
                a = b;
                b = c;
            }
        }
        
        // Sort the events by y-coordinate
        std::sort(events.begin(), events.end());
        
        // Our sweepline is an x-ordered set of pairs of edge and a boolean that
        // tells whether the edge has an odd number of edges to the left.
        typedef std::pair<Edge, bool> Elem;
        auto edgeOrderX = [&](Elem a, Elem b) {
            return segmentLeftOfAtBottom(a.first, b.first);
        };
        typedef std::set<Elem, decltype(edgeOrderX)> Sweepline;
        typedef Sweepline::iterator Iter;
        Sweepline sweepline(edgeOrderX);
        
        // Run the sweepline
        for(const Event& event : events) {
            if(event.type == Add) {
                Iter pos = sweepline.lower_bound({event.edge1, false});
                bool odd;
                if(pos == sweepline.begin()) {
                    odd = false;
                } else {
                    Iter prev = pos;
                    --prev;
                    odd = !prev->second;
                }
                pos = sweepline.insert(pos, {event.edge1, odd});
                sweepline.insert(pos, {event.edge2, !odd});
                
                // At the first Add event, we determine the correct orientation
                // for each polygon
                if(!seen[event.polyIdx]) {
                    seen[event.polyIdx] = true;
                    outer[event.polyIdx] = !odd;
                }
            }
            if(event.type == Replace) {
                Iter pos = sweepline.find({event.edge1, false});
                bool odd = pos->second;
                Iter next = pos;
                ++next;
                sweepline.erase(pos);
                sweepline.insert(next, {event.edge2, odd});
            }
            if(event.type == Remove) {
                Iter pos1 = sweepline.find({event.edge1, false});
                Iter pos2 = pos1;
                ++pos2;
                sweepline.erase(pos1);
                sweepline.erase(pos2);
            }
        }
        
        // Orient polygons
        for(std::size_t polyIdx = 0; polyIdx < boundary_.size(); ++polyIdx) {
            std::vector<Point>& poly = boundary_[polyIdx];
            
            // Find out current orientation by finding lowest vertex and looking
            // at the orientation of its neighbors
            Point a = poly[poly.size() - 2];
            Point b = poly.back();
            Point bestA = a;
            Point bestB = b;
            Point bestC = poly.front();
            for(Point c : poly) {
                if(yCoordLT(b, bestB)) {
                    bestA = a;
                    bestB = b;
                    bestC = c;
                }
                
                a = b;
                b = c;
            }
            bool ccw = isCCW(bestA, bestB, bestC);
            
            if(outer[polyIdx] != ccw) {
                std::reverse(poly.begin(), poly.end());
            }
        }
    }
    
    std::vector<std::vector<Point>> boundary_;
};

}
