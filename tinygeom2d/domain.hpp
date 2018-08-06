#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "geometry.hpp"
#include "intersection.hpp"

namespace tinygeom2d {

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
    // Time complexity: O(n log n), where n = input size.
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
        initVertexMap_();
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
    
    // Pair identifying a vertex in the boundary, consisting of the index of the
    // boundary polygon and the index of the vertex in the polygon.
    typedef std::pair<std::size_t, std::size_t> IdxPair;
    
    // Returns a mapping from the vertices of the boundary to their indices in
    // the boundary polygons. For any element (v -> (a, b)) in the map,
    // boundary()[a][b] == v.
    const std::unordered_map<Point, IdxPair>& vertexMap() const {
        return vertexMap_;
    }
    
    // Returns the indices of a boundary vertex in the boundary polygons.
    // Equivalent to vertexMap().find(vertex)->second if the vertex is a vertex
    // of the boundary.
    // Throws std::domain_error if given vertex is not a vertex of the boundary.
    IdxPair vertexID(Point vertex) const {
        auto it = vertexMap_.find(vertex);
        if(it == vertexMap_.end()) {
            throw std::domain_error("tinygeom2d::Domain::vertexID: given vertex is not in the boundary of the domain");
        }
        return it->second;
    }
    
    // Returns the position of the next vertex from given vertex on the boundary
    // polygon. 
    // Throws std::domain_error if given vertex is not a vertex of the boundary.
    Point nextVertex(Point vertex) const {
        IdxPair idxs = vertexID(vertex);
        if(idxs.second == boundary_[idxs.first].size() - 1) {
            return boundary_[idxs.first].front();
        } else {
            return boundary_[idxs.first][idxs.second + 1];
        }
    }
    
    // Returns the position of the previous vertex from given vertex on the
    // boundary polygon. 
    // Throws std::domain_error if given vertex is not a vertex of the boundary.
    Point prevVertex(Point vertex) const {
        IdxPair idxs = vertexID(vertex);
        if(idxs.second == 0) {
            return boundary_[idxs.first].back();
        } else {
            return boundary_[idxs.first][idxs.second - 1];
        }
    }
    
    // Returns true if given point is an interior point of the domain. The
    // vertices of the domain do not count as interior points.
    // Time complexity: O(n), where n = boundary size.
    bool isInteriorPoint(Point point) const {
        // Check whether the number of edges that the ray shot from point to
        // the left is odd.
        bool ret = false;
        for(const std::vector<Point>& poly : boundary_) {
            Point a = poly.back();
            if(a == point) {
                return false;
            }
            
            for(Point b : poly) {
                if(b == point) {
                    return false;
                }
                
                Point x = a;
                Point y = b;
                if(yCoordLT(y, x)) {
                    std::swap(x, y);
                }
                
                if(yCoordLT(x, point) && yCoordLT(point, y) && isCCW(x, point, y)) {
                    ret = !ret;
                }
                
                a = b;
            }
        }
        return ret;
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
            return intersection_detail::segmentLeftOfAtBottom(a.first, b.first);
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
    
    void initVertexMap_() {
        for(std::size_t polyIdx = 0; polyIdx < boundary_.size(); ++polyIdx) {
            const std::vector<Point>& poly = boundary_[polyIdx];
            for(std::size_t vertIdx = 0; vertIdx < poly.size(); ++vertIdx) {
                vertexMap_[poly[vertIdx]] = {polyIdx, vertIdx};
            }
        }
    }
    
    std::vector<std::vector<Point>> boundary_;
    std::unordered_map<Point, std::pair<std::size_t, std::size_t>> vertexMap_;
};

}
