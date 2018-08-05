#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "domain.hpp"
#include "geometry.hpp"
#include "intersection.hpp"

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
    friend std::vector<VertexVisibility> computeAllVertexVisibilities(const Domain& domain);
};

// Compute the vertices and edges in the domain visible from a given vertex of
// the domain. Throws std::domain_error if the center point is not a vertex of
// the boundary.
inline VertexVisibility computeVertexVisibility(const Domain& domain, Point center) {
    using namespace visibility_detail;
    
    VertexVisibility ret;
    ret.center_ = center;
    
    std::pair<std::size_t, std::size_t> id = domain.vertexID(center);
    
    const std::vector<Point>& centerPoly = domain.boundary()[id.first];
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
            sweepray.insert(event.edge);
        } else {
            sweepray.erase(event.edge);
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
// the domain. The implementation is faster than calling computeVertexVisibility
// for every vertex separately.
inline std::vector<VertexVisibility> computeAllVertexVisibilities(const Domain& domain) {
    // Enumerate vertices, and create doubly linked lists of them along the
    // domain boundary.
    std::vector<Point> vertPos;
    std::vector<std::size_t> prevVertIdx;
    std::vector<std::size_t> nextVertIdx;
    for(const std::vector<Point>& poly : domain.boundary()) {
        std::size_t prevIdx = vertPos.size() + poly.size() - 1;
        for(Point pos : poly) {
            std::size_t idx = vertPos.size();
            vertPos.push_back(pos);
            prevVertIdx.push_back(prevIdx);
            prevIdx = idx;
        }
    }
    std::size_t vertCount = vertPos.size();
    nextVertIdx.resize(vertCount);
    for(std::size_t v = 0; v < vertCount; ++v) {
        nextVertIdx[prevVertIdx[v]] = v;
    }
    
    // We run a simultaneous ray sweep around every vertex. We only need to
    // maintain the mapping from each vertex to the edge the ray from that
    // vertex hits first, and the only possible changes are when two rays
    // coincide or a ray moves to the next edge.
    
    // The current hit edge for each ray. The edge i means edge between vertices
    // prevVertIdx[i] and i.
    const std::size_t NoIdx = -1;
    std::vector<std::size_t> currentEdge(vertCount, NoIdx);
    
    // All rays initially point to angle zero, and we use a sweepline algorithm
    // to determine the initially visible edge for every ray
    {
        struct Event {
            // Events at same position are processed in the enum order
            enum {DelEdge, Vertex, AddEdge} type;
            
            union {
                std::size_t vertIdx;
                std::size_t edge;
            };
        };
        std::vector<Event> events;
        for(std::size_t edge = 0; edge < vertCount; ++edge) {
            // We only need upward edges
            if(yCoordLT(vertPos[prevVertIdx[edge]], vertPos[edge])) {
                events.push_back(Event{Event::AddEdge, edge});
                events.push_back(Event{Event::DelEdge, edge});
            }
        }
        for(std::size_t v = 0; v < vertCount; ++v) {
            // Only add vertices for which the ray to angle zero is not blocked
            if(angleLT(
                vertPos[v],
                vertPos[prevVertIdx[v]],
                vertPos[v],
                vertPos[nextVertIdx[v]]
            )) {
                events.push_back(Event{Event::Vertex, v});
            }
        }
        
        auto eventPos = [&](Event e) {
            switch(e.type) {
                case Event::Vertex: return vertPos[e.vertIdx];
                case Event::AddEdge: return vertPos[prevVertIdx[e.edge]];
                case Event::DelEdge: return vertPos[e.edge];
                default: return Point(); // Never happens
            }
        };
        
        std::sort(events.begin(), events.end(), [&](Event a, Event b) {
            Point aPos = eventPos(a);
            Point bPos = eventPos(b);
            
            if(aPos == bPos) {
                // Respect Event::type ordering
                return (int)a.type < (int)b.type;
            } else {
                return yCoordLT(aPos, bPos);
            }
        });
        
        auto leftOf = [&](std::size_t edge1, std::size_t edge2) {
            Point a1 = vertPos[prevVertIdx[edge1]];
            Point b1 = vertPos[edge1];
            Point a2 = vertPos[prevVertIdx[edge2]];
            Point b2 = vertPos[edge2];
            if(yCoordLT(b1, a1)) {
                std::swap(a1, b1);
            }
            if(yCoordLT(b2, a2)) {
                std::swap(a2, b2);
            }
            return segmentLeftOfAtBottom({a1, b1}, {a2, b2});
        };
        std::set<std::size_t, decltype(leftOf)> sweepline(leftOf);
        for(Event event : events) {
            switch(event.type) {
                case Event::AddEdge: {
                    sweepline.insert(event.edge);
                } break;
                
                case Event::DelEdge: {
                    sweepline.erase(event.edge);
                } break;
                
                case Event::Vertex: {
                    auto it = sweepline.lower_bound(event.vertIdx);
                    if(it != sweepline.end()) {
                        currentEdge[event.vertIdx] = *it;
                    }
                } break;
            };
        }
    }
    
    // As we simultaneously rotate the rays centered at each vertex, we maintain
    // a set of events of the ray changes ordered by the ray angle at which the
    // event happens
    struct Event {
        enum {
            AddRay, // The ray becomes unblocked
            EdgeEnds, // The ray reaches the endpoint of the current edge
            Overtake, // The ray overtakes the next ray on the same edge
        } type;
        
        // The center vertex of the event ray
        std::size_t center;
        
        // Angle of the ray at which the event occurs
        std::pair<Point, Point> angle;
        
        // It can be proven that there is only one possible event for every
        // angle
        bool operator<(const Event& o) const {
            return angleLT(angle.first, angle.second, o.angle.first, o.angle.second);
        }
    };
    std::set<Event> events;
    
    // Add all AddRay events in advance
    for(std::size_t v = 0; v < vertCount; ++v) {
        events.insert({Event::AddRay, v, {vertPos[v], vertPos[nextVertIdx[v]]}});
    }
    
    // For every edge, a doubly linked list of rays such that the edge is
    // currentEdge[ray]. Ordered by the position of the ray along the edge
    std::vector<std::size_t> rayListFirst(vertCount, NoIdx); // by edge]
    std::vector<std::size_t> rayListPrev(vertCount, NoIdx); // by vertex idx
    std::vector<std::size_t> rayListNext(vertCount, NoIdx); // by vertex idx
    
    // Fill the linked lists
    {
        std::vector<std::size_t> order(vertCount);
        for(std::size_t v = 0; v < vertCount; ++v) {
            order[v] = v;
        }
        sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return yCoordLT(vertPos[b], vertPos[a]);
        });
        for(std::size_t v : order) {
            std::size_t edge = currentEdge[v];
            if(edge != NoIdx) {
                std::size_t& first = rayListFirst[edge];
                if(first != NoIdx) {
                    rayListNext[v] = first;
                    rayListPrev[first] = v;
                }
                first = v;
            }
        }
    }
    
    // Get the next event for ray from v. Returns pair (has event, next event)
    auto nextEvent = [&](std::size_t v) -> std::pair<bool, Event> {
        if(currentEdge[v] == NoIdx) {
            return {false, {}};
        }
        if(rayListNext[v] == NoIdx) {
            return {true, Event{
                Event::EdgeEnds, v,
                {vertPos[v], vertPos[currentEdge[v]]}
            }};
        } else {
            bool hasEvent = isCCW(
                vertPos[v],
                vertPos[rayListNext[v]],
                vertPos[currentEdge[v]]
            );
            return {hasEvent, Event{
                Event::Overtake, v,
                {vertPos[v], vertPos[rayListNext[v]]}
            }};
        }
    };
    
    // We only add events that have angle larger than the previous event angle,
    // so that the rotation goes forward
    std::pair<Point, Point> prevEventAngle;
    bool hasPrevEvent = false;
    
    auto addNextEvent = [&](std::size_t v) {
        std::pair<bool, Event> p = nextEvent(v);
        if(p.first && (
            !hasPrevEvent ||
            angleLT(
                prevEventAngle.first, prevEventAngle.second,
                p.second.angle.first, p.second.angle.second
            )
        )) {
            events.insert(p.second);
        }
    };
    auto delNextEvent = [&](std::size_t v) {
        std::pair<bool, Event> p = nextEvent(v);
        if(p.first) {
            events.erase(p.second);
        }
    };
    
    // Prepopulate next events for all rays
    for(std::size_t v = 0; v < vertCount; ++v) {
        addNextEvent(v);
    }
    
    // Special case: empty domain
    if(events.empty()) {
        return std::vector<VertexVisibility>();
    }
    
    // Convenience functions for modifying ray linked lists
    auto addRayAfter = [&](std::size_t v, std::size_t x) {
        delNextEvent(x);
        currentEdge[v] = currentEdge[x];
        rayListPrev[v] = x;
        rayListNext[v] = rayListNext[x];
        rayListNext[x] = v;
        if(rayListNext[v] != NoIdx) {
            rayListPrev[rayListNext[v]] = v;
        }
        addNextEvent(x);
        addNextEvent(v);
    };
    auto addRayFirst = [&](std::size_t v, std::size_t edge) {
        currentEdge[v] = edge;
        std::size_t& first = rayListFirst[edge];
        if(first != NoIdx) {
            rayListPrev[first] = v;
        }
        rayListNext[v] = first;
        rayListPrev[v] = NoIdx;
        first = v;
        addNextEvent(v);
    };
    auto delRay = [&](std::size_t v) {
        delNextEvent(v);
        std::size_t prev = rayListPrev[v];
        if(prev == NoIdx) {
            rayListFirst[currentEdge[v]] = rayListNext[v];
            if(rayListNext[v] != NoIdx) {
                rayListPrev[rayListNext[v]] = NoIdx;
            }
        } else {
            delNextEvent(prev);
            rayListNext[prev] = rayListNext[v];
            if(rayListNext[v] != NoIdx) {
                rayListPrev[rayListNext[v]] = prev;
            }
            addNextEvent(prev);
        }
        currentEdge[v] = NoIdx;
    };
    
    // Handle events in order, and collect seen edge-vertex-pairs
    std::vector<
        std::vector<std::pair<std::size_t, std::size_t>>
    > edgeVertPairs(vertCount);
    while(!events.empty()) {
        Event event = *events.begin();
        events.erase(events.begin());
        
        prevEventAngle = event.angle;
        hasPrevEvent = true;
        
        std::size_t v = event.center;
        
        switch(event.type) {
            case Event::AddRay: {
                std::size_t x = nextVertIdx[v];
                edgeVertPairs[v].emplace_back(NoIdx, x);
                if(currentEdge[x] == NoIdx) {
                    addRayFirst(v, nextVertIdx[x]);
                } else {
                    addRayAfter(v, x);
                }
            } break;
            
            case Event::EdgeEnds: {
                std::size_t x = currentEdge[v];
                edgeVertPairs[v].emplace_back(x, x);
                delRay(v);
                if(currentEdge[x] == NoIdx) {
                    std::size_t y = nextVertIdx[x];
                    if(y != v) {
                        addRayFirst(v, y);
                    }
                } else {
                    if(x != prevVertIdx[v]) {
                        addRayAfter(v, x);
                    }
                }
            } break;
            
            case Event::Overtake: {
                std::size_t x = rayListNext[v];
                edgeVertPairs[v].emplace_back(currentEdge[v], x);
                delRay(v);
                if(x != prevVertIdx[v]) {
                    addRayFirst(v, nextVertIdx[x]);
                }
            } break;
        }
    }
    
    std::vector<VertexVisibility> ret(vertCount, VertexVisibility());
    for(std::size_t v = 0; v < vertCount; ++v) {
        ret[v].center_ = vertPos[v];
        
        // Put the values in edgeVertPairs to ret such that the edge marked by
        // NoIdx is between the end and the beginning in the verts and edges
        // arrays
        std::size_t z = 0;
        while(edgeVertPairs[v][z].first != NoIdx) {
            ++z;
        }
        
        for(std::size_t i = z; i < edgeVertPairs[v].size(); ++i) {
            ret[v].verts_.push_back(vertPos[edgeVertPairs[v][i].second]);
        }
        for(std::size_t i = 0; i < z; ++i) {
            ret[v].verts_.push_back(vertPos[edgeVertPairs[v][i].second]);
        }
        for(std::size_t i = z + 1; i < edgeVertPairs[v].size(); ++i) {
            std::size_t edge = edgeVertPairs[v][i].first;
            ret[v].edges_.emplace_back(vertPos[prevVertIdx[edge]], vertPos[edge]);
        }
        for(std::size_t i = 0; i < z; ++i) {
            std::size_t edge = edgeVertPairs[v][i].first;
            ret[v].edges_.emplace_back(vertPos[prevVertIdx[edge]], vertPos[edge]);
        }
    }
    
    return ret;
}

}
