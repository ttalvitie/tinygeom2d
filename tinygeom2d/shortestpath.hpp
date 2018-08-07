#pragma once

#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

#include "domain.hpp"
#include "geometry.hpp"
#include "intersection.hpp"
#include "visibility.hpp"

namespace tinygeom2d {

// Polygonal domain preprocessed for computing shortest paths between interior
// points of the domain.
class ShortestPathContext {
public:
    // Create context for finding shortest paths in given domain.
    // Time complexity: O(k log n) where n = domain boundary size and
    // k = output size of computeAllVertexVisibilities.
    ShortestPathContext(Domain domain)
        : domain_(std::move(domain))
    {
        // Populate vertIdx_
        vertCount_ = 0;
        for(const std::vector<Point>& poly : domain_.boundary()) {
            Point a = poly[poly.size() - 2];
            Point b = poly.back();
            for(Point c : poly) {
                if(isCCW(c, b, a)) {
                    vertIdx_[b] = vertCount_++;
                    vertPos_.push_back(b);
                }
                a = b;
                b = c;
            }
        }
        
        // Populate graph_
        graph_.resize(2 * vertCount_);
        for(const VertexVisibility& vis : computeAllVertexVisibilities(domain_)) {
            Point vPos = vis.center;
            auto vIt = vertIdx_.find(vPos);
            if(vIt != vertIdx_.end()) {
                std::size_t v = vIt->second;
                
                for(Point xPos : vis.verts) {
                    auto xIt = vertIdx_.find(xPos);
                    if(xIt != vertIdx_.end()) {
                        std::size_t x = xIt->second;
                        
                        if(isTangent_(vPos, xPos) && isTangent_(xPos, vPos)) {
                            std::size_t vDir = v;
                            if(!isLeavingCCW_(vPos, xPos)) {
                                vDir += vertCount_;
                            }
                            std::size_t xDir = x;
                            if(!isEnteringCCW_(xPos, vPos)) {
                                xDir += vertCount_;
                            }
                            graph_[vDir].push_back(xDir);
                        }
                    }
                }
            }
        }
    }
    
    // Returns the domain for which this context was created.
    const Domain& domain() const {
        return domain_;
    }
    
    // Find shortest path between interior points a and b of the domain. The
    // return value is a pair (path, length), where path is the sequence of
    // points defining the path as a polygonal line, and length is the Euclidean
    // length of the path. If there is no path between a and b, the path is
    // empty and the length is infinity. Otherwise, the path has always at least
    // one element, the first one is a and the last one is b. No element is
    // repeated on the path. Inexact double precision floating point numbers are
    // used for path lengths, but the result is still always a valid path and
    // its length should be optimal within floating point accuracy.
    // Throws std::domain_error if a or b is not an interior point of the
    // domain.
    // Time complexity: O(n log n + a log n), where n is the size of the domain
    // boundary and a is the number of visibility graph edges expanded by A*
    // (worst case bound a = O(n^2), very pessimistic in most cases).
    typedef std::pair<std::vector<Point>, double> PathResult;
    PathResult findShortestPath(Point a, Point b) const {
        const char* notInteriorMsg =
            "tinygeom2d::ShortestPathContext::findShortestPath: One of given "
            "points is not an interior point of the domain";
        const double Infinity = std::numeric_limits<double>::infinity();
        
        // Special case: two equal points
        if(a == b) {
            if(!domain_.isInteriorPoint(a)) {
                throw std::domain_error(notInteriorMsg);
            }
            return {{a}, 0.0};
        }
        
        // Special case: direct visibility
        if(isVisible(domain_, a, b)) {
            return {{a, b}, distance(a, b)};
        }
        
        // Rest of the cases: A* in the graph
        PointVisibility aVis = computePointVisibility(domain_, a);
        PointVisibility bVis = computePointVisibility(domain_, b);
        
        struct VertData {
            bool done;
            double dist;
            std::size_t parent;
            bool targetVisible;
        };
        
        const std::size_t NoIdx = -1;
        const std::size_t targetIdx = 2 * vertCount_;
        std::vector<VertData> vertData(
            2 * vertCount_ + 1,
            VertData{false, Infinity, NoIdx, false}
        );
        struct QueueElem {
            double heuristic;
            std::size_t v;
            
            bool operator<(const QueueElem& o) const {
                return heuristic > o.heuristic;
            }
        };
        std::priority_queue<QueueElem> queue;
        
        auto considerParent = [&](std::size_t v, std::size_t parent, double dist) {
            if(!vertData[v].done && dist < vertData[v].dist) {
                vertData[v].dist = dist;
                vertData[v].parent = parent;
                double heuristic = dist + distance(getVertPos_(v), b);
                queue.push({heuristic, v});
            }
        };
        
        for(Point vPos : aVis.verts) {
            auto it = vertIdx_.find(vPos);
            if(it == vertIdx_.end()) {
                continue;
            }
            std::size_t v = it->second;
            if(!isTangent_(vPos, a)) {
                continue;
            }
            if(!isEnteringCCW_(vPos, a)) {
                v += vertCount_;
            }
            considerParent(v, NoIdx, distance(a, getVertPos_(v)));
        }
        for(Point vPos : bVis.verts) {
            auto it = vertIdx_.find(vPos);
            if(it == vertIdx_.end()) {
                continue;
            }
            std::size_t v = it->second;
            if(!isTangent_(vPos, b)) {
                continue;
            }
            if(!isLeavingCCW_(vPos, b)) {
                v += vertCount_;
            }
            vertData[v].targetVisible = true;
        }
        
        while(!queue.empty()) {
            std::size_t v = queue.top().v;
            queue.pop();
            if(vertData[v].done) {
                continue;
            }
            vertData[v].done = true;
            if(v == targetIdx) {
                break;
            }
            for(std::size_t x : graph_[v]) {
                double dist = vertData[v].dist + distance(getVertPos_(v), getVertPos_(x));
                considerParent(x, v, dist);
            }
            if(vertData[v].targetVisible) {
                double dist = vertData[v].dist + distance(getVertPos_(v), b);
                considerParent(targetIdx, v, dist);
            }
        }
        
        if(vertData[targetIdx].done) {
            std::vector<Point> path;
            std::size_t v = targetIdx;
            path.push_back(b);
            while(true) {
                v = vertData[v].parent;
                if(v == NoIdx) {
                    break;
                }
                path.push_back(getVertPos_(v));
            }
            path.push_back(a);
            std::reverse(path.begin(), path.end());
            return {std::move(path), vertData[targetIdx].dist};
        } else {
            return {{}, Infinity};
        }
    }
    
private:
    Domain domain_;
    
    // Number of vertices at convex corners of the domain
    std::size_t vertCount_;
    
    // Mapping of vertices to indices used in graph_.
    std::unordered_map<Point, std::size_t> vertIdx_;
    
    // Inverse mapping of vertIdx_
    std::vector<Point> vertPos_;
    
    // Tangent visibility graph: visibility graph with only vertices that are
    // in convex corners of the domain and edges that come to both their
    // endpoints from the side. These are the only ones that can appear in
    // shortest paths. Each vertex is divided into two parts depending on the
    // orientation of the path: for vertex v, the index of the CCW facing vertex
    // is vertIdx_[v] and the index of the CW facing vertex is
    // vertIdx_[v] + vertCount_.
    std::vector<std::vector<std::size_t>> graph_;
    
    // Get vertex pos for index in graph_
    Point getVertPos_(std::size_t idx) const {
        if(idx >= vertCount_) {
            idx -= vertCount_;
        }
        return vertPos_[idx];
    }
    
    bool isTangent_(Point v, Point p) const {
        Point vNext = domain_.nextVertex(v);
        Point vPrev = domain_.prevVertex(v);
        return isCCW(vNext, v, p) || isCCW(v, vPrev, p);
    }
    bool isLeavingCCW_(Point v, Point p) const {
        return isCCW(v, domain_.prevVertex(v), p);
    }
    bool isEnteringCCW_(Point v, Point p) const {
        return isCCW(domain_.nextVertex(v), v, p);
    }
};

}
