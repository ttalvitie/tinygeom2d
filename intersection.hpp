#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <utility>
#include <vector>

#include "geometry.hpp"

namespace tinygeom2d {

// Returns true if the interiors of the segments a-b and c-d intersect, i.e.,
// there is a point that is on both segments that is not an endpoint of one of
// the segments.
inline bool intersects(Point a, Point b, Point c, Point d) {
    if(a == c || a == d || b == c || b == d) {
        if(((a == c && b == d) || (a == d && b == c)) && a != b) {
            return true;
        }
        return false;
    }
    return isCCW(a, b, c) != isCCW(a, b, d) && isCCW(c, d, a) != isCCW(c, d, b);
}

// Returns true if two of the given segments intersect.
inline bool intersects(std::vector<std::pair<Point, Point>> segments) {
    typedef std::pair<Point, Point> Segment;
    
    // Remove zero-length segments
    std::size_t s = segments.size();
    for(std::size_t i = 0; i < s; ++i) {
        if(segments[i].first == segments[i].second) {
            segments[i] = segments[--s];
        }
    }
    segments.resize(s);
    
    // Orient all segments bottom-to-top
    for(Segment& seg : segments) {
        if(yCoordLT(seg.second, seg.first)) {
            std::swap(seg.first, seg.second);
        }
    }
    
    // Create two lists of segments, one ordered by bottom y-coordinate and
    // other by top y-coordinate
    std::vector<Segment> permBottom = segments;
    std::vector<Segment> permTop = std::move(segments);
    
    std::sort(permBottom.begin(), permBottom.end(), [](Segment a, Segment b) {
        return yCoordLT(a.first, b.first);
    });
    std::sort(permTop.begin(), permTop.end(), [](Segment a, Segment b) {
        return yCoordLT(a.second, b.second);
    });
    
    // Left-to-right ordering for segments on a horizontal sweepline that
    // throws an exception if an intersection is found
    struct IntersectionFound { };
    auto orderX = [](Segment a, Segment b) {
        if(a == b) {
            return false;
        }
        if(intersects(a.first, a.second, b.first, b.second)) {
            throw IntersectionFound();
        }
        if(a.first == b.first) {
            return isCCW(a.second, b.first, b.second);
        }
        if(a.second == b.second) {
            return isCCW(a.first, b.first, b.second);
        }
        if(yCoordLT(a.first, b.first)) {
            return isCCW(a.second, a.first, b.first);
        } else {
            return isCCW(a.first, b.first, b.second);
        }
    };
    
    // Run a sweepline from bottom to top, maintaining a ordered set of
    // segments currently on the sweepline
    try {
        std::set<Segment, decltype(orderX)> sweeplineSegments(orderX);
        std::size_t bottomPos = 0;
        std::size_t topPos = 0;
        while(bottomPos != permBottom.size()) {
            // Next event might be a bottom or a top. Break ties in favor of
            // tops so we remove before adding.
            Segment nextBottom = permBottom[bottomPos];
            Segment nextTop = permTop[topPos];
            if(yCoordLT(nextBottom.first, nextTop.second)) {
                if(!sweeplineSegments.insert(nextBottom).second) {
                    // The segment was already there => intersection.
                    throw IntersectionFound();
                }
                ++bottomPos;
            } else {
                sweeplineSegments.erase(nextTop);
                ++topPos;
            }
        }
    } catch(IntersectionFound) {
        return true;
    }
    return false;
}

}
