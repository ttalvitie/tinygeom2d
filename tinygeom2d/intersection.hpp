#pragma once

#include <algorithm>
#include <cstddef>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "geometry.hpp"

namespace tinygeom2d {

namespace intersection_detail {

// Returns true if segment a is to the left of segment b in the left-to-right
// ordering for segments with intersecting y-coordinate ranges, such that the
// comparison is done at the bottommost horizontal line that intersects both
// segments. If both segments have the same bottom point, the comparison is done
// at a higher horizontal line that is infinitesimally close. The points in both
// segments must be ordered such that the bottom point is first. The segments
// must not have length zero.
inline bool segmentLeftOfAtBottom(
    std::pair<Point, Point> a,
    std::pair<Point, Point> b
) {
    if(a == b) {
        return false;
    }
    if(a.first == b.first) {
        return isCCW(a.second, b.first, b.second);
    }
    if(yCoordLT(a.first, b.first)) {
        return isCCW(a.second, a.first, b.first);
    } else {
        return isCCW(a.first, b.first, b.second);
    }
}
    
}

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

// Returns true if two segments among given segments intersect as defined for
// intersects(a, b, c, d).
// Time complexity: O(n log n), where n = input size.
inline bool intersects(std::vector<std::pair<Point, Point>> segments) {
    typedef std::pair<Point, Point> Segment;
    
    // Remove zero-length segments
    std::size_t s = segments.size();
    for(std::size_t i = 0; i < s; ++i) {
        if(segments[i].first == segments[i].second) {
            segments[i] = segments[--s];
            --i;
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
    
    // Run a sweepline from bottom to top, maintaining a ordered set of
    // segments currently on the sweepline. This way, intersecting segments
    // can be caught as neighboring segments.
    typedef std::set<Segment, bool (*)(Segment, Segment)> Sweepline;
    typedef Sweepline::iterator Iter;
    Sweepline sweepline(intersection_detail::segmentLeftOfAtBottom);
    std::size_t bottomPos = 0;
    std::size_t topPos = 0;
    while(topPos != permTop.size()) {
        // Next event might be a bottom or a top. Break ties in favor of
        // tops so we remove before adding.
        if(
            bottomPos != permBottom.size() &&
            yCoordLT(permBottom[bottomPos].first, permTop[topPos].second)
        ) {
            Iter pos;
            bool inserted;
            std::tie(pos, inserted) = sweepline.insert(permBottom[bottomPos++]);
            
            // Duplicate segment means intersection
            if(!inserted) {
                return true;
            }
            
            // Check for intersection with neighbors
            if(pos != sweepline.begin()) {
                Iter prev = pos;
                --prev;
                if(intersects(
                    prev->first, prev->second,
                    pos->first, pos->second
                )) {
                    return true;
                }
            }
            
            Iter next = pos;
            ++next;
            if(
                next != sweepline.end() &&
                intersects(
                    next->first, next->second,
                    pos->first, pos->second
                )
            ) {
                return true;
            }
        } else {
            Iter pos = sweepline.find(permTop[topPos++]);
            
            // Check for intersecting neighbors introduced by removal
            if(pos != sweepline.begin()) {
                Iter prev = pos;
                --prev;
                Iter next = pos;
                ++next;
                if(next != sweepline.end()) {
                    if(intersects(
                        prev->first, prev->second,
                        next->first, next->second)
                    ) {
                        return true;
                    }
                }
            }
            
            sweepline.erase(pos);
        }
    }
    return false;
}

}
