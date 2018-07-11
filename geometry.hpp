#pragma once

#include "common.hpp"
#include "cmpmul64.hpp"

namespace tinygeom2d {

// The range of allowed point coordinates is [MinCoord, MaxCoord]. Keeping
// coordinates in this range ensures that no overflows happen in the geometry
// functions defined in this file.
const int64_t MinCoord = -((int64_t)1 << 62);
const int64_t MaxCoord = ((int64_t)1 << 62) - 1;

// The type used for points. Points are specified using their x- and
// y-coordinates (for example, x is oriented rightwards and y is oriented
// upwards).
struct Point {
    // Construct point (0, 0).
    Point() : x(0), y(0) { }
    
    // Construct point (x, y). Throws std::out_of_range if one of the
    // coordinates are outside range [MinCoord, MaxCoord].
    Point(int64_t x, int64_t y) : x(x), y(y) {
        if(
            x < MinCoord ||
            x > MaxCoord ||
            y < MinCoord ||
            y > MaxCoord
        ) {
            throw std::out_of_range("tinygeom2d::Point coordinate out of range");
        }
    }
    
    int64_t x;
    int64_t y;
};

inline bool operator==(const Point& a, const Point& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const Point& a, const Point& b) {
    return a.x != b.x || a.y != b.y;
}

// NOTE: The library outside this file does not use the coordinates directly,
// but only using the functions below. This means that you can switch to a
// different point type just by reimplementing the functions in this file.

// NOTE: To simplify the implementation, the geometric operations below use
// a symbolic perturbation where the locations of all points are thought to be
// slightly adjusted such that the following non-degeneracy properties hold:
//   - For all points a, b, c, d such that a != b, c != d, (a, b) != (c, d)
//     and (a, b) != (d, c) it holds that the lines ab and cd are not parallel.
//     (As a special case, no three distinct points are collinear.)
//   - No two distinct points have matching y-coordinates.
// This adjustment is only symbolic in the sense that it only affects these
// degenerate cases. Any reimplementation of this file must also have these
// properties, and the results of the operations must be consistent with each
// other.

// The current symbolic perturbation works as follows:
// Let I be an injective map from points to {0, 1, 2, ...} such that
// I(x, y) < I(u, v) if (x < u || (x == u && y < v)). Now each point (x, y)
// is perturbed to location (x + t^(2 * I(x, y) + 1), y + t^(2 * I(x, y) + 2)),
// and the result of each geometric operation is obtained at the limit t -> 0
// from the positive side. In the end, all the points converge back to their
// original positions, but it can be proven that the non-degeneracy properties
// hold.

// Returns the orientation of the triangle (0, b - a, d - c): 1 if it is
// counterclockwise oriented, -1 if it is clockwise oriented and 0 if the
// points are collinear (or equivalently, two of them are equal).
inline int orientation(Point a, Point b, Point c, Point d) {
    int res = cmpMul64(b.x - a.x, d.y - c.y, b.y - a.y, d.x - c.x);
    if(res) {
        return res;
    }
    
    // Eliminate 0-cases.
    if(
        a == b ||
        c == d ||
        (a == c && b == d) ||
        (a == d && b == c)
    ) {
        return 0;
    }
    
    auto lessInI = [](Point u, Point v) {
        return make_pair(u.x, u.y) < make_pair(v.x, v.y);
    };
    
    // Normalize to the case such that b is the minimum point in I-order.
    int multiplier = 1;
    
    if(lessInI(a, b)) {
        swap(a, b);
        multiplier = -multiplier;
    }
    if(lessInI(c, d)) {
        swap(c, d);
        multiplier = -multiplier;
    }
    if(lessInI(d, b)) {
        swap(a, c);
        swap(b, d);
        multiplier = -multiplier;
    }
    
    // Make sure that b is not in {a, c, d}.
    if(b == d) {
        // We are computing orientation(b, a, c), which is the same as
        // orientation(a, b, c, a).
        d = a;
    }
    
    if(d.y != c.y) {
        return d.y > c.y ? multiplier : -multiplier;
    } else {
        return d.x < c.x ? multiplier : -multiplier;
    }
}

// Returns 1, 0 or -1 if the y-coordinate of a is greater than, equal to or
// less than the y-coordinate of b, respectively.
inline int cmpY(Point a, Point b) {
    if(a.y == b.y) {
        if(a.x == b.x) {
            return 0;
        } else {
            return a.x < b.x ? 1 : -1;
        }
    } else {
        return a.y > b.y ? 1 : -1;
    }
}

// Returns 1, 0 or -1 if b - a has angle greater than, equal to or less than
// d - c, respectively when angles are defined as standard directional angles:
//   0 degrees rightwards/positive x
//   90 degrees upwards/positive y
//   180 degrees leftwards/negative x
//   270 degrees downwards/negative y.
// The result when a = b or c = d is undefined.
inline int cmpAngle(Point a, Point b, Point c, Point d) {
    int c1 = cmpY(a, b);
    int c2 = cmpY(c, d);
    
    if(c1 == c2) {
        return orientation(c, d, a, b);
    } else {
        return c1 > c2 ? 1 : -1;
    }
}

}

// Make it possible to use unordered_set<Point>.
namespace std {
template <>
struct hash<tinygeom2d::Point> {
    size_t operator()(const tinygeom2d::Point& p) const {
        uint64_t hash = 0;
        auto add = [&](uint64_t elem) {
            const uint64_t m = UINT64_C(0xc6a4a7935bd1e995);
            elem *= m;
            elem ^= elem >> 47;
            elem *= m;
            hash ^= elem;
            hash *= m;
        };
        add(p.x);
        add(p.y);
        return (size_t)hash;
    }
};
}
