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
    Point(int64_t x, int64_t y) : x(x), y(y) { }
    
    int64_t x;
    int64_t y;
};

// NOTE: The library outside this file does not use the coordinates directly,
// but only using the functions below. This means that you can switch to a
// different point type just by reimplementing the functions in this file.

// NOTE: To simplify the implementation, the geometric operations below use
// a symbolic perturbation where the locations of all points are slightly
// adjusted such that the following non-degeneracy properties hold:
//   - No three distinct points are collinear
//   - No two distinct points have matching x-coordinates
//   - No two distinct points have matching y-coordinates
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
// original positions, but it can be proven that no three points are collinear
// anymore, and all the x- and y-coordinates are different.

// Returns the orientation of the triangle (a, b, c): 1 if it is
// counterclockwise oriented, -1 if it is clockwise oriented and 0 if the
// points are collinear (or equivalently, two of them are equal).
int orientation(Point a, Point b, Point c) {
    int res = cmpMul64(b.x - a.x, c.y - a.y, b.y - a.y, c.x - a.x);
    if(res) {
        return res;
    }
    
    // Eliminate the case where we have two equal points.
    if(
        (a.x == b.x && a.y == b.y) ||
        (a.x == c.x && a.y == c.y) ||
        (b.x == c.x && b.y == c.y)
    ) {
        return 0;
    }
    
    auto lessInI = [](Point u, Point v) {
        return make_pair(u.x, u.y) < make_pair(v.x, v.y);
    };
    
    // Rotate points such that B is the minimum in I-order.
    if(lessInI(a, c)) {
        if(lessInI(a, b)) {
            tie(a, b, c) = make_tuple(c, a, b);
        }
    } else {
        if(lessInI(c, b)) {
            tie(a, b, c) = make_tuple(b, c, a);
        }
    }
    
    if(c.y != a.y) {
        return c.y > a.y ? 1 : -1;
    } else {
        return c.x < a.x ? 1 : -1;
    }
}

// Returns 1, 0 or -1 if the x-coordinate of a is greater than, equal to or
// less than the x-coordinate of b, respectively.
int cmpX(Point a, Point b) {
    if(a.x == b.x) {
        if(a.y == b.y) {
            return 0;
        } else {
            return a.y < b.y ? 1 : -1;
        }
    } else {
        return a.x > b.x ? 1 : -1;
    }
}

// Returns 1, 0 or -1 if the y-coordinate of a is greater than, equal to or
// less than the y-coordinate of b, respectively.
int cmpY(Point a, Point b) {
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

}
