#include <iostream>
#include "tinygeom2d/shortestpath.hpp"
// BEGIN NOREADME
#include "svg.hpp"
SVG svg(tinygeom2d::Domain{});
// END NOREADME

using namespace tinygeom2d;

void showShortestPath(const ShortestPathContext& ctx, Point a, Point b) {
    // ShortestPathContext::findShortestPath returns pair (path, path length).
    // If the path is empty, there is no path between given points.
    std::pair<std::vector<Point>, double> result = ctx.findShortestPath(a, b);
    
    if(result.first.empty()) {
        std::cout << "No path between " << a << " and " << b << "\n";
// BEGIN NOREADME
        svg.drawPointWithCoords(a, colors::red);
        svg.drawPointWithCoords(b, colors::red);
// END NOREADME
    } else {
        std::cout << "Shortest path between " << a << " and " << b;
        std::cout << " of length " << result.second << ":\n";
        for(Point p : result.first) {
            std::cout << "  " << p << "\n";
        }
// BEGIN NOREADME
        bool first = true;
        Point prev;
        for(Point p : result.first) {
            bool right = true;
            bool up = true;
            if(p == Point(3, 8)) {
                right = false;
            }
            if(p == Point(9, 3)) {
                up = false;
            }
            if(p == Point(7, 12)) {
                right = false;
            }
            if(p == Point(21, 4)) {
                right = false;
                up = false;
            }
            if(p == Point(24, 4)) {
                up = false;
            }
            svg.drawPointWithCoords(p, colors::blue, right, up);
            if(!first) {
                svg.drawLine(prev, p, colors::blue);
            }
            first = false;
            prev = p;
        }
// END NOREADME
    }
    std::cout << "\n";
}

int main() {
    Domain domain({
        {{0, 0}, {10, 0}, {15, 8}, {20, 0}, {30, 0}, {30, 14}, {0, 14}},
        {{7, 12}, {11, 12}, {9, 3}},
        {{21, 4}, {24, 4}, {26, 6}, {26, 9}, {24, 11}, {21, 11}, {19, 9}, {19, 6}},
        {{12, 0}, {15, 5}, {18, 0}}
    });
// BEGIN NOREADME
    svg = SVG(domain);
// END NOREADME
    
    // Create context for computing shortest paths in domain
    ShortestPathContext ctx(domain);
    
    // Find shortest paths connecting pairs of points
    showShortestPath(ctx, {5, 5}, {27, 5});
    showShortestPath(ctx, {3, 8}, {27, 7});
    showShortestPath(ctx, {14, 1}, {22, 1});
// BEGIN NOREADME
    svg.save("shortestpath.svg");
// END NOREADME
}
