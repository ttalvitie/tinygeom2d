#include <iostream>

#include "tinygeom2d/domain.hpp"
// BEGIN NOREADME
#include "svg.hpp"
// END NOREADME

using namespace tinygeom2d;

int main() {
    // Create domain with boundary consisting of two polygons
    std::vector<std::vector<Point>> boundary = {
        {{1, 0}, {6, 2}, {10, 1}, {14, 3}, {13, 6}, {7, 7}, {0, 6}}, // Outer polygon
        {{6, 3}, {10, 3}, {7, 6}} // Triangular hole
    };
    Domain domain(boundary);
    
    // Create some points and query whether they are interior points of the domain
    Point a(2, 3);
    Point b(7, 4);
    Point c(13, 2);
    
    std::cout << a << " interior: " << (domain.isInteriorPoint(a) ? "yes" : "no") << "\n";
    std::cout << b << " interior: " << (domain.isInteriorPoint(b) ? "yes" : "no") << "\n";
    std::cout << c << " interior: " << (domain.isInteriorPoint(c) ? "yes" : "no") << "\n";
// BEGIN NOREADME
    
    SVG svg(domain);
    svg.drawPointWithCoords(a, colors::blue);
    svg.drawPointWithCoords(b, colors::red);
    svg.drawPointWithCoords(c, colors::red, true, false);
    svg.drawPointWithCoords({1, 0}, colors::black, false, true);
    svg.drawPointWithCoords({6, 2}, colors::black);
    svg.drawPointWithCoords({10, 1}, colors::black, true, false);
    svg.drawPointWithCoords({14, 3}, colors::black, false, true);
    svg.drawPointWithCoords({13, 6}, colors::black);
    svg.drawPointWithCoords({7, 7}, colors::black);
    svg.drawPointWithCoords({0, 6}, colors::black, true, false);
    svg.drawPointWithCoords({6, 3}, colors::black, false, true);
    svg.drawPointWithCoords({10, 3}, colors::black);
    svg.drawPointWithCoords({7, 6}, colors::black);
    svg.save("domain.svg");
// END NOREADME
}
