#include <iostream>

#include "tinygeom2d/domain.hpp"
#include "svg.hpp"

using namespace tinygeom2d;

int main() {
    std::vector<std::vector<Point>> boundary = {
        {{2, 0}, {5, 2}, {8, 1}, {12, 5}, {6, 7}, {0, 6}}, // Outer polygon
        {{5, 3}, {8, 3}, {6, 6}} // Triangular hole
    };
    Domain domain(boundary);
    
    Point a(2, 3);
    Point b(5, 1);
    Point c(6, 4);
    
    std::cout << a << " interior: " << (domain.isInteriorPoint(a) ? "yes" : "no") << "\n";
    std::cout << b << " interior: " << (domain.isInteriorPoint(b) ? "yes" : "no") << "\n";
    std::cout << c << " interior: " << (domain.isInteriorPoint(c) ? "yes" : "no") << "\n";
    
    SVG svg(domain);
    svg.drawPointWithCoords(a);
    svg.drawPointWithCoords(b);
    svg.drawPointWithCoords(c);
    svg.save("domain.svg");
}
