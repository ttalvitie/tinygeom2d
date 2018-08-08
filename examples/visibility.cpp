#include <iostream>
#include "tinygeom2d/visibility.hpp"
// BEGIN NOREADME
#include "svg.hpp"
// END NOREADME

using namespace tinygeom2d;

int main() {
    Domain domain({
        {{0, 5}, {10, 0}, {14, 4}, {19, 0}, {29, 5}, {25, 8}, {30, 12}, {2, 14}},
        {{10, 5}, {13, 6}, {11, 7}},
        {{6, 9}, {8, 8}, {9, 10}, {7, 12}},
        {{16, 7}, {21, 2}, {16, 11}},
    });
// BEGIN NOREADME
    SVG svg(domain);
// END NOREADME
    
    // Test visibility of pairs of points
    Point a = {21, 5};
    Point b = {23, 10};
    Point c = {25, 5};
    Point d = {27, 11};
    
    std::cout << a << " sees " << b << ": " << (isVisible(domain, a, b) ? "yes" : "no") << "\n";
    std::cout << c << " sees " << d << ": " << (isVisible(domain, c, d) ? "yes" : "no") << "\n";
    std::cout << "\n";
// BEGIN NOREADME
    
    svg.drawLine(a, b, colors::blue);
    svg.drawLine(c, d, colors::red);
    svg.drawPointWithCoords(a, colors::blue);
    svg.drawPointWithCoords(b, colors::blue, false);
    svg.drawPointWithCoords(c, colors::red);
    svg.drawPointWithCoords(d, colors::red, false);
// END NOREADME
    
    // Compute the visibility of a point
    Point p = {5, 5};
    PointVisibility vis = computePointVisibility(domain, p);
    
    std::cout << "Visible vertices from " << p << ":\n";
    for(Point v : vis.verts) {
        std::cout << "  " << v << "\n";
    }
    std::cout << "\n";
    
    std::cout << "Visibility polygon of " << p << ":\n";
// BEGIN NOREADME
    svg.drawDoublePolygon(vis.computePolygon(), colors::blue, 0.25);
    svg.drawPointWithCoords(p);
    int idx = 0;
// END NOREADME
    for(std::pair<double, double> v : vis.computePolygon()) {
        std::cout << "  (" << v.first << ", " << v.second << ")\n";
// BEGIN NOREADME
        bool right = true;
        bool up = true;
        if(idx == 0) { right = false; up = false; }
        if(idx == 1) { right = false; }
        if(idx == 2) { right = false; }
        if(idx == 3) { right = false; }
        if(idx == 6) { right = false; up = false; }
        if(idx == 7) { right = false; }
        if(idx == 11) { up = false; }
        if(idx == 12) { right = false; }
        if(idx == 13) { right = false; up = false; }
        if(idx == 14) { right = false; up = false; }
        svg.drawDoublePointWithCoords(v.first, v.second, colors::black, right, up);
        ++idx;
// END NOREADME
    }
    
// BEGIN NOREADME
    svg.save("visibility.svg");
// END NOREADME
}
