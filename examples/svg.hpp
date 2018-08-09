#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "tinygeom2d/domain.hpp"
#include "tinygeom2d/geometry.hpp"

namespace colors {
    const std::string black = "black";
    const std::string white = "white";
    const std::string lightgray = "gainsboro";
    const std::string red = "firebrick";
    const std::string blue = "mediumblue";
}

class SVG {
public:
    typedef tinygeom2d::Domain Domain;
    typedef tinygeom2d::Point Point;
    
    SVG(const Domain& domain) {
        ss << "<svg width=\"" << Width << "\" height=\"" << Height << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        ss << "  <rect x=\"0\" y=\"0\" width=\"" << Width << "\" height=\"" << Height << "\" fill=\"" << colors::lightgray << "\" />\n";
        
        const double Infinity = std::numeric_limits<double>::infinity();
        double minX = Infinity;
        double maxX = -Infinity;
        double minY = Infinity;
        double maxY = -Infinity;
        for(const std::vector<Point>& poly : domain.boundary()) {
            for(Point p : poly) {
                minX = std::min(minX, (double)p.x);
                maxX = std::max(maxX, (double)p.x);
                minY = std::min(minY, (double)p.y);
                maxY = std::max(maxY, (double)p.y);
            }
        }
        
        double adjX = 0.08 * (maxY - minY);
        double adjY = 0.08 * (maxY - minY);
        minX -= adjX;
        maxX += adjX;
        minY -= adjY;
        maxY += adjY;
        
        double aspect = Height / Width;
        if((maxY - minY) / (maxX - minX) > aspect) {
            double midX = 0.5 * (minX + maxX);
            double d = 0.5 * (maxY - minY) / aspect;
            minX = midX - d;
            maxX = midX + d;
        } else {
            double midY = 0.5 * (minY + maxY);
            double d = 0.5 * (maxX - minX) * aspect;
            minY = midY - d;
            maxY = midY + d;
        }
        
        sx = minX;
        dx = Width / (maxX - minX);
        sy = minY;
        dy = Height / (maxY - minY);
        
        ss << "  <path d=\"";
        for(const std::vector<Point>& poly : domain.boundary()) {
            bool first = true;
            for(Point p : poly) {
                ss << (first ? 'M' : 'L') << ' ';
                first = false;
                ss << mapX(p) << ' ' << mapY(p) << ' ';
            }
            ss << "Z ";
        }
        ss << "\" fill=\"" << colors::white << "\" stroke=\"black\" stroke-width=\"3\" />\n";
    }
    
    void drawDoublePointWithCoords(
        double px,
        double py,
        const std::string& color = colors::black,
        bool right = true,
        bool up = true
    ) {
        double x = mapX(px);
        double y = mapY(py);
        ss << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"7\" fill=\"" << color << "\" />\n";
        ss << "  <text x=\"" << (right ? x + 8 : x - 8) << "\" y=\"" << (up ? y - 5 : y + 15) << "\" style=\"font: normal 16px sans-serif\" fill=\"" << color << "\"" << (right ? "" : " text-anchor=\"end\"") << ">(" << coordToString(px) << ", " << coordToString(py) << ")</text>\n";
    }
    void drawPointWithCoords(
        Point p,
        const std::string& color = colors::black,
        bool right = true,
        bool up = true
    ) {
        std::pair<double, double> c = coordsAsDouble(p);
        drawDoublePointWithCoords(c.first, c.second, color, right, up);
    }
    
    void drawLine(Point a, Point b, const std::string& color = colors::black) {
        double x1 = mapX(a);
        double y1 = mapY(a);
        double x2 = mapX(b);
        double y2 = mapY(b);
        ss << "  <line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\" stroke=\"" << color << "\" stroke-width=\"4\" />\n";
    }
    void drawDoublePolygon(
        std::vector<std::pair<double, double>> poly,
        const std::string& color = colors::black,
        double opacity = 1.0
    ) {
        ss << "  <path d=\"";
        bool first = true;
        for(std::pair<double, double> p : poly) {
            ss << (first ? 'M' : 'L') << ' ';
            first = false;
            ss << mapX(p.first) << ' ' << mapY(p.second) << ' ';
        }
        ss << "Z ";
        ss << "\" fill=\"" << color << "\" opacity=\"" << opacity << "\" />\n";
    }
    
    void save(const std::string& filename) {
        std::ofstream fp;
        fp.exceptions(fp.failbit | fp.badbit);
        fp.open(filename);
        fp << ss.str() << "</svg>\n";
        fp.close();
    }
    
private:
    static constexpr double Width = 800;
    static constexpr double Height = 400;
    
    double sx;
    double sy;
    double dx;
    double dy;
    
    std::stringstream ss;
    
    double mapX(double x) {
        return dx * (x - sx);
    }
    double mapY(double y) {
        return Height - dy * (y - sy);
    }
    double mapX(Point p) {
        return mapX((double)p.x);
    }
    double mapY(Point p) {
        return mapY((double)p.y);
    }
    
    std::string coordToString(double c) {
        std::stringstream ret;
        int ci = (int)std::round(10.0 * c);
        if(ci % 10 == 0) {
            ret << ci / 10;
        } else {
            ret << ci / 10 << "." << ci % 10;
        }
        return ret.str();
    }
};
