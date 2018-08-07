#pragma once

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

#include "tinygeom2d/domain.hpp"
#include "tinygeom2d/geometry.hpp"

class SVG {
public:
    typedef tinygeom2d::Domain Domain;
    typedef tinygeom2d::Point Point;
    
    SVG(const Domain& domain) {
        ss << "<svg width=\"" << Width << "\" height=\"" << Height << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
        ss << "  <rect x=\"0\" y=\"0\" width=\"" << Width << "\" height=\"" << Height << "\" fill=\"silver\" />\n";
        
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
        
        double adjX = 0.05 * (maxX - minX);
        double adjY = 0.05 * (maxY - minY);
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
        ss << "\" fill=\"white\" stroke=\"black\" stroke-width=\"3\" />\n";
    }
    
    void drawPointWithCoords(Point p) {
        double x = mapX(p);
        double y = mapY(p);
        ss << "  <circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"7\" fill=\"red\" />\n";
        ss << "  <text x=\"" << x + 8 << "\" y=\"" << y - 5 << "\" style=\"font: normal 20px sans-serif\" fill=\"red\">(" << p.x << ", " << p.y << ")</text>\n";
    }
    
    void save(const std::string& filename) {
        std::ofstream fp;
        fp.exceptions(fp.failbit | fp.badbit);
        fp.open(filename);
        fp << ss.str() << "</svg>\n";
        fp.close();
    }
    
private:
    const double Width = 800;
    const double Height = 512;
    
    double sx;
    double sy;
    double dx;
    double dy;
    
    std::stringstream ss;
    
    double mapX(Point p) {
        return dx * ((double)p.x - sx);
    }
    double mapY(Point p) {
        return Height - dy * ((double)p.y - sy);
    }
};
