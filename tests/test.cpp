#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "tinygeom2d/int64.hpp"
#include "tinygeom2d/geometry.hpp"
#include "tinygeom2d/intersection.hpp"
#include "tinygeom2d/domain.hpp"
#include "tinygeom2d/visibility.hpp"
#include "tinygeom2d/shortestpath.hpp"

using namespace tinygeom2d;

std::mt19937 rng(std::random_device{}());

void test(bool ok, const char* str, const char* file, int line) {
    if(!ok) {
        std::cerr << "Test failed: TEST(" << str << ") @ " << file << ":" << line << std::endl;
        abort();
    }
}

#define TEST(check) test((check), #check, __FILE__, __LINE__) 

void test_int64_hpp() {
    using namespace int64_detail;
    
    // portableMulU64: random tests generated using Python 3:
    // from random import randint
    // x = randint(0, 2**64-1)
    // y = randint(0, 2**64-1)
    // print("TEST(portableMulU64(UINT64_C({}), UINT64_C({})) == std::make_pair(UINT64_C({}), UINT64_C({})));".format(x, y, (x * y) >> 64, (x * y) & ((1 << 64) - 1)))
    TEST(portableMulU64(UINT64_C(13850607296698408471), UINT64_C(1962650796196979807)) == std::make_pair(UINT64_C(1473642466662696899), UINT64_C(4227264815221106313)));
    TEST(portableMulU64(UINT64_C(10116009595631013116), UINT64_C(13460533840250291104)) == std::make_pair(UINT64_C(7381621870298187970), UINT64_C(6158501803456860544)));
    TEST(portableMulU64(UINT64_C(2440460725643442406), UINT64_C(11966028096025722649)) == std::make_pair(UINT64_C(1583077289607794060), UINT64_C(2994133882693052534)));
    TEST(portableMulU64(UINT64_C(17507301181100150215), UINT64_C(12572838332720158429)) == std::make_pair(UINT64_C(11932537607323594689), UINT64_C(4098711706989444811)));
    TEST(portableMulU64(UINT64_C(1368121625145395075), UINT64_C(18420356669026126822)) == std::make_pair(UINT64_C(1366164576311486871), UINT64_C(12310193224563368114)));
    TEST(portableMulU64(UINT64_C(10089565848189686975), UINT64_C(6093245880170678786)) == std::make_pair(UINT64_C(3332740200196728737), UINT64_C(11143210876361023358)));
    TEST(portableMulU64(UINT64_C(13146144872988902367), UINT64_C(11934585562814826332)) == std::make_pair(UINT64_C(8505229442167618003), UINT64_C(17714820330509384996)));
    TEST(portableMulU64(UINT64_C(10578583212407853803), UINT64_C(8521760899801472587)) == std::make_pair(UINT64_C(4886941372123950292), UINT64_C(14771716385615926489)));
    TEST(portableMulU64(UINT64_C(12022468742568605552), UINT64_C(11857499422431886743)) == std::make_pair(UINT64_C(7727998805728822665), UINT64_C(14325480725176820496)));
    TEST(portableMulU64(UINT64_C(11199764207885964965), UINT64_C(5255569581116580295)) == std::make_pair(UINT64_C(3190868797845636132), UINT64_C(3330522714970775363)));
    
    // portableCmpMul64: random tests generated using Python 3:
    // from random import randint
    // a = randint(-2**63, 2**63-1)
    // b = randint(-2**63, 2**63-1)
    // c = randint(-2**63, 2**63-1)
    // d = randint(-2**63, 2**63-1)
    // x = a * b
    // y = c * d
    // print("TEST(portableCmpMul64(INT64_C({}), INT64_C({}), INT64_C({}), INT64_C({})) == {});".format(a, b, c, d, int(x > y) - int(x < y)))
    TEST(portableCmpMul64(INT64_C(630025097722727616), INT64_C(4393248892236142053), INT64_C(-7724918783757161038), INT64_C(9117988473785604542)) == 1);
    TEST(portableCmpMul64(INT64_C(4627172242198388511), INT64_C(-2027138965779481110), INT64_C(-9153785664229383044), INT64_C(-5773843628854979992)) == -1);
    TEST(portableCmpMul64(INT64_C(-2935210532528881641), INT64_C(3475005425608740518), INT64_C(-4527987550838865099), INT64_C(-5931763926304303230)) == -1);
    TEST(portableCmpMul64(INT64_C(-3128268382030838515), INT64_C(-7592354795491578321), INT64_C(3773617965958904167), INT64_C(-4085160704151535788)) == 1);
    TEST(portableCmpMul64(INT64_C(-1192843127667658217), INT64_C(-744989233867644813), INT64_C(1324644132375749105), INT64_C(2592987881356638514)) == -1);
    TEST(portableCmpMul64(INT64_C(-4995491892254007254), INT64_C(-6961680915438065709), INT64_C(-2366176534911095068), INT64_C(-4912420945505275077)) == 1);
    TEST(portableCmpMul64(INT64_C(-3864321286108431747), INT64_C(7585441108130198550), INT64_C(1354054881586531315), INT64_C(-4081763376108833670)) == -1);
    TEST(portableCmpMul64(INT64_C(-8377002427053400651), INT64_C(4293282623089319352), INT64_C(6795134510553725513), INT64_C(-5427205919954543792)) == 1);
    TEST(portableCmpMul64(INT64_C(-1188104836405824677), INT64_C(8555573039862816606), INT64_C(-7932982476440084043), INT64_C(759848102779321690)) == -1);
    TEST(portableCmpMul64(INT64_C(-7282865720591720171), INT64_C(-38451447443525431), INT64_C(2183924849438275371), INT64_C(912077392808913582)) == -1);
    
    // cmpMul64: crosscheck with portableCmpMul64 in random tests
    for(int t = 0; t < 100; ++t) {
        int64_t a = std::uniform_int_distribution<int64_t>(INT64_MIN, INT64_MAX)(rng);
        int64_t b = std::uniform_int_distribution<int64_t>(INT64_MIN, INT64_MAX)(rng);
        int64_t c = std::uniform_int_distribution<int64_t>(INT64_MIN, INT64_MAX)(rng);
        int64_t d = std::uniform_int_distribution<int64_t>(INT64_MIN, INT64_MAX)(rng);
        TEST(cmpMul64(a, b, c, d) == portableCmpMul64(a, b, c, d));
    }
    
    // portableCmpMul64 and cmpMul64: equality cases
    for(int t = 0; t < 100; ++t) {
        int64_t a = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
        int64_t b = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
        int64_t c = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
        int64_t d = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
        TEST(portableCmpMul64(a * b, c * d, a * c, b * d) == 0);
        TEST(cmpMul64(a * b, c * d, a * c, b * d) == 0);
    }
    
    // portableCmpMul64 and cmpMul64: almost-equality cases
    for(int t = 0; t < 100; ++t) {
        int64_t a = 0;
        int64_t b = 0;
        int64_t c = 0;
        int64_t d = 0;
        while(!a || !b || !c || !d) {
            a = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
            b = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
            c = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
            d = std::uniform_int_distribution<int64_t>(INT32_MIN, INT32_MAX)(rng);
        }
        int64_t ab = a * b;
        int64_t cd = c * d;
        int64_t ac = a * c;
        int64_t bd = b * d;
        int64_t abSign = (ab > 0) - (ab < 0);
        int64_t cdSign = (cd > 0) - (cd < 0);
        int64_t acSign = (ac > 0) - (ac < 0);
        int64_t bdSign = (bd > 0) - (bd < 0);
        TEST(portableCmpMul64(ab + cdSign, cd, ac, bd) == 1);
        TEST(cmpMul64(ab + cdSign, cd, ac, bd) == 1);
        TEST(portableCmpMul64(ab, cd + abSign, ac, bd) == 1);
        TEST(cmpMul64(ab, cd + abSign, ac, bd) == 1);
        TEST(portableCmpMul64(ab, cd, ac - bdSign, bd) == 1);
        TEST(cmpMul64(ab, cd, ac - bdSign, bd) == 1);
        TEST(portableCmpMul64(ab, cd, ac, bd - acSign) == 1);
        TEST(cmpMul64(ab, cd, ac, bd - acSign) == 1);
        TEST(portableCmpMul64(ab - cdSign, cd, ac, bd) == -1);
        TEST(cmpMul64(ab - cdSign, cd, ac, bd) == -1);
        TEST(portableCmpMul64(ab, cd - abSign, ac, bd) == -1);
        TEST(cmpMul64(ab, cd - abSign, ac, bd) == -1);
        TEST(portableCmpMul64(ab, cd, ac + bdSign, bd) == -1);
        TEST(cmpMul64(ab, cd, ac + bdSign, bd) == -1);
        TEST(portableCmpMul64(ab, cd, ac, bd + acSign) == -1);
        TEST(cmpMul64(ab, cd, ac, bd + acSign) == -1);
    }
}

Point randomPoint() {
    std::uniform_int_distribution<int64_t> dist(MinCoord, MaxCoord);
    return Point(dist(rng), dist(rng));
}

void test_geometry_hpp() {
    // isCCW: check that sign is correct
    TEST(isCCW({0, 0}, {1, 0}, {0, 0}, {0, 1}));
    
    // isCCW: random inversion and rotation test
    for(int i = 0; i < 100; ++i) {
        Point a, b, c, d;
        while(a == b || a == c || a == d || b == c || b == d || c == d) {
            a = randomPoint();
            b = randomPoint();
            c = randomPoint();
            d = randomPoint();
        }
        
        bool o = isCCW(a, b, c, d);
        TEST(isCCW(b, a, c, d) == !o);
        TEST(isCCW(a, b, d, c) == !o);
        TEST(isCCW(b, a, d, c) == o);
        TEST(isCCW(c, d, a, b) == !o);
        TEST(isCCW(d, c, a, b) == o);
        TEST(isCCW(c, d, b, a) == o);
        TEST(isCCW(d, c, b, a) == !o);
    }
    
    // isCCW: random degenerate cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        Point c = randomPoint();
        
        TEST(!isCCW(a, a, b, c));
        TEST(!isCCW(b, c, a, a));
        TEST(!isCCW(a, b, a, b));
        TEST(!isCCW(a, b, b, a));
    }
    
    // isCCW: inversion and rotation test in small coordinates for
    // three points
    {
        for(int bx = -10; bx <= 10; ++bx) {
        for(int by = -10; by <= 10; ++by) {
        for(int cx = -10; cx <= 10; ++cx) {
        for(int cy = -10; cy <= 10; ++cy) {
            Point a(0, 0);
            Point b(bx, by);
            Point c(cx, cy);
            if(a == b || a == c || b == c) {
                continue;
            }
            
            bool o = isCCW(a, b, a, c);
            TEST(isCCW(b, c, b, a) == o);
            TEST(isCCW(c, a, c, b) == o);
            TEST(isCCW(a, c, a, b) == !o);
            TEST(isCCW(c, b, c, a) == !o);
            TEST(isCCW(b, a, b, c) == !o);
        }
        }
        }
        }
    }
    
    // isCCW: random three-point degenerate cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        TEST(!isCCW(a, a, a, b));
        TEST(!isCCW(a, b, a, a));
        TEST(!isCCW(b, a, b, a));
        TEST(!isCCW(a, a, a, a));
    }
    
    // isCCW: three-parameter version comparison in random cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        Point c = randomPoint();
        TEST(isCCW(a, b, c) == isCCW(a, b, a, c));
    }
    
    // isCCW: three-parameter version comparison in small triangles
    {
        for(int bx = -10; bx <= 10; ++bx) {
        for(int by = -10; by <= 10; ++by) {
        for(int cx = -10; cx <= 10; ++cx) {
        for(int cy = -10; cy <= 10; ++cy) {
            Point a(0, 0);
            Point b(bx, by);
            Point c(cx, cy);
            TEST(isCCW(a, b, c) == isCCW(a, b, a, c));
        }
        }
        }
        }
    }
    
    // yCoordLT: random degenerate cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        TEST(yCoordLT(a, a) == false);
    }

    // yCoordLT: random tests
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        bool o = yCoordLT(a, b);
        TEST(a == b || yCoordLT(b, a) == !o);
        if(a.y != b.y) {
            TEST(o == (a.y < b.y));
        }
    }

    // yCoordLT: total ordering for small coordinates
    {
        std::vector<Point> points;
        for(int x = -10; x <= 10; ++x) {
            for(int y = -10; y <= 10; ++y) {
                points.emplace_back(x, y);
            }
        }
        
        sort(points.begin(), points.end(), [&](Point a, Point b) {
            return yCoordLT(a, b);
        });
        
        for(int i = 0; i < (int)points.size(); ++i) {
            for(int j = 0; j < (int)points.size(); ++j) {
                TEST(yCoordLT(points[i], points[j]) == (i < j));
            }
        }
    }
    
    // angleLT: check that sign is correct
    {
        Point o;
        Point a(1, 1);
        Point b(-1, 1);
        Point c(-1, -1);
        Point d(1, -1);
        TEST(angleLT(o, a, o, b));
        TEST(angleLT(o, b, o, c));
        TEST(angleLT(o, c, o, d));
        TEST(!angleLT(o, d, o, a));
    }
    
    // angleLT: total ordering for small coordinates
    {
        Point origin(0, 0);
        
        std::vector<Point> points;
        for(int x = -10; x <= 10; ++x) {
            for(int y = -10; y <= 10; ++y) {
                if(x || y) {
                    points.emplace_back(x, y);
                }
            }
        }
        
        sort(points.begin(), points.end(), [&](Point a, Point b) {
            return angleLT(origin, a, origin, b);
        });
        
        for(int i = 0; i < (int)points.size(); ++i) {
            for(int j = 0; j < (int)points.size(); ++j) {
                TEST(angleLT(origin, points[i], origin, points[j]) == (i < j));
            }
        }
    }
    
    // normalizationFactor works for random ranges
    for(int t = 0; t < 50; ++t) {
        double minX = std::normal_distribution<double>(0.0, 1.0)(rng);
        double maxX = std::normal_distribution<double>(0.0, 1.0)(rng);
        double minY = std::normal_distribution<double>(0.0, 1.0)(rng);
        double maxY = std::normal_distribution<double>(0.0, 1.0)(rng);
        if(minX > maxX) {
            std::swap(minX, maxX);
        }
        if(minY > maxY) {
            std::swap(minY, maxY);
        }
        double factor = normalizationFactor(minX, maxX, minY, maxY);
        TEST(std::isfinite(factor));
        TEST(factor >= 0.0);
        for(int t2 = 0; t2 < 50; ++t2) {
            double x = std::uniform_real_distribution<double>(minX, maxX)(rng);
            double y = std::uniform_real_distribution<double>(minY, maxY)(rng);
            Point point(std::round(x * factor), std::round(y * factor));
        }
    }
    
    // normalizationFactor works for zero range
    {
        double factor = normalizationFactor(0.0, 0.0, 0.0, 0.0);
        TEST(std::isfinite(factor));
        TEST(factor >= 0.0);
    }
    {
        double factor = normalizationFactor(-0.0, -0.0, -0.0, -0.0);
        TEST(std::isfinite(factor));
        TEST(factor >= 0.0);
    }
    
    // Hash and comparison operators of Point
    {
        std::unordered_set<Point> x;
        x.emplace(5, 6);
        x.emplace(-5, 6);
        x.emplace(1, 3);
        x.emplace(1, 3);
        TEST(x.size() == 3);
        TEST(x.count(Point(5, 6)));
        TEST(!x.count(Point(6, 5)));
    }
}

void test_intersection_hpp() {
    // Example intersections
    {
        TEST(intersects({-1, 0}, {1, 0}, {0, -1}, {0, 1}));
        TEST(intersects({1, 0}, {-1, 0}, {0, -1}, {0, 1}));
        TEST(intersects({0, 0}, {0, 1}, {0, 0}, {0, 1}));
        TEST(intersects({0, 0}, {0, 1}, {0, 1}, {0, 0}));
    }
    
    // Example non-intersections
    {
        TEST(!intersects({-1, 0}, {0, 0}, {1, 0}, {1, 1}));
        TEST(!intersects({-2, 0}, {-1, 0}, {2, 0}, {1, 0}));
        TEST(!intersects({-2, 0}, {-1, 0}, {1, 0}, {2, 0}));
        TEST(!intersects({0, 0}, {0, -1}, {0, 0}, {0, 1}));
        TEST(!intersects({0, 0}, {0, 1}, {1, 0}, {1, 1}));
        TEST(!intersects({0, 0}, {0, 1}, {0, 0}, {1, 0}));
        TEST(!intersects({0, 0}, {0, 0}, {0, 0}, {0, 1}));
        TEST(!intersects({0, 0}, {0, 0}, {0, 0}, {0, 0}));
    }
    
    // Perturbation consistency check
    TEST(
        intersects({0, 0}, {-1, 0}, {0, -1}, {0, 1}) !=
        intersects({0, 0}, {1, 0}, {0, -1}, {0, 1})
    );
    
    // Both versions agree for small coordinates
    {
        std::vector<std::pair<Point, Point>> segments;
        segments.resize(2);
        for(int bx = -3; bx <= 3; ++bx) {
        for(int by = -3; by <= 3; ++by) {
        for(int cx = -3; cx <= 3; ++cx) {
        for(int cy = -3; cy <= 3; ++cy) {
        for(int dx = -3; dx <= 3; ++dx) {
        for(int dy = -3; dy <= 3; ++dy) {
            Point a(0, 0);
            Point b(bx, by);
            Point c(cx, cy);
            Point d(dx, dy);
            segments[0] = {a, b};
            segments[1] = {c, d};
            TEST(intersects(a, b, c, d) == intersects(segments));
        }
        }
        }
        }
        }
        }
    }
    
    // Both versions agree for random points
    for(int t = 0; t < 100; ++t) {
        Point a = randomPoint();
        Point b = randomPoint();
        Point c = randomPoint();
        Point d = randomPoint();
        std::vector<std::pair<Point, Point>> segments = {{a, b}, {c, d}};
        TEST(intersects(a, b, c, d) == intersects(segments));
    }
    
    // Random-generate non-intersecting configuration with brute-force use of
    // pairwise intersects, and verify the result of the vector version of
    // intersects.
    for(int t = 0; t < 1000; ++t) {
        std::vector<std::pair<Point, Point>> segments; 
        const std::size_t segmentCount = 10;
        while(segments.size() != segmentCount) {
            Point a = randomPoint();
            Point b = randomPoint();
            bool ok = true;
            for(std::pair<Point, Point> segment : segments) {
                if(intersects(a, b, segment.first, segment.second)) {
                    ok = false;
                    break;
                }
            }
            if(ok) {
                segments.emplace_back(a, b);
            }
        }
        std::shuffle(segments.begin(), segments.end(), rng);
        TEST(!intersects(segments));
    }
    
    // Random-generate intersecting configuration with brute-force use of
    // pairwise intersects, and verify the result of the vector version of
    // intersects.
    for(int t = 0; t < 1000; ++t) {
        std::vector<std::pair<Point, Point>> segments;
        const std::size_t segmentCount = 10;
        while(segments.size() != segmentCount) {
            Point a = randomPoint();
            Point b = randomPoint();
            bool ok = true;
            for(std::pair<Point, Point> segment : segments) {
                if(intersects(a, b, segment.first, segment.second)) {
                    ok = false;
                    break;
                }
            }
            if(segments.size() == segmentCount - 1) {
                ok = !ok;
            }
            if(ok) {
                segments.emplace_back(a, b);
            }
        }
        std::shuffle(segments.begin(), segments.end(), rng);
        if(!intersects(segments)) {
            for(auto seg : segments) {
                std::cerr << seg.first << " - " << seg.second << "\n";
            }
        }
        TEST(intersects(segments));
    }
}

// Example correctly oriented bounded domain boundary
const std::vector<std::vector<Point>> exampleBoundary = {
    {{0, 4}, {1, 1}, {3, 2}, {6, 0}, {9, 3}, {8, 8}},
    {{2, 4}, {7, 6}, {8, 3}, {5, 1}, {2, 3}},
    {{5, 2}, {6, 2}, {5, 5}, {3, 3}},
    {{4, 3}, {5, 4}, {5, 3}},
    {{7, 3}, {7, 4}, {6, 5}},
    {{1, 5}, {7, 9}, {2, 9}, {1, 7}},
    {{2, 7}, {2, 8}, {3, 8}, {4, 8}, {3, 7}}
};

// Another example boundary (not oriented correctly)
const std::vector<std::vector<Point>> exampleBoundary2 = {
    {{2, 2}, {7, 7}, {9, 1}, {14, 3}, {19, 2}, {19, 13}, {15, 16}, {13, 11}, {17, 11}, {17, 9}, {14, 8}, {10, 11}, {11, 16}, {16, 17}, {19, 19}, {1, 17}},
    {{17, 6}, {18, 10}, {15, 5}, {11, 7}, {9, 11}, {8, 10}, {10, 6}, {15, 4}},
    {{15, 7}, {16, 7}, {11, 8}},
    {{15, 12}, {17, 12}, {16, 14}},
    {{3, 6}, {7, 9}, {4, 11}, {6, 14}, {3, 15}, {6, 17}, {2, 16}, {2, 13}, {5, 14}, {3, 11}, {6, 9}},
    {{8, 12}, {9, 12}, {10, 13}, {10, 14}, {9, 16}, {8, 16}, {7, 14}, {7, 13}}
};

void test_domain_hpp() {
    // Example with randomized orientations is oriented correctly
    {
        std::vector<std::vector<Point>> correct = exampleBoundary;
        std::shuffle(correct.begin(), correct.end(), rng);
        
        for(int t = 0; t < 5; ++t) {
            std::vector<std::vector<Point>> boundary = correct;
            for(std::vector<Point>& poly : boundary) {
                if(rng() & 1) {
                    std::reverse(poly.begin(), poly.end());
                }
            }
            
            Domain domain(boundary);
            TEST(domain.boundary() == correct);
        }
    }
    
    // Example domain has the right set of interior points
    {
        Domain domain(exampleBoundary);
        std::vector<Point> interior;
        for(int y = 0; y < 10; ++y) {
            for(int x = 0; x < 10; ++x) {
                Point p(x, y);
                
                // Ignore points that are on the boundary, because their
                // interiority is dependent on the perturbation.
                if(
                    p == Point(7, 1) ||
                    p == Point(8, 2) ||
                    p == Point(4, 4) ||
                    p == Point(2, 5) ||
                    p == Point(4, 6) ||
                    p == Point(6, 7) ||
                    p == Point(4, 7) ||
                    p == Point(1, 6) ||
                    p == Point(3, 9) ||
                    p == Point(4, 9) ||
                    p == Point(5, 9) ||
                    p == Point(6, 9)
                ) {
                    continue;
                }
                
                if(domain.isInteriorPoint(p)) {
                    interior.push_back(p);
                }
            }
        }
        std::vector<Point> correct = {
            {6, 1},
            {1, 2},
            {2, 2},
            {7, 2},
            {1, 3},
            {1, 4},
            {8, 4},
            {3, 5},
            {4, 5},
            {8, 5},
            {2, 6},
            {5, 6},
            {6, 6},
            {8, 6},
            {7, 7},
            {8, 7},
            {5, 8}
        };
        TEST(interior == correct);
    }
    
    // vertexMap is correct for the example domain
    {
        std::unordered_map<Point, std::pair<std::size_t, std::size_t>> correct;
        for(std::size_t polyIdx = 0; polyIdx < exampleBoundary.size(); ++polyIdx) {
            const std::vector<Point>& poly = exampleBoundary[polyIdx];
            for(std::size_t vertIdx = 0; vertIdx < poly.size(); ++vertIdx) {
                correct[poly[vertIdx]] = {polyIdx, vertIdx};
            }
        }
        Domain domain(exampleBoundary);
        TEST(domain.vertexMap() == correct);
    }
    
    // vertexIdxPair works correctly in example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point p(x, y);
                if(domain.vertexMap().count(p)) {
                    TEST(domain.vertexIdxPair(p) == domain.vertexMap().find(p)->second);
                } else {
                    bool throws = false;
                    try {
                        domain.vertexIdxPair(p);
                    } catch(std::domain_error&) {
                        throws = true;
                    }
                    TEST(throws);
                }
            }
        }
    }
    
    // prevVertex and nextVertex work correctly in example domain
    {
        Domain domain(exampleBoundary);
        for(const std::vector<Point>& poly : domain.boundary()) {
            Point a = poly[poly.size() - 2];
            Point b = poly.back();
            for(Point c : poly) {
                TEST(domain.nextVertex(b) == c);
                TEST(domain.prevVertex(b) == a);
                a = b;
                b = c;
            }
        }
    }
    
    // Linearly mapped example is oriented correctly
    for(int t = 0; t < 100; ++t) {
        std::vector<std::vector<Point>> correct = exampleBoundary;
        int64_t a = 0;
        int64_t b = 0;
        int64_t c = 0;
        int64_t d = 0;
        while(a * d == b * c) {
            a = std::uniform_int_distribution<int64_t>(-1000000, 1000000)(rng);
            b = std::uniform_int_distribution<int64_t>(-1000000, 1000000)(rng);
            c = std::uniform_int_distribution<int64_t>(-1000000, 1000000)(rng);
            d = std::uniform_int_distribution<int64_t>(-1000000, 1000000)(rng);
        }
        for(std::vector<Point>& poly : correct) {
            for(Point& v : poly) {
                Point old = v;
                v.x = a * old.x + b * old.y;
                v.y = c * old.x + d * old.y;
            }
        }
        if(a * d < b * c) {
            for(std::vector<Point>& poly : correct) {
                std::reverse(poly.begin(), poly.end());
            }
        }
        
        {
            std::vector<std::vector<Point>> boundary = correct;
            for(std::vector<Point>& poly : boundary) {
                if(rng() & 1) {
                    std::reverse(poly.begin(), poly.end());
                }
            }
            
            Domain domain(boundary);
            TEST(domain.boundary() == correct);
        }
    }
    
    // Extract the edges of the example domain boundary
    std::vector<std::pair<Point, Point>> exampleEdges_;
    for(const std::vector<Point>& poly : exampleBoundary) {
        Point a = poly.back();
        for(Point b : poly) {
            exampleEdges_.emplace_back(a, b);
            a = b;
        }
    }
    const std::vector<std::pair<Point, Point>> exampleEdges = std::move(exampleEdges_);
    
    // Try adding triangles to the domain and check whether domaincreation fails
    // in the right cases
    for(int t = 0; t < 5000; ++t) {
        int ax = std::uniform_int_distribution<int>(-2, 11)(rng);
        int ay = std::uniform_int_distribution<int>(-2, 11)(rng);
        int bx = std::uniform_int_distribution<int>(-2, 11)(rng);
        int by = std::uniform_int_distribution<int>(-2, 11)(rng);
        int cx = std::uniform_int_distribution<int>(-2, 11)(rng);
        int cy = std::uniform_int_distribution<int>(-2, 11)(rng);
        
        Point a(ax, ay);
        Point b(bx, by);
        Point c(cx, cy);
        
        std::vector<std::vector<Point>> boundary = exampleBoundary;
        boundary.push_back({a, b, c});
        
        // See if adding the triangle breaks the boundary
        bool correct = a != b && b != c && c != a;
        for(std::pair<Point, Point> edge : exampleEdges) {
            if(
                intersects(edge.first, edge.second, a, b) ||
                intersects(edge.first, edge.second, b, c) ||
                intersects(edge.first, edge.second, c, a) ||
                edge.first == a ||
                edge.first == b ||
                edge.first == c
            ) {
                correct = false;
                break;
            }
        }
        
        // See if the call throws an error
        bool ok = true;
        try {
            Domain domain(boundary);
        } catch(std::invalid_argument) {
            ok = false;
        }
        
        TEST(correct == ok);
    }
}

void test_visibility_hpp() {
    // isDirectlyVisible throws exactly when domain.isInteriorPoint(x) is false
    // for x = a or x = b
    {
        Domain domain(exampleBoundary);
        for(int x1 = -2; x1 < 12; ++x1) {
            for(int y1 = -2; y1 < 12; ++y1) {
                for(int x2 = -2; x2 < 12; ++x2) {
                    for(int y2 = -2; y2 < 12; ++y2) {
                        Point a(x1, y1);
                        Point b(x2, y2);
                        bool throws = false;
                        try {
                            isVisible(domain, a, b);
                        } catch(std::domain_error&) {
                            throws = true;
                        }
                        TEST(throws == !domain.isInteriorPoint(a) || !domain.isInteriorPoint(b));
                    }
                }
            }
        }
    }
    
    // computePointVisibility throws exactly when domain.isInteriorPoint(center)
    // is false for the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                bool throws = false;
                try {
                    computePointVisibility(domain, center);
                } catch(std::domain_error&) {
                    throws = true;
                }
                TEST(domain.isInteriorPoint(center) != throws);
            }
        }
    }
    
    // computePointVisibility visible verts and edges vectors have the same size
    // and center is correct in example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.isInteriorPoint(center)) {
                    continue;
                }
                PointVisibility vis = computePointVisibility(domain, center);
                TEST(vis.center == center);
                TEST(vis.verts.size() == vis.edges.size());
            }
        }
    }
    
    // computePointVisibility visible verts matches brute force in example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.isInteriorPoint(center)) {
                    continue;
                }
                
                std::vector<Point> correct;
                for(const std::vector<Point>& poly : domain.boundary()) {
                    for(Point p : poly) {
                        bool ok = true;
                        for(const std::vector<Point>& poly2 : domain.boundary()) {
                            Point a = poly2.back();
                            for(Point b : poly2) {
                                if(
                                    !(a == center && b == p) &&
                                    !(a == p && b == center) &&
                                    intersects(a, b, center, p)
                                ) {
                                    ok = false;
                                }
                                a = b;
                            }
                        }
                        if(ok) {
                            correct.push_back(p);
                        }
                    }
                }
                
                std::sort(correct.begin(), correct.end(), [&](Point a, Point b) {
                    return angleLT(center, a, center, b);
                });
                
                PointVisibility vis = computePointVisibility(domain, center);
                TEST(vis.verts == correct);
            }
        }
    }
    
    // computePointVisibility visible edges are correctly oriented and between
    // correct vertices in the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.isInteriorPoint(center)) {
                    continue;
                }
                
                PointVisibility vis = computePointVisibility(domain, center);
                
                for(std::size_t i = 0; i < vis.verts.size(); ++i) {
                    Point a = vis.verts[i];
                    Point b = vis.verts[(i + 1) % vis.verts.size()];
                    Point x, y;
                    std::tie(x, y) = vis.edges[i];
                    
                    TEST(isCCW(center, x, y));
                    TEST(!isCCW(center, a, x));
                    TEST(!isCCW(center, y, b));
                }
            }
        }
    }
    
    // computeVertexVisibility throws exactly when center is not a boundary
    // vertex of the domain in the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                bool throws = false;
                try {
                    computeVertexVisibility(domain, center);
                } catch(std::domain_error&) {
                    throws = true;
                }
                TEST(domain.vertexMap().count(center) != throws);
            }
        }
    }
    
    // computeVertexVisibility verts has exactly one element more than edges
    // and center is correct in the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.vertexMap().count(center)) {
                    continue;
                }
                VertexVisibility vis = computeVertexVisibility(domain, center);
                TEST(vis.verts.size() == vis.edges.size() + 1);
                TEST(vis.center == center);
            }
        }
    }
    
    // computeVertexVisibility verts are CCW ordered in the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.vertexMap().count(center)) {
                    continue;
                }
                
                VertexVisibility vis = computeVertexVisibility(domain, center);
                for(std::size_t i = 1; i < vis.verts.size(); ++i) {
                    TEST(isCCW(center, vis.verts[i - 1], vis.verts[i]));
                }
            }
        }
    }
    
    // computeVertexVisibility visible verts matches brute force in example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.vertexMap().count(center)) {
                    continue;
                }
                
                std::pair<std::size_t, std::size_t> id = domain.vertexIdxPair(center);
                const std::vector<Point>& poly = domain.boundary()[id.first];
                Point next = poly[(id.second + 1) % poly.size()];
                Point prev = poly[(id.second + poly.size() - 1) % poly.size()];
                
                std::vector<Point> correct;
                correct.push_back(next);
                correct.push_back(prev);
                for(const std::vector<Point>& poly : domain.boundary()) {
                    for(Point vertex : poly) {
                        bool ok = true;
                        if(isCCW(center, next, prev)) {
                            if(!isCCW(center, next, vertex) || !isCCW(center, vertex, prev)) {
                                ok = false;
                            }
                        } else {
                            if(!isCCW(center, next, vertex) && !isCCW(center, vertex, prev)) {
                                ok = false;
                            }
                        }
                        for(const std::vector<Point>& poly2 : domain.boundary()) {
                            Point a = poly2.back();
                            for(Point b : poly2) {
                                if(intersects(center, vertex, a, b)) {
                                    ok = false;
                                }
                                a = b;
                            }
                        }
                        if(ok) {
                            correct.push_back(vertex);
                        }
                    }
                }
                
                VertexVisibility vis = computeVertexVisibility(domain, center);
                std::vector<Point> result = vis.verts;
                
                sort(correct.begin(), correct.end(), yCoordLT);
                sort(result.begin(), result.end(), yCoordLT);
                
                TEST(result == correct);
            }
        }
    }
    
    // computeVertexVisibility visible edges are correctly oriented and between
    // correct vertices in the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.vertexMap().count(center)) {
                    continue;
                }
                
                VertexVisibility vis = computeVertexVisibility(domain, center);
                
                for(std::size_t i = 0; i < vis.edges.size(); ++i) {
                    Point a = vis.verts[i];
                    Point b = vis.verts[i + 1];
                    Point x, y;
                    std::tie(x, y) = vis.edges[i];
                    
                    TEST(isCCW(center, x, y));
                    TEST(!isCCW(center, a, x));
                    TEST(!isCCW(center, y, b));
                }
            }
        }
    }
    
    // computeAllVertexVisibilities result matches result of multiple
    // computeVertexVisibility in the example domain
    {
        Domain domain(exampleBoundary);
        for(int x = -2; x < 12; ++x) {
            for(int y = -2; y < 12; ++y) {
                Point center(x, y);
                if(!domain.vertexMap().count(center)) {
                    continue;
                }
                
                std::vector<VertexVisibility> correct;
                for(const std::vector<Point>& poly : domain.boundary()) {
                    for(Point center : poly) {
                        correct.push_back(computeVertexVisibility(domain, center));
                    }
                }
                
                std::vector<VertexVisibility> result = computeAllVertexVisibilities(domain);
                
                auto cmp = [&](const VertexVisibility& a, const VertexVisibility& b) {
                    return yCoordLT(a.center, b.center);
                };
                sort(correct.begin(), correct.end(), cmp);
                sort(result.begin(), result.end(), cmp);
                
                TEST(result.size() == correct.size());
                for(std::size_t i = 0; i < result.size(); ++i) {
                    TEST(result[i].center == correct[i].center);
                    TEST(result[i].verts == correct[i].verts);
                    TEST(result[i].edges == correct[i].edges);
                }
            }
        }
    }
    
    // computeAllVertexVisibilities works correctly in empty domain
    {
        Domain domain;
        std::vector<VertexVisibility> result = computeAllVertexVisibilities(domain);
        TEST(result.empty());
    }
    
    // computeAllVertexVisibilities result matches result of multiple
    // computeVertexVisibility in example domain 2
    {
        Domain domain(exampleBoundary2);
        std::vector<VertexVisibility> correct;
        for(const std::vector<Point>& poly : domain.boundary()) {
            for(Point center : poly) {
                correct.push_back(computeVertexVisibility(domain, center));
            }
        }
        std::vector<VertexVisibility> result = computeAllVertexVisibilities(domain);
        
        auto cmp = [&](const VertexVisibility& a, const VertexVisibility& b) {
            return yCoordLT(a.center, b.center);
        };
        sort(correct.begin(), correct.end(), cmp);
        sort(result.begin(), result.end(), cmp);
        
        TEST(result.size() == correct.size());
        for(std::size_t i = 0; i < result.size(); ++i) {
            TEST(result[i].center == correct[i].center);
            TEST(result[i].verts == correct[i].verts);
            TEST(result[i].edges == correct[i].edges);
        }
    }
}

std::vector<Point> interiorPoints(const Domain& domain) {
    int64_t minX = MaxCoord;
    int64_t maxX = MinCoord;
    int64_t minY = MaxCoord;
    int64_t maxY = MinCoord;
    for(const std::vector<Point>& poly : domain.boundary()) {
        for(Point vertex : poly) {
            minX = std::min(minX, vertex.x);
            maxX = std::max(maxX, vertex.x);
            minY = std::min(minY, vertex.y);
            maxY = std::max(maxY, vertex.y);
        }
    }
    std::vector<Point> ret;
    for(int64_t x = minX + 1; x < maxX; ++x) {
        for(int64_t y = minY + 1; y < maxY; ++y) {
            Point p(x, y);
            if(domain.isInteriorPoint(p)) {
                ret.push_back(p);
            }
        }
    }
    return ret;
}
std::vector<Point> nonInteriorPoints(const Domain& domain) {
    int64_t minX = MaxCoord;
    int64_t maxX = MinCoord;
    int64_t minY = MaxCoord;
    int64_t maxY = MinCoord;
    for(const std::vector<Point>& poly : domain.boundary()) {
        for(Point vertex : poly) {
            minX = std::min(minX, vertex.x);
            maxX = std::max(maxX, vertex.x);
            minY = std::min(minY, vertex.y);
            maxY = std::max(maxY, vertex.y);
        }
    }
    std::vector<Point> ret;
    for(int64_t x = minX; x <= maxX; ++x) {
        for(int64_t y = minY; y <= maxY; ++y) {
            Point p(x, y);
            if(!domain.isInteriorPoint(p)) {
                ret.push_back(p);
            }
        }
    }
    return ret;
}

void test_shortestpath_hpp() {
    // Example boundaries for shortest path computations
    std::vector<std::vector<std::vector<Point>>> spExampleBoundaries;
    spExampleBoundaries.push_back(exampleBoundary);
    spExampleBoundaries.push_back(exampleBoundary2);
    for(auto& a : spExampleBoundaries) {
        for(auto& b : a) {
            for(Point& p : b) {
                p.x *= 5;
                p.y *= 5;
            }
        }
    }
    
    // findShortestPath returned path and path length match and path begins
    // ands ends in correct vertices, does not repeat any elements, and the
    // vertices between endpoints are vertices of the domain
    for(const std::vector<std::vector<Point>>& boundary : spExampleBoundaries) {
        Domain domain(boundary);
        ShortestPathContext sp(domain);
        std::vector<Point> pts = interiorPoints(domain);
        for(int t = 0; t < 100; ++t) {
            Point a = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            Point b = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            std::vector<Point> path;
            double length;
            std::tie(path, length) = sp.findShortestPath(a, b);
            if(path.empty()) {
                TEST(a != b);
                TEST(length == std::numeric_limits<double>::infinity());
            } else {
                TEST(path.front() == a);
                TEST(path.back() == b);
                TEST(path.size() > 1 || a == b);
                double correctLength = 0.0;
                for(std::size_t i = 1; i < path.size(); ++i) {
                    TEST(path[i - 1] != path[i]);
                    correctLength += distance(path[i - 1], path[i]);
                }
                TEST(std::abs(length - correctLength) < 1e-6);
                for(std::size_t i = 1; i + 1 < path.size(); ++i) {
                    TEST(domain.vertexMap().count(path[i]));
                }
            }
        }
    }
    
    // findShortestPath returns correct result when given equal interior points
    for(const std::vector<std::vector<Point>>& boundary : spExampleBoundaries) {
        Domain domain(boundary);
        ShortestPathContext sp(domain);
        std::vector<Point> pts = interiorPoints(domain);
        for(int t = 0; t < 100; ++t) {
            Point p = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            TEST(sp.findShortestPath(p, p) == std::make_pair(std::vector<Point>{p}, 0.0));
        }
    }
    
    // findShortestPath throws std::domain_error when given non-interior points
    for(const std::vector<std::vector<Point>>& boundary : spExampleBoundaries) {
        Domain domain(boundary);
        ShortestPathContext sp(domain);
        std::vector<Point> pts = interiorPoints(domain);
        std::vector<Point> nonPts = nonInteriorPoints(domain);
        for(int t = 0; t < 100; ++t) {
            Point a = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            Point b = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            int mask = std::uniform_int_distribution<int>(1, 3)(rng);
            if(mask & 1) {
                a = nonPts[std::uniform_int_distribution<std::size_t>(0, nonPts.size() - 1)(rng)];
            }
            if(mask & 2) {
                b = nonPts[std::uniform_int_distribution<std::size_t>(0, nonPts.size() - 1)(rng)];
            }
            bool throws = false;
            try {
                sp.findShortestPath(a, b);
            } catch(std::domain_error&) {
                throws = true;
            }
            TEST(throws);
        }
    }
    
    // findShortestPath returns correct path when there is direct visibility
    for(const std::vector<std::vector<Point>>& boundary : spExampleBoundaries) {
        Domain domain(boundary);
        ShortestPathContext sp(domain);
        std::vector<Point> pts = interiorPoints(domain);
        std::vector<Point> nonPts = nonInteriorPoints(domain);
        for(int t = 0; t < 100; ++t) {
            Point a, b;
            while(true) {
                a = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
                b = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
                bool ok = true;
                for(const std::vector<Point>& poly : domain.boundary()) {
                    Point x = poly.back();
                    for(Point y : poly) {
                        ok = ok && !intersects(a, b, x, y);
                        x = y;
                    }
                }
                if(ok) {
                    break;
                }
            }
            std::vector<Point> path;
            double length;
            std::tie(path, length) = sp.findShortestPath(a, b);
            if(a == b) {
                TEST(path.size() == 1);
                TEST(path.front() == a);
                TEST(length == 0.0);
            } else {
                TEST(path.size() == 2);
                TEST(path.front() == a);
                TEST(path.back() == b);
                TEST(std::abs(length - distance(a, b)) < 1e-6);
            }
        }
    }
    
    // findShortestPath path lengths satisfy triangle inequality and are symmetric
    for(const std::vector<std::vector<Point>>& boundary : spExampleBoundaries) {
        Domain domain(boundary);
        ShortestPathContext sp(domain);
        std::vector<Point> pts = interiorPoints(domain);
        for(int t = 0; t < 100; ++t) {
            Point a = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            Point b = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            Point c = pts[std::uniform_int_distribution<std::size_t>(0, pts.size() - 1)(rng)];
            double ab = sp.findShortestPath(a, b).second;
            double ba = sp.findShortestPath(b, a).second;
            double ac = sp.findShortestPath(a, c).second;
            double bc = sp.findShortestPath(b, c).second;
            TEST(std::isfinite(ab) == std::isfinite(ba));
            TEST(!std::isfinite(ab) || std::abs(ab - ba) < 1e-6);
            if(std::isfinite(ac)) {
                TEST(ac - ab - bc < 1e-6);
            } else {
                TEST(!std::isfinite(ab) || !std::isfinite(bc));
            }
        }
    }
}

int main() {
    test_int64_hpp();
    test_geometry_hpp();
    test_intersection_hpp();
    test_domain_hpp();
    test_visibility_hpp();
    test_shortestpath_hpp();
    
    std::cerr << "All tests completed successfully\n";
    
    return 0;
}
