#include <algorithm>
#include <iostream>
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

void test_domain_hpp() {
    // Example correctly oriented bounded domain boundary
    const std::vector<std::vector<Point>> example = {
        {{0, 4}, {1, 1}, {3, 2}, {6, 0}, {9, 3}, {8, 8}},
        {{2, 4}, {7, 6}, {8, 3}, {5, 1}, {2, 3}},
        {{5, 2}, {6, 2}, {5, 5}, {3, 3}},
        {{4, 3}, {5, 4}, {5, 3}},
        {{7, 3}, {7, 4}, {6, 5}},
        {{1, 5}, {7, 9}, {2, 9}, {1, 7}},
        {{2, 7}, {2, 8}, {3, 8}, {4, 8}, {3, 7}}
    };
    
    // Example with randomized orientations is oriented correctly
    {
        std::vector<std::vector<Point>> correct = example;
        std::shuffle(correct.begin(), correct.end(), rng);
        
        for(int t = 0; t < 5; ++t) {
            std::vector<std::vector<Point>> boundary = correct;
            for(std::vector<Point>& poly : boundary) {
                if(rng() & 1) {
                    std::reverse(poly.begin(), poly.end());
                }
            }
            
            Domain domain = Domain::createBounded(boundary);
            TEST(domain.boundary() == correct);
        }
        
        for(std::vector<Point>& poly : correct) {
            std::reverse(poly.begin(), poly.end());
        }
        
        for(int t = 0; t < 5; ++t) {
            std::vector<std::vector<Point>> boundary = correct;
            for(std::vector<Point>& poly : boundary) {
                if(rng() & 1) {
                    std::reverse(poly.begin(), poly.end());
                }
            }
            
            Domain domain = Domain::createUnbounded(boundary);
            TEST(domain.boundary() == correct);
        }
    }
    
    // Linearly mapped example is oriented correctly
    for(int t = 0; t < 100; ++t) {
        std::vector<std::vector<Point>> correct = example;
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
            
            Domain domain = Domain::createBounded(boundary);
            TEST(domain.boundary() == correct);
        }
        
        for(std::vector<Point>& poly : correct) {
            std::reverse(poly.begin(), poly.end());
        }
        
        {
            std::vector<std::vector<Point>> boundary = correct;
            for(std::vector<Point>& poly : boundary) {
                if(rng() & 1) {
                    std::reverse(poly.begin(), poly.end());
                }
            }
            
            Domain domain = Domain::createUnbounded(boundary);
            TEST(domain.boundary() == correct);
        }
    }
    
    // Extract the edges of the example domain boundary
    std::vector<std::pair<Point, Point>> exampleEdges_;
    for(const std::vector<Point>& poly : example) {
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
        
        std::vector<std::vector<Point>> boundary = example;
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
            if(rng() & 1) {
                Domain::createBounded(boundary);
            } else {
                Domain::createUnbounded(boundary);
            }
        } catch(std::invalid_argument) {
            ok = false;
        }
        
        TEST(correct == ok);
    }
    
    // vertex, prevVertex and nextVertex in the example
    {
        Domain domain = Domain::createBounded(example);
        for(std::size_t polyIdx = 0; polyIdx < example.size(); ++polyIdx) {
            const std::vector<Point>& poly = example[polyIdx];
            for(std::size_t vertIdx = 0; vertIdx < poly.size(); ++vertIdx) {
                size_t prevVertIdx = vertIdx ? vertIdx - 1 : poly.size() - 1;
                size_t nextVertIdx = vertIdx == poly.size() - 1 ? 0 : vertIdx + 1;
                
                VertexID id = {polyIdx, vertIdx};
                TEST(domain.vertex(id) == poly[vertIdx]);
                TEST(domain.prevVertex(id) == poly[prevVertIdx]);
                TEST(domain.nextVertex(id) == poly[nextVertIdx]);
            }
        }
    }
}

int main() {
    test_int64_hpp();
    test_geometry_hpp();
    test_intersection_hpp();
    test_domain_hpp();
    
    std::cerr << "All tests completed successfully\n";
    
    return 0;
}
