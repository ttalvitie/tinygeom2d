#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

#include "int64.hpp"
#include "geometry.hpp"

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
    // orientation: check that sign is correct
    TEST(orientation({0, 0}, {1, 0}, {0, 0}, {0, 1}) == 1);
    
    // orientation: random inversion and rotation test
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        Point c = randomPoint();
        Point d = randomPoint();
        int o = orientation(a, b, c, d);
        TEST(o >= -1 && o <= 1);
        TEST(orientation(b, a, c, d) == -o);
        TEST(orientation(a, b, d, c) == -o);
        TEST(orientation(b, a, d, c) == o);
        TEST(orientation(c, d, a, b) == -o);
        TEST(orientation(d, c, a, b) == o);
        TEST(orientation(c, d, b, a) == o);
        TEST(orientation(d, c, b, a) == -o);
    }
    
    // orientation: random degenerate cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        Point c = randomPoint();
        TEST(orientation(a, a, b, c) == 0);
        TEST(orientation(b, c, a, a) == 0);
        TEST(orientation(a, b, a, b) == 0);
        TEST(orientation(a, b, b, a) == 0);
    }
    
    // orientation: inversion and rotation test in small coordinates for
    // three points
    {
        for(int bx = -5; bx <= 5; ++bx) {
        for(int by = -5; by <= 5; ++by) {
        for(int cx = -5; cx <= 5; ++cx) {
        for(int cy = -5; cy <= 5; ++cy) {
            Point a(0, 0);
            Point b(bx, by);
            Point c(cx, cy);
            int o = orientation(a, b, a, c);
            TEST(o >= -1 && o <= 1);
            TEST(orientation(b, c, b, a) == o);
            TEST(orientation(c, a, c, b) == o);
            TEST(orientation(a, c, a, b) == -o);
            TEST(orientation(c, b, c, a) == -o);
            TEST(orientation(b, a, b, c) == -o);
        }
        }
        }
        }
    }
    
    // orientation(a, b, c): random three-point degenerate cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        TEST(orientation(a, a, a, b) == 0);
        TEST(orientation(a, b, a, a) == 0);
        TEST(orientation(b, a, b, a) == 0);
        TEST(orientation(a, a, a, a) == 0);
    }
    
    // cmpY: random degenerate cases
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        TEST(cmpY(a, a) == 0);
    }

    // cmpY: random tests
    for(int i = 0; i < 100; ++i) {
        Point a = randomPoint();
        Point b = randomPoint();
        int o = cmpY(a, b);
        TEST(o >= -1 && o <= 1);
        TEST(cmpY(b, a) == -o);
        if(a.y != b.y) {
            TEST(o == (a.y > b.y) - (a.y < b.y));
        }
    }

    // cmpY: total ordering for small coordinates
    {
        std::vector<Point> points;
        for(int x = -10; x <= 10; ++x) {
            for(int y = -10; y <= 10; ++y) {
                points.emplace_back(x, y);
            }
        }
        
        sort(points.begin(), points.end(), [&](Point a, Point b) {
            return cmpY(a, b) == -1;
        });
        
        for(int i = 0; i < (int)points.size(); ++i) {
            for(int j = 0; j < (int)points.size(); ++j) {
                TEST(cmpY(points[i], points[j]) == (i > j) - (i < j));
            }
        }
    }
    
    // cmpAngle: check that sign is correct
    {
        Point o;
        Point a(1, 1);
        Point b(-1, 1);
        Point c(-1, -1);
        Point d(1, -1);
        TEST(cmpAngle(o, a, o, b) == -1);
        TEST(cmpAngle(o, b, o, c) == -1);
        TEST(cmpAngle(o, c, o, d) == -1);
        TEST(cmpAngle(o, d, o, a) == 1);
    }
    
    // cmpAngle: total ordering for small coordinates
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
            return cmpAngle(origin, a, origin, b) == -1;
        });
        
        for(int i = 0; i < (int)points.size(); ++i) {
            for(int j = 0; j < (int)points.size(); ++j) {
                TEST(cmpAngle(origin, points[i], origin, points[j]) == (i > j) - (i < j));
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

int main() {
    test_int64_hpp();
    test_geometry_hpp();
    
    std::cerr << "All tests completed successfully\n";
    
    return 0;
}
