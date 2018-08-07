# tinygeom2d
Tiny 2D geometry library

## API reference
### geometry.hpp: Geometry primitives
```c++
namespace tinygeom2d {

// The range of allowed point coordinates is [MinCoord, MaxCoord]. Keeping
// coordinates in this range ensures that no overflows happen in the geometry
// functions defined in this file.
const int64_t MinCoord = -((int64_t)1 << 62);
const int64_t MaxCoord = ((int64_t)1 << 62) - 1;

// The type used for points. Points are specified using their x- and
// y-coordinates (positive x-axis points rightwards and positive y-axis points
// upwards).
struct Point {
    // Construct point (0, 0).
    Point();

    // Construct point (x, y). Throws std::out_of_range if one of the
    // coordinates are outside range [MinCoord, MaxCoord].
    Point(int64_t x, int64_t y);

    // Coordinates of the point, must be in range [MinCoord, MaxCoord].
    int64_t x;
    int64_t y;
};

// Make it possible to use Point with std::cout/std::cerr.
std::ostream& operator<<(std::ostream& out, Point p);

// Point comparisons.
bool operator==(const Point& a, const Point& b);
bool operator!=(const Point& a, const Point& b);

// Returns true if the triangle (0, b - a, d - c) is strictly counterclockwise
// oriented (not collinear).
bool isCCW(Point a, Point b, Point c, Point d);

// Returns true if the triangle (a, b, c) is strictly counterclockwise
// oriented (not collinear).
bool isCCW(Point a, Point b, Point c);

// Returns true if the y-coordinate of a is strictly less than the y-coordinate
// of b.
bool yCoordLT(Point a, Point b);

// Returns true if b - a has angle strictly less than d - c, when angles are
// defined as standard directional angles:
//   0 degrees rightwards/positive x
//   90 degrees upwards/positive y
//   180 degrees leftwards/negative x
//   270 degrees downwards/negative y.
// The result is undefined if a = b or c = d.
bool angleLT(Point a, Point b, Point c, Point d);

// Returns Euclidean distance between two points as a double-precision floating
// point number. Used for shortest path lengths because comparing sums of
// Euclidean distances is complicated and expensive.
double distance(Point a, Point b);

// Convenience function that returns a safe multiplicative factor for converting
// floating point coordinates to the correct range for Points, if the original
// x-coordinates are in range [minX, maxX] and the original y-coordinates are in
// range [minY, maxY]. Using the result factor, original point (x, y) in range
// [minX, maxX] x [minY, maxY] should be converted to a Point as follows:
// Point point(std::round(x * factor), std::round(y * factor)). The inverse
// conversion back to the original coordinate space is obtained by
// (point.x / factor, point.y / factor). Note that the conversion is lossy, and
// may change the geometry or map different points to the same point, but in
// most practical settings, the precision should suffice.
double normalizationFactor(
    double minX, double maxX,
    double minY, double maxY
);

}

// Make it possible to use std::unordered_set<Point>.
namespace std {
template <>
struct hash<tinygeom2d::Point> {
    std::size_t operator()(const tinygeom2d::Point& p) const;
};
}
```
### intersection.hpp: Intersection detection
```c++
namespace tinygeom2d {

// Returns true if the interiors of the segments a-b and c-d intersect, i.e.,
// there is a point that is on both segments that is not an endpoint of one of
// the segments.
bool intersects(Point a, Point b, Point c, Point d);

// Returns true if two segments among given segments intersect as defined for
// intersects(a, b, c, d).
// Time complexity: O(n log n), where n = input size.
bool intersects(std::vector<std::pair<Point, Point>> segments);

}
```
### domain.hpp: Polygonal domains
```c++
namespace tinygeom2d {

// A bounded domain in the plane with polygonal boundaries. The domain may be
// disconnected and contain holes.
class Domain {
public:
    // Creates an empty domain.
    Domain();

    // Given boundary polygons, create a bounded domain, i.e. points inside an
    // odd number of polygons are inside the domain. The orientations of the
    // boundary polygons do not matter. No vertex may appear twice in the
    // boundary, and boundary edges may not intersect each other. Each boundary
    // polygon must contain at least three vertices.
    // Throws std::invalid_argument if the boundary is invalid.
    // Time complexity: O(n log n), where n = input size.
    Domain(std::vector<std::vector<Point>> boundary);

    // Returns the boundary polygons of the domain. The polygons are oriented
    // such that outer boundaries are counterclockwise oriented and inner
    // boundaries (holes in the domain) are clockwise oriented. The ordering
    // of the boundary polygons is guaranteed to be the same as in the argument
    // given to the constructor, but the ordering of the vertices in each
    // polygon may have been reversed.
    const std::vector<std::vector<Point>>& boundary() const;

    // Pair identifying a vertex in the boundary, consisting of the index of the
    // boundary polygon and the index of the vertex in the polygon.
    typedef std::pair<std::size_t, std::size_t> IdxPair;

    // Returns a mapping from the vertices of the boundary to their indices in
    // the boundary polygons. For any element (v -> (a, b)) in the map,
    // boundary()[a][b] == v.
    const std::unordered_map<Point, IdxPair>& vertexMap() const;

    // Returns the indices of a boundary vertex in the boundary polygons.
    // Equivalent to vertexMap().find(vertex)->second if the vertex is a vertex
    // of the boundary.
    // Throws std::domain_error if given vertex is not a vertex of the boundary.
    IdxPair vertexIdxPair(Point vertex) const;

    // Returns the position of the next vertex from given vertex on the boundary
    // polygon. 
    // Throws std::domain_error if given vertex is not a vertex of the boundary.
    Point nextVertex(Point vertex) const;

    // Returns the position of the previous vertex from given vertex on the
    // boundary polygon. 
    // Throws std::domain_error if given vertex is not a vertex of the boundary.
    Point prevVertex(Point vertex) const;

    // Returns true if given point is an interior point of the domain. The
    // vertices of the domain do not count as interior points.
    // Time complexity: O(n), where n = boundary size.
    bool isInteriorPoint(Point point) const;
};

}
```
### visibility.hpp: Visibility computations
```c++
namespace tinygeom2d {

// Returns true if interior point a is directly visible from interior point b in
// the domain.
// Throws std::domain_error if a or b is not an interior point of domain.
// Time complexity: O(n), where n = domain boundary size.
bool isVisible(const Domain& domain, Point a, Point b);

// The visibile vertices and edges of a domain from a center point inside the
// domain, computed by computePointVisibility.
struct PointVisibility {
    // The center point, i.e. the point for which the visible vertices and edges
    // are given in this structure.
    Point center;

    // The vertices visible from the center, ordered by angle (in the sense of
    // angleLT).
    std::vector<Point> verts;

    // The parts of edges visible from the center. The edges are ordered such
    // that the visible part of edges()[i] is between vertices verts()[i] and
    // verts()[i + 1] (or verts()[0] if i + 1 == verts.size()) when looking from
    // center(). Note that the endpoints of a visible edge are not always
    // visible - the edge may be behind the limiting vertices verts()[i]
    // and verts()[i + 1]. Each edge (a, b) is ordered such that (center, a, b)
    // is CCW oriented, and a appears before b in the boundary of the domain.
    std::vector<std::pair<Point, Point>> edges;
};

// Compute the visible vertices and edges in the domain from given center point
// in the interior of the domain, i.e. domain.isInteriorPoint(center) returns
// true.
// Throws std::domain_error if the center point is not in the interior.
// Time complexity: O(n log n), where n = domain boundary size.
PointVisibility computePointVisibility(const Domain& domain, Point center);

// The visibile vertices and edges of a domain from a center vertex in the
// boundary of the domain, computed by computeVertexVisibility or
// computeAllVertexVisibilities.
struct VertexVisibility {
    // The center vertex, i.e. the vertex for which the visible vertices and
    // edges are given in this structure.
    Point center;

    // The vertices visible from the center, in CCW order around the center.
    // The first and last vertices are always the next and previous vertices
    // from the center in the boundary polygon, respectively. The center vertex
    // is not listed as visible.
    std::vector<Point> verts;

    // The parts of edges visible from the center. The edges are ordered such
    // that the visible part of edges()[i] is between vertices verts()[i] and
    // verts()[i + 1] when looking from center(). Note that the endpoints of a
    // visible edge are not always visible - the edge may be behind the limiting
    // vertices verts()[i] and verts()[i + 1]. Each edge (a, b) is ordered such
    // that (center, a, b) is CCW oriented, and a appears before b in the
    // boundary of the domain. There is always one less elements in edges than
    // in verts.
    std::vector<std::pair<Point, Point>> edges;
};

// Compute the vertices and edges in the domain visible from a given vertex of
// the domain. Throws std::domain_error if the center point is not a vertex of
// the domain boundary.
// Time complexity: O(n log n), where n = domain boundary size.
VertexVisibility computeVertexVisibility(const Domain& domain, Point center);

// Returns the vector of results of computeVertexVisibility for every vertex of
// the domain. The implementation is faster than calling computeVertexVisibility
// for every vertex separately.
// Time complexity: O(k log n), where k = output size, n = domain boundary size.
std::vector<VertexVisibility> computeAllVertexVisibilities(const Domain& domain);

}
```
### shortestpath.hpp: Shortest path computations
```c++
namespace tinygeom2d {

// Polygonal domain preprocessed for computing shortest paths between interior
// points of the domain.
class ShortestPathContext {
public:
    // Create context for finding shortest paths in given domain.
    // Time complexity: O(k log n) where n = domain boundary size and
    // k = output size of computeAllVertexVisibilities.
    ShortestPathContext(Domain domain);

    // Returns the domain for which this context was created.
    const Domain& domain() const;

    // Find shortest path between interior points a and b of the domain. The
    // return value is a pair (path, length), where path is the sequence of
    // points defining the path as a polygonal line, and length is the Euclidean
    // length of the path. If there is no path between a and b, the path is
    // empty and the length is infinity. Otherwise, the path has always at least
    // one element, the first one is a and the last one is b. No element is
    // repeated on the path. Inexact double precision floating point numbers are
    // used for path lengths, but the result is still always a valid path and
    // its length should be optimal within floating point accuracy.
    // Throws std::domain_error if a or b is not an interior point of the
    // domain.
    // Time complexity: O(n log n + a log n), where n is the size of the domain
    // boundary and a is the number of visibility graph edges expanded by A*
    // (worst case bound a = O(n^2), very pessimistic in most cases).
    typedef std::pair<std::vector<Point>, double> PathResult;
    PathResult findShortestPath(Point a, Point b) const;
};

}
```
### int64.hpp: Portable 64-bit integer multiplication comparison
```c++
namespace tinygeom2d {

using std::int64_t;
using std::uint64_t;

// Returns 1, 0 or -1, if a * b is greater than, equal to or less than c * d,
// respectively.
int cmpMul64(int64_t a, int64_t b, int64_t c, int64_t d);

}
```
