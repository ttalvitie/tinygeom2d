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
// y-coordinates (for example, x is oriented rightwards and y is oriented
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

// Returns true if the y-coordinate of a is strictly less than, equal to or
// less than the y-coordinate of b, respectively.
bool yCoordLT(Point a, Point b);

// Returns true if b - a has angle strictly less than d - c, when angles are
// defined as standard directional angles:
//   0 degrees rightwards/positive x
//   90 degrees upwards/positive y
//   180 degrees leftwards/negative x
//   270 degrees downwards/negative y.
// The result is undefined if a = b or c = d.
bool angleLT(Point a, Point b, Point c, Point d);

}

// Make it possible to use std::unordered_set<Point>.
namespace std {
template <>
struct hash<tinygeom2d::Point> {
    std::size_t operator()(const tinygeom2d::Point& p) const;
};
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
    Domain(std::vector<std::vector<Point>> boundary);

    // Returns the boundary polygons of the domain. The polygons are oriented
    // such that outer boundaries are counterclockwise oriented and inner
    // boundaries (holes in the domain) are clockwise oriented. The ordering
    // of the boundary polygons is guaranteed to be the same as in the argument
    // given to the constructor, but the ordering of the vertices in each
    // polygon may have been reversed.
    const std::vector<std::vector<Point>>& boundary() const;

    // Returns a mapping from the vertices of the boundary to their locations
    // in the boundary expressed as pairs (a, b), where a is the index of the
    // polygon and b is the index of the vertex in the polygon. Thus, for any
    // element (v -> (a, b)) in the map, boundary()[a][b] == v.
    const std::unordered_map<Point, std::pair<std::size_t, std::size_t>>& vertexMap() const;

    // Returns a pair consisting of the index of the polygon boundary and the
    // index of the vertex in that polygon for given boundary vertex. Equivalent
    // to vertexMap().find(vertex)->second if the vertex is a vertex of the boundary.
    // Throws std::domain_error if vertex is not a vertex of the boundary.
    std::pair<std::size_t, std::size_t> vertexID(Point vertex) const;

    // Returns true if given point is an interior point of the domain. The
    // vertices of the domain do not count as interior points.
    bool isInteriorPoint(Point point) const;

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

// Returns true if segment a is to the left of segment b in the left-to-right
// ordering for segments with intersecting y-coordinate ranges, such that the
// comparison is done at the bottommost horizontal line that intersects both
// segments. If both segments have the same bottom point, the comparison is done
// at a higher horizontal line that is infinitesimally close. The points in both
// segments must be ordered such that the bottom point is first. The segments
// must not have length zero.
bool segmentLeftOfAtBottom(
    std::pair<Point, Point> a,
    std::pair<Point, Point> b
);

// Returns true if two of the given segments intersect.
bool intersects(std::vector<std::pair<Point, Point>> segments);

}
```
### visibility.hpp: Visibility computations
```c++
namespace tinygeom2d {

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
// true.  Throws std::domain_error if the center point is not in the interior.
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
// the boundary.
VertexVisibility computeVertexVisibility(const Domain& domain, Point center);

// Returns the vector of results of computeVertexVisibility for every vertex of
// the domain. The implementation is faster than calling computeVertexVisibility
// for every vertex separately.
std::vector<VertexVisibility> computeAllVertexVisibilities(const Domain& domain);

}
```
