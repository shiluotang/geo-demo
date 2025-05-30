#include <cmath>
#include <cassert>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <valarray>
#include <algorithm>

#include <gtest/gtest.h>

#ifndef USE_UNSTABLE_GEOS_CPP_API
#   define USE_UNSTABLE_GEOS_CPP_API
#   include <geos/geom/Geometry.h>
#   include <geos/geom/GeometryFactory.h>
#   include <geos/geom/CoordinateArraySequenceFactory.h>
#   include <geos/geom/CoordinateSequence.h>
#   include <geos/geom/Triangle.h>
#   include <geos/geom/Triangle.h>
#endif

namespace geographic {

struct degrees {

    static
    double deg2rad(double deg) {
        static const double coeff = std::asin(1) / 90.0;
        return deg * coeff;
    }

    static
    double rad2deg(double rad) {
        static const double coeff = 90.0 / std::asin(1);
        return rad * coeff;
    }
};

struct coord3 {
    coord3(double ax, double ay, double az) :x(ax), y(ay), z(az) { }
    coord3() : x(0), y(0), z(0) { }

    coord3& operator+=(coord3 const &other) { return this->add_self(other); }
    coord3& operator-=(coord3 const &other) { return this->add_self(other.scale(-1)); }
    coord3& operator*=(double n) { return this->scale_self(n); }
    coord3& operator/=(double n) { return this->scale_self(1.0 / n); }
    coord3 operator+(coord3 const &other) const { return this->add(other); }
    coord3 operator-(coord3 const &other) const { return this->add(other.scale(-1)); }
    coord3 operator*(double n) const { return this->scale(n); }
    coord3 operator/(double n) const { return this->scale(1.0 / n); }
    coord3 operator-() const { return this->scale(-1); }
    double operator*(coord3 const &other) const { return this->dot(other); }

    static
    coord3 of_degrees(double b, double l, double h) {
        return coord3(degrees::deg2rad(b), degrees::deg2rad(l), h);
    }

    double length() const {
        return std::sqrt(this->dot(*this));
    }

    coord3& add_self(coord3 const &other) {
        this->x += other.x;
        this->y += other.y;
        this->z += other.z;
        return *this;
    }

    coord3 add(coord3 const &other) const {
        return coord3(*this).add_self(other);
    }

    coord3& scale_self(double n) {
        this->x *= n;
        this->y *= n;
        this->z *= n;
        return *this;
    }

    coord3 scale(double n) const { return coord3(*this).scale_self(n); }

    double dot(coord3 const &other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    void print(std::ostream &out) const {
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<double>::digits10);
        oss << "[" << x << ", " << y << ", " << z << "]";
        out << oss.str();
    }

    double x;
    double y;
    double z;
};

struct blh_deg_view {
    blh_deg_view(coord3 const &c) :_M_ref(c) { }

    double x() const { return degrees::rad2deg(_M_ref.x); }
    double y() const { return degrees::rad2deg(_M_ref.y); }
    double z() const { return _M_ref.z; }

    void print(std::ostream &out) const {
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<double>::digits10);
        oss << "[" << x() << ", " << y() << ", " << z() << "]";
        out << oss.str();
    }

    coord3 const & _M_ref;
};

std::ostream& operator<<(std::ostream &os, blh_deg_view const &v) {
    v.print(os);
    return os;
}

std::ostream& operator<<(std::ostream &os, coord3 const &c) {
    c.print(os);
    return os;
}

struct matrix {
    matrix(
        int nrows,
        int ncols,
        std::valarray<double> li = {})
        : _M_nrows(nrows)
        , _M_ncols(ncols)
        , _M_data(nrows * ncols)
    {
        // FIXME _SCL_SECURE_NO_WARNINGS
        // std::copy_n(std::begin(li), nrows * ncols, std::begin(_M_data));
        for (std::size_t i = 0, n = std::min<size_t>(li.size(), nrows * ncols);
                i < n; ++i)
            _M_data[i] = li[i];
    }

    double elem(int i, int j) const { return _M_data[i * _M_ncols + j]; }
    double& elem(int i, int j) { return _M_data[i * _M_ncols + j]; }
    double operator()(int i, int j) const { return elem(i, j); }
    double& operator()(int i, int j) { return elem(i, j); }
    matrix& operator+=(matrix const &other) { return this->add_self(other); }
    matrix& operator-=(matrix const &other) { return this->add_self(other.scale(-1)); }
    matrix& operator*=(double n) { return this->scale_self(n); }
    matrix& operator/=(double n) { return this->scale_self(1.0 / n); }
    matrix operator+(matrix const &other) const { return this->add(other); }
    matrix operator-(matrix const &other) const { return this->add(other.scale(-1)); }
    matrix operator*(double n) const { return this->scale(n); }
    matrix operator/(double n) const { return this->scale(1.0 / n); }
    matrix operator-() const { return this->scale(-1); }
    matrix operator*(matrix const &other) const { return this->multiple(other); }

    int nrows() const { return _M_nrows; }
    int ncols() const { return _M_ncols; }

    matrix transpose() const {
        matrix r(_M_ncols, _M_nrows);
        for (int i = 0; i < _M_nrows; ++i)
            for (int j = 0; j < _M_ncols; ++j)
                r.elem(j, i) = elem(i, j);
        return r;
    }

    matrix& add_self(matrix const &other) {
        for (std::size_t i = 0, n = _M_data.size(); i < n; ++i)
            _M_data[i] += other._M_data[i];
        return *this;
    }

    matrix add(matrix const &other) const {
        return matrix(*this).add_self(other);
    }

    matrix& scale_self(double n) {
        for (std::size_t i = 0, n = _M_data.size(); i < n; ++i)
            _M_data[i] *= n;
        return *this;
    }

    matrix scale(double n) const {
        return matrix(*this).scale_self(n);
    }

    matrix multiple(matrix const &other) const {
        assert(_M_ncols == other._M_nrows);
        matrix r(_M_nrows, other._M_ncols);
        for (int i = 0, n = this->_M_nrows; i < n; ++i) {
            for (int j = 0, m = other._M_ncols; j < m; ++j) {
                double sum = 0.0;
                for (int k = 0, p = this->_M_ncols; k < p; ++k)
                    sum += this->operator()(i, k) * other(k, j);
                r(i, j) = sum;
            }
        }
        return r;
    }

    void print_row(std::ostream &out, int i) const {
        int width = 13;
        std::ostringstream oss;
        oss << "\n\t[";
        oss << std::right << std::setprecision(7);
        if (_M_ncols > 0) {
            oss << std::setw(width) << this->elem(i, 0);
            for (int j = 1; j < _M_ncols; ++j)
                oss << ", " << std::setw(width) << this->elem(i, j);
        }
        oss << "]";
        out << oss.str();
    }

    void print(std::ostream &out) const {
        out << "[";
        if (_M_nrows > 0) {
            print_row(out, 0);
            for (int i = 1; i < _M_nrows; ++i) {
                out << ", ";
                print_row(out, i);
            }
        }
        out << "\n]";
    }

    static matrix eye(int n) {
        matrix m(n, n);
        for (int i = 0; i < n; ++i)
            m.elem(i, i) = 1.0;
        return m;
    }

    int _M_nrows;
    int _M_ncols;
    std::valarray<double> _M_data;
};

struct matrix4 : public matrix {
    explicit
    matrix4(std::valarray<double> li = {}) : matrix(4, 4, li) { }

    static matrix translate(double dx, double dy, double dz) {
        matrix m = matrix::eye(4);
        m.elem(0, 3) = dx;
        m.elem(1, 3) = dy;
        m.elem(2, 3) = dz;
        return m;
    }

    static matrix rotate(int axis, double theta) {
        matrix m = matrix::eye(4);
        double c = std::cos(theta);
        double s = std::sin(theta);
        // (c s -s c) transform coordinate system
        // (c -s s c) transform object coordinate
        double r[] = {c, s, -s, c};
        int idx = 0;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                m.elem((i + axis + 1) % 3, (j + axis + 1) % 3) = r[idx++];
        return m;
    }

    static matrix rotatex(double theta) { return matrix4::rotate(0, theta); }
    static matrix rotatey(double theta) { return matrix4::rotate(1, theta); }
    static matrix rotatez(double theta) { return matrix4::rotate(2, theta); }

    static matrix translate(coord3 const &delta) {
        return matrix4::translate(delta.x, delta.y, delta.z);
    }

    static matrix hcoord(coord3 const &c) {
        matrix m = matrix(4, 1);
        m.elem(0, 0) = c.x;
        m.elem(1, 0) = c.y;
        m.elem(2, 0) = c.z;
        m.elem(3, 0) = 1;
        return m;
    }
};

matrix operator*(double n, matrix const &m) {
    return m.scale(n);
}

std::ostream& operator<<(std::ostream& os, matrix const &m) {
    m.print(os);
    return os;
}

struct coordsys {
    virtual coord3 blh2xyz(coord3 const &blh) const {
        throw std::runtime_error("not implemented!");
    }

    virtual coord3 xyz2blh(coord3 const &xyz) const {
        throw std::runtime_error("not implemented!");
    }

    virtual coord3 xyz2enu(coord3 const &xyz, coord3 const &orig) const {
        throw std::runtime_error("not implemented!");
    }

    virtual coord3 enu2xyz(coord3 const &enu, coord3 const &orig) const {
        throw std::runtime_error("not implemented!");
    }

    virtual void print(std::ostream &out) const {}
};

std::ostream& operator<<(std::ostream &os, coordsys const &c) {
    c.print(os);
    return os;
}

struct coordsys_ellipsoid : public coordsys {

    coordsys_ellipsoid(double aA, double af) :A(aA), f(af) { }
    coordsys_ellipsoid() :A(0), f(0) { }

    virtual void print(std::ostream &out) const override {
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<double>::digits10)
            << "{A = " << this->A << ", f = " << this->f << "}";
        out << oss.str();
    }

    virtual coord3 blh2xyz(coord3 const &blh) const override {
        double B = blh.x;
        double L = blh.y;
        double H = blh.z;
        double cosB = std::cos(B);
        double sinB = std::sin(B);
        double cosL = std::cos(L);
        double sinL = std::sin(L);
        double e_sq = this->f * (2 - this->f);
        double N = this->A / std::sqrt(1.0 - e_sq * sinB * sinB);
        coord3 xyz;
        xyz.x = (N + H) * cosB * cosL;
        xyz.y = (N + H) * cosB * sinL;
        xyz.z = (N * (1 - e_sq) + H) * sinB;
        return xyz;
    }

    double A;
    double f;

    static coordsys_ellipsoid WGS84;
};

struct coordsys_sphere : public coordsys {
    coordsys_sphere(double ar) : r(ar) {}
    coordsys_sphere() :r(0) { }

    virtual void print(std::ostream &out) const override {
        std::ostringstream oss;
        oss << std::setprecision(std::numeric_limits<double>::digits10)
            << "{r = " << this->r << "}";
        out << oss.str();
    }

    virtual coord3 blh2xyz(coord3 const &blh) const override {
        double B = blh.x;
        double L = blh.y;
        double H = blh.z;
        double cosB = std::cos(B);
        double sinB = std::sin(B);
        double cosL = std::cos(L);
        double sinL = std::sin(L);
        coord3 xyz;
        xyz.x = (this->r + H) * cosB * cosL;
        xyz.y = (this->r + H) * cosB * sinL;
        xyz.z = (this->r + H) * sinB;
        return xyz;
    }

    virtual coord3 xyz2blh(coord3 const &xyz) const override {
        coord3 blh;
        double N = xyz.length();
        blh.x = std::asin(xyz.z / N);
        blh.y = std::atan2(xyz.y, xyz.x);
        blh.z = N - this->r;
        return blh;
    }

    virtual coord3 xyz2enu(coord3 const &xyz, coord3 const &orig) const override {
        double B0 = orig.x;
        double L0 = orig.y;
        double H0 = orig.z;
        double cosB0 = std::cos(B0);
        double sinB0 = std::sin(B0);
        double cosL0 = std::cos(L0);
        double sinL0 = std::sin(L0);
        coord3 oxyz = this->blh2xyz(orig);
        matrix rotate = {
            4, 4,
            {
                -sinL0,				cosL0,				0,		0,
                -sinB0 * cosL0,		-sinB0 * sinL0,		cosB0,	0,
                cosB0 * cosL0,		cosB0 * sinL0,		sinB0,	0,
                0,					0,					0,		1,
            }
        };
        matrix trans = matrix4::translate(-oxyz);
        matrix coord = matrix4::hcoord(xyz);
        matrix v = rotate * trans * coord;
        coord3 enu;
        enu.x = v.elem(0, 0);
        enu.y = v.elem(1, 0);
        enu.z = v.elem(2, 0);
        return enu;
    }

    virtual coord3 enu2xyz(coord3 const &enu, coord3 const &orig) const override {
        double B0 = orig.x;
        double L0 = orig.y;
        double H0 = orig.z;
        double cosB0 = std::cos(B0);
        double sinB0 = std::sin(B0);
        double cosL0 = std::cos(L0);
        double sinL0 = std::sin(L0);
        coord3 oxyz = this->blh2xyz(orig);
        matrix rotate = {
            4, 4,
            {
                -sinL0, -sinB0 * cosL0, cosB0 * cosL0,	0,
                cosL0,  -sinB0 * sinL0, cosB0 * sinL0,	0,
                0,		cosB0,			sinB0,			0,
                0,		0,				0,				1
            }
        };
        matrix trans = matrix4::translate(oxyz);
        matrix coord = matrix4::hcoord(enu);
        matrix v = trans * rotate * coord;
        coord3 xyz;
        xyz.x = v.elem(0, 0);
        xyz.y = v.elem(1, 0);
        xyz.z = v.elem(2, 0);
        return xyz;
    }

    double r;

    static coordsys_sphere WGS84;
};

coordsys_ellipsoid coordsys_ellipsoid::WGS84(6378137.0, 1 / 298.257223563);
coordsys_sphere coordsys_sphere::WGS84(6378137.0);

struct lubound {
    lubound() {
        reset();
    }

    void reset() {
        _M_lbound = std::numeric_limits<double>::infinity();
        _M_ubound = -std::numeric_limits<double>::infinity();
    }

    bool update(double const &value) {
        bool updated = false;
        if (value < _M_lbound) {
            _M_lbound = value;
            updated = true;
        }
        if (value > _M_ubound) {
            _M_ubound = value;
            updated = true;
        }
        return updated;
    }

    double middle() const {
        return (_M_lbound + _M_ubound) * 0.5;
    }

    double _M_lbound;
    double _M_ubound;
};

struct bbox {
    bbox(int n) :_M_boundaries(n) {
    }

    bool expand_dimension(int dimension, double const &value) {
        return _M_boundaries[dimension].update(value);
    }

    void reset() {
        for (lubound &b : _M_boundaries)
            b.reset();
    }

    std::vector<double> center() const {
        std::vector<double> r;
        for (lubound const &b : _M_boundaries)
            r.push_back(b.middle());
        return r;
    }

    std::vector<lubound> _M_boundaries;
};

struct bbox3 : public bbox {
    bbox3() :bbox(3) {}

    bool expand(double x, double y, double z) {
        bool updated = false;
        updated = _M_boundaries[0].update(x) || updated;
        updated = _M_boundaries[1].update(y) || updated;
        updated = _M_boundaries[2].update(z) || updated;
        return updated;
    }

    coord3
    center_coord3() const {
        std::vector<double> const &c = center();
        return coord3(c[0], c[1], c[2]);
    }
};

} // namespace geographic

namespace {

struct buffer_analysis {
    struct parameters {
        parameters()
            : distance(0)
            , quad_segments(8)
            , end_cap_style(1)
        {
        }

        int distance;
        int quad_segments;
        int end_cap_style;
    };

    bool create_rectangular_buffer_area(
            std::vector<geographic::coord3> const &icoords,
            std::vector<geographic::coord3> &ocoords) {
        using namespace geographic;
        using geos::geom::DefaultCoordinateSequenceFactory;
        using geos::geom::CoordinateSequence;
        using geos::geom::Coordinate;
        using geos::geom::CoordinateArraySequenceFactory;
        using geos::geom::CoordinateArraySequence;
        using geos::geom::GeometryFactory;
        using geos::geom::Geometry;
        using geos::geom::LineString;

        // create buffer area
        CoordinateSequence::Ptr coords = CoordinateArraySequenceFactory
            ::instance()->create();
        {
            CoordinateArraySequence &seq =
                dynamic_cast<CoordinateArraySequence&>(*coords);
            for (coord3 const &c : icoords)
                seq.add(Coordinate(c.x, c.y, c.z));
        }
        LineString::Ptr ls = GeometryFactory
            ::getDefaultInstance()->createLineString(std::move(coords));
        if (!ls)
            return false;
        Geometry::Ptr buffer = ls->buffer(
                _M_params.distance,
                _M_params.quad_segments,
                _M_params.end_cap_style);
        if (!buffer)
            return false;
        ocoords.clear();
        coords = buffer->getCoordinates();
        for (int i = 0, n = coords->getSize(); i < n; ++i) {
            Coordinate const &c = coords->getAt(i);
            ocoords.push_back(coord3(c.x, c.y, std::isnan(c.z) ? 0 : c.z));
        }
        return true;
    }

    bool create_geo_buffer_area(
            std::vector<geographic::coord3> const &icoords_blh,
            std::vector<geographic::coord3> &ocoords_blh) {
        using namespace geographic;
        coordsys const &rf = coordsys_sphere::WGS84;

        bbox3 bb3;
        for (coord3 const &c : icoords_blh)
            bb3.expand(c.x, c.y, c.z);
        coord3 icoords_blh_center = bb3.center_coord3();
        std::vector<coord3> icoords_xyz;
        for (coord3 const &c : icoords_blh)
            icoords_xyz.push_back(rf.blh2xyz(c));
        std::vector<coord3> icoords_enu;
        std::vector<coord3> ocoords_enu;
        for (coord3 const &c : icoords_xyz)
            icoords_enu.push_back(rf.xyz2enu(c, icoords_blh_center));
        if (!create_rectangular_buffer_area(icoords_enu, ocoords_enu))
            return false;
        std::vector<coord3> ocoords_xyz;
        for (coord3 const &c : ocoords_enu)
            ocoords_xyz.push_back(rf.enu2xyz(c, icoords_blh_center));
        for (coord3 const &c : ocoords_xyz)
            ocoords_blh.push_back(rf.xyz2blh(c));
        return true;
    }

    parameters _M_params;
};

} // namespace anonymous

TEST(test_geo, test_geo_buffer) {
    using namespace geographic;

    coordsys const &rf = coordsys_sphere::WGS84;
    double distance = 100000;
    std::vector<coord3> icoords;
    std::vector<coord3> ocoords;
    icoords.push_back(coord3::of_degrees(31, 121, 100));
    icoords.push_back(coord3::of_degrees(32, 122, 100));
    icoords.push_back(coord3::of_degrees(34, 123, 100));
    icoords.push_back(coord3::of_degrees(32, 124, 100));
    icoords.push_back(coord3::of_degrees(31, 125, 100));
    icoords.push_back(coord3::of_degrees(30, 126, 100));
    buffer_analysis analysis;
    analysis._M_params.distance = distance;
    analysis.create_geo_buffer_area(icoords, ocoords);
    for (coord3 const &c : icoords)
        std::cout << "icoord = " << blh_deg_view(c) << std::endl;
    for (coord3 const &c : ocoords)
        std::cout << "ocoord = " << blh_deg_view(c) << std::endl;
}
