#include <cmath>
#include <cassert>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <valarray>
#include <algorithm>

#include <gtest/gtest.h>


namespace geographic {

double deg2rad(double deg) {
    static const double coeff = std::asin(1) * 2 / 180.0;
    return deg * coeff;
}

double rad2deg(double rad) {
    static const double coeff = 180.0 / std::asin(1);
    return rad * coeff;
}

struct coord3 {
    coord3(double ax, double ay, double az) :x(ax), y(ay), z(az) { }
    coord3() : x(0), y(0), z(0) {}

    coord3& operator+=(coord3 const &other) { return this->add_self(other); }
    coord3& operator-=(coord3 const &other) { return this->add_self(other.scale(-1)); }
    coord3& operator*=(double n) { return this->scale_self(n); }
    coord3& operator/=(double n) { return this->scale_self(1.0 / n); }
    coord3 operator+(coord3 const &other) const { return this->add(other); }
    coord3 operator-(coord3 const &other) const { return this->add(other.scale(-1)); }
    coord3 operator*(double n) const { return this->scale(n); }
    coord3 operator/(double n) { return this->scale(1.0 / n); }
    coord3 operator-() const { return this->scale(-1); }

    static
    coord3 of_degrees(double b, double l, double h) {
        return coord3(deg2rad(b), deg2rad(l), h);
    }

    double length() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    coord3& add_self(coord3 const &other) {
        this->x += other.x;
        this->y += other.y;
        this->z += other.z;
        return *this;
    }

    coord3 add(coord3 const &other) const {
        coord3 copy(*this);
        copy.add_self(other);
        return copy;
    }

    coord3& scale_self(double n) {
        this->x *= n;
        this->y *= n;
        this->z *= n;
        return *this;
    }

    coord3 scale(double n) const {
        coord3 copy(*this);
        copy.scale_self(n);
        return copy;
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
        for (std::size_t i = 0, n = std::min<size_t>(li.size(), nrows * ncols); i < n; ++i)
            _M_data[i] = li[i];
    }

    double elem(int i, int j) const { return _M_data[i * _M_ncols + j]; }
    double& elem(int i, int j) { return _M_data[i * _M_ncols + j]; }
    double operator()(int i, int j) const { return elem(i, j); }
    double& operator()(int i, int j) { return elem(i, j); }
    matrix operator*(matrix const &other) const { return this->multiple(other); }
    matrix& operator+=(matrix const &other) { return this->add_self(other);  }
    matrix& operator*=(double n) { return this->scale_self(n);  }
    matrix& operator-=(matrix const &other) { return this->add_self(other.scale(-1));  }
    matrix& operator/=(double n) { return this->scale_self(1.0 / n);  }
    matrix operator+(matrix const &other) { return this->add(other);  }
    matrix operator*(double n) { return this->scale(n);  }
    matrix operator-(matrix const &other) { return this->add(other.scale(-1));  }
    matrix operator/(double n) { return this->scale(1.0 / n);  }
    matrix operator-() const { return this->scale(-1);  }

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
        matrix copy(*this);
        copy.add_self(other);
        return copy;
    }

    matrix& scale_self(double n) {
        for (std::size_t i = 0, n = _M_data.size(); i < n; ++i)
            _M_data[i] *= n;
        return *this;
    }

    matrix scale(double n) const {
        matrix copy(*this);
        copy.scale_self(n);
        return copy;
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
        out << "{ A =" << this->A << ", f = " << this->f << "}";
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
        out << "{ r =" << this->r << "}";
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

} // namespace geographic


int main(int argc, char* argv[]) {
    return EXIT_SUCCESS;
}
