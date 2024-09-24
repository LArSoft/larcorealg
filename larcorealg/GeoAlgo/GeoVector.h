/**
 * \file GeoVector.h
 *
 * \ingroup GeoAlgo
 *
 * \brief Class def header for a class Point and Vector
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOVECTOR_H
#define BASICTOOL_GEOVECTOR_H

#include "larcorealg/GeoAlgo/GeoAlgoConstants.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <functional>
#include <ostream>
#include <stddef.h>
#include <vector>

namespace geoalgo {

  /**
     \class Vector
     This class represents an n-dimensional vector
  */
  class Vector : public std::vector<double> {
    friend class Trajectory;
    friend class HalfLine;
    friend class LineSegment;
    friend class Sphere;
    friend class GeoAlgo;

  public:
    /// Default ctor
    Vector() : std::vector<double>() {}

    /// Ctor to instantiate with invalid value
    Vector(size_t n) : std::vector<double>(n, kINVALID_DOUBLE) {}

    /// Default ctor w/ a bare std::vector<double>
    Vector(const std::vector<double>& obj) : std::vector<double>(obj) {}

    Vector(double const x, double const y);                 ///< ctor w/ x & y
    Vector(double const x, double const y, double const z); ///< ctor w/ x, y & z
    Vector(TVector3 const& pt);                             ///< ctor w/ TVector3
    Vector(TLorentzVector const& pt);                       ///< ctor w/ TLorentzVector

    void Normalize(); ///< Normalize itself

    bool IsValid() const;    ///< Check if point is valid
    double SqLength() const; ///< Compute the squared length of the vector
    double Length() const;   ///< Compute the length of the vector
    Vector Dir() const;      ///< Return a direction unit vector
    double Phi() const;      ///< Compute the angle Phi
    double Theta() const;    ///< Compute the angle theta

    double SqDist(Vector const& obj) const; ///< Compute the squared distance to another vector
    double Dist(Vector const& obj) const;   ///< Compute the distance to another vector
    double Dot(Vector const& obj) const;    /// Compute a dot product of two vectors
    Vector Cross(Vector const& obj) const;  /// Compute a cross product of two vectors
    double Angle(Vector const& obj) const;  /// Compute an opening angle w.r.t. the given vector

    TLorentzVector ToTLorentzVector()
      const; ///< Convert geovector to TLorentzVector (with 4th element set equal to 0)

    /// Dimensional check for a compatibility
    void compat(Vector const& obj) const;

    /// rotation operations
    void RotateX(double const& theta);
    void RotateY(double const& theta);
    void RotateZ(double const& theta);

  protected:
    /// Compute the squared-distance to another vector w/o dimension check
    double _SqDist_(Vector const& obj) const;
    /// Compute the distance to another vector w/o dimension check
    double _Dist_(Vector const& obj) const;
    /// Compute a dot product w/o dimention check.
    double _Dot_(Vector const& obj) const;
    /// Compute a cross product w/o dimension check.
    Vector _Cross_(Vector const& obj) const;
    /// Compute the angle in degrees between 2 vectors w/o dimension check.
    double _Angle_(Vector const& obj) const;

  public:
    //
    // binary/uniry operators
    //
    inline Vector& operator+=(Vector const& rhs)
    {
      for (size_t i = 0; i < size(); ++i)
        (*this)[i] += rhs[i];
      return *this;
    }

    inline Vector& operator-=(Vector const& rhs)
    {
      for (size_t i = 0; i < size(); ++i)
        (*this)[i] -= rhs[i];
      return *this;
    }

    inline Vector& operator*=(double const rhs)
    {
      for (auto& v : *this)
        v *= rhs;
      return *this;
    }

    inline Vector& operator/=(double const rhs)
    {
      for (auto& v : *this)
        v /= rhs;
      return *this;
    }

    inline Vector& operator=(Vector const& rhs)
    {
      this->resize(rhs.size());
      for (size_t i = 0; i < rhs.size(); ++i)
        (*this)[i] = rhs[i];
      return (*this);
    }

    inline Vector operator+(Vector const& rhs) const
    {
      Vector res((*this));
      res += rhs;
      return res;
    }

    inline Vector operator-(Vector const& rhs) const
    {
      Vector res((*this));
      res -= rhs;
      return res;
    }

    inline double operator*(Vector const& rhs) const
    {
      double res = 0;
      for (size_t i = 0; i < size(); ++i)
        res += (*this)[i] * rhs[i];
      return res;
    }

    inline Vector operator*(double const& rhs) const
    {
      Vector res((*this));
      res *= rhs;
      return res;
    }

    inline Vector operator/(double const& rhs) const
    {
      Vector res((*this));
      res /= rhs;
      return res;
    }

    inline bool operator<(Vector const& rhs) const
    {
      compat(rhs);
      for (size_t i = 0; i < size(); ++i)
        if ((*this)[i] < rhs[i]) return true;
      return false;
    }

    inline bool operator<(double const& rhs) const { return Length() < rhs; }

    inline bool operator==(Vector const& rhs) const
    {
      compat(rhs);
      for (size_t i = 0; i < size(); ++i)
        if ((*this)[i] != rhs[i]) return false;
      return true;
    }

    inline bool operator!=(Vector const& rhs) const { return !(*this == rhs); }

/// Streamer
#ifndef __CINT__
    friend std::ostream& operator<<(std::ostream& o, ::geoalgo::Vector const& a)
    {
      o << "Vector (";
      for (auto const& v : a)
        o << v << " ";
      o << ")";
      return o;
    }
#endif
  };

  /// Point has same feature as Vector
  typedef Vector Vector_t;
  typedef Vector Point_t;
}

// Define a pointer comparison
namespace std {
  template <>
  class less<geoalgo::Vector*> {
  public:
    bool operator()(const geoalgo::Vector* lhs, const geoalgo::Vector* rhs)
    {
      return (*lhs) < (*rhs);
    }
  };
}

#endif
/** @} */ // end of doxygen group
