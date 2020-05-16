#ifndef _chi_meshvector_h
#define _chi_meshvector_h
#include<cmath>
#include <sstream>

//=============================================== General 3D vector structure
/**General vertex structure*/
struct chi_mesh::Vector3
{
  double x;
  double y;
  double z;

  Vector3()  {
    x=0.0; y=0.0; z=0.0;
  }

  Vector3(double a){
    x=0.0; y=0.0; z=a;
  }

  Vector3(double a, double b){
    x=a; y=b; z=0.0;
  }

  Vector3(double a, double b, double c){
    x=a; y=b; z=c;
  }

  Vector3 operator+(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->x + that.x;
    newVector.y = this->y + that.y;
    newVector.z = this->z + that.z;

    return newVector;
  }

  Vector3 operator-(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->x - that.x;
    newVector.y = this->y - that.y;
    newVector.z = this->z - that.z;

    return newVector;
  }

  Vector3& operator=(const Vector3& that)
  {
    this->x = that.x;
    this->y = that.y;
    this->z = that.z;

    return *this;
  }

  Vector3 operator*(double value) const
  {
    Vector3 newVector;
    newVector.x = this->x*value;
    newVector.y = this->y*value;
    newVector.z = this->z*value;

    return newVector;
  }

  Vector3 operator/(double value) const
  {
    Vector3 newVector;
    newVector.x = this->x/value;
    newVector.y = this->y/value;
    newVector.z = this->z/value;

    return newVector;
  }

  Vector3 Cross(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->y*that.z - this->z*that.y;
    newVector.y = this->z*that.x - this->x*that.z;
    newVector.z = this->x*that.y - this->y*that.x;

    return newVector;
  }

  double Dot(const Vector3& that) const
  {
    double value = 0.0;
    value += this->x*that.x;
    value += this->y*that.y;
    value += this->z*that.z;

    return value;
  }

  void Normalize()
  {
    double norm = this->Norm();

    x /= norm;
    y /= norm;
    z /= norm;
  }

  Vector3 Normalized() const
  {
    double norm = this->Norm();

    Vector3 newVector;
    newVector.x = this->x/norm;
    newVector.y = this->y/norm;
    newVector.z = this->z/norm;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that it is greater than the given tolerance, otherwise it is zeroed.*/
  Vector3 InverseZeroIfGreater(double tol) const
  {
    Vector3 newVector;
    newVector.x = (std::fabs(this->x)>tol)? 1.0/this->x : 0.0;
    newVector.y = (std::fabs(this->y)>tol)? 1.0/this->y : 0.0;
    newVector.z = (std::fabs(this->z)>tol)? 1.0/this->z : 0.0;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that the inversion is not infinite, otherwise it is zeroed.*/
  Vector3 InverseZeroIfInf() const
  {
    Vector3 newVector;
    double dx_inv = 1.0/this->x;
    double dy_inv = 1.0/this->y;
    double dz_inv = 1.0/this->z;

    newVector.x = (std::isinf(dx_inv))? dx_inv : 0.0;
    newVector.y = (std::isinf(dy_inv))? dy_inv : 0.0;
    newVector.z = (std::isinf(dz_inv))? dz_inv : 0.0;

    return newVector;
  }

  double Norm() const
  {
    double value = 0.0;
    value += this->x*this->x;
    value += this->y*this->y;
    value += this->z*this->z;

    value = sqrt(value);

    return value;
  }

  double NormSquare()
  {
    double value = 0.0;
    value += this->x*this->x;
    value += this->y*this->y;
    value += this->z*this->z;

    return value;
  }

  void Print()
  {
    std::cout<<this->x << " ";
    std::cout<<this->y << " ";
    std::cout<<this->z;
  }

  friend std::ostream & operator<< (std::ostream& out, Vector3& v)
  {
    out << "[" << v.x << " " << v.y << " " << v.z << "]";

    return out;
  }

  std::string PrintS() const
  {
    std::stringstream out;
    out << "[" << x << " " << y << " " << z << "]";

    return out.str();
  }
};

chi_mesh::Vector3 operator*(double value,const chi_mesh::Vector3& a);
//{
//  chi_mesh::Vector3 newVector;
//  newVector.x = a.x*value;
//  newVector.y = a.y*value;
//  newVector.z = a.z*value;
//
//  return newVector;
//}

#endif
