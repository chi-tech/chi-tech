#ifndef _chi_meshvector_h
#define _chi_meshvector_h
#include<cmath>
#include <sstream>
//#include"../chi_mesh.h"
//=============================================== Vertex structure
struct chi_mesh::Vector
{
  double x;
  double y;
  double z;

  Vector()  {
    x=0.0; y=0.0; z=0.0;
  }

  Vector(double a, double b){
    x=a; y=b; z=0.0;
  }

  Vector(double a, double b, double c){
    x=a; y=b; z=c;
  }

  Vector operator+(const Vector& that)
  {
    Vector newVector;
    newVector.x = this->x + that.x;
    newVector.y = this->y + that.y;
    newVector.z = this->z + that.z;

    return newVector;
  }

  Vector operator-(const Vector& that)
  {
    Vector newVector;
    newVector.x = this->x - that.x;
    newVector.y = this->y - that.y;
    newVector.z = this->z - that.z;

    return newVector;
  }

  Vector& operator=(const Vector& that)
  {
    this->x = that.x;
    this->y = that.y;
    this->z = that.z;

    return *this;
  }

  Vector operator*(double value)
  {
    Vector newVector;
    newVector.x = this->x*value;
    newVector.y = this->y*value;
    newVector.z = this->z*value;

    return newVector;
  }

  Vector operator/(double value)
  {
    Vector newVector;
    newVector.x = this->x/value;
    newVector.y = this->y/value;
    newVector.z = this->z/value;

    return newVector;
  }

  Vector Cross(const Vector& that)
  {
    Vector newVector;
    newVector.x = this->y*that.z - this->z*that.y;
    newVector.y = this->z*that.x - this->x*that.z;
    newVector.z = this->x*that.y - this->y*that.x;

    return newVector;
  }

  double Dot(const Vector& that)
  {
    double value = 0.0;
    value += this->x*that.x;
    value += this->y*that.y;
    value += this->z*that.z;

    return value;
  }

  double Norm()
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

  friend std::ostream & operator<< (std::ostream& out, Vector& v)
  {
    out << "[" << v.x << " " << v.y << " " << v.z << "]";

    return out;
  }

  std::string PrintS()
  {
    std::stringstream out;
    out << "[" << x << " " << y << " " << z << "]";

    return out.str();
  }
};

#endif
