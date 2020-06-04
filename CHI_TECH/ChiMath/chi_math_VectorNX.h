#ifndef chi_math_VectorNX_h
#define chi_math_VectorNX_h

#include <iostream>

#include <vector>
#include <cmath>
#include <sstream>

#include<type_traits>

namespace chi_math{
    /**  
     * \author JerryW
    */
    template <int N>
    struct VectorNX{
        std::vector<NumberFormat> elements;
        //default constructor
        VectorNX()
        {
            static_assert(std::is_floating_point<NumberFormat>::value,
            "Only floating point number formats are supported for VectorNX." );
            elements.resize(N,NumberFormat());
        }

        //constructor with value
        VectorNX(NumberFormat value)
        {
            static_assert(std::is_floating_point<NumberFormat>::value,
            "Only floating point number formats are supported for VectorNX." );            
            elements.resize(N,value);
        }

        /**Returns a copy of the value at the given index.*/
        NumberFormat operator[] (int i)
        {
            return elements[i];
        }

        /**Returns a reference of the value at the given index.*/
        NumberFormat& operator ()(int i)
        {
            return elements[i];
        }

        /**Component-wise addition.*/
        VectorNX operator+(const VectorNX& rhs) const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N;++i)
                newVector.elements[i] = elements[i] + rhs.elements[i];
            return newVector;
        }

        /**Component-wise subtraction.*/
        VectorNX operator-(const VectorNX& rhs) const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N;++i)
                newVector.elements[i] = elements[i] - rhs.elements[i];
            return newVector;
        }

        /**Component-wise copy.*/
        VectorNX& operator=(const VectorNX& rhs) 
        {
            elements = rhs.elements;
            return *this;
        }

        /**Vector component-wise multiplication by scalar.*/
        VectorNX operator*(NumberFormat value) const
        {
            VectorNX<N> newVector;
            for (int i = 0;i<N;++i)
                newVector.elements[i] = elements[i] * value;
            return newVector;
        }

        /**Vector component-wise multiplication.*/
        VectorNX operator*(const VectorNX& rhs) const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N ; ++i)
                newVector.elements[i] = elements[i] * rhs.elements[i];
            return newVector;
        }

        /**Vector component-wise division by scalar.*/
        VectorNX operator/(NumberFormat value) const
        {
            VectorNX<N> newVector;
            for (int i = 0;i<N;++i)
                newVector.elements[i] = elements[i] / value;
            return newVector;
        }

        /**Vector component-wise division.*/
        VectorNX operator/(const VectorNX& rhs) const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N ; ++i)
                newVector.elements[i] = elements[i] / rhs.elements[i];
            return newVector;
        }

        //cross product

        /**Vector dot-product.*/
        NumberFormat Dot(const VectorNX& rhs) const
        {
            NumberFormat value = 0.0;
            for (int i = 0; i<N; ++i)
                value += elements[i]*rhs.elements[i];
            return value;
        }        
        /**Computes the L2-norm of the vector. Otherwise known as the length of
        * a 3D vector.*/
        NumberFormat Norm() const
        {
            NumberFormat value = 0.0;
            for (int i = 0; i<N;++i)
                value += elements[i]*elements[i];
            value = sqrt(value);
            return value;
        }

        /**Computes the square of the L2-norm of the vector. This eliminates the
        * usage of the square root and is therefore less expensive that a proper
        * L2-norm. Useful if only comparing distances.*/
        NumberFormat NormSquare() const
        {
            NumberFormat value = 0.0;
            for (int i = 0; i<N;++i)
                value += elements[i]*elements[i];

            return value;            
        }

        /**Normalizes the vector in-place.*/
        void Normalize()
        {
            NumberFormat norm = this->Norm();
            for (int i = 0;i<N;++i)
                elements[i] = elements[i]/norm;
        }

        /**Returns a normalized version of the vector.*/
        VectorNX Normalized() const
        {
            NumberFormat norm = this->Norm();
            VectorNX<N> newVector;
            for (int i = 0;i<N;++i)
                newVector.elements[i] = elements[i]/norm;
            return newVector;
        }

        /**Returns a vector v^* where each element is inverted provided
        * that it is greater than the given tolerance, otherwise the offending entry
        * is set to 0.0.*/
        VectorNX InverseZeroIfSmaller(NumberFormat tol) const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N;++i)
                newVector.elements[i] = (std::fabs(elements[i])>tol) ? 1.0/elements[i] : 0.0;
            
            return newVector;
        }

        /**Returns a vector v^* where each element is inverted provided
        * that it is greater than the given tolerance, otherwise the offending entry
        * is set to 1.0.*/
        VectorNX InverseOneIfSmaller(NumberFormat tol) const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N;++i)
                newVector.elements[i] = (std::fabs(elements[i])>tol) ? 1.0/elements[i] : 1.0;
            return newVector;
        }

        /**Returns a vector v^* where each element is inverted provided
        * that the inversion is not infinite, otherwise it is zeroed.*/
        VectorNX InverseZeroIfInf() const
        {
            VectorNX<N> newVector;
            for (int i = 0; i<N; ++i){
                NumberFormat dn_inv = 1.0/elements[i];
                newVector.elements[i] = (std::isinf(dn_inv))? dn_inv : 0.0;
            }

            return newVector;
        }
        
        /**prints the vector to standard cout*/
        void Print() const
        {
            for (int i = 0; i<N-1; i++)
                std::cout<<elements[i]<<" ";
            std::cout<<elements[N-1];
        }

        //overloading <<
        

        /**Prints the vector to a string and then returns the string.*/
        std::string PrintS() const
        {
            std::stringstream out;
            out<<"[";
            for (int i = 0; i<N-1 ; ++i)
                out<<elements[i]<<" ";
            out<<elements[N-1]<<"]";

            return out.str();

        }

    };
}

#endif