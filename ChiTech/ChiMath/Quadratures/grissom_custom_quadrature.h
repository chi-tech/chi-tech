//
// Created by john on 1/5/21.
//

#ifndef CHITECH_GRISSOM_CUSTOM_QUADRATURE_H
#define CHITECH_GRISSOM_CUSTOM_QUADRATURE_H

#include "angular_quadrature_base.h"

namespace chi_math
{
    class GrissomCustomQuadrature : public AngularQuadrature
    {
    private:
        const std::string filename;
    public:
        explicit
        GrissomCustomQuadrature(std::string& in_filename) :
        AngularQuadrature(AngularQuadratureType::Arbitrary),
        filename(in_filename)
        {}

        void BuildDiscreteToMomentOperator(int scatt_order, bool oneD) override;
        void BuildMomentToDiscreteOperator(int scatt_order, bool oneD) override;
    };
}

#endif //CHITECH_GRISSOM_CUSTOM_QUADRATURE_H
