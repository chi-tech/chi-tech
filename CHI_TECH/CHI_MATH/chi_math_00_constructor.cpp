#include"chi_math.h"

CHI_MATH::CHI_MATH()
{
	this->IVP_coeff1.SetSizeZero(6, 4);
	//ai
	IVP_coeff1.sij(1, 1, 16.0 / 135.0);
	IVP_coeff1.sij(2, 1, 0.0);
	IVP_coeff1.sij(3, 1, 6656.0/12825.0);
	IVP_coeff1.sij(4, 1, 28561.0/56430.0);
	IVP_coeff1.sij(5, 1, -9.0/50.0);
	IVP_coeff1.sij(6, 1, 2.0/55.0);

	//ai-bi
	IVP_coeff1.sij(1, 2,1.0/360.0);
	IVP_coeff1.sij(2, 2,0.0);
	IVP_coeff1.sij(3, 2,-128.0/4275.0);
	IVP_coeff1.sij(4, 2,-2197.0/75240.0);
	IVP_coeff1.sij(5, 2,1.0/50.0);
	IVP_coeff1.sij(6, 2,2.0/55.0);

	//bi
	IVP_coeff1.sij(1, 3, IVP_coeff1.ij(1, 1) - IVP_coeff1.ij(1, 2));
	IVP_coeff1.sij(2, 3, IVP_coeff1.ij(2, 1) - IVP_coeff1.ij(2, 2));
	IVP_coeff1.sij(3, 3, IVP_coeff1.ij(3, 1) - IVP_coeff1.ij(3, 2));
	IVP_coeff1.sij(4, 3, IVP_coeff1.ij(4, 1) - IVP_coeff1.ij(4, 2));
	IVP_coeff1.sij(5, 3, IVP_coeff1.ij(5, 1) - IVP_coeff1.ij(5, 2));
	IVP_coeff1.sij(6, 3, IVP_coeff1.ij(6, 1) - IVP_coeff1.ij(6, 2));

	//ci
	IVP_coeff1.sij(1, 4, 0.0);
	IVP_coeff1.sij(2, 4, 0.25);
	IVP_coeff1.sij(3, 4, 3.0/8.0);
	IVP_coeff1.sij(4, 4, 12.0/13.0);
	IVP_coeff1.sij(5, 4, 1.0);
	IVP_coeff1.sij(6, 4, 0.5);


	this->IVP_coeff2.SetSizeZero(6, 5);
	//row 2
	IVP_coeff2.sij(2, 1, 0.25);
	//row 3
	IVP_coeff2.sij(3, 1, 3.0/32.0);
	IVP_coeff2.sij(3, 2, 9.0/32.0);
	//row 4
	IVP_coeff2.sij(4, 1, 1932.0/2197.0);
	IVP_coeff2.sij(4, 2, -7200.0/2197.0);
	IVP_coeff2.sij(4, 3, 7296.0/2197.0);
	//row 5
	IVP_coeff2.sij(5, 1,439.0/216.0);
	IVP_coeff2.sij(5, 2,-8.0);
	IVP_coeff2.sij(5, 3,3680.0/513.0);
	IVP_coeff2.sij(5, 4,-845.0/4104.0);
	//row 6
	IVP_coeff2.sij(6, 1,-8.0/27.0);
	IVP_coeff2.sij(6, 2,2.0);
	IVP_coeff2.sij(6, 3,-3544.0/2565.0);
	IVP_coeff2.sij(6, 4,1859.0/4104.0);
	IVP_coeff2.sij(6, 5,-11.0/40.0);
}