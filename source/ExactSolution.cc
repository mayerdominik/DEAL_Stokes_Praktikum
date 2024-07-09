#include "ExactSolution.h"

namespace project {
	using namespace dealii;
	
	template <int dim>
	void ExactSolution<dim>::vector_value(const Point<dim> &p,
                                          Vector<double>   &values) const
	{
		const double R_x = p[0];
		const double R_y = p[1];

		constexpr double pi  = numbers::PI;
		constexpr double pi2 = numbers::PI * numbers::PI;

		/*values[0] = -std::exp(R_x * (-std::sqrt(25.0 + 4 * pi2) + 5.0)) *
			  std::cos(2 * R_y * pi) + 1;
			  
		values[1] = (1.0L / 2.0L) * (-std::sqrt(25.0 + 4 * pi2) + 5.0) *
			std::exp(R_x * (-std::sqrt(25.0 + 4 * pi2) + 5.0)) *
			std::sin(2 * R_y * pi) / pi;*/
			
			//(-2 * (1-4*(coordsX*coordsX + coordsY*coordsY)) * 8 * coordsY )- velocity_X
			//
			  
		values[0] = -2 * (1 - 4 * (R_x * R_x + R_y * R_y)) * 8 * R_y;
		values[1] = 2 * (1 - 4 * (R_x * R_x + R_y * R_y)) * 8 * R_x;
		
		// x = d/dy Phi
		// y = -d/dx Phi
		
		// Phi = (1-4*(x^2 +y^2) ^2
		// x = - 2 * (1-4*(x*x +y*y)) * 8 * y = -16 * y + 16 * 4 * y * x * x + 16 * 4 * y * y * y 
		// y = 2 * (1-4*(x*x +y*y)) * 8 * x = 16 * x  - 64 * x * x * x - 64 * y * y * x
			  
		//Paraview Calculator: 
		// X-Coordinate: -exp(coordsX*(-sqrt(25.0 + 4 * 3.1415*3.1415) + 5.0)) * cos(2*coordsY * 3.1415) + 1
		//Y-Coordinate: 0.5 * (-sqrt(25.0 + 4 * 3.1415 * 3.1415) + 5.0) *exp(coordsX * (-sqrt(25.0 + 4 * 3.1415 * 3.1415) + 5.0)) * sin(2 * coordsY * 3.1415) / 3.1415


		
		/*values[dim] =
		-1.0L / 2.0L * std::exp(R_x * (-2 * std::sqrt(25.0 + 4 * pi2) + 10.0)) -
		2.0 *
		(-6538034.74494422 +
		 0.0134758939981709 * std::exp(4 * std::sqrt(25.0 + 4 * pi2))) /
		(-80.0 * std::exp(3 * std::sqrt(25.0 + 4 * pi2)) +
		 16.0 * std::sqrt(25.0 + 4 * pi2) *
		   std::exp(3 * std::sqrt(25.0 + 4 * pi2))) -
		1634508.68623606 * std::exp(-3.0 * std::sqrt(25.0 + 4 * pi2)) /
		(-10.0 + 2.0 * std::sqrt(25.0 + 4 * pi2)) +
		(-0.00673794699908547 * std::exp(std::sqrt(25.0 + 4 * pi2)) +
		3269017.37247211 * std::exp(-3 * std::sqrt(25.0 + 4 * pi2))) /
		(-8 * std::sqrt(25.0 + 4 * pi2) + 40.0) +
		0.00336897349954273 * std::exp(1.0 * std::sqrt(25.0 + 4 * pi2)) /
		(-10.0 + 2.0 * std::sqrt(25.0 + 4 * pi2));*/

	}


template class ExactSolution<2>;
template class ExactSolution<3>;
 }
