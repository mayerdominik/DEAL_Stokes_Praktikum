#include "RightHandSide.h"

namespace project{
    using namespace dealii;

    //template <int dim>
    //Tensor<1, dim> RightHandSide<dim>::value(const Point<dim> & /*p*/) const
    //{
    //    return Tensor<1, dim>();
    //}


   template <int dim>
    void RightHandSide<dim>::value_list(const std::vector<Point<dim>> &vp,
                                        std::vector<Tensor<1, dim>> &values) const
    {
        for (unsigned int c = 0; c < vp.size(); ++c)
        {
            values[c] = RightHandSide<dim>::value(vp[c]);
        }
    }
    
  template <int dim>
  Tensor<1, dim> RightHandSide<dim>::value(const Point<dim> &p) const
  {
    Tensor<1, dim>  values;
    const double R_x = p[0];
    const double R_y = p[1];
 
    constexpr double pi  = numbers::PI;
    constexpr double pi2 = numbers::PI * numbers::PI;
 
    values[0] = -1.0L / 2.0L * (-2 * std::sqrt(25.0 + 4 * pi2) + 10.0) *
                  std::exp(R_x * (-2 * std::sqrt(25.0 + 4 * pi2) + 10.0)) -
                0.4 * pi2 * std::exp(R_x * (-std::sqrt(25.0 + 4 * pi2) + 5.0)) *
                  std::cos(2 * R_y * pi) +
                0.1 *
                  Utilities::fixed_power<2>(-std::sqrt(25.0 + 4 * pi2) + 5.0) *
                  std::exp(R_x * (-std::sqrt(25.0 + 4 * pi2) + 5.0)) *
                  std::cos(2 * R_y * pi);
    values[1] = 0.2 * pi * (-std::sqrt(25.0 + 4 * pi2) + 5.0) *
                  std::exp(R_x * (-std::sqrt(25.0 + 4 * pi2) + 5.0)) *
                  std::sin(2 * R_y * pi) -
                0.05 *
                  Utilities::fixed_power<3>(-std::sqrt(25.0 + 4 * pi2) + 5.0) *
                  std::exp(R_x * (-std::sqrt(25.0 + 4 * pi2) + 5.0)) *
                  std::sin(2 * R_y * pi) / pi;
 
   // values[dim - 1] = 0;
    return values;
  }
 

    template class RightHandSide<2>;
    template class RightHandSide<3>;

}
