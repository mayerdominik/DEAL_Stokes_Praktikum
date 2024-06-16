#include "RightHandSide.h"

namespace project{
    using namespace dealii;

    template <int dim>
    Tensor<1, dim> RightHandSide<dim>::value(const Point<dim> & /*p*/) const
    {
        return Tensor<1, dim>();
    }


    template <int dim>
    void RightHandSide<dim>::value_list(const std::vector<Point<dim>> &vp,
                                        std::vector<Tensor<1, dim>> &values) const
    {
        for (unsigned int c = 0; c < vp.size(); ++c)
        {
            values[c] = RightHandSide<dim>::value(vp[c]);
        }
    }

    template class RightHandSide<2>;
    template class RightHandSide<3>;

}