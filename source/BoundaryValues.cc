#include "BoundaryValues.h"


namespace project{
    using namespace dealii;
    template <int dim>
    double BoundaryValues<dim>::value(const Point<dim> & p,
                                      const unsigned int component) const
    {
        Assert(component < this->n_components,
               ExcIndexRange(component, 0, this->n_components));

        if (component == 0)
        {
            if (p[0] < 0)
                return -1;
            else if (p[0] > 0)
                return 1;
            else
                return 0;
        }

        return 0;
    }


    template <int dim>
    void BoundaryValues<dim>::vector_value(const Point<dim> &p,
                                           Vector<double> &  values) const
    {
        for (unsigned int c = 0; c < this->n_components; ++c)
            values(c) = BoundaryValues<dim>::value(p, c);
    }
template class BoundaryValues<2>;
template class BoundaryValues<3>;
}