/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/visualization/velocity_inclination.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VelocityInclination<dim>::
      VelocityInclination ()
        :
        DataPostprocessorScalar<dim> ("velocity_inclination",
                                      update_values | update_q_points)
      {}



      template <int dim>
      void
      VelocityInclination<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian, spherical and ellipsoidal.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
        if (dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()) != 0)
          {
            coordinate_system = Utilities::Coordinates::cartesian;
          }
        else
          {
            coordinate_system = Utilities::Coordinates::spherical;
          }

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // extract the primal variables
            Tensor<1,dim> v;
            for (unsigned int d=0; d<dim; ++d)
              v[d] = input_data.solution_values[q][d];

            computed_quantities[q](0) = Utilities::compute_vector_inclination_wrt_horizontal(input_data.evaluation_points[q],v,coordinate_system);
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VelocityInclination,
                                                  "velocity inclination",
                                                  "A visualization output object that generates output "
                                                  "for the norm of the strain rate, i.e., for the quantity "
                                                  "$\\sqrt{\\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)}$ "
                                                  "in the incompressible case and "
                                                  "$\\sqrt{[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]}$ "
                                                  "in the compressible case.")
    }
  }
}
