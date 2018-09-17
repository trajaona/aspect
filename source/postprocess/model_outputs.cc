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



#include <aspect/postprocess/model_outputs.h>

#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ModelOutputs<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
							   update_gradients |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<Tensor<1,dim> > velocity_values(n_q_points);
      std::vector<Tensor<2,dim> > velocity_gradients (n_q_points, Tensor<2,dim>());

      // Write time, position, velocities, velocity gradient tensor
      // into "model_outputs.txt".
      const std::string filename = (this->get_output_directory() +
                                    "model_ouputs.txt");
      std::ofstream f (filename.c_str());
      f << ("# time  "
            "X  "
            "Y  "
            "Z  ")
        << ("Vx  "
           "Vy  "
           "Vz  ")
        << ("D11  "
		    "D12  "
		    "D13  "
		    "D21  "
            "D22  "
            "D23  "
            "D31  "
            "D32  "
            "D33  ");
       f << '\n';

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_solution(),
                                                                                           velocity_gradients);
            const std::vector<Point<dim> > &position_point = fe_values.get_quadrature_points();
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                f << this->get_time()
                  << ' '
                  << position_point[q]
				  << ' ';

                for (unsigned int i = 0; i < dim; ++i)
				f << velocity_values[q][i] *  year_in_seconds<< ' ';
                for (unsigned int i = 0; i < Tensor<2,dim>::n_independent_components ; ++i)
                f << velocity_gradients[q][Tensor<2,dim>::unrolled_to_component_indices(i)] << ' ';
                f << '\n';
              }
           }

      return std::make_pair (std::string ("Writing model outputs:"),
              filename);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ModelOutputs,
                                  "model outputs",
                                  "A postprocessor that output time, position, velocities, "
                                  "and velocity gradient tensor "
                                  "at all quadrature points")
  }
}
