/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/heat_flux_map.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      template <int dim>
      std::vector<std::vector<std::pair<double, double> > >
      compute_heat_flux_through_boundary_faces (const SimulatorAccess<dim> &simulator_access)
      {
        std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area(simulator_access.get_triangulation().n_active_cells(),
                                                                                 std::vector<std::pair<double, double> >(GeometryInfo<dim>::faces_per_cell,
                                                                                     std::pair<double,double>()));

        // create a quadrature formula based on the temperature element alone.
        const QGauss<dim-1> quadrature_formula (simulator_access.get_fe().base_element(simulator_access.introspection().base_elements.temperature).degree+1);
        FEFaceValues<dim> fe_face_values (simulator_access.get_mapping(),
                                          simulator_access.get_fe(),
                                          quadrature_formula,
                                          update_gradients      | update_values |
                                          update_normal_vectors |
                                          update_q_points       | update_JxW_values);

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_face_values.n_quadrature_points, simulator_access.n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_face_values.n_quadrature_points, simulator_access.n_compositional_fields());

        std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());

        // loop over all of the surface cells and evaluate the heat flux
        typename DoFHandler<dim>::active_cell_iterator
        cell = simulator_access.get_dof_handler().begin_active(),
        endc = simulator_access.get_dof_handler().end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->at_boundary(f))
                {
                  fe_face_values.reinit (cell, f);
                  in.reinit(fe_face_values, cell, simulator_access.introspection(), simulator_access.get_solution(), false);
                  simulator_access.get_material_model().evaluate(in, out);

                  fe_face_values[simulator_access.introspection().extractors.temperature].get_function_gradients (simulator_access.get_solution(),
                      temperature_gradients);

                  // Calculate the normal conductive heat flux given by the formula
                  //   j = - k * n . grad T
                  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                    {
                      const double thermal_conductivity
                        = out.thermal_conductivities[q];

                      heat_flux_and_area[cell->active_cell_index()][f].first += -thermal_conductivity *
                                                                                (temperature_gradients[q] * fe_face_values.normal_vector(q)) *
                                                                                fe_face_values.JxW(q);
                      heat_flux_and_area[cell->active_cell_index()][f].second += fe_face_values.JxW(q);
                    }
                }
        return heat_flux_and_area;
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    HeatFluxMap<dim>::execute (TableHandler &)
    {
      std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area =
        internal::compute_heat_flux_through_boundary_faces (*this);

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output;
      std::vector<std::pair<Point<dim>,double> > stored_values;

      // loop over all of the surface cells and evaluate the heat flux
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f) &&
                (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top" ||
                 this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "bottom"))
              {
                // evaluate position of heat flow to write into output file
                const Point<dim> midpoint_at_surface = cell->face(f)->center();

                const double flux_density = heat_flux_and_area[cell->active_cell_index()][f].first /
                                            heat_flux_and_area[cell->active_cell_index()][f].second;

                // store final position and heat flow
                stored_values.push_back (std::make_pair(midpoint_at_surface, flux_density));
              }


      // Write the solution to an output stream
      for (unsigned int i=0; i<stored_values.size(); ++i)
        {
          output << stored_values[i].first
                 << ' '
                 << stored_values[i].second
                 << std::endl;
        }



      std::string filename = this->get_output_directory() +
                             "heat_flux." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 567;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((dim==2)? "x y" : "x y z")
               << " heat flux" << std::endl;

          // first write out the data we have created locally
          file << output.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                        this->get_mpi_communicator(), &status);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          MPI_Send (&output.str()[0], output.str().size()+1, MPI_CHAR, 0, mpi_tag,
                    this->get_mpi_communicator());
        }

      return std::pair<std::string,std::string>("Writing heat flux map:",
                                                filename);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(HeatFluxMap,
                                  "heat flux map",
                                  "A postprocessor that computes the (conductive) heat flux "
                                  "density across each boundary. The heat density flux is computed in "
                                  "outward direction, i.e., from the domain to the outside, "
                                  "using the formula $-k \\nabla T \\cdot \\mathbf n$, where "
                                  "$k$ is the thermal conductivity as reported by the material "
                                  "model, $T$ is the temperature, and $\\mathbf n$ is the outward "
                                  "normal. Note that the quantity so computed does not include "
                                  "any energy transported across the boundary by material "
                                  "transport in cases where $\\mathbf u \\cdot \\mathbf n \\neq 0$. "
                                  "The integrated heat flux for each boundary can be obtained "
                                  "from the heat flux statistics postprocessor.")
  }
}
