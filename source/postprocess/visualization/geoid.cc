/*
  Copyright (C) 2016 - 2018 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/postprocess/visualization/geoid.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      Geoid<dim>::execute() const
      {
        AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim> >(this->get_geometry_model()),
                     ExcMessage("The geoid postprocessor is currently only implemented for "
                                "the spherical shell geometry model."));

        std::pair<std::string, Vector<float> *>
        return_value ("geoid",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        const Postprocess::Geoid<dim> &geoid =
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Geoid<dim> >();

        // loop over all of the surface cells and if one less than h/3 away from
        // the top or bottom surface
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              (*return_value.second)(cell->active_cell_index()) = 0.0;
              if (cell->at_boundary())
                {
                  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                    {
                      if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                        {
                          // Get the location of the cell for the expansion
                          const Point<dim> p = cell->face(f)->center();
                          (*return_value.second)(cell->active_cell_index())  = geoid.evaluate(p);

                        }
                    }
                }
            }

        return return_value;
      }

      template <int dim>
      std::list<std::string>
      Geoid<dim>::required_other_postprocessors() const
      {
        return std::list<std::string> (1, "geoid");
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Geoid,
                                                  "geoid",
                                                  "Visualization for the geoid solution. The geoid is given "
                                                  "by the equivalent water column height due to a gravity perturbation. "
                                                  "(Units: m)")
    }
  }
}
