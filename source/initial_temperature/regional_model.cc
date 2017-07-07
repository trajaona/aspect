/*
 * lithosphere_geotherm.cc
 *
 *  Created on: May 26, 2016
 *      Author: trajaona
 */

/*
   Copyright (C) 2016 by the authors of the ASPECT code.

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


#include <aspect/utilities.h>
#include <fstream>
#include <iostream>
#include <deal.II/base/std_cxx11/array.h>
#include "../../include/aspect/initial_temperature/regional_model.h"


namespace aspect
{
  namespace InitialTemperature
  {
    namespace internal
    {
      namespace Tomography
      {

        class VsPerturbationLookup
        {
          public:

            VsPerturbationLookup (const std::string &filename,
                                  const double &h,
                                  const double &v,
                                  const MPI_Comm &comm)
            {
              h_delta = h;
              v_delta = v;
              /**
               * Read data from disk and distribute among processe
               */
              std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

              /**
               * Reading data lines.
               */
              double latitude_Vs, longitude_Vs, depth_Vs, perturbation ;
              while (in >> latitude_Vs >> longitude_Vs >> depth_Vs >> perturbation)
                {
                  latitudes_Vs.push_back(latitude_Vs);
                  longitudes_Vs.push_back(longitude_Vs);
                  depths_Vs.push_back(depth_Vs*1000.);
                  perturbations.push_back(perturbation);
                }
            }

            bool
            is_in_tomo_domain(const std_cxx11::array<double,3> &wgscoord) const
            {
              for (unsigned int i = 0; i <= latitudes_Vs.size();)
                {
                  if (std::fabs(wgscoord[2] - latitudes_Vs[i]) <= (h_delta/2) &&
                      std::fabs(wgscoord[1] - longitudes_Vs[i]) <= (h_delta/2) &&
                      std::fabs(wgscoord[0] - depths_Vs[i])  <= (v_delta/2))

                    return true;
                  else
                    i++;
                }
              return false;
            }

            double get_perturbation (const  std_cxx11::array<double,3> &wgscoord) const
            {
              for (unsigned int i = 0; i <= latitudes_Vs.size();)
                {
                  if (std::fabs(wgscoord[2] - latitudes_Vs[i]) <= h_delta/2 &&
                      std::fabs(wgscoord[1] - longitudes_Vs[i]) <= h_delta/2 &&
                      std::fabs(wgscoord[0] - depths_Vs[i])  <= v_delta/2)

                    return perturbations[i];
                  else
                    i++;
                }
              Assert (false, ExcInternalError());
              return 0;
            }

          private:
            std::vector<double>  latitudes_Vs;
            std::vector<double>  longitudes_Vs;
            std::vector<double>  depths_Vs;
            std::vector<double>  perturbations;
            double h_delta;
            double v_delta;
        };
      }
    }

    template <int dim>
    void
    RegionalModel<dim>::initialize()
    {
      perturbation_lookup.reset(new internal::Tomography::VsPerturbationLookup(data_directory+perturbation_file_name,
                                                                               h_delta,
                                                                               v_delta,
                                                                               this->get_mpi_communicator()));
    }

    template <>
    double
    RegionalModel<2>::initial_temperature (const Point<2> &) const
    {
      AssertThrow (false, ExcMessage ("The 'regional geotherm' initial temperature plugin is only"
                                      "implemented for 3d cases."));
      return 0;
    }

    template <int dim>
    double
    RegionalModel<dim>::initial_temperature (const Point<dim> &position) const
    {
      const double background_temperature = this->get_adiabatic_conditions().temperature(position);

      const double depth                   = this->get_geometry_model().depth(position);
      std_cxx11::array<double,dim> wcoord  = Utilities::Coordinates::WGS84_coordinates(position);

      wcoord[0] = depth;
      double perturbation;
      if (perturbation_lookup->is_in_tomo_domain(wcoord)==true)
        perturbation = (perturbation_lookup->get_perturbation(wcoord));
      else
        perturbation = 0.0;
      const double density_perturbation = dlnvs_to_dlnrho * perturbation;

      const double temperature_perturbation =  (-1./thermal_alpha) * density_perturbation;
      // std::cout<<depth<<"  "<<background_temperature<<"  "<<perturbation<<"  "<<temperature_perturbation<<std::endl;

      return background_temperature+temperature_perturbation ;
    }

    template <int dim>
    void
    RegionalModel<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection("Regional model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-conditions/adiabatic-boundary/",
                             Patterns::DirectoryName (),
                             "The path to the Vs perturbation data files");
          prm.declare_entry ("Horizontal resolution", "1.0",
                             Patterns::Double (0),
                             "horizontal resolution of the 3D tomography: in degree");
          prm.declare_entry ("Vs perturbation filename",
                             "dlnvs.africa.txt",
                             Patterns::FileName (),
                             "File from which the Vs perturbation data is read.");
          prm.declare_entry ("Vs to density scaling", "0.25",
                             Patterns::Double (0),
                             "This parameter specifies how the perturbation in shear wave velocity "
                             "as prescribed by SAVANI is scaled into a density perturbation. "
                             "See the general description of this model for more detailed information.");
          prm.declare_entry ("Thermal expansion coefficient in initial temperature scaling", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Vertical resolution", "25000",
                             Patterns::Double (0),
                             "resolution along depth of the 3D tomography: $m$");
          prm.declare_entry ("Reference temperature", "1600.0",
                             Patterns::Double (0),
                             "The reference temperature that is perturbed by the Vs perturbation"
                             ".Only used in incompressible models.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    RegionalModel<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Regional model");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
          perturbation_file_name = prm.get("Vs perturbation filename");
          h_delta               = prm.get_double ("Horizontal resolution");
          v_delta               = prm.get_double ("Vertical resolution");
          dlnvs_to_dlnrho          = prm.get_double ("Vs to density scaling");
          thermal_alpha          = prm.get_double ("Thermal expansion coefficient in initial temperature scaling");
          reference_temperature   = prm.get_double ("Reference temperature");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }

}


namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL (RegionalModel,
                                        "regional model",
                                        "An initial temperature condition that allows to combine Global and regional tomography  "
                                        "Both models use the same scaling technics")
  }
}




