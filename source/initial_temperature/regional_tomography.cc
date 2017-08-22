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
#include <aspect/initial_temperature/regional_tomography.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {
    namespace internal
    {
      namespace Tomography
      {

        class ShearWaveVelocityLookup
        {
          public:


            ShearWaveVelocityLookup (const std::string &filename,
                                     const double &h,
                                     const double &v,
                                     const MPI_Comm &comm)
            {
              h_delta = h;
              v_delta = v;
              /**
               * Read data from disk and distribute among processes
               */
              std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

              /**
               * Reading data lines.
               */
              double latitude_Vs, longitude_Vs, depth_Vs, shear_wave_velocity ;
              while (in >> latitude_Vs >> longitude_Vs >> depth_Vs >> shear_wave_velocity)
                {
                  latitudes_Vs.push_back(latitude_Vs);
                  longitudes_Vs.push_back(longitude_Vs);
                  depths_Vs.push_back(depth_Vs*1000.);
                  shear_wave_velocities.push_back(shear_wave_velocity);
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

            double get_shear_wave_velocity (const  std_cxx11::array<double,3> &wgscoord) const
            {
              for (unsigned int i = 0; i <= latitudes_Vs.size();)
                {
                  if (std::fabs(wgscoord[2] - latitudes_Vs[i]) <= h_delta/2 &&
                      std::fabs(wgscoord[1] - longitudes_Vs[i]) <= h_delta/2 &&
                      std::fabs(wgscoord[0] - depths_Vs[i])  <= v_delta/2)

                    return shear_wave_velocities[i];
                  else
                    i++;
                }
              std::cout<<wgscoord[0]<<"   "<<wgscoord[1]<<"   "<<wgscoord[2]<<std::endl;
              Assert (false, ExcInternalError());
              return 0;
            }

          private:
            std::vector<double>  latitudes_Vs;
            std::vector<double>  longitudes_Vs;
            std::vector<double>  depths_Vs;
            std::vector<double>  shear_wave_velocities;
            double h_delta;
            double v_delta;
        };
      }
    }


    template <int dim>
    void
    RegionalTomography<dim>::initialize()
    {
      Vs_look_up.reset(new internal::Tomography::ShearWaveVelocityLookup(data_directory+file_name,
                                                                         h_delta,
                                                                         v_delta,
                                                                         this->get_mpi_communicator()));
    }


    template <int dim>
    double
    RegionalTomography<dim>::vs_to_temperature (const double &pressure,
                                                const double &shear_wave_velocity,
                                                const double &depth) const
    {
      // Calculate Vs*: removed non-active part of the pressure dependence of Vs.
      const double shear_wave_velocity_star = shear_wave_velocity / (1.0 + bv * (depth));
      //std::cout<<depth<<"  "<<shear_wave_velocity<<"  "<<shear_wave_velocity_star<<std::endl;
      // Initial temperature for Newton iteration*
      double temperature =  (shear_wave_velocity_star > 4400)
                            ?
                            ((shear_wave_velocity_star - c)/m)
                            : 1000.0;
      double shear_wave_velocity_residual = tolerance+0.001;
      double shear_wave_velocity_deriv;
      double temperature_iteration = 0.0;
      int it = 0;
      while (fabs(shear_wave_velocity_residual) > tolerance && it <= maximum_iteration)
        {
          shear_wave_velocity_residual = - shear_wave_velocity_star  +
                                         m * temperature +
                                         c + frequency_factor *
                                         exp (- (activation_energy + pressure * activation_volume) /
                                              (gaz_constant * (temperature +273.15)) );

          shear_wave_velocity_deriv    =  m +
                                          frequency_factor *
                                          ((activation_energy + pressure * activation_volume) /
                                           (gaz_constant * (temperature+273.15) * (temperature+273.15)) *
                                           exp (- (activation_energy + pressure * activation_volume) /
                                                (gaz_constant * (temperature+273.15) )));

          temperature  -= shear_wave_velocity_residual / shear_wave_velocity_deriv;
          it++;
        }
      return temperature;
    }


    template <>
    double
    RegionalTomography<2>::initial_temperature (const Point<2> &) const
    {
      AssertThrow (false, ExcMessage ("The 'regional tomography' initial temperature plugin is only"
                                      "implemented for 3d cases."));
      return 0;
    }

    template <int dim>
    double
    RegionalTomography<dim>::initial_temperature (const Point<dim> &position) const
    {
      const double pressure                = this->get_adiabatic_conditions().pressure(position);
      std_cxx11::array<double,dim> wcoord  = Utilities::Coordinates::WGS84_coordinates(position);
      const double depth                   = this->get_geometry_model().depth(position);
      double temp;
      const double depth_min = 75000;
      const double depth_max = 350000;
      if (depth < depth_min)
        {
          wcoord[0] = depth_min;
          const double a = (vs_to_temperature(pressure, Vs_look_up->get_shear_wave_velocity(wcoord), depth_min) )/ depth_min;
          const double b = 300.0; //in Kelvin
          temp = a * depth + b;
        }
      else if (depth > depth_max)
      {
    	  wcoord[0] = depth_max;
    	  temp = (vs_to_temperature(pressure, Vs_look_up->get_shear_wave_velocity(wcoord), depth_max)) + (depth_max/depth)*0.0005;
      }
      else
    	  wcoord[0] = depth;
    	  temp = vs_to_temperature(pressure, Vs_look_up->get_shear_wave_velocity(wcoord), depth);
        return temp;
    }

    template <int dim>
    void
    RegionalTomography<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-conditions/adiabatic-boundary/",
                             Patterns::DirectoryName (),
                             "The path to the isotherm depth data file");
          prm.declare_entry ("Tomography data filename",
                             "adiabatic_boundary.txt",
                             Patterns::FileName (),
                             "File from which the isotherm depth data is read.");
          prm.declare_entry ("Tomography vertical resolution", "20000.0",
                             Patterns::Double (0),
                             "resolution along depth of the 3D tomography: $m$");
          prm.declare_entry ("Tomography horizontal resolution", "3.0",
                             Patterns::Double (0),
                             "resolution along depth of the 3D tomography: $degree$");
          prm.declare_entry ("Constant one", "-2.8e-1",
                             Patterns::Double (),
                             "Linear slope constant of Vs to temperature relationship");
          prm.declare_entry ("Constant two", "4700.00",
                             Patterns::Double (0),
                             "Linear intercept of Vs to temperature relationship");
          prm.declare_entry ("Constant three", "3.84e-7",
                             Patterns::Double (0),
                             "Empirical constant to remove non-active part of the pressure on Vs");
          prm.declare_entry ("Frequency factor", "-1.8e16",
                             Patterns::Double (),
                             "Exponential factor  of Vs to temperature relation");
          prm.declare_entry ("Activation Energy", "409e3",
                             Patterns::Double (0),
                             "The activation energy");
          prm.declare_entry ("Activation volume", "10e-6",
                             Patterns::Double (0),
                             "Activation volume of Vs to temperature relation");
          prm.declare_entry ("Gaz constant", "8.3144",
                             Patterns::Double (0),
                             "Gaz constant");
          prm.declare_entry ("Error tolerance", "0.1",
                             Patterns::Double (0),
                             "Newton-Rapthson Error tolerance");
          prm.declare_entry ("Maximum iteration", "200",
                             Patterns::Integer (),
                             "Maximum iteration allowed for Newton-Rapthson method to solve for temperature");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    RegionalTomography<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
          data_directory              = Utilities::replace_in_string(prm.get ("Data directory"),
                                                                     "$ASPECT_SOURCE_DIR",
                                                                     ASPECT_SOURCE_DIR);
          file_name                   = prm.get("Tomography data filename");
          m                           = prm.get_double ("Constant one");
          bv                          = prm.get_double ("Constant three");
          c                           = prm.get_double ("Constant two");
          frequency_factor            = prm.get_double ("Frequency factor");
          activation_energy           = prm.get_double ("Activation Energy");
          activation_volume           = prm.get_double ("Activation volume");
          gaz_constant                = prm.get_double ("Gaz constant");
          tolerance                   = prm.get_double ("Error tolerance");
          maximum_iteration           = prm.get_integer ("Maximum iteration");
          h_delta                     = prm.get_double ("Tomography horizontal resolution");
          v_delta                     = prm.get_double ("Tomography vertical resolution");
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
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL (RegionalTomography,
                                        "regional tomography",
                                        "An initial temperature condition that allows to convert 3D Vsv data into  "
                                        "temperature field.")
  }
}




