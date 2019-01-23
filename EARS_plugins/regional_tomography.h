/*
   Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_temperature/adiabatic_boundary.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;


    template <int dim>
    class RegionalTomography:  public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        RegionalTomography ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * Return Vs as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        ascii_grid_vs (const Point<dim> &position) const;

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters that this class needs.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters above from the parameter file.
         */
        virtual
        void parse_parameters (ParameterHandler &prm);

      private:

        /*
         * Sublithospheric geotherm parameters
         */
        double max_grid_depth;
        double smoothing_length_scale;
        double thermal_alpha;
        double vs_to_density;
        double lab_isotherm_temperature;

        Utilities::AsciiDataBoundary<dim> ascii_data_lab;

        types::boundary_id surface_boundary_id;

        /*
         * Continetal lithosphere geotherm parameters and function
         */

        /*
         * function that calculate temperature within the lithopshere as a function of depth
         */
        double lithosphere_geotherm (const double &z /*depth*/,
        		                                  const double &lithosphere_thickness,
        		                                  const double &upper_crust_thickness,
									              const double &middle_crust_thickness,
									              const double &lower_crust_thickness,
									              const double &lab_temp) const;

        /*double atmospheric_thermal_gradient;
         *vector<double> heat_flux;
         *vector<double> layers_temperature;
         *vector<double> heat_productions;
         *vector<double>
         */
    };
  }
}

