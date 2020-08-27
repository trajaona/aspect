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

#ifndef _aspect_initial_temperature_continental_geotherm_h
#define _aspect_initial_temperature_continental_geotherm_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;


    template <int dim>
    class ContinentalGeotherm:  public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        ContinentalGeotherm ();

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
        double thermal_alpha;
        double vs_to_density;
        double lab_isotherm_temperature;
        double surface_temperature;
        double thermal_conductivity_of_lithosphere;
        std::vector<double> heat_productions;
        double maximum_rift_thickness;
        double minimum_craton_thickness;
        double rift_surface_heat_flow;
        double transition_surface_heat_flow;
        double craton_surface_heat_flow; 

        Utilities::AsciiDataBoundary<dim> ascii_data_lab;

        types::boundary_id surface_boundary_id;
    

        /* Parameters for reference profile of shear wave velocity */
        std::vector<double> spline_depths;
        std::vector<double> reference_Vs;
        std::string reference_profile_filename;
        bool use_reference_profile;
        bool  add_temperature_perturbation;
        std::string data_directory;
        double heat_flow_at_lab; 
        std::string plate_boundaries_file_name;

        std::vector<Point<2>> boundaries_point_lists;



        /*
         * Continetal lithosphere geotherm parameters and function
         */

        /*
         * function that calculate temperature within the lithopshere as a function of depth
         */
       double continental_geotherm_method2 (const double &depth /*depth*/,
        		                    const double &lithosphere_bottom,
        		                    const double &upper_crust_bottom,
				            const double &middle_crust_bottom,
				            const double &lower_crust_bottom) const;
        
       double continental_geotherm_method1 (const bool &rift,
                                            const double &z, 
                                            const double &lithospheric_bottom,
                                            const double &upper_crust_bottom, 
                                            const double &middle_crust_bottom,
                                            const double &lower_crust_bottom) const ;

       bool   use_thick_craton;
       bool   use_uniform_crustal_thicknesses;
       bool   use_uniform_LAB;

     
    };
  }
}

#endif
