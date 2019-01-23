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


#include "regional_tomography.h"
#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {

   template <int dim>
   RegionalTomography<dim>::RegionalTomography()
   :
   surface_boundary_id(1)
   {}


   template <int dim>
   void
   RegionalTomography<dim>::initialize ()
   {
    Utilities::AsciiDataInitial<dim>::initialize(2);
    // Find the boundary indicator that represents the surface
    surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

    std::set<types::boundary_id> surface_boundary_set;
    surface_boundary_set.insert(surface_boundary_id);

    // The input ascii table contains two components, the crust depth and the LAB depth
    ascii_data_lab.initialize(surface_boundary_set,8);
   }


   template <int dim>
   double
   RegionalTomography<dim>::
   ascii_grid_vs (const Point<dim> &position) const
   {
    const double vs_perturbation = Utilities::AsciiDataInitial<dim>::get_data_component(position,1);
    return (vs_perturbation)*0.01;
   }


   template <int dim>
   double
   RegionalTomography<dim>::lithosphere_geotherm (const double &z,
           		                                  const double &lithosphere_bottom,
           		                                  const double &upper_crust_bottom,
   									              const double &middle_crust_bottom,
   									              const double &lower_crust_bottom,
												  const double &lab_temp) const
   {

	   std::vector<double> A{0, 0.4e-7, 0.4e-7, 0.1e-7};  /* Layers heat production */
	   std::vector<double> T(5, 0.0);  /* Boundaries bottom temperature */
	   std::vector<double> Q{0.08, 0., 0., 0., 0.};  /* Boundaries bottom heat flow */
	   double k = 3.0;  /* Thermal conductivity W/(m*K) */
	   T[4] = lab_temp;

	   /*
        * Calculate layers thickness
        */
	   std::vector<double> dz(5, 0.0); /* Layers thickness */
	   dz[0] = upper_crust_bottom;
	   dz[1] = middle_crust_bottom - upper_crust_bottom;
	   dz[2] = lower_crust_bottom  - middle_crust_bottom;
	   dz[3] = lithosphere_bottom  - lower_crust_bottom;
	   dz[4] = lithosphere_bottom;

	  /* Calculate radiogenic heat production in the upper crust */
      // A [0] = ( ( 2*k ) / ( dz[0] * ( 2*dz[4] - dz[0] ) ) ) *
      // 		   ( T[0] - T[4] + ( ( Q[0]*dz[4] - A[1]*dz[1] * ( dz[2]+dz[3] ) - A[2]*dz[2]*dz[3] ) / k ) -
      //	       ( A[1]*dz[1]*dz[1] + A[2]*dz[2]*dz[2] + A[3]*dz[3]*dz[3]) / ( 2*k ) );

	  /* Calculate surface heat flows at boundaries between layer */
      // Q[1] = Q[0] - A[0]*dz[0];
      // Q[2] = Q[1] - A[1]*dz[1];
      // Q[3] = Q[2] - A[2]*dz[2];
      // Q[4] = Q[3] - A[3]*dz[3];

      /* Calculate temperatures at top of layers */
      // T [1] = T [0] + (Q[0]/k)*dz[0] - (A[0]*dz[0]*dz[0]) / (2.0*k);
      // T [2] = T [1] + (Q[1]/k)*dz[1] - (A[1]*dz[1]*dz[1]) / (2.0*k);
      // T [3] = T [2] + (Q[2]/k)*dz[2] - (A[2]*dz[2]*dz[2]) / (2.0*k);

	   /* Calculate temperature at the top of the mantle lithosphere */
	   T[3] = T[4] - (Q[4]/k)*dz[3] - (A[3]*dz[3]*dz[3])/(2*k);

	   /* Calculate heat flow at the top of the mantle lithosphere */
	   Q[3] = Q[4] + (A[3]*dz[3]);

	   /* Calculate temperature at that the top of the lower crust */
	   T[2] = T[3] - (Q[3]/k)*dz[2] - (A[2]*dz[2]*dz[2])/(2.0*k);

	   /* Calculate heat flow at the top of the lower crust */
	   Q[2] = Q[3] + (A[2]*dz[2]);

	   /* Calculate heat flow at the top of the middle crust */
	   T[1] = T[2] - (Q[2]/k)*dz[1] - (A[1]*dz[1]*dz[1]) / (2.0*k);

	   /* Calculate radiogenic heat concentration at the top of the upper crust */
	   A[0] = (T[1] -T[0] - (Q[1]/k)*dz[0]) * (2.0*k) / (dz[0]*dz[0]);

	   /* Calculate heat flow at the top of the upper crust */
	   Q[0] = Q[1] + (A[0]*dz[0]);

       /* calculate temperature as function of depth */
	   double temperature;
       if ( z  <= upper_crust_bottom )
    	  temperature =  T[0] + (Q[0]/k)*z - ( A[0]*(z*z)) / (2*k);
       else if ( z > upper_crust_bottom && z <= middle_crust_bottom )
    	  temperature = T[1] + (Q[1]/k)*(z - dz[0]) - (A[1]*((z - dz[0])*(z - dz[0]))) / (2*k);
       else if ( z > middle_crust_bottom && z <= lower_crust_bottom )
    	  temperature =  T[2] + (Q[2]/k)*(z - dz[0]-dz[1]) - (A[2]*((z - dz[0]-dz[1])*(z - dz[0]-dz[1])))/ (2*k);
       else
    	  temperature = T[3] + (Q[3]/k)*(z-dz[0]-dz[1]-dz[2]) - (A[3]*((z-dz[0]-dz[1]-dz[2])*(z-dz[0]-dz[1]-dz[2])))/(2*k) ;

       return temperature;
   }


    template <int dim>
    double
    RegionalTomography<dim>::initial_temperature (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      std::array<double,dim> wcoord      = Utilities::Coordinates::WGS84_coordinates(position);
      const double radius = 6378137.;

      double vs_perturbation;
      if (depth <= max_grid_depth - smoothing_length_scale)
      {
    	  vs_perturbation = ascii_grid_vs(position);
      }
      //add smoothing between the two models
      else if (depth > max_grid_depth - smoothing_length_scale && depth < max_grid_depth)
      {
    	   const double scale_factor = (depth-(max_grid_depth-smoothing_length_scale))/smoothing_length_scale;
    	   vs_perturbation = 0.0 *(scale_factor) + ascii_grid_vs(position)*(1.0-scale_factor);
      }
      else
      {
    	   vs_perturbation = 0.0;
      }
      const double density_perturbation =  vs_to_density * vs_perturbation;
      const double upper_crust_bottom   =  ascii_data_lab.get_data_component(surface_boundary_id, position, 3);
      const double middle_crust_bottom  =  ascii_data_lab.get_data_component(surface_boundary_id, position, 2);
      const double lower_crust_bottom   =  ascii_data_lab.get_data_component(surface_boundary_id, position, 1);
      const double isotherm_depth       =  ascii_data_lab.get_data_component(surface_boundary_id, position, 0);

      double temperature_perturbation;
      if (depth < max_grid_depth && depth > isotherm_depth)
      {
    	  // scale the density perturbation into a temperature perturbation
    	  temperature_perturbation =  -1./thermal_alpha* density_perturbation;
      }
      else
    	  // set heterogeneity to zero down to a specified depth
    	  temperature_perturbation = 0.0;

      wcoord[0] = radius - isotherm_depth;

      Point<dim> lab_position =  Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(wcoord);;

      double lab_temperature = lab_isotherm_temperature + ascii_grid_vs(lab_position);

      double lithosphere_geotherm_temp =  lithosphere_geotherm (depth,
      		                                                    isotherm_depth,
	  													        upper_crust_bottom,
	  													        middle_crust_bottom,
	 													        lower_crust_bottom,
														        lab_temperature);
      std::cout<<upper_crust_bottom<<"  "
         	   <<middle_crust_bottom<<"  "
			   <<lower_crust_bottom<<"   "
			   <<lab_temperature<<"    "
			   <<temperature_perturbation<<std::endl;

     if (depth > isotherm_depth)
            return  lab_temperature + (depth - isotherm_depth) * 0.0005;
     else
    	    return lithosphere_geotherm_temp;
    }

    template <int dim>
    void
    RegionalTomography<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
            prm.declare_entry ("Maximum grid depth", "350000.0",
          	                 Patterns::Double (0),
          	                 "The maximum depth of the Vs ascii grid. The model will read in  "
          	                 "Vs from golbal model or use temperature from adiabatic boundary below this depth.");
            prm.declare_entry ("Thermal expansion coefficient", "2e-5",
          	                 Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
            prm.declare_entry ("Smoothing length scale", "10000.0",
          	                 Patterns::Double (0),
          	                 "The depth range (above maximum grid depth) over which to smooth. "
          	                 "The boundary is smoothed using a depth weighted combination of Vs.");
            prm.declare_entry ("Vs to density", "0.25",
                      	     Patterns::Double (0),
                      	     "Vs to density scaling factor");
            prm.declare_entry ("Surface heat flow", "0.04",
                             Patterns::Double (0),
                             "Heat flow at the surface. "
							 "Units: W.m^-2");
            prm.declare_entry ("LAB isotherm temperature", "1600",
                             Patterns::Double (0),
                             "Isother temperatue at the LAB. "
							 "Units: K");
            Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
                                                              "regional_tomography_3d.txt");
        }
      prm.leave_subsection();
      }
      prm.leave_subsection();

      Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
         	    	     	    	         	        "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
         	    	     	    	         	        "kenya.lithosphere.txt");
    }

    template <int dim>
    void
    RegionalTomography<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
          max_grid_depth           = prm.get_double ("Maximum grid depth");
          smoothing_length_scale   = prm.get_double ("Smoothing length scale");
          thermal_alpha            = prm.get_double ("Thermal expansion coefficient");
          vs_to_density            = prm.get_double ("Vs to density");
          lab_isotherm_temperature = prm.get_double ("LAB isotherm temperature");

          Utilities::AsciiDataBase<dim>::parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      ascii_data_lab.initialize_simulator (this->get_simulator());

      // Note: parse_parameters will call initialize for us
      ascii_data_lab.parse_parameters(prm);

    }

  }

}


namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(RegionalTomography,
                                              "regional tomography",
                                              "An initial temperature condition that allows for discretizing "
                                              "is designed specifically for the ellipsoidal chunk geometry model.")
  }
}
