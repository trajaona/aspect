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

#include <aspect/initial_temperature/regional_tomography.h>
#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {

   template <int dim>
   RegionalTomography<dim>::RegionalTomography()
   :
   surface_boundary_id(5)
   {}
  
   template <int dim>
   void
   RegionalTomography<dim>::initialize ()
   {
    // initialize AsciiDataInitial to read tomography file
    Utilities::AsciiDataInitial<dim>::initialize(1);
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
       return  Utilities::AsciiDataInitial<dim>::get_data_component(position,0);
   }

   template <int dim>
   double
   RegionalTomography<dim>::continental_geotherm_method1 (const double &z,
                                                          const double &lithosphere_bottom,
                                                          const double &upper_crust_bottom,
                                                          const double &middle_crust_bottom,
                                                          const double &lower_crust_bottom) const
   {
      std::vector<double> dz(5);
      std::vector<double> Q(5);
      std::vector<double> A(4);
      std::vector<double> T(5);
      double k = thermal_conductivity_of_lithosphere;

     // Define surface heat flow for different function of regions.
      if (lithosphere_bottom < 99000.0)
      {
         Q[0] = rift_surface_heat_flow; 
         A = {0.0,0.4e-6,0.4e-6,0.1e-6};
      }
      else if (lithosphere_bottom > 140000.0)
      {
         Q[0] = craton_surface_heat_flow;
         A = {0.0,0.4e-6,0.2e-6,0.1e-6}; 
      }
      else
      { 
         Q[0] = transition_surface_heat_flow; 
         A = {0.0,0.4e-6,0.4e-6,0.1e-6};
      }
      
      // Surface temperature 
      T[0] = surface_temperature;
      
      // LAB isotherm temperature
      T[4] = lab_isotherm_temperature; 

      // Heat productions for each layer. 
     A = heat_productions;
      
      // Calculate layers thicknesses 
      dz[0] = upper_crust_bottom;
      dz[1] = middle_crust_bottom - upper_crust_bottom;
      dz[2] = lower_crust_bottom  - middle_crust_bottom;
      dz[3] = lithosphere_bottom  - lower_crust_bottom;
      dz[4] = lithosphere_bottom;

      // Calculate radiogenic heat production in the upper crust
       A[0] = ( ( 2*k ) / ( dz[0] * ( 2*dz[4] - dz[0] ) ) )
              * ( T[0] - T[4] + ( ( Q[0]*dz[4] - A[1]*dz[1] * ( dz[2]+dz[3] )- A[2]*dz[2]*dz[3] ) / k ) 
              - ( A[1]*dz[1]*dz[1] + A[2]*dz[2]*dz[2] + A[3]*dz[3]*dz[3]  ) / ( 2*k ) );

      // Surface heat flows at boudaries between layers
       Q[1] = Q[0] - A[0]*dz[0];
       Q[2] = Q[1] - A[1]*dz[1];
       Q[3] = Q[2] - A[2]*dz[2];
       Q[4] = Q[3] - A[3]*dz[3];

      // Calculate temperatures at top of layers
       T[1] = T[0] + (Q[0]/k)*dz[0] - (A[0]*dz[0]*dz[0])/(2.*k);
       T[2] = T[1] + (Q[1]/k)*dz[1] - (A[1]*dz[1]*dz[1])/(2.*k);
       T[3] = T[2] + (Q[2]/k)*dz[2] - (A[2]*dz[2]*dz[2])/(2.*k);
    
    double temp;
     if  ( z  <= upper_crust_bottom )
        temp = T[0] + (Q[0]/k)*z - (A[0]*(z*z))/(2*k);
     else if ( z > upper_crust_bottom && z <= middle_crust_bottom )
        temp  = T[1] + (Q[1]/k)*(z - dz[0]) - (A[1]*((z - dz[0])*(z - dz[0])))/(2*k);
     else if ( z > middle_crust_bottom && z <= lower_crust_bottom )
        temp  = T[2] + (Q[2]/k)*(z-dz[0]-dz[1]) - (A[2]*((z-dz[0]-dz[1])*(z-dz[0]-dz[1])))/(2*k);
     else 
        temp  = T[3] + (Q[3]/k)*(z-dz[0]-dz[1]-dz[2]) - (A[3]*((z-dz[0]-dz[1]-dz[2])*(z-dz[0]-dz[1]-dz[2])))/(2*k);

    return temp;
   }


  template <int dim>
  double
  RegionalTomography<dim>::continental_geotherm_method2 (const double &depth,
         		                                  const double &lithosphere_bottom,
           		                                  const double &upper_crust_bottom,
   			                                  const double &middle_crust_bottom,
   			                                  const double &lower_crust_bottom) const
   {
           std::vector<double> dz(5);
           std::vector<double> Q(5);
           std::vector<double> A(4);
           std::vector<double> T(5);
           double k;
           double z = depth;
	   k = thermal_conductivity_of_lithosphere;  /* Thermal conductivity W/(m*K) */
	   Q[4] =  heat_flow_at_lab;  /* Boundaries bottom heat flow */
	   T[4] = lab_isotherm_temperature;
	   A[0] = surface_temperature;

	   
           /* Generate layers thickness withing the lithosphere. */
	   dz[0] = upper_crust_bottom;
	   dz[1] = middle_crust_bottom - upper_crust_bottom;
	   dz[2] = lower_crust_bottom  - middle_crust_bottom;
	   dz[3] = lithosphere_bottom  - lower_crust_bottom;
	   dz[4] = lithosphere_bottom;
 
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

	   /* Calculate heat flow at the top of the middle crust*/
           Q[1] = Q[2] + (A[1]*dz[1]);

	   /* Calculate radiogenic heat concentration at the top of the upper crust */
	   A[0] = (T[1] -T[0] - (Q[1]/k)*dz[0]) * (2.0*k) / (dz[0]*dz[0]);

	   /* Calculate heat flow at the top of the upper crust */
	   Q[0] = Q[1] + (A[0]*dz[0]);

       /* calculate temperature as function of depth */
	   double temp;
       if ( z  <= upper_crust_bottom )
    	  temp =  T[0] + (Q[0]/k)*z - ( A[0]*(z*z)) / (2*k);
       else if ( z > upper_crust_bottom && z <= middle_crust_bottom )
    	  temp = T[1] + (Q[1]/k)*(z - dz[0]) - (A[1]*((z - dz[0])*(z - dz[0]))) / (2*k);
       else if ( z > middle_crust_bottom && z <= lower_crust_bottom )
    	  temp =  T[2] + (Q[2]/k)*(z - dz[0]-dz[1]) - (A[2]*((z - dz[0]-dz[1])*(z - dz[0]-dz[1])))/ (2*k);
       else
    	  temp = T[3] + (Q[3]/k)*(z-dz[0]-dz[1]-dz[2]) - (A[3]*((z-dz[0]-dz[1]-dz[2])*(z-dz[0]-dz[1]-dz[2])))/(2*k);
      
       return temp;
   }
    
    template <int dim>
    double
    RegionalTomography<dim>::initial_temperature (const Point<dim> &position) const
    {
        const double depth = this->get_geometry_model().depth(position);
        double vs_perturbation = 0.0;
        std::array<double,dim> wcoord      = Utilities::Coordinates::WGS84_coordinates(position);
        const Point<2> wpoint (wcoord[1], wcoord[2]);       

       //   Read ascii data file containing depth of lithospheric layers. 
        const double rho_uc =   ascii_data_lab.get_data_component(surface_boundary_id, position, 6); 
        double isotherm_depth;
        double upper_crust_bottom;
        double middle_crust_bottom; 
        double lower_crust_bottom;
 
        if (use_uniform_crustal_thicknesses)
        {
           upper_crust_bottom  = 10000.0;
           middle_crust_bottom = 25000.0;
           lower_crust_bottom  = 40000.0; 
        }
        else
        {
           upper_crust_bottom   =    ascii_data_lab.get_data_component(surface_boundary_id, position, 3);
           middle_crust_bottom  =    ascii_data_lab.get_data_component(surface_boundary_id, position, 2);
           lower_crust_bottom   =      ascii_data_lab.get_data_component(surface_boundary_id, position, 1);
        }
        
        isotherm_depth =    ascii_data_lab.get_data_component(surface_boundary_id, position, 0);
        isotherm_depth = (isotherm_depth > 99000.0 && isotherm_depth < 120010.0)? 110000.0:isotherm_depth;
        
       isotherm_depth = (ascii_data_lab.get_data_component(surface_boundary_id, position, 7)==1)? rift_thickness  :isotherm_depth;
        
        // lithosphere. More realistic would be a half space cooling model. TO DO List. 
        double temperature_perturbation = 0.0;
        double lab_temperature;
        aspect::Utilities::tk::spline s;
        s.set_points(spline_depths, reference_Vs);
         
        if (use_reference_profile)
        {
               vs_perturbation = (ascii_grid_vs(position) - s(depth))/s(depth);
        }
        else
               vs_perturbation = ascii_grid_vs(position)*0.01;

    	 // scale the density perturbation into a temperature perturbation. Velocity perturbation make sense only above 400Km depth. Do not apply perturbation below that
    	 // Set a maximum perturbation to 300 K and minimum -300 K. DO TO: Find a more resonnable perturbation value. 
        double max_grid_depth = 410000.0;
        if (depth > isotherm_depth)
           {
               // scale the density perturbation into a temperature perturbation
              temperature_perturbation = -1./thermal_alpha * vs_to_density * vs_perturbation;
              temperature_perturbation = std::min(std::max(temperature_perturbation, -300.0),300.0);
           }
         else
           {
              // set heterogeneity to zero down to a specified depth
               temperature_perturbation = 0.0;
           }
    
           if (add_temperature_perturbation == false)
              temperature_perturbation = 0.0;
        
         double temperature;
        if (depth < isotherm_depth)
          temperature =  continental_geotherm_method1 (depth, isotherm_depth, upper_crust_bottom, middle_crust_bottom, lower_crust_bottom);
        else 
          temperature =  lab_isotherm_temperature + (depth - isotherm_depth) * 0.0005 + temperature_perturbation;


        return temperature;
    }

    template <int dim>
    void
    RegionalTomography<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
             prm.enter_subsection("LAB file");
             {
                 Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                            "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
           //                                                 "Emry.litho.aspect.input.txt");
                                                             "ears_synthetic_litho_updated.txt");
             }
            prm.leave_subsection();
  
            prm.enter_subsection("Tomography file");
             {
                Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
                                                              "Emry_shear_wave_velocity.txt");
             }
           prm.leave_subsection();

            prm.declare_entry ("Thermal expansion coefficient", "4e-5",
          	                 Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
            prm.declare_entry ("Vs to density", "0.1",
                      	     Patterns::Double (0),
                      	     "Vs to density scaling factor");
            prm.declare_entry ("Surface heat flow", "0.04",
                             Patterns::Double (0),
                             "Heat flow at the surface. "
							 "Units: W.m^-2");
            prm.declare_entry ("Surface temperature", "293",
                             Patterns::Double (0),
                             "Temperature at the surface or top boundary of the model. "
							 "Units: K");
            prm.declare_entry ("LAB isotherm temperature", "1600",
                              Patterns::Double (0),
                              "Isother temperatue at the LAB. "
							  "Units: K");
            prm.declare_entry ("Heat productions in the lithosphere", "0.4e-7",
                               Patterns::List(Patterns::Double(0)),
                               "Heat productions in the upper crust, middle crust, lower crust and mantle lithosphere. "
                               "Unites: ");
            prm.declare_entry ("Thermal conductivity of the lithosphere", "3.0",
                               Patterns::List(Patterns::Double(0)),
                               "Thermal conductivy of the lithosphere. It is used to calculate continetal geotherm. "
                               "Unites: W/(m*K)");
            prm.declare_entry ("Heat flow at LAB", "0.03",
                               Patterns::List(Patterns::Double(0)),
                               "Thermal conductivy of the lithosphere. It is used to calculate continetal geotherm. "
                               "Unites: (W/m^2)");
            prm.declare_entry ("Heat production in the lithosphere", "0.1",
                               Patterns::List(Patterns::Double(0)),
                               "List of heat production values in the upper crust, middle crust, "
                               "lower crust and mantle lithosphere respectively.  Units: W.m^3");
            prm.declare_entry ("Data dir", "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
                             Patterns::DirectoryName (),
                             "Directory where the reference profile file is located.");
            prm.declare_entry ("Reference profile filename",
                             "AK135F_AVG.txt",
                             Patterns::FileName (),
                             "File from which the isotherm depth data is read.");

            prm.declare_entry ("Use reference profile", "false",
                             Patterns::Bool (),
                             "Option to take the thermal expansion coefficient from the "
                             "material model instead of from what is specified in this "
                             "section.");
            prm.declare_entry ("Use temperature perturbation", "false",
                             Patterns::Bool (),
                             "Option to take the thermal expansion coefficient from the "
                             "material model instead of from what is specified in this "
                             "section.");
            prm.declare_entry ("Maximum rift thickness", "75000.0",
                             Patterns::Double (0),
                             "Maximum lithospheric thickness below which region is considered as a rift ");
            prm.declare_entry ("Minimum craton thickness", "150000.0",
                              Patterns::Double (0),
                             "Minimum lithospheric thickness above which region is considered as a craton ");
            prm.declare_entry ("Surface heat flow at rift regions", "0.08",
                              Patterns::Double (0),
                             "Surface heat flow at rift regions");
            prm.declare_entry ("Surface heat flow at transition", "0.06",
                              Patterns::Double (0),
                             "Surface heat flow at normal lithospheric regions or transition");
            prm.declare_entry ("Surface heat flow at cratonic regions", "0.05",
                              Patterns::Double (0),
                             "Surface heat flow at cratonic regions");
            prm.declare_entry ("Craton filename",
                              "african_plate_boundary.txt",
                              Patterns::FileName (),
                              "File from which the isotherm depth data is read."); 
            prm.declare_entry ("Use uniform crustal layers thicknesses", "true",
                              Patterns::Bool (),
                             "Option to use uniform crustal thicknesses instead of reading it from file. ");
            prm.declare_entry ("Use uniform lithosphere thickness", "true",
                              Patterns::Bool (),
                             "Option to use uniform lithospherick thicknesses instead of reading it from file. ");
            prm.declare_entry ("Add thick craton", "false",
                              Patterns::Bool (),
                             "Option to add thick craton at region define in the ascii data file. ");

        }
      prm.leave_subsection();
      }
      prm.leave_subsection();

    }

    template <int dim>
    void
    RegionalTomography<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Regional tomography");
        {
         prm.enter_subsection("LAB file");
          {
             ascii_data_lab.initialize_simulator (this->get_simulator());
            ascii_data_lab.parse_parameters(prm);
          }
          prm.leave_subsection();
          
         prm.enter_subsection("Tomography file");
          {
               Utilities::AsciiDataBase<dim>::parse_parameters(prm);
          }
          prm.leave_subsection();

  
          vs_to_density                         = prm.get_double ("Vs to density");
          thermal_alpha                         = prm.get_double ("Thermal expansion coefficient");
          lab_isotherm_temperature              = prm.get_double ("LAB isotherm temperature");
          maximum_rift_thickness                = prm.get_double ("Maximum rift thickness"); 
          minimum_craton_thickness              = prm.get_double ("Minimum craton thickness");
          rift_surface_heat_flow                = prm.get_double ("Surface heat flow at rift regions"); 
          transition_surface_heat_flow          = prm.get_double ("Surface heat flow at transition");
          craton_surface_heat_flow              = prm.get_double ("Surface heat flow at cratonic regions");
          surface_temperature                   = prm.get_double ("Surface temperature");
          thermal_conductivity_of_lithosphere   = prm.get_double ("Thermal conductivity of the lithosphere");
          heat_flow_at_lab                      = prm.get_double ("Heat flow at LAB");
          reference_profile_filename            = prm.get ("Reference profile filename");
          data_dir                             = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data dir")); 
          use_reference_profile                 = prm.get_bool ("Use reference profile");
          add_temperature_perturbation          = prm.get_bool ("Use temperature perturbation");
          plate_boundaries_file_name            = prm.get("Craton filename");
          use_thick_craton                      = prm.get_bool ("Add thick craton");
          use_uniform_crustal_thicknesses       = prm.get_bool ("Use uniform lithosphere thickness");
          use_uniform_LAB                       = prm.get_bool ("Use uniform crustal layers thicknesses");
     
          heat_productions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat production in the lithosphere"))),
                                                                     4,
                                                                     "Heat production in the lithosphere");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
     
      if (use_reference_profile)
      { 
      const std::string filename = data_dir+reference_profile_filename;
      
      /**
 *        * Read data from disk and distribute among processes
 *               */
      std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

      /**
       * Reading data lines.
       */
      double depths, ref_vs;
      while (in >> depths >> ref_vs)
        {
          spline_depths.push_back(depths);
          reference_Vs.push_back(ref_vs);
        }
  
      }
    
       /*Read craton file */
      const std::string filename = data_dir+plate_boundaries_file_name;

      /**
       * Read data from disk and distribute among processes
       */
       std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

      /**
       * Reading data lines
       */
       double longitude, latitude;

      while (in >> longitude >> latitude)
      {
          Point<2> point (longitude,latitude);
          boundaries_point_lists.push_back(point);
      }

    
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
