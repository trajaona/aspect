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


#include <aspect/initial_composition/rift.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>
#include <fstream>
#include <iostream>


namespace aspect
{
  namespace InitialComposition
  {

    template <int dim>
    Rift<dim>::Rift ()
    :
	surface_boundary_id(5)
    {}

    template <int dim>
    void
	Rift<dim>::initialize ()
    {
        // Find the boundary indicator that represents the surface
        surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

        std::set<types::boundary_id> surface_boundary_set;
        surface_boundary_set.insert(surface_boundary_id);

        // The input ascii table contains two components, the crust depth and the LAB depth
        Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,8); //change to 8
        ascii_data_topo.initialize(surface_boundary_set,1);
    }

    template <int dim>
    double
    Rift<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      const double depth = this->get_geometry_model().depth(position);
      std::array<double,dim> wcoord      = Utilities::Coordinates::WGS84_coordinates(position);
      wcoord[0] = depth;
      const Point<2> wpoint (wcoord[1], wcoord[2]);
      const double topo = ascii_data_topo.get_data_component(surface_boundary_id, position, 0);
      double oceanic_region =  Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 6);
      double rho_uc; 
      double rho_mc;
      double rho_lc;
      double base_of_upper_crust;
      double base_of_middle_crust;
      double base_of_lower_crust;
      double lab;

      if (use_uniform_density)
      {
       rho_uc = 2700.0;
       rho_mc = 2825.0;
       rho_lc = 3000.0;
      }
      else 
      {
       rho_uc =  Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 6);
       rho_mc =  Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 5);
       rho_lc =  Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 4);
      }
   
      if (use_uniform_crustal_thicknesses)
      {
       base_of_upper_crust    = 10000.0;
       base_of_middle_crust   = 25000.0;
       base_of_lower_crust    = 40000.0;
      }
      else
      {
       base_of_upper_crust    = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 3);
       base_of_middle_crust   = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 2);
       base_of_lower_crust    = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 1);
      }

        lab  = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 0);                
      
     double  base_of_lithosphere =  lab;
     base_of_lithosphere = (base_of_lithosphere > 99000.0 && base_of_lithosphere < 120010.0)? 110000.0:base_of_lithosphere;
     base_of_lithosphere = (base_of_lithosphere > 85000.0 && base_of_lithosphere < 95000.0)? 80000.0:base_of_lithosphere;
/*        double rift_width;
        double rift_thickness;
        if (wcoord[2] > 6.0 &&  wcoord[1] < 55)
        {
            rift_width = 0.0;
            rift_thickness =  70000.0;
        }
        else if (wcoord[1] > 55)
        {
            rift_width = 30.0;
            rift_thickness = 70000.0;
        }
        else
        {
            rift_width = 0.0;
            rift_thickness = 70000.0;
        }
*/
     //  double base_of_lithosphere = lab;
  /*    base_of_lithosphere = ( Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 7) ==1)?rift_thickness: base_of_lithosphere;
      base_of_lithosphere =  (111*std::abs(Utilities::signed_distance_to_polygon<2>(boundaries_point_lists, wpoint )) < rift_width) ?rift_thickness :base_of_lithosphere;
      bool western_branch  = (wcoord[1] < 32.0) ? true:false;
      base_of_lithosphere = (lab < 90000.0 && western_branch==true)?(rift_thickness + 5):base_of_lithosphere;
      
      base_of_lithosphere = (lab > 100000.0)?150000.0:base_of_lithosphere;
      base_of_lithosphere = (base_of_lithosphere > 99999.0 && base_of_lithosphere < 100005.0) ? 120000.0: base_of_lithosphere;      
      double rift  = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 0);
       
       */
      // calculate  layers thickness. 
      const double h_uc = base_of_upper_crust; 
      const double h_mc = base_of_middle_crust - base_of_upper_crust;
      const double h_lc = base_of_lower_crust -  base_of_middle_crust;
      const double h_ml = base_of_lithosphere - base_of_lower_crust;
      

      // Indexes of lithospheric compositonal fields.
      const unsigned int upper_crust_idx = this->introspection().compositional_index_for_name("upper_crust");
      const unsigned int middle_crust_idx = this->introspection().compositional_index_for_name("middle_crust");
      const unsigned int lower_crust_idx = this->introspection().compositional_index_for_name("lower_crust");
      const unsigned int mantle_lithosphere_idx = this->introspection().compositional_index_for_name("mantle_lithosphere");
      
      const unsigned int upper_crust_dens_idx = this->introspection().compositional_index_for_name("upper_crust_density");
      const unsigned int middle_crust_dens_idx = this->introspection().compositional_index_for_name("middle_crust_density");
      const unsigned int lower_crust_dens_idx = this->introspection().compositional_index_for_name("lower_crust_density");
      const unsigned int mantle_lithosphere_dens_idx = this->introspection().compositional_index_for_name("mantle_lithosphere_density");
      const unsigned int plastic_strain_idx = this->introspection().compositional_index_for_name("plastic_strain");
      //const unsigned int rifts_idx = this->introspection().compositional_index_for_name("rifts");
      
      // distance to polygone in km
      /*
      double distance_to_plate_boundary = 111.0 * std::abs(Utilities::signed_distance_to_polygon<2>(boundaries_point_lists, wpoint));
      const  Point<2> p1 (31.466, 2.420); 
      const  Point<2> p2 (33.1, 4.280);
      const std::array<Point<2>,2> l_1 = {p1,p2};
      const  Point<2> p3 (33.756, -8.545);
      const  Point<2> p4 (35.269, -7.096);
      const std::array<Point<2>,2> l_2 = {p3,p4};
      const  Point<2> p5 (35.369, -7.096);
      const  Point<2> p6 (35.389, -4.328);
      const std::array<Point<2>,2> l_3 = {p5,p6};

      double seg1 =  111.0 * std::abs(Utilities::distance_to_line(l_1, wpoint));    
      double seg2 =  111.0 * std::abs(Utilities::distance_to_line(l_2, wpoint));    
      double seg3 =  111.0 * std::abs(Utilities::distance_to_line(l_3, wpoint));    
      bool strong = (seg1 < 170.0 || seg2 < 150.0 || seg3 < 175.0) ? true:false;
       //base_of_lithosphere = 100000.0; 
       */
      // Lithospheric compositional fields
      if  (depth <= base_of_upper_crust &&  n_comp == upper_crust_idx) // uppper crust
          return 1.;
      else if (depth > base_of_upper_crust && depth <= base_of_middle_crust && n_comp == middle_crust_idx) // middle_crust
          return 1.;
      else if (depth > base_of_middle_crust && depth <= base_of_lower_crust && n_comp == lower_crust_idx) // lower_crust
          return 1.;
      else if (depth > base_of_lower_crust && depth <  base_of_lithosphere   && n_comp == mantle_lithosphere_idx)  // mantle lithosphere
          return 1;
      else if  (depth <= base_of_upper_crust &&  n_comp == upper_crust_dens_idx)// uppper crust
          return rho_uc*(topo + h_uc)/h_uc;
      else if (depth > base_of_upper_crust && depth <= base_of_middle_crust && n_comp == middle_crust_dens_idx) // middle_crust
          return rho_mc;
      else if (depth > base_of_middle_crust && depth <= base_of_lower_crust && n_comp == lower_crust_dens_idx) // lower_crust
          return rho_lc;
      else if (depth > base_of_lower_crust && depth < base_of_lithosphere  && n_comp == mantle_lithosphere_dens_idx) // lower_crust
      {
         double  w1 =  rho_uc * topo;
         double  w2 =  rho_uc * h_uc; 
         double  w3 =  rho_mc * h_mc;
         double  w4 =  rho_lc * h_lc; 
         
         double  w_crust = w1 + w2 + w3 + w4;
         double L = 150000.0;
         double rho_ref =  3300.0;
         
         double  w_mm    = (L - base_of_lithosphere)*3300.0;  
         double h_mtl = base_of_lithosphere - base_of_lower_crust;
         
         double w_ref = rho_ref * L;
         double w_mtl = w_ref - (w_crust + w_mm);
         double density_compensation =  w_mtl / h_mtl;
         return std::min(density_compensation,4316.0); 
       // return rho_ref = 3300.0;
     } 
       //else if (depth <   base_of_lithosphere   &&  base_of_lithosphere  < 90000.0 && strong == false && n_comp == plastic_strain_idx)
       else if (depth < base_of_lithosphere   &&  base_of_lithosphere  < 95000.0 && n_comp == plastic_strain_idx)
       {
        return 1.0; 
       }
      else
       return 0.;
    }


    template <int dim>
    void
    Rift<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
 
      Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                        "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
//                                			"Emry.litho.aspect.input.txt");
                                                        "ears_synthetic_litho_updated.txt");
    
           prm.enter_subsection("Rift");
           {

             prm.enter_subsection("Topography");
             {
              Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                             "$ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/",
                                                             "ears_topo.txt"); 
              }
          prm.leave_subsection();
      
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
        	                 Patterns::DirectoryName (),
        	    	         "The path to the isotherm depth data file");

          prm.declare_entry ("Plate boundary file name",
        	             "african_plate_boundary.txt",
        	             Patterns::FileName (),
        	             "File from which the isotherm depth data is read.");
          prm.declare_entry ("Isostatic compenstion depth", "100000.0",
                             Patterns::Double (0),
                             "The depth at which pressure is equal. Units: $m");
           prm.declare_entry ("Reference column density", "3195.8568",
                             Patterns::Double (0),
                             "The density of the reference column for isostatic compensation. Units: $kg/m3.");
           prm.declare_entry ("Rift thickness", "75000.0",
                             Patterns::Double (0),
                             "Lithospheric thickness below which region is considered as a rift ");
           prm.declare_entry ("Use uniform crustal thicknesses", "true",
                             Patterns::Bool (),
                             "Option to use uniform crustal thicknesses rather than ascii file. ");
           prm.declare_entry ("Use uniform crustal densities", "true",
                             Patterns::Bool (),
                             "Option to use uniform crustal thicknesses rather than ascii file. ");
           prm.declare_entry ("Use uniform lithospheric thickness", "true",
                             Patterns::Bool (),
                             "Option to use uniform crustal thicknesses rather than ascii file. ");
           prm.declare_entry ("No topography", "false",
                             Patterns::Bool (),
                             "Option to use uniform crustal thicknesses rather than ascii file. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Rift<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
    	Utilities::AsciiDataBase<dim>::parse_parameters(prm);
        prm.enter_subsection("Rift");
        {
          prm.enter_subsection("Topography");
          {
             ascii_data_topo.initialize_simulator (this->get_simulator()); 
             ascii_data_topo.parse_parameters(prm);
          }
           prm.leave_subsection();

          data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));

          plate_boundaries_file_name         = prm.get("Plate boundary file name");
          z_comp                             = prm.get_double ("Isostatic compenstion depth");          
          rho_refc                           = prm.get_double ("Reference column density");          
          rift_thickness                     = prm.get_double ("Rift thickness");
          use_uniform_crustal_thicknesses    = prm.get_bool ("Use uniform crustal thicknesses");
          use_uniform_density                = prm.get_bool ("Use uniform crustal densities");
          use_uniform_LAB                    = prm.get_bool ("Use uniform lithospheric thickness");
          no_topography                      = prm.get_bool ("No topography");
          
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      const std::string filename = data_directory+plate_boundaries_file_name;

      /**
       * Read data from disk and distribute among processes
       */
       std::istringstream in(Utilities::read_and_distribute_file_content(filename, this->get_mpi_communicator()));

      /**
       * Reading data lines.
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

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Rift,
                                              "rift",
                                              "Specify the composition in terms of an explicit formula. The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
