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


#include <aspect/initial_composition/lithosphere.h>
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
    Lithosphere<dim>::Lithosphere ()
    :
	surface_boundary_id(5)
    {}

    template <int dim>
    void
	Lithosphere<dim>::initialize ()
    {
        // Find the boundary indicator that represents the surface
        surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

        std::set<types::boundary_id> surface_boundary_set;
        surface_boundary_set.insert(surface_boundary_id);

        // The input ascii table contains two components, the crust depth and the LAB depth
        Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                      7);
        ascii_data_topo.initialize(surface_boundary_set,1);
    }

    template <int dim>
    double
    Lithosphere<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      const double depth = this->get_geometry_model().depth(position);
      std::array<double,dim> wcoord      = Utilities::Coordinates::WGS84_coordinates(position);
      wcoord[0] = depth;

      const Point<2> wpoint (wcoord[1], wcoord[2]);

     // const double base_of_upper_crust              = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 3);
     // const double base_of_middle_crust             = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 2);
     // const double base_of_lower_crust              = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 1);
     //const double base_of_lithosphere              = std::max(Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 0),100000.0);
         const double base_of_lithosphere              = 100000.0;
        const double base_of_upper_crust = 20000.0;
        const double base_of_middle_crust = 30000.0;
        const double base_of_lower_crust = 35000.0; 
      // Indexes of lithospheric compositonal fields.
      const unsigned int upper_crust_idx = this->introspection().compositional_index_for_name("upper_crust");
      const unsigned int middle_crust_idx = this->introspection().compositional_index_for_name("middle_crust");
      const unsigned int lower_crust_idx = this->introspection().compositional_index_for_name("lower_crust");
      const unsigned int mantle_lithosphere_idx = this->introspection().compositional_index_for_name("mantle_lithosphere");
      
      const unsigned int upper_crust_dens_idx = this->introspection().compositional_index_for_name("upper_crust_density");
     const unsigned int middle_crust_dens_idx = this->introspection().compositional_index_for_name("middle_crust_density");
      const unsigned int lower_crust_dens_idx = this->introspection().compositional_index_for_name("lower_crust_density");
  
     // Lithospheric compositional fields
      if  (depth <= base_of_upper_crust &&  n_comp == upper_crust_idx) // uppper crust
          return 1.;
      else if (depth > base_of_upper_crust && depth <= base_of_middle_crust && n_comp == middle_crust_idx) // middle_crust
          return 1.;
      else if (depth > base_of_middle_crust && depth <= base_of_lower_crust && n_comp == lower_crust_idx) // lower_crust
          return 1.;
      else if (depth > base_of_lower_crust && depth <= base_of_lithosphere && n_comp == mantle_lithosphere_idx)  // mantle lithosphere
          return 1;
      else if  (depth <= base_of_upper_crust &&  n_comp == upper_crust_dens_idx) // uppper crust
         return Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 6);
      else if (depth > base_of_upper_crust && depth <= base_of_middle_crust && n_comp == middle_crust_dens_idx) // middle_crust
          return Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 5);
      else if (depth > base_of_middle_crust && depth <= base_of_lower_crust && n_comp == lower_crust_dens_idx) // lower_crust
          return Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 4);
       else
        return 0.;
    }


    template <int dim>
    void
    Lithosphere<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
  
    	Utilities::AsciiDataBase<dim>::declare_parameters(prm,
    	     	    	         	                  "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
       	     	    	         	                  "ears_synthetic_litho_with_tanzania_craton.txt");

        prm.enter_subsection("Lithosphere");
        {
           prm.enter_subsection("Topography");
           {
            Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                             "$ASPECT_SOURCE_DIR/data/geometry-model/initial-topography-model/ascii-data/",
                                                             "ears_topography.txt"); 
           }
          prm.leave_subsection();
      
          Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                             "$ASPECT_SOURCE_DIR/data/data/geometry-model/initial-topography-model/ascii-data/",
                                                             "ears_topography.txt");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-composition/ascii-data/",
        	                 Patterns::DirectoryName (),
        	    	         "The path to the isotherm depth data file");

          prm.declare_entry ("Craton filename",
        	                 "plate_boundaries.txt",
        	                 Patterns::FileName (),
        	                 "File from which the isotherm depth data is read.");
          prm.declare_entry ("Moho", "30000.0",
                             Patterns::Double (0),
                             "Moho depth. Units: $m.");
          prm.declare_entry ("LAB isotherm", "1673",
                             Patterns::Double (0),
                             "Temperature at the base of the lithosphere. Units: $K.");
                   prm.declare_entry ("Topography factor", "1.0",
                             Patterns::Double (0),
                             "Temperature at the base of the lithosphere. Units: $K.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Lithosphere<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
    	Utilities::AsciiDataBase<dim>::parse_parameters(prm);
        prm.enter_subsection("Lithosphere");
        {
          prm.enter_subsection("Topography");
          {
             ascii_data_topo.initialize_simulator (this->get_simulator()); 
             ascii_data_topo.parse_parameters(prm);
          }
          prm.leave_subsection();

          data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));

          plate_boundaries_file_name   = prm.get("Craton filename");
          moho                            = prm.get_double ("Moho");
          LAB_isotherm                    = prm.get_double ("LAB isotherm");
          topo_fac                    = prm.get_double ("Topography factor");
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
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Lithosphere,
                                              "lithosphere",
                                              "Specify the composition in terms of an explicit formula. The format of these "
                                              "functions follows the syntax understood by the "
                                              "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
