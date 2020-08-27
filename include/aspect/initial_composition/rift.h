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


#ifndef _aspect_initial_composition_rift_h
#define _aspect_initial_composition_rift_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Rift : public Interface<dim>,
	  public Utilities::AsciiDataBoundary<dim>
    {
      public:
        /**
         * Constructor.
         */
        Rift ();

        void
        initialize ();

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        
        double z_comp;

        double rho_refc;
        double rift_thickness;

        bool  use_uniform_crustal_thicknesses;
        bool  use_uniform_density;
        bool  use_uniform_LAB;
        bool no_topography;

        std::string plate_boundaries_file_name;
       
        std::vector<Point<2>> boundaries_point_lists;

        std::string data_directory;

        types::boundary_id surface_boundary_id;
         Utilities::AsciiDataBoundary<dim> ascii_data_topo; 
    };
  }
}


#endif
