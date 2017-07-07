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
#ifndef __aspect__initial_conditions_regional_model_h
#define __aspect__initial_conditions_regional_model_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <deal.II/base/std_cxx11/array.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    namespace internal
    {
      namespace Tomography
      {
        class VsPerturbationLookup;
      }
    }


    template <int dim>
    class RegionalModel : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

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

        /**
         * Initialization function. Loads data and sets up
         * pointers.
         */
        void initialize();

      private:

        std::string perturbation_file_name;
        double h_delta;
        double v_delta;
        double dlnvs_to_dlnrho;
        double  thermal_alpha;
        double reference_temperature;
        std::string data_directory;

        /*
         * Pointer to an object that reads the shear wave velocity at a
         * given point for 3D regional model.
        */
        std_cxx11::shared_ptr<internal::Tomography::VsPerturbationLookup> perturbation_lookup;
    };
  }
}


#endif


