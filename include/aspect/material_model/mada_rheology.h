/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__mada_rheology_h
#define __aspect__mada_rheology_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_temperature/adiabatic_boundary.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /*
    * A material model that consists of globally constant values for all
    * material parameters except density and viscosity.
    *
    * The model is considered incompressible, following the definition
    * described in Interface::is_compressible. This is essentially the
    * material model used in the step-32 tutorial program.
    *
    * @ingroup MaterialModels
    */
    template <int dim>
    class MadaRheology : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

        virtual bool is_non_linear() const;

//TODO: should we make this a virtual function as well? where is it used?
        double reference_thermal_diffusivity () const;

        double reference_cp () const;
        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


      private:
        std::vector<double>  latitudes_iso;
        std::vector<double>  longitudes_iso;
        std::vector<double>  depths_iso;
        std::string crustal_file_name;
        std::string data_directory;
        double h_delta;



        double mantle_density;
        double bottom_rho;
        double crust_density;
        double reference_rho;
        double reference_T;
        double eta;
        double composition_viscosity_prefactor;
        double thermal_viscosity_exponent;
        double thermal_alpha;
        double reference_specific_heat;
        bool density_no_temperature;
        bool non_linear;
        /*
         * Rheology parameters
         */
        double grain_size;
        double activation_energie_diffusion;
        double activation_volume_diffusion;
        double stress_exponent_diffusion;
        double grain_size_exponent_diffusion;
        double activation_energie_dislocation;
        double activation_volume_dislocation;
        double stress_exponent_dislocation;
        double prefactor_diffusion;
        double prefactor_dislocation;
        double min_strain_rate;
        double min_visc;
        double max_visc;
        double C_OH;
        double  lab_depth;


        /*
         * Parameters for file reading
         */
        std::string crustal_file;
        std::vector<double>  latitudes_crust;
        std::vector<double>  longitudes_crust;
        std::vector<double>  crustal_depths;

        /*
         * The Moho data resolution.
         */
        double delta;

        /**
         * The thermal conductivity.
         */
        double k_value;

        double compositional_delta_rho;

        double diffusion_creep (const double &pressure,
                                const double &temperature) const;
        double dislocation_creep (const double &pressure,
                                  const double &temperature,
                                  const SymmetricTensor<2,dim> &strain_rate) const;

        /*
         * Return true if the given coordinates is located in the crust.
         * Return false otherwise.
         */
        bool crustal_region (const Point<dim> &position) const;
    };

  }
}

#endif
