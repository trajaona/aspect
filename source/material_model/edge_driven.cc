/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include "edge_driven.h"

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    namespace internal
    {
      namespace Boundary
      {

        class BoundaryLookup
        {
          public:

            BoundaryLookup (const std::string &filename,
                            const double &h,
                            const MPI_Comm &comm)
            {
              h_delta = h;
              /**
               * Read data from disk and distribute among processes
               */
              std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

              /**
               * Reading data lines.
               */
              double latitude_Vs, longitude_Vs, boundary_depth;
              while (in >> latitude_Vs >> longitude_Vs >> boundary_depth)
                {
                  latitudes_Vs.push_back(latitude_Vs);
                  longitudes_Vs.push_back(longitude_Vs);
                  /*
                   * Convert depth to meter
                   */
                  boundary_depths.push_back(boundary_depth * 1000.0);
                }
            }

            bool
            is_in_tomo_domain(const std_cxx11::array<double,3> &wgscoord) const
            {
              for (unsigned int i = 0; i <= latitudes_Vs.size();)
                {
                  if (std::fabs(wgscoord[2] - latitudes_Vs[i]) <= (h_delta) &&
                      std::fabs(wgscoord[1] - longitudes_Vs[i]) <= (h_delta))

                    return true;
                  else
                    i++;
                }
              return false;
            }

            double get_boundary_depth (const  std_cxx11::array<double,3> &wgscoord) const
            {
              for (unsigned int i = 0; i <= latitudes_Vs.size();)
                {
                  if (std::fabs(wgscoord[2] - latitudes_Vs[i]) <= h_delta &&
                      std::fabs(wgscoord[1] - longitudes_Vs[i]) <= h_delta)

                    return boundary_depths[i];
                  else
                    i++;
                }
              Assert (false, ExcInternalError());
              return 0;
            }

          private:
            std::vector<double>  latitudes_Vs;
            std::vector<double>  longitudes_Vs;
            std::vector<double>  boundary_depths;
            double h_delta;
        };
      }
    }

    template <int dim>
    void
    EdgeDriven<dim>::initialize()
    {
      litho_depth_lookup.reset(new internal::Boundary::BoundaryLookup(data_directory+isotherm_file_name,
                                                                      h_delta,
                                                                      this->get_mpi_communicator()));
    }
    
    template <int dim>
    double
    EdgeDriven<dim>::
    boundary_smoothing (const double &depth,
                        const double &boundary_depth,
                        const double &smoothing_thickness) const
    {
       return 0.5 + 0.5 * tanh((depth - boundary_depth) / smoothing_thickness);   
    }

    template <int dim>
    double
    EdgeDriven<dim>::
    diffusion_creep (const double &pressure,
                     const double &temperature) const
    {
      //Diffusion creep flow laws following Karato and Wu, 1993.
      const double shear_modulus = 80e9;
      const double b             = 0.5e-9;
      const double B_diff = std::pow(prefactor_diffusion / shear_modulus, -1) * std::pow(b / grain_size, -grain_size_exponent_diffusion );
      return 0.5 * B_diff * std::exp((activation_energie_diffusion + pressure * activation_volume_diffusion) / 
    		 (constants::gas_constant * temperature));
    }

    template <int dim>
    double
    EdgeDriven<dim>::
    dislocation_creep (const double &pressure,
                       const double &temperature,
                       const SymmetricTensor<2,dim> &strain_rate) const
    {
      const double p_disl = 0;                // Dislocation creep grain size exponent
      const double r_disl = 1.2;              // Dislocation creep water exponent (Freed et al., 2012 and ref. therein)
      const double strain_o = 1e-18;          // Set small strain rate for which limit viscosity calculation
      const double edot_ii_dislocation = second_invariant(strain_rate) + strain_o; 
      const double B_disl = std::pow(prefactor_dislocation * (std::pow(grain_size,(-1.0 * p_disl))) * 
     		            (std::pow(C_OH,(r_disl))) , 1 / stress_exponent_dislocation);
      return 0.5 * B_disl *
    		 std::pow(edot_ii_dislocation,(1 / stress_exponent_dislocation - 1)) * 
    		 std::exp((activation_energie_dislocation + pressure * activation_volume_dislocation) / 
    		 (stress_exponent_dislocation * constants::gas_constant * temperature));
    }
    
    
    template <>
    void
    EdgeDriven<2>::evaluate(const MaterialModel::MaterialModelInputs<2> &in,
                            MaterialModel::MaterialModelOutputs<2> &out) const
    {
      Assert (false, ExcNotImplemented());
      //return 0;
    }

    template <int dim>
    void
    EdgeDriven<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // evaluate the base model to get the crustal viscosity
      MaterialModel::MaterialModelOutputs<dim> crust_out = out;
      base_model->evaluate(in, crust_out);

      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const Point<3> pos = in.position[i];
          const std_cxx11::array<double,dim> wcoord      = Utilities::Coordinates::WGS84_coordinates(pos);
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double delta_temp = in.temperature[i]-reference_T;
          const double c = (in.composition[i].size()>0)
                           ?
                           std::max(0.0, in.composition[i][0])
                           :
                           0.0;
          // Specify density of the crust and density of the lithosphere
          if (crustal_region(in.position[i])==true)
              out.densities[i] = crust_density * (1 - thermal_alpha * (in.temperature[i] - reference_T));
          else
              out.densities[i] = mantle_density * (1 - thermal_alpha * (in.temperature[i] - reference_T));
         
          // The viscosity of the crust is determined by the base model.
          // The rheology of the lithosphere is governed by disolcation creep flow law 
          // and diffusion creep for of the sublithospheric mantle. 
          double s = 0.0;
          double viscosity;
          if (depth <  moho)
          {
            viscosity = crust_out.viscosities[i];
          }
          else if (depth >= moho && depth < litho_depth_lookup->get_boundary_depth(wcoord))
          {
             s = boundary_smoothing (depth, moho, smoothing_thickness);
             viscosity = (1-s) * crust_viscosity + 
                         s * dislocation_creep (pressure, temperature, in.strain_rate[i]);
          }
          else
          { 
             s = boundary_smoothing (depth, litho_depth_lookup->get_boundary_depth(wcoord), smoothing_thickness);
             viscosity = (1-s) * dislocation_creep (pressure, temperature, in.strain_rate[i]) +
                         s + diffusion_creep(pressure, temperature);
          }
          out.viscosities[i] = std::min(std::max(viscosity, min_visc), max_visc);
          out.thermal_expansion_coefficients[i] = thermal_alpha;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    double
    EdgeDriven<dim>::
    reference_viscosity () const
    {
      return crust_viscosity;
    }
    template <int dim>
    double
    EdgeDriven<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    EdgeDriven<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    EdgeDriven<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    EdgeDriven<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    bool
    EdgeDriven<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    bool
    EdgeDriven<dim>::
    is_litho_dislocaion_creep () const
    {
      return  disl_creep_litho;
    }

    template <>
    bool
    EdgeDriven<2>::crustal_region (const  Point<2> &) const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }

    template <int dim>
    bool
    EdgeDriven<dim>::crustal_region (const Point<dim> &position) const
    {
      const double depth                       = this->get_geometry_model().depth(position);
      if (depth <= moho)
        return true;
      else
        return false;
    }


    template <int dim>
    void
    EdgeDriven<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Edge driven");
        {
          prm.declare_entry ("Crust density", "2700",
                             Patterns::Double (0),
                             "Reference density $\\rho_c$. Units: $kg/m^3$.");
          prm.declare_entry  ("Reference strain rate","1.0e-15",Patterns::Double(0),
                              "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Moho", "30000.0",
                             Patterns::Double (0),
                             "Moho depth. Units: $m.");
           prm.declare_entry ("Smoothing thickness", "20000.0",
                             Patterns::Double (0),
                             "Thickness of the smoothed region. Units: $m.");
          prm.declare_entry ("Reference density", "3400",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Mantle lithosphere density", "3400",
                             Patterns::Double (0),
                             "Density of the lithosheric mantle $\\rho_ML$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Crust viscosity", "1e25",
                             Patterns::Double (0),
                             "The viscosity of the crust which can be  modified later by a visco-plastic rheology "
                             "modified by only temperature. Units: $kg/m/s$.");
          prm.declare_entry ("Base model","simpler",
                             Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                             "The name of a material model that will be used for the crustal viscosity. "
                             "Valid values for this parameter "
                             "are the names of models that are also valid for the "
                             "``Material models/Model name'' parameter. See the documentation for "
                             "that for more information.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $C_p$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Density differential for compositional field 1", "0",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                             "change in composition}$.");
          prm.declare_entry ("Denisity adiabatic condition","true",
                             Patterns::Bool (),
                             "Use adiabatic condition for density");
          prm.declare_entry ("Dislocation creep for lithosphere","true",
                             Patterns::Bool (),
                             "Use dislocation creep flow law for the lithospheric region");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-temperature/adiabatic-boundary/",
                             Patterns::DirectoryName (),
                             "The path to the isotherm depth data file");
          prm.declare_entry ("Isotherm depth filename",
                             "litho.africa.txt",
                             Patterns::FileName (),
                             "File from which the isotherm depth data is read.");
          prm.declare_entry ("Horizontal resolution", "3.0",
                             Patterns::Double (0),
                             "horizontal resolution of the 3D tomography: in degree");
          //rheology parameters
          prm.declare_entry ("Grain size", "1e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Minimum strain rate", "1.4e-18", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Maximum viscosity", "1e34", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("LAB depth", "150000",
                             Patterns::Double (0),
                             "Lithosphere Asthenosphere boundary depth");
          //difusion creep parameters
          prm.declare_entry ("Activation energie for diffusion creep", "300e3",
                             Patterns::List(Patterns::Double(0)),
                             "Aactivation energies, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Activation volume for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponent for diffusion creep", "1",
                             Patterns::List(Patterns::Double(0)),
                             "$n_dislocation$, for background mantle and compositional fields, "
                             "Units: None");
          prm.declare_entry ("Grain size exponent for diffusion creep", "2.5",
                             Patterns::List(Patterns::Double(0)),
                             "$m_diffusion$, for background mantle "
                             "Units: None");
          prm.declare_entry ("Prefactor for diffusion creep", "0.11e-15",
                             Patterns::List(Patterns::Double(0)),
                             "Viscosity prefactor, $A$, for background mantle,  "
                             "Units: $Pa^{-n_{diffusion}} m^{n_{diffusion}/m_{diffusion}} s^{-1}$");
          
          //dislocation creep parameters
          prm.declare_entry ("Activation energie for dislocation creep", "430e3",
                             Patterns::List(Patterns::Double(0)),
                             "Aactivation energies, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Activation volume for dislocation creep", "10e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponent for dislocation creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "Stress exponent, $n_dislocation$, for background mantle, "
                             "Units: None");
          prm.declare_entry ("Prefactor for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0)),
                             "Viscosity prefactor, $A$, for background mantle, "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Water content", "1000",
                             Patterns::List(Patterns::Double(0)),
                             "Water conten for dry olivine, $E_a$, for background mantle Units: $J / mol$");
          
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    EdgeDriven<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Edge driven");
        {
          data_directory                  = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
          isotherm_file_name              = prm.get("Isotherm depth filename");
          moho                            = prm.get_double ("Moho");
          smoothing_thickness             = prm.get_double ("Smoothing thickness");
          crust_density                   = prm.get_double ("Crust density");
          reference_rho                   = prm.get_double ("Reference density");
          mantle_density                  = prm.get_double ("Mantle lithosphere density");
          reference_T                     = prm.get_double ("Reference temperature");
          crust_viscosity                 = prm.get_double ("Crust viscosity");
          // Crustal viscosity material model
          AssertThrow( prm.get("Base model") != "edge driven",
                      ExcMessage("You may not use ``edge driven'' as the base model for "
                                  "a this material model.") );
          // create the crustal viscosity model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
             sim->initialize_simulator (this->get_simulator());

          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent      = prm.get_double ("Thermal viscosity exponent");
          k_value                         = prm.get_double ("Thermal conductivity");
          reference_specific_heat         = prm.get_double ("Reference specific heat");
          thermal_alpha                   = prm.get_double ("Thermal expansion coefficient");
          compositional_delta_rho         = prm.get_double ("Density differential for compositional field 1");
          density_no_temperature          = prm.get_bool ("Denisity adiabatic condition");
          disl_creep_litho                = prm.get_bool ("Dislocation creep for lithosphere");
          h_delta                         = prm.get_double ("Horizontal resolution");
          //rheology parameters
          grain_size                      = prm.get_double("Grain size");
          min_strain_rate                 = prm.get_double("Minimum strain rate");
          max_visc                        = prm.get_double ("Maximum viscosity");
          min_visc                        = prm.get_double ("Minimum viscosity");

          //diffusion creep parameters
          activation_energie_diffusion    = prm.get_double ("Activation energie for diffusion creep");
          activation_volume_diffusion     = prm.get_double ("Activation volume for diffusion creep");
          stress_exponent_diffusion       = prm.get_double ("Stress exponent for diffusion creep");
          grain_size_exponent_diffusion   = prm.get_double ("Grain size exponent for diffusion creep");
          prefactor_diffusion             = prm.get_double ("Prefactor for diffusion creep");
          C_OH                            = prm.get_double ("Water content");
          
          //diffusion creep parameters;
          activation_energie_dislocation  = prm.get_double ("Activation energie for dislocation creep");
          activation_volume_dislocation   = prm.get_double ("Activation volume for dislocation creep");
          stress_exponent_dislocation     = prm.get_double ("Stress exponent for dislocation creep");
          prefactor_dislocation           = prm.get_double ("Prefactor for dislocation creep");
          C_OH                            = prm.get_double ("Water content");
          lab_depth                       = prm.get_double("LAB depth");
          
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::none;

      this->model_dependence.viscosity |= NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate;
      if (composition_viscosity_prefactor != 1.0)
        this->model_dependence.viscosity |= NonlinearDependence::compositional_fields;

      if (thermal_alpha != 0)
        this->model_dependence.density |=NonlinearDependence::temperature;
      if (compositional_delta_rho != 0)
        this->model_dependence.density |=NonlinearDependence::compositional_fields;

      // Parse the parameters of the base model
      base_model->parse_parameters(prm);
      this->model_dependence.viscosity = base_model->get_model_dependence().viscosity;
    }
  }
}
// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EdgeDriven,
                                   "edge driven",
                                   "A material model for edge driven convection model. There reference density, "
                                   "of the crust is 2700 $kg/m^3$ and 3400 $kg/m^3$ for the mantle. "
                                   "The rheology of the lithosphere and the sublithospheric mantle "
			           "folows dislocation creep and diffusion creep flow laws respectively (Karato and Wu, 1993),"
	                           "We apply a smoothing function at the interface of the lithosphere and sublithospheric mantle.")
  }
}
