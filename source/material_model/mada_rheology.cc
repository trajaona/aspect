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


#include <aspect/material_model/mada_rheology.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    double
    MadaRheology<dim>::
    diffusion_creep (const double &pressure,
                     const double &temperature) const
    {
      // If strain rate is zero (like during the first time step) set it to some very small number
      // to prevent a division-by-zero, and a floating point exception.
      // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
      // strain rate (often simplified as epsilondot_ii)
      //const double edot_ii_diffusion = std::max(std::fabs(second_invariant(strain_rate)),
      //                                 min_strain_rate * min_strain_rate);

      //For diffusion creep, viscosity is grain size dependent
      //double prefactor_stress_diffusion = std::exp((activation_energie_diffusion + pressure*activation_volume_diffusion)/
      //                                  (stress_exponent_diffusion*constants::gas_constant*temperature)) /
      //                                 std::pow(prefactor_diffusion *
      //                                       std::pow(grain_size, -1/grain_size_exponent_diffusion)  * C_OH, 1/stress_exponent_diffusion);

      //double strain_factor = 0.5*pow(edot_ii_diffusion, ((1 - stress_exponent_diffusion) / stress_exponent_diffusion));
      //return prefactor_stress_diffusion * strain_factor;
      const double B_diff = prefactor_diffusion;
      return 0.5 * B_diff * std::exp((activation_energie_diffusion+pressure*activation_volume_diffusion)/(constants::gas_constant*temperature));
    }

    template <int dim>
    double
    MadaRheology<dim>::
    dislocation_creep (const double &pressure,
                       const double &temperature,
                       const SymmetricTensor<2,dim> &strain_rate) const
    {
      // If strain rate is zero (like during the first time step) set it to some very small number
      // to prevent a division-by-zero, and a floating point exception.
      // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
      // strain rate (often simplified as epsilondot_ii)
      //const double edot_ii_dislocation = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
      // min_strain_rate * min_strain_rate);
      //const double r_disl = 1.2;
      // For dislocation creep, viscosity is independent
      //double prefactor_stress_dislocation = std::pow(prefactor_dislocation, -1/stress_exponent_dislocation) *
      //                                   std::pow(C_OH, -r_disl/stress_exponent_dislocation) *
      //                  std::exp((activation_energie_dislocation + pressure*activation_volume_dislocation)/
      //                                         (stress_exponent_dislocation*constants::gas_constant*temperature));
      //double strain_factor = 0.5*pow(edot_ii_dislocation, ((1 - stress_exponent_dislocation) / stress_exponent_dislocation));
      //return prefactor_stress_dislocation * strain_factor;
      const double R = 8.3144;                  // gas constant J/K.mol
      const double n_disl = 3.5;              // n for dislocation creep (Freed et al., 2012 and ref. therein)
      const double E_disl = 480000;           // Activation energy J/mol (Freed et al., 2012 and ref. therein)
      const double V_disl = 0.00000000001;      // Activation volume m^3/mol (Freed et al., 2012 and ref. therein)
      const double A_disl = 30000000;         // Material parameter Pa^-n s^-1
      const double d = 0.003;               // Grain size m
      const double p_disl = 0;              // Dislocation creep grain size exponent
      const double C_OH = 1000;             // Olivine water content H/10^6Si (Freed et al., 2012 and ref. therein)
      const double r_disl = 1.2;              // Dislocation creep water exponent (Freed et al., 2012 and ref. therein)
      const double strain_o = 1e-18;          // Set small strain rate for which limit viscosity calculation

      const double B_disl = A_disl * (std::pow(d,(-1.0*p_disl))) * (std::pow(C_OH,(r_disl)));
      const double exp_factor = std::exp((E_disl + pressure*V_disl)/(n_disl*R*temperature));
      const double n_factor = 1/n_disl;
      const double exp_B_disl = std::pow(B_disl,n_factor);
      const double secInvStrainRate = second_invariant(strain_rate) + strain_o;
      const double strain_factor = 0.5 * std::pow(secInvStrainRate,(n_factor-1));
      return strain_factor*exp_factor*exp_B_disl;
    }

    template <>
    void
    MadaRheology<2>::evaluate(const MaterialModel::MaterialModelInputs<2> &in,
                              MaterialModel::MaterialModelOutputs<2> &out) const
    {
      Assert (false, ExcNotImplemented());
      //return 0;
    }

    template <int dim>
    void
    MadaRheology<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double delta_temp = in.temperature[i]-reference_T;
          const double c = (in.composition[i].size()>0)
                           ?
                           std::max(0.0, in.composition[i][0])
                           :
                           0.0;
          // Define density
          double rho;
          if (crustal_region(in.position[i])==true)
            out.densities[i]= crust_density;
          else
            out.densities[i] = mantle_density * (1 - thermal_alpha * (in.temperature[i] - reference_T))
                               + compositional_delta_rho * c;

          
          // tanngent hyperbolic
          //  s = 0.5 + 0.5*tanh((depth - (litho_depth_lookup->get_boundary_depth(wcoord)))/30000.0);
          //Here 'eta' is the viscosity of the lithosphere
          
          // Define viscosity
          double viscosity = 0.0;
          if ( depth <= lab_depth)
            viscosity = eta;
          else
            viscosity = diffusion_creep(pressure, temperature);
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
    MadaRheology<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    MadaRheology<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    MadaRheology<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    MadaRheology<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    MadaRheology<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    bool
    MadaRheology<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    bool
    MadaRheology<dim>::
    is_non_linear () const
    {
      return non_linear;
    }

    template <>
    bool
    MadaRheology<2>::crustal_region (const  Point<2> &) const
    {
      Assert (false, ExcNotImplemented());
      return 0;
    }

    template <int dim>
    bool
    MadaRheology<dim>::crustal_region (const Point<dim> &position) const
    {
      const double depth                       = this->get_geometry_model().depth(position);
      if (depth <= 30000.0)
        return true;
      else
        return false;
    }


    template <int dim>
    void
    MadaRheology<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Mada rheology");
        {
          prm.declare_entry ("Horizontal resolution", "2.0",
                             Patterns::Double (0),
                             "Horizontal resolution of crust thickness");
          prm.declare_entry ("Crust density", "2700",
                             Patterns::Double (0),
                             "Reference density $\\rho_c$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference density", "3400",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Mantle lithosphere density", "3400",
                             Patterns::Double (0),
                             "Density of the lithosheric mantle $\\rho_ML$. Units: $kg/m^3$.");
          prm.declare_entry ("Bottom density", "4400",
                             Patterns::Double (0),
                             "Density at the base of the model $\\rho_B$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("LAB depth", "150000",
                             Patterns::Double (0),
                             "Lithosphere Asthenosphere boundary depth");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$. This viscosity may be "
                             "modified by both temperature and compositional dependencies. Units: $kg/m/s$.");
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
          prm.declare_entry ("Non linear rheology","true",
                             Patterns::Bool (),
                             "Use non linear rheology");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/initial-conditions/adiabatic-boundary/",
                             Patterns::DirectoryName (),
                             "The path to the isotherm depth data file");
          //rheology parameters
          prm.declare_entry ("Grain size", "1e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Minimum strain rate", "1.4e-18", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");

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
    MadaRheology<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Mada rheology");
        {
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
          crust_density                  = prm.get_double ("Crust density");
          reference_rho                  = prm.get_double ("Reference density");
          bottom_rho                     = prm.get_double ("Bottom density");
          mantle_density                 = prm.get_double ("Mantle lithosphere density");
          reference_T                    = prm.get_double ("Reference temperature");
          eta                            = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent      = prm.get_double ("Thermal viscosity exponent");
          k_value                         = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");
          density_no_temperature     = prm.get_bool ("Denisity adiabatic condition");
          non_linear                 = prm.get_bool ("Non linear rheology");
          h_delta                    = prm.get_double ("Horizontal resolution");
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
          activation_energie_dislocation  = prm.get_double ("Activation energie for dislocation creep");
          activation_volume_dislocation   = prm.get_double ("Activation volume for dislocation creep");
          stress_exponent_dislocation     = prm.get_double ("Stress exponent for dislocation creep");
          prefactor_dislocation           = prm.get_double ("Prefactor for dislocation creep");
          C_OH                            = prm.get_double ("Water content");
          lab_depth                      = prm.get_double("LAB depth");


          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model simple with Thermal viscosity exponent can not have reference_T=0."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
      /* After parsing the parameters for depth dependent, it is essential to parse
        parameters related to the base model. */

      // Declare dependencies on solution variables
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::none;

      this->model_dependence.viscosity |= NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate;
      if (composition_viscosity_prefactor != 1.0)
        this->model_dependence.viscosity |= NonlinearDependence::compositional_fields;

      if (thermal_alpha != 0)
        this->model_dependence.density |=NonlinearDependence::temperature;
      if (compositional_delta_rho != 0)
        this->model_dependence.density |=NonlinearDependence::compositional_fields;

    }
  }
}
// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MadaRheology,
                                   "mada rheology",
                                   "A material model that has constant values "
                                   "for all coefficients but the density and viscosity. The defaults for all "
                                   "coefficients are chosen to be similar to what is believed to be correct "
                                   "for Earth's mantle. All of the values that define this model are read "
                                   "from a section ``Material model/Simple model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Simple_20model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant: "
                                   "\\begin{align}"
                                   "  \\eta(p,T,\\mathfrak c) &= \\tau(T) \\zeta(\\mathfrak c) \\eta_0, \\\\"
                                   "  \\rho(p,T,\\mathfrak c) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0 + \\Delta\\rho \\; c_0,"
                                   "\\end{align}"
                                   "where $c_0$ is the first component of the compositional vector "
                                   "$\\mathfrak c$ if the model uses compositional fields, or zero otherwise. "
                                   "\n\n"
                                   "The temperature pre-factor for the viscosity formula above is "
                                   "defined as "
                                   "\\begin{align}"
                                   "  \\tau(T) &= H\\left(e^{-\\beta (T-T_0)/T_0}\\right),"
                                   "  \\qquad\\qquad H(x) = \\begin{cases}"
                                   "                            10^{-2} & \\text{if}\\; x<10^{-2}, \\\\"
                                   "                            x & \\text{if}\\; 10^{-2}\\le x \\le 10^2, \\\\"
                                   "                            10^{2} & \\text{if}\\; x>10^{2}, \\\\"
                                   "                         \\end{cases}"
                                   "\\end{align} "
                                   "where $\\beta$ corresponds to the input parameter ``Thermal viscosity exponent'' "
                                   "and $T_0$ to the parameter ``Reference temperature''. If you set $T_0=0$ "
                                   "in the input file, the thermal pre-factor $\\tau(T)=1$."
                                   "\n\n"
                                   "The compositional pre-factor for the viscosity is defined as "
                                   "\\begin{align}"
                                   "  \\zeta(\\mathfrak c) &= \\xi^{c_0}"
                                   "\\end{align} "
                                   "if the model has compositional fields and equals one otherwise. $\\xi$ "
                                   "corresponds to the parameter ``Composition viscosity prefactor'' in the "
                                   "input file."
                                   "\n\n"
                                   "Finally, in the formula for the density, $\\alpha$ corresponds to the "
                                   "``Thermal expansion coefficient'' and "
                                   "$\\Delta\\rho$ "
                                   "corresponds to the parameter ``Density differential for compositional field 1''."
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\alpha(T-T_{\\text{ref}}))$. "
                                   "\n\n"
                                   "\\note{Despite its name, this material model is not exactly ``simple'', "
                                   "as indicated by the formulas above. While it was originally intended "
                                   "to be simple, it has over time acquired all sorts of temperature "
                                   "and compositional dependencies that weren't initially intended. "
                                   "Consequently, there is now a ``simpler'' material model that now fills "
                                   "the role the current model was originally intended to fill.}")
  }
}
