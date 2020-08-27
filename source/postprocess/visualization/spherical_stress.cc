/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/postprocess/visualization/spherical_stress.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      void
      SphericalStress<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == SymmetricTensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // Compute the viscosity...
        this->get_material_model().evaluate(in, out);

        // ...and use it to compute the stresses
        for (unsigned int q=0; q<n_quadrature_points; ++q)
        {
            const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            const double eta = out.viscosities[q];
       
            const SymmetricTensor<2,dim> stress = 2*eta*compressible_strain_rate + in.pressure[q] * unit_symmetric_tensor<dim>();
            const SymmetricTensor<2,dim> dev =  stress - 1./3 * trace(stress) * unit_symmetric_tensor<dim>();
           
            std::array<std::array<double,dim>,dim> stress_tensor;
            stress_tensor[0][0] = dev[dev.unrolled_to_component_indices(0)];
            stress_tensor[1][1] = dev[dev.unrolled_to_component_indices(1)];
            stress_tensor[2][2] = dev[dev.unrolled_to_component_indices(2)];
            stress_tensor[0][1] = dev[dev.unrolled_to_component_indices(3)];
            stress_tensor[1][0] = dev[0][1];
            stress_tensor[0][2] = dev[dev.unrolled_to_component_indices(4)];
            stress_tensor[2][0] = stress_tensor[0][2];
            stress_tensor[1][2] = stress[stress.unrolled_to_component_indices(5)];
            stress_tensor[2][1] = stress_tensor[1][2];        

    
            std::array<double,dim> scoord = Utilities::Coordinates::cartesian_to_spherical_coordinates(input_data.evaluation_points[q]);            
            std::array<std::array<double,dim>,dim> n;
            std::array<std::array<double,dim>,dim> nT;
            std::array<std::array<double,dim>,dim> s;
            std::array<std::array<double,dim>,dim> sp;
             
            switch (dim)
            { 
            case 2:
            {
            n[0][0] = std::cos(scoord[1]);
            n[1][1] = -std::cos(scoord[1]);
            n[0][1] = std::sin(scoord[1]);
            n[1][0] = -std::sin(scoord[1]);
            break;
            }
            case 3:
            {
            n[2][0] = std::sin(scoord[2]) * std::cos(scoord[1]);
            n[2][1] = std::sin(scoord[2]) * std::sin(scoord[1]);
            n[2][2] = std::cos(scoord[2]);
            
            n[1][1] = std::cos(scoord[2]) * std::sin(scoord[1]);
            n[1][0] = std::cos(scoord[2]) * std::cos(scoord[1]);
            n[1][2] = -std::sin(scoord[2]);
            
            n[0][0] = -std::sin(scoord[1]);
            n[0][1] = std::cos(scoord[1]);
            n[0][2] = 0.0;
            break;
            }
            default:
            {
            AssertThrow(false,ExcNotImplemented());
            }
            }
           
            for (unsigned int i = 0; i <= dim-1; ++i)
                   for (unsigned int j = 0; j <= dim-1; ++j)
                         nT[j][i] = n[i][j];           
 
            for (unsigned int i = 0; i <= dim-1; ++i)
               for (unsigned int j = 0; j <= dim-1; ++j)
                   for (unsigned int k = 0; k <= dim-1; ++k)
                         s[i][j] += stress_tensor[i][k]*nT[k][j]; 

            for (unsigned int i = 0; i <= dim-1; ++i)
               for (unsigned int j = 0; j <= dim-1; ++j)
                   for (unsigned int k = 0; k <= dim-1; ++k)
                         sp[i][j] += n[i][k]*s[k][j];
            
              
               computed_quantities[q](0) = sp[0][0]; 
               computed_quantities[q](1) = sp[1][1]; 
               computed_quantities[q](2) = sp[2][2]; 
               computed_quantities[q](3) = sp[0][1]; 
               computed_quantities[q](4) = sp[0][2]; 
               computed_quantities[q](5) = sp[1][2]; 

        }
      }


      template <int dim>
      std::vector<std::string>
      SphericalStress<dim>::get_names () const
      {
        std::vector<std::string> names;
        switch (dim)
          {
            case 2:
              names.emplace_back("stress_xx");
              names.emplace_back("stress_yy");
              names.emplace_back("stress_xy");
              break;

            case 3:
              names.emplace_back("stress_ee");
              names.emplace_back("stress_nn");
              names.emplace_back("stress_uu");
              names.emplace_back("stress_en");
              names.emplace_back("stress_eu");
              names.emplace_back("stress_nu");
          break;

            default:
              Assert (false, ExcNotImplemented());
          }

        return names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      SphericalStress<dim>::get_data_component_interpretation () const
      {
        return
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          (SymmetricTensor<2,dim>::n_independent_components,
           DataComponentInterpretation::component_is_scalar);
      }



      template <int dim>
      UpdateFlags
      SphericalStress<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values | update_quadrature_points;
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SphericalStress,
                                                  "spherical stress",
                                                  "A visualization output object that generates output "
                                                  "for the 3 (in 2d) or 6 (in 3d) components of the stress "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "$2\\eta\\varepsilon(\\mathbf u)+pI$ "
                                                  "in the incompressible case and "
                                                  "$2\\eta\\left[\\varepsilon(\\mathbf u)-"
                                                  "\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I\\right]+pI$ "
                                                  "in the compressible case.")
    }
  }
}
