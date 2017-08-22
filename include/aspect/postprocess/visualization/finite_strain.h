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


#ifndef _aspect_postprocess_visualization_finite_strain_h
#define _aspect_postprocess_visualization_finite_strain_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>
#include <boost/numeric/ublas/io.hpp>
#include </Users/trajaona/software/aspect/cookbooks/polar_decomposition/polar-decomposition/polar_decomposition.hpp>
#include </Users/trajaona/software/eigen/Eigen/Eigenvalues>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the strain rate at every
       * point. The scalar strain rate is defined as $\sqrt{ (\varepsilon -
       * \tfrac 13 \textrm{trace}\ \varepsilon \mathbf 1) : \varepsilon -
       * \tfrac 13 \textrm{trace}\ \varepsilon \mathbf 1}$.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class FiniteStrain
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          FiniteStrain ();

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

          /**
           * Save the state of this object.
           */
          virtual
          void save (std::map<std::string, std::string> &status_strings) const;

          /**
           * Restore the state of the object.
           */
          virtual
          void load (const std::map<std::string, std::string> &status_strings);

          /**
           * Serialize the contents of this class as far as they are not read
           * from input parameter files.
           */
          template <class Archive>
          void serialize (Archive &ar, const unsigned int version);

          virtual
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const;

        private:
		  std::vector<Point<dim> >                                       evaluation_points;
	      std::vector<std_cxx11::array<double,dim> >                     spherical_evaluation_points;
		  std::vector<std_cxx11::array<double,dim> >                     ellipsoidal_evaluation_points;
		  std::vector<std::pair<double, std::vector<Vector<double> > > > point_values;

		  /**
		   * The coordinate representation to evaluate the function. Possible
		   * choices are depth, cartesian and spherical.
		   */
		  Utilities::Coordinates::CoordinateSystem coordinate_system;

      };
    }
  }
}

#endif
