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


#include <aspect/postprocess/visualization/finite_strain.h>
#include <boost/numeric/ublas/io.hpp>
#include </Users/trajaona/software/eigen/Eigen/Eigenvalues>
#include <aspect/global.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/std_cxx11/array.h>
#include <cstddef>
#include <cmath>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      FiniteStrain<dim>::
      FiniteStrain ()
        :
        DataPostprocessorScalar<dim> ("finite_strain",
                                      update_gradients | update_q_points)
      {}



      template <int dim>
      void
      FiniteStrain<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

        Eigen::MatrixXf stretching_tensor (dim, dim);
        Eigen::MatrixXf eigen_values (dim, dim);
        Eigen::MatrixXf eigen_vectors (dim, dim);
        boost::numeric::ublas::matrix<double> strain(dim,dim);
        boost::numeric::ublas::matrix<double> L;
        boost::numeric::ublas::matrix<double> Q;
        std::vector<std::vector<double>> principal_stretch (n_quadrature_points, std::vector<double>(dim));
        std::vector<double>  lambda(dim);
        std::vector<Vector<double> >
        current_point_values (evaluation_points.size(),
                              Vector<double> (this->introspection().n_components));


        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // Get velocity gradient at the Q points.
            Tensor<2,dim> velocity_gradients;
            for (unsigned int d=0; d<dim; ++d)
              velocity_gradients[d] = input_data.solution_gradients[q][d];

            for (unsigned i = 0; i < strain.size1 (); ++ i)
                   for (unsigned j = 0; j < strain.size2 (); ++ j)
                	   strain (i, j) = velocity_gradients[i][j];

            polar::polar_decomposition(strain, Q, L);

            for (unsigned i = 0; i < L.size1 (); ++ i)
                   for (unsigned j = 0; j < L.size2 (); ++ j)
                        stretching_tensor (i, j) = L (i, j);

            Eigen::EigenSolver<Eigen::MatrixXf> es (stretching_tensor);
            eigen_values = es.pseudoEigenvalueMatrix();
            eigen_vectors = es.pseudoEigenvectors();

            for (unsigned int k = 0; k < dim; k++)
            	lambda[k] = std::fabs(eigen_values(k,k));

            std::vector<unsigned int> idx(dim);
            idx = aspect::Utilities::get_sorted_indexes(lambda);

            for (unsigned int d = 0; d < dim; d++)
            	principal_stretch[q][d] = eigen_vectors(d, idx.front());


            //std::cout<<"Original matrix" << std::endl;
            //std::cout<< stretching_tensor << std::endl;
            //std::cout<< "Eigen values"<< std::endl;
            //std::cout<< eigen_values << std::endl;
            //std::cout<< "Eigen vectors"<< std::endl;
            //std::cout<< eigen_vectors << std::endl;
            //std::cout<< " Index"<< std::endl;
            //std::cout<< idx[0] << idx[1] << idx[2] << std::endl;
            //std::cout<<"Principale strain"<<std::endl;
            //std::cout<< maximun_principale_strain[q] <<std::endl;

            computed_quantities[q](0) = std::log(lambda[idx.front()] / lambda[idx.back()]);
          }


        for (unsigned int p=0; p<evaluation_points.size(); ++p)
		{
		  // try to evaluate the solution at this point. in parallel, the point
		  // will be on only one processor's owned cells, so the others are
		  // going to throw an exception. make sure at least one processor
		  // finds the given point
		  bool point_found = false;
		  try
			{
			  VectorTools::point_value(this->get_mapping(),
									   this->get_dof_handler(),
									   principal_stretch,
									   evaluation_points[p],
									   current_point_values[p]);
			  point_found = true;
			}
		  catch (const VectorTools::ExcPointNotAvailableHere &)
			{
			  // ignore
			}

		  // ensure that at least one processor found things
		  const int n_procs = Utilities::MPI::sum (point_found ? 1 : 0, this->get_mpi_communicator());
		  AssertThrow (n_procs > 0,
					   ExcMessage ("While trying to evaluate the solution at point " +
								   Utilities::to_string(evaluation_points[p][0]) + ", " +
								   Utilities::to_string(evaluation_points[p][1]) +
								   (dim == 3
									?
									", " + Utilities::to_string(evaluation_points[p][2])
									:
									"") + "), " +
								   "no processors reported that the point lies inside the " +
								   "set of cells they own. Are you trying to evaluate the " +
								   "solution at a point that lies outside of the domain?"
								  ));

		  // Reduce all collected values into local Vector
		  Utilities::MPI::sum (current_point_values[p], this->get_mpi_communicator(),
							   current_point_values[p]);

		  // Normalize in cases where points are claimed by multiple processors
		  if (n_procs > 1)
			current_point_values[p] /= n_procs;
		}

      }

      template <int dim>
      void
      FiniteStrain<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Point values");
          {
            prm.declare_entry("Evaluation points", "",
                              // a list of points, separated by semicolons; each point has
                              // exactly 'dim' components/coordinates, separated by commas
                              Patterns::List (Patterns::List (Patterns::Double(), dim, dim, ","),
                                              0, Patterns::List::max_int_value, ";"),
                              "The list of points at which the solution should be evaluated. "
                              "Points need to be separated by semicolons, and coordinates of "
                              "each point need to be separated by commas.");

            prm.declare_entry ("Coordinate system", "cartesian",
                               Patterns::Selection ("cartesian|spherical|depth|ellipsoidal"),
                               "A selection that determines the assumed coordinate "
                               "system for the function variables. Allowed values "
                               "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                               "are interpreted as r,phi or r,phi,theta in 2D/3D "
                               "respectively with theta being the polar angle. `depth' "
                               "will create a function, in which only the first "
                               "parameter is non-zero, which is interpreted to "
                               "be the depth of the point.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      FiniteStrain<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Point values");
          {
            const std::vector<std::string> point_list
              = Utilities::split_string_list(prm.get("Evaluation points"), ';');
            for (unsigned int p=0; p<point_list.size(); ++p)
              {
                const std::vector<std::string> coordinates
                  = Utilities::split_string_list(point_list[p], ',');
                AssertThrow (coordinates.size() == dim,
                             ExcMessage ("In setting up the list of evaluation points for the <Point values> "
                                         "postprocessor, one of the evaluation points reads <"
                                         + point_list[p] +
                                         ">, but this does not correspond to a list of numbers with "
                                         "as many coordinates as you run your simulation in."));
                coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));

                Point<dim> point;
                if (coordinate_system ==  Utilities::Coordinates::CoordinateSystem::cartesian)
                {
              	  for (unsigned int d=0; d<dim; ++d)
              	  point[d] = Utilities::string_to_double (coordinates[d]);
                    evaluation_points.push_back (point);
                }
                else if (coordinate_system ==  Utilities::Coordinates::CoordinateSystem::ellipsoidal)
                {
              	  const double semi_major_axis_a = 6378137.0;
              	  const double eccentricity = 8.1819190842622e-2;
              	  const double degree_to_radian = numbers::PI/180.;
              	  std_cxx11::array<double,dim> phi_theta_d;
              	  for (unsigned int d=0; d<dim; ++d)
              	  phi_theta_d[d] = Utilities::string_to_double (coordinates[d]);
              	  phi_theta_d[0] *= degree_to_radian;
              	  phi_theta_d[1] *= degree_to_radian;
              	  point = Utilities::Coordinates::ellipsoidal_to_cartesian_coordinates<dim>(phi_theta_d, semi_major_axis_a, eccentricity);
  		          evaluation_points.push_back (point);
              	  ellipsoidal_evaluation_points.push_back(phi_theta_d);
                }
                else if (coordinate_system ==  Utilities::Coordinates::CoordinateSystem::spherical)
                {
              	  const double degree_to_radian = numbers::PI/180.;
  				  std_cxx11::array<double,dim> scoord;
  				  for (unsigned int d=0; d<dim; ++d)
  				  scoord[d] = (Utilities::string_to_double (coordinates[d]));
  				  scoord[0] *= degree_to_radian;
  				  scoord[1] *= degree_to_radian;
  				  point = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(scoord);
  				  evaluation_points.push_back (point);
  				  spherical_evaluation_points.push_back(scoord);
  			  }
                else
              	  AssertThrow(false, ExcNotImplemented());

              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      template <class Archive>
      void FiniteStrain<dim>::serialize (Archive &ar, const unsigned int)
      {
        ar &evaluation_points
        & point_values;
      }


      template <int dim>
      void
      FiniteStrain<dim>::save (std::map<std::string, std::string> &status_strings) const
      {
        std::ostringstream os;
        aspect::oarchive oa (os);
        oa << (*this);

        status_strings["FiniteStrain"] = os.str();
      }


      template <int dim>
      void
      FiniteStrain<dim>::load (const std::map<std::string, std::string> &status_strings)
      {
        // see if something was saved
        if (status_strings.find("FiniteStrain") != status_strings.end())
          {
            std::istringstream is (status_strings.find("FiniteStrain")->second);
            aspect::iarchive ia (is);
            ia >> (*this);
          }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(FiniteStrain,
                                                  "finite strain",
                                                  "Calculate stretching and output azimuth of the principal axes of the stretching tensor ")
    }
  }
}
