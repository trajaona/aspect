#include <aspect/initial_conditions/interface.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /** Here is class that compute initial temperature field. Below the lithosphere
     * Temperature increase adiabatically while above the lithosphere the temperature
     * is a linear interpolation from 0 degree at the surface to 1600 degree at the base of
     * the lithosphere. This initial condition apply to the EllipsoidalChunk Geometry

     template <int dim>
     class AdiabaticBoundary : public Interface<dim>
     {
       public:
    	  /** Return the initial temperature as a function of position
    	  *
    	  */
    	  virtual
		  double initial_temperature (const Point<dim> &position) const;

    	  /**
    	  * declare the parameters that this class needs
    	  */
          static
		  void
		  declare_parameters (ParameterHandler &prm);
          /**
           * Read the parameters above from the parameter file
           */
          virtual
		  void
		  parse_parameters (ParameterHandler &prm);


       private:
          std::vector<double>  latitudes_iso;
		  std::vector<double>  longitudes_iso;
		  std::vector<double>  depths_iso;
		  std::string litho_isotherm_file_name;
		  std::string line;
          double delta;
          int litho_flag;
          int number_coords_litho;


          /**
          * A function that read lithosphere isotherm file and return the value of the depth for each position
          */
          double
          get_lithosphere_isotherm (const double latitude,
                                    const double longitude) const;


          /**
           * Return latitude and longitude from cartesian x,y and z (wgs84)
           */
          std::pair<double, double>
          lat_long_from_xyz_wgs84(const Point<3> &pos) const;

          double
          /**
           * Return distance from center of the WGS84 to a point on the surface
           */
          radius_wgs84(const double &theta) const;

     };
  }
}

