// Think about adding these header files when needed

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /** Here is class that compute initial temperature field. Below the lithosphere
     * Temperature increase adiabatically while above the lithosphere the temperature
     * is a linear interpolation from 0 degree at the surface to 1600 degree at the base of
     * the lithosphere. This initial condition apply to the EllipsoidalChunk Geometry

     // Name of the class is temporary, could be changed
      *
      */
     template <int dim>
     class GeothermalModel : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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


       /*
        * The private member of the class
        */
       private:
           /*
            * We need the read the lithosphere isotherm depth in function of coordinates
            * latitude and longitude
            */
          std::vector<double>  latitudes_iso;
		  std::vector<double>  longitudes_iso;
		  std::vector<double>  depths_iso;
		  std::string litho_isotherm_file_name;
		  std::string line;
          double delta;
          int litho_flag;
          int number_coords_litho;


          /**
          * a function that read lithosphere isotherm file and return the value of the depth for each position
          */
          double
          get_lithosphere_isotherm (const double latitude,
                                    const double longitude) const;


          /**
           * Return latitude and longitude from cartesian x,y and z (wgs84)
           */
          std::pair<double, double>
          lat_long_from_xyz_wgs84(const Point<3> &pos) const;


     };
  }
}

