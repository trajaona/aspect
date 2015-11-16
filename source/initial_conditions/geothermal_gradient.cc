// Geothermal gradient beneath madagascar using the ellipsoid model

#include <aspect/initial_conditions/interface.h>
#include <aspect/initial_conditions/geothermal_gradient.h>
#include <fstream>
#include <iostream>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * Function that return latitude and longitude from cartesian
     * coordinates that account for ellipsoidal shape of the Earth
     */

    template <int dim>
    std::pair<double, double>
    GeothermalModel<dim>::
	lat_long_from_xyz_wgs84(const Point<3> &pos) const
    {
    	    /* WGS84 ellipsoid constants */
    	    const double radius = 6378137;
    	    const double ellipticity = 8.1819190842622e-2;

    	    /* calculations */
    	    const double b = std::sqrt(radius * radius
    	                               * (1 - ellipticity * ellipticity));
    	    const double ep = std::sqrt((radius * radius - b * b) / (b * b));
    	    const double p = std::sqrt(pos(0) * pos(0) + pos(1) * pos(1));
    	    const double th = std::atan2(radius * pos(2), b * p);
    	    const double lon = std::atan2(pos(1), pos(0));
    	    const double lat = std::atan2((pos(2) + ep * ep * b * std::sin(th)
    	                                   * std::sin(th) * std::sin(th)),
    	                                  (p - (ellipticity
    	                                        * ellipticity
    	                                        * radius
    	                                        * (std::cos(th) * std::cos(th)
    	                                        * std::cos(th)))));

    	    /* convert to degrees */
    	    const double lon_degrees = lon * (180 / numbers::PI);
    	    const double lat_degrees = lat * (180 / numbers::PI);

    	    /* set all longitudes between [0,360] */
    	    if (lon_degrees < 0)
    	      return std::make_pair(lat_degrees, lon_degrees + 360);
    	    else if (lon_degrees > 360)
    	      return std::make_pair(lat_degrees, lon_degrees - 360);
    	    else
    	      return std::make_pair(lat_degrees, lon_degrees);

    }


    template <>
    double
    GeothermalModel<2>::get_lithosphere_isotherm(const double,
                                              const double) const
    {       abort();
        	return 0;
    }

    template <int dim>
    double
    GeothermalModel<dim>::get_lithosphere_isotherm(const double latitude,
                                                  const double longitude) const
    {
          // loop over the entire array and see if we find a point
          // that's within delta of what we're looking for. the data
          // is arranged in a way that keeps the latitude constant
          // while running over all longitudes, and when we're done
          // with that changes the latitude by one step. so if the
          // latitude is wrong, we can definitely skip ahead a whole
          // set of longitudes. The number of values to skip is calculated.
          for (unsigned int i = 0; i <= latitudes_iso.size();)
            if (std::fabs(latitude - latitudes_iso[i]) <= delta)
              {
                if (std::fabs(longitude - longitudes_iso[i]) <= delta)
                  return -depths_iso[i]*1000;
                else
                  ++i;
              }
            else
              i += number_coords_litho;

          Assert(false, ExcInternalError());
          return 0;
     }//End definition of get_lithosphere


    template <int dim>
    double
    GeothermalModel<dim>::initial_temperature (const Point<dim> &) const
    {
          Assert (false, ExcNotImplemented());
          return 0;
    }

    template <>
    double
    GeothermalModel<3>::initial_temperature (const Point<3> &position) const
    {
    // get the depth of the lithosphere isotherm for the current lat/long
          // position
          const std::pair<double, double> lat_long = lat_long_from_xyz_wgs84(position);

          //TODO: this is the depth with respect to the sphere; need WGS84 here
          const double radius = 6378137;
          const double depth = position.norm() - radius;

          // if above the isotherm, use a linear behavior up to the surface
          // below the isotherm, use an increase of .5 degrees per kilometer
          // (0.0005 degrees per meter)
          const double isotherm_temp = 1673.15;

          const double isotherm_depth = this->get_lithosphere_isotherm(lat_long.first, lat_long.second);
          if (depth < isotherm_depth)
            return isotherm_temp - (depth-isotherm_depth) * 0.0005;
          else
            return 273.15+depth/isotherm_depth*1400;
     }//End the function initialize_temperature


    template <int dim>
    void
    GeothermalModel<dim>::declare_parameters(ParameterHandler &prm)
    {
         prm.enter_subsection("Initial conditions");
         {
           prm.enter_subsection("Geothermal model");
           {
             prm.declare_entry("Lithosphere Isotherm filename",
                               "thickness.txt",
                               Patterns::FileName(),
                               "Surface coordinates and depths to the 1673.15 K isotherm. Units: degrees and kilometers.");
           }
           prm.leave_subsection();
         }
         prm.leave_subsection();
    }//End the void function declare_parameters

    template <int dim>
    void
    GeothermalModel<dim>::parse_parameters(ParameterHandler &prm)
    {
        prm.enter_subsection("Initial conditions");
        {
          prm.enter_subsection("Geothermal model");
          {
            litho_isotherm_file_name = prm.get("Lithosphere Isotherm filename");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        std::ifstream input1(litho_isotherm_file_name.c_str());
        AssertThrow (input1.is_open(),
                     ExcMessage (std::string("Can't read from file <") + litho_isotherm_file_name + ">"));


        // Start loop with int count to calculate delta below for lithospheric thickness
        int count = 0;
        while (true)
          { 
            double latitude_iso, longitude_iso, depth_iso;
            input1 >> latitude_iso >>longitude_iso >> depth_iso;
            if (input1.eof())
            break;
            latitudes_iso.push_back(latitude_iso);
            longitudes_iso.push_back(longitude_iso);
            depths_iso.push_back(depth_iso);

            /** Find first 2 numbers that are different to use in
             * calculating half the difference between each position as delta.
             */
            if (count < 2 )
              {
                count ++;
              }
            if (count == 2)
              {
                if ((std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9) && (std::fabs(longitudes_iso[0] - longitudes_iso[1]) > 1e-9))
                  {
                    // Stop program if file formatted incorrectly.
                    std::cout << ""<< std::endl;
                    throw std::ios_base::failure("Lithospheric thickness file not formatted correctly. " + litho_isotherm_file_name + "Make sure you have lat, lon, value with lat. or lon. varying.");
                  }
                if ((std::fabs(latitudes_iso[0] - latitudes_iso[1]) < 1e-9) && (std::fabs(longitudes_iso[0] - longitudes_iso[1]) < 1e-9))
                  {
                    // Stop program if file formatted incorrectly.
                    std::cout << ""<< std::endl;
                    throw std::ios_base::failure("Lithospheric thickness file not formatted correctly. " + litho_isotherm_file_name + "Make sure you have lat, lon, value with lat. or lon. varying.");
                  }

                if (std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9)
                  {
                    // Calculate delta as half the distance between points.
                    delta = std::fabs((0.5)*(latitudes_iso[0] - latitudes_iso[1]));
                    // If flag is 0 then longitudes grouped and we calculate delta from latitudes
                    litho_flag = 0;
                  }
                else
                  {
                    // Calculate delta as half the distance between points.
                    delta = std::fabs((0.5)*(longitudes_iso[0] - longitudes_iso[1]));
                    // If flag is 1 then latitudes are grouped and we calculate delta from longitudes
                    litho_flag = 1;
                  }
                std::cout << ""<< std::endl;
                std::cout<<"Lithosphere thickness delta = "<< delta << std::endl;
                std::cout << ""<< std::endl;
                std::cout<<"Resolution of input lithosphere thickness in meters is approximately = "<< delta*111*2 << std::endl;
                std::cout << ""<< std::endl;
                count++;
              }
          } // End of loop for calculating delta for the lithosphere thickness file

        //Calculate the number of unique longitudes or latitudes from the lithosphere isotherm file and crustal thickness file.
        double c,d;
        int count3;
        count3 = 0;
        c = 0;
        d = 0;

        if ( litho_flag == 1 )
          {
            c = latitudes_iso[0];
            d = latitudes_iso[1];
            count3 = 2;

            while (c-d < 1e-9)
              {
                c = d;
                d = latitudes_iso[count3];
                count3++;
              }
          }

        if ( litho_flag == 0 )
          {
            c = longitudes_iso[0];
            d = longitudes_iso[1];
            count3 = 2;

            while (c-d < 1e-9)
              {
                c = d;
                d = longitudes_iso[count3];
                count3++;
              }
          }

        std::cout << ""<< std::endl;
        std::cout<<"number of unique latitudes or longitudes in lithosphere thickness file= "<< count3 - 1 << std::endl;
        number_coords_litho = count3-1;
        std::cout << ""<< std::endl;
      }
    
  



  }
}


namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(GeothermalModel,
                                       "geothermal model",
                                       "In subsection lithosphere isotherm we define "
                                       "the extent of the conductive heat "
                                       "equation is assumed to be 1400 C (1673.15 K) "
                                       "as previously used by Bird et al., 2008 "
                                       "and Stamps et al. (in prep). This assumption "
                                       "is consistent with the Schubert 2001 "
                                       "definition of the mechanical lithosphere.")
  }
}


