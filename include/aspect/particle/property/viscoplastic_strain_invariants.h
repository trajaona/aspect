/*
 Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_property_viscoplastic_strain_invariant_h
#define _aspect_particle_property_viscoplastic_strain_invariant_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/visco_plastic.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that calculates the finite strain invariant a particle has
       * experienced and assigns it to either the plastic and/or viscous strain field based
       * on whether the material is plastically yielding, or the total strain field
       * used in the visco plastic material model. The implementation of this property
       * is equivalent to the implementation for compositional fields that is located in
       * the plugin <code>benchmarks/buiter_et_al_2008_jgr/plugin/finite_strain_invariant.cc</code>,
       * and is effectively the same as what the visco plastic material model uses for compositional fields.
       * @ingroup ParticleProperties
       */
      template <int dim>
      class ViscoPlasticStrainInvariant : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
          * Constructor.
          */
          ViscoPlasticStrainInvariant ();

          virtual void initialize ();


          /**
          * @copydoc aspect::Particle::Property::Interface::initialize_one_particle_property()
          **/
          virtual
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const;

          /**
          * @copydoc aspect::Particle::Property::Interface::update_one_particle_property()
          **/
          virtual
          void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim> &position,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim> > &gradients,
                                        const ArrayView<double> &particle_properties) const;

          /**
          * @copydoc aspect::Particle::Property::Interface::need_update()
          **/
          UpdateTimeFlags
          need_update () const;

          /**
          * @copydoc aspect::Particle::Property::Interface::get_needed_update_flags()
          **/
          virtual
          UpdateFlags
          get_needed_update_flags () const;

          /**
          * @copydoc aspect::Particle::Property::Interface::get_property_information()
          **/
          virtual
          std::vector<std::pair<std::string, unsigned int> >
          get_property_information() const;

        private:
          unsigned int n_components;

      };
    }
  }
}

#endif
