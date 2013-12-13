// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2010-2013 by Andreas Lauser
  Copyright (C) 2011 by Markus Wolff

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \ingroup NcpModel
 *
 * \brief Defines default values for the properties required for the
 *        NCP compositional multi-phase model.
 */
#ifndef EWOMS_NCP_PROPERTY_DEFAULTS_HH
#define EWOMS_NCP_PROPERTY_DEFAULTS_HH

#include "ncpmodel.hh"
#include "ncplocalresidual.hh"
#include "ncpfluxvariables.hh"
#include "ncpprimaryvariables.hh"
#include "ncpboundaryratevector.hh"
#include "ncpratevector.hh"
#include "ncpvolumevariables.hh"
#include "ncpnewtonmethod.hh"
#include "ncpindices.hh"
#include "ncpproperties.hh"

#include <ewoms/models/common/multiphasebaseproblem.hh>
#include <ewoms/models/common/velocity.hh>

#include <opm/material/constraintsolvers/CompositionFromFugacities.hpp>
#include <opm/material/heatconduction/DummyHeatConductionLaw.hpp>

namespace Opm {
namespace Properties {
/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system.
 */
SET_INT_PROP(NcpModel, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_INT_PROP(NcpModel, NumPhases, GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases);

/*!
 * \brief Set the property for the number of equations and primary variables.
 */
SET_INT_PROP(NcpModel, NumEq,GET_PROP_TYPE(TypeTag, Indices)::numEq);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(NcpModel, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

/*!
 * \brief Extract the type parameter objects for the heat conduction law
 *        from the law itself
 */
SET_TYPE_PROP(NcpModel,
              HeatConductionLaw,
              Opm::DummyHeatConductionLaw<typename GET_PROP_TYPE(TypeTag, Scalar)>);

/*!
 * \brief Extract the type parameter objects for the heat conduction law
 *        from the law itself.
 */
SET_TYPE_PROP(NcpModel,
              HeatConductionLawParams,
              typename GET_PROP_TYPE(TypeTag, HeatConductionLaw)::Params);

/*!
 * \brief Set the themodynamic constraint solver which calculates the
 *        composition of any phase given all component fugacities.
 */
SET_PROP(NcpModel, NcpCompositionFromFugacitiesSolver)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    typedef Opm::CompositionFromFugacities<Scalar, FluidSystem> type;
};

//! Use the Ncp local jacobian operator for the compositional NCP model
SET_TYPE_PROP(NcpModel,
              LocalResidual,
              Ewoms::NcpLocalResidual<TypeTag>);

//! Use the Darcy relation by default
SET_TYPE_PROP(NcpModel, VelocityModule, Ewoms::DarcyVelocityModule<TypeTag>);

//! Use the Ncp specific newton method for the compositional NCP model
SET_TYPE_PROP(NcpModel, NewtonMethod, Ewoms::NcpNewtonMethod<TypeTag>);

//! the Model property
SET_TYPE_PROP(NcpModel, Model, Ewoms::NcpModel<TypeTag>);

//! The type of the base base class for actual problems
SET_TYPE_PROP(NcpModel, BaseProblem, Ewoms::MultiPhaseBaseProblem<TypeTag>);

//! use an isothermal model by default
SET_BOOL_PROP(NcpModel, EnableEnergy, false);

//! disable diffusion by default
SET_BOOL_PROP(NcpModel, EnableDiffusion, false);

//! do not use smooth upwinding by default
SET_BOOL_PROP(NcpModel, EnableSmoothUpwinding, false);

//! the RateVector property
SET_TYPE_PROP(NcpModel, RateVector, Ewoms::NcpRateVector<TypeTag>);

//! the BoundaryRateVector property
SET_TYPE_PROP(NcpModel, BoundaryRateVector, Ewoms::NcpBoundaryRateVector<TypeTag>);

//! the PrimaryVariables property
SET_TYPE_PROP(NcpModel, PrimaryVariables, Ewoms::NcpPrimaryVariables<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(NcpModel, VolumeVariables, Ewoms::NcpVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(NcpModel, FluxVariables, Ewoms::NcpFluxVariables<TypeTag>);

//! truncate the newton update for the first 3 iterations of a time step
SET_INT_PROP(NcpModel, NcpNewtonNumChoppedIterations, 3);

//! The indices required by the compositional NCP model
SET_TYPE_PROP(NcpModel, Indices, Ewoms::NcpIndices<TypeTag, 0>);

//! The unmodified weight for the pressure primary variable
SET_SCALAR_PROP(NcpModel, NcpPressureBaseWeight, 1.0);
//! The weight for the saturation primary variables
SET_SCALAR_PROP(NcpModel, NcpSaturationsBaseWeight, 1.0);
//! The unmodified weight for the fugacity primary variables
SET_SCALAR_PROP(NcpModel, NcpFugacitiesBaseWeight, 1.0);

}} // namespace Properties, Opm

#endif
