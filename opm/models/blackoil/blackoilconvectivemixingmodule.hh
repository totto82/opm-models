// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Classes required for molecular diffusion.
 */
#ifndef EWOMS_CONVECTIVEMIXING_MODULE_HH
#define EWOMS_CONVECTIVEMIXING_MODULE_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <stdexcept>
#include <iostream>

namespace Opm {

/*!
 * \copydoc Opm::BlackOilConvectiveMixingModule
 * \brief Provides the requiered methods for dynamic convective mixing.
 */
template <class TypeTag>
class BlackOilConvectiveMixingModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { dimWorld = GridView::dimensionworld };

public:

    #if HAVE_ECL_INPUT
    /*!
     * \brief Initialize all internal data structures needed by the convective mixing module
     */
    static void initFromState(const EclipseState& eclState, const Schedule& schedule)
    {
        std::size_t numRegions = eclState.runspec().tabdims().getNumPVTTables();
        int episodeIdx = 0; // we dont allow for DRSDTCON to change during iterations. TODO: Add check
        const auto& control = schedule[episodeIdx].oilvap();
        active_ = control.drsdtConvective();
        Xhi_.resize(numRegions); 
        Smo_.resize(numRegions); 
        for (size_t i = 0; i < numRegions; ++i ) {
            Xhi_[i] = control.getMaxDRSDT(i);
            Smo_[i] = control.getSmo(i);
        }
    }

    static void beginEpisode(const EclipseState& eclState, const Schedule& schedule, const int episodeIdx)
    {
        // check that Xhi and Smo didn't change
        std::size_t numRegions = eclState.runspec().tabdims().getNumPVTTables();
        const auto& control = schedule[episodeIdx].oilvap();
        active_ = control.drsdtConvective();
        Xhi_.resize(numRegions); 
        Smo_.resize(numRegions); 
        bool changed = false;
        for (size_t i = 0; i < numRegions; ++i ) {
            if (Xhi_[i] != control.getMaxDRSDT(i) || Smo_[i] != control.getSmo(i)) {
                changed = true;
            }
        }
        if (changed) {
            throw("SMO changed");
        }
    }
    #endif

    template <class Context>
    static void addConvectiveMixingFlux(RateVector& flux, 
                                        const Context& elemCtx,
                                        unsigned scvfIdx, 
                                        unsigned timeIdx)
    {

        // need for dary flux calculation
        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        unsigned interiorDofIdx = scvf.interiorIndex();
        unsigned exteriorDofIdx = scvf.exteriorIndex();
        assert(interiorDofIdx != exteriorDofIdx);

        const auto& globalIndexIn = stencil.globalSpaceIndex(interiorDofIdx);
        const auto& globalIndexEx = stencil.globalSpaceIndex(exteriorDofIdx);
        Scalar trans = problem.transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar faceArea = scvf.area();
        const Scalar g = problem.gravity()[dimWorld - 1];
        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);
        const Scalar zIn = problem.dofCenterDepth(elemCtx, interiorDofIdx, timeIdx);
        const Scalar zEx = problem.dofCenterDepth(elemCtx, exteriorDofIdx, timeIdx);
        const Scalar distZ = zIn - zEx;
        addConvectiveMixingFlux(flux,
                                intQuantsIn,
                                intQuantsEx,
                                globalIndexIn,
                                globalIndexEx,
                                distZ * g,
                                trans,
                                faceArea);
    }





    /*!
     * \brief Adds the diffusive mass flux flux to the flux vector over a flux
     *        integration point.
      */
    static void addConvectiveMixingFlux(RateVector& flux,
                            const IntensiveQuantities& intQuantsIn,
                            const IntensiveQuantities& intQuantsEx,
                            const unsigned globalIndexIn,
                            const unsigned globalIndexEx,
                            const Scalar distZg, 
                            const Scalar trans,
                            const Scalar faceArea)
    {


        if (!active_) {
            return;
        }

        const Scalar rs_zero_in = 0.0;
        const auto& t_in = Opm::getValue(intQuantsIn.fluidState().temperature(FluidSystem::oilPhaseIdx));
        const auto& p_in = Opm::getValue(intQuantsIn.fluidState().pressure(FluidSystem::oilPhaseIdx));
        const auto& rssat_in = Opm::getValue(intQuantsIn.fluidState().RsSat()); 
        const auto rho_in
            = FluidSystem::oilPvt().inverseFormationVolumeFactor(intQuantsIn.pvtRegionIndex(), t_in, p_in, rs_zero_in)
            * FluidSystem::oilPvt().oilReferenceDensity(intQuantsIn.pvtRegionIndex());
        const auto rho_sat_in = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(intQuantsIn.pvtRegionIndex(), t_in, p_in)
            * (FluidSystem::oilPvt().oilReferenceDensity(intQuantsIn.pvtRegionIndex())
            + rssat_in * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, intQuantsIn.pvtRegionIndex()));
        
        //exteriour

        const Scalar rs_zero_ex = 0.0;
        const auto& t_ex = Opm::getValue(intQuantsEx.fluidState().temperature(FluidSystem::oilPhaseIdx));
        const auto& p_ex = Opm::getValue(intQuantsEx.fluidState().pressure(FluidSystem::oilPhaseIdx));
        const auto& rssat_ex = Opm::getValue(intQuantsEx.fluidState().RsSat());
        const auto rho_ex
            = FluidSystem::oilPvt().inverseFormationVolumeFactor(intQuantsEx.pvtRegionIndex(), t_ex, p_ex, rs_zero_ex)
            * FluidSystem::oilPvt().oilReferenceDensity(intQuantsEx.pvtRegionIndex());
        const auto rho_sat_ex = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(intQuantsEx.pvtRegionIndex(), t_ex, p_ex)
            * (FluidSystem::oilPvt().oilReferenceDensity(intQuantsEx.pvtRegionIndex())
            + rssat_ex * FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, intQuantsEx.pvtRegionIndex()));
        
        //rho difference approximation
        const auto delta_rho = (rho_sat_ex + rho_sat_in - rho_in -rho_ex)/2;
        const auto pressure_difference_convective_mixing =  delta_rho * distZg;

        //if change in pressure
        if (Opm::abs(pressure_difference_convective_mixing) > 1e-12){

            // find new upstream direction
            short interiorDofIdx = 0; // NB
            short exteriorDofIdx = 1; // NB
            unsigned upIdx = 0;//scvf.interiorIndex();
            unsigned downIdx = 1;//scvf.exteriorIndex();

            if (pressure_difference_convective_mixing > 0) {
                upIdx = exteriorDofIdx;//scvf.exteriorIndex();
                downIdx = interiorDofIdx;//scvf.interiorIndex();
            }

            const auto& up = (upIdx == interiorDofIdx) ? intQuantsIn : intQuantsEx;
            unsigned globalUpIndex = (upIdx == interiorDofIdx) ? globalIndexIn : globalIndexEx;

            const auto& down = (downIdx == interiorDofIdx) ? intQuantsIn : intQuantsEx;
            unsigned globalDownIndex = (downIdx == interiorDofIdx) ? globalIndexIn : globalIndexEx;

            
            const auto& Rs =  up.fluidState().Rs();
            //const Evaluation SoMax = 0.0;
            const auto& RsSat = up.fluidState().RsSat();
            const Evaluation& transMult = up.rockCompTransMultiplier();
            
            Evaluation sg = up.fluidState().saturation(FluidSystem::gasPhaseIdx);
            Evaluation S = (Rs - RsSat * sg) / (RsSat * ( 1.0 - sg));
            if ( (S > Smo_[up.pvtRegionIndex()] || down.fluidState().Rs() > 0) ) {
                const auto& invB = up.fluidState().invB(FluidSystem::oilPhaseIdx);
                const auto& visc = up.fluidState().viscosity(FluidSystem::oilPhaseIdx);
                // what will be the flux when muliplied with trans_mob
                const auto convectiveFlux = -trans*transMult*Xhi_[up.pvtRegionIndex()]*invB*pressure_difference_convective_mixing*Rs/(visc*faceArea);
                unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);                   
                if (globalUpIndex == globalIndexIn)
                    flux[conti0EqIdx + activeGasCompIdx] += convectiveFlux;
                else
                    flux[conti0EqIdx + activeGasCompIdx] += Opm::getValue(convectiveFlux);
            }
        } 
    };

    private:
    static bool active_;
    static std::vector<double> Xhi_;
    static std::vector<double> Smo_;
};

template <typename TypeTag>
bool
BlackOilConvectiveMixingModule<TypeTag>::active_;

template <typename TypeTag>
std::vector<double>
BlackOilConvectiveMixingModule<TypeTag>::Xhi_;

template <typename TypeTag>
std::vector<double>
BlackOilConvectiveMixingModule<TypeTag>::Smo_;

}
#endif



