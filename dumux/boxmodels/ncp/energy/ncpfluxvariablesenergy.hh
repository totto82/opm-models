// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008,2009 by Andreas Lauser                               *
 *   Copyright (C) 2008,2009 by Melanie Darcis                               *
 *                                                                           *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Contains the quantities to calculate the energy flux in the
 *        Ncp box model.
 */
#ifndef DUMUX_NCP_ENERGY_FLUX_VARIABLES_HH
#define DUMUX_NCP_ENERGY_FLUX_VARIABLES_HH

#include <dumux/boxmodels/ncp/ncpproperties.hh>

#include <dumux/common/spline.hh>

#include <dune/common/fvector.hh>

namespace Dumux
{

template <class TypeTag, bool enableEnergy/*=false*/>
class NcpFluxVariablesEnergy
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

public:
    NcpFluxVariablesEnergy()
    {}

    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {}
};

template <class TypeTag>
class NcpFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GridView::ctype CoordScalar;

    enum {
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    NcpFluxVariablesEnergy()
    {}

    void update(const ElementContext &elemCtx, int scvfIdx, int timeIdx)
    {
        const auto &fvElemGeom = elemCtx.fvElemGeom(timeIdx);
        const auto &scvf = fvElemGeom.subContVolFace[scvfIdx];

        // calculate temperature gradient using finite element
        // gradients
        DimVector tmp(0.0);
        DimVector temperatureGrad(0.);
        for (int scvIdx = 0; scvIdx < elemCtx.numScv(); scvIdx++)
        {
            tmp = fvElemGeom.subContVolFace[scvfIdx].grad[scvIdx];
            tmp *= elemCtx.volVars(scvIdx, timeIdx).fluidState().temperature(/*phaseIdx=*/0);
            temperatureGrad += tmp;
        }

        // scalar product of temperature gradient and scvf normal
        temperatureGradNormal_ = 0.0;
        for (int i = 0; i < dimWorld; ++ i)
            temperatureGradNormal_ += scvf.normal[i]*temperatureGrad[i];

        const auto &volVarsInside = elemCtx.volVars(scvf.i, timeIdx);
        const auto &volVarsOutside = elemCtx.volVars(scvf.j, timeIdx);

        // arithmetic mean
        heatConductivity_ = 
            0.5 * (volVarsInside.heatConductivity()
                   +
                   volVarsOutside.heatConductivity());
        Valgrind::CheckDefined(heatConductivity_);
    }

    /*!
     * \brief The lumped / average conductivity of solid plus phases \f$[W/mK]\f$.
     */
    Scalar heatConductivity() const
    { return heatConductivity_; }

    /*!
     * \brief The normal of the gradient of temperature .
     */
    Scalar temperatureGradNormal() const
    { return temperatureGradNormal_; }

private:
    Scalar heatConductivity_;
    Scalar temperatureGradNormal_;
};

} // end namepace

#endif
