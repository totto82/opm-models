/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser
 *   Institute of Hydraulic Engineering                                      *
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
 * \ingroup Components
 * \brief A component that only throws exceptions.
 *
 * Its main purpose is to make things compile and give runtime errors
 * if it is actually used.
 */
#ifndef DUMUX_NULL_COMPONENT_HH
#define DUMUX_NULL_COMPONENT_HH

#include "component.hh"

namespace Dumux
{
/*!
 * \ingroup Components
 *
 * \brief A component that only throws exceptions.
 *
 * Its main purpose is to make things compile and give runtime errors
 * if it is actually used.
 */
template <class Scalar>
class NullComponent : public Component<Scalar, NullComponent<Scalar> >
{
};

} // end namepace

#endif
