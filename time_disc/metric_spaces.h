/*
 * Copyright (c) 2014-2017:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
 *
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef METRIC_SPACES_H_
#define METRIC_SPACES_H_

// c++ includes
#include <vector>
#include <cmath>

// ug includes
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/operator/debug_writer.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#include "lib_disc/function_spaces/integrate.h"


/// Abstract base class for (algebraic) vectors
template<typename TVector>
class IBanachSpace
{

public:
	virtual ~IBanachSpace() {}

	/// euclidean norm (default)
	virtual double norm(TVector &x) const
	{ return x.norm(); }

	virtual double distance(TVector& x, TVector& y) const
	{
		SmartPtr<TVector> delta = x.clone();
		*delta -= y;
		return norm(*delta);
	}
};



/// Abstract base class for grid functions
template<typename TGridFunction>
class IGridFunctionSpace  : public IBanachSpace<typename TGridFunction::vector_type>
{

public:
	typedef typename TGridFunction::vector_type vector_type;
	typedef TGridFunction grid_function_type;

	/// DTOR
	virtual ~IGridFunctionSpace() {}

	/// norm (for grid functions)
	virtual double norm(TGridFunction& x) const
	{ return x.norm(); }

	/// distance (for grid functions)
	virtual double distance(TGridFunction& x, TGridFunction& y) const
	{
		SmartPtr<TGridFunction> delta = x.clone();
		*delta -= y;
		return norm(*delta);
	}

	/// norm (for vectors)
	virtual double norm(vector_type &x) const
	{
		TGridFunction* gfX=dynamic_cast< TGridFunction*>(&x);
		UG_ASSERT(gfX!=NULL, "Huhh: GridFunction required!");
		return norm(*gfX);
	}

	/// distance (for vectors)
	virtual double distance(vector_type &x, vector_type &y) const
	{ return distance(static_cast<TGridFunction &>(x), static_cast<TGridFunction &>(y)); }




};


#endif
