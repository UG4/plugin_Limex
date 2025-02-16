/*
 * Copyright (c) 2025: Goethe University Frankfurt
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

#ifndef __H__LIMEX__SIMPLE_INTEGRATOR_HPP__
#define __H__LIMEX__SIMPLE_INTEGRATOR_HPP__

#include "lib_algebra/operator/debug_writer.h"

namespace ug {

/// Integrate (a non-linear problem) over a given time interval
template<class TDomain, class TAlgebra>
class SimpleTimeIntegrator :
		public INonlinearTimeIntegrator<TDomain, TAlgebra>,
		public ITimeDiscDependentObject<TAlgebra>,
		public DebugWritingObject<TAlgebra>
{
protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;

public:
	typedef INonlinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef ITimeDiscretization<TAlgebra> time_disc_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef IGridFunctionSpace<grid_function_type> grid_function_space_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

	// constructor
	SimpleTimeIntegrator (SmartPtr<time_disc_type> tDisc)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc),
	  m_spBanachSpace(new AlgebraicSpace<grid_function_type>() ),
	  m_spDerivative(SPNULL), m_initial_consistency_error(0.0)

	{}

	SimpleTimeIntegrator
	(
		SmartPtr<time_disc_type> tDisc,
		SmartPtr<grid_function_space_type> spSpace
	)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc),
	  m_spBanachSpace(spSpace),
	  m_spDerivative(SPNULL), m_initial_consistency_error(0.0)
	{}


	bool apply
	(
		SmartPtr<grid_function_type> u1,
		number t1,
		ConstSmartPtr<grid_function_type> u0,
		number t0
	)
	{
		time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
		
		if (tdisc.num_stages() == 1)
			return apply_single_stage(u1,t1,u0,t0);
		else
			return apply_multi_stage(u1,t1,u0,t0);
	}

	void set_derivative(SmartPtr<grid_function_type> udot)
	{ m_spDerivative = udot; }

	SmartPtr<grid_function_type> get_derivative()
	{ return m_spDerivative; }

	number get_consistency_error() const
	{ return m_initial_consistency_error; }

	void set_banach_space(SmartPtr<IGridFunctionSpace<grid_function_type> > spSpace)
	{ m_spBanachSpace = spSpace; }

protected:

	bool apply_single_stage
	(
		SmartPtr<grid_function_type> u1,
		number t1,
		ConstSmartPtr<grid_function_type> u0,
		number t0
	);
	
	bool apply_multi_stage
	(
		SmartPtr<grid_function_type> u1,
		number t1,
		ConstSmartPtr<grid_function_type> u0,
		number t0
	);

	inline bool hasTerminated(double tCurrent, double tStart, double tFinal) const
	{
	 	/*return (! ((tCurrent < tFinal) && (tFinal-tCurrent > base_type::m_precisionBound)));*/
		return tCurrent >= tFinal || tFinal-tCurrent < (tFinal-tStart)*base_type::m_precisionBound;
	}

	/// metric
	SmartPtr<IGridFunctionSpace<grid_function_type> > m_spBanachSpace;

	SmartPtr<grid_function_type> m_spDerivative;

	number m_initial_consistency_error;
};

} // ug

#include "simple_integrator_impl.hpp"

#endif // __H__LIMEX__SIMPLE_INTEGRATOR_HPP__
