/*
 * Copyright (c) 2014-2016:  G-CSC, Goethe University Frankfurt
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

#ifndef LIMEX_H_
#define LIMEX_H_

// extern libraries
#include <vector>
#include <cmath>

// ug libraries
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"

#include "lib_algebra/algebra_template_define_helper.h"
#include "lib_algebra/operator/debug_writer.h"

namespace ug{


/// linear implicit time stepping scheme
/**
 * This time stepping scheme discretizes equations of the form
 * \f[
 * 	M \partial_t u(t) = f(t)
 * \f]
 * as
 * \f[
 * 	(M - \Delta t J) \left( u(t^{k+1}) - u(t^k) \right)  =  \Delta t \cdot f(t^{k})
 * \f]
 *
 * Thus, for \f$\theta = 1 \f$ this is the Backward-Euler time stepping.
 */
template <class TAlgebra>
class LinearImplicitEuler
: public ITimeDiscretization<TAlgebra>,
  public DebugWritingObject<TAlgebra>

{
public:
/// Type of base class
	typedef ITimeDiscretization<TAlgebra> base_type;

/// Type of algebra
	typedef TAlgebra algebra_type;

/// Type of algebra matrix
	typedef typename algebra_type::matrix_type matrix_type;

/// Type of algebra vector
	typedef typename algebra_type::vector_type vector_type;

/// Domain Discretization type
	typedef IDomainDiscretization<algebra_type>	domain_discretization_type;

public:
	/// CTOR
	LinearImplicitEuler(SmartPtr<IDomainDiscretization<algebra_type> > spDD)
		: ITimeDiscretization<TAlgebra>(spDD),
		  m_pPrevSol(NULL),
		  m_dt(0.0),
		  m_futureTime(0.0),
		  m_spMatrixJDisc(spDD), m_spMatrixJOp(SPNULL), m_bMatrixJNeedsUpdate(true), m_useLinearMode(false),
		  m_spGammaDisc(SPNULL), m_spGammaOp(SPNULL), m_bGammaNeedsUpdate(true),
		  m_useCachedMatrices(true)
	{}

	/// CTOR
	LinearImplicitEuler(SmartPtr<IDomainDiscretization<algebra_type> > spDefectDisc,
						SmartPtr<IDomainDiscretization<algebra_type> > spMatrixJDisc)
		: ITimeDiscretization<TAlgebra>(spDefectDisc),
		  m_pPrevSol(NULL),
		  m_dt(0.0),
		  m_futureTime(0.0),
		  m_spMatrixJDisc(spMatrixJDisc), m_spMatrixJOp(SPNULL), m_bMatrixJNeedsUpdate(true), m_useLinearMode(false), 
		  m_spGammaDisc(SPNULL), m_spGammaOp(SPNULL), m_bGammaNeedsUpdate(true),
		  m_useCachedMatrices(true)
	{}

	/// CTOR
	LinearImplicitEuler(SmartPtr<IDomainDiscretization<algebra_type> > spDefectDisc,
						SmartPtr<IDomainDiscretization<algebra_type> > spMatrixJDisc,
						SmartPtr<IDomainDiscretization<algebra_type> > spGammaDisc)
		: ITimeDiscretization<TAlgebra>(spDefectDisc),
		  m_pPrevSol(NULL),
		  m_dt(0.0),
		  m_futureTime(0.0),
		  m_spMatrixJDisc(spMatrixJDisc), m_spMatrixJOp(SPNULL), m_bMatrixJNeedsUpdate(true), m_useLinearMode(false), 
		  m_spGammaDisc(spGammaDisc), m_spGammaOp(SPNULL), m_bGammaNeedsUpdate(true),
		  m_useCachedMatrices(true)
	{}

	/// DTOR
	virtual ~LinearImplicitEuler(){};

public:
	/// \copydoc ITimeDiscretization::num_prev_steps()
	virtual size_t num_prev_steps() const {return m_prevSteps;}

	///	\copydoc ITimeDiscretization::prepare_step()
	virtual void prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
	                          number dt);

	///	\copydoc ITimeDiscretization::prepare_step_elem()
	virtual void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
	                               number dt, const GridLevel& gl);

	///	\copydoc ITimeDiscretization::finish_step()
	virtual void finish_step(SmartPtr<VectorTimeSeries<vector_type> > currSol) {};

	///	\copydoc ITimeDiscretization::finish_step_elem()
	virtual void finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
	                              const GridLevel& gl);

	virtual number future_time() const {return m_futureTime;}

public:

	/// Meant to assemble J(u) c = d(u), but abused here...  (u not used!)
	void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl);

	/// Meant to assemble d(u), but abused here... (u not used!)
	void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl);

	/// Should return (M+tau A) delta = tau f
	void assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl);

	void assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl);

	void assemble_rhs(vector_type& b, const GridLevel& gl);

	void adjust_solution(vector_type& u, const GridLevel& gl);

	 virtual size_t num_stages() const {return 1;};
	 virtual void set_stage(size_t stage) {};


	/// Some simplifications for linear systems. (In this case, the mass matrix is not re-assembled.)
	void enable_linear_mode() { m_useLinearMode = true; }
	void disable_linear_mode() { m_useLinearMode = false; }
	bool use_linear_mode() const { return m_useLinearMode; }

	void enable_matrix_cache() { m_useCachedMatrices = true; }
	void disable_matrix_cache() { m_useCachedMatrices = false; }
	void set_matrix_cache(bool useCache) { m_useCachedMatrices = useCache; }

protected:

	virtual number update_scaling(std::vector<number>& vSM,
			                              std::vector<number>& vSA,
			                              number dt, number currentTime,
			                              ConstSmartPtr<VectorTimeSeries<vector_type> > prevSol)

	{
		//	resize scaling factors
		vSM.resize(1);
		vSM[0] = 1.0;

		vSA.resize(1);
		vSA[0] = dt;

		return currentTime + dt;
	}

	static const size_t m_prevSteps=1;			///< number of previous steps needed.
	std::vector<number> m_vScaleMass;			///< Scaling for mass part
	std::vector<number> m_vScaleStiff;			///< Scaling for stiffness part

	SmartPtr<VectorTimeSeries<vector_type> > m_pPrevSol;		///< Previous solutions
	number m_dt; 								///< Time Step size
	number m_futureTime;						///< Future Time

	// discretization for defect
	using base_type::m_spDomDisc;


	// constant matrix $$ \tau J$$
	SmartPtr<IDomainDiscretization<algebra_type> > m_spMatrixJDisc;
	SmartPtr<AssembledLinearOperator<algebra_type> > m_spMatrixJOp;		///< Operator
	bool m_bMatrixJNeedsUpdate;
	bool m_useLinearMode;

	// Matrix: $\Gamma[u0, u0']$
	SmartPtr<IDomainDiscretization<algebra_type> > m_spGammaDisc;		///< Gamma disc
	SmartPtr<AssembledLinearOperator<algebra_type> > m_spGammaOp;	    ///< Gamma operator
	bool m_bGammaNeedsUpdate;

	SmartPtr<matrix_type> m_spMatrixCacheMk;

	bool m_useCachedMatrices;

public:

	/// Invalidate all cached operators
	void invalidate()
	{
		m_spMatrixCacheMk = SPNULL;
		m_bMatrixJNeedsUpdate = true;
		invalidate_gamma();
	}

	/// Invalidate Gamma operator
	void invalidate_gamma()
	{ m_spGammaOp = SPNULL; m_bGammaNeedsUpdate = true; }

	void set_gamma_disc(SmartPtr<IDomainDiscretization<algebra_type> > spGammaDisc)
	{m_spGammaDisc = spGammaDisc;}

	using DebugWritingObject<TAlgebra>::set_debug;

protected:
	using DebugWritingObject<TAlgebra>::write_debug;
	using DebugWritingObject<TAlgebra>::debug_writer;

};



}  // namespace ug
#endif
