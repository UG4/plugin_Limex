
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


// ug4 headers
#include "common/assert.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"
#include "lib_disc/function_spaces/grid_function_util.h"

// plugin headers
#include "linear_implicit_timestep.h"
#include "../limex_tools.h"

namespace ug{


template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
             number dt)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:prepare_step" << std::endl);

	PROFILE_BEGIN_GROUP(LinearImplicitEuler_prepare_step, "discretization LinearImplicitEuler");
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("LinearImplicitEuler::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol->size() << " passed.\n");

//	remember old values
	m_pPrevSol = prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              m_pPrevSol);

	std::cout << "PREP: "<< m_vScaleMass[0] <<", " << m_vScaleStiff[0] << ", " <<m_dt << ", " << m_pPrevSol->time(0) << std::endl;

//	prepare time step (elemDisc-wise)
	try
	{
		this->m_spDomDisc->prepare_timestep(m_pPrevSol, m_futureTime);
		this->m_spMatrixJDisc->prepare_timestep(m_pPrevSol, m_futureTime);

		if (m_spGammaDisc.valid())
		{ m_spGammaDisc->prepare_timestep(m_pPrevSol, m_futureTime);}
	}
	UG_CATCH_THROW("LinearlyImplicitEuler: Cannot prepare time step.");

	// create matrix J
	if (m_spMatrixJOp.invalid())
	{ m_spMatrixJOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spMatrixJDisc)); }

	// create matrix Gamma
	if (m_spGammaDisc != SPNULL && m_spGammaOp == SPNULL)
	{ m_spGammaOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spGammaDisc)); }

	/*{
		m_spMatrixJOp->init(*m_pPrevSol->oldest());
	}
*/
}

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
                  number dt, const GridLevel& gl)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:prepare_step_elem" << std::endl);

	PROFILE_BEGIN_GROUP(LinearImplicitEuler_step_elem, "discretization LinearImplicitEuler");
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("LinearImplicitEuler::prepare_step_elem:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol->size() << " passed.\n");

//	remember old values
	m_pPrevSol = prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              m_pPrevSol);
// 	prepare timestep
	try{
		this->m_spDomDisc->prepare_timestep(m_pPrevSol, m_futureTime, gl);
		this->m_spMatrixJDisc->prepare_timestep(m_pPrevSol, m_futureTime, gl);
	} UG_CATCH_THROW("LinearImplicitEuler: Cannot prepare timestep.");

	// Aux linear operator
	if (m_spMatrixJOp.invalid()) {
		m_spMatrixJOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spMatrixJDisc));
	}


	if (m_spGammaDisc.valid() && m_spGammaOp == SPNULL) {
		m_spGammaOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spGammaDisc));
	}

	// m_spMatrixJOp->init(*m_pPrevSol->oldest());
	//std::cout << "PREPELEM: "<< m_vScaleMass[0] <<", " << m_vScaleStiff[0] << ", " <<m_dt << ", " << m_pPrevSol->time(0) << std::endl;

}

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
adjust_solution(vector_type& u, const GridLevel& gl)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:adjust_solution" << &u << std::endl);
	PROFILE_BEGIN_GROUP(LinearImplicitEuler_adjust_solution, "discretization LinearImplicitEuler");

	//	adjust solution
	try{
		this->m_spDomDisc->adjust_solution(u, m_futureTime, gl);
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot adjust solution.");

}



template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
                 const GridLevel& gl)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:finish_step_elem" << std::endl);
	
//	perform checks whether 'currSol' is a solutionTimeSeries only with the new values
	if(currSol->time(0) != m_futureTime)
		UG_THROW("LinearImplicitEuler::finish_step_elem:"
				" The solution of the SolutionTimeSeries used in this function"
				" does not coincide with the current solution! ");

	// 	finish timestep using the current solution
	try{
		this->m_spDomDisc->finish_timestep(currSol, gl);
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot finish timestep.");
}


/*
 *
 *  Non-linear system
 *
 *	\partial m(u) = f(u)
 *
 *  */


/** WARNING: This function is abused
 * Must return: $Mk + \tau J$
 * */
template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_jacobian(matrix_type& J_limex, const vector_type& u, const GridLevel& gl)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:assemble_jacobian" << std::endl);

	PROFILE_BEGIN_GROUP(LinearImplicitEuler_assemble_jacobian, "discretization LinearImplicitEuler");

	//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("LinearImplicitEuler::assemble_jacobian:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

	//	push unknown solution to solution time series
	//	ATTENTION: Here, we must cast away the constness of the solution, but note,
	//			   that we pass pPrevSol as a const object in assemble_... Thus,
	//			   the solution will not be changed there and we pop it from the
	//			   Solution list afterwards, such that nothing happens to u
		// \todo: avoid this hack, use smart ptr properly
		int DummyRefCount = 2;
		SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
		m_pPrevSol->push(pU, m_futureTime);

	// assemble "Jacobian" (M_{k-1} + \tau J) using current iterate
	try{

		// check assertions
		UG_ASSERT(m_spMatrixJDisc.valid(), "Huhh: Invalid matrix discretization")

		if (!m_useCachedMatrices)
		{
			// un-cached, (i.e., approximate Newton)
			this->m_spMatrixJDisc->assemble_jacobian(J_limex, m_pPrevSol, m_dt, gl);
			UG_DLOG(LIB_LIMEX, 5, "Computed J_limex =" << J_limex << " for dt=" << m_dt << std::endl);
			write_debug(J_limex, "myMatrixAssembled.mat");
		}
		else
		{
			/* Assemble (Mk + \tau J)
			 *
			 * Note: Mk has Dirichlet rows (e.g., for inout bnd cond)
			 */
			// J_limex.set(0.0);
			this->m_spMatrixJDisc->assemble_jacobian(J_limex, m_pPrevSol, 0.0, gl);
			UG_DLOG(LIB_LIMEX, 3, "> Computed Mk (" << &J_limex << " : "<<J_limex <<
						" at " << m_pPrevSol->oldest_time() << ", " << GetNNZs(J_limex) << " nonzeros)" << std::endl);
			write_debug(J_limex, "myMk.mat");

			const double mydt = 1.0;    // use time step 1.0 for assembly
			matrix_type &J_stiff = m_spMatrixJOp->get_matrix();

			if (m_bMatrixJNeedsUpdate)
			{
				// First part of J: -df/du
				this->m_spMatrixJDisc->assemble_jacobian(J_stiff, m_pPrevSol, mydt, gl);
				UG_DLOG(LIB_LIMEX, 3, "> Cached  J_0 = -df/du (" << &J_stiff << ":"<< J_stiff<< ")"<< std::endl);
				write_debug(J_stiff, "myStiff0.mat");

				// subtracting mass matrix yields J
				MatAddNonDirichlet<matrix_type>(J_stiff, 1.0, J_stiff, -1.0, J_limex);
				write_debug(J_stiff, "myStiff1.mat");

				// Second part of J: Gamma update (if applicable)
				if (m_spGammaDisc.valid())
				{
					UG_ASSERT(m_spGammaOp != SPNULL, "Huhh: No operator??? ");
					matrix_type &myGamma = m_spGammaOp->get_matrix();
					this->m_spGammaDisc->assemble_jacobian(myGamma, m_pPrevSol, mydt, gl);
					UG_DLOG(LIB_LIMEX, 3, "Assembled Gamma0 ("<< myGamma);
					UG_DLOG(LIB_LIMEX, 3,  " at " << m_pPrevSol->oldest_time() <<", " << GetNNZs(myGamma) << " nnz)" << std::endl);
					write_debug(myGamma, "myGamma.mat");

					MatAddNonDirichlet<matrix_type>(J_stiff, 1.0, J_stiff, 1.0, myGamma);
					write_debug(J_stiff, "myStiff2.mat");
					UG_DLOG(LIB_LIMEX, 3, "Added Gamma0"<< std::endl)
				}

				// do not need to recompute
				m_bMatrixJNeedsUpdate = false;
			}

			write_debug(J_stiff, "myStiffX.mat");

			// Updating Jtotal = Mk+ \tau Ja (for non-dirichlet)
			UG_ASSERT (J_limex.num_rows() == J_stiff.num_rows(), "Huhh: Row dimension does not match: "
						<< J_limex.num_rows() <<"!=" << J_stiff.num_rows());
			UG_ASSERT (J_limex.num_cols() == J_stiff.num_cols(), "Huhh: Col dimension does not match:"
						<< J_limex.num_cols() <<"!=" << J_stiff.num_cols());

			// Note: J_limex has Dirichlet values
			MatAddNonDirichlet<matrix_type>(J_limex, 1.0, J_limex, m_dt, J_stiff);
			UG_DLOG(LIB_LIMEX, 3, "Computed $J=Mk+ tau J_a$ (" << J_limex << " with tau=" << m_dt);
			UG_DLOG(LIB_LIMEX, 3,  " at " << m_pPrevSol->oldest_time() <<", " << GetNNZs(J_limex) << " nonzeros)" << std::endl);

			write_debug(J_limex, "myMatrixCached.mat");
		} // if(! m_useCachedMatrices)

	}
	UG_CATCH_THROW("LinearlyImplicitEuler: Cannot assemble jacobian(s).");

	//UG_ASSERT(0, "Crashing")

	//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

/** WARNING: This function is abused
 * Must return : d_A(k-1):= tau * F(k-1) - A(k-1) u(k-1) */
template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:assemble_defect" << std::endl);

	PROFILE_BEGIN_GROUP(LinearImplicitEuler_assemble_defect, "discretization LinearImplicitEuler");

	//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("LinearImplicitEuler::assemble_defect:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

	//	push unknown solution to solution time series
	//	ATTENTION: Here, we must cast away the constness of the solution, but note,
	//			   that we pass pPrevSol as a const object in assemble_... Thus,
	//			   the solution will not be changed there and we pop it from the
	//			   Solution list afterwards, such that nothing happens to u
    // \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
	m_pPrevSol->push(pU, m_futureTime);

// 	future solution part
	try{

		// d:=\tau {f(k-1) - A(k-1) u(k-1)}
		std::vector<number> vScaleMass(1, 0.0); // 0
		std::vector<number> vScaleStiff(1, m_dt);
		this->m_spDomDisc->assemble_defect(d, m_pPrevSol, vScaleMass, vScaleStiff, gl);

		// d := d + (M - \tau J) (u(k-1)-u)
		/*
		vector_type deltau = *m_pPrevSol->oldest()->clone();
		deltau -= u;

		this->m_spDomDisc->assemble_jacobian(m_spMatrixJOp->get_matrix(), m_pPrevSol, m_dt, gl);
		m_spMatrixJOp->apply_sub(d, deltau);
		*/
	}UG_CATCH_THROW("LinearlyImplicitEuler: Cannot assemble defect.");

	//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();

}


/* RHS (non-linear system) */
template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl)
{
	UG_ASSERT(0, "Really wanna use me???");
	PROFILE_BEGIN_GROUP(LinearImplicitEuler_assemble_rhs, "discretization LinearImplicitEuler");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("LinearImplicitEuler::assemble_linear:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_rhs(b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble jacobian.");

}

/*
 *
 *  Linear system
 *
 *
 */

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_linear(matrix_type& J, vector_type& b, const GridLevel& gl)
{
	UG_DLOG(LIB_LIMEX, 5, "LinearlyImplicitEuler:assemble_linear" << std::endl);

	PROFILE_BEGIN_GROUP(LinearImplicitEuler_assemble_linear, "discretization LinearImplicitEuler");

//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
	UG_THROW("LinearImplicitEuler::assemble_defect:"
			" Number of previous solutions must be at least "<<
			m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

	// 	future solution part
		try{
			UG_ASSERT(m_pPrevSol->size()<2, "Only need one solution,,,");

			// J(u_{k+1} - u_k) = -\tau f_k
			// NOTE: both routines have invoked the contraints, so constraint adjustment are assumed to be correct!
			this->assemble_jacobian(J, *m_pPrevSol->oldest(), gl);
			this->assemble_defect(b, *m_pPrevSol->oldest(), gl);

			// J u_{k+1} = -\tau f_k + J u_k
			//AxpyCommonSparseMatrix(J, b, 0.0, *m_pPrevSol->oldest(), +1.0, *m_pPrevSol->oldest());
			J.axpy(b, -1.0, b, 1.0, *m_pPrevSol->oldest());
		}UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble defect.");
}



/* RHS (linear system) */
template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_rhs(vector_type& b, const GridLevel& gl)
{
	UG_ASSERT(0, "Really wanna use me???");
	PROFILE_BEGIN_GROUP(LinearImplicitEuler_assemble_rhs, "discretization LinearImplicitEuler");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("LinearImplicitEuler::assemble_rhs:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_rhs(b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble jacobian.");

}


////////////////////////////////////////////////////////////////////////
//	template instantiations for all current algebra types.

UG_ALGEBRA_CPP_TEMPLATE_DEFINE_ALL(LinearImplicitEuler)


}; // namespace ug
