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

#include "linear_implicit_timestep.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"

namespace ug{

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
             number dt)
{
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
		this->m_spMatrixDisc->prepare_timestep(m_pPrevSol, m_futureTime);

		this->m_spGammaDisc->prepare_timestep(m_pPrevSol, m_futureTime);
	}
	UG_CATCH_THROW("ThetaTimeStep: Cannot prepare time step.");

	//m_JLinOp = make_sp(new AssembledLinearOperator<TAlgebra>(m_spJAss));

	m_JLinOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spMatrixDisc));

	if (m_spGammaDisc != SPNULL && m_spGammaOp == SPNULL)
	{ m_spGammaOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spGammaDisc)); }
	/*{
		m_JLinOp->init(*m_pPrevSol->oldest());
	}
*/
}

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
                  number dt, const GridLevel& gl)
{
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
		this->m_spMatrixDisc->prepare_timestep(m_pPrevSol, m_futureTime, gl);
	} UG_CATCH_THROW("LinearImplicitEuler: Cannot prepare timestep.");

	// Aux linear operator
	m_JLinOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spMatrixDisc));


	if (m_spGammaDisc != SPNULL && m_spGammaOp == SPNULL)
	{ m_spGammaOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spGammaDisc)); }

	// m_JLinOp->init(*m_pPrevSol->oldest());
	//std::cout << "PREPELEM: "<< m_vScaleMass[0] <<", " << m_vScaleStiff[0] << ", " <<m_dt << ", " << m_pPrevSol->time(0) << std::endl;

}

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
adjust_solution(vector_type& u, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(LinearImplicitEuler_adjust_solution, "discretization LinearImplicitEuler");
//	adjust solution
	try{
		this->m_spDomDisc->adjust_solution(u, m_futureTime, gl);
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot adjust solution.");
	//std::cout << "ADJUST: "<< m_futureTime << std::endl;

}



template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
                 const GridLevel& gl)
{
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
 *  */


/** WARNING: This function is abused
 * Must return: $Mk + \tau J$
 * */
template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl)
{
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

	// assemble "jacobian"  using current iterate
	try{
		// (M_{k-1} + \tau J)
		// WARNING: This would only work for constant mass matrix !!!
		this->m_spDomDisc->assemble_jacobian(J, m_pPrevSol, m_dt, gl);


		/* Add (M_k-M_{k-1})
		matrix_type M;
		M=J;

		this->m_spDomDisc->assemble_mass_matrix(M, u, gl);
		MatAdd(J, 1.0, J, 1.0, M);

		this->m_spDomDisc->assemble_mass_matrix(M, *m_pPrevSol->oldest(), gl);
		MatAdd(J, 1.0, J, -1.0, M);
		 */

		// (Re-)assemble Gamma
		if (m_bGammaNeedsUpdate == true)
		{
			UG_LOG("Assembling GAMMA")
			this->m_spGammaDisc->assemble_jacobian(m_spGammaOp->get_matrix(), m_pPrevSol, m_dt, gl);
			m_bGammaNeedsUpdate = false;
		}

		//
		if (m_spGammaDisc != SPNULL)
		{
			UG_ASSERT(m_spGammaOp != SPNULL, "Huhh: No operator??? ");
			UG_LOG("Adding GAMMA")
			MatAdd(J, 1.0, J, 1.0, m_spGammaOp->get_matrix());
		}


	} UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble jacobian.");

	//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

/** WARNING: This function is abused
 * Must return : d_A(k-1):= tau * F(k-1) - A(k-1) u(k-1) */
template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl)
{
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
		/*	vector_type deltau = *m_pPrevSol->oldest()->clone();
		deltau -= u;

		this->m_spDomDisc->assemble_jacobian(m_JLinOp->get_matrix(), m_pPrevSol, m_dt, gl);
		m_JLinOp->apply_sub(d, deltau);
		*/
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble defect.");

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
 *  */

template <typename TAlgebra>
void LinearImplicitEuler<TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl)
{

	PROFILE_BEGIN_GROUP(LinearImplicitEuler_assemble_linear, "discretization LinearImplicitEuler");

//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
	UG_THROW("LinearImplicitEuler::assemble_defect:"
			" Number of previous solutions must be at least "<<
			m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

	// 	future solution part
		try{
			// A:= (M(k-1) + tau A(k-1) )
			this->m_spDomDisc->assemble_jacobian(A, m_pPrevSol, m_dt, gl);

			std::vector<number> vScaleMass(1, 0.0);
			std::vector<number> vScaleStiff(1, m_dt);
			this->m_spDomDisc->assemble_defect(b, m_pPrevSol, vScaleMass, vScaleStiff, gl);

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
		UG_THROW("LinearImplicitEuler::assemble_linear:"
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
