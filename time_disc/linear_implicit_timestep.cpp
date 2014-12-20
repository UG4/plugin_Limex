/*
 * linear_implicit_timestep.cpp
 *
 *  Created on: 20.05.2014
 *      Author: anaegel
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

	//m_JLinOp = make_sp(new AssembledLinearOperator<TAlgebra>(m_spJAss));
	m_JLinOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spDomDisc));
	m_JLinOp->init(*m_pPrevSol->oldest());

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

		this->m_spDomDisc->prepare_timestep(m_pPrevSol, gl);
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot prepare timestep.");

// Aux linear operator
	m_JLinOp = make_sp(new AssembledLinearOperator<TAlgebra>(this->m_spDomDisc));
	m_JLinOp->init(*m_pPrevSol->oldest());
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
	}UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble jacobian.");

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

// 	future solution part
	try{

		// d:=\tau {f(k-1) - A(k-1) u(k-1)}
		std::vector<number> vScaleMass(1, 0.0); // 0
		std::vector<number> vScaleStiff(1, m_dt);
		this->m_spDomDisc->assemble_defect(d, m_pPrevSol, vScaleMass, vScaleStiff, gl);

		// d:=d+(M + \tau J) (u(k-1)-u)
		vector_type deltau = *m_pPrevSol->oldest()->clone();
		deltau -= u;

		this->m_spDomDisc->assemble_jacobian(m_JLinOp->get_matrix(), m_pPrevSol, m_dt, gl);
		m_JLinOp->apply_sub(d, deltau);

	}UG_CATCH_THROW("LinearImplicitEuler: Cannot assemble defect.");

	//std::cout << "DEFECT: "<< m_pPrevSol->time(0) << std::endl;

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

#ifdef UG_CPU_1
template class LinearImplicitEuler<CPUAlgebra>;
#endif
#ifdef UG_CPU_2
template class LinearImplicitEuler<CPUBlockAlgebra<2> >;
#endif
#ifdef UG_CPU_3
template class LinearImplicitEuler<CPUBlockAlgebra<3> >;
#endif
#ifdef UG_CPU_4
template class LinearImplicitEuler<CPUBlockAlgebra<4> >;
#endif

}; // namespace ug
