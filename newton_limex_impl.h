/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit,
 *         but largely copied from ugcore Newton implementation by Andreas Vogel
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

#include <sstream>

#include "lib_disc/function_spaces/grid_function_util.h"
//#include "common/util/string_util.h"
#include "newton_limex.h"


#define PROFILE_NEWTON
#ifdef PROFILE_NEWTON
	#define NEWTON_PROFILE_FUNC()		PROFILE_FUNC_GROUP("Newton")
	#define NEWTON_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "Newton")
	#define NEWTON_PROFILE_END()		PROFILE_END()
#else
	#define NEWTON_PROFILE_FUNC()
	#define NEWTON_PROFILE_BEGIN(name)
	#define NEWTON_PROFILE_END()
#endif


namespace ug{

template <typename TAlgebra>
LimexNewtonSolver<TAlgebra>::
LimexNewtonSolver()
: m_spLinearSolver(NULL),
  m_N(NULL),
  m_J(NULL),
  m_spAss(NULL),
  m_linSolverSteps(0),
  m_linSolverRate(0.0)
{}


template <typename TAlgebra>
LimexNewtonSolver<TAlgebra>::
LimexNewtonSolver(SmartPtr<IOperator<vector_type> > N)
: m_spLinearSolver(NULL),
  m_N(NULL),
  m_J(NULL),
  m_spAss(NULL),
  m_linSolverSteps(0),
  m_linSolverRate(0.0)
{
	init(N);
}


template <typename TAlgebra>
LimexNewtonSolver<TAlgebra>::
LimexNewtonSolver(SmartPtr<IAssemble<TAlgebra> > spAss)
: m_spLinearSolver(NULL),
  m_N(new AssembledOperator<TAlgebra>(spAss)),
  m_J(NULL),
  m_spAss(spAss),
  m_linSolverSteps(0),
  m_linSolverRate(0.0)
{}


template <typename TAlgebra>
bool LimexNewtonSolver<TAlgebra>::init(SmartPtr<IOperator<vector_type> > N)
{
	NEWTON_PROFILE_BEGIN(NewtonLimexSolver_init);
	m_N = N.template cast_dynamic<AssembledOperator<TAlgebra> >();
	if (m_N.invalid())
		UG_THROW("NewtonLimexSolver: currently only works for AssembledDiscreteOperator.");

	m_spAss = m_N->discretization();
	return true;
}


template <typename TAlgebra>
bool LimexNewtonSolver<TAlgebra>::prepare(vector_type& u)
{
	return true;
}


template <typename TAlgebra>
bool LimexNewtonSolver<TAlgebra>::apply(vector_type& u)
{
	NEWTON_PROFILE_BEGIN(NewtonLimexSolver_apply);
	// check for linear solver
	if (m_spLinearSolver.invalid())
		UG_THROW("NewtonLimexSolver::apply: Linear solver not set.");

	// prepare Jacobian
	if (m_J.invalid() || m_J->discretization() != m_spAss)
		m_J = make_sp(new AssembledLinearOperator<TAlgebra>(m_spAss));

	m_J->set_level(m_N->level());

	// create tmp vectors
	SmartPtr<vector_type> spD = u.clone_without_values();
	SmartPtr<vector_type> spC = u.clone_without_values();

	// set Dirichlet values
	try {m_N->prepare(u);}
	UG_CATCH_THROW("NewtonLimexSolver::prepare: Operator preparation failed.");

	// compute defect
	try
	{
		NEWTON_PROFILE_BEGIN(NewtonComputeDefect1);
		m_N->apply(*spD, u);
		NEWTON_PROFILE_END();
	}
	UG_CATCH_THROW("NewtonLimexSolver::apply: Defect computation failed.");

	// increase offset of output for linear solver
	const int stdLinOffset = m_spLinearSolver->standard_offset();
	m_spLinearSolver->convergence_check()->set_offset(stdLinOffset + 3);

// perform exactly one Newton step
	// set c = 0
	NEWTON_PROFILE_BEGIN(NewtonSetCorretionZero);
	spC->set(0.0);
	NEWTON_PROFILE_END();

	// compute Jacobian
	try
	{
		NEWTON_PROFILE_BEGIN(NewtonComputeJacobian);
		m_J->init(u);
		NEWTON_PROFILE_END();
	}
	UG_CATCH_THROW("NewtonLimexSolver::apply: Jacobian initialization failed.");

	// init Jacobian inverse
	try
	{
		NEWTON_PROFILE_BEGIN(NewtonPrepareLinSolver);
		if (!m_spLinearSolver->init(m_J, u))
		{
			UG_LOGN("ERROR in 'NewtonLimexSolver::apply': Cannot init inverse linear "
				   "operator for Jacobi operator.");
			return false;
		}
		NEWTON_PROFILE_END();
	}
	UG_CATCH_THROW("NewtonLimexSolver::apply: Initialization of Linear Solver failed.");

	// solve linearized system
	try
	{
		NEWTON_PROFILE_BEGIN(NewtonApplyLinSolver);
		if (!m_spLinearSolver->apply(*spC, *spD))
		{
			UG_LOGN("ERROR in 'NewtonLimexSolver::apply': Cannot apply inverse linear "
					"operator for Jacobi operator.");
			return false;
		}
		NEWTON_PROFILE_END();
	}
	UG_CATCH_THROW("NewtonLimexSolver::apply: Application of Linear Solver failed.");

	// store convergence history
	m_linSolverSteps = m_spLinearSolver->step();
	m_linSolverRate = m_spLinearSolver->convergence_check()->avg_rate();

	// update solution
	u -= *spC;

	// apply constraints
	m_N->prepare(u);

	// reset offset of output for linear solver to previous value
	m_spLinearSolver->convergence_check()->set_offset(stdLinOffset);

	return true;
}


template <typename TAlgebra>
number LimexNewtonSolver<TAlgebra>::linear_solver_rate() const
{
	return m_linSolverRate;
}

template <typename TAlgebra>
int LimexNewtonSolver<TAlgebra>::linear_solver_steps() const
{
	return m_linSolverSteps;
}

template <typename TAlgebra>
std::string LimexNewtonSolver<TAlgebra>::config_string() const
{
	std::stringstream ss;
	ss << "NewtonLimexSolver\n";

	ss << " LinearSolver: ";
	if (m_spLinearSolver.valid())	ss << ConfigShift(m_spLinearSolver->config_string()) << "\n";
	else							ss << " NOT SET!\n";

	return ss.str();
}


}
