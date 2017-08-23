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

#ifndef UG_PLUGINS__LIMEX__NEWTON_LIMEX
#define UG_PLUGINS__LIMEX__NEWTON_LIMEX

#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/debug_writer.h"

#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"


namespace ug {

/// Newton solver for assembling-based discretizations solved using Limex
template <typename TAlgebra>
class LimexNewtonSolver
	: 	public IOperatorInverse<typename TAlgebra::vector_type>
		//public DebugWritingObject<TAlgebra>
{
	public:
		/// algebra type
		typedef TAlgebra algebra_type;

		/// vector type
		typedef typename TAlgebra::vector_type vector_type;

		/// matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
		/// default constructor
		LimexNewtonSolver();

		/// constructor setting operator
		LimexNewtonSolver(SmartPtr<IOperator<vector_type> > N);

		/// constructor using assembling
		LimexNewtonSolver(SmartPtr<IAssemble<TAlgebra> > spAss);

		/// sets the linear solver
		void set_linear_solver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver)
		{m_spLinearSolver = LinearSolver;}

		/// This operator inverts the operator N: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > N);

		/// prepare operator
		virtual bool prepare(vector_type& u);

		/// apply operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

		/**
		 * @brief Returns information about configuration parameters.
		 * This should return necessary information about parameters and possibly
		 * calling config_string of subcomponents.
		 *
		 * @returns std::string  necessary information about configuration parameters
		 */
		virtual std::string config_string() const;

		/// prints average linear solver convergence
		number linear_solver_rate() const;

		/// information on linear solver convergence
		int linear_solver_steps() const;

	private:
		/// help functions for debug output
		/// @{
		void write_debug(const vector_type& vec, const char* filename);
		void write_debug(const matrix_type& mat, const char* filename);
		/// @}

	private:
		/// linear solver
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spLinearSolver;


		/// assembling routine
		SmartPtr<AssembledOperator<algebra_type> > m_N;

		/// jacobi operator
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J;

		/// assembling
		SmartPtr<IAssemble<TAlgebra> > m_spAss;


		/// convergence history of linear solver
		/// @{
		size_t m_linSolverSteps;
		number m_linSolverRate;
		/// @}
};

}	// namespace ug

#include "newton_limex_impl.h"

#endif  // UG_PLUGINS__LIMEX__NEWTON_LIMEX
