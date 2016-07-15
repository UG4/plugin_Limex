/*
 * linear_implicit_timestep.h
 *
 *  Created on: 20.05.2014
 *      Author: anaegel
 *
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
: public ITimeDiscretization<TAlgebra>
	//	: public MultiStepTimeDiscretization<TAlgebra>
{
public:
/// Type of algebra
	typedef TAlgebra algebra_type;

/// Type of algebra matrix
	typedef typename algebra_type::matrix_type matrix_type;

/// Type of algebra vector
	typedef typename algebra_type::vector_type vector_type;

/// Domain Discretization type
	typedef IDomainDiscretization<algebra_type>	domain_discretization_type;

public:
/// constructor
	LinearImplicitEuler(SmartPtr<IDomainDiscretization<algebra_type> > spDD)
		: ITimeDiscretization<TAlgebra>(spDD),
		  m_pPrevSol(NULL)
	{}

	virtual ~LinearImplicitEuler(){};

/// \copydoc ITimeDiscretization::num_prev_steps()
	virtual size_t num_prev_steps() const {return m_prevSteps;}

///	\copydoc ITimeDiscretization::prepare_step()
	virtual void prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
	                          number dt);

///	\copydoc ITimeDiscretization::prepare_step_elem()
	virtual void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
	                               number dt, const GridLevel& gl);

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

protected:

	virtual number update_scaling(std::vector<number>& vSM,
			                              std::vector<number>& vSA,
			                              number dt, number currentTime,
			                              ConstSmartPtr<VectorTimeSeries<vector_type> > prevSol)

	{
		// const number theta = 1.0;

		//	resize scaling factors
		vSM.resize(1);
		vSM[0] = 1.0;
		//vSM[1] = -1.0;

		vSA.resize(1);
		vSA[0] = dt;
		//vSA[1] = (1.0-theta) * dt;
		return currentTime + dt;
	}


	static const size_t m_prevSteps=1;			///< number of previous steps needed.
	std::vector<number> m_vScaleMass;			///< Scaling for mass part
	std::vector<number> m_vScaleStiff;			///< Scaling for stiffness part

	SmartPtr<VectorTimeSeries<vector_type> > m_pPrevSol;		///< Previous solutions
	SmartPtr<AssembledLinearOperator<algebra_type> > m_JLinOp;	///< Operator

	number m_dt; 								///< Time Step size
	number m_futureTime;						///< Future Time


};



}  // namespace ug
#endif
