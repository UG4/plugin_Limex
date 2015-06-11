/*
 * time_integrator.hpp
 *
 *  Created on: 15.08.2014
 *      Author: anaegel
 */

#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/assemble_interface.h" // TODO: missing IAssemble in following file:
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"

// own headers
#include "time_extrapolation.h"
#include "../limex_tools.h"

namespace ug {


/// Integrates over a given time interval
template<class TDomain, class TAlgebra>
class ITimeIntegrator : public IOperator< GridFunction<TDomain, TAlgebra> >
{
	public:

		typedef TAlgebra algebra_type;
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename  TAlgebra::matrix_type matrix_type;


		typedef ITimeDiscretization<TAlgebra> time_disc_base_type;
		typedef ThetaTimeStep<TAlgebra> time_disc_type;


		typedef TDomain domain_type;
		typedef GridFunction<TDomain, TAlgebra> grid_function_type;
		typedef DomainDiscretization<TDomain, TAlgebra> domain_disc_type;

	protected:
		SmartPtr<domain_disc_type> m_spDomainDisc;

		double m_lower_tim;
		double m_upper_tim;
		double m_dt;
		number m_theta;

		SmartPtr<time_disc_type> m_spTimeDisc;        // set during init



	public:
		// constructor
		ITimeIntegrator(SmartPtr<domain_disc_type> domDisc)
		: m_spDomainDisc(domDisc),
		  m_lower_tim(0.0), m_upper_tim(0.0), m_dt(1.0),
		  m_theta(1.0)
		 {}

		/// virtual	destructor
		virtual ~ITimeIntegrator() {};

	///	init operator depending on a function u
	/**
	 * This method initializes the operator. Once initialized the 'apply'-method
	 * can be called. The function u is passed here, since the linear operator
	 * may be the linearization of some non-linear operator. Thus, the operator
	 * depends on the linearization point.
	 * If the operator is not a linearization, this method can be implemented
	 * by simply calling init() and forgetting about the linearization point.
	 *
	 * \param[in]	u		function (linearization point)
	 * \returns 	bool	success flag
	 */
	 virtual void init(grid_function_type const& u)
	 {
	 	UG_ASSERT(m_spDomainDisc.valid(), "TimeIntegrator<TDomain, TAlgebra>::init: m_spDomainDisc invalid.");
	 	m_spTimeDisc = make_sp(new time_disc_type(m_spDomainDisc, m_theta));
	 }


	///	init operator
	/**
	 * This method initializes the operator. Once initialized the 'apply'-method
	 * can be called.
	 * \returns 	bool	success flag
	 */
	  void init()
	  { UG_THROW("Please init with grid function!"); }

	///	prepares functions for application

	/*	This method is used to prepare the in- and output vector used later in apply.
	 * It can be called, after the init method has been called at least once.
	 * The prepare method is e.g. used to set dirichlet values.
	 */
	void prepare(grid_function_type& u) {}

	//! Apply operator
	/*! This method applies the operator, i.e, advances the time step*/
	void apply(grid_function_type& u1, const grid_function_type& u0)
	{
	  UG_THROW("Fix interfaces!");
	  // apply(u1, m_upper_tim, u0, m_lower_tim); 
	}

    virtual void apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0) = 0;

	//! Set initial time step
	void set_time_step(double dt)
	{ m_dt=dt;}

	void set_theta(double theta)
	{ m_theta=theta;}

	SmartPtr<time_disc_type> get_time_disc() {return m_spTimeDisc;}

};



/// integration of linear systems
template<class TDomain, class TAlgebra>
class ILinearTimeIntegrator : public ITimeIntegrator<TDomain, TAlgebra>
{

public:
	typedef ITimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;
	typedef typename base_type::vector_type vector_type;
	typedef IPreconditionedLinearOperatorInverse<vector_type> linear_solver_type;
	typedef AssembledLinearOperator<TAlgebra> assembled_operator_type;

	// forward constructor
	ILinearTimeIntegrator(SmartPtr<domain_disc_type> domDisc)
	: base_type(domDisc) {}

	void set_linear_solver(SmartPtr<linear_solver_type> lSolver)
	{ m_spLinearSolver=lSolver;}

protected:
	SmartPtr<linear_solver_type> m_spLinearSolver;

};




/// Integrate over a given time interval (for a linear problem)
template<class TDomain, class TAlgebra>
class LinearTimeIntegrator :
	public ILinearTimeIntegrator<TDomain, TAlgebra>
{
private:

protected:

public:
	typedef ILinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

	// constructor
	LinearTimeIntegrator (SmartPtr< typename base_type::domain_disc_type> domDisc)
	: base_type(domDisc) {}

	//void init(grid_function_type const& u);
	void apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
};



template<typename TDomain, typename TAlgebra>
void LinearTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{

	LIMEX_PROFILE_FUNC()

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	typename base_type::time_disc_type &tdisc = *base_type::m_spTimeDisc;

	// create solution vector & right hand side
	SmartPtr<typename base_type::vector_type> uold= u0->clone();
	SmartPtr<typename base_type::vector_type> b= u0->clone_without_values();

	// solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->clear();
	m_spSolTimeSeries->push(uold, t0);

	// create matrix operator
	SmartPtr<typename base_type::assembled_operator_type> spAssOp=make_sp(new typename base_type::assembled_operator_type(base_type::m_spTimeDisc, gl));

	// integrate
	UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"]\n");
	 double t = t0;
	 number dt_assembled = -1.0;   // invalid
	 int step = 1;

	 number currdt = base_type::m_dt;

	 while(1e-10*t1 < (t1-t))
	 {
		 UG_LOG("+++ Timestep +++" << step++ << "\n");
		 // determine step size
		 number dt = std::min(currdt, t1-t);

		 // prepare step
		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 if (fabs(dt-dt_assembled) > 1e-10*dt)
		 {
			 // re-assemble operator
			 UG_LOG("+++ Reassemble (t=" << t << ", dt=" << dt <<")\n");
			 tdisc.assemble_linear(*spAssOp, *b, gl);
			 (base_type::m_spLinearSolver)->init(spAssOp, *u1);
			 dt_assembled = dt;
		 }
		 else
		 {
			 // keep old operator
			 tdisc.assemble_rhs(*b, gl);
		 }

		 // execute step
		 if (base_type::m_spLinearSolver->apply(*u1, *b))
		 {
			 // ACCEPTING:
			 // push updated solution into time series
			 t += dt;
			 SmartPtr<typename base_type::vector_type> tmp = m_spSolTimeSeries->oldest();
			 VecAssign(*tmp, *u1);
			 m_spSolTimeSeries->push_discard_oldest(tmp, t);
		 }
		 else
		 {
			 // DISCARDING
			 currdt *= 0.5;
		 }

	 }

};


/// Integrate over a given time interval (for a linear problem)
template<class TDomain, class TAlgebra>
class ConstStepLinearTimeIntegrator :
	public ILinearTimeIntegrator<TDomain, TAlgebra>
{

public:
	typedef ILinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

private:

protected:

	int m_numSteps;
public:

	// constructor
	ConstStepLinearTimeIntegrator (SmartPtr< typename base_type::domain_disc_type> domDisc)
	: base_type(domDisc), m_numSteps(1) {}

	void apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
	void set_num_steps(int steps) {m_numSteps = steps;}
};



template<typename TDomain, typename TAlgebra>
void ConstStepLinearTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	LIMEX_PROFILE_FUNC()
	// short-cuts
	GridLevel const &gl = u0->grid_level();
	typename base_type::time_disc_type &tdisc = *base_type::m_spTimeDisc;

	// create solution vector & right hand side
	SmartPtr<typename base_type::vector_type> uold= u0->clone();
	SmartPtr<typename base_type::vector_type> b= u0->clone_without_values();

	// solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->push(uold, t0);

	SmartPtr<typename base_type::assembled_operator_type> spAssOp = SPNULL;

	// select number of steps
	 double t = t0;
	 int numSteps = round((t1-t0) / base_type::m_dt);
	 number currdt = (t1-t0) / numSteps;
	
	 //std::cerr << "+++ Integrating: ["<< t0 <<", "<< t1 <<"] with dt=" << currdt << "("<< numSteps<< " iters)\n";
	 UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"] with dt=" << currdt << "("<< numSteps<< " iters)");
	 
	 // integrate
	 for(int step = 1; step<=numSteps; ++step)
	 {
	         // determine step size
	         // number dt = std::min(currdt, t1-t);
		 const number dt = currdt;
		 UG_LOG("+++ Const timestep +++" << step<< "t=" << t << "->" << dt);// std::endl;
		
		 // prepare step
		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 if (spAssOp==SPNULL)
		 {
			 // assemble operator
			 UG_LOG("+++ Assemble (t=" << t << ", dt=" << dt <<")\n");
			 spAssOp=make_sp(new typename base_type::assembled_operator_type(base_type::m_spTimeDisc, gl));
			 tdisc.assemble_linear(*spAssOp, *b, gl);
			 (base_type::m_spLinearSolver)->init(spAssOp, *u1);
		 }
		 else
		 {
			 // std::cerr << "Recycling timestep " << step << "\n";
			 // keep old operator
			 tdisc.assemble_rhs(*b, gl);
		 }

		 // execute step
		 if (base_type::m_spLinearSolver->apply(*u1, *b))
		 {
			 // ACCEPTING:
			 // push updated solution into time series
			 t += dt;
			 SmartPtr<typename base_type::vector_type> tmp = m_spSolTimeSeries->oldest();
			 VecAssign(*tmp, *u1);
			 m_spSolTimeSeries->push_discard_oldest(tmp, t);
		 }
		 else
		 {
			UG_THROW("Const Time step failed!!!")
		 }

	 }

};



/// Integrate over a given time interval (for a linear problem)
template<class TDomain, class TAlgebra>
class TimeIntegratorLinearAdaptive :
	public ILinearTimeIntegrator<TDomain, TAlgebra>
{
public:
	typedef ILinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

private:


protected:

	SmartPtr<typename base_type::time_disc_type> m_spTimeDisc2;        // set during init

  double m_tol, m_dtmin, m_dtmax;


public:
  TimeIntegratorLinearAdaptive (SmartPtr< typename base_type::domain_disc_type> domDisc) : base_type(domDisc), m_tol(1e-2), m_dtmin(-1.0), m_dtmax(-1.0)
	{
		m_spTimeDisc2 = make_sp(new typename base_type::time_disc_type(base_type::m_spDomainDisc, base_type::m_theta));
	}

	void init(grid_function_type const& u);
	void apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

    void set_tol(double tol) {m_tol = tol;}
    void set_time_step_min(number dt) {m_dtmin = dt;}
    void set_time_step_max(number dt) {m_dtmax = dt;}
};


template<typename TDomain, typename TAlgebra>
void TimeIntegratorLinearAdaptive<TDomain, TAlgebra>::init(grid_function_type const& u)
{
	// call base
	base_type::init(u);
	m_spTimeDisc2 = make_sp(new typename base_type::time_disc_type(base_type::m_spDomainDisc, base_type::m_theta));
}

template<typename TDomain, typename TAlgebra>
void TimeIntegratorLinearAdaptive<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	// short-cuts
	GridLevel const &gl = u0->grid_level();
	typename base_type::time_disc_type &tdisc = *base_type::m_spTimeDisc;
	typename base_type::time_disc_type &tdisc2 = *m_spTimeDisc2;

	// create solution vector & right hand side
	SmartPtr<typename base_type::vector_type> uold= u0->clone();
	SmartPtr<typename base_type::vector_type> b= u0->clone_without_values();

	// create matrix operator
	SmartPtr<typename base_type::assembled_operator_type> spAssOp=make_sp(new typename base_type::assembled_operator_type(base_type::m_spTimeDisc, gl));

	// create additional solutions
	SmartPtr<typename base_type::grid_function_type>  u2old = u0->clone();
	SmartPtr<typename base_type::grid_function_type>  u2 = u0->clone();

	// solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->push(uold, t0);

	SmartPtr<vector_time_series_type> m_spSolTimeSeries2;
	m_spSolTimeSeries2=make_sp(new vector_time_series_type());
	m_spSolTimeSeries2->push(u2old, t0);
	
	// automatic selection of min/max time step
	if (m_dtmin<0) m_dtmin= (t1-10)/1e+5; 
	if (m_dtmax<0) m_dtmax= (t1-10)/10; 

	// Aitken Neville extrapolation
	const size_t tsteps[2] = {1,2};
	std::vector<size_t> vsteps (tsteps, tsteps+2);
	AitkenNevilleTimex<typename base_type::vector_type> timex(vsteps);

	// integrate
	UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"]\n");
	 double t = t0;
	 int step = 0;

	 number dt = base_type::m_dt;
	 while(t<t1)
	 {
	   // step: t -> t+dt
	   bool bSuccess = false;
	   while (!bSuccess){
		
	     // determine step size
	     if (dt<m_dtmin){
	       UG_LOG("Step size below minimum")
		 }
		 dt = std::min(dt, t1-t);
	   
		 // basic step
		 UG_LOG("+++ Timestep: " << ++step << "\n");
		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 tdisc.assemble_linear(*spAssOp, *b, gl);
		 base_type::m_spLinearSolver->init(spAssOp, *u1);
		 base_type::m_spLinearSolver->apply(*u1, *b);
		 
		 
		 // control 1/2
		 UG_LOG("+++ Control: " << step << "\n");
		 tdisc2.prepare_step(m_spSolTimeSeries, 0.5*dt);
		 tdisc2.assemble_linear(*spAssOp, *b, gl);
		 base_type::m_spLinearSolver->init(spAssOp, *u2);
		 base_type::m_spLinearSolver->apply(*u2, *b);
		 
		 // push back solution
		 SmartPtr<typename base_type::vector_type> tmp2 =  m_spSolTimeSeries2->oldest();
		 (*tmp2)=*(u2.template cast_static<typename base_type::vector_type>());
		 m_spSolTimeSeries2->push_discard_oldest(tmp2, t + 0.5*dt);
		 
		 // control 2/2
		 tdisc2.prepare_step(m_spSolTimeSeries2, 0.5*dt);
		 tdisc2.assemble_linear(*spAssOp, *b, gl);
		 base_type::m_spLinearSolver->init(spAssOp, *u2);
		 base_type::m_spLinearSolver->apply(*u2, *b);

		 // obtain extrapolated solution
		 timex.set_solution(u1, 0);
		 timex.set_solution(u2, 1);
		 timex.apply();

		 // predict (subsequent) time step
		 number eps = timex.get_error_estimate();
		 number lambda = std::pow(0.8* m_tol/eps, 0.5); 

		 number dtEst= dt*lambda; 
		 dtEst = std::min(dtEst, 1.5*dt); 
		 dtEst = std::min(dtEst, m_dtmax); 
		 
		 dt = dtEst;
		 if (eps <= m_tol)
		 {
		   // ACCEPT STEP (and thus solution u2)
		   UG_LOG("ACCEPTING solution, dtnew=" << dt);
		   bSuccess = true;
		 }
		 else
		 {
		   // DISCARD step
		   UG_LOG("DISCARDING solution, dtnew=" << dt);

		   // => reset solutions
		   VecAssign(*u1, *uold);
		   
		   // => timeSeries2 has been updated...
		   SmartPtr<typename base_type::vector_type> tmp = m_spSolTimeSeries2->oldest();
		   VecAssign(*tmp, *uold);
		   m_spSolTimeSeries2->push_discard_oldest(tmp, t);

		 }

	   } 

	 
	   // prepare next loop
	   t += dt;

	   // push updated solution into time series (and continue)
	   SmartPtr<typename base_type::vector_type> tmp = m_spSolTimeSeries->oldest();
	   VecAssign(*tmp, *u2);
	   m_spSolTimeSeries->push_discard_oldest(tmp, t);

	 }

};

/// integration of non-linear systems
template<class TDomain, class TAlgebra>
class INonlinearTimeIntegrator : public ITimeIntegrator<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;
	typedef typename base_type::vector_type vector_type;
	typedef IOperatorInverse<vector_type> solver_type;
	typedef AssembledOperator<TAlgebra> assembled_operator_type;

	INonlinearTimeIntegrator(SmartPtr<domain_disc_type> domDisc)
	: base_type(domDisc) {}

	void set_solver(SmartPtr<solver_type> solver)
	{ m_spSolver=solver;}
protected:
	SmartPtr<solver_type> m_spSolver;
};

/// integration of non-linear systems with steo
template<class TDomain, class TAlgebra>
class INonlinearTimeIntegratorWithBounds : public INonlinearTimeIntegrator<TDomain, TAlgebra>
{
public:
	typedef INonlinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;

	INonlinearTimeIntegratorWithBounds(SmartPtr<domain_disc_type> domDisc)
	: base_type(domDisc), m_dtMin(0.0), m_dtMax(std::numeric_limits<double>::max()),  m_redFac(1.0), m_incFac(1.0) {}

	void set_dt_min(double min) { m_dtMin = min; }
	void set_dt_max(double max) { m_dtMax = max; }
	void set_decrease_factor(double dec) { m_redFac = dec; }
	void set_increase_factor(double inc) { m_incFac = inc; }

protected:
	double m_dtMin, m_dtMax;
	double m_redFac, m_incFac;
};



/// Integrate (a non-linear problem) over a given time interval
template<class TDomain, class TAlgebra>
class SimpleTimeIntegrator :
		public INonlinearTimeIntegratorWithBounds<TDomain, TAlgebra>
{

public:
	typedef INonlinearTimeIntegratorWithBounds<TDomain, TAlgebra> base_type;
	typedef typename base_type::domain_disc_type domain_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

	// constructor
	SimpleTimeIntegrator (SmartPtr< typename base_type::domain_disc_type> domDisc)
	: base_type(domDisc) {}

	void apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

};


template<typename TDomain, typename TAlgebra>
void SimpleTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{

	LIMEX_PROFILE_FUNC()

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	typename base_type::time_disc_type &tdisc = *base_type::m_spTimeDisc;
	typename base_type::solver_type &solver = *base_type::m_spSolver;

	// create solution vector & right hand side
	SmartPtr<typename base_type::vector_type> uold= u0->clone();

	// init solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->clear();
	m_spSolTimeSeries->push(uold, t0);

	// init solver (and matrix operator)
	SmartPtr<typename base_type::assembled_operator_type> spAssOp=make_sp(new typename base_type::assembled_operator_type(base_type::m_spTimeDisc, gl));
	solver.init(spAssOp);

	// integrate
	double t = t0;
	number currdt = base_type::m_dt;
	int step = 1;

	 UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"] with " << currdt <<"\n");


	 bool failed = false;
	 while(1e-10*t1 < (t1-t) && !failed)
	 {
		 UG_LOG("+++ Timestep +++" << step++ << "\n");

		 // determine step size
		 UG_COND_THROW(currdt < base_type::m_dtMin, "Time step size below minimum. ABORTING!")
		 number dt = std::min(currdt, t1-t);

		 // init step
		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 if (solver.prepare(*u1) == false)
		 {
			 UG_LOG("Initialzation failed! RETRY");
			 currdt *= base_type::m_redFac;
			 continue;
		 };

		 // execute step
		 if (solver.apply(*u1))
		 {
			 // ACCEPT step
			 // update time
			 t += dt;

			 // push updated solution into time series (and continue)
			 SmartPtr<typename base_type::vector_type> tmp = m_spSolTimeSeries->oldest();
			 VecAssign(*tmp, *u1);
			 m_spSolTimeSeries->push_discard_oldest(tmp, t);

		 }
		 else
		 {
			 // REJECT step
			 UG_LOG("Solution failed! RETRY");
			 currdt *= base_type::m_redFac;
			 continue;
		 };




		// tdisc.finish_step_elem(m_spSolTimeSeries, dt);



	 }

};


} // namespace ug

#endif /* TIME_INTEGRATOR_HPP_ */
