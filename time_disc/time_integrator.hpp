/*
 * Copyright (c) 2014-2020:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel, Andreas Kreienbuehl
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

#ifndef __H__LIMEX__TIME_INTEGRATOR_HPP__
#define __H__LIMEX__TIME_INTEGRATOR_HPP__

#if __cplusplus >= 201103L
#define OVERRIDE override
#else
#define OVERRIDE
#endif

// std headers.
#include <string>

// UG4 headers
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
#include "lib_disc/time_disc/time_integrator_subject.hpp"
#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "lib_disc/function_spaces/grid_function_util.h" // SaveVectorForConnectionViewer
#include "lib_disc/function_spaces/interpolate.h" //Interpolate
#include "lib_disc/function_spaces/integrate.h" //Integral
#include "lib_disc/function_spaces/grid_function.h" //GridFunction


// Plugin headers
#include "time_extrapolation.h"
#include "../limex_tools.h"

namespace ug {

/// Integrates over a given time interval [a,b] with step size dt
template<class TDomain, class TAlgebra>
class ITimeIntegrator
		:	public IOperator< GridFunction<TDomain, TAlgebra> >,
			public TimeIntegratorSubject<TDomain, TAlgebra>
{
	public:

		typedef TAlgebra algebra_type;
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;


		typedef TDomain domain_type;
		typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	protected:
		double m_dt;
		double m_lower_tim;
		double m_upper_tim;

		double m_precisionBound;

		bool m_bNoLogOut;


	public:
		// constructor
		ITimeIntegrator()
		: m_dt(1.0), m_lower_tim(0.0), m_upper_tim(0.0), m_precisionBound(1e-10), m_bNoLogOut(false)
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
	 //	UG_ASSERT(m_spDomainDisc.valid(), "TimeIntegrator<TDomain, TAlgebra>::init: m_spDomainDisc invalid.");
	 //	m_spTimeDisc = make_sp(new time_disc_type(m_spDomainDisc, m_theta));
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
	{ UG_THROW("Fix interfaces!"); }

    virtual bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0) = 0;

	//! Set initial time step
	void set_time_step(double dt)
	{ m_dt = dt; }

	double get_time_step()
	{ return m_dt; }

	void set_precision_bound(double precisionBound)
	{ m_precisionBound = precisionBound; return; }

	void set_no_log_out(bool bNoLogOut)
	{ m_bNoLogOut = bNoLogOut; return; }

};


/// Integration of linear systems
template<class TDomain, class TAlgebra>
class ILinearTimeIntegrator : public ITimeIntegrator<TDomain, TAlgebra>
{

public:
	typedef ITimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::vector_type vector_type;
	typedef ILinearOperatorInverse<vector_type> linear_solver_type;
	typedef AssembledLinearOperator<TAlgebra> assembled_operator_type;

	// forward constructor
	ILinearTimeIntegrator()
	: base_type() {}

	ILinearTimeIntegrator(SmartPtr<linear_solver_type> lSolver)
	: base_type(),  m_spLinearSolver(lSolver)
	{}


	void set_linear_solver(SmartPtr<linear_solver_type> lSolver)
	{ m_spLinearSolver=lSolver;}

protected:
	SmartPtr<linear_solver_type> m_spLinearSolver;

};


/// ITimeDiscDependentObject
template<class TAlgebra>
class ITimeDiscDependentObject
{
public:
	typedef ITimeDiscretization<TAlgebra> time_disc_type;

	ITimeDiscDependentObject(SmartPtr<time_disc_type> spTimeDisc) :
		m_spTimeDisc(spTimeDisc)
	{}

	SmartPtr<time_disc_type> get_time_disc() {return m_spTimeDisc;}
protected:
	SmartPtr<time_disc_type> m_spTimeDisc;
};


/// Integrate over a given time interval (for a linear problem)
template<class TDomain, class TAlgebra>
class LinearTimeIntegrator :
	public ILinearTimeIntegrator<TDomain, TAlgebra>,
	public ITimeDiscDependentObject<TAlgebra>

{
private:

protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;
public:
	typedef ILinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef ITimeDiscretization<TAlgebra> time_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

	// constructor
	LinearTimeIntegrator (SmartPtr< time_disc_type> tDisc)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc) {}

	LinearTimeIntegrator (SmartPtr< time_disc_type> tDisc, 	SmartPtr<typename base_type::linear_solver_type> lSolver)
	: base_type(lSolver), ITimeDiscDependentObject<TAlgebra>(tDisc) {}

	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
};


/// Integrate over a given time interval (for a linear problem)
template<class TDomain, class TAlgebra>
class ConstStepLinearTimeIntegrator :
	public ILinearTimeIntegrator<TDomain, TAlgebra>,
	public ITimeDiscDependentObject<TAlgebra>
{
protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;
public:
	typedef ILinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef ITimeDiscretization<TAlgebra> time_disc_type;
	typedef typename base_type::linear_solver_type linear_solver_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

private:

protected:

	int m_numSteps;
public:

	// constructor
	ConstStepLinearTimeIntegrator (SmartPtr<time_disc_type> tDisc)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc), m_numSteps(1) {}

	ConstStepLinearTimeIntegrator (SmartPtr<time_disc_type> tDisc, SmartPtr<typename base_type::linear_solver_type> lSolver)
	: base_type(lSolver), ITimeDiscDependentObject<TAlgebra>(tDisc), m_numSteps(1) {}

	void set_num_steps(int steps) {m_numSteps = steps;}

	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
};


/// Integrate over a given time interval (for a linear problem)
template<class TDomain, class TAlgebra>
class TimeIntegratorLinearAdaptive :
	public ILinearTimeIntegrator<TDomain, TAlgebra>,
	public ITimeDiscDependentObject<TAlgebra>
{
protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;
	
public:
	typedef ILinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef ITimeDiscretization<TAlgebra> time_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

protected:

	SmartPtr<time_disc_type> m_spTimeDisc2;        // set during init

	double m_tol, m_dtmin, m_dtmax;


public:
  TimeIntegratorLinearAdaptive (SmartPtr<time_disc_type> tDisc1, SmartPtr<time_disc_type> tDisc2)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc1), m_tol(1e-2), m_dtmin(-1.0), m_dtmax(-1.0)
	{
		m_spTimeDisc2 = tDisc2;
	}

	void init(grid_function_type const& u)
	{
		// call base
		base_type::init(u);
		//m_spTimeDisc2 = make_sp(new typename base_type::time_disc_type(base_type::m_spDomainDisc, base_type::m_theta));
	}
	
	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

    void set_tol(double tol) {m_tol = tol;}
    void set_time_step_min(number dt) {m_dtmin = dt;}
    void set_time_step_max(number dt) {m_dtmax = dt;}
};


class TimeStepBounds
{
public:
	TimeStepBounds()
	: m_dtMin(0.0), m_dtMax(std::numeric_limits<double>::max()),  m_redFac(1.0), m_incFac(1.0)
	{}

	void set_dt_min(double min) { m_dtMin = min; }
	double get_dt_min() { return m_dtMin; }

	void set_dt_max(double max) { m_dtMax = max; }
	double get_dt_max() { return m_dtMax; }

	void set_reduction_factor(double dec) { m_redFac = dec; }
	double get_reduction_factor() { return m_redFac; }

	void set_increase_factor(double inc) { m_incFac = inc; }
	double get_increase_factor() { return m_incFac; }

	void rescale(double alpha)
	{ m_dtMin*= alpha; m_dtMax*= alpha;}

protected:
	double m_dtMin, m_dtMax;
	double m_redFac, m_incFac;
};

/// Integration of non-linear systems (with bounds on dt)
template<class TDomain, class TAlgebra>
class INonlinearTimeIntegrator
: public ITimeIntegrator<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::vector_type vector_type;
	typedef IOperatorInverse<vector_type> solver_type;
	typedef AssembledOperator<TAlgebra> assembled_operator_type;

	INonlinearTimeIntegrator()
	: m_dtBounds() {}

	void set_solver(SmartPtr<solver_type> solver)
	{ m_spSolver=solver;}

	ConstSmartPtr<solver_type> get_solver() const
	{ return m_spSolver;}

	SmartPtr<solver_type> get_solver()
	{ return m_spSolver;}

	void set_dt_min(double min) { m_dtBounds.set_dt_min(min); }
	double get_dt_min() { return m_dtBounds. get_dt_min(); }

	void set_dt_max(double max) { m_dtBounds.set_dt_max(max); }
	double get_dt_max() { return m_dtBounds.get_dt_max(); }

	void set_reduction_factor(double dec) { m_dtBounds.set_reduction_factor(dec); }
	double get_reduction_factor() { return m_dtBounds.get_reduction_factor(); }

	void set_increase_factor(double inc) { m_dtBounds.set_increase_factor(inc); }
	double get_increase_factor() { return m_dtBounds.get_increase_factor(); }

protected:
	SmartPtr<solver_type> m_spSolver;
	TimeStepBounds m_dtBounds;

};


//! This class integrates (t0, t1] with stops at intermediate points tk.
template<class TDomain, class TAlgebra>
class DiscontinuityIntegrator :
		public INonlinearTimeIntegrator<TDomain, TAlgebra>
{
public:

	typedef INonlinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::grid_function_type grid_function_type;

	DiscontinuityIntegrator(SmartPtr<base_type> baseIntegrator) :
		base_type(), m_wrappedIntegrator(baseIntegrator), m_timePoints() {};

	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
	{
		int dstep = 0;
		auto tpoint = m_timePoints.begin();
		double tcurr = (tpoint == m_timePoints.end()) ? t1 : (*tpoint);
		double eps = 1e-8;

		// Perform first step.
		//this->notify_init_step(u0, dstep, t0,  tcurr-t0);
		bool status = m_wrappedIntegrator->apply(u1, tcurr*(1.0-eps), u0, t0);
		this->notify_finalize_step(u1, dstep++, tcurr, tcurr-t0);

		// Repeat for all intermediate points.
		while (tpoint != m_timePoints.end())
		{
			tpoint++;
			double tnext = (tpoint == m_timePoints.end()) ? t1 : (*tpoint);

			// Perform step.
			//this->notify_init_step(u1, dstep, tcurr, tnext-tcurr);
			status = status && m_wrappedIntegrator->apply(u1, tnext*(1.0-eps), u1, tcurr);
			this->notify_finalize_step(u1, dstep++, tnext, tnext-tcurr);

			tcurr = tnext;
		}
		this->notify_end(u1, dstep, t1, 0.0);
		return status;
	}

	void insert_points (std::vector<double> points)
	{
		m_timePoints = points;
	}
protected:
	SmartPtr<base_type> m_wrappedIntegrator;
	std::vector<double> m_timePoints;
};

} // namespace ug

#include "time_integrator_impl.hpp"

#endif /* __H__LIMEX__TIME_INTEGRATOR_HPP__ */
