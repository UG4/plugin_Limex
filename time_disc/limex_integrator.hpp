/*
 * Copyright (c) 2014 - 2020:  G-CSC, Goethe University Frankfurt
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

#ifndef LIMEX_INTEGRATOR_HPP_
#define LIMEX_INTEGRATOR_HPP_
/*
#define XMTHREAD_BOOST
#ifdef XMTHREAD_BOOST
#include <boost/thread/thread.hpp>
#endif
*/
#include <string>

#include "common/stopwatch.h"

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/operator/debug_writer.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/assemble_interface.h" // TODO: missing IAssemble in following file:
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/function_spaces/grid_function_util.h" // SaveVectorForConnectionViewer
#include "lib_disc/function_spaces/metric_spaces.h"
#include "lib_disc/io/vtkoutput.h"

#include "lib_grid/refinement/refiner_interface.h"


// own headers
#include "time_extrapolation.h"
#include "time_integrator.hpp"
#include "../limex_tools.h"
//#include "../multi_thread_tools.h"

#undef LIMEX_MULTI_THREAD

namespace ug {


/*
template<class TI>
class TimeIntegratorThread : public TI
{
	TimeIntegratorThread(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
	{

	}

	// thre
	void apply()
	{
		TI::apply(u1, t1, u0, t0);
	}


};*/

static void MyPrintError(UGError &err)
{
	for(size_t i=0;i<err.num_msg();++i)
	{
		UG_LOG("MYERROR "<<i<<":"<<err.get_msg(i)<<std::endl);
		UG_LOG("     [at "<<err.get_file(i)<<
		       ", line "<<err.get_line(i)<<"]\n");
	}
}

class ILimexRefiner
{
	virtual ~ILimexRefiner(){};

protected:
	SmartPtr<IRefiner> m_spRefiner;
};

//! Abstract class for the cost of a limex stage.
/*! The LimexTimeIntegrator requires a model for computing the cost per stage.  */
class ILimexCostStrategy
{
public:
	/// destructor
	virtual ~ILimexCostStrategy(){}

	/// provides the cost for all 0...nstages stages.
	virtual void update_cost(std::vector<number> &costA, const std::vector<size_t> &vSteps, const size_t nstages) = 0;
};


/// Cost is identical to (summation over) number of steps
class LimexDefaultCost : public ILimexCostStrategy
{
public:
	LimexDefaultCost(){};
	void update_cost(std::vector<number> &m_costA, const std::vector<size_t> &m_vSteps, const size_t nstages)
	{
		UG_ASSERT(m_costA.size() >= nstages, "Huhh: Vectors match in size:" << m_costA.size() << "vs." << nstages);

		UG_LOG("A_0="<< m_vSteps[0] << std::endl);
		m_costA[0] = (1.0)*m_vSteps[0];
		for (size_t i=1; i<nstages; ++i)
		{
			m_costA[i] = m_costA[i-1] + (1.0)*m_vSteps[i];
			UG_LOG("A_i="<< m_vSteps[i] << std::endl);
		}
	}
};

/// For
class LimexNonlinearCost : public ILimexCostStrategy
{
public:
	LimexNonlinearCost() :
	m_cAssemble (1.0), m_cMatAdd(0.0), m_cSolution(1.0), m_useGamma(1)
	{}

	void update_cost(std::vector<number> &m_costA, const std::vector<size_t> &nSteps, const size_t nstages)
	{
		UG_ASSERT(m_costA.size() >= nstages, "Huhh: Vectors match in size:" << m_costA.size() << "vs." << nstages);

		// 3 assemblies (M0, J0, Gamma0)
		m_costA[0] = (2.0+m_useGamma)*m_cAssemble + nSteps[0]*((1.0+m_useGamma)*m_cMatAdd + m_cSolution);
		for (size_t i=1; i<=nstages; ++i)
		{
			// n-1 assemblies for Mk, 2n MatAdds, n solutions
			m_costA[i] = m_costA[i-1] + (nSteps[i]-1) * m_cAssemble + nSteps[i]*((1.0+m_useGamma)*m_cMatAdd + m_cSolution);
		}
	}

protected:
	double m_cAssemble;  ///! Cost for matrix assembly
	double m_cMatAdd;	 ///! Cost for J=A+B
	double m_cSolution;  ///! Cost for solving Jx=b

	int m_useGamma;

};



//! Base class for LIMEX time integrator
template<class TDomain, class TAlgebra>
class LimexTimeIntegrator
: public INonlinearTimeIntegrator<TDomain, TAlgebra>,
  public VectorDebugWritingObject<typename TAlgebra::vector_type>
{


public:
		typedef TAlgebra algebra_type;
		typedef typename algebra_type::matrix_type matrix_type;
		typedef typename algebra_type::vector_type vector_type;
		typedef GridFunction<TDomain, TAlgebra> grid_function_type;
		typedef INonlinearTimeIntegrator<TDomain, TAlgebra> base_type;
		typedef typename base_type::solver_type solver_type;

		typedef IDomainDiscretization<algebra_type>	domain_discretization_type;
		typedef LinearImplicitEuler<algebra_type> timestep_type;
		typedef AitkenNevilleTimex<vector_type> timex_type;
		typedef INonlinearTimeIntegrator<TDomain, TAlgebra> itime_integrator_type;
		typedef SimpleTimeIntegrator<TDomain, TAlgebra> time_integrator_type;
		typedef ISubDiagErrorEst<vector_type> error_estim_type;

		//! Contains all data for parallel execution of time steps
		class ThreadData
		{
			//typedef boost::thread thread_type;
		public:

			ThreadData(SmartPtr<timestep_type> spTimeStep)
			: m_stepper(spTimeStep)
			{}

			ThreadData(){}
			SmartPtr<timestep_type> get_time_stepper()
			{ return m_stepper; }


			void set_solver(SmartPtr<solver_type> solver)
			{ m_solver = solver;}

			SmartPtr<solver_type> get_solver()
			{ return m_solver;}

			void set_error(int e)
			{ m_error=e; }


			void get_error()
			{ /*return m_error;*/ }


			void set_solution(SmartPtr<grid_function_type> sol)
			{ m_sol = sol;}

			SmartPtr<grid_function_type> get_solution()
			{ return m_sol;}

			void set_derivative(SmartPtr<grid_function_type> sol)
			{ m_dot = sol;}

			SmartPtr<grid_function_type> get_derivative()
			{ return m_dot;}

		protected:
			 // includes time step series
			SmartPtr<timestep_type> m_stepper;
			SmartPtr<solver_type> m_solver;

			SmartPtr<grid_function_type> m_sol;
			SmartPtr<grid_function_type> m_dot;
			int m_error;
		};

		typedef std::vector<SmartPtr<ThreadData> > thread_vector_type;

public:

	/// forward debug info to time integrators
	void set_debug_for_timestepper(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
	{
		for(size_t i=0; i<m_vThreadData.size(); ++i)
		{
				m_vThreadData[i].get_time_stepper()->set_debug(spDebugWriter);
		}
		UG_LOG("set_debug:" << m_vThreadData.size());
	}

	using VectorDebugWritingObject<vector_type>::set_debug;
protected:
	using VectorDebugWritingObject<vector_type>::write_debug;



public:
		LimexTimeIntegrator(int nstages)
		: m_tol(0.01),
		  m_rhoSafety(0.8),
		  m_sigmaReduction(0.5),
		  m_nstages(nstages),
		  m_gamma(m_nstages+1),
		  m_costA(m_nstages+1),
		  m_monitor(((m_nstages)*(m_nstages))), // TODO: wasting memory here!
		  m_workload(m_nstages),
		  m_lambda(m_nstages), 
		  m_num_reductions(m_nstages, 0),
		  m_max_reductions(2),
		  m_asymptotic_order(1000),
		  m_consistency_error(m_nstages),
		  m_greedyOrderIncrease(0.0),
		  m_useCachedMatrices(false),
		  m_spCostStrategy(make_sp<LimexDefaultCost>(new LimexDefaultCost())),
		  m_spBanachSpace(new IGridFunctionSpace<grid_function_type>()),              // default algebraic space
		  m_bInterrupt(false),
		  m_limex_step(1)
		{
			m_vThreadData.reserve(m_nstages);
			m_vSteps.reserve(m_nstages);

			// init exponents (i.e. k+1, k, 2k+1, ...)
			init_gamma();
		}


		/// tolerance
		void set_tolerance(double tol) { m_tol = tol;}
		void set_stepsize_safety_factor(double rho) { m_rhoSafety = rho;}
		void set_stepsize_reduction_factor(double sigma) { m_sigmaReduction = sigma;}
		void set_stepsize_greedy_order_factor(double sigma) { m_greedyOrderIncrease = sigma;}

		void set_max_reductions(size_t nred) { m_max_reductions = nred;}
		void set_asymptotic_order(size_t q) { m_asymptotic_order = q;}

		/// add an error estimator
		void add_error_estimator(SmartPtr<error_estim_type> spErrorEstim)
		{ m_spErrorEstimator = spErrorEstim; }

		//! add a new stage (at end of list)
		void add_stage_base(size_t nsteps, SmartPtr<solver_type> solver, SmartPtr<domain_discretization_type> spDD, SmartPtr<domain_discretization_type> spGamma=SPNULL)
		{
			UG_ASSERT(m_vThreadData.size() == m_vSteps.size(), "ERROR: m_vThreadData and m_vSteps differ in size!");

			UG_ASSERT(m_vThreadData.empty() || m_vSteps.back()<nsteps, "ERROR: Sequence of steps must be increasing." );

			// all entries have been set
			if (m_vThreadData.size() == m_nstages)
			{ return; }


			// a) set number of steps
			m_vSteps.push_back(nsteps);

			// b) set time-stepper
			SmartPtr<timestep_type> limexStepSingleton;
#ifndef LIMEX_MULTI_THREAD

			if(m_vThreadData.size()>0) {
				// re-use time-stepper (if applicable)
				limexStepSingleton = m_vThreadData.back().get_time_stepper();
			}
			else
			{
				// create time-stepper
#endif
				// for mult-threading, each info object has own time-stepper
				if (spGamma.invalid()) {
					limexStepSingleton = make_sp(new timestep_type(spDD));
				} else {
					limexStepSingleton = make_sp(new timestep_type(spDD, spDD, spGamma));
				}
				UG_ASSERT(limexStepSingleton.valid(), "Huhh: Invalid pointer")
#ifndef LIMEX_MULTI_THREAD
			}
#endif
			// propagate debug info
			limexStepSingleton->set_debug(VectorDebugWritingObject<vector_type>::vector_debug_writer());
			m_vThreadData.push_back(ThreadData(limexStepSingleton));

			// c) set solver
			UG_ASSERT(solver.valid(), "Huhh: Need to supply solver!");
			m_vThreadData.back().set_solver(solver);
			UG_ASSERT(m_vThreadData.back().get_solver().valid(), "Huhh: Need to supply solver!");
		}

		void add_stage(size_t nsteps, SmartPtr<solver_type> solver, SmartPtr<domain_discretization_type> spDD)
		{ add_stage_base(nsteps, solver, spDD); }

		void add_stage_ext(size_t nsteps, SmartPtr<solver_type> solver, SmartPtr<domain_discretization_type> spDD, SmartPtr<domain_discretization_type> spGamma)
		{ add_stage_base(nsteps, solver, spDD, spGamma); }

		///! TODO: remove this function!
		void add_stage(size_t i, size_t nsteps, SmartPtr<domain_discretization_type> spDD, SmartPtr<solver_type> solver)
		{
			UG_LOG("WARNING: add_stage(i, nsteps ,...) is deprecated. Please use 'add_stage(nsteps ,...) instead!'");
			add_stage(nsteps, solver, spDD);
		}

		//!
		SmartPtr<timestep_type> get_time_stepper(size_t i)
		{
			UG_ASSERT(i<m_vThreadData.size(), "Huhh: Invalid entry");
			return m_vThreadData[i].get_time_stepper();
		}




protected:
		//! Initialize integrator threads (w/ solutions)
		void init_integrator_threads(ConstSmartPtr<grid_function_type> u);

		//! (Tentatively) apply integrators
		int apply_integrator_threads(number dtcurr, ConstSmartPtr<grid_function_type> u0, number t0, size_t nstages);

		//! e.g. wait for all threads to complete
		void join_integrator_threads();

		//! Override thread-wise solutions with common solution
		void update_integrator_threads(ConstSmartPtr<grid_function_type> ucommon, number t);

		//! Dispose integrator threads (w/ solutions)
		void dispose_integrator_threads();

public:
		//! Integrating from t0 -> t1
		bool apply(SmartPtr<grid_function_type> u, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

		number 	get_cost(size_t i) { return m_costA[i]; }
		number 	get_gamma(size_t i) { return m_gamma[i]; }
		number 	get_workload(size_t i) { return m_workload[i]; }


protected:
		number& monitor(size_t k, size_t q) {
		  UG_ASSERT(k<m_nstages, "Huhh: k mismatch");
		  UG_ASSERT(q<m_nstages, "Huhh: q mismatch");
		  return m_monitor[k+(m_nstages)*q]; 
		}

		/// aux: compute exponents gamma_k (for roots)
		void init_gamma()
		{
			for (size_t k=0; k<=m_nstages; ++k)
			{
				m_gamma[k] = 2.0+k;
			}
		}

		/// Updating workloads A_i for computing T_ii
		//  (depends on m_vSteps, which must have been initialized!)
		void update_cost()
		{
			m_spCostStrategy->update_cost(m_costA, m_vSteps, m_nstages);
		}

		/// convergence monitor
		// (depends on cost, which must have been initialized!)
		void update_monitor()
		{
		  UG_ASSERT(m_costA.size()>= m_nstages, "Cost vector too small!")
		  for (size_t k=0; k<m_nstages-1; ++k)
		  {
		    UG_LOG("k= "<< k<< ", A[k]=" << m_costA[k] << ", gamma[k]=" << m_gamma[k] << "\t");
		    for (size_t q=0; q<m_nstages-1; ++q)
		      {
			// Deuflhard: Order and stepsize, ... eq. (3.7)
			double gamma = (m_costA[k+1] - m_costA[0] + 1.0)/(m_costA[q+1] - m_costA[0] + 1.0);
			double alpha = pow(m_tol, gamma);
			
			// for fixed order q, the monitor indicates the performance penalty compared to a strategy using k stages only
			// cf. eq. (4.6)
			monitor(k,q) = pow(alpha/(m_tol * m_rhoSafety), 1.0/m_gamma[k]);
			UG_LOG(monitor(k,q) << "[" << pow(alpha/(m_tol), 1.0/m_gamma[k]) <<"," << gamma<< ","<< m_costA[k+1] << "," << m_costA[q+1]<< ","<< alpha << "]" << "\t");
			// UG_LOG(  << "\t");
			
		      }
		    UG_LOG(std::endl);
		  }

		}

		// Find row k=1, ..., ntest-1 minimizing (estimated) error eps[kmin]
		// Also: predict column q with minimal workload W_{k+1,k} = A_{k+1} * lambda_{k+1}
		size_t find_optimal_solution(const std::vector<number>& eps, size_t ntest, /*size_t &kf,*/ size_t &qpred);

public:
		/// setter for time derivative info (optional for setting $\Gamma$)
		void set_time_derivative(SmartPtr<grid_function_type> udot) {m_spDtSol = udot;}

		/// getter for time derivative info (optional for setting $\Gamma$)
		SmartPtr<grid_function_type> get_time_derivative() {return m_spDtSol;}

		/// status for time derivative info (optional for setting $\Gamma$)
		bool has_time_derivative() {return m_spDtSol!=SPNULL;}


		void enable_matrix_cache() { m_useCachedMatrices = true; }  	///< Select classic LIMEX
		void disable_matrix_cache() { m_useCachedMatrices = false; }    ///< Select approximate Newton (default)

		void select_cost_strategy(SmartPtr<ILimexCostStrategy> cost)
		{ m_spCostStrategy = cost;}


		/// set banach space (e.g. for computing consistency error)
		void set_space(SmartPtr<IGridFunctionSpace<grid_function_type> > spSpace)
		{ 
		  m_spBanachSpace = spSpace;
		  UG_LOG("set_space:" << m_spBanachSpace->config_string()); 
		}

		/// interrupt execution of apply() by external call via observer
		void interrupt() {m_bInterrupt = true;}

protected:

		double m_tol;
		double m_rhoSafety;
		double m_sigmaReduction;
		SmartPtr<error_estim_type> m_spErrorEstimator;     // (smart ptr for) error estimator

		unsigned int m_nstages; 				///< Number of Aitken-Neville stages
		std::vector<size_t> m_vSteps;			///< generating sequence for extrapolation
		std::vector<ThreadData> m_vThreadData;	///< vector with thread information

		std::vector<number> m_gamma;			///< gamma_i: exponent

		std::vector<number> m_costA;			///< Cost A_i (for completing stage i)
		std::vector<number> m_monitor;			///< Convergence monitor \alpha

		std::vector<number> m_workload;
		std::vector<number> m_lambda;

		std::vector<size_t> m_num_reductions;		///< history of reductions
		size_t m_max_reductions;
		size_t m_asymptotic_order;  				///< For PDEs, we may apply an symptotic order reduction

		std::vector<number> m_consistency_error;///<Consistency error
  		double m_greedyOrderIncrease;

		SmartPtr<grid_function_type> m_spDtSol;   ///< Time derivative
		bool m_useCachedMatrices;

		SmartPtr<ILimexCostStrategy> m_spCostStrategy;

		/// metric space
		SmartPtr<IGridFunctionSpace<grid_function_type> > m_spBanachSpace;

		bool m_bInterrupt;

		int m_limex_step;
};


/*! Create private solutions for each thread */
template<class TDomain, class TAlgebra>
void LimexTimeIntegrator<TDomain,TAlgebra>::init_integrator_threads(ConstSmartPtr<grid_function_type> u)
{
	PROFILE_FUNC_GROUP("limex");
	const int nstages = m_vThreadData.size()-1;
	for (int i=nstages; i>=0; --i)
	{
		m_vThreadData[i].set_solution(u->clone());
		m_vThreadData[i].set_derivative(u->clone());
		m_vThreadData[i].get_time_stepper()->set_matrix_cache(m_useCachedMatrices);
	}
}

/*! Create private solutions for each thread */
template<class TDomain, class TAlgebra>
void LimexTimeIntegrator<TDomain,TAlgebra>::dispose_integrator_threads()
{
	PROFILE_FUNC_GROUP("limex");
	const int nstages = m_vThreadData.size()-1;
	for (int i=nstages; i>=0; --i)
	{
		m_vThreadData[i].set_solution(SPNULL);
		m_vThreadData[i].set_derivative(SPNULL);
		// m_vThreadData[i].get_time_stepper()->set_matrix_cache(m_useCachedMatrices);
	}
}


// create (& execute) threads
/*boost::thread_group g;
		typename thread_vector_type::reverse_iterator rit=m_vThreadData.rbegin();
		for (rit++; rit!= m_vThreadData.rend(); ++rit)
		{

			boost::thread *t =new boost::thread(boost::bind(&ThreadSafeTimeIntegrator::apply, *rit));
			//g.add_thread(t);

			g.create_thread(boost::bind(&ThreadSafeTimeIntegrator::apply, *rit));

		}*/


/* TODO: PARALLEL execution?*/
template<class TDomain, class TAlgebra>
int LimexTimeIntegrator<TDomain,TAlgebra>::apply_integrator_threads(number dtcurr, ConstSmartPtr<grid_function_type> u0, number t0, size_t nstages)
{
	PROFILE_FUNC_GROUP("limex");
	update_cost();		// compute cost A_i (alternative: measure times?)
	update_monitor();	// convergence monitor

	/*
			int tn = omp_get_thread_num();
			int nt = omp_get_num_threads();
			omp_set_num_threads(nstages);
	 */
	int error = 0;
	//const int nstages = m_vThreadData.size()-1;
	//	#pragma omp for private(i) // shared (nstages, u1) schedule(static)
	for (int i=0; i<=(int)nstages; ++i)
	{
		/*
				std::cerr << "I am " << tn << " of " << nt << " ("<< i<< "/" << nstages<<")!" << std::endl;
				UGMultiThreadEnvironment mt_env;
		 */

		// copy data to private structure (one-to-many)
		//m_vThreadData[i].set_solution(u1->clone());

		// switch to "child" comm
		// mt_env.begin();

		// integrate (t0, t0+dtcurr)
		time_integrator_type integrator(m_vThreadData[i].get_time_stepper());
		integrator.set_time_step(dtcurr/m_vSteps[i]);
		integrator.set_dt_min(dtcurr/m_vSteps[i]);
		integrator.set_dt_max(dtcurr/m_vSteps[i]);
		integrator.set_reduction_factor(0.0);                 // quit immediately, if step fails
		integrator.set_solver(m_vThreadData[i].get_solver());
		integrator.set_derivative(m_vThreadData[i].get_derivative());

		UG_ASSERT(m_spBanachSpace.valid(), "Huhh: Need valid (default) banach space");
		integrator.set_banach_space(m_spBanachSpace);
		UG_LOG("Set space:" << m_spBanachSpace->config_string());

		bool exec = true;
		try
		{
			exec = integrator.apply(m_vThreadData[i].get_solution(), t0+dtcurr, u0, t0);
			m_consistency_error[i] = integrator.get_consistency_error();
		}
		catch(ug::UGError& err)
		{

			exec = false;
			error += (1 << i);
			UG_LOGN("Step "<< i<< " failed on stage " << i << ": " << err.get_msg());
			MyPrintError(err);

		}

		if (!exec)
		{
			// Additional actions at failure
		}

		// switch to "parent" comm
		//mt_env.end();
	} /*for all stages loop*/



	return error;
}


template<class TDomain, class TAlgebra>
void LimexTimeIntegrator<TDomain,TAlgebra>::join_integrator_threads()
{
	// join all threads
	// g.join_all();
	const int nstages = m_vThreadData.size()-1;
	for (int i=nstages; i>=0; --i)
	{
		m_vThreadData[i].get_time_stepper()->invalidate();
	}
}


template<class TDomain, class TAlgebra>
void LimexTimeIntegrator<TDomain,TAlgebra>::update_integrator_threads(ConstSmartPtr<grid_function_type> ucommon, number t)
{
	const int nstages = m_vThreadData.size()-1;
	for (int i=nstages; i>=0; --i)
	{
		UG_ASSERT(m_vThreadData[i].get_solution()->size()==ucommon->size(), "LIMEX: Vectors must match in size!")
		*m_vThreadData[i].get_solution() = *ucommon;
	}
}

template<class TDomain, class TAlgebra>
size_t LimexTimeIntegrator<TDomain,TAlgebra>::
find_optimal_solution(const std::vector<number>& eps, size_t ntest, /*size_t &kf,*/ size_t &qpred)
{

	const size_t qold=qpred;

	size_t jbest = 1;
	qpred = 1;

	size_t j=1;
	size_t k=j-1;

	m_lambda[k] = pow(m_rhoSafety*m_tol/eps[j], 1.0/m_gamma[k]);   // 1/epsilon(k)
	m_workload[k] = m_costA[j]/m_lambda[k];
	UG_LOG("j=" << j << ": eps=" << eps[j]  << ", lambda(j)=" <<m_lambda[k]  << ", epsilon(j)=" <<1.0/m_lambda[k] << "<= alpha(k, qcurr)=" << monitor(k,qold-1) << "< alpha(k, qcurr+1)=" << monitor(k,qold) <<", A="<< m_costA[j] << ", W="<< m_workload[k] <<std::endl);

	for (j=2; j<ntest; ++j)
	{
	  k = j-1;
	  m_lambda[k] = pow(m_rhoSafety*m_tol/eps[j], 1.0/m_gamma[k]);
	  m_workload[k] = m_costA[j]/m_lambda[k];
	  UG_LOG("j=" << j << ": eps=" << eps[j]  << ", lambda(j)=" <<m_lambda[k]  << ", epsilon(j)=" <<1.0/m_lambda[k] << "<= alpha(k, qcurr)=" << monitor(k,qold-1) << "< alpha(k, qcurr+1)=" << monitor(k,qold) <<", A="<< m_costA[j] << ", W="<< m_workload[k] <<std::endl);


	  // TODO: Convergence monitor

	  qpred = (m_workload[qpred-1] > m_workload[k]) ? j : qpred;
	  jbest = (eps[jbest] > eps [j]) ? j : jbest;
	}

	return jbest;
}

template<class TDomain, class TAlgebra>
bool LimexTimeIntegrator<TDomain,TAlgebra>::
apply(SmartPtr<grid_function_type> u, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	PROFILE_FUNC_GROUP("limex");
#ifdef UG_OPENMP
	// create multi-threading environment
	//int nt = std::min(omp_get_max_threads(), m_nstages);

#endif

	// NOTE: we use u as common storage for future (and intermediate) solution(s)
	if (u.get() != u0.get())    // only copy if not already identical, otherwise: PST_UNDEFINED!
		*u = *u0;

	// initialize integrator threads
	// (w/ solutions)
	init_integrator_threads(u0);


	// write_debug
	for (unsigned int i=0; i<m_vThreadData.size(); ++i)
	{
		std::ostringstream ossName;
		ossName << std::setfill('0') << std::setw(4);
		ossName << "Limex_Init_iter" << 0 << "_stage" << i;
		write_debug(*m_vThreadData[i].get_solution(), ossName.str().c_str());
	}

	number t = t0;
	double dtcurr = ITimeIntegrator<TDomain, TAlgebra>::get_time_step();

	size_t kmax = m_vThreadData.size();   		// maximum number of stages
	size_t qpred = kmax-1;  	 						// predicted optimal order
	size_t qcurr = qpred;								// current order

	// for Gustafsson/lundh/Soederlind type PID controller
	/*size_t qlast = 0;
	double dtlast = 0.0;
	std::vector<double> epslast(kmax, m_rhoSafety*m_tol);
*/

	// time integration loop
	SmartPtr<grid_function_type> ubest = SPNULL;
	size_t limex_total = 1;
	size_t limex_success = 0;
	size_t ntest;    ///< active number of stages <= kmax
	size_t jbest;

	m_bInterrupt = false;
	//bool bProbation = false;
	bool bAsymptoticReduction = false;

	const size_t nSwitchHistory=16;
	const size_t nSwitchLookBack=5;
	int vSwitchHistory[nSwitchHistory] ={ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	timex_type timex(m_vSteps);
	while ((t < t1) && ((t1-t) > base_type::m_precisionBound))
	{
		int err = 0;

		//UG_DLOG(LIB_LIMEX, 5, "+++ LimexTimestep +++" << limex_step << "\n");
		UG_LOG("+++ LimexTimestep +++" << m_limex_step << "\n");


		// save time stamp for limex step start
		Stopwatch stopwatch;
		stopwatch.start();


		// determine step size
		number dt = std::min(dtcurr, t1-t);
		UG_COND_THROW(dt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!");


		// Notify init observers. (NOTE: u = u(t))¸
		itime_integrator_type::notify_init_step(u, m_limex_step, t, dt);


		// determine number of stages to investigate
		qcurr = qpred;
		ntest = std::min(kmax, qcurr+1);
		UG_LOG("ntest="<< ntest << std::endl);

		// checks
		UG_ASSERT(m_vSteps.size() >= ntest, "Huhh: sizes do not match: " << m_vSteps.size() << "<"<<ntest);
		UG_ASSERT(m_vThreadData.size() >= ntest, "Huhh: sizes do not match: " << m_vThreadData.size() << "< "<< ntest);

		///////////////////////////////////////
		// PARALLEL EXECUTION: BEGIN

		// write_debug
		for (size_t i=0; i<ntest; ++i)
		{
			std::ostringstream ossName;
			ossName << std::setfill('0') << std::setw(4);
			ossName << "Limex_BeforeSerial_iter" << m_limex_step << "_stage" << i << "_total" << limex_total;
			write_debug(*m_vThreadData[i].get_solution(), ossName.str().c_str());
		}

		// integrate: t -> t+dt
		err = apply_integrator_threads(dt, u, t, ntest-1);

		// write_debug
		for (size_t i=0; i<ntest; ++i)
		{
			std::ostringstream ossName;
			ossName << std::setfill('0') << std::setw(4);
			ossName << "Limex_AfterSerial_iter" << m_limex_step << "_stage" << i << "_total" << limex_total;
			write_debug(*m_vThreadData[i].get_solution(), ossName.str().c_str());
		}

		join_integrator_threads();
		// PARALLEL EXECUTION: END
		///////////////////////////////////////


		///////////////////////////////////////
		// SERIAL EXECUTION: BEGIN

		// sanity checks
		UG_ASSERT(m_spErrorEstimator.valid(), "Huhh: Invalid Error estimator?");

		double epsmin = 0.0;

		bool limexConverged = false;
		if (err==0)
		{	// Compute extrapolation at t+dtcurr (SERIAL)
			// Consistency check.
			for (unsigned int i=1; i<ntest; ++i)
			{
				UG_LOG("Checking consistency:" <<m_consistency_error[i] <<"/" << m_consistency_error[i-1] << "="<< m_consistency_error[i]/m_consistency_error[i-1] << std::endl);
			}

			// Extrapolation.
			timex.set_error_estimate(m_spErrorEstimator);
			for (unsigned int i=0; i<ntest; ++i)
			{
				UG_ASSERT(m_vThreadData[i].get_solution().valid(), "Huhh: no valid solution?");
				timex.set_solution(m_vThreadData[i].get_solution(), i);
			}
			timex.apply(ntest);

			// write_debug
			for (size_t i=0; i<ntest; ++i)
			{
				std::ostringstream ossName;
				ossName << std::setfill('0') << std::setw(4);
				ossName << "Limex_Extrapolates_iter" << m_limex_step << "_stage" << i << "_total" << limex_total;
				write_debug(*m_vThreadData[i].get_solution(), ossName.str().c_str());
			}
			limex_total++;

			// Obtain sub-diagonal error estimates.
			const std::vector<number>& eps = timex.get_error_estimates();
			UG_ASSERT(ntest<=eps.size(), "Huhh: Not enough solutions?");

			// select optimal solution (w.r.t error) AND
			// predict optimal order (w.r.t. workload) for next step
			jbest = find_optimal_solution(eps, ntest, qpred);
			UG_ASSERT(jbest < ntest, "Huhh: Not enough solutions?");

			// best solution
			ubest  = timex.get_solution(jbest).template cast_dynamic<grid_function_type>();
			epsmin = eps[jbest];

			// check for convergence
			limexConverged = (epsmin <= m_tol) ;


			if (limex_success>3)
			{

				vSwitchHistory[m_limex_step%nSwitchHistory] = (qpred - qcurr);
				UG_DLOG(LIB_LIMEX, 5, "LIMEX-ASYMPTOTIC-ORDER switch:  = " <<  (qpred - qcurr)<< std::endl);

				size_t nSwitches=0;
				for (int s=nSwitchLookBack-1; s>=0; s--)
				{
					nSwitches += std::abs(vSwitchHistory[(m_limex_step-s)%nSwitchHistory]);
					UG_DLOG(LIB_LIMEX, 6, "LIMEX-ASYMPTOTIC-ORDER: s[" << s<< "] = " <<  vSwitchHistory[(m_limex_step-s)%nSwitchHistory] << std::endl);
				}
				UG_DLOG(LIB_LIMEX, 5,"LIMEX-ASYMPTOTIC-ORDER: nSwitches = " <<  nSwitches << std::endl);


			/*
			if (bProbation)
			{

				if ((qpred < qcurr)	|| (!limexConverged))	// consecutive (follow-up) order decrease
				{
					bProbation = true;
					m_num_reductions[0]++;
					UG_LOG("LIMEX-ASYMPTOTIC-ORDER: Decrease on parole detected: "<< qcurr << " -> " << qpred << "("<< m_num_reductions[0]<<")"<<std::endl);

				}
				else
				{
					// reset all counters
					bProbation = false;
					UG_LOG("LIMEX-ASYMPTOTIC-ORDER: Probation off!"<<std::endl);

				}
			}
			else
			{
				if ((qpred < qcurr)	|| (!limexConverged))// 1st order decrease
				{
					bProbation = true;
					m_num_reductions[0]++;
					UG_LOG("LIMEX-ASYMPTOTIC-ORDER:  Decrease detected: "<< qcurr << " -> " << qpred << "("<< m_num_reductions[0]<<")"<<std::endl);
				}
				else
				{
					m_num_reductions[0] = 0;
					UG_LOG("LIMEX-ASYMPTOTIC-ORDER: I am totally free!"<<std::endl);
				}


			}
*/

			// bAsymptoticReduction = (m_num_reductions[0] >= m_max_reductions) || bAsymptoticReduction;

			// bAsymptoticReduction = (nSwitches >= m_max_reductions) || bAsymptoticReduction;

			if (nSwitches >= m_max_reductions) {


				m_asymptotic_order = std::min<size_t>(m_asymptotic_order-1, 2);
				bAsymptoticReduction = true;

				for (int s=nSwitchLookBack-1; s>=0; s--)
				{ vSwitchHistory[(m_limex_step-s)%nSwitchHistory]=0; }


			}



			// asymptotic order reduction
			//if (m_num_reductions[0] >= m_max_reductions)
			if	(bAsymptoticReduction)
			{
				kmax = (kmax >= m_asymptotic_order ) ? m_asymptotic_order : kmax;
				qpred = kmax - 1;
				UG_DLOG(LIB_LIMEX, 5, "LIMEX-ASYMPTOTIC-ORDER: Reduction: "<< qpred );
				UG_DLOG(LIB_LIMEX, 5, "(kmax=" << kmax << ", asymptotic"<<  m_asymptotic_order <<") after "<<m_max_reductions<< std::endl);
			}

			} // (limex_success>3)

			/*
			// adjust for order bound history
			if ((m_num_reductions[qpred] > m_max_reductions) && (qpred > 1))
			{
					UG_LOG("prohibited  q:"<< qpred << "(" << m_num_reductions[qpred] << ")"<<std::endl);
					//qpred--;
					qpred = kmax = 2;
			} else
			{
				UG_LOG("keeping  q:"<< qpred << "(" << m_num_reductions[qpred] << ")"<<std::endl);
			}
			*/

			//double pid = 1.0;

			// constant order?
			/*if (qcurr == qlast)
			{
				double dtRatio =dtcurr/dtlast;
				// Dd, Band 3 (DE), p. 360
				int k = ntest-1;
				while (k--) {
					double epsRatio =eps[k+1]/epslast[k+1];

					double qopt = log(epsRatio) / log(dtRatio);  // take dtcurr here, as dt may be smaller

					UG_LOG("effective p["<< k << "]="<< qopt-1.0 <<","<< eps[k+1] <<","<< epslast[k+1] <<","  <<  epsRatio << "("<<  log(epsRatio)  << ","<<  pow(epsRatio, 1.0/m_gamma[k])  <<")"<< dtcurr <<","<< dtlast   << "," << dtRatio << "," <<log(dtRatio) << std::endl);
				}


				double epsRatio = eps[qcurr]/epslast[qcurr];
				pid = dtRatio/pow(epsRatio, 1.0/m_gamma[qcurr-1]);

				UG_LOG("pid=" << pid << std::endl);

			}*/


			// select (predicted) order for next step
			double dtpred = dtcurr*std::min(m_lambda[qpred-1], itime_integrator_type::get_increase_factor());
			//double dtpred = dtcurr*m_lambda[qpred-1];
			UG_LOG("+++++\nget_increase_factor() gives "<<itime_integrator_type::get_increase_factor()<<" \n+++++++")
			UG_LOG("koptim=\t" << jbest << ",\t eps(k)=" << epsmin << ",\t q=\t" << qpred<< "("<<  ntest << "), lambda(q)=" << m_lambda[qpred-1] << ", alpha(q-1,q)=" << monitor(qpred-1, qpred) << "dt(q)=" << dtpred<< std::endl);

			// EXTENSIONS: convergence model
			if (limexConverged)
			{
				limex_success++;

				// a) aim for order increase in next step
			  if ((qpred+1==ntest)  /* increase by one possible? */
			      //    && (m_lambda[qpred-1]>1.0)
				  // && (m_num_reductions[qpred+1] <= m_max_reductions)
			      && (kmax>ntest)) /* still below max? */
				{
					const double alpha = monitor(qpred-1, qpred);
					UG_LOG("CHECKING for order increase: "<< m_costA[qpred] << "*" << alpha << ">" << m_costA[qpred+1]);
					// check, whether further increase could still be efficient
					if (m_costA[qpred] * alpha > m_costA[qpred+1])
					{
						qpred++;    			// go for higher order
						if (m_greedyOrderIncrease >0.0) {
						  dtpred *= (1.0-m_greedyOrderIncrease) + m_greedyOrderIncrease*alpha;		// & adapt time step  // TODO: check required!
						}
						UG_LOG("... yes.\n")

						// update history
						vSwitchHistory[m_limex_step%nSwitchHistory] = (qpred - qcurr);
						UG_LOG("LIMEX-ASYMPTOTIC-ORDER switch update:  = " <<  (qpred - qcurr)<< std::endl);

					} else {
						UG_LOG("... nope.\n")
					}

				}


				// b) monitor convergence (a-priori check!)

			}

			// parameters for subsequent step
			dtcurr = std::min(dtpred, itime_integrator_type::get_dt_max());


		}
		else
		{
			// solver failed -> cut time step
			dtcurr *= m_sigmaReduction;
		}

		// output compute time for Limex step
		number watchTime = stopwatch.ms()/1000.0;
		UG_LOGN("Time: " << std::setprecision(3) << watchTime << "s");

		if ((err==0) && limexConverged)
		{
			// ACCEPT time step
			UG_LOG("+++ LimexTimestep +++" << m_limex_step << " ACCEPTED"<< std::endl);
			UG_LOG("               :\t time \t dt (success) \t dt (pred) \tq=\t order (curr)" << qcurr+1 << std::endl);
			UG_LOG("LIMEX-ACCEPTING:\t" << t <<"\t"<< dt << "\t" << dtcurr << "\tq=\t" << qcurr+1 << std::endl);


			// update PID controller
			/*qlast = qcurr;
			epslast =  timex.get_error_estimates();  // double check ???
			dtlast = dt;
*/
			// Compute time derivatives (by extrapolation)
			if (this->has_time_derivative())
			{
				UG_LOG("Computing derivative" << std::endl);
				grid_function_type &udot = *get_time_derivative();

				for (size_t i=0; i<=jbest; ++i)
				{
					timex.set_solution(m_vThreadData[i].get_derivative(), i);
				}
				timex.apply(jbest+1, false); // do not compute error

				udot = *timex.get_solution(jbest).template cast_dynamic<grid_function_type>();

				std::ostringstream ossName;
				ossName << std::setfill('0');
				ossName << "Limex_Derivative_iter" << m_limex_step << "_total" << limex_total;
				write_debug(udot, ossName.str().c_str());
			}


			// Copy best solution.
			UG_ASSERT(ubest.valid(), "Huhh: Invalid error estimate?");
			*u = *ubest;
			t += dt;

			// Initiate post process.
			itime_integrator_type::notify_finalize_step(u, m_limex_step++, t, dt);

			// make sure that all threads continue
			// with identical initial value u(t)
			// update_integrator_threads(ubest, t);


			// working on last row => increase order
			//if (ntest == q+1) ntest++;
		}
		else
		{
			// DISCARD time step
			UG_LOG("+++ LimexTimestep +++" << m_limex_step << " FAILED" << std::endl);
			UG_LOG("LIMEX-REJECTING:\t" << t <<"\t"<< dt << "\t" << dtcurr << std::endl);

			itime_integrator_type::notify_rewind_step(ubest, m_limex_step, t+dt, dt);

		}

		// interrupt if interruption flag set
		if (m_bInterrupt)
		{
			UG_LOGN("Limex interrupted by external command.");
			break;
		}

		// SERIAL EXECUTION: END
		///////////////////////////////////////
		update_integrator_threads(u, t);


		// SOLVE

		// ESTIMATE

		// MARK

		// REFINE


	} // time integration loop

	dispose_integrator_threads();

	return true;
} // apply


} // namespace ug

#endif /* TIME_INTEGRATOR_HPP_ */
