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

#ifndef LIMEX_INTEGRATOR_HPP_
#define LIMEX_INTEGRATOR_HPP_
/*
#define XMTHREAD_BOOST
#ifdef XMTHREAD_BOOST
#include <boost/thread/thread.hpp>
#endif
*/
#include <string>

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
#include "lib_disc/io/vtkoutput.h"

#include "lib_grid/refinement/refiner_interface.h"


// own headers
#include "time_extrapolation.h"
#include "time_integrator.hpp"
#include "../limex_tools.h"
//#include "../multi_thread_tools.h"


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

			SmartPtr<timestep_type> get_time_stepper()
			{ return m_stepper; }

			void set_solution(SmartPtr<grid_function_type> sol)
			{ m_sol = sol;}

			SmartPtr<grid_function_type> get_solution()
			{ return m_sol;}

			void set_solver(SmartPtr<solver_type> solver)
			{ m_solver = solver;}

			SmartPtr<solver_type> get_solver()
			{ return m_solver;}

			void set_error(int e)
			{ m_error=e; }

			void get_error()
			{ return m_error; }


		protected:
			 // includes time step series
			SmartPtr<timestep_type> m_stepper;
			SmartPtr<grid_function_type> m_sol;
			SmartPtr<solver_type> m_solver;
			int m_error;

		};

		typedef std::vector<SmartPtr<ThreadData> > thread_vector_type;

public:
	using VectorDebugWritingObject<vector_type>::set_debug;
protected:
	using VectorDebugWritingObject<vector_type>::write_debug;



public:
		LimexTimeIntegrator(int nstages)
		: m_tol(0.01),
		  m_rhoSafety(0.8),
		  m_sigmaReduction(0.5),
		  m_nstages(nstages),
		  m_costA(m_nstages+1),
		  m_gamma(m_nstages+1),
		  m_monitor(((m_nstages+1)*(m_nstages+1))), // TODO: wasting memory here!
		  m_workload(m_nstages+1),
		  m_lambda(m_nstages+1), 
		  m_greedyOrderIncrease(0.0)
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


		void add_error_estimator(SmartPtr<error_estim_type> spErrorEstim)
		{ m_spErrorEstimator = spErrorEstim; }

		//! add a new stage (at end of list)
		void add_stage(size_t nsteps, SmartPtr<domain_discretization_type> spDD, SmartPtr<solver_type> solver)
		{
			UG_ASSERT(m_vThreadData.size() == m_vSteps.size(), "ERROR: m_vThreadData and m_vSteps differ in size!");

			UG_ASSERT(m_vThreadData.empty() || m_vSteps.back()<nsteps, "ERROR: Sequence of steps must be increasing." );

			if (m_vThreadData.size() == m_nstages)
			{
				return;
			}

			m_vSteps.push_back(nsteps);

			m_vThreadData.push_back(ThreadData(make_sp(new timestep_type(spDD))));
			m_vThreadData.back().set_solver(solver);

			UG_ASSERT(solver.valid(), "Huhh: Need to supply solver!");
			UG_ASSERT(m_vThreadData.back().get_solver().valid(), "Huhh: Need to supply solver!");
		}

		///! TODO: remove this function!
		void add_stage(size_t i, size_t nsteps, SmartPtr<domain_discretization_type> spDD, SmartPtr<solver_type> solver)
		{
			UG_LOG("WARNING: add_stage(i, nsteps ,...) is deprecated. Please use 'add_stage(nsteps ,...) instead!'");
			add_stage(nsteps, spDD, solver);
		}



protected:

		//! Initialize integrator threads (w/ solutions)
		/*! Create private solutions for each thread */
		void init_integrator_threads(ConstSmartPtr<grid_function_type> u)
		{
			const int nstages = m_vThreadData.size()-1;
			for (int i=nstages; i>=0; --i)
			{
				m_vThreadData[i].set_solution(u->clone());
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

		//! (Tentatively) apply integrators
		// TODO: PARALLEL execution?
		int apply_integrator_threads(number dtcurr, ConstSmartPtr<grid_function_type> u0, number t0, size_t nstages)
		{

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
			for (int i=nstages; i>=0; --i)
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

				bool exec = true;
				try
				{
					exec = integrator.apply(m_vThreadData[i].get_solution(), t0+dtcurr, u0, t0);
				}
				catch(ug::UGError& err)
				{
					exec = false;
					error += (1 << i);
					UG_LOG("Step "<< i<< " failed: " << error << " "<< (1<< i) << ":");
					MyPrintError(err);

				}

				if (!exec)
				{

				}

				// switch to "parent" comm
				//mt_env.end();
			} /*for-loop*/



			return error;
		}

		//! e.g. wait for all threads to complete
		void join_integrator_threads()
		{
				// join all threads
				// g.join_all();
		}

		//! Override thread-wise solutions with common solution
		void update_integrator_threads(ConstSmartPtr<grid_function_type> ucommon, number t)
		{
			const int nstages = m_vThreadData.size()-1;
			for (int i=nstages; i>=0; --i)
			{
				UG_ASSERT(m_vThreadData[i].get_solution()->size()==ucommon->size(), "LIMEX: Vectors must match in size!")
				*m_vThreadData[i].get_solution() = *ucommon;
			}
		}



public:
		//! integrating from t0 -> t1
		bool apply(SmartPtr<grid_function_type> u, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
		{
#ifdef UG_OPENMP
			// create multi-threading environment
			//int nt = std::min(omp_get_max_threads(), m_nstages);

#endif

			// NOTE: we use u as common storage for future (and intermediate) solution(s)
			*u = *u0;

			// initialize integrator threads
			// (w/ solutions)
			init_integrator_threads(u0);


			// write_debug
			char name[40];
			for (unsigned int i=0; i<m_vThreadData.size(); ++i)
			{
				sprintf(name, "Limex_Init_iter%03d_stage%03d", 0, i);
				write_debug(*m_vThreadData[i].get_solution(), name);
			}




			number t = t0;
			double dtcurr = ITimeIntegrator<TDomain, TAlgebra>::get_time_step();

			const size_t kmax = m_vThreadData.size();   // maximum number of stages
			size_t qpred = kmax-1;  	 						// predicted optimal order
			size_t qcurr = qpred;

			// double lambda=1.0;   						// step length increase/decrease

			// time integration loop
			SmartPtr<grid_function_type> ubest = SPNULL;
			int limex_step = 1;
			size_t limex_total = 1;
			size_t ntest;    // active number of stages <= kmax
			while ((t < t1) && ((t1-t) > base_type::m_precisionBound))
			{
				int err = 0;

				UG_DLOG(LIB_LIMEX, 5, "+++ LimexTimestep +++" << limex_step << "\n");

				// determine step size
				number dt = std::min(dtcurr, t1-t);
				UG_COND_THROW(dt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!");

				// number of stages to investigate
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
						sprintf(name, "Limex_BeforeSerial_iter%03d_stage%03lu_total%04lu", limex_step, i, limex_total);
						write_debug(*m_vThreadData[i].get_solution(), name);
				}

				// integrate: t -> t+dt
				err = apply_integrator_threads(dt, u, t, ntest-1);

				// write_debug
				for (size_t i=0; i<ntest; ++i)
				{
					sprintf(name, "Limex_AfterSerial_iter%03d_stage%03lu_total%04lu", limex_step, i, limex_total);
					write_debug(*m_vThreadData[i].get_solution(), name);
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
				{
					// compute extrapolation at t+dtcurr (SERIAL)
					timex_type timex(m_vSteps);
					timex.set_error_estimate(m_spErrorEstimator);
					for (unsigned int i=0; i<ntest; ++i)
					{
						timex.set_solution(m_vThreadData[i].get_solution(), i);
					}
					timex.apply(ntest);

					// write_debug
					for (size_t i=0; i<ntest; ++i)
					{
						sprintf(name, "Limex_Extrapolates_iter%03d_stage%03lu_total%04lu", limex_step, i, limex_total);
						write_debug(*m_vThreadData[i].get_solution(), name);
					}
					limex_total++;

					// obtain sub-diagonal error estimates
					const std::vector<number>& eps = timex.get_error_estimates();
					UG_ASSERT(ntest<=eps.size(), "Huhh: Not enough solutions?");

					// select optimal solution (w.r.t error) AND
					// predict optimal order (w.r.t. workload) for next step
					size_t kbest = find_optimal_solution(eps, ntest, qpred);
					UG_ASSERT(kbest < ntest, "Huhh: Not enough solutions?");


					// best solution
					ubest  = timex.get_solution(kbest).template cast_dynamic<grid_function_type>();
					epsmin = eps[kbest]; /*were: kbest*/

					// check for convergence
					limexConverged = (epsmin <= m_tol);

					// select predicted order for next step
					double dtpred = m_lambda[qpred]*dtcurr;
					UG_LOG("koptim=\t" << kbest << ",\t eps(k)=" << epsmin << ",\t q=\t" << qpred<< "("<<  ntest << "), lambda(q)=" << m_lambda[qpred] << ", alpha(q,q)=" << monitor(qpred, qpred) << "dt(q)=" << dtpred<< std::endl);

					// EXTENSIONS: convergence model
					if (limexConverged)
					{
						// a) aim for order increase in next step
						if ((qpred+1==ntest)  /* increase by one possible? */
							 && (kmax>ntest)) /* still below max? */
						{
							const double alpha = monitor(qpred, qpred+1);
							UG_LOG("CHECKING for order increase: "<< m_costA[qpred] << "*" << alpha << ">" << m_costA[qpred+1]);
							// check, whether further increase could still be efficient
							if (m_costA[qpred] * alpha > m_costA[qpred+1])
							{
								qpred++;    			// go for higher order
								if (m_greedyOrderIncrease >0.0) {
								  dtpred *= m_greedyOrderIncrease*alpha;		// & adapt time step  // TODO: check required!
								}
								UG_LOG("... yes.\n")

							} else {
								UG_LOG("... nope.\n")
							}

						}


						// b) monitor convergence (a-priori check!)

					}

					// parameters for subsequent step
					// step length increase/reduction
					/*double dtpred = std::min(dtcurr*lambda,
							dtcurr*itime_integrator_type::get_increase_factor());*/
					dtcurr = std::min(dtpred, itime_integrator_type::get_dt_max());


				}
				else
				{
					// solver failed -> cut time step
					dtcurr *= m_sigmaReduction;
				}


				if ((err==0) && limexConverged)
				{
					// ACCEPT time step
					UG_LOG("+++ LimexTimestep +++" << limex_step << " ACCEPTED"<< std::endl);
					UG_LOG("LIMEX-ACCEPTING:\t" << t <<"\t"<< dt << "\t" << dtcurr << "\tq=\t" << qcurr+1 << std::endl);

					// copy best solution
					UG_ASSERT(ubest.valid(), "Huhh: Invalid error estimate?");
					*u = *ubest;
					t += dt;

					// make sure that all threads continue
					// with identical initial value u(t)
					// update_integrator_threads(ubest, t);


					// working on last row => increase order
					//if (ntest == q+1) ntest++;

					// post process
					itime_integrator_type::notify_step_postprocess(ubest, limex_step++, t, dt);
				}
				else
				{
					// DISCARD time step
					UG_LOG("+++ LimexTimestep +++" << limex_step << " FAILED" << std::endl);
					UG_LOG("LIMEX-REJECTING:\t" << t <<"\t"<< dt << "\t" << dtcurr << std::endl);

				}

				// SERIAL EXECUTION: END
				///////////////////////////////////////
				update_integrator_threads(u, t);



				// SOLVE

				// ESTIMATE

				// MARK



				// REFINE


			} // time integration loop


			return true;
		} // apply


		number 	get_cost(size_t i) { return m_costA[i]; }
		number 	get_gamma(size_t i) { return m_gamma[i]; }
		number 	get_workload(size_t i) { return m_workload[i]; }


	protected:
		number& monitor(size_t k, size_t q) { return m_monitor[k+(m_nstages+1)*q]; }

		/// aux: compute exponents gamma_k (for roots)
		void init_gamma()
		{
			for (size_t k=0; k<=m_nstages; ++k)
			{
				m_gamma[k] = k+1;
			}
		}

  /// Updating workloads A_i for computing T_ii
  //  (depends on m_vSteps, which must have been initialized!)
  void update_cost()
  {
    //UG_LOG("A_0="<< m_vSteps[0] << std::endl);
    m_costA[0] = (1.0)*m_vSteps[0];
    for (size_t i=1; i<=m_nstages; ++i)
      {
	m_costA[i] = m_costA[i-1] + (1.0)*m_vSteps[i];
	//UG_LOG("A_i="<< m_vSteps[i] << std::endl);
      }
  }

			/// convergence monitor
			// (depends on cost, which must have been initialized!)
			void update_monitor()
			{


				for (size_t k=0; k<=m_nstages; ++k)
				{
					UG_LOG("A[k]=" << m_costA[k] << ", gamma[k]=" << m_gamma[k] << "\t");
					for (size_t q=0; q<=m_nstages; ++q)
					{
						// Deuflhard: order and stepsize, ... eq. (3.7)
						double gamma = (m_costA[k] - m_costA[0] + 1.0)/(m_costA[q] - m_costA[0] + 1.0);
						double alpha = pow(m_tol, gamma);

						// for fixed order q, the monitor indicates the performance penalty compared to a strategy using k stages only
						// cf. eq. (4.6)
						monitor(k,q) = pow(alpha/(m_tol*m_rhoSafety), 1.0/m_gamma[k]);
						UG_LOG(monitor(k,q) << "[" << pow(alpha/(m_tol), 1.0/m_gamma[k]) << "]" << "\t");
						// UG_LOG(  << "\t");



					}
					UG_LOG(std::endl);
				}


			}

			// Find row k=1, ..., ntest-1 minimizing (estimated) error eps[kmin]
			// Also: predict column q with minimal workload W_{k+1,k} = A_{k+1} * lambda_{k+1}
			size_t find_optimal_solution(const std::vector<number>& eps, size_t ntest, /*size_t &kf,*/ size_t &qpred)
			{

				const size_t qold=qpred;

				size_t kbest = 1;
				qpred = 1;

				size_t k=1;
				m_lambda[k] = pow(m_rhoSafety*m_tol/eps[k], 1.0/m_gamma[k]);   // 1/epsilon(k)
				m_workload[k] = m_costA[k]/m_lambda[k];
				UG_LOG("k=" << k << ": eps=" << eps[k]  << ", lambda(k)=" <<m_lambda[k]  << ", epsilon(k)=" <<1.0/m_lambda[k] << "<= alpha(k, qcurr)=" << monitor(k-1,qold) << "< alpha(k, qcurr+1)=" << monitor(k-1,qold+1) <<", A="<< m_costA[k] << ", W="<< m_workload[k] <<std::endl);

				for (k=2; k<ntest; ++k)
				{
					m_lambda[k] = pow(m_rhoSafety*m_tol/eps[k], 1.0/m_gamma[k]);
					m_workload[k] = m_costA[k]/m_lambda[k];
					UG_LOG("k=" << k << ": eps=" << eps[k]  << ", lambda(k)=" <<m_lambda[k]   << ", epsilon(k)=" <<1.0/m_lambda[k] << "<= alpha(k, qcurr)=" << monitor(k-1,qold) << "< alpha(k, qcurr+1)=" << monitor(k-1,qold+1) <<", A="<< m_costA[k] << ", W="<< m_workload[k] <<std::endl);

					// TODO: Convergence monitor

					qpred = (m_workload[qpred] > m_workload[k]) ? k : qpred;
					kbest = (eps[kbest] > eps [k]) ? k : kbest;
				}

				return kbest;
			}


protected:

		double m_tol;
		double m_rhoSafety;
		double m_sigmaReduction;
		SmartPtr<error_estim_type> m_spErrorEstimator;     // (smart ptr for) error estimator


		unsigned int m_nstages; 							// stages in Aitken-Neville
		std::vector<size_t> m_vSteps;						// generating sequence for extrapolation
		std::vector<ThreadData> m_vThreadData;				// vector with thread information




		std::vector<number> m_gamma;			/// gamma_i: exponent

		std::vector<number> m_costA;			// A_i: cost (for completing stage i)
		std::vector<number> m_monitor;			/// convergence monitor \alpha

		std::vector<number> m_workload;
		std::vector<number> m_lambda;

  double m_greedyOrderIncrease;

};

} // namespace ug

#endif /* TIME_INTEGRATOR_HPP_ */
