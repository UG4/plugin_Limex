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
		: m_nstages(nstages),
		  m_tol(0.01),
		  m_rhoSafety(0.8),
		  m_costA(m_nstages),
		  m_gamma(m_nstages),
		  m_alpha(m_nstages*m_nstages),
		  m_lambda(m_nstages),
		  m_workload(m_nstages)
		{
			m_vThreadData.reserve(m_nstages);
			m_vSteps.reserve(m_nstages);

			// create integrators
			/*for (unsigned int i=0; i<m_nstages; ++i){

				 m_vSteps.push_back(i+1);
			}*/


			init_gamma();



		}

		/// tolerance
		void set_tolerance(double tol) { m_tol = tol;}

		void add_error_estimator(SmartPtr<error_estim_type> spErrorEstim)
		{ m_spErrorEstimator = spErrorEstim; }

		void add_stage(int i, int nsteps, SmartPtr<domain_discretization_type> spDD, SmartPtr<solver_type> solver)
		{
			if (i>=m_nstages) return;

			m_vSteps.push_back(nsteps);

			m_vThreadData.push_back(ThreadData(make_sp(new timestep_type(spDD))));
			m_vThreadData.back().set_solver(solver);

			UG_ASSERT(solver.valid(), "Huhh: Need to supply solver!");
			UG_ASSERT(m_vThreadData.back().get_solver().valid(), "Huhh: Need to supply solver!");


		}

protected:
		//! initialize integrator threads (w/ solutions)
		void init_integrator_threads(SmartPtr<grid_function_type> u1)
		{
			const int nstages = m_vThreadData.size()-1;
			for (int i=nstages; i>=0; --i)
			{
				m_vThreadData[i].set_solution(u1->clone());
			}
		}


		//! (tentatively) apply integrators
		// TODO: PARALLEL execution?
		int apply_integrator_threads(number dtcurr, SmartPtr<grid_function_type> u0, number t0, size_t nstages)
		{
			// compute cost A_i (alternative: measure times?)
			init_cost();

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
					UG_LOG("Step "<< i<< " failed: " << error << " "<< (1<< i));

				}

				if (!exec)
				{

				}

				// switch to "parent" comm
				//mt_env.end();
			} /*for-loop*/



			return error;
		}

		void update_integrator_threads(ConstSmartPtr<grid_function_type> u, number t)
		{

		}

		void sync_integrator_threads()
		{
			// join all threads
			// g.join_all();
		}

public:
		//! integrating from t0 -> t1
		bool apply(SmartPtr<grid_function_type> u, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
		{


#ifdef UG_OPENMP
			// create multi-threading environment
			//int nt = std::min(omp_get_max_threads(), m_nstages);

#endif

			// IMPORTANT NOTE: we use u as common storage for future (and intermediate) solution(s)
			*u = *u0;


			// initialize integrator threads
			// (w/ solutions)
			init_integrator_threads(u);


			// write_debug
			char name[40];
			for (unsigned int i=0; i<m_vThreadData.size(); ++i)
			{
				sprintf(name, "Limex_Init_iter%03d_stage%03d", 0, i);
				write_debug(*m_vThreadData[i].get_solution(), name);
			}

			// initialize integrator threads
			// (w/ solutions)

			/*m_vThreadData[k--].solution() = u1;
			while(k>=0)
			{
				m_vThreadData[k--].solution() = u1->clone();
			}
*/
			// create (& execute) threads
			/*boost::thread_group g;
			typename thread_vector_type::reverse_iterator rit=m_vThreadData.rbegin();
			for (rit++; rit!= m_vThreadData.rend(); ++rit)
			{

				boost::thread *t =new boost::thread(boost::bind(&ThreadSafeTimeIntegrator::apply, *rit));
				//g.add_thread(t);

				g.create_thread(boost::bind(&ThreadSafeTimeIntegrator::apply, *rit));

			}*/


			number t = t0;
			double dtcurr = ITimeIntegrator<TDomain, TAlgebra>::get_time_step();

			const size_t kmax = m_vThreadData.size();   // maximum number of stages
			size_t q = 0;    // current order
			size_t kf = 2;   // upper end

			// time integration loop
			int limex_step = 1;
			while ((t < t1) && ((t1-t) > base_type::m_precisionBound))
			{
				int err = 0;

				UG_DLOG(LIB_LIMEX, 5, "+++ LimexTimestep +++" << limex_step << "\n");

				// determine step size
				number dt = std::min(dtcurr, t1-t);
				UG_COND_THROW(dt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!");

				// number of stages to investigate
				kf = std::min(kmax, q+2);
				UG_LOG("kf="<< kf << std::endl);

				// checks
				UG_ASSERT(m_vSteps.size() >= kf, "Huhh: sizes do not match: " << m_vSteps.size() << "<"<<kf);
				UG_ASSERT(m_vThreadData.size() >= kf, "Huhh: sizes do not match: " << m_vThreadData.size() << "< "<< kf);

				// write_debug
				for (unsigned int i=0; i<kf; ++i)
				{
						sprintf(name, "Limex_BeforeSerial_iter%03d_stage%03d", limex_step, i);
						write_debug(*m_vThreadData[i].get_solution(), name);
				}

				// integrate: t -> t+dt
				err = apply_integrator_threads(dt, u, t, kf-1);
				sync_integrator_threads();


				// write_debug
				for (unsigned int i=0; i<kf; ++i)
				{
					sprintf(name, "Limex_AfterSerial_iter%03d_stage%03d", limex_step, i);
					write_debug(*m_vThreadData[i].get_solution(), name);
				}

				// post-cond checks
				UG_ASSERT(m_spErrorEstimator.valid(), "Huhh: Invalid Error estimator?");

				double epsq = 0.0;
				SmartPtr<grid_function_type> ubest = SPNULL;

				if (err==0)
				{
					// compute extrapolation at t+dtcurr (SERIAL)
					timex_type timex(m_vSteps);
					timex.set_error_estimate(m_spErrorEstimator);
					for (unsigned int i=0; i<kf; ++i)
					{
						timex.set_solution(m_vThreadData[i].get_solution(), i);
					}
					timex.apply(kf);

					// write_debug
					for (unsigned int i=0; i<kf; ++i)
					{
						sprintf(name, "Limex_Extra_iter%03d_stage%03d", limex_step, i);
						write_debug(*m_vThreadData[i].get_solution(), name);
					}

					// obtain sub-diagonal error estimates
					const std::vector<number>& eps = timex.get_error_estimates();

					UG_ASSERT(kf<=eps.size(), "Huhh: Not enough solutions?");

					// detect optimal order (w.r.t. workload)
					q = compute_optimal_order(eps, kf);
					ubest = timex.get_solution(q).template cast_dynamic<grid_function_type>();

					double lambda = m_lambda[q]; 				// step length increase/reduction
					double dtq = std::min(dtcurr*lambda,
							dtcurr*itime_integrator_type::get_increase_factor());
					epsq= eps[q];
					UG_LOG("order=" << q<< "("<<  kf <<"), eps(q)=" << epsq << ", lambda0(q)=" << m_rhoSafety*m_tol/epsq << ", lambda(q)=" << lambda << "dtest=" << dtq<< std::endl);
					dtcurr = std::min(dtq, itime_integrator_type::get_dt_max());

				}
				else
				{
					// solver failed -> cut time step
					dtcurr *=0.5;
				}
				// err == 0


				if ((err==0) && (epsq <= m_tol))
				{
					// ACCEPT time step
					UG_LOG("+++ LimexTimestep +++" << limex_step << " ACCEPTED"<< std::endl);
					UG_LOG("LIMEX-ACCEPTING:\t" << t <<"\t"<< dt << "\t" << dtcurr << std::endl);

					// copy best solution
					UG_ASSERT(ubest.valid(), "Huhh: Invalid error estimate?");
					*u = *ubest;
					t += dt;

					// working on last row => increase order
					// if (kf == q+1) kf++;

					// post process
					itime_integrator_type::notify_step_postprocess(u, limex_step++, t, dt);
				}
				else
				{
					// DISCARD time step
					UG_LOG("+++ LimexTimestep +++" << limex_step << " FAILED" << std::endl);
					UG_LOG("LIMEX-REJECTING:\t" << t <<"\t"<< dt << "\t" << dtcurr << std::endl);

				}

			} // time integration loop



			// notify_step_postprocess(u, 1, 10.0*dtcurr, dtcurr);

			return true;
		} // apply


		number get_cost(size_t i) { return m_costA[i]; }
		number get_gamma(size_t i) { return m_gamma[i]; }
		number get_alpha(size_t k, size_t q) { return m_alpha[k]; }
		number get_workload(size_t i) { return m_workload[i]; }


	protected:
			/// aux: compute workloads A_i for computing T_ii
			void init_cost()
			{
				//UG_LOG("A_0="<< m_vSteps[0] << std::endl);
				m_costA[0] = (1.0)*m_vSteps[0];
				for (size_t i=1; i<m_nstages; ++i)
				{
					m_costA[i] = m_costA[i-1] + (1.0)*m_vSteps[i];
					//UG_LOG("A_i="<< m_vSteps[i] << std::endl);
				}
			}

			/// aux: compute exponents gamma_k (for roots)
			void init_gamma()
			{
				for (size_t k=0; k<m_nstages; ++k)
				{
					m_gamma[k] = k+1;
				}
			}

			// Find column k=0, ..., kf-1
			// minimizing W_{k+1,k} = A_{k+1} * lambda_{k+1},
			size_t compute_optimal_order(const std::vector<number>& eps, size_t kf)
			{


				m_lambda[0] = pow(m_rhoSafety*m_tol/eps[0], 1.0/m_gamma[0]);
				m_workload[0] = m_costA[0]/m_lambda[0];
				UG_LOG("k=0, A="<< m_costA[0] << ", W="<< m_workload[0] << ", \lambda=" <<m_lambda[0] <<std::endl);

				size_t q = 0;
				number Wq = m_workload[0];

				for (size_t k=1; k<kf; ++k)
				{
					m_lambda[k] = pow(m_rhoSafety*m_tol/eps[k], 1.0/m_gamma[k]);
					m_workload[k] = m_costA[k]/m_lambda[k];
					UG_LOG("k=" << k <<", A="<< m_costA[k] << ", W="<< m_workload[k] << ", \lambda=" <<m_lambda[k] <<std::endl);

					if (Wq > m_workload[k]) {
						q = k;
						Wq = m_workload[k];
					}
				}
				return q;
			}


protected:
		// SmartPtr<domain_discretization_type> m_spDD;		// underlying domain disc

		unsigned int m_nstages; 							// stages in Aitken-Neville
		std::vector<size_t> m_vSteps;						// generating sequence for extrapolation
		std::vector<ThreadData> m_vThreadData;				// vector with thread information



		std::vector<number> m_costA;			// A_i: cost (for completing stage)
		std::vector<number> m_workload;


		std::vector<number> m_gamma;			/// gamma_i: exponent
		std::vector<number> m_alpha;
		std::vector<number> m_lambda;


		//TimeStepBounds m_dtBounds;
		double m_tol;
		double m_rhoSafety;
		SmartPtr<error_estim_type> m_spErrorEstimator;     // (smart ptr for) error estimator
};

} // namespace ug

#endif /* TIME_INTEGRATOR_HPP_ */
