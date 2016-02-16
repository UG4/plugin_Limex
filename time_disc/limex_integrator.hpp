/*
 * time_integrator.hpp
 *
 *  Created on: 15.08.2014
 *      Author: anaegel
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
: public INonlinearTimeIntegrator<TDomain, TAlgebra>
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
		typedef SimpleTimeIntegrator<TDomain, TAlgebra> time_integrator_type;

		class ThreadData

		{
			//typedef boost::thread thread_type;
		public:

			ThreadData(SmartPtr<timestep_type> spTimeStep)
			: m_stepper(spTimeStep)
			{

			}

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



		protected:
			 // includes time step series
			SmartPtr<timestep_type> m_stepper;
			SmartPtr<grid_function_type> m_sol;
			SmartPtr<solver_type> m_solver;

		};

		typedef std::vector<SmartPtr<ThreadData> > thread_vector_type;


		LimexTimeIntegrator(int nstages)
		: m_nstages(nstages), m_tol(0.01), m_rhoSafety(0.8)
		{
			m_vThreadData.reserve(m_nstages);
			m_vSteps.reserve(m_nstages);
			// create integrators
			for (unsigned int i=0; i<m_nstages; ++i){

				m_vSteps.push_back(i+1);
			}
		}

		void set_tolerance(double tol)
		{m_tol = tol;}

		void add_stage(int i, int nsteps, SmartPtr<domain_discretization_type> spDD, SmartPtr<solver_type> solver)
		{
			if (i>=m_nstages) return;

			m_vSteps.push_back(nsteps);
			m_vThreadData.push_back(ThreadData(make_sp(new timestep_type(spDD))));

			m_vThreadData.back().set_solver(solver);

			UG_ASSERT(solver.valid(), "Huhh: Need to supply solver!")
			UG_ASSERT(m_vThreadData.back().get_solver().valid(), "Huhh: Need to supply solver!")

		}

		void apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
		{

#ifdef UG_OPENMP
			// create multi-threading environment
			//int nt = std::min(omp_get_max_threads(), m_nstages);

#endif


			double dtcurr = ITimeIntegrator<TDomain, TAlgebra>::get_time_step();
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

			//unsigned int k = m_vThreadData.size()-1;
			const int nstages = m_vThreadData.size()-1;

			for (int i=nstages; i>=0; --i)
			{
				m_vThreadData[i].set_solution(u1->clone());
			}

//#pragma omp parallel shared (u1)
			{
				// PARALLEL execution

/*
				int tn = omp_get_thread_num();
				int nt = omp_get_num_threads();

				omp_set_num_threads(nstages);
*/
				int i;
			//	#pragma omp for private(i) // shared (nstages, u1) schedule(static)
				for (i=nstages; i>=0; --i)
				{
					/*
					std::cerr << "I am " << tn << " of " << nt << " ("<< i<< "/" << nstages<<")!" << std::endl;

					UGMultiThreadEnvironment mt_env;
					 */

					// copy data to private structure (one-to-many)
					//m_vThreadData[i].set_solution(u1->clone());

					// switch to "child" comm
					// mt_env.begin();

					// integrate
					time_integrator_type integrator(m_vThreadData[i].get_time_stepper());
					integrator.set_time_step(dtcurr/m_vSteps[i]);
					integrator.set_solver(m_vThreadData[i].get_solver());
					integrator.apply(m_vThreadData[i].get_solution() , t1, u0, t0);

					// switch to "parent" comm
					//mt_env.end();
				} /*for-loop*/

			} /*omp parallel */
			// join all threads
			// g.join_all();


			// SERIAL EXECUTION:
			// compute extrapolation
			UG_ASSERT(m_vSteps.size()==m_vThreadData.size(), "Huhh: sizes do not match!");
			timex_type timex(m_vSteps);
			for (unsigned int i=0; i<m_vThreadData.size(); ++i)
			{
				timex.set_solution(m_vThreadData[i].get_solution(), i);
			}
			timex.apply();
			double eps = timex.get_error_estimate();
			double lambda=pow(m_rhoSafety*m_tol/eps, 0.5);

			std::cerr << "eps=" << eps << "lambda=" << lambda << std::endl;


			if (eps<= m_tol)
			{
				// ACCEPT time step
				double dtest = std::min(dtcurr*lambda, dtcurr*m_dtBounds.get_increase_factor());
				dtcurr = std::min(dtest, m_dtBounds.get_dt_max());

			} else
			{
				// DISCARD time step
				double dtest = std::min(dtcurr*lambda, dtcurr*m_dtBounds.get_reduction_factor());
				UG_ASSERT (eps<= m_tol, "Time step failed!")
			}



			// copy solution
			*u1 = *(m_vThreadData[nstages].get_solution());

			// notify_step_postprocess(u1, 1, 10.0*dtcurr, dtcurr);

		}

protected:
		// SmartPtr<domain_discretization_type> m_spDD;		// underlying domain disc

		unsigned int m_nstages; 							// stages in Aitken-Neville
		std::vector<size_t> m_vSteps;						// generating sequence for extrapolation
		std::vector<ThreadData> m_vThreadData;				// vector with thread information

		TimeStepBounds m_dtBounds;
		double m_tol;
		double m_rhoSafety;

};

} // namespace ug

#endif /* TIME_INTEGRATOR_HPP_ */
