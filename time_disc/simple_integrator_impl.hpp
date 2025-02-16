/*
 * Copyright (c) 2025: Goethe University Frankfurt
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

#ifndef __H__LIMEX__SIMPLE_INTEGRATOR_IMPL_HPP__
#define __H__LIMEX__SIMPLE_INTEGRATOR_IMPL_HPP__

namespace ug {

template<typename TDomain, typename TAlgebra>
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_single_stage
(
	SmartPtr<grid_function_type> u1,
	number t1,
	ConstSmartPtr<grid_function_type> u0,
	number t0
)
{
	LIMEX_PROFILE_FUNC()

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
	typename base_type::solver_type &solver = *base_type::m_spSolver;

	// create solution vector & right hand side
	SmartPtr<grid_function_type> uold;

	// init solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;        ///< contains all solutions compute so far
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->clear();
	m_spSolTimeSeries->push(u0->clone(), t0);

	// init solver (and matrix operator)
	SmartPtr<typename base_type::assembled_operator_type> spAssOp;
	spAssOp = make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));
	solver.init(spAssOp);

	// integrate
	double t = t0;
	number currdt = base_type::m_dt;
	int step = 1;

	double final_dt = base_type::m_dt;

	if(!base_type::m_bNoLogOut)
		UG_LOG("+++ Integrating: [\t"<< t0 <<"\t, \t"<< t1 <<"\t] with\t" << currdt <<"\n");

	while(!hasTerminated(t, t0, t1)) {
		if(!base_type::m_bNoLogOut)
			UG_LOG("+++ Timestep +++" << step << "\n");
			
		if (this->debug_writer_valid())
		{
			char debug_name_ext[16]; snprintf(debug_name_ext, 16, "%04d", step);
			this->enter_debug_writer_section(std::string("SimpleTimeIntegrator_step") + debug_name_ext);
		}

		// determine step size
		UG_COND_THROW(currdt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!")
			number dt = std::min(currdt, t1-t);
		final_dt = dt;

		// prepare step
		tdisc.prepare_step(m_spSolTimeSeries, dt);
		if (solver.prepare(*u1) == false)
		{
			if(!base_type::m_bNoLogOut)
				UG_LOG("Initialization failed! RETRY");

			currdt *= base_type::get_reduction_factor();
			continue;
		}
		//UG_LOG("m_spSolTimeSeries.size="<< m_spSolTimeSeries->size());

		// execute step
		if (solver.apply(*u1))
		{
			//
			// ACCEPT step
			//

			// post prcess (e.g. physics)
			if(!base_type::m_bNoLogOut)
			{
				// m_spSolTimeSeries->oldest() actually holds a pointer to a grid function
				// but as the time series does not know this, we have to cast ourselves
				// SmartPtr<grid_function_type> tmpOld = m_spSolTimeSeries->oldest().template cast_static<grid_function_type>();
				this->notify_finalize_step(u1, step, t+dt, dt);
			}

			// update time
			t += dt;

			// push updated solution into time series (and continue)
			//SmartPtr<typename base_type::vector_type> utmp = m_spSolTimeSeries->oldest();
			//VecAssign(*utmp, static_cast<typename base_type::vector_type> (*u1) );
			uold = m_spSolTimeSeries->push_discard_oldest(u1->clone(), t).template cast_static<grid_function_type>();
			}
		else
		{
			//
			// REJECT step
			//
			UG_LOG("Solution failed! RETRY");
			currdt *= base_type::get_reduction_factor();
			continue;
		}

		// consistency check
		if (step == 1 && m_spDerivative.valid())
		{
			UG_LOG("Computing consistency error: "<< std::endl);
			UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*u0),
					"Huhh: Different vectors required!");
			*m_spDerivative = *u0;
			m_initial_consistency_error = m_spBanachSpace->distance(*m_spDerivative, *u1);
		}

		step++;
		// tdisc.finish_step_elem(m_spSolTimeSeries, dt);

		this->leave_debug_writer_section();
	}

	if(base_type::m_bNoLogOut)
	{
		this->notify_finalize_step(u1, /*uold,*/ step, t, final_dt);
	}


	if (m_spDerivative.valid())
	{

		//
		// approximate derivative (by forward difference)
		//
		UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*uold),
				  "Huhh: Different vectors required!");

		VecScaleAdd(static_cast<typename TAlgebra::vector_type&>(*m_spDerivative),
				1.0/final_dt, static_cast<typename TAlgebra::vector_type&>(*u1),
				-1.0/final_dt, static_cast<typename TAlgebra::vector_type&>(*uold));

	}

	m_spSolTimeSeries->clear();

	return true;

};

template<typename TDomain, typename TAlgebra>
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_multi_stage
(
	SmartPtr<grid_function_type> u1,
	number t1,
	ConstSmartPtr<grid_function_type> u0,
	number t0
)
{

	LIMEX_PROFILE_FUNC()

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
	typename base_type::solver_type &solver = *base_type::m_spSolver;

	//using TimeIntegratorSubject<TDomain,TAlgebra>::notify_step_postprocess;

	// create solution vector & right hand side
	SmartPtr<grid_function_type> uold = u0->clone();

	// init solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->clear();
	m_spSolTimeSeries->push(u0->clone(), t0);

	// init solver (and matrix operator)
	SmartPtr<typename base_type::assembled_operator_type> spAssOp=make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));
	solver.init(spAssOp);

	// integrate
	double t = t0;
	number currdt = base_type::m_dt;
	int step = 1;

	double final_dt = base_type::m_dt;

	if(!base_type::m_bNoLogOut)
		UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"] with " << currdt <<"\n");

	while(!hasTerminated(t, t0, t1))
	{
		if(!base_type::m_bNoLogOut)
			UG_LOG("++++++ TIMESTEP " << step++ << " BEGIN (current time: " << t << ") ++++++\n");

		// determine step size
		UG_COND_THROW(currdt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!")
		number dt = std::min(currdt, t1-t);

		final_dt = dt;

		double told = t;

		const int num_stages = tdisc.num_stages();
		int s=1;
		do // 'for-each-stage' loop
		{ 

			// for (int s=1; s<=num_stages; ++s)
			if(!base_type::m_bNoLogOut)
				UG_LOG("+++ STAGE "<< s << " BEGIN +++\n");

			// set stage
			tdisc.set_stage(s);

			// prepare step
			tdisc.prepare_step(m_spSolTimeSeries, dt);
			if (solver.prepare(*u1) == false) break;

			// execute step
			if (!solver.apply(*u1)) break;

			// stage was successful:
			// a. update (intermediate) time
			t = tdisc.future_time();

			// b. push updated solution into time series (and continue)
			SmartPtr<typename base_type::vector_type> oldest= m_spSolTimeSeries->oldest();
			VecAssign(*oldest, static_cast<typename base_type::vector_type> (*u1));
			m_spSolTimeSeries->push_discard_oldest(oldest, t);

			// c. output
			if(!base_type::m_bNoLogOut)
				UG_LOG("+++ STAGE "<< s << " END +++\n");

		} while ((++s) <=num_stages);


		if (s<=num_stages)
		{
			// REJECT time step
			if(!base_type::m_bNoLogOut)
				UG_LOG("Solution failed! RETRY");

			currdt *= this->get_reduction_factor();
			t = told;
			continue;
		}
		else
		{
			if(!base_type::m_bNoLogOut)
			{
				this->notify_finalize_step(u1, /*uold, */step, t, dt);
			}

			// ACCEPT time step
			if (!hasTerminated(t, t0, t1))
				*uold = *u1;   // save solution (but not in last step)
			// tdisc.finish_step_elem(m_spSolTimeSeries, dt);

			if(!base_type::m_bNoLogOut)
				UG_LOG("++++++ TIMESTEP " << step++ << " END   (current time: " << t << ") ++++++\n");
		}
	}

	if(base_type::m_bNoLogOut)
	{
		this->notify_finalize_step(u1, /*uold,*/ step, t, final_dt);
	}
	return true;
};

} // ug

#endif // __H__LIMEX__SIMPLE_INTEGRATOR_IMPL_HPP__
