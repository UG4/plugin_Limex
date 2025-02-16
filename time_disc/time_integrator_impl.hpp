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

#ifndef __H__LIMEX__TIME_INTEGRATOR_IMPL_HPP__
#define __H__LIMEX__TIME_INTEGRATOR_IMPL_HPP__

namespace ug {

template<typename TDomain, typename TAlgebra>
bool LinearTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{

	LIMEX_PROFILE_FUNC()
	UG_COND_THROW(!base_type::m_spLinearSolver.valid(), "Linear solver invalid");

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;

	// create solution vector & right hand side
	SmartPtr<typename base_type::vector_type> uold= u0->clone();
	SmartPtr<typename base_type::vector_type> b= u0->clone_without_values();

	// solution time series
	SmartPtr<vector_time_series_type> m_spSolTimeSeries;
	m_spSolTimeSeries=make_sp(new vector_time_series_type());
	m_spSolTimeSeries->clear();
	m_spSolTimeSeries->push(uold, t0);

	// create matrix operator
	SmartPtr<typename base_type::assembled_operator_type> spAssOp=make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));

	// integrate
	if(!base_type::m_bNoLogOut)
		UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"]\n");

	double t = t0;
	number dt_assembled = -1.0;   // invalid
	int step = 1;

	number currdt = base_type::m_dt;

	while((t < t1) && (t1-t > base_type::m_precisionBound))
	{

		if(!base_type::m_bNoLogOut)
		{	UG_LOG("+++ Timestep +++" << step << "\n"); }

		// determine step size
		number dt = std::min(currdt, t1-t);

		// prepare step
		tdisc.prepare_step(m_spSolTimeSeries, dt);
		if (fabs(dt-dt_assembled) > base_type::m_precisionBound)
		{
			// re-assemble operator
			if(!base_type::m_bNoLogOut)
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
			VecAssign(*tmp, *u1.template cast_dynamic<typename base_type::vector_type>());
			m_spSolTimeSeries->push_discard_oldest(tmp, t);

			this->notify_finalize_step(u1, step++, t+dt, dt);
		}
		else
		{
			// DISCARDING
			currdt *= 0.5;
		}

	}

	this->notify_end(u1, step, t1, currdt);

	return true;
};


template<typename TDomain, typename TAlgebra>
bool ConstStepLinearTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	LIMEX_PROFILE_FUNC()
	UG_COND_THROW(!base_type::m_spLinearSolver.valid(), "Linear solver invalid");

	// short-cuts
	GridLevel const &gl = u0->grid_level();
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;

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
	if(!base_type::m_bNoLogOut)
	{
		UG_LOG("+++ Integrating: [\t"<< t0 <<"\t, \t"<< t1 <<"\t] with dt=\t" << currdt << "("<< numSteps<< " iters)" << std::endl);
	}
	 
	 // integrate
	 for(int step = 1; step<=numSteps; ++step)
	 {
	     // determine step size
	     // number dt = std::min(currdt, t1-t);
		 const number dt = currdt;

		if(!base_type::m_bNoLogOut)
		{
			UG_LOG("+++ Const timestep +++" << step<< "(t=" << t << ", dt=" << dt << ")"<< std::endl);
		}
		this->notify_init_step(u1, step, t, dt);
		
		// prepare step
		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 if (spAssOp==SPNULL)
		 {
			 // Assemble operator.
			if(!base_type::m_bNoLogOut) UG_LOG("+++ Assemble (t=" << t << ", dt=" << dt <<")" << std::endl);

			 spAssOp=make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));
			 tdisc.assemble_linear(*spAssOp, *b, gl);
			 (base_type::m_spLinearSolver)->init(spAssOp, *u1);
		 }
		 else
		 {
			 // Recycle existing operator.
			 // std::cerr << "Recycling timestep " << step << "\n";
			 tdisc.assemble_rhs(*b, gl);
		 }

		 // execute step
		 if (base_type::m_spLinearSolver->apply(*u1, *b))
		 {
			 // ACCEPTING: push updated solution into time series
			 t += dt;
			 SmartPtr<typename base_type::vector_type> tmp = m_spSolTimeSeries->oldest();
			 VecAssign(*tmp,  *u1.template cast_dynamic<typename base_type::vector_type>());
			 m_spSolTimeSeries->push_discard_oldest(tmp, t);
			 this->notify_finalize_step(u1, step, t, dt);  // t updated already!
		 }
		 else
		 {
			UG_THROW("ConstStepLinearTimeIntegrator::apply failed!!!");
			this->notify_rewind_step(u1, step, t+dt, dt);
		 }

	 }

	 this->notify_end(u1, numSteps, t1, currdt);

	 return true;
};

template<typename TDomain, typename TAlgebra>
bool TimeIntegratorLinearAdaptive<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	// short-cuts
	GridLevel const &gl = u0->grid_level();
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
	time_disc_type &tdisc2 = *m_spTimeDisc2;

	// create solution vector & right hand side
	SmartPtr<typename base_type::vector_type> uold= u0->clone();
	SmartPtr<typename base_type::vector_type> b= u0->clone_without_values();

	// create matrix operator
	SmartPtr<typename base_type::assembled_operator_type> spAssOp=make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));

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
	if (m_dtmin <= 0.0) { m_dtmin = (t1 - t0)/1.0e+5; }
	if (m_dtmax <= 0.0) { m_dtmax = (t1 - t0)/10.0; }

	// Aitken Neville extrapolation
	const size_t tsteps[2] = {1,2};
	std::vector<size_t> vsteps (tsteps, tsteps+2);
	AitkenNevilleTimex<typename base_type::vector_type> timex(vsteps);

	// integrate
	if(!base_type::m_bNoLogOut)
		UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"]\n");

	 double t = t0;
	 int step = 0;

	 number dt = base_type::m_dt;
	 while((t < t1) && (t1-t > base_type::m_precisionBound))
	 {
	   // step: t -> t+dt
	   bool bSuccess = false;
	   while (!bSuccess){
		
	     // determine step size
	     if (dt<m_dtmin){
				if(!base_type::m_bNoLogOut)
	       UG_LOG("Step size below minimum")
		 }
		 dt = std::min(dt, t1-t);
	   
		 // basic step
			if(!base_type::m_bNoLogOut)
			 UG_LOG("+++ Timestep: " << ++step << "\n");

		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 tdisc.assemble_linear(*spAssOp, *b, gl);
		 base_type::m_spLinearSolver->init(spAssOp, *u1);
		 base_type::m_spLinearSolver->apply(*u1, *b);
		 
		 
		 // control 1/2
			if(!base_type::m_bNoLogOut)
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
		 number eps = timex.get_error_estimates()[0];
		 number lambda = std::pow(0.8* m_tol/eps, 0.5); 

		 number dtEst= dt*lambda; 
		 dtEst = std::min(dtEst, 1.5*dt); 
		 dtEst = std::min(dtEst, m_dtmax); 
		 
		 dt = dtEst;
		 if (eps <= m_tol)
		 {
		   // ACCEPT STEP (and thus solution u2)
			if(!base_type::m_bNoLogOut)
		   UG_LOG("ACCEPTING solution, dtnew=" << dt);

		   bSuccess = true;
		 }
		 else
		 {
		   // DISCARD step
			if(!base_type::m_bNoLogOut)
		   UG_LOG("DISCARDING solution, dtnew=" << dt);

		   // => reset solutions
		   VecAssign(*u1.template cast_dynamic<typename base_type::vector_type>(), *uold);
		   
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
	   VecAssign(*tmp, static_cast<typename base_type::vector_type> (*u2));
	   m_spSolTimeSeries->push_discard_oldest(tmp, t);

	 }

	 return true;
};

} // namespace ug

#endif /* __H__LIMEX__TIME_INTEGRATOR_IMPL_HPP__ */
