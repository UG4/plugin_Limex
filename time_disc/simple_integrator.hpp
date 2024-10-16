#include "lib_algebra/operator/debug_writer.h"

namespace ug {

/// Integrate (a non-linear problem) over a given time interval
template<class TDomain, class TAlgebra>
class SimpleTimeIntegrator :
		public INonlinearTimeIntegrator<TDomain, TAlgebra>,
		public ITimeDiscDependentObject<TAlgebra>,
		public DebugWritingObject<TAlgebra>
{
protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;

public:
	typedef INonlinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef ITimeDiscretization<TAlgebra> time_disc_type;
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef IGridFunctionSpace<grid_function_type> grid_function_space_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

	// constructor
	SimpleTimeIntegrator (SmartPtr<time_disc_type> tDisc)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc),
	  m_spBanachSpace(new AlgebraicSpace<grid_function_type>() ),
	  m_spDerivative(SPNULL), m_initial_consistency_error(0.0)

	{}

	SimpleTimeIntegrator (SmartPtr<time_disc_type> tDisc, SmartPtr<grid_function_space_type> spSpace)
		: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc),
		  m_spBanachSpace(spSpace),
		  m_spDerivative(SPNULL), m_initial_consistency_error(0.0)
	{}


	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
	{
		time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
		if (tdisc.num_stages() == 1)
			return apply_single_stage(u1,t1,u0,t0);
		else
			return apply_multi_stage(u1,t1,u0,t0);
	}

	void set_derivative(SmartPtr<grid_function_type> udot)
	{ m_spDerivative = udot; }

	SmartPtr<grid_function_type> get_derivative()
	{ return m_spDerivative; }

	number get_consistency_error() const
	{ return m_initial_consistency_error; }

	void set_banach_space(SmartPtr<IGridFunctionSpace<grid_function_type> > spSpace)
	{ m_spBanachSpace = spSpace; }

protected:
	bool apply_single_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
	bool apply_multi_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

	inline bool hasTerminated(double tCurrent, double tStart, double tFinal) const
	{
	 	/*return (! ((tCurrent < tFinal) && (tFinal-tCurrent > base_type::m_precisionBound)));*/
		return tCurrent >= tFinal || tFinal-tCurrent < (tFinal-tStart)*base_type::m_precisionBound;
	}

	/// metric
	SmartPtr<IGridFunctionSpace<grid_function_type> > m_spBanachSpace;

	SmartPtr<grid_function_type> m_spDerivative;

	number m_initial_consistency_error;
};


template<typename TDomain, typename TAlgebra>
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_single_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
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
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_multi_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
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
