
namespace ug {

/// Integrate (a non-linear problem) over a given time interval
template<class TDomain, class TAlgebra>
class SimpleTimeIntegrator :
		public INonlinearTimeIntegrator<TDomain, TAlgebra>,
		public ITimeDiscDependentObject<TAlgebra>
{
protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;

	using TimeIntegratorSubject<TDomain, TAlgebra>::notify_init_step;
	using TimeIntegratorSubject<TDomain, TAlgebra>::notify_preprocess_step;
	using TimeIntegratorSubject<TDomain, TAlgebra>::notify_postprocess_step;

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
	  m_spBanachSpace(new IGridFunctionSpace<grid_function_type>() ),
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
		if (tdisc.num_stages() == 1) {
			return apply_single_stage(u1,t1,u0,t0);
		}else{ untested();
			return apply_multi_stage(u1,t1,u0,t0);
		}
	}

	void set_derivative(SmartPtr<grid_function_type> udot)
	{ m_spDerivative = udot; }

	SmartPtr<grid_function_type> get_derivative()
	{ return m_spDerivative; }

	number get_consistency_error() const
	{ return m_initial_consistency_error; }

	void set_banach_space(SmartPtr<IGridFunctionSpace<grid_function_type> > spSpace)
	{ m_spBanachSpace = spSpace; }

	void set_finished_tester(SmartPtr<FinishedTester> p) {
		_finished_tester = p;
	}

	void set_b_finish_time_step(bool x) {
		_bFinishTimeStep = x;
	}

	void set_output(SmartPtr<ITimeIntegratorObserver<TDomain, TAlgebra>> x) {
		_output = x;
	}
protected:
	bool apply_single_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
	bool apply_multi_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

	bool hasTerminated(double tCurrent, double tStart, double tFinal, unsigned step=0) const {
		if(_finished_tester){
			return _finished_tester->is_finished(tCurrent, step);
		}else{ untested();
			/*return (! ((tCurrent < tFinal) && (tFinal-tCurrent > base_type::m_precisionBound)));*/
			return tCurrent >= tFinal || tFinal-tCurrent < (tFinal-tStart)*base_type::m_precisionBound;
		}
	}

	/// metric
	SmartPtr<IGridFunctionSpace<grid_function_type> > m_spBanachSpace;

	SmartPtr<grid_function_type> m_spDerivative;

	number m_initial_consistency_error;

private: // to base class?
	mutable /*BUG*/ SmartPtr<FinishedTester> _finished_tester;
	bool _bFinishTimeStep{false};
	SmartPtr<ITimeIntegratorObserver<TDomain, TAlgebra>> _output;
}; // SimpleTimeIntegrator


template<typename TDomain, typename TAlgebra>
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_single_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	LIMEX_PROFILE_FUNC()

	// short-cuts
	assert(u0);
	GridLevel const &gl = u0->grid_level();
	assert(tdisc_dep_type::m_spTimeDisc);
	time_disc_type &tdisc = *tdisc_dep_type::m_spTimeDisc;
	assert(base_type::m_spSolver);
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

	if(!base_type::m_bNoLogOut) {
		UG_LOG("+++ Integrating: [\t"<< t0 <<"\t, \t"<< t1 <<"\t] with\t" << currdt <<"\n");
	}else{
	}

	while(!hasTerminated(t, t0, t1, step)) {
		if(!base_type::m_bNoLogOut){
			UG_LOG("++++++ TIMESTEP " << step << " BEGIN (current time: " << lua_like_string(t) << ") ++++++\n");
		}else{
		}

		// determine step size
		UG_COND_THROW(currdt < base_type::get_dt_min(),
				"Time step size " + std::to_string(currdt) + " below minimum. ABORTING!");
		number dt = std::min(currdt, t1-t);
		final_dt = dt;

		UG_LOG("++++++ Time step size: " + lua_like_string(dt) + "\n");

		bool initres = notify_init_step(u1, step, t, dt);
		if(initres){
		}else{ untested();
			incomplete();//?
			UG_LOG("++++++ Preparation of the time step failed.\n");
			break;
		}

		// foreach stage { ?

		//// from SolveNonlinearProblem
		bool pp_res = notify_preprocess_step(u1, step, t, dt);
		if(pp_res){
		}else{ untested();
			UG_LOG("++++++ PreProcess failed.\n");
			// newtonSuccess = false; incomplete.
			break;
		}
		/////

		// TODO: from SolveNonlinearProblem
		// if ((_endTime-(_time+currdt))/_maxStepSize < _relPrecisionBound){ untested()
		// }
		tdisc.set_stage(1);

		// prepare step
		tdisc.prepare_step(m_spSolTimeSeries, dt);
		if (solver.prepare(*u1) == false) { untested();
			if(!base_type::m_bNoLogOut)
				UG_LOG("Initialization failed! RETRY");

			currdt *= base_type::get_reduction_factor();
			continue;
		}else{
		}
		//UG_LOG("m_spSolTimeSeries.size="<< m_spSolTimeSeries->size());

		// .. much more (obsolete?) code in SolveNonlinearProblem
		// is this obsolete (and part of newton solver?)

		// execute step
		if (solver.apply(*u1)) {

			pp_res = notify_postprocess_step(u1, step, tdisc.future_time(), dt);
			if(pp_res){
			}else{
				incomplete();
				UG_LOG("++++++ PostProcess failed.\n");
				//						newtonSuccess = false;
				break;
			}
			//
			// ACCEPT step
			//

			// post prcess (e.g. physics)
			if(base_type::m_bNoLogOut){ untested();
			}else{
				// m_spSolTimeSeries->oldest() actually holds a pointer to a grid function
				// but as the time series does not know this, we have to cast ourselves
				// SmartPtr<grid_function_type> tmpOld = m_spSolTimeSeries->oldest().template cast_static<grid_function_type>();
#if 1 // this works with gkn77, test: total_nitrat.dat
				this->notify_finalize_step(u1, step, t, dt);
#else // this does not
				this->notify_finalize_step(u1, step, t+dt, dt);
#endif
			}

			// update time
			t += dt;

			// push updated solution into time series (and continue)
			//SmartPtr<typename base_type::vector_type> utmp = m_spSolTimeSeries->oldest();
			//VecAssign(*utmp, static_cast<typename base_type::vector_type> (*u1) );
			uold = m_spSolTimeSeries->push_discard_oldest(u1->clone(), t).template cast_static<grid_function_type>();
			}else{ untested();
				//
				// REJECT step
				//
				UG_LOG("Solution failed! RETRY");
				currdt *= base_type::get_reduction_factor();
				continue;
			}

		// consistency check
		if (step == 1 && m_spDerivative.valid())
		{ untested();
			UG_LOG("Computing consistency error: "<< std::endl);
			UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*u0),
					"Huhh: Different vectors required!");
			*m_spDerivative = *u0;
			m_initial_consistency_error = m_spBanachSpace->distance(*m_spDerivative, *u1);
		}else{
		}

		// } stage loop?


		if(_bFinishTimeStep){
			tdisc.finish_step_elem(m_spSolTimeSeries, dt);
		}else{
		}

		if(_output){
			_output->step_process(u1, step, tdisc.future_time(), 0.);
		}else{
		}

		UG_LOG("++++++ TIMESTEP " + std::to_string(step)
				+ " END   (current time: " + lua_like_string(t) + ") ++++++\n");

		step++;

	} // time loop

	if(base_type::m_bNoLogOut) { untested();
		this->notify_finalize_step(u1, /*uold,*/ step, t, final_dt);
	}else{
	}


	if (m_spDerivative.valid())
	{ untested();

		//
		// approximate derivative (by forward difference)
		//
		UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*uold),
		          "Huhh: Different vectors required!");

		VecScaleAdd(static_cast<typename TAlgebra::vector_type&>(*m_spDerivative),
				1.0/final_dt, static_cast<typename TAlgebra::vector_type&>(*u1),
				-1.0/final_dt, static_cast<typename TAlgebra::vector_type&>(*uold));

	}else{
	}

	m_spSolTimeSeries->clear();

	return true;

} // apply_single_stage

template<typename TDomain, typename TAlgebra>
bool SimpleTimeIntegrator<TDomain, TAlgebra>::apply_multi_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{ untested();

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

	if(!base_type::m_bNoLogOut){ untested();
		UG_LOG("+++ Integrating: ["<< t0 <<", "<< t1 <<"] with " << currdt <<"\n");
	}else{ untested();
	}

	// similar to SolveLinearTimeProblem in ugcore/.../time_step_util.lua
	// while not finishedTester:is_finished(time, step) do ...
	while(!hasTerminated(t, t0, t1))
	{ untested();
		if(!base_type::m_bNoLogOut)
			UG_LOG("++++++ TIMESTEP " << step++ << " BEGIN (current time: " << lua_like_string(t) << ") ++++++\n");

		// determine step size
		UG_COND_THROW(currdt < base_type::get_dt_min(), "Time step size below minimum. ABORTING!")
		number dt = std::min(currdt, t1-t);

		final_dt = dt;

		double told = t;

		const int num_stages = tdisc.num_stages();
		int s=1;
		do // 'for-each-stage' loop
		{ untested();

			// for (int s=1; s<=num_stages; ++s)
			if(!base_type::m_bNoLogOut)
				UG_LOG("+++ STAGE "<< s << " BEGIN +++\n");

			// set stage
			tdisc.set_stage(s);

			// prepare step
			tdisc.prepare_step(m_spSolTimeSeries, dt);
			if (solver.prepare(*u1) == false) { untested();
				break;
			}else{ untested();
			}

			// execute step
			if (!solver.apply(*u1)) { untested();
				break;
			}else{ untested();
			}

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
		{ untested();
			// REJECT time step
			if(!base_type::m_bNoLogOut)
				UG_LOG("Solution failed! RETRY");

			currdt *= this->get_reduction_factor();
			t = told;
			continue;
		}
		else
		{ untested();
			if(!base_type::m_bNoLogOut)
			{ untested();
				this->notify_finalize_step(u1, /*uold, */step, t, final_dt);
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
	{ untested();
		this->notify_finalize_step(u1, /*uold,*/ step, t, final_dt);
	}
	return true;
}

} // ug
