/*
 * Copyright (c) 2014-2016:  G-CSC, Goethe University Frankfurt
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
#include "lib_disc/function_spaces/grid_function_util.h" // SaveVectorForConnectionViewer
#include "lib_disc/function_spaces/interpolate.h" //Interpolate
#include "lib_disc/function_spaces/integrate.h" //Integral
#include "lib_disc/io/vtkoutput.h"

// std headers
#include <string>

// own headers
#include "time_extrapolation.h"
#include "../limex_tools.h"

namespace ug {

/// Base class for time integration observer
template<class TDomain, class TAlgebra>
class ITimeIntegratorObserver
{
public:
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	virtual ~ITimeIntegratorObserver() {}
	//virtual void init(){}
	//virtual void finish(){}

	virtual void step_preprocess(SmartPtr<grid_function_type> u, int step, number time, number dt) {}
	virtual void step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt) {}

};

/// Sample class for integration observer: Output to VTK
template<class TDomain, class TAlgebra>
class VTKOutputObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef VTKOutput<TDomain::dim> vtk_type;

	VTKOutputObserver()
	:  m_sp_vtk(SPNULL), m_filename("0000"), m_plotStep(0.0) {}

	VTKOutputObserver(const char *filename, SmartPtr<vtk_type> vtk)
	: m_sp_vtk(vtk), m_filename(filename), m_plotStep(0.0) {}

	VTKOutputObserver(const char *filename, SmartPtr<vtk_type> vtk, number plotStep)
	: m_sp_vtk(vtk), m_filename(filename), m_plotStep(plotStep) {}

	virtual ~VTKOutputObserver()
	{ m_sp_vtk = SPNULL; }

	virtual void step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt)
	{
		if (!m_sp_vtk.valid())
			return;

		if (m_plotStep == 0.0)
		{
			m_sp_vtk->print(m_filename.c_str(), *uNew, step, time);
			return;
		}

		// otherwise, only plot data at multiples of given time step (interpolate if necessary)
		number rem = fmod(time-dt, m_plotStep);
		number curTime = time - dt - rem + m_plotStep;
		int curStep = (int) ((curTime + 0.5*m_plotStep) / m_plotStep);

		if (curTime > time)
			return;

		SmartPtr<grid_function_type> uCur = uNew->clone_without_values();
		while (curTime <= time)
		{
			number alpha = (time-curTime) / dt;
			VecScaleAdd(static_cast<typename TAlgebra::vector_type&>(*uCur),
				alpha, static_cast<typename TAlgebra::vector_type&>(*uOld),
				1.0 - alpha, static_cast<typename TAlgebra::vector_type&>(*uNew));

			m_sp_vtk->print(m_filename.c_str(), *uCur, curStep, curTime);

			curTime = (++curStep) * m_plotStep;
		}
	}

	void close(SmartPtr<grid_function_type> u)
	{
		if (m_sp_vtk.valid())
			m_sp_vtk->write_time_pvd(m_filename, u);
	}
protected:
	SmartPtr<vtk_type> m_sp_vtk;
	std::string m_filename;
	number m_plotStep;
};



/// Sample class for integration observer: Output to VTK
template<class TDomain, class TAlgebra>
class ConnectionViewerOutputObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	ConnectionViewerOutputObserver(const char *filename)
	: m_filename(filename), m_outputTime(-1.0) {}

	ConnectionViewerOutputObserver(const char *filename, number t_out)
	: m_filename(filename), m_outputTime(t_out) {}

	virtual ~ConnectionViewerOutputObserver()
	{}

	virtual void step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt)
	{
		// quit, if time does not match
		if (m_outputTime >=0.0 && time != m_outputTime) return;

		SaveVectorForConnectionViewer<grid_function_type>(*uNew, m_filename.c_str());
	}

protected:
	std::string m_filename;
	number m_outputTime;

};


#if 0
template <typename TData, typename TDataIn1, typename TDataIn2>
class LuaFunction2 // : public IFunction<TData, TDataIn1, typename TDataIn2>
{
	public:
	///	constructor
	LuaFunction2();
		virtual ~LuaFunction2() {};

	///	sets the Lua function used to compute the data
	/**
	 * This function sets the lua callback. The name of the function is
	 * passed as a string. Make sure, that the function name is defined
	 * when executing the script.
	 */
		void set_lua_callback(const char* luaCallback, size_t numArgs);

	///	evaluates the data
		virtual void operator() (TData& out, int numArgs1, int numArgs2,...);

	protected:
	///	callback name as string
		std::string m_cbValueName;

	///	reference to lua function
		int m_cbValueRef;

	///	lua state
		lua_State*	m_L;

	///	number of arguments to use
		size_t m_numArgs;
};


template <typename TData, typename TDataIn1, typename TDataIn2>
LuaFunction2<TData,TDataIn1,TDataIn2>::LuaFunction2() : m_numArgs(0)
{
	m_L = ug::script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
}

template <typename TData, typename TDataIn1, typename TDataIn2>
void LuaFunction2<TData,TDataIn1,TDataIn2>::set_lua_callback(const char* luaCallback, size_t numArgs)
{
//	store name (string) of callback
	m_cbValueName = luaCallback;

//	obtain a reference
	lua_getglobal(m_L, m_cbValueName.c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW("LuaFunction::set_lua_callback(...):"
				"Specified lua callback does not exist: " << m_cbValueName);
	}

//	store reference to lua function
	m_cbValueRef = luaL_ref(m_L, LUA_REGISTRYINDEX);

//	remember number of arguments to be used
	m_numArgs = numArgs;
}
#endif

/*
SmartUserDataWrapper* CreateNewUserData(lua_State* L, const SmartPtr<void>& ptr,
											  const char* metatableName)
{
//	create the userdata
	SmartUserDataWrapper* udata = (SmartUserDataWrapper*)lua_newuserdata(L,
											sizeof(SmartUserDataWrapper));
	new(udata) SmartUserDataWrapper;

//	associate the object with the userdata.
	udata->smartPtr = ptr;
	udata->type = SMART_POINTER;

//	associate the metatable (userdata is already on the stack)
	luaL_getmetatable(L, metatableName);
	lua_setmetatable(L, -2);

	return udata;
}
*/
/*
template <typename TData, typename TDataIn1, typename TDataIn2>
void LuaFunction2<TData,TDataIn1,TDataIn2>::operator() (TData& out, int numArgs1, SmartPtr<TDataIn1> valsArgs1[],
														int numArgs2, ...)
{
	PROFILE_CALLBACK_BEGIN(operatorBracket);

		UG_ASSERT((numArgs1+numArgs2) == (int)m_numArgs, "Number of arguments mismatched.");

	//	push the callback function on the stack
		lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_cbValueRef);

	//	get list of arguments
		va_list ap;
		va_start(ap, numArgs2);

	//	read all arguments and push them to the lua stack
		for(int i = 0; i < numArgs1; ++i)
		{

			CreateNewUserData(m_L, &valArgs1[i], "");

		}
		for(int i = 0; i < numArgs2; ++i)
		{
			TDataIn2 val = va_arg(ap, TDataIn2);
			lua_traits<TDataIn2>::push(m_L, val);
		}


	//	end read in of parameters
		va_end(ap);

	//	compute total args size
		size_t argSize = lua_traits<TDataIn1>::size * numArgs1;
		argSize += lua_traits<TDataIn2>::size * numArgs2;

	//	compute total return size
		size_t retSize = lua_traits<TData>::size;

	//	call lua function
		if(lua_pcall(m_L, argSize, retSize, 0) != 0)
			UG_THROW("LuaFunction::operator(...): Error while "
						"running callback '" << m_cbValueName << "',"
						" lua message: "<< lua_tostring(m_L, -1));

		try{
		//	read return value
			lua_traits<TData>::read(m_L, out);
			UG_COND_THROW(IsFiniteAndNotTooBig(out)==false, out);
		}
		UG_CATCH_THROW("LuaFunction::operator(...): Error while running "
						"callback '" << m_cbValueName << "'");

	//	pop values
		lua_pop(m_L, retSize);

	    PROFILE_CALLBACK_END();
}
*/
/// Integration observer: Output using Lua callback
/**!
 * TODO: Calls a LUA function with signature
 * func (SmartPtr<G> u, int step, number dt, number t)
 */
template<class TDomain, class TAlgebra>
class LuaCallbackObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef LuaFunction<number, number> lua_function_type;

	LuaCallbackObserver()
	{}

	virtual ~LuaCallbackObserver()
	{}

	// TODO: replace by call 'func (SmartPtr<G> u, int step, number dt, number t)'
	virtual void step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt)
	{
		number dummy_return;
		m_u = uNew;
		m_callback(dummy_return, numArgs2, (number) step, time, dt);

	}
	void set_callback(const char* luaCallback)
	{
		m_callback.set_lua_callback(luaCallback, numArgs2);
	};

	SmartPtr<grid_function_type> get_current_solution()
	{ return m_u; }

protected:
	// TODO: replace by appropriate call-back
	LuaFunction<number, number> m_callback;
	const static size_t numArgs1=0;  // num SmartPtr
	const static size_t numArgs2=3;
	SmartPtr<grid_function_type> m_u;
};


template<class TDomain, class TAlgebra>
class PlotRefOutputObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef LuaFunction<number, number> lua_function_type;
	typedef VTKOutput<TDomain::dim> vtk_type;

	PlotRefOutputObserver(SmartPtr<UserData<number, grid_function_type::dim> > spExactSol)
	{ m_spReference = spExactSol; }

#ifdef UG_FOR_LUA
	PlotRefOutputObserver(const char *ExactSol)
	: m_sp_vtk(SPNULL)
	{ m_spReference = make_sp(new LuaUserData<number, grid_function_type::dim>(ExactSol)); }

	PlotRefOutputObserver(const char *ExactSol, SmartPtr<vtk_type> vtk)
	: m_sp_vtk(vtk)
	{ m_spReference = make_sp(new LuaUserData<number, grid_function_type::dim>(ExactSol)); }

#endif

	virtual ~PlotRefOutputObserver()
	{}

	// TODO: replace by call 'func (SmartPtr<G> u, int step, number dt, number t)'
	virtual void step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt)
	{
		UG_LOG("L2Error(\t"<< time << "\t) = \t" << L2Error(m_spReference, uNew, "c", time, 4) << std::endl);
		if (m_sp_vtk.valid())
		{
			SmartPtr<grid_function_type> ref = uNew->clone();
			Interpolate<grid_function_type> (m_spReference, ref, "c", time);
			m_sp_vtk->print("MyReference", *ref, step, time);
		}



	}

protected:
	// TODO: replace by appropriate call-back
	SmartPtr<UserData<number, grid_function_type::dim> > m_spReference;
	SmartPtr<vtk_type> m_sp_vtk;
};

/// Integration observer: Output using Lua callback
/**!
 * TODO: should be replaced by LUA observer!
 */
template<class TDomain, class TAlgebra>
class IntegrationOutputObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef LuaFunction<number, number> lua_function_type;
	typedef VTKOutput<TDomain::dim> vtk_type;
protected:
	struct IntegralSpecs
	{
		IntegralSpecs(const char* cmp, const char* subsets, int quadOrder, const char *idString) :
			m_cmp(cmp), m_subsets(subsets), m_quadOrder(quadOrder), m_idString(idString)
		{};
		std::string m_cmp;
		std::string m_subsets;
		int m_quadOrder;
		std::string m_idString;
	};

public:
	IntegrationOutputObserver() : m_vIntegralData()
	{}

	virtual ~IntegrationOutputObserver()
	{}

	// TODO: replace by call 'func (SmartPtr<G> u, int step, number dt, number t)'
	virtual void step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt)
	{

		for (typename std::vector<IntegralSpecs>::iterator it = m_vIntegralData.begin();
			 it!=m_vIntegralData.end(); ++it)
		{
			number value=Integral(uNew, it->m_cmp.c_str(), it->m_subsets.c_str(), it->m_quadOrder);
			UG_LOG("Integral(\t"<< it->m_idString << "\t"<< time << "\t)=\t" << value << std::endl);
		}


	}


	void add_integral_specs(const char* cmp, const char* subsets, int quadOrder, const char* idString)
	{
		m_vIntegralData.push_back(IntegralSpecs(cmp, subsets, quadOrder, idString));
	}

protected:

	std::vector<IntegralSpecs> m_vIntegralData;
	//const char* cmp,
	//                const char* subsets,
	// int quadOrder
};



/// Base class for observer attachment
template<class TDomain, class TAlgebra>
class TimeIntegratorSubject
{
public:
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> process_observer_type;
	typedef typename std::vector<SmartPtr<process_observer_type> > process_observer_container_type;

protected:
	process_observer_container_type m_vProcessObservers;

public:
	//! register observer
	void attach_observer(SmartPtr<process_observer_type> obs)
	{m_vProcessObservers.push_back(obs);}

	//! notify all observers that time step evolution starts
	void notify_step_preprocess(SmartPtr<grid_function_type> u, int step, number time, number dt)
	{
		for (typename process_observer_container_type::iterator it = m_vProcessObservers.begin(); it!= m_vProcessObservers.end(); ++it)
		{(*it)->step_preprocess(u, step, time, dt); }
	}

	//! notify all observers that time step has been evolved (successfully)
	void notify_step_postprocess(SmartPtr<grid_function_type> uNew, SmartPtr<grid_function_type> uOld, int step, number time, number dt)
	{
		for (typename process_observer_container_type::iterator it = m_vProcessObservers.begin(); it!= m_vProcessObservers.end(); ++it)
		{(*it)->step_postprocess(uNew, uOld, step, time, dt); }
	}

};


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
		//SmartPtr<domain_disc_type> m_spDomainDisc;


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



/// integration of linear systems
template<class TDomain, class TAlgebra>
class ILinearTimeIntegrator : public ITimeIntegrator<TDomain, TAlgebra>
{

public:
	typedef ITimeIntegrator<TDomain, TAlgebra> base_type;
	typedef typename base_type::vector_type vector_type;
	typedef IPreconditionedLinearOperatorInverse<vector_type> linear_solver_type;
	typedef AssembledLinearOperator<TAlgebra> assembled_operator_type;

	// forward constructor
	ILinearTimeIntegrator()
	: base_type() {}

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

	//void init(grid_function_type const& u);
	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
};



template<typename TDomain, typename TAlgebra>
bool LinearTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{

	LIMEX_PROFILE_FUNC()

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
			UG_LOG("+++ Timestep +++" << step++ << "\n");

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
		}
		else
		{
			// DISCARDING
			currdt *= 0.5;
		}

	}

	return true;
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
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

private:

protected:

	int m_numSteps;
public:

	// constructor
	ConstStepLinearTimeIntegrator (SmartPtr<time_disc_type> tDisc)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc), m_numSteps(1) {}

	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
	void set_num_steps(int steps) {m_numSteps = steps;}
};



template<typename TDomain, typename TAlgebra>
bool ConstStepLinearTimeIntegrator<TDomain, TAlgebra>::apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0)
{
	LIMEX_PROFILE_FUNC()
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
		UG_LOG("+++ Integrating: [\t"<< t0 <<"\t, \t"<< t1 <<"\t] with dt=\t" << currdt << "("<< numSteps<< " iters)");
	 
	 // integrate
	 for(int step = 1; step<=numSteps; ++step)
	 {
	         // determine step size
	         // number dt = std::min(currdt, t1-t);
		 const number dt = currdt;

			if(!base_type::m_bNoLogOut)
				UG_LOG("+++ Const timestep +++" << step<< "t=" << t << "->" << dt);// std::endl;
		
		 // prepare step
		 tdisc.prepare_step(m_spSolTimeSeries, dt);
		 if (spAssOp==SPNULL)
		 {
			 // assemble operator
			if(!base_type::m_bNoLogOut)
			 UG_LOG("+++ Assemble (t=" << t << ", dt=" << dt <<")\n");

			 spAssOp=make_sp(new typename base_type::assembled_operator_type(tdisc_dep_type::m_spTimeDisc, gl));
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
			 VecAssign(*tmp,  *u1.template cast_dynamic<typename base_type::vector_type>());
			 m_spSolTimeSeries->push_discard_oldest(tmp, t);
		 }
		 else
		 {
			UG_THROW("Const Time step failed!!!")
		 }

	 }

	 return true;
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

private:


protected:

	SmartPtr<time_disc_type> m_spTimeDisc2;        // set during init

  double m_tol, m_dtmin, m_dtmax;


public:
  TimeIntegratorLinearAdaptive (SmartPtr<time_disc_type> tDisc1, SmartPtr<time_disc_type> tDisc2)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc1), m_tol(1e-2), m_dtmin(-1.0), m_dtmax(-1.0)
	{
		m_spTimeDisc2 = tDisc2;
	}

	void init(grid_function_type const& u);
	bool apply(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

    void set_tol(double tol) {m_tol = tol;}
    void set_time_step_min(number dt) {m_dtmin = dt;}
    void set_time_step_max(number dt) {m_dtmax = dt;}
};


template<typename TDomain, typename TAlgebra>
void TimeIntegratorLinearAdaptive<TDomain, TAlgebra>::init(grid_function_type const& u)
{
	// call base
	base_type::init(u);
	//m_spTimeDisc2 = make_sp(new typename base_type::time_disc_type(base_type::m_spDomainDisc, base_type::m_theta));
}

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
	if (m_dtmin<0.0) m_dtmin = (t1-10.0)/1.0e+5; 
	if (m_dtmax<0.0) m_dtmax = (t1-10.0)/10.0; 

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

/// integration of non-linear systems (with bounds on dt)
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




/// Integrate (a non-linear problem) over a given time interval
template<class TDomain, class TAlgebra>
class SimpleTimeIntegrator :
		public INonlinearTimeIntegrator<TDomain, TAlgebra>,
		public ITimeDiscDependentObject<TAlgebra>
{
protected:
	typedef ITimeDiscDependentObject<TAlgebra> tdisc_dep_type;


public:
	typedef INonlinearTimeIntegrator<TDomain, TAlgebra> base_type;
	typedef ITimeDiscretization<TAlgebra> time_disc_type;
	typedef typename base_type::grid_function_type grid_function_type;
	typedef VectorTimeSeries<typename base_type::vector_type> vector_time_series_type;

	// constructor
	SimpleTimeIntegrator (SmartPtr<time_disc_type> tDisc)
	: base_type(), ITimeDiscDependentObject<TAlgebra>(tDisc), m_spDerivative(SPNULL), m_initial_consistency_error(0.0)
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

protected:
	bool apply_single_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);
	bool apply_multi_stage(SmartPtr<grid_function_type> u1, number t1, ConstSmartPtr<grid_function_type> u0, number t0);

	inline bool hasTerminated(double tCurrent, double tStart, double tFinal) const
	{
	 	/*return (! ((tCurrent < tFinal) && (tFinal-tCurrent > base_type::m_precisionBound)));*/
		return ((tCurrent >= tFinal) ||
				(tFinal-tCurrent < (tFinal-tStart)*base_type::m_precisionBound + base_type::m_precisionBound));
	
	}

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

	 while(!hasTerminated(t, t0, t1))
	 {
		if(!base_type::m_bNoLogOut)
			 UG_LOG("+++ Timestep +++" << step << "\n");

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

		 UG_LOG("m_spSolTimeSeries.size="<< m_spSolTimeSeries->size());
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
					SmartPtr<grid_function_type> tmpOld = m_spSolTimeSeries->oldest().template cast_static<grid_function_type>();
					this->notify_step_postprocess(u1, tmpOld, step, t, dt);
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
		 UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*u0),
		 					 				  "Huhh: Different vectors required!");
		 		VecScaleAdd((typename base_type::vector_type&) *m_spDerivative, 1.0, *u1, -1.0, *u0);
		 		m_initial_consistency_error = m_spDerivative->norm();

		}

		 step++;
		// tdisc.finish_step_elem(m_spSolTimeSeries, dt);

	 }

	if(base_type::m_bNoLogOut)
	{
		this->notify_step_postprocess(u1, uold, step, t, final_dt);
	}


	if (m_spDerivative.valid())
	{
		/*
		UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*u0),
					 				  "Huhh: Different vectors required!");
		VecScaleAdd((typename base_type::vector_type&) *m_spDerivative, 1.0, *u1, -1.0, *u0);
		m_initial_consistency_error = m_spDerivative->norm();
*/

		//
		// approximate derivative (by forward difference)
		//
		UG_ASSERT(static_cast<typename base_type::vector_type*> (&*u1) != &(*uold),
				  "Huhh: Different vectors required!");
		VecScaleAdd((typename base_type::vector_type&) *m_spDerivative, 1.0/final_dt, *u1, -1.0/final_dt, *uold);
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
				this->notify_step_postprocess(u1, uold, step, t, dt);
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
		this->notify_step_postprocess(u1, uold, step, t, final_dt);
	}
	return true;
};


} // namespace ug

#endif /* TIME_INTEGRATOR_HPP_ */
