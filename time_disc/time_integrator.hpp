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

#ifndef TIME_INTEGRATOR_HPP_
#define TIME_INTEGRATOR_HPP_

#if __cplusplus >= 201103L
#define OVERRIDE override
#else
#define OVERRIDE
#endif

// std headers.
#include <string>

// UG4 headers.
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
#include "lib_disc/io/vtkoutput.h"



// My headers.
#include "time_extrapolation.h"
#include "../limex_tools.h"

namespace ug {


/// Sample class for integration observer: Output to VTK
template<class TDomain, class TAlgebra>
class VTKOutputObserver
: public ITimeIntegratorObserver<TDomain, TAlgebra>,
  public ITimeIntegratorStageObserver_start<TDomain, TAlgebra>,
  public ITimeIntegratorStageObserver_finalize<TDomain, TAlgebra>,
  public ITimeIntegratorStageObserver_end<TDomain, TAlgebra>
{
public:
	typedef ITimeIntegratorObserver<TDomain, TAlgebra> base_type;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef VTKOutput<TDomain::dim> vtk_type;

	VTKOutputObserver()
	:  m_sp_vtk(SPNULL), m_filename("0000"), m_startTime(0.0), m_plotStep(0.0) {}

	VTKOutputObserver(const char *filename, SmartPtr<vtk_type> vtk)
	: m_sp_vtk(vtk), m_filename(filename), m_startTime(0.0), m_plotStep(0.0) {}

	VTKOutputObserver(const char *filename, SmartPtr<vtk_type> vtk, number plotStep)
	: m_sp_vtk(vtk), m_filename(filename), m_startTime(0.0), m_plotStep(plotStep) {}

	virtual ~VTKOutputObserver()
	{ m_sp_vtk = SPNULL; }


	void set_output_scales(const std::vector<number>& vScales)
	{
		m_vOutputScales = vScales;
	}


	bool step_process(SmartPtr<grid_function_type> uNew, int step, number time, number dt) OVERRIDE
	{
		return true;
	}


	bool start_action(SmartPtr<grid_function_type> u, int step, number time, number dt) OVERRIDE
	{
		if (!m_sp_vtk.valid())
			return false;

		writeToFile(u, step, time);

		m_startTime = time;
		m_uOld = u->clone();

		return true;
	}


	bool finalize_action(SmartPtr<grid_function_type> uNew, int step, number time, number dt) OVERRIDE
	{
		if (!m_sp_vtk.valid())
			return false;

		if (m_plotStep == 0.0)
		{
			writeToFile(uNew, step, time);
			return true;
		}

		// otherwise, only plot data at multiples of given time step (interpolate if necessary)
		number rem = fmod(time - dt - m_startTime, m_plotStep);
		number nextPlotTimePt = time - dt - rem + m_plotStep;
		int nextStep = (int) ((nextPlotTimePt - m_startTime + 0.5 * m_plotStep) / m_plotStep);

		if (nextPlotTimePt > time)
		{
			m_uOld = uNew->clone();
			return true;
		}

		SmartPtr<grid_function_type> uCur = uNew->clone_without_values();
		while (nextPlotTimePt <= time)
		{
			number alpha = (time - nextPlotTimePt) / dt;
			VecScaleAdd(static_cast<typename TAlgebra::vector_type&>(*uCur),
				alpha, static_cast<typename TAlgebra::vector_type&>(*m_uOld),
				1.0 - alpha, static_cast<typename TAlgebra::vector_type&>(*uNew));

			writeToFile(uCur, nextStep, nextPlotTimePt);

			nextPlotTimePt = (++nextStep) * m_plotStep;
		}

		m_uOld = uNew->clone();
		return true;
	}


	bool end_action(SmartPtr<grid_function_type> u, int step, number time, number dt) OVERRIDE
	{
		if (!m_sp_vtk.valid())
			return false;

		m_sp_vtk->write_time_pvd(m_filename.c_str(), *u);
		return true;
	}


private:
	void writeToFile(SmartPtr<grid_function_type> u, int step, number time)
	{
		if (m_vOutputScales.size())
		{
			SmartPtr<grid_function_type> uTmp = u->clone_without_values();
			ScaleGF<grid_function_type>(uTmp, u, m_vOutputScales);
			m_sp_vtk->print(m_filename.c_str(), *uTmp, step, time);
		}
		else
			m_sp_vtk->print(m_filename.c_str(), *u, step, time);
	}

protected:
	SmartPtr<vtk_type> m_sp_vtk;
	SmartPtr<grid_function_type> m_uOld;
	std::string m_filename;
	number m_startTime;
	number m_plotStep;
	std::vector<number> m_vOutputScales;
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

	bool step_process(SmartPtr<grid_function_type> uNew, /*SmartPtr<grid_function_type> uOld, */int step, number time, number dt) OVERRIDE
	{
		// quit, if time does not match
		if (m_outputTime >=0.0 && time != m_outputTime) return true;

		SaveVectorForConnectionViewer<grid_function_type>(*uNew, m_filename.c_str());
		return true;
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
	bool step_process(SmartPtr<grid_function_type> uNew, /* SmartPtr<grid_function_type> uOld, */ int step, number time, number dt) OVERRIDE
	{
		UG_LOG("L2Error(\t"<< time << "\t) = \t" << L2Error(m_spReference, uNew, "c", time, 4) << std::endl);
		if (m_sp_vtk.valid())
		{
			SmartPtr<grid_function_type> ref = uNew->clone();
			Interpolate<grid_function_type> (m_spReference, ref, "c", time);
			m_sp_vtk->print("MyReference", *ref, step, time);
			return true;
		}

		return false;

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
	bool step_process(SmartPtr<grid_function_type> uNew, /*SmartPtr<grid_function_type> uOld,*/ int step, number time, number dt) OVERRIDE
	{

		for (typename std::vector<IntegralSpecs>::iterator it = m_vIntegralData.begin();
			 it!=m_vIntegralData.end(); ++it)
		{
			number value=Integral(uNew, it->m_cmp.c_str(), it->m_subsets.c_str(), it->m_quadOrder);
			UG_LOG("Integral(\t"<< it->m_idString << "\t"<< time << "\t)=\t" << value << std::endl);
		}

		return true;


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



/// integration of linear systems
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


} // ug

#include "simple_integrator.hpp"

namespace ug {

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

#endif /* TIME_INTEGRATOR_HPP_ */
