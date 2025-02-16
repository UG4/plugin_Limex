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

#ifndef __H__LIMEX__AUX_OUTPUT_OBSERVER_HPP__
#define __H__LIMEX__AUX_OUTPUT_OBSERVER_HPP__

#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "lib_disc/io/vtkoutput.h"

namespace ug {

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

} // namespace ug

#endif /* __H__LIMEX__AUX_OUTPUT_OBSERVER_HPP__ */
