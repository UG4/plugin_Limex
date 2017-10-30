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


#include <string>

// bridge
#include "bridge/util.h"
#include "bridge/util_overloaded.h"
#include "bridge/util_domain_algebra_dependent.h"
// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time


// ug
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "common/math/math_vector_matrix/math_vector.h"


// plugin
#include "time_disc/time_integrator.hpp"
#include "time_disc/time_extrapolation.h"
#include "time_disc/linear_implicit_timestep.h"
#include "time_disc/limex_integrator.hpp"
#include "newton_limex.h"




// include for plugins
using namespace std;
using namespace ug::bridge;

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif

namespace ug{

/* Define LIB_LIMEX DebugID */
DebugID LIB_LIMEX("LIB_LIMEX");

namespace Limex{



/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra dependent parts
 * of the plugin. All Functions and Classes depending on both Domain and Algebra
 * are to be placed here when registering. The method is called for all
 * available Domain and Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain, typename TAlgebra>
static void DomainAlgebra(Registry& reg, string grp)
{
	//	some defines/typedefs
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();
	typedef MultiStepTimeDiscretization<TAlgebra> TMultiStepTimeDisc;
	typedef ITimeDiscretization<TAlgebra> TTimeDisc;
	typedef typename TAlgebra::vector_type TVector;
	typedef GridFunction<TDomain,TAlgebra> TGridFunction;

	{
		/// GridFunctionEstimator (REPLACED, sub-diagonal)
		typedef ISubDiagErrorEst<TVector> TBase;
		typedef GridFunctionEstimator<TDomain, TAlgebra> T;
		string name = string("GridFunctionEstimator").append(suffix);

		reg.add_class_<T, TBase>(name, grp)
					   .template add_constructor<void (*)() >("")
					   .template add_constructor<void (*)(number) >("")
					   .add_method("set_reference_norm", &T::set_reference_norm)
					   .add_method("add", &T::add)
					   .add_method("config_string", &T::config_string)
					   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "GridFunctionEstimator", tag);
	}



	{
		// SupErrorEvaluator
		typedef SupErrorEvaluator<TGridFunction> T;
		typedef IComponentSpace<TGridFunction> TBase;

		string name = string("SupErrorEvaluator").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)(const char *) >("fctNames")
		   .template add_constructor<void (*)(const char *, number) >("fctNames, scale")
		   .template add_constructor<void (*)(const char *, const char *, number) >("fctNames, subsetNames, scale")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SupErrorEvaluator", tag);
	}

	{
			// UserDataEvaluatorNumber
			typedef UserDataEvaluator<TGridFunction, number > T;
			typedef IComponentSpace<TGridFunction> TBase;

			string name = string("UserDataEvaluatorNumber").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
			   .template add_constructor<void (*)(const char *) >("fctNames")
			   .template add_constructor<void (*)(const char *, int) >("fctNames, scale")
			   .template add_constructor<void (*)(const char *, int, number) >("fctNames, subsetNames, scale")
			   .add_method("set_user_data", &T::set_user_data)
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "UserDataEvaluatorNumber", tag);
	}

	{
			// UserDataEvaluatorVector
			typedef UserDataEvaluator<TGridFunction, MathVector<TGridFunction::dim> > T;
			typedef IComponentSpace<TGridFunction> TBase;

			string name = string("UserDataEvaluatorVector").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.template add_constructor<void (*)(const char *) >("fctNames")
				.template add_constructor<void (*)(const char *, int) >("fctNames, scale")
				.template add_constructor<void (*)(const char *, int, number) >("fctNames, subsetNames, scale")
				.add_method("set_user_data", &T::set_user_data)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "UserDataEvaluatorVector", tag);
	}


	{
		// ScaledGridFunctionEstimator (sub-diagonal)
		typedef ISubDiagErrorEst<TVector> TBase;
		typedef ScaledGridFunctionEstimator<TDomain, TAlgebra> T;
		string name = string("ScaledGridFunctionEstimator").append(suffix);

		reg.add_class_<T, TBase>(name, grp)
		   .template add_constructor<void (*)() >("Default constructor")
		   .add_method("add", &T::add)
		   .add_method("config_string", &T::config_string)
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ScaledGridFunctionEstimator", tag);
	}


	{
		// ITimeIntegratorObserver (virtual base class)
		typedef ITimeIntegratorObserver<TDomain, TAlgebra> T;
		string name = string("ITimeIntegratorObserver").append(suffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, "ITimeIntegratorObserver", tag);
	}

	{
		// VTKOutputObserver
		typedef VTKOutputObserver<TDomain, TAlgebra> T;

		string name = string("VTKOutputObserver").append(suffix);
		reg.add_class_<T, typename T::base_type>(name, grp)
			.template add_constructor<void (*)(const char*, SmartPtr<typename T::vtk_type>) >("")
			.template add_constructor<void (*)(const char*, SmartPtr<typename T::vtk_type>, number) >("")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "VTKOutputObserver", tag);
	}

	{
		// ConnectionViewerOutputObserver
		typedef ConnectionViewerOutputObserver<TDomain, TAlgebra> T;

		string name = string("ConnectionViewerOutputObserver").append(suffix);
		reg.add_class_<T, typename T::base_type>(name, grp)
		   .template add_constructor<void (*)(const char*) >("")
		   .template add_constructor<void (*)(const char*, number) >("")
		   .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConnectionViewerOutputObserver", tag);
	}

	{
			// LuaOutputObserver
			typedef LuaCallbackObserver<TDomain, TAlgebra> T;

			string name = string("LuaCallbackObserver").append(suffix);
			reg.add_class_<T, typename T::base_type>(name, grp)
				.template add_constructor<void (*)() >("")
		 	 	.add_method("set_callback", &T::set_callback)
				.add_method("get_current_solution", &T::get_current_solution)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "LuaCallbackObserver", tag);
	}

	{
				// PlotRefOutputObserver
				typedef PlotRefOutputObserver<TDomain, TAlgebra> T;

				string name = string("PlotRefOutputObserver").append(suffix);
				reg.add_class_<T, typename T::base_type>(name, grp)
				    .template add_constructor<void (*)(const char*) >("")
					.template add_constructor<void (*)(const char*, SmartPtr<typename T::vtk_type>) >("")
					.set_construct_as_smart_pointer(true);
				reg.add_class_to_group(name, "PlotRefOutputObserver", tag);
		}


	{
			// IntegrationOutputObserver
			typedef IntegrationOutputObserver<TDomain, TAlgebra> T;

			string name = string("IntegrationOutputObserver").append(suffix);
			reg.add_class_<T, typename T::base_type>(name, grp)
				    .template add_constructor<void (*)() >("")
					 .add_method("add_integral_specs", &T::add_integral_specs)
					.set_construct_as_smart_pointer(true);
				reg.add_class_to_group(name, "IntegrationOutputObserver", tag);
		}

	{
		// ITimeIntegrator (virtual base class)
		typedef ITimeIntegrator<TDomain, TAlgebra> T;
		string name = string("ITimeIntegrator").append(suffix);
		reg.add_class_<T>(name, grp)
				  .add_method("set_time_step", &T::set_time_step)
				  .add_method("set_precision_bound", &T::set_precision_bound)
				  .add_method("set_no_log_out", &T::set_no_log_out)
				  .add_method("init", (void (T::*)(TGridFunction const&u) ) &T::init, "","")
				  .add_method("attach_observer", &T::attach_observer);
		reg.add_class_to_group(name, "ITimeIntegrator", tag);
	}

	{
		// ILinearTimeIntegrator
		typedef ITimeIntegrator<TDomain, TAlgebra> TBase;
		typedef ILinearTimeIntegrator<TDomain, TAlgebra> T;
		string name = string("ILinearTimeIntegrator").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
					  .add_method("set_linear_solver", &T::set_linear_solver);
		reg.add_class_to_group(name, "ILinearTimeIntegrator", tag);
	}

	{
		// LinearTimeIntegrator
		// (e.g., implicit Euler for linear problem)
		typedef ILinearTimeIntegrator<TDomain, TAlgebra> TBase;
		typedef LinearTimeIntegrator<TDomain, TAlgebra> T;
		//typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;

		string name = string("LinearTimeIntegrator").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
				  .template add_constructor<void (*)(SmartPtr<TTimeDisc>) >("")
				  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")
				  .add_method("get_time_disc", &T::get_time_disc)
				  .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearTimeIntegrator", tag);

	}

	{
		// LinearTimeIntegrator
		// (e.g., implicit Euler for linear problem)

		typedef ILinearTimeIntegrator<TDomain, TAlgebra> TBase;
		typedef ConstStepLinearTimeIntegrator<TDomain, TAlgebra> T;
		//typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;
		//typedef MultiStepTimeDiscretization<TAlgebra> TTimeDisc;

		string name = string("ConstStepLinearTimeIntegrator").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
					  .template add_constructor<void (*)(SmartPtr<TTimeDisc>) >("")
					  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")
					  .add_method("get_time_disc", &T::get_time_disc)
					  .add_method("set_num_steps", &T::set_num_steps)
					  .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstStepLinearTimeIntegrator", tag);

	}

	{
		// Adaptive LinearTimeIntegrator
		typedef ILinearTimeIntegrator<TDomain, TAlgebra> TBase;
		typedef TimeIntegratorLinearAdaptive<TDomain, TAlgebra> T;
		//typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;
		//typedef MultiStepTimeDiscretization<TAlgebra> TTimeDisc;

		string name = string("TimeIntegratorLinearAdaptive").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
				  .template add_constructor<void (*)(SmartPtr<TTimeDisc>, SmartPtr<TTimeDisc>) >("tdisc1, tdisc2")
				  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")
				  .add_method("get_time_disc", &T::get_time_disc)
				  .add_method("set_tol", &T::set_tol)
				  .add_method("set_time_step_max", &T::set_time_step_max)
				  .add_method("set_time_step_min", &T::set_time_step_min)
				  .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TimeIntegratorLinearAdaptive", tag);
	}

	{
			// INonlinearTimeIntegrator
			typedef ITimeIntegrator<TDomain, TAlgebra> TBase;
			typedef INonlinearTimeIntegrator<TDomain, TAlgebra> T;
			string name = string("INonlinearTimeIntegrator").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
						  .add_method("set_solver", &T::set_solver)
						  .add_method("set_dt_min", &T::set_dt_min)
						  .add_method("set_dt_max", &T::set_dt_max)
						  .add_method("set_reduction_factor", &T::set_reduction_factor)
						  .add_method("set_increase_factor", &T::set_increase_factor);
			reg.add_class_to_group(name, "INonlinearTimeIntegrator", tag);
	}

	/*{
		// INonlinearTimeIntegratorWithBounds

		typedef INonlinearTimeIntegrator<TDomain, TAlgebra> TBase;
		typedef INonlinearTimeIntegratorWithBounds<TDomain, TAlgebra> T;
		string name = string("ITimeIntegratorWithBounds").append(suffix);
		reg.add_class_<T, TBase>(name, grp)

		reg.add_class_to_group(name, "ITimeIntegratorWithBounds", tag);
	}*/

	{
		// SimpleTimeIntegrator
		// (e.g., implicit Euler for linear problem)
		typedef INonlinearTimeIntegrator<TDomain, TAlgebra> TBase;
		typedef SimpleTimeIntegrator<TDomain, TAlgebra> T;
		//typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;
		//typedef MultiStepTimeDiscretization<TAlgebra> TTimeDisc;

		string name = string("SimpleTimeIntegrator").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
							  .template add_constructor<void (*)(SmartPtr<TTimeDisc>) >("")
							  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")
							  .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "SimpleTimeIntegrator", tag);
	}

	{
			// LimexTimeIntegrator
			typedef INonlinearTimeIntegrator<TDomain, TAlgebra> TBase;
			typedef LimexTimeIntegrator<TDomain, TAlgebra> T;
			typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;

			string name = string("LimexTimeIntegrator").append(suffix);
			reg.add_class_<T,TBase>(name, grp)

			  //.template add_constructor<void (*)() >("")
			  //.ADD_CONSTRUCTOR( (SmartPtr<TDomainDisc>, int) ) ("Domain disc|number of steps (vector)")
			  .ADD_CONSTRUCTOR( (int) ) ("number of stages")
			  .add_method("set_tolerance", &T::set_tolerance)
			  .add_method("set_stepsize_safety_factor",&T::set_stepsize_safety_factor)
			  .add_method("set_stepsize_reduction_factor", &T::set_stepsize_reduction_factor)
			  .add_method("set_stepsize_greedy_order_factor", &T::set_stepsize_greedy_order_factor)
			  .add_method("add_error_estimator", &T::add_error_estimator)
			  .add_method("add_stage", (void (T::*)(size_t, size_t,  SmartPtr<typename T::domain_discretization_type>, SmartPtr<typename T::solver_type>) ) &T::add_stage)
			  .add_method("add_stage", (void (T::*)(size_t, SmartPtr<typename T::solver_type>, SmartPtr<typename T::domain_discretization_type>) ) &T::add_stage)
			  .add_method("add_stage", (void (T::*)(size_t, SmartPtr<typename T::solver_type>, SmartPtr<typename T::domain_discretization_type>, SmartPtr<typename T::domain_discretization_type>) ) &T::add_stage_ext)
			  .add_method("set_debug", &T::set_debug)
			  .add_method("set_debug_for_timestepper", &T::set_debug_for_timestepper)
			  .add_method("has_time_derivative", &T::has_time_derivative)
			  .add_method("get_time_derivative", &T::get_time_derivative)
			  .add_method("set_time_derivative", &T::set_time_derivative)
			  .add_method("enable_matrix_cache", &T::enable_matrix_cache)
			  .add_method("disable_matrix_cache", &T::disable_matrix_cache)
			 // .add_method("enable_linear_mode", &T::enable_linear_mode)
			 // .add_method("disable_linear_mode", &T::disable_linear_mode)
			  .add_method("select_cost_strategy", &T::select_cost_strategy)
			  .add_method("set_space", &T::set_space)
			  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")
			  .add_method("interrupt", &T::interrupt, "", "", "interrupt execution of apply()")
			  .set_construct_as_smart_pointer(true);

			reg.add_class_to_group(name, "LimexTimeIntegrator", tag);
	}

}

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

}

/**
 * Function called for the registration of Dimension dependent parts
 * of the plugin. All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
//	useful defines
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

}

/**
 * Function called for the registration of Algebra dependent parts
 * of the plugin. All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string parentGroup)
{

	//	typedefs for Vector and Matrix
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;

	//	useful defines
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

/*
	{
		// IErrorEvaluator (abstract base class)
		typedef IBanachSpace<vector_type> T;

		string name = string("IBanachSpace").append(suffix);
		reg.add_class_<T>(name, parentGroup)
							.add_method("norm", &T::norm)
							.add_method("distance", &T::distance);
		reg.add_class_to_group(name, "IBanachSpace", tag);
	}
*/
	//	LinearlyImplicitEuler
	{
		std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
		typedef ITimeDiscretization<TAlgebra> TBase;
		typedef LinearImplicitEuler<TAlgebra> T;
		string name = string("LinearImplicitEuler").append(suffix);
		reg.add_class_<T, TBase>(name, grp)
								.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("LinearImplicitEuler")
								.add_method("set_gamma_disc", &T::set_gamma_disc)
								.add_method("enable_matrix_cache", &T::enable_matrix_cache)
								.add_method("disable_matrix_cache", &T::disable_matrix_cache)
								.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearImplicitEuler", tag);
	}

	//	Time extrapolation
		{
			std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
			typedef typename TAlgebra::vector_type VT;
			typedef ISubDiagErrorEst<VT> ET;
			typedef AitkenNevilleTimex<VT> T;
			string name = string("AitkenNevilleTimex").append(suffix);
			reg.add_class_<T>(name, grp)
						.ADD_CONSTRUCTOR( (std::vector<size_t> nsteps) ) ("number of steps (vector)")
						.ADD_CONSTRUCTOR( (std::vector<size_t> nsteps, SmartPtr<ET> ) ) ("number of steps (vector)")
						.add_method("set_solution", &T::set_solution)
						.add_method("get_solution", &T::get_solution)
						.add_method("apply", (void (T::*)()) &T::apply)
						.add_method("apply", (void (T::*)(size_t, bool)) &T::apply)
						//.add_method("get_error_estimate", &T::get_error_estimate())
						//.add_method("get_error_estimate", &T::get_error_estimate(int))
						.add_method("get_error_estimate",  (double (T::*)(void)) &T::get_error_estimate, "","")
					    .add_method("get_error_estimate",  (double (T::*)(int)) &T::get_error_estimate, "","")
						.add_method("set_error_estimate", &T::set_error_estimate)
						.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "AitkenNevilleTimex", tag);
		}

		// interface for error estimator for time extrapolation
		{
				std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
				typedef ISubDiagErrorEst<vector_type> T;
				string name = string("ISubDiagErrorEst").append(suffix);
				reg.add_class_<T>(name, grp)
				   .add_method("config_string", &T::config_string);
				reg.add_class_to_group(name, "ISubDiagErrorEst", tag);
		}

		// L2 norm error estimator for time extrapolation
		{
			std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
			typedef Norm2Estimator<vector_type> T;
			typedef ISubDiagErrorEst<vector_type> TBase;
			string name = string("Norm2Estimator").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
			   .ADD_CONSTRUCTOR( (void) ) ("")
			   .add_method("set_stride", &T::set_stride)
			   .add_method("set_offset", &T::set_offset)
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "Norm2Estimator", tag);
		}

		// Inf norm error estimator for time extrapolation
		{
			std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
			typedef NormInfEstimator<vector_type> T;
			typedef ISubDiagErrorEst<vector_type> TBase;
			string name = string("NormInfEstimator").append(suffix);
			reg.add_class_<T,TBase>(name, grp)
			   .ADD_CONSTRUCTOR( (void) ) ("")
			   .add_method("set_stride", &T::set_stride)
			   .add_method("set_offset", &T::set_offset)
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NormInfEstimator", tag);
		}

		// Rel norm error estimator for time extrapolation
		{
			std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
			typedef NormRelEstimator<vector_type> T;
			typedef ISubDiagErrorEst<vector_type> TBase;
			string name = string("NormRelEstimator").append(suffix);
			reg.add_class_<T,TBase>(name, grp)
			   .ADD_CONSTRUCTOR( (void) ) ("")
			   .set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "NormRelEstimator", tag);
		}

		// LimexNewton
		{
			std::string grp = parentGroup;
			grp.append("/Discretization/Nonlinear");
			typedef LimexNewtonSolver<TAlgebra> T;
			typedef IOperatorInverse<vector_type> TBase;
			string name = string("LimexNewtonSolver").append(suffix);
			reg.add_class_<T, TBase>(name, grp)
				.add_constructor()
				.template add_constructor<void (*)(SmartPtr<IOperator<vector_type> >)>("Operator")
				.template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> >)>("AssemblingRoutine")
				.add_method("set_linear_solver", &T::set_linear_solver, "", "linear solver")
				.add_method("init", &T::init, "success", "op")
				.add_method("prepare", &T::prepare, "success", "u")
				.add_method("apply", &T::apply, "success", "u")
				.add_method("linear_solver_rate", &T::linear_solver_rate)
				.add_method("linear_solver_steps", &T::linear_solver_steps)
				.add_method("config_string", &T::config_string)
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "LimexNewtonSolver", tag);
		}

}

/**
 * Function called for the registration of Domain and Algebra independent parts
 * of the plugin. All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string parentGroup)
{

	std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
	typedef ILimexCostStrategy TBase;
	reg.add_class_<TBase>(std::string("ILimexCostStrategy"), grp);

	{
		typedef LimexDefaultCost T;
		string name = string("LimexDefaultCost");
		reg.add_class_<T, TBase>(name, grp)
		   .add_constructor()
		   .set_construct_as_smart_pointer(true);
	}

	{
			typedef LimexNonlinearCost T;
			string name = string("LimexNonlinearCost");
			reg.add_class_<T, TBase>(name, grp)
			   .add_constructor()
			   .set_construct_as_smart_pointer(true);
	}
}

}; // end Functionality

// end group sample_plugin
/// \}

}// end of namespace Sample


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_Limex(Registry* reg, string grp)
{
	grp.append("/Limex");
	typedef Limex::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(*reg,grp);
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
		RegisterAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

extern "C" UG_API void
FinalizeUGPlugin_Limex()
{
}

}//	end of namespace ug
