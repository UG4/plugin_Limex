// created by Arne Naegel

#include <string>

// bridge
#include "bridge/util.h"
#include "bridge/util_overloaded.h"
#include "bridge/util_domain_algebra_dependent.h"
// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time


// ug
#include "lib_disc/function_spaces/grid_function.h"
#include "common/ug_config.h"
#include "common/error.h"


// plugin
#include "time_disc/time_integrator.hpp"
#include "time_disc/time_extrapolation.h"
#include "time_disc/linear_implicit_timestep.h"


// include for plugins
using namespace std;
using namespace ug::bridge;

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif

namespace ug{
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
	//	useful defines
	string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
	string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();

	typedef ApproximationSpace<TDomain> TApproximationSpace;
	typedef GridFunction<TDomain,TAlgebra> TGridFunction;

	{
		// ITimeIntegrator
		typedef ITimeIntegrator<TDomain, TAlgebra> T;
		string name = string("ITimeIntegrator").append(suffix);
		reg.add_class_<T>(name, grp)
		  .add_method("set_time_step", &T::set_time_step)
		  .add_method("set_theta", &T::set_theta)
		  .add_method("init", (void (T::*)(TGridFunction const&u) ) &T::init, "","");
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
		typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;

		string name = string("LinearTimeIntegrator").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
		  .template add_constructor<void (*)(SmartPtr<TDomainDisc>) >("")
		  // .add_method("apply", &T::apply)
		  //.add_method("apply", (void (T::*)(TGridFunction &u, TGridFunction const &u0) ) &T::apply, "","")
		  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")
		  .add_method("get_time_disc", &T::get_time_disc)
		  .set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "LinearTimeIntegrator", tag);

	}

	{
			// LinearTimeIntegrator
			// (e.g., implicit Euler for linear problem)
			typedef ITimeIntegrator<TDomain, TAlgebra> TBase2;
			typedef ILinearTimeIntegrator<TDomain, TAlgebra> TBase;
			typedef ConstStepLinearTimeIntegrator<TDomain, TAlgebra> T;
			typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;

			string name = string("ConstStepLinearTimeIntegrator").append(suffix);
			reg.add_class_<T,TBase,TBase2>(name, grp)
			  .template add_constructor<void (*)(SmartPtr<TDomainDisc>) >("")
			  // .add_method("apply", &T::apply)
			  //.add_method("apply", (void (T::*)(TGridFunction &u, TGridFunction const &u0) ) &T::apply, "","")
			  //.add_method("apply", (void (T::*)(TGridFunction &u, number time, TGridFunction const &u0, number time0) ) &T::apply, "","")
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
		typedef DomainDiscretization<TDomain, TAlgebra> TDomainDisc;

		string name = string("TimeIntegratorLinearAdaptive").append(suffix);
		reg.add_class_<T,TBase>(name, grp)
		  .template add_constructor<void (*)(SmartPtr<TDomainDisc>) >("")
		  //.add_method("apply", (void (T::*)(TGridFunction &u, TGridFunction const &u0) ) &T::apply, "","")
		  //.add_method("apply", (void (T::*)(TGridFunction &u, number time, TGridFunction const &u0, number time0) ) &T::apply, "","")
		  .add_method("apply", (void (T::*)(SmartPtr<TGridFunction> u, number time, ConstSmartPtr<TGridFunction> u0, number time0) ) &T::apply, "","")						  
		  .add_method("get_time_disc", &T::get_time_disc)
		  .add_method("set_tol", &T::set_tol)
		  .add_method("set_time_step_max", &T::set_time_step_max)
		  .add_method("set_time_step_min", &T::set_time_step_min)
		.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "TimeIntegratorLinearAdaptive", tag);
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


	//	LinearImplicitEuler
		{
				std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
				typedef ITimeDiscretization<TAlgebra> TBase;
				typedef LinearImplicitEuler<TAlgebra> T;
				string name = string("LinearImplicitEuler").append(suffix);
				reg.add_class_<T, TBase>(name, grp)
						.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("LinearImplicitEuler")
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
						.add_method("apply", &T::apply)
						.add_method("get_error_estimate", &T::get_error_estimate)
						.add_method("set_error_estimate", &T::set_error_estimate)
						.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, "AitkenNevilleTimex", tag);
		}

		// interface for error estimator for time extrapolation
		{
				std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
				typedef ISubDiagErrorEst<vector_type> T;
				string name = string("ISubDiagErrorEst").append(suffix);
				reg.add_class_<T>(name, grp);
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



}

/**
 * Function called for the registration of Domain and Algebra independent parts
 * of the plugin. All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{

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
