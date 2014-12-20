// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include "bridge/util.h"

// replace this with util_domain_dependent.h or util_algebra_dependent.h
// to speed up compilation time
#include "bridge/util_domain_algebra_dependent.h"
#include "bridge/util_overloaded.h"

#include "time_disc/time_extrapolation.h"
#include "time_disc/linear_implicit_timestep.h"

#include "common/ug_config.h"
#include "common/error.h"
#include <string>

using namespace std;
using namespace ug::bridge;

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif

namespace ug{
namespace Sample{

/** 
 *  \defgroup sample_plugin Sample
 *  \ingroup plugins_experimental
 *  This is just a sample plugin.
 *  \{
 */

void PluginSaysHello()
{
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
	cout << "Hello, I'm your plugin on proc " <<
				pcl::ProcRank() << "." << endl;
	pcl::SynchronizeProcesses();
#else
	UG_LOG("Hello, I'm your personal plugin in serial environment!\n");
#endif
}

void CrashFct(const string& reason)
{
	UG_THROW("I Crash because: "<< reason);
}

void CrashFctFatal(const string& reason)
{
	UG_THROW("I Crash fatal because: "<< reason);
}

void PluginCrashes()
{
	try{
		CrashFct("Some funny reason");
	}
	catch(bad_alloc& err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashes");
}

void PluginCrashesFatal()
{
	try{
		CrashFctFatal("Some fatal reason");
	}
	catch(bad_alloc& err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashesFatal");
}

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
						.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >)>("Domain Discretization")
					//	.template add_constructor<void (*)(SmartPtr<IDomainDiscretization<TAlgebra> >,int)>("Domain Discretization#Order")
						//.add_method("set_order", &T::set_order, "", "Order")
						//.add_method("set_stage", &T::set_stage, "", "Stage")
						.set_construct_as_smart_pointer(true);
				reg.add_class_to_group(name, "LinearImplicitEuler", tag);
		}

	//	Time extrapolation
		{
			std::string grp = parentGroup; grp.append("/Discretization/TimeDisc");
			typedef AitkenNevilleTimex<typename TAlgebra::vector_type> T;
			string name = string("AitkenNevilleTimex").append(suffix);
			reg.add_class_<T>(name, grp)
						.ADD_CONSTRUCTOR( (std::vector<size_t> nsteps) ) ("number of steps (vector)")
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
	reg.add_function("PluginSaysHello", &PluginSaysHello, grp)
		.add_function("PluginCrashes", &PluginCrashes, grp)
		.add_function("PluginCrashesFatal", &PluginCrashesFatal, grp);
}

}; // end Functionality

// end group sample_plugin
/// \}

}// end of namespace Sample


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_Sample(Registry* reg, string grp)
{
	grp.append("/Sample");
	typedef Sample::Functionality Functionality;

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
FinalizeUGPlugin_Sample()
{
}

}//	end of namespace ug
