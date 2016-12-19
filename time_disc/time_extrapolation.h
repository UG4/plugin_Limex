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

#ifndef TIME_EXTRAPOLATION_H_
#define TIME_EXTRAPOLATION_H_

#include <vector>
#include <cmath>

// ug libraries
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_algebra/lib_algebra.h"
//#include "lib_algebra/parallelization/parallel_vector.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#include "lib_disc/function_spaces/integrate.h"

#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/linker/scale_add_linker.h"

namespace ug{

/*
 std::vector<size_t> steps {1, 2, 4};
 timex = new AitkenNevilleTimex(steps);

 //
 timex.set_solution(sol0, 0);  // single step
 timex.set_solution(sol1, 1);  // double step
 timex.set_solution(sol2, 2);  // four step

 timex.set_estimator(L2NormEstimator())
 timex.set_estimator(InfNormEstimator())
 timex.set_estimator(RelativeEstimator())
 //
 timex.apply()                // updating sol1 and sol2

 timex.get_error_estimate()

 */

namespace tools {


	//! calculates vUpdate = vUpdate + alpha2*vFine + alpha3*vCoarse. for doubles
	inline void VecScaleAddWithNormRel(double &vUpdate, double alpha2, const double &vFine, double alpha3, const double &vCoarse, double &norm)
	{
		const double update = alpha2*vFine + alpha3*vCoarse;
		vUpdate = vUpdate + update;
		norm = std::max(norm, 0.5*fabs(update)/(1.0+fabs(vFine)+fabs(vCoarse)));
	}

	//! calculates vUpdate = vUpdate + alpha2*vFine + alpha3*vCoarse
	template<typename vector_t, template <class T> class TE_VEC>
	inline void VecScaleAddWithNormRel(TE_VEC<vector_t> &vUpdate, double alpha2, const TE_VEC<vector_t> &vFine, double alpha3, const TE_VEC<vector_t> &vCoarse, double &norm)
	{
		for(size_t i=0; i<vUpdate.size(); i++)
			VecScaleAddWithNormRel(vUpdate[i], alpha2, vFine[i], alpha3, vCoarse[i], norm);
	}

	//! calculates vUpdate = vUpdate + alpha2*vFine + alpha3*vCoarse. for doubles
	inline void VecScaleAddWithNormInf(double &vUpdate, double alpha2, const double &vFine, double alpha3, const double &vCoarse, double &norm)
	{
		const double update = alpha2*vFine + alpha3*vCoarse;
		vUpdate = vUpdate + update;
		norm = std::max(norm, fabs(update));
	}

	//! calculates vUpdate = vUpdate + alpha2*vFine + alpha3*vCoarse
	template<typename vector_t, template <class T> class TE_VEC>
	inline void VecScaleAddWithNormInf(TE_VEC<vector_t> &vUpdate, double alpha2, const TE_VEC<vector_t> &vFine, double alpha3, const TE_VEC<vector_t> &vCoarse, double &norm, const int delta=1, const int offset=0)
	{
		// std::cerr << norm << " "<< delta << offset << std::endl;
		for(size_t i=offset; i<vUpdate.size(); i+=delta)
			VecScaleAddWithNormInf(vUpdate[i], alpha2, vFine[i], alpha3, vCoarse[i], norm);

		// std::cerr << norm << std::endl;
	}


	//! calculates vUpdate = vUpdate + alpha2*vFine + alpha3*vCoarse. for doubles
	inline void VecScaleAddWithNorm2(double &vUpdate, double alpha2, const double &vFine, double alpha3, const double &vCoarse, double &norm)
	{
		const double update = alpha2*vFine + alpha3*vCoarse;
		vUpdate = vUpdate+ update;
		norm += update*update;
	}

	//! calculates vUpdate = vUpdate + alpha2*vFine + alpha3*vCoarse
	template<typename vector_t, template <class T> class TE_VEC>
	inline void VecScaleAddWithNorm2(TE_VEC<vector_t> &vUpdate, double alpha2, const TE_VEC<vector_t> &vFine, double alpha3, const TE_VEC<vector_t> &vCoarse, double &norm, const int delta=1, const int offset=0)
	{
		for(size_t i=offset; i<vUpdate.size(); i+=delta)
			VecScaleAddWithNorm2(vUpdate[i], alpha2, vFine[i], alpha3, vCoarse[i], norm);
	}



}

/// Interface for error estimator
template <class TVector>
class ISubDiagErrorEst
{
public:
	// constructor
	ISubDiagErrorEst() : m_est(0.0) {};

	// vUpdateructor
	virtual ~ISubDiagErrorEst() {};

	// evaluate
	virtual bool update(SmartPtr<TVector> vUpdate, number alpha2, SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse) = 0;

	// get estimate
	number get_current_estimate() {return m_est; };
	void reset_estimate() {  m_est=0.0; };


protected:
	number m_est;
};


/// Evaluate using (algebraic) infinity norm
template <class TVector>
class NormInfEstimator : public ISubDiagErrorEst<TVector>
{
protected:
	typedef ISubDiagErrorEst<TVector> base_type;
	int m_stride;
	int m_offset;

public:
	// constructor
	NormInfEstimator() :
	ISubDiagErrorEst<TVector>(), m_stride(1), m_offset(0) {};

	// destructor
	//~NormInfEstimator() {}

	// apply
	bool update(SmartPtr<TVector> vUpdate, number alpha2, SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		const int delta = m_stride;
		const int offset = m_offset;
		base_type::m_est=0.0;
		tools::VecScaleAddWithNormInf(*vUpdate, alpha2, *vFine, -alpha2, *vCoarse, base_type::m_est, delta, offset);
		return true;
	}

	// offset (e.g., component to work on for systems w/ CPU1)
		void set_offset(int offset) {m_offset=offset;}

		// delta (e.g., total number components for systems w/ CPU1)
		void set_stride(int delta) {m_stride=delta;}

};


/// Evaluate using (algebraic) L2 norm
template <class TVector>
class Norm2Estimator : public ISubDiagErrorEst<TVector>
{
protected:
	typedef ISubDiagErrorEst<TVector> base_type;
	int m_stride;
	int m_offset;

public:
	// constructor
	Norm2Estimator() :
		ISubDiagErrorEst<TVector>(), m_stride(1), m_offset(0) {};
	Norm2Estimator(int stride) :
		ISubDiagErrorEst<TVector>(), m_stride(stride), m_offset(0) {};
	Norm2Estimator(int delta, int offset) :
		ISubDiagErrorEst<TVector>(), m_stride(delta), m_offset (offset)  {};

	// apply
	bool update(SmartPtr<TVector> vUpdate,  number alpha2,  SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		const int delta = m_stride;
		const int offset = m_offset;
		base_type::m_est=0.0;
		tools::VecScaleAddWithNorm2(*vUpdate, alpha2, *vFine, -alpha2, *vCoarse, base_type::m_est, delta, offset);
#ifdef UG_PARALLEL
		double locEst = base_type::m_est;
		double globEst =0.0;
		vUpdate->layouts()->proc_comm().allreduce(&locEst, &globEst, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
		base_type::m_est = globEst;
#endif
		base_type::m_est = sqrt(base_type::m_est);
		return true;
	}

	// offset (e.g., component to work on for systems w/ CPU1)
	void set_offset(int offset) { m_offset=offset;}

	// delta (e.g., total number components for systems w/ CPU1)
	void set_stride(int delta) { m_stride=delta;}

};


/// Evaluate using (algebraic) L2 norm
template <class TVector>
class NormRelEstimator : public ISubDiagErrorEst<TVector>
{
protected:
	typedef ISubDiagErrorEst<TVector> base_type;

public:
	// constructor
	NormRelEstimator() : ISubDiagErrorEst<TVector>() {};

	// apply w/ rel norm
	bool update(SmartPtr<TVector> vUpdate, number alpha2,  SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		base_type::m_est=0;
		tools::VecScaleAddWithNormRel(*vUpdate, alpha2, *vFine, -alpha2, *vCoarse, base_type::m_est);
		return true;
	}
};


/// Evaluate using (algebraic) L2 norm
template <class TDomain, class TAlgebra>
class GridFunctionEstimator : public ISubDiagErrorEst<typename TAlgebra::vector_type>
{
protected:
	typedef typename TAlgebra::vector_type TVector;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;

	int m_quadorder;
	number m_refNormValue;

	struct GridFunctionEvaluator
	{
		std::string m_fctNames;
		int m_quadorder;
		int m_type;
		double m_scale;

		GridFunctionEvaluator(const GridFunctionEvaluator &val)
		{
			this->m_fctNames = val.m_fctNames;
			this->m_quadorder = val.m_quadorder;
			this->m_type = val.m_type;
			this->m_scale = val.m_scale;
		}

		GridFunctionEvaluator(const char *fctNames) :
		m_fctNames(fctNames), m_quadorder(3), m_type(0), m_scale(1.0) {}

		GridFunctionEvaluator(const char *fctNames, int order) :
		m_fctNames(fctNames), m_quadorder(order), m_type(0), m_scale(1.0) {}

		GridFunctionEvaluator(const char *fctNames, int order, int type) :
		m_fctNames(fctNames), m_quadorder(order), m_type(type), m_scale(1.0) {}

		GridFunctionEvaluator(const char *fctNames, int order, int type, double scale) :
		m_fctNames(fctNames), m_quadorder(order), m_type(type), m_scale(scale) {}

		double compute_norm(SmartPtr<grid_function_type> uFine) const
		{
			return (m_type ==1) ? H1SemiNorm<grid_function_type>(uFine, m_fctNames.c_str(), m_quadorder)
			: L2Norm(uFine, m_fctNames.c_str(), m_quadorder);
		};

		double compute_error(SmartPtr<grid_function_type> uFine, SmartPtr<grid_function_type> uCoarse) const
		{
			return (m_type ==1) ? H1Error<grid_function_type>(uFine, m_fctNames.c_str(), uCoarse, m_fctNames.c_str() ,m_quadorder)
			: L2Error(uFine, m_fctNames.c_str(), uCoarse, m_fctNames.c_str() ,m_quadorder);
		};

	};

	std::vector<GridFunctionEvaluator> m_evaluators;

public:
	typedef ISubDiagErrorEst<TVector> base_type;

	// constructor
	GridFunctionEstimator(const char *fctNames) :
	ISubDiagErrorEst<TVector>(), m_refNormValue(0.0)
	{ this->add(fctNames); }

	GridFunctionEstimator(const char *fctNames, int order) :
	ISubDiagErrorEst<TVector>(), m_refNormValue(0.0)
	{ this->add(fctNames, order); }

	GridFunctionEstimator(const char *fctNames, int order, double ref) :
	ISubDiagErrorEst<TVector>(), m_refNormValue(ref)
	{ this->add(fctNames, order); }

	void add(const char *fctNames)
	{
		m_evaluators.push_back(GridFunctionEvaluator(fctNames));
	}

	void add(const char *fctNames, int order)
	{
		m_evaluators.push_back(GridFunctionEvaluator(fctNames, order));
	}

	void add4(const char *fctNames, int order, int type, double scale)
	{
		m_evaluators.push_back(GridFunctionEvaluator(fctNames, order, type, scale));
	}

	// apply w/ rel norm
	bool update(SmartPtr<TVector> vUpdate, number alpha,  SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		typedef ScaleAddLinker<number, TDomain::dim, number> linker_type;
		typedef GridFunction<TDomain, TAlgebra> grid_function_type;
		typedef GridFunctionNumberData<grid_function_type> TNumberData;

		// try upcast
		SmartPtr<grid_function_type> uFine = vFine.template cast_dynamic<grid_function_type>();
		SmartPtr<grid_function_type> uCoarse = vCoarse.template cast_dynamic<grid_function_type>();
		if (uFine.invalid() || uCoarse.invalid()) return false;


		// error estimate
		if (m_refNormValue<=0.0)
		{
			// relative error estimator
			//number unorm = L2Norm(uFine, m_fctNames.c_str(), m_quadorder);
			//number enorm = alpha*L2Error(uFine, m_fctNames.c_str(), uCoarse, m_fctNames.c_str() ,m_quadorder);
			number unorm = 0.0;
			number enorm = 0.0;
			for (typename std::vector<GridFunctionEvaluator>::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
			{
				unorm +=  it->compute_norm(uFine);
				enorm +=  alpha * it->compute_error(uFine, uCoarse);
			}

			base_type::m_est = enorm/unorm;

			std::cerr << "unorm=" << unorm << "enorm=" << enorm << "eps="<< base_type::m_est << std::endl;
		}
		else
		{
			// weighted error estimator
		//	number enorm = alpha*L2Error(uFine, m_fctNames.c_str(), uCoarse, m_fctNames.c_str() ,m_quadorder);
			number enorm = 0.0;
			for (typename std::vector<GridFunctionEvaluator>::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
			{

				enorm +=  alpha * it->compute_error(uFine, uCoarse);
			}
			base_type::m_est = enorm/m_refNormValue;

			std::cerr << "unorm (FIXED)=" << m_refNormValue << "enorm=" << enorm << "eps="<< base_type::m_est << std::endl;

		}

		// update
		VecScaleAdd(*vUpdate, 1.0+alpha, *vFine, -alpha, *vCoarse);
		return true;
	}

	void set_reference_norm(number norm)
	{m_refNormValue = norm; }
};


template <typename TVector>
class AitkenNevilleTimex
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:

		/** Aitken Neville scheme with a given number */
		AitkenNevilleTimex(std::vector<size_t> nsteps)
		: m_num_steps(nsteps),
		  m_stepsize(0.0),
		  m_subdiag (make_sp(new Norm2Estimator<TVector>())),
		  m_solution(nsteps.size()),
		  m_subdiag_error_est(nsteps.size(), INFINITY)
		{};

		AitkenNevilleTimex(std::vector<size_t> nsteps, SmartPtr<ISubDiagErrorEst<vector_type> > error)
		: m_num_steps(nsteps),
		  m_stepsize(0.0),
		  m_subdiag(error),
		  m_solution(nsteps.size()),
		  m_subdiag_error_est(nsteps.size(), INFINITY)
		{};

		virtual ~AitkenNevilleTimex() {}

		void set_global_stepsize(number H) {m_stepsize=H;}
		number get_global_stepsize() {return m_stepsize;}

		//! set solution (for stage i)
		 void set_solution(SmartPtr<vector_type> soli, int i)
		 { m_solution[i] = soli; }

		//! get solution (on stage i)
		SmartPtr<vector_type> get_solution(size_t i)
		{ return m_solution[i]; }

		//! set error estimator
		void set_error_estimate(SmartPtr<ISubDiagErrorEst<vector_type> > subdiag)
		{
			m_subdiag = subdiag;
		}


		const std::vector<number>& get_error_estimates() const
		{return m_subdiag_error_est; }

		//! error estimate on stage k
		number get_error_estimate(int k) const
		{ return m_subdiag_error_est[k];}


		//! best error estimate
		/*number get_error_estimate()
		{
			std::vector<number>::iterator best = std::min_element(m_subdiag_error_est.begin(), m_subdiag_error_est.end());
			return *best;
		}

		//! return best index
		int get_best_index() const
		{
			std::vector<number>::iterator best = std::min_element(m_subdiag_error_est.begin(), m_subdiag_error_est.end());
			//std::cout << "min element at: " << std::distance(std::begin(m_subdiag_error_est), best);
			return std::distance(m_subdiag_error_est.begin(), best);
		}

		 */

		/**
		 * Triangular Aitken Neville extrapolation:
		 *
		 * T_11
		 * T_21 T_22
		 *
		 * T_N1        T_NN
		 *
		 * for T_ik
 		 * */
		void apply(size_t nstages)
		{
			UG_ASSERT(nstages <= m_solution.size(),
					 "Dimensions do not match:"  << nstages << ">" << m_solution.size());

			// clear (for safety reasons...)
			for (size_t k=1; k<m_solution.size(); ++k)
			{ m_subdiag_error_est[k] = INFINITY; }

			//m_subdiag_error_est[0] = ;
			// process columns (left to right)
			for (size_t k=1; k<nstages; ++k)
			{

				// process rows (bottom up, allows recycling memory)
				for (size_t i=nstages-1; i>=k; --i)
				{
					UG_ASSERT(m_solution[i].valid(), "Invalid SmarPtr!");
					UG_ASSERT(m_solution[i-1].valid(), "Invalid SmarPtr!");

					SmartPtr<vector_type> solcoarse = m_solution[i-1];
					SmartPtr<vector_type> solfine = m_solution[i];

					// (2^p -1)
					// m_solution[i] += (1.0/scal)*(m_solution[i]- m_solution[i-1]);
					const number scaling = ((1.0*m_num_steps[i])/(1.0*m_num_steps[i-k])-1.0);
					UG_LOG("scaling="<<i << ","<< k <<
							": ns["<<i<<"]="<< m_num_steps[i] <<
							"ns2["<<i-k<<"]=" <<  m_num_steps[i-k] <<"=" << scaling << std::endl);

					if (i==k)
					{
						// compute subdiagonal error estimate
						m_subdiag->update(solfine, (1.0/scaling), solfine,  solcoarse);
						number subdiag_error_est=m_subdiag->get_current_estimate();

						//m_subdiag_error_est[k] = sqrt(subdiag_error_est*scaling);
						m_subdiag_error_est[k]=subdiag_error_est; //*scaling;

						UG_LOG(" ErrorEst["<< k<<"]=" << m_subdiag_error_est[k] << ";" << std::endl);
					}
					else
					{
						// standard case
						VecScaleAdd(*solfine, (1.0+1.0/scaling), *solfine,
										      -(1.0/scaling), *solcoarse);

					}
				} // rows i

			} // columns k

		}

		/// apply for all stages
		void apply()
		{
			//UG_ASSERT(m_num_steps.size() == m_solution.size(), "Dimensions do not match");
			const size_t nstages = m_num_steps.size();
			apply(nstages);
		}



protected:
		number substep(size_t i) {return m_stepsize/m_num_steps[i];}

private:
		/** time step */
		number m_stepsize;
		static const int m_order=1;

		/** evaluation on last stage */
		SmartPtr<ISubDiagErrorEst<vector_type> > m_subdiag;

		/** n_i: number of intermediate steps (per stage)*/
		std::vector<size_t> m_num_steps;

		/** vector of solutions (per stage)*/
		std::vector<SmartPtr<vector_type> > m_solution;

		/** sub-diagonal error estimate (per stage)*/
		std::vector<number> m_subdiag_error_est;


};



/** Estimate*/
/*
class TimeStepEstimator{

public:
	TimeStepEstimator(double tol) :
	m_rho(0.9), m_gamma(2.0), m_tol(tol)
	{

	}

	bool check(double eps, double &factor)
	{
		if (eps>=tol) return false;

		double val = (m_rho*m_tol)/eps;
		factor = pow(val, m_gamma);

	}


private:
	double m_rho;
	double m_gamma;
	double m_tol;
};
*/

}
#endif
