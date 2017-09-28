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
#include "lib_algebra/operator/debug_writer.h"

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

//! Interface for sub-diagonal error estimator (w.r.t time in Aitken-Neville scheme)
/*! Given u_fine and u_coarse,
 *
 * */
template <class TVector>
class ISubDiagErrorEst
{
public:
	// constructor
	ISubDiagErrorEst() : m_est(0.0) {};

	// destructor
	virtual ~ISubDiagErrorEst() {};

	//! compute update vUpdate = vFine + alpha * (vFine- vCoarse) AND estimate error \| alpha * (vFine- vCoarse)\|
	virtual bool update(SmartPtr<TVector> vUpdate, number alpha2, SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse) = 0;

	//! get estimate
	number get_current_estimate() {return m_est; };
	void reset_estimate() {  m_est=0.0; };

	virtual std::string config_string() const {return "ISubDiagErrorEst";}

protected:
	number m_est;
};


//! Evaluate using (algebraic) infinity norm
/*! Optional: provide component offset o and stride s, such that vector is evaluated at vec[o+k*s]
 * */
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


//! Evaluate using (algebraic) L2 norm
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

/*
template<class TVector>
class VectorDebugWritingEstimator
		: public VectorDebugWritingObject<TVector>
{
public:
///	type of vector
	typedef TVector vector_type;

public:
	VectorDebugWritingEstimator()
	: VectorDebugWritingObject<vector_type>() {}

	VectorDebugWritingEstimator(SmartPtr<IVectorDebugWriter<vector_type> > spDebugWriter)
	: VectorDebugWritingObject<vector_type>(spDebugWriter) {}

	int get_call_id() { return m_dgbCall; }
	void inc_call_id() { m_dgbCall++; }

protected:
	///	call counter
	int m_dgbCall;
};
*/



//! Estimate the error (based on the difference between two grid functions)
template <typename TGridFunction>
class IErrorEvaluator
{
protected:
	std::string m_fctNames;
	const char* m_ssNames;
	int m_quadorder;
	number m_scale;

public:
	IErrorEvaluator(const IErrorEvaluator<TGridFunction> &val)
	{
		this->m_fctNames = val.m_fctNames;
		this->m_ssNames = val.m_ssNames;
		this->m_quadorder = val.m_quadorder;
		this->m_scale = val.m_scale;
	}

	IErrorEvaluator(const char *fctNames) :
		m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(3), m_scale(1.0) {}

	IErrorEvaluator(const char *fctNames, int order) :
		m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(order), m_scale(1.0) {}

	IErrorEvaluator(const char *fctNames, int order, number scale) :
		m_fctNames(fctNames), m_ssNames(NULL), m_quadorder(order), m_scale(scale) {}

	IErrorEvaluator(const char *fctNames, const char* ssNames, int order, number scale) :
		m_fctNames(fctNames), m_ssNames(ssNames), m_quadorder(order), m_scale(scale) {}

	virtual ~IErrorEvaluator() {};

	virtual double compute_norm(SmartPtr<TGridFunction> u) const
	{ return 0.0; }

	virtual double compute_error(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{ return 0.0; }

	/// print config string
	std::string config_string() const
	{
		std::stringstream ss;

		if (this->m_ssNames)
		ss << this->m_fctNames << ", " << this->m_ssNames << ", " << this->m_quadorder
			<< ", type=" <<", scale=" << this->m_scale << std::endl;

		else
			ss << this->m_fctNames << ",  (no name), " << this->m_quadorder
						<< ", type=" <<", scale=" << this->m_scale << std::endl;
		return ss.str();
	}

};

/** Evaluates difference between two grid functions in L2 norm */
template <typename TGridFunction>
class L2ErrorEvaluator :
		public IErrorEvaluator<TGridFunction>
{
public:
	typedef IErrorEvaluator<TGridFunction> base_type;

	L2ErrorEvaluator(const char *fctNames) : base_type(fctNames) {};
	L2ErrorEvaluator(const char *fctNames, int order) : base_type(fctNames, order) {};
	L2ErrorEvaluator(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	L2ErrorEvaluator(const char *fctNames, const char* ssNames, int order, number scale)
		: base_type(fctNames, ssNames, order, scale) {};
	~L2ErrorEvaluator() {};

	double compute_norm(SmartPtr<TGridFunction> uFine) const
	{ return L2Norm(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, base_type::m_ssNames); }

	double compute_error(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{ return L2Error(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(),
		base_type::m_quadorder, base_type::m_ssNames);}
};


/** Evaluates distance between two grid functions in H1 semi-norm */
template <typename TGridFunction>
class H1SemiErrorEvaluator :
		public IErrorEvaluator<TGridFunction>
{
public:
	typedef IErrorEvaluator<TGridFunction> base_type;
	typedef typename H1SemiDistIntegrand<TGridFunction>::weight_type weight_type;

	H1SemiErrorEvaluator(const char *fctNames) : base_type(fctNames) {};
	H1SemiErrorEvaluator(const char *fctNames, int order) : base_type(fctNames, order) {};
	H1SemiErrorEvaluator(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	H1SemiErrorEvaluator(const char *fctNames, int order, number scale, SmartPtr<weight_type> spWeight)
	: base_type(fctNames, order, scale),  m_spWeight(spWeight) {};

	~H1SemiErrorEvaluator() {};

	double compute_norm(SmartPtr<TGridFunction> uFine) const
	{ return H1SemiNorm<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder, m_spWeight); }

	double compute_error(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{ return H1SemiDistance<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder, m_spWeight); }

	void set_weight(SmartPtr<weight_type> spWeight)
	{ m_spWeight = spWeight; }

protected:
	SmartPtr<weight_type> m_spWeight;
};

/** Evaluates difference between two grid functions in H1 semi-norm */
template <typename TGridFunction>
class H1ErrorEvaluator :
		public IErrorEvaluator<TGridFunction>
{
public:
	typedef IErrorEvaluator<TGridFunction> base_type;

	H1ErrorEvaluator(const char *fctNames) : base_type(fctNames) {};
	H1ErrorEvaluator(const char *fctNames, int order) : base_type(fctNames, order) {};
	H1ErrorEvaluator(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	~H1ErrorEvaluator() {};

	double compute_norm(SmartPtr<TGridFunction> uFine) const
	{ return H1Norm<TGridFunction>(uFine, base_type::m_fctNames.c_str(), base_type::m_quadorder); }

	double compute_error(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{ return H1Error<TGridFunction>(uFine, base_type::m_fctNames.c_str(), uCoarse, base_type::m_fctNames.c_str(), base_type::m_quadorder); }

};

/** Evaluates difference between two grid functions in L_inf norm */
template <typename TGridFunction>
class SupErrorEvaluator
: public IErrorEvaluator<TGridFunction>
{
	public:
		typedef IErrorEvaluator<TGridFunction> base_type;

		SupErrorEvaluator(const char *fctNames) : base_type(fctNames) {};
		SupErrorEvaluator(const char *fctNames, number scale) : base_type(fctNames, 1, scale) {};
		SupErrorEvaluator(const char *fctNames, const char* ssNames, number scale)
			: base_type(fctNames, ssNames, 1, scale) {};
		~SupErrorEvaluator() {};


		double compute_norm(SmartPtr<TGridFunction> uFine) const
		{
			// gather subsets in group
			SubsetGroup ssGrp(uFine->domain()->subset_handler());
			if (base_type::m_ssNames != NULL)
				ssGrp.add(TokenizeString(base_type::m_ssNames));
			else
				ssGrp.add_all();

			double maxVal = 0.0;

			// loop subsets
			for (size_t i = 0; i < ssGrp.size(); ++i)
			{
				// get subset index
				const int si = ssGrp[i];

				// loop elements of subset and dim
				switch (ssGrp.dim(i))
				{
					case DIM_SUBSET_EMPTY_GRID: break;
					case 0: maxVal = std::max(maxVal, findFctMaxOnSubset<Vertex>(uFine, si)); break;
					case 1: maxVal = std::max(maxVal, findFctMaxOnSubset<Edge>(uFine, si)); break;
					case 2: maxVal = std::max(maxVal, findFctMaxOnSubset<Face>(uFine, si)); break;
					case 3: maxVal = std::max(maxVal, findFctMaxOnSubset<Volume>(uFine, si)); break;
					default: UG_THROW("IntegrateSubsets: Dimension " << ssGrp.dim(i) << " not supported.");
				}
			}

			#ifdef UG_PARALLEL
				// max over processes
				if (pcl::NumProcs() > 1)
				{
					pcl::ProcessCommunicator com;
					number local = maxVal;
					com.allreduce(&local, &maxVal, 1, PCL_DT_DOUBLE, PCL_RO_MAX);
				}
			#endif

			// return the result
			return maxVal;
		}

		double compute_error(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
		{
			UG_COND_THROW(uFine->dof_distribution().get() != uCoarse->dof_distribution().get(),
				"Coarse and fine solutions do not have the same underlying dof distro.");

			SmartPtr<TGridFunction> uErr = uCoarse->clone();
			uErr->operator-=(*uFine);
			return compute_norm(uErr);
		}

	protected:
		template <typename TBaseElem>
		number findFctMaxOnSubset(ConstSmartPtr<TGridFunction> u, int si) const
		{
			ConstSmartPtr<DoFDistribution> dd = u->dof_distribution();

			size_t fct;
			try {fct = dd->fct_id_by_name(base_type::m_fctNames.c_str());}
			UG_CATCH_THROW("Function index could not be determined.\n"
				"Bear in mind that only one function can be evaluated in this error evaluator.");

			number maxVal = 0.0;
			typename DoFDistribution::traits<TBaseElem>::const_iterator it = dd->begin<TBaseElem>(si);
			typename DoFDistribution::traits<TBaseElem>::const_iterator itEnd = dd->end<TBaseElem>(si);
			for (; it != itEnd; ++it)
			{
				std::vector<DoFIndex> vInd;

				// we compare against all indices on the element
				// which means most indices will be compared against quite often
				// but as this is a sup norm, this is not a problem (only in terms of performance)
				size_t nInd = dd->dof_indices(*it, fct, vInd, false, false);
				for (size_t i = 0; i < nInd; ++i)
					maxVal = std::max(maxVal, fabs(DoFRef(*u, vInd[i])));
			}

			return maxVal;
		}
};



/*!
 * For arbitrary $\rho$, this class defines the integrand $|\rho(u_1)- \rho(u_2)|^2$.
 * If the grid function $u_2$ is not provided, the class returns $|\rho(u_1)|^2$
 * */
template <typename TDataIn, typename TGridFunction>
class DeltaSquareIntegrand
	: public StdIntegrand<number, TGridFunction::dim, DeltaSquareIntegrand<TDataIn, TGridFunction> >
{
	public:
	//	world dimension of grid function
		static const int worldDim = TGridFunction::dim;

	//	data type
		typedef TDataIn data_type;

	private:

		static number product(const number &x, const number &y)
		{return x*y;}

		static number product(const MathVector<worldDim> &x, const MathVector<worldDim> &y)
		{return VecDot(x,y);}


	//  data to integrate
		SmartPtr<UserData<TDataIn, worldDim> > m_spData;

	// 	grid function
		SmartPtr<TGridFunction> m_spGridFct;
		SmartPtr<TGridFunction> m_spGridFct2;

	//	time
		number m_time;

	public:
	/// constructor
		DeltaSquareIntegrand(SmartPtr<UserData<TDataIn, worldDim> > spData,
		                SmartPtr<TGridFunction> spGridFct,
						SmartPtr<TGridFunction> spGridFct2,
		                number time)
		: m_spData(spData), m_spGridFct(spGridFct), m_spGridFct2(spGridFct2), m_time(time)
		{
			m_spData->set_function_pattern(spGridFct->function_pattern());
		};

	/// constructor
		DeltaSquareIntegrand(SmartPtr<UserData<TDataIn, worldDim> > spData,
							  number time)
		: m_spData(spData), m_spGridFct(NULL), m_spGridFct2(NULL), m_time(time)
		{
			if(m_spData->requires_grid_fct())
				UG_THROW("UserDataDeltaIntegrand: Missing GridFunction, but "
						" data requires grid function.")
		};




		template <int elemDim>
		void get_values(TDataIn vValue[],
					  ConstSmartPtr<UserData<TDataIn, worldDim> > spData,
					  ConstSmartPtr<TGridFunction> spGridFct,
		              const MathVector<worldDim> vGlobIP[],
		              GridObject* pElem,
		              const MathVector<worldDim> vCornerCoords[],
		              const MathVector<elemDim> vLocIP[],
		              const MathMatrix<elemDim, worldDim> vJT[],
		              const size_t numIP)
		{
		//	get local solution if needed
			if(spData->requires_grid_fct())
			{
			//	create storage
				LocalIndices ind;
				LocalVector u;

			// 	get global indices
				spGridFct->indices(pElem, ind);

			// 	adapt local algebra
				u.resize(ind);

			// 	read local values of u
				GetLocalVector(u, *spGridFct);
				std::cout << u << std::endl;

			//	compute data
				try{
					(*spData)(vValue, vGlobIP, m_time, this->m_si, pElem,
								vCornerCoords, vLocIP, numIP, &u, &vJT[0]);
				}
				UG_CATCH_THROW("UserDataDeltaIntegrand: Cannot evaluate data.");
			}
			else
			{
			//	compute data
				try{
					(*spData)(vValue, vGlobIP, m_time, this->m_si, numIP);
				}
				UG_CATCH_THROW("UserDataDeltaIntegrand: Cannot evaluate data.");
			}

		};


		/// \copydoc IIntegrand::values
			template <int elemDim>
			void evaluate(number vValue[],
			              const MathVector<worldDim> vGlobIP[],
			              GridObject* pElem,
			              const MathVector<worldDim> vCornerCoords[],
			              const MathVector<elemDim> vLocIP[],
			              const MathMatrix<elemDim, worldDim> vJT[],
			              const size_t numIP)
			{
				std::vector<TDataIn> v1(numIP);

				get_values<elemDim>(&v1[0], m_spData, m_spGridFct, vGlobIP, pElem, vCornerCoords, vLocIP, vJT, numIP);
				std::cout << "--- got v1!" << std::endl;

				if (m_spGridFct2.valid())
				{
					std::vector<TDataIn> v2(numIP);
				/*	m_spGridFct->set(0.5);
					m_spGridFct2->set(0.5);*/
					get_values<elemDim>(&v2[0], m_spData, m_spGridFct2, vGlobIP, pElem, vCornerCoords, vLocIP, vJT, numIP);
					std::cout << "--- got v2!" << std::endl;

					for (size_t ip=0; ip<numIP; ++ip)
					{

						std::cout << std::setprecision(12) << v1[ip] <<" "<<  std::setprecision(12) << v2[ip] << std::endl;
						v1[ip] -= v2[ip];
						vValue[ip] = this->product(v1[ip], v1[ip]);

					}
				} else {

					for (size_t ip=0; ip<numIP; ++ip)
					{ vValue[ip] = this->product(v1[ip], v1[ip]); }

				}
			};
};

//! Evaluate the difference for a (dependent) UserData type induced by different grid functions
/*! UserData maybe of type TDataInput, i.e., number/vector/matrix/... */
template <typename TGridFunction, typename TDataInput>
class UserDataEvaluator :
		public IErrorEvaluator<TGridFunction>
{
public:
	typedef IErrorEvaluator<TGridFunction> base_type;
	//typedef MathVector<TGridFunction::dim> TDataInput;
	typedef UserData<TDataInput, TGridFunction::dim> input_user_data_type;

	UserDataEvaluator(const char *fctNames) : base_type(fctNames) {};
	UserDataEvaluator(const char *fctNames, int order) : base_type(fctNames, order) {};
	UserDataEvaluator(const char *fctNames, int order, number scale) : base_type(fctNames, order, scale) {};
	~UserDataEvaluator() {};

	void set_user_data(SmartPtr<input_user_data_type> spData)
	{
		m_userData = spData;
	}


	double compute_norm(SmartPtr<TGridFunction> uFine) const
	{
		SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= make_sp(new DeltaSquareIntegrand<TDataInput, TGridFunction> (m_userData, uFine, SPNULL, 0.0));

		return sqrt(IntegrateSubsets(spIntegrand, uFine, base_type::m_ssNames, base_type::m_quadorder, "best"));
	}

	double compute_error(SmartPtr<TGridFunction> uFine, SmartPtr<TGridFunction> uCoarse) const
	{
		SmartPtr<IIntegrand<number, TGridFunction::dim> > spIntegrand
		= make_sp(new DeltaSquareIntegrand<TDataInput, TGridFunction> (m_userData, uFine, uCoarse, 0.0));

		std::cout << "uFine="<<(void*) &(*uFine) << ", uCoarse="<< (void*)&(*uCoarse) << std::endl;

		return sqrt(IntegrateSubsets(spIntegrand, uFine, base_type::m_ssNames, base_type::m_quadorder, "best"));
	}

protected:
	SmartPtr<input_user_data_type> m_userData;
};


/*
/// Evaluate difference between two functions (w.r.t various norms)
template <class TDomain, class TAlgebra>
class GridFunctionEstimator :
		public ISubDiagErrorEst<typename TAlgebra::vector_type>
{
protected:
	typedef typename TAlgebra::vector_type TVector;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef IErrorEvaluator<grid_function_type> evaluator_type;

	std::vector<SmartPtr<evaluator_type> > m_evaluators;
	number m_refNormValue;

public:
	typedef ISubDiagErrorEst<TVector> base_type;

	// constructor
	ScaledGridFunctionEstimator() : m_refNormValue(0.0) {}
	ScaledGridFunctionEstimator(number ref) : m_refNormValue(ref) {}

	void add(SmartPtr<evaluator_type> eval)
	{
		m_evaluators.push_back(eval);
	}

	// apply w/ rel norm
	bool update(SmartPtr<TVector> vUpdate, number alpha,  SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		// typedef ScaleAddLinker<number, TDomain::dim, number> linker_type;
		typedef GridFunctionNumberData<TGridFunction> TNumberData;

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
			double est = 0.0;
			for (typename std::vector<evaluator_type>::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
			{
				enorm =  alpha * it->compute_error(uFine, uCoarse);
				unorm =  std::max(it->compute_norm(uFine), 1e-10);
				est += (enorm*enorm)/(unorm*unorm);
				std::cerr << "unorm=" << unorm << "enorm=" << enorm << "est="<<est << std::endl;
			}

			base_type::m_est = sqrt(est)/m_evaluators.size();
			std::cerr << "eps="<< base_type::m_est << std::endl;

		}
		else
		{
			// weighted error estimator
			number enorm = 0.0;
			for (typename std::vector<evaluator_type>::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
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


	/// print config string
	std::string config_string() const
	{
		std::stringstream ss;
		ss << "GridFunctionEstimator:\n";
		for (typename std::vector<GridFunctionEvaluator>::const_iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
		{
			ss << it->config_string();
		}
		return ss.str();
	}
};
*/

/// Evaluate using continuous norm (DEPRECATED!)
template <class TDomain, class TAlgebra>
class GridFunctionEstimator :
		public ISubDiagErrorEst<typename TAlgebra::vector_type>
{
protected:
	typedef typename TAlgebra::vector_type TVector;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef IErrorEvaluator<grid_function_type> evaluator_type;

	number m_refNormValue;

	std::vector<evaluator_type> m_evaluators;


public:
	typedef ISubDiagErrorEst<TVector> base_type;

	GridFunctionEstimator(const char *fctNames) :
		ISubDiagErrorEst<TVector>(), m_refNormValue(0.0)
	{
		this->add(fctNames);
	}

	GridFunctionEstimator(const char *fctNames, int order) :
		ISubDiagErrorEst<TVector>(), m_refNormValue(0.0)
	{
		this->add(fctNames, order);
	}

	GridFunctionEstimator(const char *fctNames, int order, double ref) :
		ISubDiagErrorEst<TVector>(), m_refNormValue(ref)
	{
		this->add(fctNames, order);
	}

	//! Add an L2 error for fctNames
	void add(const char *fctNames)
	{
		m_evaluators.push_back(L2ErrorEvaluator<grid_function_type>(fctNames));
	}

	void add(const char *fctNames, int order)
	{
		m_evaluators.push_back(L2ErrorEvaluator<grid_function_type>(fctNames, order));
	}

	void add4(const char *fctNames, int order, int type, double scale)
	{
		if (type == 0) {
			m_evaluators.push_back(L2ErrorEvaluator<grid_function_type>(fctNames, order, scale));
		} else if (type == 1) {
			m_evaluators.push_back(H1SemiErrorEvaluator<grid_function_type>(fctNames, order, scale));
		} else { UG_LOG(""); }
	}

	// apply w/ rel norm
	bool update(SmartPtr<TVector> vUpdate, number alpha,  SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		// typedefs
		typedef ScaleAddLinker<number, TDomain::dim, number> linker_type;
		typedef GridFunction<TDomain, TAlgebra> grid_function_type;

		//  upcast to GridFunction
		SmartPtr<grid_function_type> uFine = vFine.template cast_dynamic<grid_function_type>();
		SmartPtr<grid_function_type> uCoarse = vCoarse.template cast_dynamic<grid_function_type>();
		if (uFine.invalid() || uCoarse.invalid()) return false;

		// error estimate
		if (m_refNormValue<=0.0)
		{
			// relative error estimator
			number unorm = 0.0;
			number enorm = 0.0;
			for (typename std::vector<evaluator_type>::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
			{
				unorm +=  it->compute_norm(uFine);
				enorm +=  alpha * it->compute_error(uFine, uCoarse);
			}

			base_type::m_est = enorm/unorm;

			std::cerr << "unorm=" << unorm << "\tenorm=" << enorm << "\teps="<< base_type::m_est << std::endl;
		}
		else
		{
			// weighted error estimator
			number enorm = 0.0;
			for (typename std::vector<evaluator_type>::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
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


	/// print config string
	std::string config_string() const
	{
		std::stringstream ss;
		ss << "GridFunctionEstimator:\n";
		for (typename std::vector<evaluator_type>::const_iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
		{
			ss << it->config_string();
		}
		return ss.str();
	}
};

/// Evaluate difference between two functions (w.r.t various norms)
template <class TDomain, class TAlgebra>
class ScaledGridFunctionEstimator :
		public ISubDiagErrorEst<typename TAlgebra::vector_type>
{
protected:
	typedef typename TAlgebra::vector_type TVector;
	typedef GridFunction<TDomain, TAlgebra> grid_function_type;
	typedef IErrorEvaluator<grid_function_type> evaluator_type;

	std::vector<SmartPtr<evaluator_type> > m_evaluators;

public:
	typedef ISubDiagErrorEst<TVector> base_type;

	// constructor
	ScaledGridFunctionEstimator() : base_type() {}

	void add(SmartPtr<evaluator_type> eval)
	{ m_evaluators.push_back(eval); }

	// apply w/ rel norm
	bool update(SmartPtr<TVector> vUpdate, number alpha,  SmartPtr<TVector> vFine, SmartPtr<TVector> vCoarse)
	{
		// typedef GridFunctionNumberData<TGridFunction> TNumberData;

		// upcast
		SmartPtr<grid_function_type> uFine = vFine.template cast_dynamic<grid_function_type>();
		SmartPtr<grid_function_type> uCoarse = vCoarse.template cast_dynamic<grid_function_type>();
		if (uFine.invalid() || uCoarse.invalid()) return false;

		// error estimate
		number est = 0.0;
		for (typename std::vector<SmartPtr<evaluator_type> >::iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
		{
			// use sub-diagonal error estimator (i.e. multiply with alpha)
			double enorm =  alpha * (*it)->compute_error(uFine, uCoarse);
			double unorm =  std::max((*it)->compute_norm(uFine), 1e-10*enorm);
			est += (enorm*enorm)/(unorm*unorm);
			UG_LOGN("unorm=" << unorm << "\tenorm=" << enorm << "\tratio2="<< (enorm*enorm)/(unorm*unorm));
		}

		base_type::m_est = sqrt(est)/m_evaluators.size();
		UG_LOGN("eps="<< base_type::m_est);

		// update
		VecScaleAdd(*vUpdate, 1.0+alpha, *vFine, -alpha, *vCoarse);
		return true;
	}


	/// print config string
	std::string config_string() const
	{
		std::stringstream ss;
		ss << "ScaledGridFunctionEstimator:\n";
		for (typename std::vector<SmartPtr<evaluator_type> >::const_iterator it = m_evaluators.begin(); it!= m_evaluators.end(); ++it)
		{
			ss << (*it)->config_string();
		}
		return ss.str();
	}
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
		: m_stepsize(0.0),
		  m_subdiag (make_sp(new Norm2Estimator<TVector>())),
		  m_num_steps(nsteps),
		  m_solution(nsteps.size()),
		  m_subdiag_error_est(nsteps.size(), INFINITY)
		{};

		AitkenNevilleTimex(std::vector<size_t> nsteps, SmartPtr<ISubDiagErrorEst<vector_type> > error)
		: m_stepsize(0.0),
		  m_subdiag(error),
		  m_num_steps(nsteps),
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
		void apply(size_t nstages, bool with_error=true)
		{
			UG_ASSERT(nstages <= m_solution.size(),
					 "Dimensions do not match:"  << nstages << ">" << m_solution.size());

			if (with_error)
			{
				// reset (for safety reasons...)
				for (size_t k=1; k<m_solution.size(); ++k)
				{ m_subdiag_error_est[k] = INFINITY; }
			}

			//m_subdiag_error_est[0] = ;
			// process columns (left to right)
			for (size_t k=1; k<nstages; ++k)
			{

				// process rows (bottom-up -> recycling memory)
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

					if (with_error && (i==k))
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
