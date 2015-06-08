// external libraries

#ifndef TIME_EXTRAPOLATION_H_
#define TIME_EXTRAPOLATION_H_

#include <vector>
#include <cmath>

// ug libraries
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_algebra/lib_algebra.h"

#include "lib_disc/time_disc/time_disc_interface.h"

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


	//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3. for doubles
	inline void VecScaleAddWithNormRel(double &dest, double alpha1, const double &v1, double alpha2, const double &v2, double alpha3, const double &v3, double &norm)
	{
		const double update = alpha2*v2 + alpha3*v3;
		dest = alpha1*v1 + update;
		norm = std::max(norm, 0.5*fabs(update)/(1.0+fabs(v2)+fabs(v3)));
	}

	//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
	template<typename vector_t, template <class T> class TE_VEC>
	inline void VecScaleAddWithNormRel(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1, double alpha2, const TE_VEC<vector_t> &v2, double alpha3, const TE_VEC<vector_t> &v3, double &norm)
	{
		for(size_t i=0; i<dest.size(); i++)
			VecScaleAddWithNormRel(dest[i], alpha1, v1[i], alpha2, v2[i], alpha3, v3[i], norm);
	}

	//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3. for doubles
	inline void VecScaleAddWithNormInf(double &dest, double alpha1, const double &v1, double alpha2, const double &v2, double alpha3, const double &v3, double &norm)
	{
		const double update = alpha2*v2 + alpha3*v3;
		dest = alpha1*v1 + update;
		norm = std::max(norm, fabs(update));
	}

	//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
	template<typename vector_t, template <class T> class TE_VEC>
	inline void VecScaleAddWithNormInf(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1, double alpha2, const TE_VEC<vector_t> &v2, double alpha3, const TE_VEC<vector_t> &v3, double &norm, const int delta=1, const int offset=0)
	{
		std::cerr << norm << " "<< delta << offset << std::endl;
		for(size_t i=offset; i<dest.size(); i+=delta)
			VecScaleAddWithNormInf(dest[i], alpha1, v1[i], alpha2, v2[i], alpha3, v3[i], norm);

		std::cerr << norm << std::endl;
	}


	//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3. for doubles
	inline void VecScaleAddWithNorm2(double &dest, double alpha1, const double &v1, double alpha2, const double &v2, double alpha3, const double &v3, double &norm)
	{
		const double update = alpha2*v2 + alpha3*v3;
		dest = alpha1*v1 + update;
		norm += update*update;
	}

	//! calculates dest = alpha1*v1 + alpha2*v2 + alpha3*v3
	template<typename vector_t, template <class T> class TE_VEC>
	inline void VecScaleAddWithNorm2(TE_VEC<vector_t> &dest, double alpha1, const TE_VEC<vector_t> &v1, double alpha2, const TE_VEC<vector_t> &v2, double alpha3, const TE_VEC<vector_t> &v3, double &norm, const int delta=1, const int offset=0)
	{
		for(size_t i=offset; i<dest.size(); i+=delta)
			VecScaleAddWithNorm2(dest[i], alpha1, v1[i], alpha2, v2[i], alpha3, v3[i], norm);
	}



}

/// Interface for error estimator
template <class TVector>
class ISubDiagErrorEst
{
public:
	// constructor
	ISubDiagErrorEst() : m_est(0.0) {};

	// destructor
	virtual ~ISubDiagErrorEst() {};

	// evaluate
	virtual bool update(TVector &dest, number alpha1, TVector &v1, number alpha2, TVector &v2, number alpha3, TVector &v3) = 0;

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
	bool update(TVector &dest, number alpha1, TVector &v1, number alpha2, TVector &v2, number alpha3, TVector &v3)
	{
		const int delta = m_stride;
		const int offset = m_offset;
		base_type::m_est=0.0;
		tools::VecScaleAddWithNormInf(dest, alpha1, v1, alpha2, v2, alpha3, v3, base_type::m_est, delta, offset);
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
	bool update(TVector &dest, number alpha1, TVector &v1, number alpha2, TVector &v2, number alpha3, TVector &v3)
	{
		const int delta = m_stride;
		const int offset = m_offset;
		base_type::m_est=0.0;
		tools::VecScaleAddWithNorm2(dest, alpha1, v1, alpha2, v2, alpha3, v3, base_type::m_est, delta, offset);
#ifdef UG_PARALLEL
		double locEst = base_type::m_est;
		double globEst =0.0;
		dest.layouts()->proc_comm().allreduce(&locEst, &globEst, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
		base_type::m_est = globEst;
#endif
		base_type::m_est = sqrt(base_type::m_est);
		return true;
	}

	// offset (e.g., component to work on for systems w/ CPU1)
	void set_offset(int offset) {m_offset=offset;}

	// delta (e.g., total number components for systems w/ CPU1)
	void set_stride(int delta) {m_stride=delta;}

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
	bool update(TVector &dest, number alpha1, TVector &v1, number alpha2, TVector &v2, number alpha3, TVector &v3)
	{
		base_type::m_est=0;
		tools::VecScaleAddWithNormRel(dest, alpha1, v1, alpha2, v2, alpha3, v3, base_type::m_est);
		return true;
	}
};


/// Evaluate using (algebraic) L2 norm
template <class TVector>
class GridFunctionNormEstimator : public ISubDiagErrorEst<TVector>
{
protected:
	typedef ISubDiagErrorEst<TVector> base_type;

public:
	// constructor
	GridFunctionNormEstimator() : ISubDiagErrorEst<TVector>() {};

	// apply w/ rel norm
	bool update(TVector &dest, number alpha1, TVector &v1, number alpha2, TVector &v2, number alpha3, TVector &v3)
	{

		UG_THROW("Please implement me!");
		SmartPtr<TVector> delta = v2.clone();
		delta = alpha2*v2 + alpha3*v3;		// now do something with delta

		dest = alpha1*v1+ delta;
		return false;
	}
};


template <typename TVector>
class AitkenNevilleTimex
{
	public:
	///	vector type of solutions
		typedef TVector vector_type;

	public:

		/** Aitken Neville scheme with a given umber */
		AitkenNevilleTimex(std::vector<size_t> nsteps)
		: m_num_steps(nsteps),
		  m_stepsize(0.0),
		  m_solution(nsteps.size()),
		  m_subdiag_error_est(nsteps.size(), INFINITY),
		  m_subdiag (make_sp(new Norm2Estimator<TVector>())){};

		AitkenNevilleTimex(std::vector<size_t> nsteps, SmartPtr<ISubDiagErrorEst<vector_type> > error)
		: m_num_steps(nsteps),
		  m_stepsize(0.0),
		  m_solution(nsteps.size()),
		  m_subdiag_error_est(nsteps.size(), INFINITY),
		  m_subdiag(error){};

		virtual ~AitkenNevilleTimex() {}

		void set_global_stepsize(number H) {m_stepsize=H;}
		number get_global_stepsize() {return m_stepsize;}


		 void set_solution(SmartPtr<vector_type> soli, int i)
		 { m_solution[i] = soli; }

		SmartPtr<vector_type> get_solution(size_t i)
		{ return m_solution[i]; }

		number get_error_estimate()
		{ return m_subdiag_error_est.back();}

		void set_error_estimate(SmartPtr<ISubDiagErrorEst<vector_type> > subdiag)
		{ m_subdiag = subdiag; }

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
		void apply()
		{
			UG_ASSERT(  m_num_steps.size() == m_solution.size(), "Dimensions do not match");
			const size_t N = m_num_steps.size();

			//m_subdiag_error_est[0] = ;
			// process columns (left to right)
			for (size_t k=1; k<N; ++k)
			{

				// process rows (bottom up, allows recycling memory)
				for (size_t i=N-1; i>k; --i)
				{
					UG_ASSERT(m_solution[i].valid(), "Invalid SmarPtr!");
					UG_ASSERT(m_solution[i-1].valid(), "Invalid SmarPtr!");

					vector_type &solcoarse = *m_solution[i-1];
					vector_type &solfine = *m_solution[i];

					// (2^p -1)
					// m_solution[i] += (1.0/scal)*(m_solution[i]- m_solution[i-1]);
					const number scaling = (m_num_steps[i]/m_num_steps[i-k]-1.0);
					VecScaleAdd(solfine, (1.0+1.0/scaling), solfine,
										-(1.0/scaling), solcoarse);

				}

				// subdiagonal error estimate
				{
					UG_ASSERT(m_solution[k].valid(), "Invalid SmarPtr!");
					UG_ASSERT(m_solution[k-1].valid(), "Invalid SmarPtr!");

					vector_type &solcoarse = *m_solution[k-1];
					vector_type &solfine = *m_solution[k];


					const number scaling = ((1.0*m_num_steps[k])/m_num_steps[k-1]-1.0);

					m_subdiag->update(solfine, 1.0, solfine, (1.0/scaling), solfine,  -(1.0/scaling), solcoarse);
					number subdiag_error_est=m_subdiag->get_current_estimate();

					//m_subdiag_error_est[k] = sqrt(subdiag_error_est*scaling);
					m_subdiag_error_est[k]=subdiag_error_est*scaling;

					UG_LOG(" ErrorEst["<< k<<"]=" << m_subdiag_error_est[k] << ";" << std::endl);
				}

			}

		}

protected:
		number substep(size_t i) {return m_stepsize/m_num_steps[i];}

private:
		/** time step */
		number m_stepsize;
		static const int m_order=1;


		/** number of intermediate steps (per stage)*/
		std::vector<size_t> m_num_steps;

		/** vector of solutions (per stage)*/
		std::vector<SmartPtr<vector_type> > m_solution;

		/** sub-diagonal error estimate (per stage)*/
		std::vector<number> m_subdiag_error_est;

		/** evaluation on last stage */
		SmartPtr<ISubDiagErrorEst<vector_type> > m_subdiag;

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
