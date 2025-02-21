/*
 * SPDX-FileCopyrightText: Copyright (c) 2014-2025:  G-CSC, Goethe University Frankfurt
 * SPDX-License-Identifier: LicenseRef-UG4-LGPL-3.0
 *
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

#ifndef __H__LIMEX__DATA_OUTPUT_OBSERVER_HPP__
#define __H__LIMEX__DATA_OUTPUT_OBSERVER_HPP__

#include <string>
#include <vector>

#include "lib_disc/time_disc/time_integrator_observers/time_integrator_observer_interface.h"
#include "lib_disc/io/vtkoutput.h"

namespace ug {

/// Class for integration observer: Output to VTK
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

/// Class for integration observer: Output to ConnectionViewer
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

} // namespace ug

#endif /* __H__LIMEX__DATA_OUTPUT_OBSERVER_HPP__ */
