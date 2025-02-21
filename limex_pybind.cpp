/*
 * SPDX-FileCopyrightText: Copyright (c) 2014-2025:  Goethe University Frankfurt
 * SPDX-License-Identifier: LicenseRef-UG4-LGPL-3.0
 *
 * Author: Arne Naegel
 *
 */

#ifdef UG_USE_PYBIND11

#include "limex_plugin.h"

PYBIND11_MODULE(pylimex, m)
{
	m.doc() = "Limex module";
	m.attr("__name__") = "ug4py.limex";

	ug::pybind::Registry registry(m);
	std::string name("Limex");

	ug::Limex::InitUGPlugin(&registry, name);
}
#endif
