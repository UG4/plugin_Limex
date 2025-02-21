/*
 * SPDX-FileCopyrightText: Copyright (c) 2014-2025:  Goethe University Frankfurt
 * SPDX-License-Identifier: LicenseRef-UG4-LGPL-3.0
 *
 * Author: Arne Naegel
 *
 */
#pragma once

#include "bridge/util.h"

extern "C" void InitUGPlugin_Limex(ug::bridge::Registry* reg, std::string grp);

#ifdef UG_USE_PYBIND11

#include "bindings/pybind/ug_pybind.h"

namespace ug {
namespace Limex{
	void InitUGPlugin(ug::pybind::Registry* reg, std::string grp);
}
}
#endif
