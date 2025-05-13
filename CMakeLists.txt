################################################################################
# SAMPLE
#
# created by Sebastian Reiter
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name, the path (relative to ug's root path) and
# a list of your source files.
################################################################################

# SPDX-FileCopyrightText:  Copyright (c) 2025:  G-CSC, Goethe University Frankfurt
# SPDX-License-Identifier: LicenseRef-UG4-LGPL-3.0

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	Limex)
set(SOURCES		limex_plugin.cpp
				time_disc/time_extrapolation.cpp
				time_disc/linear_implicit_timestep.cpp
)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 3.0...3.15)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

################################################################################
# Classic binding (static or dynamic plugin).
################################################################################
if(NOT USE_PYBIND11)
if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
	set_target_properties(${pluginName} PROPERTIES INSTALL_RPATH "$ORIGIN/../../lib/")
	install(TARGETS ${pluginName} LIBRARY DESTINATION bin/plugins COMPONENT plugins)
endif(buildEmbeddedPlugins)
endif(NOT USE_PYBIND11)

################################################################################
# Python binding (static plugin, dynamic python interface).
################################################################################
if(USE_PYBIND11)
	# Create a STATIC library for UG4 stuff from the sources
	add_library(${pluginName} STATIC ${SOURCES})
	
	SET(myPluginSources limex_pybind.cpp)
	SET(myLibraries ug4 ${pluginName})
	
	# First argument must match module name in PYBIND11_MODULE call
	python_add_library(pylimex MODULE  ${myPluginSources} ${SOURCES} WITH_SOABI)
	#ug4pybind_add_module(pylimex ${myPluginSources} ${myLibraries})
	target_link_libraries (pylimex PRIVATE ${pluginName} ug4_s)
	set_target_properties(pylimex PROPERTIES INSTALL_RPATH "$ORIGIN/..:$ORIGIN/../../../lib")
	install(TARGETS pylimex LIBRARY DESTINATION ug4py COMPONENT pymodules)
endif(USE_PYBIND11)
