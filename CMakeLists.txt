################################################################################
# SAMPLE
#
# created by Sebastian Reiter
#
# This script shows how easy it is to add a build script for your plugins.
# Simply specify the plugins name, the path (relative to ug's root path) and
# a list of your source files.
################################################################################

# Name of your plugin, path to your plugin (relative to ug's root path) and sources.
set(pluginName	Limex)
set(SOURCES		limex_plugin.cpp
				time_disc/time_extrapolation.cpp
				time_disc/linear_implicit_timestep.cpp
)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 3.0)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

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

################################################################################
# Python binding
################################################################################

if(USE_PYBIND11)
	SET(PYBIND_SOURCES 	limex_pybind.cpp)
	SET(myLibraries ug4 ${pluginName})
	# First argument must match module name in PYBIND11_MODULE call
	ug4pybind_add_module(pylimex ${PYBIND_SOURCES} ${myLibraries})
	target_link_libraries (pylimex PRIVATE ${pluginName})
	set_target_properties(pylimex PROPERTIES INSTALL_RPATH "$ORIGIN/..:$ORIGIN/../../../lib")
	install(TARGETS pylimex LIBRARY DESTINATION bin/plugins/ug4 COMPONENT plugins)
endif(USE_PYBIND11)
