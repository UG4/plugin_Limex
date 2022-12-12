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
cmake_minimum_required(VERSION 2.8)

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
endif(buildEmbeddedPlugins)

################################################################################
# Python binding
################################################################################
if(USE_PYBIND11)
	SET(pyPluginName pylimex) # Must match module name after PYBIND11_MODULE call
	pybind11_add_module(${pyPluginName} limex_pybind.cpp)
	set_target_properties(${pyPluginName} PROPERTIES	
  		CXX_STANDARD 11 
  		CXX_STANDARD_REQUIRED YES
   		CXX_EXTENSIONS NO
   		OUTPUT_NAME ug4py/${pyPluginName})
   	# We need to link against 'this' plugin and UG4 library.
	target_link_libraries (${pyPluginName} PRIVATE ug4 ${pluginName})
endif(USE_PYBIND11)
