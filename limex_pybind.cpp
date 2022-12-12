#include "limex_plugin.h"

#ifdef UG_USE_PYBIND11
PYBIND11_MODULE(pylimex, m)
{
	m.doc() = "Limex module";
	m.attr("__name__") = "ug4py.limex";

	ug::pybind::Registry registry(m);
	std::string name("Limex");

	ug::Limex::InitUGPlugin(&registry, name);
}
#endif
