#include <array>
#include <ios>
#include <iterator>
#include <memory>
#include <pulchra.hpp>
#include <sstream>
#include <sstream> // __str__
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <functional>
#include <string>

#ifndef BINDER_PYBIND11_TYPE_CASTER
	#define BINDER_PYBIND11_TYPE_CASTER
	PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
	PYBIND11_DECLARE_HOLDER_TYPE(T, T*)
	PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

void bind_pulchra(std::function< pybind11::module &(std::string const &namespace_) > &M)
{
	{ // pulchra::PulchraResult file:pulchra.hpp line:50
		pybind11::class_<pulchra::PulchraResult, std::shared_ptr<pulchra::PulchraResult>> cl(M("pulchra"), "PulchraResult", "");
		cl.def( pybind11::init<const std::string &, const std::string &>(), pybind11::arg("_pdb_str"), pybind11::arg("_traj_pdb_str") );

		cl.def( pybind11::init( [](pulchra::PulchraResult const &o){ return new pulchra::PulchraResult(o); } ) );
		cl.def_readwrite("pdb_str", &pulchra::PulchraResult::pdb_str);
		cl.def_readwrite("traj_pdb_str", &pulchra::PulchraResult::traj_pdb_str);
		cl.def("assign", (struct pulchra::PulchraResult & (pulchra::PulchraResult::*)(const struct pulchra::PulchraResult &)) &pulchra::PulchraResult::operator=, "pulchra.PulchraResult.operator=(const struct pulchra.PulchraResult &) --> struct pulchra.PulchraResult &", pybind11::return_value_policy::automatic, pybind11::arg(""));
	}
	// pulchra::run_from_pdb_str(const std::string &) file:pulchra.hpp line:3295
	M("pulchra").def("run_from_pdb_str", (struct pulchra::PulchraResult (*)(const std::string &)) &pulchra::run_from_pdb_str, "pulchra.run_from_pdb_str(const std.string &) --> struct pulchra.PulchraResult", pybind11::arg("pdb_str"));
	M("").attr("run_from_pdb_str") = M("pulchra").attr("run_from_pdb_str");
}
