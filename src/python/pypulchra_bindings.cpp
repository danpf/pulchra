#include <algorithm>
#include <array>
#include <functional>
#include <ios>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>

#include <pulchra.hpp>

typedef std::function<pybind11::module &(std::string const &)> ModuleGetter;
#ifndef BINDER_PYBIND11_TYPE_CASTER
#define BINDER_PYBIND11_TYPE_CASTER
PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)
PYBIND11_DECLARE_HOLDER_TYPE(T, T *)
PYBIND11_MAKE_OPAQUE(std::shared_ptr<void>)
#endif

PYBIND11_MODULE(_pypulchra, root_module) {
  root_module.doc() = "_pypulchra module";

  std::map<std::string, pybind11::module> modules;
  ModuleGetter M = [&](std::string const &namespace_) -> pybind11::module & {
    auto it = modules.find(namespace_);
    if (it == modules.end())
      throw std::runtime_error(
          "Attempt to access pybind11::module for namespace " + namespace_ +
          " before it was created!!!");
    return it->second;
  };

  modules[""] = root_module;

  static std::vector<std::string> const reserved_python_words{
      "nonlocal",
      "global",
  };

  auto mangle_namespace_name([](std::string const &ns) -> std::string {
    if (std::find(reserved_python_words.begin(), reserved_python_words.end(),
                  ns) == reserved_python_words.end())
      return ns;
    else
      return ns + '_';
  });

  std::vector<std::pair<std::string, std::string>> const sub_modules{
      {"", "pulchra"},
  };
  for (auto const &p : sub_modules)
    modules[p.first.size() ? p.first + "::" + p.second : p.second] =
        modules[p.first].def_submodule(
            mangle_namespace_name(p.second).c_str(),
            ("Bindings for " + p.first + "::" + p.second + " namespace")
                .c_str());

  // pybind11::class_<std::shared_ptr<void>>(M(""), "_encapsulated_data_");

  { // pulchra::PulchraResult file:pulchra.hpp line:50
    pybind11::class_<pulchra::PulchraResult,
                     std::shared_ptr<pulchra::PulchraResult>>
        cl(M("pulchra"), "PulchraResult", "");
    cl.def(pybind11::init<const std::string &, const std::string &>(),
           pybind11::arg("_pdb_str"), pybind11::arg("_traj_pdb_str"));

    cl.def(pybind11::init([](pulchra::PulchraResult const &o) {
      return new pulchra::PulchraResult(o);
    }));
    cl.def_readwrite("pdb_str", &pulchra::PulchraResult::pdb_str);
    cl.def_readwrite("traj_pdb_str", &pulchra::PulchraResult::traj_pdb_str);
    cl.def(
        "assign",
        (struct pulchra::PulchraResult &
         (pulchra::PulchraResult::*)(const struct pulchra::PulchraResult &)) &
            pulchra::PulchraResult::operator=,
        "pulchra.PulchraResult.operator=(const struct pulchra.PulchraResult &) "
        "--> struct pulchra.PulchraResult &",
        pybind11::return_value_policy::automatic, pybind11::arg(""));
  }
  // pulchra::run_from_pdb_str(const std::string &) file:pulchra.hpp line:3295
  M("pulchra").def("run_from_pdb_str",
                   (struct pulchra::PulchraResult(*)(const std::string &)) &
                       pulchra::run_from_pdb_str,
                   "pulchra.run_from_pdb_str(const std.string &) --> struct "
                   "pulchra.PulchraResult",
                   pybind11::arg("pdb_str"));
  M("").attr("run_from_pdb_str") = M("pulchra").attr("run_from_pdb_str");
}
