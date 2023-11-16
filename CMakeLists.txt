cmake_minimum_required(VERSION 3.15...3.26)

project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.11.1)
FetchContent_MakeAvailable(pybind11)
find_package(pybind11 CONFIG REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

add_library(pulchra src/cpp/pulchra.cpp)
target_compile_features(pulchra PRIVATE cxx_std_11)
set_property(TARGET pulchra PROPERTY POSITION_INDEPENDENT_CODE ON)
install(TARGETS pulchra DESTINATION pulchra)

python_add_library(_pypulchra MODULE src/python/pypulchra_bindings.cpp WITH_SOABI)
target_include_directories(_pypulchra PRIVATE src/cpp)
target_link_libraries(_pypulchra PRIVATE pybind11::headers pulchra)
target_compile_definitions(_pypulchra PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_compile_features(_pypulchra PRIVATE cxx_std_11)

install(TARGETS _pypulchra DESTINATION pypulchra)
