#!/usr/bin/env python3

import glob
import os
import sys
import shutil
import subprocess
from distutils.sysconfig import get_python_inc
from pathlib import Path
import argparse
from typing import List


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binder-executable", required=True, help="Where binder and clang are, must be built and in the same directory"
    )
    parser.add_argument(
        "--output-directory", required=True, help="directory to build/output the cmake bindings and executables"
    )
    parser.add_argument(
        "--module-name", required=True, help="what you would like to call this project (ie- import module-name)"
    )
    parser.add_argument("--project-source", required=True, help="the location of the project's source files")
    parser.add_argument("--pybind11-source", required=True, help="the location of the pybind11 source directory")

    parser.add_argument(
        "--source-directories-to-include",
        nargs="+",
        required=False,
        default=[],
        help="If you require any extra source directories to build your project, list them here",
    )

    parser.add_argument("--config-file", help="A config file for you to use")
    parser.add_argument("--extra-binder-flags", default="", help="Extra binder flags, for debugging typically you use: --trace --annotate-includes")

    return parser.parse_args()


def get_all_project_source_files(project_source: str) -> List[str]:
    ps_pth = Path(project_source)
    lf = lambda x: list(ps_pth.rglob(x))
    all_source_files = lf("*.hpp") + lf("*.cpp") + lf("*.h") + lf("*.hh") + lf("*.cc") + lf("*.c")
    return [str(x) for x in all_source_files]

def make_and_write_all_includes(all_project_source_files: List[str], out_all_includes_fn: str) -> None:
    all_includes = []
    for filename in all_project_source_files:
        with open(filename, "r") as fh:
            for line in fh:
                if line.startswith("#include"):
                    line = line.strip()
                    # if '"' in line:
                    #     line = line.replace('"', "<")[:-1] + ">"
                    all_includes.append(line)
    all_includes = list(set(all_includes))
    # This is to ensure that the list is always the same and doesn't
    # depend on the filesystem state.  Not technically necessary, but
    # will cause inconsistent errors without it.
    all_includes.sort()
    with open(out_all_includes_fn, "w") as fh:
        for include in all_includes:
            fh.write(f"{include}\n")


def make_bindings_code(
        binder_executable: str,
        all_includes_fn: str,
        bindings_dir: str,
        python_module_name: str,
        extra_source_directories: List[str],
        extra_binder_flags: str,
    config_file: str,
):
    includes = " ".join([f"-I{x}" for x in extra_source_directories])
    if config_file:
        config_file = f" --config {config_file}"
    else:
        config_file = ""
    command = (
        f"{binder_executable} --root-module {python_module_name}"
        f" --prefix {bindings_dir} {extra_binder_flags}"
        f" {config_file}"
        f" {all_includes_fn} -- -std=c++11"
        f" {includes} -DNDEBUG -v"
    )
    print("running command", command)
    ret = subprocess.run(command.split())
    if ret.returncode != 0:
        raise RuntimeError(f"Bad command return {command}")
    sources_to_compile = []
    with open(os.path.join(bindings_dir, f"{python_module_name}.sources"), "r") as fh:
        for line in fh:
            l = line.strip()
            if l in sources_to_compile:
                raise RuntimeError("WARNING - DUPLICATED SOURCE - DO NOT NAME YOUR MODULE THE SAME AS ONE OF YOUR NAMESPACES/CLASSES")
            sources_to_compile.append(l)
    return sources_to_compile


def compile_sources(sources_to_compile: List[str], bindings_dir: str, module_name: str, directories_to_include: List[str], pybind11_source: str,
        all_project_source_files: List[str]):
    lines_to_write = []
    lines_to_write.append(f"cmake_minimum_required(VERSION 3.4...3.18)")
    lines_to_write.append(f"project({module_name})")
    lines_to_write.append(f"add_subdirectory(\"{pybind11_source}\" \"${{CMAKE_CURRENT_BINARY_DIR}}/pybind11_build\")")
    lines_to_write.append("")
    # lines_to_write.append(f"add_subdirectory(\"{pybind11" "${CMAKE_CURRENT_BINARY_DIR}/testlib_build"))
    tolink = []
    for source_fn in all_project_source_files:
        if "c" in Path(source_fn).suffix:
            source_lib_name = str(source_fn).replace("/", "_").replace(".","_")
            lines_to_write.append(f"add_library({source_lib_name} STATIC ${{CMAKE_SOURCE_DIR}}/{source_fn})")
            lines_to_write.append(f"set_target_properties({source_lib_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)")
            lines_to_write.append(f"target_include_directories({source_lib_name} PRIVATE {' '.join(directories_to_include)})")
            tolink.append(source_lib_name)

    sources_to_compile_for_cmake = ["${CMAKE_CURRENT_BINARY_DIR}/"+x for x in sources_to_compile]
    lines_to_write.append(f"pybind11_add_module({module_name} MODULE {' '.join(sources_to_compile_for_cmake)})")
    lines_to_write.append(f"target_include_directories({module_name} PRIVATE {' '.join(directories_to_include)})")
    lines_to_write.append(f"set_target_properties({module_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)")

    lines_to_write.append(f"target_link_libraries({module_name} PRIVATE {' '.join(tolink)})")


    with open("CMakeLists.txt", "w") as f:
        for line in lines_to_write:
            f.write(f"{line}\n")

    # Done making CMakeLists.txt
    subprocess.run("cmake -G Ninja -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang ..".split(), cwd=bindings_dir)
    subprocess.run("ninja -v".split(), cwd=bindings_dir)
    if bindings_dir not in sys.path:
        sys.path.append(bindings_dir)

    sys.stdout.flush()
    print("Testing Python lib...")
    import pypulchra

    print(dir(pypulchra))
    print(pypulchra)


def main():
    args = parseargs()
    extra_source_directories = [args.project_source, get_python_inc(), "pybind11/include"] + args.source_directories_to_include
    extra_source_directories = [str(Path(x).resolve()) for x in extra_source_directories]

    shutil.rmtree(args.output_directory, ignore_errors=True)
    Path(args.output_directory).mkdir(exist_ok=True, parents=True)
    all_includes_fn = str((Path(args.output_directory) / "all_includes.hpp").resolve())

    all_project_source_files = get_all_project_source_files(args.project_source)
    make_and_write_all_includes(all_project_source_files, all_includes_fn)
    sources_to_compile = make_bindings_code(
        args.binder_executable,
        all_includes_fn,
        args.output_directory,
        args.module_name,
        extra_source_directories,
        args.extra_binder_flags,
        args.config_file,
        
    )
    compile_sources(sources_to_compile, args.output_directory, args.module_name, extra_source_directories, args.pybind11_source,
        all_project_source_files)


if __name__ == "__main__":
    main()
