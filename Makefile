# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/jul/.local/lib/python3.10/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/jul/.local/lib/python3.10/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jul/Documents/SimulationOfMat/git/DEAL_Stokes_Praktikum

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jul/Documents/SimulationOfMat/git/DEAL_Stokes_Praktikum

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/home/jul/.local/lib/python3.10/site-packages/cmake/data/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test
.PHONY : test/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/home/jul/.local/lib/python3.10/site-packages/cmake/data/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/home/jul/.local/lib/python3.10/site-packages/cmake/data/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jul/Documents/SimulationOfMat/git/DEAL_Stokes_Praktikum/CMakeFiles /home/jul/Documents/SimulationOfMat/git/DEAL_Stokes_Praktikum//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/jul/Documents/SimulationOfMat/git/DEAL_Stokes_Praktikum/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named stokes.g

# Build rule for target.
stokes.g: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 stokes.g
.PHONY : stokes.g

# fast build rule for target.
stokes.g/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/build
.PHONY : stokes.g/fast

#=============================================================================
# Target rules for targets named stokes_2d.g

# Build rule for target.
stokes_2d.g: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 stokes_2d.g
.PHONY : stokes_2d.g

# fast build rule for target.
stokes_2d.g/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_2d.g.dir/build.make CMakeFiles/stokes_2d.g.dir/build
.PHONY : stokes_2d.g/fast

#=============================================================================
# Target rules for targets named stokes_3d.g

# Build rule for target.
stokes_3d.g: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 stokes_3d.g
.PHONY : stokes_3d.g

# fast build rule for target.
stokes_3d.g/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_3d.g.dir/build.make CMakeFiles/stokes_3d.g.dir/build
.PHONY : stokes_3d.g/fast

#=============================================================================
# Target rules for targets named indent

# Build rule for target.
indent: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 indent
.PHONY : indent

# fast build rule for target.
indent/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/indent.dir/build.make CMakeFiles/indent.dir/build
.PHONY : indent/fast

#=============================================================================
# Target rules for targets named compile_test_executables

# Build rule for target.
compile_test_executables: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 compile_test_executables
.PHONY : compile_test_executables

# fast build rule for target.
compile_test_executables/fast:
	$(MAKE) $(MAKESILENT) -f tests/CMakeFiles/compile_test_executables.dir/build.make tests/CMakeFiles/compile_test_executables.dir/build
.PHONY : compile_test_executables/fast

#=============================================================================
# Target rules for targets named tests.template.debug

# Build rule for target.
tests.template.debug: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 tests.template.debug
.PHONY : tests.template.debug

# fast build rule for target.
tests.template.debug/fast:
	$(MAKE) $(MAKESILENT) -f tests/CMakeFiles/tests.template.debug.dir/build.make tests/CMakeFiles/tests.template.debug.dir/build
.PHONY : tests.template.debug/fast

#=============================================================================
# Target rules for targets named tests.template.debug.test

# Build rule for target.
tests.template.debug.test: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 tests.template.debug.test
.PHONY : tests.template.debug.test

# fast build rule for target.
tests.template.debug.test/fast:
	$(MAKE) $(MAKESILENT) -f tests/CMakeFiles/tests.template.debug.test.dir/build.make tests/CMakeFiles/tests.template.debug.test.dir/build
.PHONY : tests.template.debug.test/fast

source/BoundaryValues.o: source/BoundaryValues.cc.o
.PHONY : source/BoundaryValues.o

# target to build an object file
source/BoundaryValues.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/BoundaryValues.cc.o
.PHONY : source/BoundaryValues.cc.o

source/BoundaryValues.i: source/BoundaryValues.cc.i
.PHONY : source/BoundaryValues.i

# target to preprocess a source file
source/BoundaryValues.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/BoundaryValues.cc.i
.PHONY : source/BoundaryValues.cc.i

source/BoundaryValues.s: source/BoundaryValues.cc.s
.PHONY : source/BoundaryValues.s

# target to generate assembly for a file
source/BoundaryValues.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/BoundaryValues.cc.s
.PHONY : source/BoundaryValues.cc.s

source/InverseMatrix.o: source/InverseMatrix.cc.o
.PHONY : source/InverseMatrix.o

# target to build an object file
source/InverseMatrix.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/InverseMatrix.cc.o
.PHONY : source/InverseMatrix.cc.o

source/InverseMatrix.i: source/InverseMatrix.cc.i
.PHONY : source/InverseMatrix.i

# target to preprocess a source file
source/InverseMatrix.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/InverseMatrix.cc.i
.PHONY : source/InverseMatrix.cc.i

source/InverseMatrix.s: source/InverseMatrix.cc.s
.PHONY : source/InverseMatrix.s

# target to generate assembly for a file
source/InverseMatrix.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/InverseMatrix.cc.s
.PHONY : source/InverseMatrix.cc.s

source/RightHandSide.o: source/RightHandSide.cc.o
.PHONY : source/RightHandSide.o

# target to build an object file
source/RightHandSide.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/RightHandSide.cc.o
.PHONY : source/RightHandSide.cc.o

source/RightHandSide.i: source/RightHandSide.cc.i
.PHONY : source/RightHandSide.i

# target to preprocess a source file
source/RightHandSide.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/RightHandSide.cc.i
.PHONY : source/RightHandSide.cc.i

source/RightHandSide.s: source/RightHandSide.cc.s
.PHONY : source/RightHandSide.s

# target to generate assembly for a file
source/RightHandSide.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/RightHandSide.cc.s
.PHONY : source/RightHandSide.cc.s

source/SchurComplement.o: source/SchurComplement.cc.o
.PHONY : source/SchurComplement.o

# target to build an object file
source/SchurComplement.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/SchurComplement.cc.o
.PHONY : source/SchurComplement.cc.o

source/SchurComplement.i: source/SchurComplement.cc.i
.PHONY : source/SchurComplement.i

# target to preprocess a source file
source/SchurComplement.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/SchurComplement.cc.i
.PHONY : source/SchurComplement.cc.i

source/SchurComplement.s: source/SchurComplement.cc.s
.PHONY : source/SchurComplement.s

# target to generate assembly for a file
source/SchurComplement.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/SchurComplement.cc.s
.PHONY : source/SchurComplement.cc.s

source/StokesProblem.o: source/StokesProblem.cc.o
.PHONY : source/StokesProblem.o

# target to build an object file
source/StokesProblem.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/StokesProblem.cc.o
.PHONY : source/StokesProblem.cc.o

source/StokesProblem.i: source/StokesProblem.cc.i
.PHONY : source/StokesProblem.i

# target to preprocess a source file
source/StokesProblem.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/StokesProblem.cc.i
.PHONY : source/StokesProblem.cc.i

source/StokesProblem.s: source/StokesProblem.cc.s
.PHONY : source/StokesProblem.s

# target to generate assembly for a file
source/StokesProblem.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes.g.dir/build.make CMakeFiles/stokes.g.dir/source/StokesProblem.cc.s
.PHONY : source/StokesProblem.cc.s

source/main.o: source/main.cc.o
.PHONY : source/main.o

# target to build an object file
source/main.cc.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_2d.g.dir/build.make CMakeFiles/stokes_2d.g.dir/source/main.cc.o
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_3d.g.dir/build.make CMakeFiles/stokes_3d.g.dir/source/main.cc.o
.PHONY : source/main.cc.o

source/main.i: source/main.cc.i
.PHONY : source/main.i

# target to preprocess a source file
source/main.cc.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_2d.g.dir/build.make CMakeFiles/stokes_2d.g.dir/source/main.cc.i
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_3d.g.dir/build.make CMakeFiles/stokes_3d.g.dir/source/main.cc.i
.PHONY : source/main.cc.i

source/main.s: source/main.cc.s
.PHONY : source/main.s

# target to generate assembly for a file
source/main.cc.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_2d.g.dir/build.make CMakeFiles/stokes_2d.g.dir/source/main.cc.s
	$(MAKE) $(MAKESILENT) -f CMakeFiles/stokes_3d.g.dir/build.make CMakeFiles/stokes_3d.g.dir/source/main.cc.s
.PHONY : source/main.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... test"
	@echo "... compile_test_executables"
	@echo "... indent"
	@echo "... tests.template.debug.test"
	@echo "... stokes.g"
	@echo "... stokes_2d.g"
	@echo "... stokes_3d.g"
	@echo "... tests.template.debug"
	@echo "... source/BoundaryValues.o"
	@echo "... source/BoundaryValues.i"
	@echo "... source/BoundaryValues.s"
	@echo "... source/InverseMatrix.o"
	@echo "... source/InverseMatrix.i"
	@echo "... source/InverseMatrix.s"
	@echo "... source/RightHandSide.o"
	@echo "... source/RightHandSide.i"
	@echo "... source/RightHandSide.s"
	@echo "... source/SchurComplement.o"
	@echo "... source/SchurComplement.i"
	@echo "... source/SchurComplement.s"
	@echo "... source/StokesProblem.o"
	@echo "... source/StokesProblem.i"
	@echo "... source/StokesProblem.s"
	@echo "... source/main.o"
	@echo "... source/main.i"
	@echo "... source/main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

