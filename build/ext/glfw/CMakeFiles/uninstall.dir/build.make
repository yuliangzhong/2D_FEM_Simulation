# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ubuntu/CMM/a4-yuliangzhong

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ubuntu/CMM/a4-yuliangzhong/build

# Utility rule file for uninstall.

# Include the progress variables for this target.
include ext/glfw/CMakeFiles/uninstall.dir/progress.make

ext/glfw/CMakeFiles/uninstall:
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/glfw && /usr/local/bin/cmake -P /home/ubuntu/CMM/a4-yuliangzhong/build/ext/glfw/cmake_uninstall.cmake

uninstall: ext/glfw/CMakeFiles/uninstall
uninstall: ext/glfw/CMakeFiles/uninstall.dir/build.make

.PHONY : uninstall

# Rule to build all files generated by this target.
ext/glfw/CMakeFiles/uninstall.dir/build: uninstall

.PHONY : ext/glfw/CMakeFiles/uninstall.dir/build

ext/glfw/CMakeFiles/uninstall.dir/clean:
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/glfw && $(CMAKE_COMMAND) -P CMakeFiles/uninstall.dir/cmake_clean.cmake
.PHONY : ext/glfw/CMakeFiles/uninstall.dir/clean

ext/glfw/CMakeFiles/uninstall.dir/depend:
	cd /home/ubuntu/CMM/a4-yuliangzhong/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/CMM/a4-yuliangzhong /home/ubuntu/CMM/a4-yuliangzhong/ext/glfw /home/ubuntu/CMM/a4-yuliangzhong/build /home/ubuntu/CMM/a4-yuliangzhong/build/ext/glfw /home/ubuntu/CMM/a4-yuliangzhong/build/ext/glfw/CMakeFiles/uninstall.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/glfw/CMakeFiles/uninstall.dir/depend
