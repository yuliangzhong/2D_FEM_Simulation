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

# Include any dependencies generated for this target.
include ext/imgui/CMakeFiles/imgui.dir/depend.make

# Include the progress variables for this target.
include ext/imgui/CMakeFiles/imgui.dir/progress.make

# Include the compile flags for this target's objects.
include ext/imgui/CMakeFiles/imgui.dir/flags.make

ext/imgui/CMakeFiles/imgui.dir/imgui.cpp.o: ext/imgui/CMakeFiles/imgui.dir/flags.make
ext/imgui/CMakeFiles/imgui.dir/imgui.cpp.o: ../ext/imgui/imgui.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a4-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ext/imgui/CMakeFiles/imgui.dir/imgui.cpp.o"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui.cpp.o -c /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui.cpp

ext/imgui/CMakeFiles/imgui.dir/imgui.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui.cpp.i"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui.cpp > CMakeFiles/imgui.dir/imgui.cpp.i

ext/imgui/CMakeFiles/imgui.dir/imgui.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui.cpp.s"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui.cpp -o CMakeFiles/imgui.dir/imgui.cpp.s

ext/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o: ext/imgui/CMakeFiles/imgui.dir/flags.make
ext/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o: ../ext/imgui/imgui_demo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a4-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object ext/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui_demo.cpp.o -c /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_demo.cpp

ext/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui_demo.cpp.i"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_demo.cpp > CMakeFiles/imgui.dir/imgui_demo.cpp.i

ext/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui_demo.cpp.s"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_demo.cpp -o CMakeFiles/imgui.dir/imgui_demo.cpp.s

ext/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o: ext/imgui/CMakeFiles/imgui.dir/flags.make
ext/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o: ../ext/imgui/imgui_draw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a4-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object ext/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui_draw.cpp.o -c /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_draw.cpp

ext/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui_draw.cpp.i"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_draw.cpp > CMakeFiles/imgui.dir/imgui_draw.cpp.i

ext/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui_draw.cpp.s"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_draw.cpp -o CMakeFiles/imgui.dir/imgui_draw.cpp.s

ext/imgui/CMakeFiles/imgui.dir/imgui_widgets.cpp.o: ext/imgui/CMakeFiles/imgui.dir/flags.make
ext/imgui/CMakeFiles/imgui.dir/imgui_widgets.cpp.o: ../ext/imgui/imgui_widgets.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ubuntu/CMM/a4-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object ext/imgui/CMakeFiles/imgui.dir/imgui_widgets.cpp.o"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imgui.dir/imgui_widgets.cpp.o -c /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_widgets.cpp

ext/imgui/CMakeFiles/imgui.dir/imgui_widgets.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgui.dir/imgui_widgets.cpp.i"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_widgets.cpp > CMakeFiles/imgui.dir/imgui_widgets.cpp.i

ext/imgui/CMakeFiles/imgui.dir/imgui_widgets.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgui.dir/imgui_widgets.cpp.s"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui/imgui_widgets.cpp -o CMakeFiles/imgui.dir/imgui_widgets.cpp.s

# Object files for target imgui
imgui_OBJECTS = \
"CMakeFiles/imgui.dir/imgui.cpp.o" \
"CMakeFiles/imgui.dir/imgui_demo.cpp.o" \
"CMakeFiles/imgui.dir/imgui_draw.cpp.o" \
"CMakeFiles/imgui.dir/imgui_widgets.cpp.o"

# External object files for target imgui
imgui_EXTERNAL_OBJECTS =

ext/imgui/libimgui.a: ext/imgui/CMakeFiles/imgui.dir/imgui.cpp.o
ext/imgui/libimgui.a: ext/imgui/CMakeFiles/imgui.dir/imgui_demo.cpp.o
ext/imgui/libimgui.a: ext/imgui/CMakeFiles/imgui.dir/imgui_draw.cpp.o
ext/imgui/libimgui.a: ext/imgui/CMakeFiles/imgui.dir/imgui_widgets.cpp.o
ext/imgui/libimgui.a: ext/imgui/CMakeFiles/imgui.dir/build.make
ext/imgui/libimgui.a: ext/imgui/CMakeFiles/imgui.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ubuntu/CMM/a4-yuliangzhong/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libimgui.a"
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean_target.cmake
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imgui.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ext/imgui/CMakeFiles/imgui.dir/build: ext/imgui/libimgui.a

.PHONY : ext/imgui/CMakeFiles/imgui.dir/build

ext/imgui/CMakeFiles/imgui.dir/clean:
	cd /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui && $(CMAKE_COMMAND) -P CMakeFiles/imgui.dir/cmake_clean.cmake
.PHONY : ext/imgui/CMakeFiles/imgui.dir/clean

ext/imgui/CMakeFiles/imgui.dir/depend:
	cd /home/ubuntu/CMM/a4-yuliangzhong/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ubuntu/CMM/a4-yuliangzhong /home/ubuntu/CMM/a4-yuliangzhong/ext/imgui /home/ubuntu/CMM/a4-yuliangzhong/build /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui /home/ubuntu/CMM/a4-yuliangzhong/build/ext/imgui/CMakeFiles/imgui.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ext/imgui/CMakeFiles/imgui.dir/depend

