# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/alic/src/hpstr

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alic/src/hpstr/build

# Include any dependencies generated for this target.
include processors/CMakeFiles/processors.dir/depend.make

# Include the progress variables for this target.
include processors/CMakeFiles/processors.dir/progress.make

# Include the compile flags for this target's objects.
include processors/CMakeFiles/processors.dir/flags.make

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o: ../processors/src/ClusterOnTrackAnaProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/ClusterOnTrackAnaProcessor.cxx

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/ClusterOnTrackAnaProcessor.cxx > CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/ClusterOnTrackAnaProcessor.cxx -o CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o: ../processors/src/ECalDataProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/ECalDataProcessor.cxx

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/ECalDataProcessor.cxx > CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/ECalDataProcessor.cxx -o CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o: ../processors/src/EventProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/EventProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/EventProcessor.cxx

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/EventProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/EventProcessor.cxx > CMakeFiles/processors.dir/src/EventProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/EventProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/EventProcessor.cxx -o CMakeFiles/processors.dir/src/EventProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o: ../processors/src/HPSEventProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/HPSEventProcessor.cxx

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/HPSEventProcessor.cxx > CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/HPSEventProcessor.cxx -o CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o: ../processors/src/ParticleProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/ParticleProcessor.cxx

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/ParticleProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/ParticleProcessor.cxx > CMakeFiles/processors.dir/src/ParticleProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/ParticleProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/ParticleProcessor.cxx -o CMakeFiles/processors.dir/src/ParticleProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o: ../processors/src/RefittedTracksProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/RefittedTracksProcessor.cxx

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/RefittedTracksProcessor.cxx > CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/RefittedTracksProcessor.cxx -o CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o: ../processors/src/SvtCondAnaProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/SvtCondAnaProcessor.cxx

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/SvtCondAnaProcessor.cxx > CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/SvtCondAnaProcessor.cxx -o CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o: ../processors/src/SvtDataProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/SvtDataProcessor.cxx

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/SvtDataProcessor.cxx > CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/SvtDataProcessor.cxx -o CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o: ../processors/src/SvtRawDataProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/SvtRawDataProcessor.cxx

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/SvtRawDataProcessor.cxx > CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/SvtRawDataProcessor.cxx -o CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o: ../processors/src/TrackingProcessor.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o -c /home/alic/src/hpstr/processors/src/TrackingProcessor.cxx

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/TrackingProcessor.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/TrackingProcessor.cxx > CMakeFiles/processors.dir/src/TrackingProcessor.cxx.i

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/TrackingProcessor.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/TrackingProcessor.cxx -o CMakeFiles/processors.dir/src/TrackingProcessor.cxx.s

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.requires

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.provides: processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.provides

processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o


processors/CMakeFiles/processors.dir/src/utilities.cxx.o: processors/CMakeFiles/processors.dir/flags.make
processors/CMakeFiles/processors.dir/src/utilities.cxx.o: ../processors/src/utilities.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object processors/CMakeFiles/processors.dir/src/utilities.cxx.o"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/processors.dir/src/utilities.cxx.o -c /home/alic/src/hpstr/processors/src/utilities.cxx

processors/CMakeFiles/processors.dir/src/utilities.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/processors.dir/src/utilities.cxx.i"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alic/src/hpstr/processors/src/utilities.cxx > CMakeFiles/processors.dir/src/utilities.cxx.i

processors/CMakeFiles/processors.dir/src/utilities.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/processors.dir/src/utilities.cxx.s"
	cd /home/alic/src/hpstr/build/processors && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alic/src/hpstr/processors/src/utilities.cxx -o CMakeFiles/processors.dir/src/utilities.cxx.s

processors/CMakeFiles/processors.dir/src/utilities.cxx.o.requires:

.PHONY : processors/CMakeFiles/processors.dir/src/utilities.cxx.o.requires

processors/CMakeFiles/processors.dir/src/utilities.cxx.o.provides: processors/CMakeFiles/processors.dir/src/utilities.cxx.o.requires
	$(MAKE) -f processors/CMakeFiles/processors.dir/build.make processors/CMakeFiles/processors.dir/src/utilities.cxx.o.provides.build
.PHONY : processors/CMakeFiles/processors.dir/src/utilities.cxx.o.provides

processors/CMakeFiles/processors.dir/src/utilities.cxx.o.provides.build: processors/CMakeFiles/processors.dir/src/utilities.cxx.o


# Object files for target processors
processors_OBJECTS = \
"CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/EventProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o" \
"CMakeFiles/processors.dir/src/utilities.cxx.o"

# External object files for target processors
processors_EXTERNAL_OBJECTS =

processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/src/utilities.cxx.o
processors/libprocessors.so: processors/CMakeFiles/processors.dir/build.make
processors/libprocessors.so: processing/libprocessing.so
processors/libprocessors.so: analysis/libanalysis.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libCore.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libImt.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libRIO.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libNet.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libHist.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGraf.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGraf3d.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGpad.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libROOTDataFrame.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libTree.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libTreePlayer.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libRint.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libPostscript.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libMatrix.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libPhysics.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libMathCore.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libThread.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libMultiProc.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libPyROOT.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGeom.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libEve.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGui.so
processors/libprocessors.so: /home/alic/work/hps/hpsgit/LCIO/install/lib/liblcio.so
processors/libprocessors.so: /usr/lib/x86_64-linux-gnu/libpython2.7.so
processors/libprocessors.so: event/libevent.so
processors/libprocessors.so: /home/alic/work/hps/hpsgit/LCIO/install/lib/liblcio.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libCore.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libImt.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libRIO.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libNet.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libHist.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGraf.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGraf3d.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGpad.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libROOTDataFrame.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libTree.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libTreePlayer.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libRint.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libPostscript.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libMatrix.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libPhysics.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libMathCore.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libThread.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libMultiProc.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libPyROOT.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGeom.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libEve.so
processors/libprocessors.so: /home/alic/src/root-6.18.04/buildV61804/lib/libGui.so
processors/libprocessors.so: processors/CMakeFiles/processors.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/alic/src/hpstr/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX shared library libprocessors.so"
	cd /home/alic/src/hpstr/build/processors && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/processors.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
processors/CMakeFiles/processors.dir/build: processors/libprocessors.so

.PHONY : processors/CMakeFiles/processors.dir/build

processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/ClusterOnTrackAnaProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/ECalDataProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/EventProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/HPSEventProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/ParticleProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/RefittedTracksProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/SvtCondAnaProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/SvtDataProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/SvtRawDataProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/TrackingProcessor.cxx.o.requires
processors/CMakeFiles/processors.dir/requires: processors/CMakeFiles/processors.dir/src/utilities.cxx.o.requires

.PHONY : processors/CMakeFiles/processors.dir/requires

processors/CMakeFiles/processors.dir/clean:
	cd /home/alic/src/hpstr/build/processors && $(CMAKE_COMMAND) -P CMakeFiles/processors.dir/cmake_clean.cmake
.PHONY : processors/CMakeFiles/processors.dir/clean

processors/CMakeFiles/processors.dir/depend:
	cd /home/alic/src/hpstr/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alic/src/hpstr /home/alic/src/hpstr/processors /home/alic/src/hpstr/build /home/alic/src/hpstr/build/processors /home/alic/src/hpstr/build/processors/CMakeFiles/processors.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : processors/CMakeFiles/processors.dir/depend

