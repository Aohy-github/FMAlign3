cmake_minimum_required(VERSION 3.14)

project(Crossfire CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_definitions(-mavx -mavx2 -mlzcnt -mbmi -DDEBUG=1)
set(TBB_TEST OFF CACHE BOOL "Enable testing" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "Disable building oneTBB examples" FORCE)

add_subdirectory(submodule/WFA2-lib)
add_subdirectory(submodule/oneTBB)
add_subdirectory(submodule/plog)
add_subdirectory(submodule/spdlog)



find_package(Boost 1.71.0 COMPONENTS program_options REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)


find_package(OpenMP REQUIRED)



add_executable(Crossfire
        main.cpp
        src/fcgr.cpp
        src/readfile.cpp
        src/cluster.cpp
        src/align_seq.cpp
        submodule/WFA2-lib/bindings/cpp/WFAligner.cpp
        include/bwt.cpp
        include/ssw.c
        include/ssw_cpp.cpp
)


target_link_libraries(Crossfire PRIVATE
        TBB::tbb
        OpenMP::OpenMP_CXX
        Threads::Threads
        wfa2::wfa2
        plog::plog
        ${Boost_LIBRARIES}
        spdlog::spdlog
)


target_include_directories(Crossfire PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)
