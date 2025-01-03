cmake_minimum_required(VERSION 3.14)
# project(Crossfire LANGUAGES CUDA CXX)
project(Crossfire CXX C CUDA)

# 指定使用 GCC 12 和 G++ 12
# set(CMAKE_C_COMPILER gcc-12)
# set(CMAKE_CXX_COMPILER g++-12)

set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_definitions(-mavx -mavx2 -mlzcnt -mbmi -DDEBUG=1)


add_subdirectory(PSA/WFA2-lib)

# 查找 TBB
find_package(TBB REQUIRED)

# 包含 Boost 库
find_package(Boost 1.71.0 COMPONENTS program_options REQUIRED)

# 包含 Threads 库
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# 查找并包含 OpenMP 库
find_package(OpenMP REQUIRED)


# 添加可执行文件时只包含源文件
add_executable(Crossfire
        main.cpp
        src/fcgr.cpp
        src/readfile.cpp
        src/cluster.cpp
        src/align_seq.cpp
        src/fmindex.cpp
        PSA/WFA2-lib/bindings/cpp/WFAligner.cpp
        include/bwt.cpp
        include/ssw.c
        include/ssw_cpp.cpp
        # include/ssw.cpp
)


# 启用 CUDA 可分离编译
set_target_properties(Crossfire PROPERTIES CUDA_SEPARABLE_COMPILATION ON)



# 链接 库
target_link_libraries(Crossfire PRIVATE
        TBB::tbb
        OpenMP::OpenMP_CXX
        Threads::Threads
        wfa2
        ${Boost_LIBRARIES}
)

set_target_properties(Crossfire PROPERTIES
    CUDA_ARCHITECTURES 89  # 根据您的 GPU 计算能力设置
)

# 包含头文件目录
target_include_directories(Crossfire PRIVATE
        ${PROJECT_SOURCE_DIR}/plog/include
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/PSA/WFA2-lib/bindings/cpp
)
