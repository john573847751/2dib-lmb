cmake_minimum_required(VERSION 3.10)
project(LBESolver VERSION 1.0.0 LANGUAGES CXX)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置构建类型
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 添加编译选项
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -O3 -march=native)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g -O0)
    endif()
elseif(MSVC)
    add_compile_options(/W4)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(/Od /Zi)
    endif()
endif()

# 查找OpenMP
find_package(OpenMP REQUIRED)

# 包含目录
include_directories(${CMAKE_SOURCE_DIR}/include)

# 源文件
set(SOURCES
    src/main.cpp
    src/LBESolver.cpp
    src/VTKWriter.cpp
    src/ImmersedBoundary.cpp
)

# 创建可执行文件
add_executable(lbe_solver ${SOURCES})

# 链接OpenMP
target_link_libraries(lbe_solver PRIVATE OpenMP::OpenMP_CXX)

# 安装规则
install(TARGETS lbe_solver DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)

# 打印配置信息
message(STATUS "")
message(STATUS "LBE Solver Configuration:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "  OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "  Binary directory: ${CMAKE_BINARY_DIR}")
message(STATUS "")
