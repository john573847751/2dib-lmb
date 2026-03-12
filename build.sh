#!/bin/bash

# LBE Solver 构建脚本

set -e

echo "==================================="
echo "LBE Solver Build Script"
echo "==================================="

# 创建构建目录
if [ ! -d "build" ]; then
    echo "Creating build directory..."
    mkdir build
fi

cd build

echo "Configuring with CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

echo "Building..."
make -j$(nproc)

echo "==================================="
echo "Build completed successfully!"
echo "==================================="
echo ""
echo "To run the solver:"
echo "  cd build && ./lbe_solver"
echo ""
echo "VTK files will be saved in: ./vtk_output/"
echo "Final results will be saved to: final_results.csv"
echo ""
