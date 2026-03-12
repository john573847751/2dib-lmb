#!/bin/bash

# LBE Solver 运行脚本

set -e

echo "==================================="
echo "Running LBE Solver"
echo "==================================="

# 检查可执行文件是否存在
if [ ! -f "build/lbe_solver" ]; then
    echo "Error: lbe_solver not found. Please build first using ./build.sh"
    exit 1
fi

# 创建vtk输出目录（如果不存在）
mkdir -p vtk_output

# 运行求解器
cd build
./lbe_solver

echo "==================================="
echo "Simulation completed!"
echo "==================================="
echo ""
echo "Results:"
echo "  - VTK files: vtk_output/"
echo "  - CSV results: final_results.csv"
echo ""
echo "To visualize with ParaView:"
echo "  paraview vtk_output/solution_*.vtk"
echo ""
