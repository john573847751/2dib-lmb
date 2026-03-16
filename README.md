# LBE Solver - 格子玻尔兹曼求解器

基于浸入边界法的二维格子玻尔兹曼求解器（D2Q9模型）

## 项目结构

```
newib-lbm/
├── CMakeLists.txt          # CMake配置文件
├── README.md               # 项目说明
├── include/                # 头文件目录
│   ├── Constants.hpp       # 常量和参数定义
│   ├── Types.hpp           # 类型定义和数据结构
│   ├── LBESolver.hpp       # 求解器类定义
│   └── VTKWriter.hpp       # VTK文件写入器
└── src/                    # 源文件目录
    ├── main.cpp            # 主程序
    ├── LBESolver.cpp       # 求解器实现
    └── VTKWriter.cpp       # VTK写入器实现
```

## 特性

- **现代C++17**: 使用智能指针、constexpr等现代C++特性
- **OpenMP并行**: 支持多线程并行计算
- **VTK可视化**: 每200步自动保存VTK格式的结果
- **浸入边界法**: 支持复杂边界条件
- **模块化设计**: 清晰的代码结构和接口

## 编译与运行

### 1. 创建构建目录
```bash
mkdir build && cd build
```

### 2. 使用CMake配置项目
```bash
cmake ..
```

### 3. 编译
```bash
make -j
```

### 4. 运行
```bash
./lbe_solver
```

## 可视化结果

VTK文件保存在 `./vtk_output/` 目录中。可以使用以下软件查看：

- **ParaView**: 免费开源的VTK可视化软件
  ```bash
  paraview vtk_output/solution_000000.vtk
  ```

- **VisIt**: 高性能科学可视化工具

## 参数配置

在 `src/main.cpp` 中可以修改以下参数：

```cpp
const int nx = 100;          // X方向网格数
const int ny = 100;          // Y方向网格数
const double tau = 0.6;      // 弛豫时间
const double dt = 1.0;       // 时间步长
const int num_steps = 1000;  // 总步数
const int output_interval = 200;  // VTK输出间隔
```

## 代码说明

### 主要类

- **LBESolver**: 主求解器类，封装了LBE的所有计算逻辑
- **VTKWriter**: VTK文件写入器，支持标量场和向量场的可视化
- **BoundaryPoint**: 边界点数据结构，存储边界信息

### 算法流程

1. 初始化平衡分布函数
2. 对每个时间步：
   - 计算平衡分布
   - 碰撞步骤
   - 传播步骤
   - 计算宏观量
   - 边界速度修正
   - 定期保存VTK结果

### 物理模型

- 使用D2Q9离散速度模型
- 采用Bhatnagar-Gross-Krook (BGK) 碰撞算子
- 浸入边界法处理复杂边界

## 注意事项

- 需要支持C++17的编译器
- 需要安装OpenMP库
- VTK结果需要外部软件查看

## 性能优化

- 启用了O3优化和原生指令集
- 使用OpenMP实现并行化
- 采用现代C++特性提高代码效率

## 参考资料

- [The Lattice Boltzmann Method: Principles and Practice](https://doi.org/10.1007/978-3-319-44649-3)

