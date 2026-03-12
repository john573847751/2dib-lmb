#pragma once

#include <vector>
#include <memory>
#include "Constants.hpp"

namespace lbm {

// 前向声明
class LBESolver;

// 类型别名
template<typename T>
using Vector2D = std::vector<std::vector<T>>;

template<typename T>
using Vector3D = std::vector<std::vector<std::vector<T>>>;

using DistributionFunction = Vector3D<double>;
using ScalarField = Vector2D<double>;
using VectorField = std::pair<ScalarField, ScalarField>;

// 边界点信息结构体
struct BoundaryPoint {
    double x{0.0};      // x坐标
    double y{0.0};      // y坐标
    double ux{0.0};     // x方向速度
    double uy{0.0};     // y方向速度
    double ds{0.0};     // 边界段长度
    
    BoundaryPoint() = default;
    BoundaryPoint(double _x, double _y, double _ux = 0.0, double _uy = 0.0, double _ds = 0.1)
        : x(_x), y(_y), ux(_ux), uy(_uy), ds(_ds) {}
};

// 网格尺寸结构体
struct GridSize {
    int nx{0};
    int ny{0};
    
    GridSize(int _nx, int _ny) : nx(_nx), ny(_ny) {}
};

// 模拟参数结构体
struct SimulationParams {
    double tau{0.6};        // 弛豫时间
    double dt{1.0};         // 时间步长
    int num_steps{1000};    // 总步数
    int output_interval{200}; // 输出间隔
    
    SimulationParams(double _tau = 0.6, double _dt = 1.0, 
                     int _steps = 1000, int _interval = 200)
        : tau(_tau), dt(_dt), num_steps(_steps), output_interval(_interval) {}
};

} // namespace lbm
