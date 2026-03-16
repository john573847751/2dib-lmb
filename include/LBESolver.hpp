#pragma once

#include "Constants.hpp"
#include "Types.hpp"
#include "VTKWriter.hpp"
#include "ImmersedBoundary.hpp"
#include <vector>
#include <string>
#include <memory>

namespace lbm {

class LBESolver {
public:
    LBESolver(const GridSize& grid_size, const SimulationParams& params);
    
    // 删除拷贝构造函数和赋值操作符
    LBESolver(const LBESolver&) = delete;
    LBESolver& operator=(const LBESolver&) = delete;
    
    // 移动构造函数和赋值操作符
    LBESolver(LBESolver&&) noexcept = default;
    LBESolver& operator=(LBESolver&&) noexcept = default;
    
    ~LBESolver() = default;

    // 主求解函数
    void solve();
    
    // 设置边界条件（使用浸没边界方法）
    void setupBoundary(ImmersedBoundary& ib);

    // 清除边界条件
    void clearBoundary();
    
    // 获取结果
    const ScalarField& getDensity() const { return rho_; }
    const ScalarField& getVelocityX() const { return ux_; }
    const ScalarField& getVelocityY() const { return uy_; }
    
    // 输出到CSV（可选）
    void outputToCSV(const std::string& filename) const;

private:
    // 网格和参数
    GridSize grid_size_;
    SimulationParams params_;
    
    // 物理场
    ScalarField rho_;        // 密度
    ScalarField ux_;         // x方向速度
    ScalarField uy_;         // y方向速度
    ScalarField u_star_x_;   // 中间x速度
    ScalarField u_star_y_;   // 中间y速度
    ScalarField delta_u_x_;  // 速度修正x分量
    ScalarField delta_u_y_;  // 速度修正y分量
    
    // 分布函数
    DistributionFunction f_;       // 分布函数
    DistributionFunction feq_;     // 平衡分布函数
    DistributionFunction f_temp_;  // 临时分布函数
    
    // 浸没边界方法
    std::unique_ptr<ImmersedBoundary> ib_method_;
    
    // VTK写入器
    std::unique_ptr<VTKWriter> vtk_writer_;
    
    // 私有成员函数
    void initializeEquilibrium();
    void computeEquilibrium();
    void collisionStep();
    void streamingStep();
    void computeMacroscopic();
    void computeBoundaryVelocityCorrection();
    void computeForceDensity();
    void outputVTK(int step);
    void applyInletBoundary(double inlet_velocity);
    void applyOutletBoundary();
};

} // namespace lbm
