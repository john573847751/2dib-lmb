#pragma once

#include "Types.hpp"
#include <vector>
#include <string>

namespace lbm {

// 前向声明
class ImmersedBoundary;

// 几何形状基类
class Geometry {
public:
    virtual ~Geometry() = default;

    // 获取边界点
    virtual std::vector<BoundaryPoint> getBoundaryPoints() const = 0;

    // 获取几何中心
    virtual double getCenterX() const = 0;
    virtual double getCenterY() const = 0;
};

// 圆形几何形状
class CircleGeometry : public Geometry {
public:
    CircleGeometry(double radius, double center_x, double center_y,
                  const std::vector<double>& velocity = {0.0, 0.0}, int num_points = 100);

    std::vector<BoundaryPoint> getBoundaryPoints() const override;

    double getCenterX() const override { return center_x_; }
    double getCenterY() const override { return center_y_; }
    double getRadius() const { return radius_; }

private:
    double radius_;
    double center_x_;
    double center_y_;
    double velocity_x_;
    double velocity_y_;
    int num_points_;
};

// 浸没边界方法类
class ImmersedBoundary {
public:
    ImmersedBoundary();

    // 添加几何形状
    void addGeometry(std::unique_ptr<Geometry> geometry);

    // 初始化（构建矩阵等）
    void initialize(const GridSize& grid_size);

    // 检查是否已初始化
    bool isInitialized() const { return initialized_; }

    // 获取边界点总数
    int getNumBoundaryPoints() const;

    // 计算边界速度修正
    void computeVelocityCorrection(const ScalarField& u_star_x, const ScalarField& u_star_y,
                                  ScalarField& delta_u_x, ScalarField& delta_u_y,
                                  const GridSize& grid_size) const;

    // 清除所有几何形状
    void clear();

private:
    // Delta函数
    double deltaFunction(double r) const;

    // 插值函数D
    double interpolationD(int i, int j, double x_B, double y_B) const;

    // 构建矩阵A
    void buildMatrixA(const GridSize& grid_size);

    // 计算逆矩阵
    void computeInverseMatrix();

    // 几何形状列表
    std::vector<std::unique_ptr<Geometry>> geometries_;

    // 边界点列表（所有几何形状的边界点合并）
    std::vector<BoundaryPoint> boundary_points_;

    // 矩阵A及其逆矩阵
    std::vector<std::vector<double>> A_matrix_;
    std::vector<std::vector<double>> A_inv_;

    // 初始化标志
    bool initialized_;
};

} // namespace lbm
