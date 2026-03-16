#include "ImmersedBoundary.hpp"
#include "Constants.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>

namespace lbm {

// ============ CircleGeometry 实现 ============

CircleGeometry::CircleGeometry(double radius, double center_x, double center_y,
                              const std::vector<double>& velocity, int num_points)
    : radius_(radius), center_x_(center_x), center_y_(center_y),
      velocity_x_(0.0), velocity_y_(0.0), num_points_(num_points) {
    // 参数验证
    if (radius <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }
    if (num_points <= 0) {
        throw std::invalid_argument("Number of points must be positive");
    }
    
    // 安全地设置速度分量
    if (velocity.size() > 0) {
        velocity_x_ = velocity[0];
        if (velocity.size() > 1) {
            velocity_y_ = velocity[1];
        }
    }
}


std::vector<BoundaryPoint> CircleGeometry::getBoundaryPoints() const {
    std::vector<BoundaryPoint> points;
    points.reserve(num_points_);

    const double ds = 2.0 * M_PI * radius_ / num_points_;

    for (int l = 0; l < num_points_; l++) {
        const double angle = 2.0 * M_PI * l / num_points_;
        const double x = center_x_ + radius_ * std::cos(angle);
        const double y = center_y_ + radius_ * std::sin(angle);

        points.emplace_back(x, y, velocity_x_, velocity_y_, ds);
    }

    return points;
}

// ============ ImmersedBoundary 实现 ============

ImmersedBoundary::ImmersedBoundary() : initialized_(false) {}

void ImmersedBoundary::addGeometry(std::unique_ptr<Geometry> geometry) {
    geometries_.push_back(std::move(geometry));
    initialized_ = false;  // 需要重新初始化
}

void ImmersedBoundary::initialize(const GridSize& grid_size) {
    // 合并所有几何形状的边界点
    boundary_points_.clear();
    for (const auto& geometry : geometries_) {
        auto points = geometry->getBoundaryPoints();
        boundary_points_.insert(boundary_points_.end(),
                               std::make_move_iterator(points.begin()),
                               std::make_move_iterator(points.end()));
    }

    std::cout << "ImmersedBoundary initialized with " << boundary_points_.size()
              << " boundary points from " << geometries_.size() << " geometries." << std::endl;

    // 构建矩阵
    buildMatrixA(grid_size);
    computeInverseMatrix();

    initialized_ = true;
}

/**
 * @brief 获取边界点的数量
 * @return 边界点的数量
 */
int ImmersedBoundary::getNumBoundaryPoints() const {
    return static_cast<int>(boundary_points_.size());
}

void ImmersedBoundary::computeVelocityCorrection(const ScalarField& u_star_x,
                                                  const ScalarField& u_star_y,
                                                  ScalarField& delta_u_x,
                                                  ScalarField& delta_u_y,
                                                  const GridSize& grid_size) const {
    if (!initialized_) {
        std::cerr << "Warning: ImmersedBoundary not initialized!" << std::endl;
        return;
    }

    const int n_boundary = static_cast<int>(boundary_points_.size());
    const int b_size = n_boundary * 2;

    // 计算右端项 b = U_B - sum(u* * D)
    std::vector<double> b(b_size, 0.0);

    for (int l = 0; l < n_boundary; l++) {
        double sum_ux = 0.0, sum_uy = 0.0;

        const int x_start = std::max(0, static_cast<int>(boundary_points_[l].x - 3));
        const int x_end = std::min(grid_size.nx, static_cast<int>(boundary_points_[l].x + 4));
        const int y_start = std::max(0, static_cast<int>(boundary_points_[l].y - 3));
        const int y_end = std::min(grid_size.ny, static_cast<int>(boundary_points_[l].y + 4));

        for (int i = x_start; i < x_end; i++) {
            for (int j = y_start; j < y_end; j++) {
                const double D_il = interpolationD(i, j, boundary_points_[l].x, boundary_points_[l].y);
                sum_ux += u_star_x[i][j] * D_il;
                sum_uy += u_star_y[i][j] * D_il;
            }
        }

        b[2 * l] = boundary_points_[l].ux - sum_ux;
        b[2 * l + 1] = boundary_points_[l].uy - sum_uy;
    }

    // 求解线性系统 AX = b
    std::vector<double> X(b_size, 0.0);
    for (int i = 0; i < b_size; i++) {
        for (int j = 0; j < b_size; j++) {
            X[i] += A_inv_[i][j] * b[j];
        }
    }

    // 将边界速度修正应用到欧拉网格
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < grid_size.nx; i++) {
        for (int j = 0; j < grid_size.ny; j++) {
            delta_u_x[i][j] = 0.0;
            delta_u_y[i][j] = 0.0;

            for (int l = 0; l < n_boundary; l++) {
                const double D_ij = interpolationD(i, j, boundary_points_[l].x, boundary_points_[l].y);
                delta_u_x[i][j] += X[2 * l] * D_ij * boundary_points_[l].ds;
                delta_u_y[i][j] += X[2 * l + 1] * D_ij * boundary_points_[l].ds;
            }
        }
    }
}

/**
 * @brief 清除浸入边界方法的所有数据
 * 
 * 重置所有内部状态，包括几何体、边界点、矩阵数据和初始化标志。
 * 调用此函数后，需要重新初始化才能进行计算。
 */
void ImmersedBoundary::clear() {
    geometries_.clear();
    boundary_points_.clear();
    A_matrix_.clear();
    A_inv_.clear();
    initialized_ = false;
}

/**
 * @brief 计算浸入边界方法中的 delta 函数值
 * 
 * 该函数实现了一个余弦型 delta 函数，用于在浸入边界方法中
 * 将边界力从拉格朗日网格插值到欧拉网格。
 * 
 * @param r 距离参数（无量纲）
 * @return double delta 函数值，当 |r| > 2 时返回 0
 */
double ImmersedBoundary::deltaFunction(double r) const {
    const double abs_r = std::abs(r);
    if (abs_r <= 2.0) {
        return 0.25 * (1.0 + std::cos(M_PI * abs_r / 2.0));
    }
    return 0.0;
}

/**
 * @brief 计算离散 delta 函数的插值权重
 * @param i 网格点的 x 方向索引
 * @param j 网格点的 y 方向索引
 * @param x_B 边界点的 x 坐标
 * @param y_B 边界点的 y 坐标
 * @return 插值权重值，为 x 和 y 方向 delta 函数的乘积
 */
double ImmersedBoundary::interpolationD(int i, int j, double x_B, double y_B) const {
    const double dx = static_cast<double>(i) - x_B;
    const double dy = static_cast<double>(j) - y_B;
    return deltaFunction(dx) * deltaFunction(dy);
}


void ImmersedBoundary::buildMatrixA(const GridSize& grid_size) {
    const int n_boundary = static_cast<int>(boundary_points_.size());
    const int matrix_size = n_boundary * 2;

    A_matrix_.assign(matrix_size, std::vector<double>(matrix_size, 0.0));

    for (int l = 0; l < n_boundary; l++) {
        for (int m = 0; m < n_boundary; m++) {
            double sum_x = 0.0, sum_y = 0.0;

            const int x_start = std::max(0, static_cast<int>(boundary_points_[l].x - 3));
            const int x_end = std::min(grid_size.nx, static_cast<int>(boundary_points_[l].x + 4));
            const int y_start = std::max(0, static_cast<int>(boundary_points_[l].y - 3));
            const int y_end = std::min(grid_size.ny, static_cast<int>(boundary_points_[l].y + 4));

            for (int i = x_start; i < x_end; i++) {
                for (int j = y_start; j < y_end; j++) {
                    const double D_il = interpolationD(i, j, boundary_points_[l].x, boundary_points_[l].y);
                    const double D_im = interpolationD(i, j, boundary_points_[m].x, boundary_points_[m].y);

                    sum_x += D_il * D_im * boundary_points_[m].ds;
                    sum_y += D_il * D_im * boundary_points_[m].ds;
                }
            }

            // A矩阵的x分量
            A_matrix_[2 * l][2 * m] = sum_x;
            A_matrix_[2 * l][2 * m + 1] = 0.0;
            // A矩阵的y分量
            A_matrix_[2 * l + 1][2 * m] = 0.0;
            A_matrix_[2 * l + 1][2 * m + 1] = sum_y;
        }
    }
}

void ImmersedBoundary::computeInverseMatrix() {
    const int n = static_cast<int>(A_matrix_.size());
    A_inv_.assign(n, std::vector<double>(n, 0.0));

    // 单位矩阵
    for (int i = 0; i < n; i++) {
        A_inv_[i][i] = 1.0;
    }

    // 复制A矩阵用于高斯-约旦消元
    std::vector<std::vector<double>> A = A_matrix_;

    // 高斯-约旦消元法
    for (int i = 0; i < n; i++) {
        // 寻找主元
        int max_row = i;
        double max_val = std::abs(A[i][i]);

        for (int k = i + 1; k < n; k++) {
            const double val = std::abs(A[k][i]);
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }

        // 交换行
        if (max_row != i) {
            std::swap(A[i], A[max_row]);
            std::swap(A_inv_[i], A_inv_[max_row]);
        }

        // 检查奇异矩阵
        if (std::abs(A[i][i]) < std::numeric_limits<double>::epsilon()) {
            std::cerr << "Error: Matrix is singular!" << std::endl;
            return;
        }

        // 归一化
        const double pivot = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= pivot;
            A_inv_[i][j] /= pivot;
        }

        // 消元
        for (int k = 0; k < n; k++) {
            if (k != i) {
                const double factor = A[k][i];
                for (int j = 0; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                    A_inv_[k][j] -= factor * A_inv_[i][j];
                }
            }
        }
    }
}

} // namespace lbm
