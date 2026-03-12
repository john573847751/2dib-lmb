#include "LBESolver.hpp"
#include "Constants.hpp"
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <filesystem>

namespace lbm {

LBESolver::LBESolver(const GridSize& grid_size, const SimulationParams& params)
    : grid_size_(grid_size), params_(params) {
    
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    // 初始化数组
    rho_ = ScalarField(nx, std::vector<double>(ny, 1.0));
    ux_ = ScalarField(nx, std::vector<double>(ny, 0.0));
    uy_ = ScalarField(nx, std::vector<double>(ny, 0.0));
    u_star_x_ = ScalarField(nx, std::vector<double>(ny, 0.0));
    u_star_y_ = ScalarField(nx, std::vector<double>(ny, 0.0));
    delta_u_x_ = ScalarField(nx, std::vector<double>(ny, 0.0));
    delta_u_y_ = ScalarField(nx, std::vector<double>(ny, 0.0));
    
    f_ = DistributionFunction(nx, ScalarField(ny, std::vector<double>(Q, 0.0)));
    feq_ = DistributionFunction(nx, ScalarField(ny, std::vector<double>(Q, 0.0)));
    f_temp_ = DistributionFunction(nx, ScalarField(ny, std::vector<double>(Q, 0.0)));
    
    // 初始化VTK写入器
    vtk_writer_ = std::make_unique<VTKWriter>();
    // 使用绝对路径确保目录创建在正确位置
    const std::string vtk_dir = std::filesystem::absolute("../vtk_output").string();
    vtk_writer_->setOutputDirectory(vtk_dir);
    
    // 初始化分布函数为平衡态
    initializeEquilibrium();
}

void LBESolver::initializeEquilibrium() {
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const double rho_val = 1.0;
            const double ux_val = 0.0;
            const double uy_val = 0.0;
            
            for (int alpha = 0; alpha < Q; alpha++) {
                const double eu = EX[alpha] * ux_val + EY[alpha] * uy_val;
                const double uu = ux_val * ux_val + uy_val * uy_val;
                
                feq_[i][j][alpha] = OMEGA[alpha] * rho_val *
                                  (1.0 + eu / CS2 + 0.5 * (eu * eu) / (CS2 * CS2) - 0.5 * uu / CS2);
                f_[i][j][alpha] = feq_[i][j][alpha];
            }
        }
    }
}

void LBESolver::setupCircularBoundary(double radius, double center_x, double center_y,
                                      const std::vector<double>& boundary_velocity) {
    const int num_points = 100;
    boundary_points_.clear();
    boundary_points_.reserve(num_points);
    
    const double ds = 2.0 * M_PI * radius / num_points;
    
    for (int l = 0; l < num_points; l++) {
        const double angle = 2.0 * M_PI * l / num_points;
        const double x = center_x + radius * std::cos(angle);
        const double y = center_y + radius * std::sin(angle);
        
        boundary_points_.emplace_back(x, y, boundary_velocity[0], boundary_velocity[1], ds);
    }
    
    buildMatrixA();
    computeInverseMatrix();
}

double LBESolver::deltaFunction(double r) const {
    const double abs_r = std::abs(r);
    if (abs_r <= 2.0) {
        return 0.25 * (1.0 + std::cos(M_PI * abs_r / 2.0));
    }
    return 0.0;
}

double LBESolver::interpolationD(int i, int j, double x_B, double y_B) const {
    const double dx = static_cast<double>(i) - x_B;
    const double dy = static_cast<double>(j) - y_B;
    return deltaFunction(dx) * deltaFunction(dy);
}

void LBESolver::buildMatrixA() {
    const int n_boundary = static_cast<int>(boundary_points_.size());
    const int matrix_size = n_boundary * 2;
    
    A_matrix_.assign(matrix_size, std::vector<double>(matrix_size, 0.0));
    
    for (int l = 0; l < n_boundary; l++) {
        for (int m = 0; m < n_boundary; m++) {
            double sum_x = 0.0, sum_y = 0.0;
            
            const int x_start = std::max(0, static_cast<int>(boundary_points_[l].x - 3));
            const int x_end = std::min(grid_size_.nx, static_cast<int>(boundary_points_[l].x + 4));
            const int y_start = std::max(0, static_cast<int>(boundary_points_[l].y - 3));
            const int y_end = std::min(grid_size_.ny, static_cast<int>(boundary_points_[l].y + 4));
            
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

void LBESolver::computeInverseMatrix() {
    const int n = static_cast<int>(A_matrix_.size());
    A_inv_.assign(n, std::vector<double>(n, 0.0));
    
    // 单位矩阵
    for (int i = 0; i < n; i++) {
        A_inv_[i][i] = 1.0;
    }
    
    // 高斯-约旦消元法
    for (int i = 0; i < n; i++) {
        // 寻找主元
        int max_row = i;
        double max_val = std::abs(A_matrix_[i][i]);
        
        for (int k = i + 1; k < n; k++) {
            const double val = std::abs(A_matrix_[k][i]);
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }
        
        // 交换行
        if (max_row != i) {
            std::swap(A_matrix_[i], A_matrix_[max_row]);
            std::swap(A_inv_[i], A_inv_[max_row]);
        }
        
        // 检查奇异矩阵
        if (std::abs(A_matrix_[i][i]) < std::numeric_limits<double>::epsilon()) {
            std::cerr << "Error: Matrix is singular!" << std::endl;
            return;
        }
        
        // 归一化
        const double pivot = A_matrix_[i][i];
        for (int j = 0; j < n; j++) {
            A_matrix_[i][j] /= pivot;
            A_inv_[i][j] /= pivot;
        }
        
        // 消元
        for (int k = 0; k < n; k++) {
            if (k != i) {
                const double factor = A_matrix_[k][i];
                for (int j = 0; j < n; j++) {
                    A_matrix_[k][j] -= factor * A_matrix_[i][j];
                    A_inv_[k][j] -= factor * A_inv_[i][j];
                }
            }
        }
    }
}

void LBESolver::computeEquilibrium() {
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const double rho_val = rho_[i][j];
            const double ux_val = ux_[i][j];
            const double uy_val = uy_[i][j];
            
            for (int alpha = 0; alpha < Q; alpha++) {
                const double eu = EX[alpha] * ux_val + EY[alpha] * uy_val;
                const double uu = ux_val * ux_val + uy_val * uy_val;
                
                feq_[i][j][alpha] = OMEGA[alpha] * rho_val *
                                  (1.0 + eu / CS2 + 0.5 * (eu * eu) / (CS2 * CS2) - 0.5 * uu / CS2);
            }
        }
    }
}

void LBESolver::collisionStep() {
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    const double inv_tau = 1.0 / params_.tau;
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int alpha = 0; alpha < Q; alpha++) {
                f_temp_[i][j][alpha] = f_[i][j][alpha] -
                                     inv_tau * (f_[i][j][alpha] - feq_[i][j][alpha]);
            }
        }
    }
    
    // 更新分布函数
    std::swap(f_, f_temp_);
}

void LBESolver::streamingStep() {
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    // 创建临时数组
    DistributionFunction f_new = f_;
    
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            for (int alpha = 0; alpha < Q; alpha++) {
                const int ip = i - EX[alpha];
                const int jp = j - EY[alpha];
                
                if (ip >= 0 && ip < nx && jp >= 0 && jp < ny) {
                    f_new[i][j][alpha] = f_[ip][jp][alpha];
                }
            }
        }
    }
    
    // 应用边界条件（反弹边界）
    for (int i = 0; i < nx; i++) {
        for (int alpha = 0; alpha < Q; alpha++) {
            f_new[i][0][alpha] = f_[i][0][alpha];  // 底部边界
            f_new[i][ny - 1][alpha] = f_[i][ny - 1][alpha];  // 顶部边界
        }
    }
    
    for (int j = 0; j < ny; j++) {
        for (int alpha = 0; alpha < Q; alpha++) {
            f_new[0][j][alpha] = f_[0][j][alpha];  // 左边界
            f_new[nx - 1][j][alpha] = f_[nx - 1][j][alpha];  // 右边界
        }
    }
    
    f_ = std::move(f_new);
}

void LBESolver::computeMacroscopic() {
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            // 计算密度
            double rho_val = 0.0;
            double sum_ex_f = 0.0, sum_ey_f = 0.0;
            
            for (int alpha = 0; alpha < Q; alpha++) {
                rho_val += f_[i][j][alpha];
                sum_ex_f += EX[alpha] * f_[i][j][alpha];
                sum_ey_f += EY[alpha] * f_[i][j][alpha];
            }
            
            rho_[i][j] = rho_val;
            u_star_x_[i][j] = sum_ex_f / rho_val;
            u_star_y_[i][j] = sum_ey_f / rho_val;
            
            // 添加速度修正
            ux_[i][j] = u_star_x_[i][j] + delta_u_x_[i][j];
            uy_[i][j] = u_star_y_[i][j] + delta_u_y_[i][j];
        }
    }
}

void LBESolver::computeBoundaryVelocityCorrection() {
    const int n_boundary = static_cast<int>(boundary_points_.size());
    const int b_size = n_boundary * 2;
    
    // 计算右端项 b = U_B - sum(u* * D)
    std::vector<double> b(b_size, 0.0);
    
    for (int l = 0; l < n_boundary; l++) {
        double sum_ux = 0.0, sum_uy = 0.0;
        
        const int x_start = std::max(0, static_cast<int>(boundary_points_[l].x - 3));
        const int x_end = std::min(grid_size_.nx, static_cast<int>(boundary_points_[l].x + 4));
        const int y_start = std::max(0, static_cast<int>(boundary_points_[l].y - 3));
        const int y_end = std::min(grid_size_.ny, static_cast<int>(boundary_points_[l].y + 4));
        
        for (int i = x_start; i < x_end; i++) {
            for (int j = y_start; j < y_end; j++) {
                const double D_il = interpolationD(i, j, boundary_points_[l].x, boundary_points_[l].y);
                sum_ux += u_star_x_[i][j] * D_il;
                sum_uy += u_star_y_[i][j] * D_il;
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
    for (int i = 0; i < grid_size_.nx; i++) {
        for (int j = 0; j < grid_size_.ny; j++) {
            delta_u_x_[i][j] = 0.0;
            delta_u_y_[i][j] = 0.0;
            
            for (int l = 0; l < n_boundary; l++) {
                const double D_ij = interpolationD(i, j, boundary_points_[l].x, boundary_points_[l].y);
                delta_u_x_[i][j] += X[2 * l] * D_ij * boundary_points_[l].ds;
                delta_u_y_[i][j] += X[2 * l + 1] * D_ij * boundary_points_[l].ds;
            }
        }
    }
}

void LBESolver::computeForceDensity() {
    // 这里简化处理，实际应用中需要根据具体物理模型计算
    // F_alpha = (1-1/(2*tau)) * omega_alpha * [(e_alpha-u)/cs2 + (e_alpha*u)/cs4 * e_alpha] * f
}

void LBESolver::applyInletBoundary(double inlet_velocity) {
    // 在左边界(x=0)设置均匀入口速度
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    for (int j = 0; j < ny; j++) {
        // 设置入口速度
        ux_[0][j] = inlet_velocity;
        uy_[0][j] = 0.0;
        rho_[0][j] = 1.0;
        
        // 根据入口速度计算平衡分布函数
        for (int alpha = 0; alpha < Q; alpha++) {
            const double eu = EX[alpha] * inlet_velocity;
            const double uu = inlet_velocity * inlet_velocity;
            
            f_[0][j][alpha] = OMEGA[alpha] * rho_[0][j] *
                             (1.0 + eu / CS2 + 0.5 * (eu * eu) / (CS2 * CS2) - 0.5 * uu / CS2);
        }
    }
}

void LBESolver::applyOutletBoundary() {
    // 在右边界(x=nx-1)设置零梯度出口条件
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;
    
    for (int j = 0; j < ny; j++) {
        // 零梯度条件：出口值等于相邻内部节点值
        rho_[nx - 1][j] = rho_[nx - 2][j];
        ux_[nx - 1][j] = ux_[nx - 2][j];
        uy_[nx - 1][j] = uy_[nx - 2][j];
    }
}

void LBESolver::outputVTK(int step) {
    const std::string filename = "solution";
    vtk_writer_->writeLBESolution(filename, rho_, ux_, uy_, grid_size_.nx, grid_size_.ny, step);
}

void LBESolver::solve() {
    // 设置圆形边界（直径=10，半径=5），位置在(50, 50)
    const double radius = 5.0;
    const double center_x = 50.0;
    const double center_y = 50.0;
    setupCircularBoundary(radius, center_x, center_y);
    
    // 计算Re=50对应的入口速度
    const double viscosity = (params_.tau - 0.5) * CS2;
    const double characteristic_length = 2.0 * radius;  // 直径
    const double reynolds_number = 50.0;
    const double inlet_velocity = reynolds_number * viscosity / characteristic_length;
    
    std::cout << "Reynolds number: " << reynolds_number << std::endl;
    std::cout << "Inlet velocity: " << inlet_velocity << std::endl;
    std::cout << "Viscosity: " << viscosity << std::endl;
    std::cout << "=======================================" << std::endl;
    
    for (int step = 0; step < params_.num_steps; step++) {
        // LBE演化 (初始F_alpha=0)
        computeEquilibrium();
        collisionStep();
        streamingStep();
        computeMacroscopic();
        
        // 应用入口和出口边界条件
        applyInletBoundary(inlet_velocity);
        applyOutletBoundary();
        
        // 计算边界速度修正
        computeBoundaryVelocityCorrection();
        
        // 修正速度并计算力密度
        computeMacroscopic();
        computeForceDensity();
        
        // 输出VTK（每output_interval步）
        if (step % params_.output_interval == 0) {
            outputVTK(step);
            
            // 输出进度
            const int mid_x = grid_size_.nx / 2;
            const int mid_y = grid_size_.ny / 2;
            const double max_vel = std::sqrt(ux_[mid_x][mid_y] * ux_[mid_x][mid_y] +
                                             uy_[mid_x][mid_y] * uy_[mid_x][mid_y]);
            std::cout << "Step: " << step << ", Max velocity: " << max_vel << std::endl;
        }
    }
    
    // 输出最终结果
    outputVTK(params_.num_steps);
}

void LBESolver::outputToCSV(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << "x,y,rho,ux,uy\n";
    
    for (int i = 0; i < grid_size_.nx; i++) {
        for (int j = 0; j < grid_size_.ny; j++) {
            file << i << "," << j << "," << rho_[i][j] << ","
                 << ux_[i][j] << "," << uy_[i][j] << "\n";
        }
    }
    
    file.close();
}

} // namespace lbm
