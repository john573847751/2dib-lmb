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

void LBESolver::setupBoundary(ImmersedBoundary& ib) {
    ib_method_ = std::make_unique<ImmersedBoundary>(std::move(ib));
    ib_method_->initialize(grid_size_);
}

void LBESolver::clearBoundary() {
    ib_method_.reset();
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
    if (ib_method_ && ib_method_->isInitialized()) {
        ib_method_->computeVelocityCorrection(u_star_x_, u_star_y_,
                                            delta_u_x_, delta_u_y_,
                                            grid_size_);
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
    // 自由出口边界条件（Free outflow boundary condition）
    // 使用二阶外推：出口值等于内部节点的外推值
    // f_out = 2*f_inner - f_more_inner (二阶精度)
    const int nx = grid_size_.nx;
    const int ny = grid_size_.ny;

    for (int j = 0; j < ny; j++) {
        // 二阶零梯度条件（出口速度等于相邻内部节点的外推值）
        ux_[nx - 1][j] = 2.0 * ux_[nx - 2][j] - ux_[nx - 3][j];
        uy_[nx - 1][j] = 2.0 * uy_[nx - 2][j] - uy_[nx - 3][j];
        rho_[nx - 1][j] = 2.0 * rho_[nx - 2][j] - rho_[nx - 3][j];

        // 同时更新分布函数（使用外推）
        for (int alpha = 0; alpha < Q; alpha++) {
            // 只对向外流动的方向进行外推
            if (EX[alpha] > 0) {  // 流出方向
                f_[nx - 1][j][alpha] = 2.0 * f_[nx - 2][j][alpha] - f_[nx - 3][j][alpha];
            } else if (EX[alpha] == 0) {
                // 零速度方向，保持不变
            }
            // 对于向内流动的方向，使用平衡态
        }
    }

    // 对右上角和右下角进行角点处理
    // 右边界顶部角点
    int j = ny - 1;
    ux_[nx - 1][j] = 2.0 * ux_[nx - 2][j] - ux_[nx - 3][j];
    uy_[nx - 1][j] = 2.0 * uy_[nx - 2][j] - uy_[nx - 3][j];
    rho_[nx - 1][j] = 2.0 * rho_[nx - 2][j] - rho_[nx - 3][j];

    // 右边界底部角点
    j = 0;
    ux_[nx - 1][j] = 2.0 * ux_[nx - 2][j] - ux_[nx - 3][j];
    uy_[nx - 1][j] = 2.0 * uy_[nx - 2][j] - uy_[nx - 3][j];
    rho_[nx - 1][j] = 2.0 * rho_[nx - 2][j] - rho_[nx - 3][j];
}

void LBESolver::outputVTK(int step) {
    const std::string filename = "solution";
    vtk_writer_->writeLBESolution(filename, rho_, ux_, uy_, grid_size_.nx, grid_size_.ny, step);
}

void LBESolver::solve() {
    // 计算Re=50对应的入口速度
    // 注意：入口速度现在需要在main.cpp中根据具体几何形状计算
    const double viscosity = (params_.tau - 0.5) * CS2;
    const double reynolds_number = 50.0;
    const double characteristic_length = 10.0;  // 默认特征长度
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
