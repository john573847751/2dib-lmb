#include "LBESolver.hpp"
#include "ImmersedBoundary.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <omp.h>

int main(int argc, char* argv[]) {
    try {
        // 默认线程数
        int num_threads = 4;

        // 解析命令行参数 -t 或 --threads 指定线程数
        for (int i = 1; i < argc; i++) {
            std::string arg = argv[i];
            if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
            } else if (arg == "-h" || arg == "--help") {
                std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
                std::cout << "Options:" << std::endl;
                std::cout << "  -t, --threads N  Set number of OpenMP threads (default: 4)" << std::endl;
                std::cout << "  -h, --help       Show this help message" << std::endl;
                return 0;
            }
        }

        // 设置OpenMP线程数
        omp_set_num_threads(num_threads);
        std::cout << "Using " << num_threads << " OpenMP threads" << std::endl;

        // 网格尺寸
        const int nx = 300;
        const int ny = 100;
        lbm::GridSize grid_size(nx, ny);

        // 模拟参数
        const double tau = 0.6;      // 弛豫时间
        const double dt = 1.0;       // 时间步长
        const int num_steps = 80000;  // 总步数
        const int output_interval = 2000;  // VTK输出间隔

        lbm::SimulationParams params(tau, dt, num_steps, output_interval);

        // 创建LBE求解器
        auto solver = std::make_unique<lbm::LBESolver>(grid_size, params);

        // 创建浸没边界方法
        lbm::ImmersedBoundary ib;

        // 添加多个圆柱体（可在此处添加更多几何形状）
        // 圆柱体1：半径=5，中心在(50, 50)
        std::vector<double> vel1 = {0.0, 0.0};
        auto cylinder1 = std::make_unique<lbm::CircleGeometry>(
            5.0, 50.0, 50.0, vel1, 100);
        ib.addGeometry(std::move(cylinder1));

        // 圆柱体2：半径=3，中心在(100, 30) - 示例：添加第二个圆柱
        // auto cylinder2 = std::make_unique<lbm::CircleGeometry>(
        //     3.0, 100.0, 30.0, {0.0, 0.0}, 80);
        // ib.addGeometry(std::move(cylinder2));

        // 初始化浸没边界
        solver->setupBoundary(ib);

        // 计算Re=50对应的入口速度
        const double viscosity = (params.tau - 0.5) * lbm::CS2;
        const double characteristic_length = 10.0;  // 圆柱直径
        const double reynolds_number = 50.0;
        const double inlet_velocity = reynolds_number * viscosity / characteristic_length;

        std::cout << "Starting LBE simulation..." << std::endl;
        std::cout << "Grid size: " << nx << " x " << ny << std::endl;
        std::cout << "Total steps: " << num_steps << std::endl;
        std::cout << "Output interval: " << output_interval << " steps" << std::endl;
        std::cout << "Relaxation time: " << tau << std::endl;
        std::cout << "Reynolds number: " << reynolds_number << std::endl;
        std::cout << "Inlet velocity: " << inlet_velocity << std::endl;
        std::cout << "Viscosity: " << viscosity << std::endl;
        std::cout << "=======================================" << std::endl;

        // 运行求解
        solver->solve();

        // 输出最终结果到CSV
        //solver->outputToCSV("final_results.csv");

        std::cout << "=======================================" << std::endl;
        std::cout << "LBE simulation completed successfully!" << std::endl;
        std::cout << "VTK files are saved in ./vtk_output/ directory" << std::endl;
        std::cout << "Final results saved to final_results.csv" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred!" << std::endl;
        return 1;
    }
}
