#include "LBESolver.hpp"
#include <iostream>
#include <memory>

int main(int argc, char* argv[]) {
    try {
    // 网格尺寸
    const int nx = 300;
    const int ny = 100;
    lbm::GridSize grid_size(nx, ny);
        
        // 模拟参数
        const double tau = 0.6;      // 弛豫时间
        const double dt = 1.0;       // 时间步长
        const int num_steps = 1000;  // 总步数
        const int output_interval = 200;  // VTK输出间隔
        
        lbm::SimulationParams params(tau, dt, num_steps, output_interval);
        
        // 创建LBE求解器
        auto solver = std::make_unique<lbm::LBESolver>(grid_size, params);
        
        std::cout << "Starting LBE simulation..." << std::endl;
        std::cout << "Grid size: " << nx << " x " << ny << std::endl;
        std::cout << "Total steps: " << num_steps << std::endl;
        std::cout << "Output interval: " << output_interval << " steps" << std::endl;
        std::cout << "Relaxation time: " << tau << std::endl;
        std::cout << "=======================================" << std::endl;
        
        // 运行求解
        solver->solve();
        
        // 输出最终结果到CSV
        solver->outputToCSV("final_results.csv");
        
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
