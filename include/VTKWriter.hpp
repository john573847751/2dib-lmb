#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <filesystem>

namespace lbm {

// 简单的VTK写入器 - 输出为ASCII格式
class VTKWriter {
public:
    VTKWriter() = default;
    
    // 设置输出目录
    void setOutputDirectory(const std::string& dir);
    
    // 写入二维标量场
    template<typename T>
    bool writeScalarField2D(const std::string& filename, 
                            const std::vector<std::vector<T>>& data,
                            int nx, int ny,
                            const std::string& field_name = "scalar");
    
    // 写入二维向量场
    template<typename T>
    bool writeVectorField2D(const std::string& filename,
                            const std::vector<std::vector<T>>& ux,
                            const std::vector<std::vector<T>>& uy,
                            int nx, int ny,
                            const std::string& field_name = "velocity");
    
    // 写入完整的LBE数据（密度+速度）
    template<typename T>
    bool writeLBESolution(const std::string& filename,
                          const std::vector<std::vector<T>>& rho,
                          const std::vector<std::vector<T>>& ux,
                          const std::vector<std::vector<T>>& uy,
                          int nx, int ny,
                          int step = 0);

private:
    std::string output_dir_ = ".";
    
    // 确保目录存在
    bool ensureDirectoryExists(const std::string& path);
    
    // 生成带时间步的文件名
    std::string generateFilename(const std::string& base_name, int step);
};

} // namespace lbm
