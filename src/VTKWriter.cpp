#include "VTKWriter.hpp"
#include <filesystem>

namespace fs = std::filesystem;

namespace lbm {

void VTKWriter::setOutputDirectory(const std::string& dir) {
    output_dir_ = dir;
    ensureDirectoryExists(output_dir_);
}

bool VTKWriter::ensureDirectoryExists(const std::string& path) {
    try {
        if (!fs::exists(path)) {
            fs::create_directories(path);
        }
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error creating directory: " << e.what() << std::endl;
        return false;
    }
}

std::string VTKWriter::generateFilename(const std::string& base_name, int step) {
    std::ostringstream oss;
    oss << output_dir_ << "/" << base_name << "_" 
        << std::setw(6) << std::setfill('0') << step << ".vtk";
    return oss.str();
}

template<typename T>
bool VTKWriter::writeScalarField2D(const std::string& filename, 
                                   const std::vector<std::vector<T>>& data,
                                   int nx, int ny,
                                   const std::string& field_name) {
    std::string full_path = generateFilename(filename, 0);
    std::ofstream file(full_path);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << full_path << std::endl;
        return false;
    }
    
    // VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "LBE Scalar Field Output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "POINTS " << nx * ny << " float\n";
    
    // 写入网格点坐标
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << i << " " << j << " 0.0\n";
        }
    }
    
    // 写入标量数据
    file << "POINT_DATA " << nx * ny << "\n";
    file << "SCALARS " << field_name << " float 1\n";
    file << "LOOKUP_TABLE default\n";
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << data[i][j] << "\n";
        }
    }
    
    file.close();
    return true;
}

template<typename T>
bool VTKWriter::writeVectorField2D(const std::string& filename,
                                   const std::vector<std::vector<T>>& ux,
                                   const std::vector<std::vector<T>>& uy,
                                   int nx, int ny,
                                   const std::string& field_name) {
    std::string full_path = generateFilename(filename, 0);
    std::ofstream file(full_path);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << full_path << std::endl;
        return false;
    }
    
    // VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "LBE Vector Field Output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "POINTS " << nx * ny << " float\n";
    
    // 写入网格点坐标
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << i << " " << j << " 0.0\n";
        }
    }
    
    // 写入向量数据
    file << "POINT_DATA " << nx * ny << "\n";
    file << "VECTORS " << field_name << " float\n";
    
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << ux[i][j] << " " << uy[i][j] << " 0.0\n";
        }
    }
    
    file.close();
    return true;
}

template<typename T>
bool VTKWriter::writeLBESolution(const std::string& filename,
                                 const std::vector<std::vector<T>>& rho,
                                 const std::vector<std::vector<T>>& ux,
                                 const std::vector<std::vector<T>>& uy,
                                 int nx, int ny,
                                 int step) {
    std::string full_path = generateFilename(filename, step);
    std::ofstream file(full_path);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << full_path << std::endl;
        return false;
    }
    
    // VTK文件头
    file << "# vtk DataFile Version 3.0\n";
    file << "LBE Simulation Solution\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_GRID\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "POINTS " << nx * ny << " float\n";
    
    // 写入网格点坐标
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << i << " " << j << " 0.0\n";
        }
    }
    
    // 写入数据
    file << "POINT_DATA " << nx * ny << "\n";
    
    // 密度场
    file << "SCALARS density float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << rho[i][j] << "\n";
        }
    }
    
    // 速度场
    file << "VECTORS velocity float\n";
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            file << ux[i][j] << " " << uy[i][j] << " 0.0\n";
        }
    }
    
    file.close();
    return true;
}

// 显式实例化
template bool VTKWriter::writeScalarField2D<double>(
    const std::string&, const std::vector<std::vector<double>>&, int, int, const std::string&);
    
template bool VTKWriter::writeVectorField2D<double>(
    const std::string&, const std::vector<std::vector<double>>&,
    const std::vector<std::vector<double>>&, int, int, const std::string&);
    
template bool VTKWriter::writeLBESolution<double>(
    const std::string&, const std::vector<std::vector<double>>&,
    const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&,
    int, int, int);

} // namespace lbm
