#pragma once

#include <array>
#include <cmath>

namespace lbm {

// D2Q9 模型参数
inline constexpr int Q = 9;
inline constexpr int DIM = 2;

// 权重系数
inline constexpr std::array<double, Q> OMEGA = {
    4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};

// 离散速度方向
inline constexpr std::array<int, Q> EX = {0, 1, 0, -1, 0, 1, -1, -1, 1};
inline constexpr std::array<int, Q> EY = {0, 0, 1, 0, -1, 1, 1, -1, -1};

// 物理常数
inline constexpr double CS2 = 1.0 / 3.0;  // 声速平方
inline constexpr double CS4 = CS2 * CS2;

} // namespace lbm
