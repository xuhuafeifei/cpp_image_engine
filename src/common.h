//
// Created by 25080 on 2025/1/27.
//
#pragma once
#ifndef RMRENDERER_COMMON_H
#define RMRENDERER_COMMON_H

#include "../lib/raylib/raylib.h"
#include "../lib/glm/glm.hpp"

typedef glm::vec4 sphere;

float toDegree(float t) {
    return t * 180.0f / PI;
}


// 命名空间Math
// my math
namespace mm {

    // 计算两个浮点数的最大值
    float max(float a, float b) {
        return (a > b) ? a : b;
    }

    // 计算两个浮点数的最小值
    float min(float a, float b) {
        return (a < b) ? a : b;
    }

    // 计算三维向量的长度
    float length(const glm::vec3& v) {
        return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    // 计算两个三维向量对应分量的最大值
    glm::vec3 max(const glm::vec3& a, const glm::vec3& b) {
        return glm::vec3(mm::max(a.x, b.x), mm::max(a.y, b.y), mm::max(a.z, b.z));
    }

    // 计算两个三维向量对应分量的最小值
    glm::vec3 min(const glm::vec3& a, const glm::vec3& b) {
        return glm::vec3(mm::min(a.x, b.x), mm::min(a.y, b.y), mm::min(a.z, b.z));
    }
}

#endif //RMRENDERER_COMMON_H