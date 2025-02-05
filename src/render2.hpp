//
// Created by 25080 on 2025/2/4.
//

#ifndef RMRENDERER_RENDER2_H
#define RMRENDERER_RENDER2_H

#pragma once

//#include "raylib.h"
//#include "glm/glm.hpp"
#include "common.h"

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;

// glm::vec3 light = glm::vec3(-10, 15, 2);

const glm::vec4 circle = glm::vec4({0, 0, 7, 2});

glm::vec3 getCamera()
{
    return glm::vec3(0., 0, -1.5);
}

glm::vec2 fixUV(int x, int y)
{
    float m = std::min(width, height);
    float u = (2.0 * x - width) / m;
    float v = (2.0 * y - height) / m;
    return {u, v};
}

unsigned char convert(float v)
{
    if (v < 0)
        return 0;
    if (v > 255)
        return 255;
    return (unsigned char)v;
}

float sdfSphere(glm::vec3 p, glm::vec3 o, float r)
{
    return glm::length(p - o) - r;
}

float sdfGround(glm::vec3 p)
{
    // 很操蛋, 我的这个y轴是向下为正
    // 所以这里应该是地面坐标 - p.y
    return 2 - p.y;
}

// quadratic polynomial
float smin( float a, float b, float k )
{
    float h = glm::clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return glm::mix(b, a, h) - k * h * (1. - h);
}

float map(glm::vec3 p)
{
    // return smin(sdfSphere(p, {circle.x, circle.y, circle.z}, circle.w), sdfGround(p), 1.5);
    return 1e10;
}

// 环境光影响
float calcAO(glm::vec3 p, glm::vec3 n)
{
    float occ = 0.;
    float sca = 1.;
    for (int i = 0; i < 5; ++i)
    {
        float h = 0.01 + 0.03 * float(i);
        float d = map(p + n * h);
        occ += (h - d) * sca;
        sca *= 0.95;
        if (occ > 0.35)
        {
            break;
        }
    }
    return glm::clamp(1. - 3. * occ, 0., 1.);
}

float calcAO2(glm::vec3 n)
{
    // 我的y轴向下为正, 所以需要取符号
    return 0.5 + 0.5 * (-n.y);
}

Color fromVec(glm::vec3 v)
{
    return Color{
            convert(v.x * 255),
            convert(v.y * 255),
            convert(v.z * 255),
            255};
}

float rayMarch(glm::vec3 ro, glm::vec3 rd)
{
    float t = tMin;

    for (int i = 0; i < maxMarchTime && t < tMax; i++)
    {
        glm::vec3 p = ro + rd * t;
        float d = map(p);
        if (d < delta)
            break;
        t += d;
    }

    return t;
}

glm::vec3 calcNormal(glm::vec3 p)
{
    const glm::vec3 v1(1, -1, -1);
    const glm::vec3 v2(-1, -1, 1);
    const glm::vec3 v3(-1, 1, -1);
    const glm::vec3 v4(1);
    return glm::normalize(
            glm::vec3(
                    map(p + v1 * delta) * v1 +
                    map(p + v2 * delta) * v2 +
                    map(p + v3 * delta) * v3 +
                    map(p + v4 * delta) * v4));
}

// imporve: 阴影没有波纹, 但w需要小于0.5, 否则渲染不好
float softshadow( glm::vec3 ro, glm::vec3 rd, float mint, float maxt, float w )
{
    float res = 1.0;
    float ph = 1e20;
    float t = mint;
    for( int i=0; i<256 && t<maxt; i++ )
    {
        auto k = ro + rd*t;
        float h = map(k);
        if( h<0.001 )
            return 0.0;
        float y = h*h/(2.0*ph);
        float d = std::sqrt(h*h-y*y);
        res = std::min( res, d/(w*std::max(0.0f,t-y)) );
        ph = h;
        t += h;
    }
    return res;
}

/**
 * 看陈老师写的, 但不太好嵌入到当前项目当中
 */
glm::vec3 sky_draw(glm::vec3 rd)
{
    // Sky gradient
    float y_ = -rd.y; // Invert y for easier gradient calculation
    glm::vec3 color = glm::vec3(0.3, 0.5, 0.85) - glm::vec3(y_ * y_ * 0.5); // Base sky color
    color = glm::mix(color, glm::vec3(0.85) * glm::vec3(0.7, 0.75, 0.85), std::pow(1. - std::max(y_, 0.0f), 4.)); // Horizon blend

    // Sun direction and glow
    glm::vec3 sunLight = glm::normalize(glm::vec3(0., 0.4, -1.0)); // Sun direction
    float sundot = glm::clamp(glm::dot(rd, sunLight), 0.f, 1.f); // Dot product for sun position

    // Sun glow effect
    color += glm::vec3(0.25) * glm::vec3(1., 0.7, 0.4) * glm::vec3(glm::pow(sundot, 5.0)); // Outer glow
    color += glm::vec3(0.25) * glm::vec3(1., 0.8, 0.6) * glm::vec3(glm::pow(sundot, 64.0)); // Mid glow
    color += glm::vec3(0.2) * glm::vec3(1., 0.8, 0.6) * glm::vec3(glm::pow(sundot, 512.0)); // Inner glow (brightest)

    return color;
}

/**
 * 二次函数
 */
float smooth_func(float f) {
    // a控制锐边程度, a越大, 锐边越明显
    const float a = 0.6f;
    const float c = 0.f;
    return  a * a * f + c;
}

glm::vec3 sky_draw_my(glm::vec2 uv)
{
    float y = glm::clamp(std::abs(uv.y), 0.f, 0.5f);
    // 背景色
    glm::vec3 color = glm::vec3(0.3, 0.5, 0.85);
    // 添加太阳(x, y)
    glm::vec2 sun = glm::vec2(0., 0.);
    // 到太阳的距离
    float dist = glm::length(sun - uv);
    // 距离越小, 越亮
    // 通过pow函数, 让反比例函数迅速衰减, 大的值迅速膨胀, 小的值迅速衰减, 抑制太阳光晕
    // alpha控制衰减速率
    float alpha = 1.2;
    glm::vec3 sunLight = glm::vec3(glm::pow(1 / (dist + alpha), 2));
    // 增加渐变(中间亮, 上下暗)
    color = color + glm::vec3(sunLight);
    // 增加环境光(整体调暗)
    color = color - glm::vec3(0.13);
    // 增加y轴渐变
    auto c = smooth_func(y);
    color = color - c;
    return color;
}

Color render(int x, int y)
{
    // Convert pixel coordinates to UV coordinates
    glm::vec2 uv = fixUV(x, y);

    // Camera position and ray direction
    glm::vec3 ro = getCamera();
    glm::vec3 rd = glm::normalize(glm::vec3(uv, 0) - ro);

    // auto color = sky_draw(rd);
    auto color = sky_draw_my(uv);

    // Ray marching for scene objects
    glm::vec3 light = glm::vec3(10, -15, -5); // Light position
    float t = rayMarch(ro, rd); // Ray march to find intersection

    if (t < tMax) // If the ray hits an object
    {
        glm::vec3 p = ro + rd * t; // Intersection point
        glm::vec3 n = calcNormal(p); // Surface normal
        float diff = glm::dot(glm::normalize(light - p), n); // Diffuse lighting
        diff *= softshadow(p, glm::normalize(light - p), 0.01, 100, 0.15); // Soft shadows
        color = glm::vec3(1) * diff; // Apply lighting to object

        // Ambient occlusion (AO)
        float ao = calcAO(p, n); // Calculate AO
        ao = calcAO2(n); // Optional: Additional AO calculation
        color *= ao; // Apply AO to object color
    }

    // Convert final color to output format
    return fromVec(color);
}


#endif //RMRENDERER_RENDER2_H
