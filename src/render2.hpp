//
// Created by 25080 on 2025/2/4.
//

/**
 * 绘制天空
 */
#ifndef RMRENDERER_RENDER2_H
#define RMRENDERER_RENDER2_H

#pragma once

//#include "raylib.h"
//#include "glm/glm.hpp"
#include "common.h"
#include "camera2.h"

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
float smooth_func(float f)
{
    // a控制锐边程度, a越大, 锐边越明显
    const float a = 0.42f;
    const float c = 0.f;
    return  a * a * f + c;
}

/**
 * 分段函数
 */
float seg_water_ground_func(float y, float p)
{
    float res;
    float c = 0.2;
    if (y < p)
    {
        res = c / (1 + p - y);
    }
    else
    {
        // exp用于控制海面颜色, exp越大, 海面越深
        float exp = 10;
        res = glm::pow(y - p, exp) + c;
    }
    return res;
}

/**
 * 增加海平面
 */
float water_ground(float y)
{
    return  seg_water_ground_func(y, 0.25);
}

glm::vec3 sky_draw_my(glm::vec2 uv)
{
    float y = glm::clamp(std::abs(uv.y), 0.f, 0.5f);
    // 背景色
    glm::vec3 color = glm::vec3(150./255, 195./255, 247./255);
    // 添加太阳(x, y)
    glm::vec2 sun = glm::vec2(-1., -1);
    // 到太阳的距离
    float dist = glm::length(sun - uv);
    // 距离越小, 越亮
    // 通过pow函数, 让反比例函数迅速衰减, 大的值迅速膨胀, 小的值迅速衰减, 抑制太阳光晕
    // alpha控制衰减速率. alpha越小, 越亮
    float alpha = 1.1;
    glm::vec3 sunLight = glm::vec3(glm::pow(1 / (dist + 1 * alpha), 2));
    // 增加渐变(中间亮, 上下暗)
    color = color + glm::vec3(sunLight);
    // 增加环境光(整体调暗)
    float k = 0.50; // k越大, 越暗
    color = color - glm::vec3(k);
    // 增加y轴渐变
    auto c = smooth_func(y);
    color = color - c;
    color += water_ground(uv.y);
//    color += glm::vec3(water_ground(uv.y), 0., -0.12);
    return color;
}

glm::mat2 mat = glm::mat2(.6, -0.8, 0.8, 0.6);

float random(float x)
{
    return glm::fract(23345.6754 * sin(x));
}

/**
 * 原本的sin函数会生成无理数, 在某些情况下会不连续
 */
float hash12(glm::vec2 p)
{
    glm::vec3 p3 = glm::fract(glm::vec3(p.x, p.y, 0) * glm::vec3(.1031));
    p3 += glm::dot(p3, glm::vec3(p3.y, p3.z, p3.x) + glm::vec3(33.33));
    return glm::fract((p3.x + p3.y) * p3.z);
}

float random(glm::vec2 pos)
{

    /*
    auto f = 5820.34280 * sin(glm::dot(pos, glm::vec2(24.432, 409.43)));
    f = glm::fract(f);
    return glm::abs(f);
     */
    return hash12(pos);
}

//glm::vec3 noise(glm::vec2 pos)
//{
//    auto i = glm::floor(pos);
//    auto f = glm::fract(pos);
//    auto u = glm::smoothstep(0.f, 1.f, f);
//
//    auto a = random(i);
//    auto b = random(i + glm::vec2(1, 0));
//    auto c = random(i + glm::vec2(0, 1));
//    auto d = random(i + glm::vec2(1, 1));
//
//    auto u_x = u.x;
//    auto u_y = u.y;
//    auto u_x_inv = 1.0f - u_x;
//    auto u_y_inv = 1.0f - u_y;
//
//    auto x1 = a * u_x_inv + b * u_x;
//    auto x2 = c * u_x_inv + d * u_x;
//    auto y1 = x1 * u_y_inv + x2 * u_y;
//
//    return glm::vec3(y1);
//}
//
//float fbm(glm::vec2 p)
//{
////    auto a = 0.f;
////    auto fac = 0.45f;
////    for (int i = 0; i < 4; ++i)
////    {
////        a += fac * noise(p).x;
////        p = glm::vec2(1.5) * mat * p;
////        fac *= 0.12;
////    }
////    return a;
//    return noise(p).x;
//}

glm::vec3 noise(glm::vec2 pos)
{
    auto i = glm::floor(pos);
    auto f = glm::fract(pos);
    auto u = glm::smoothstep(0.f, 1.f, f);

    auto a = random(i);
    auto b = random(i + glm::vec2(1, 0));
    auto c = random(i + glm::vec2(0, 1));
    auto d = random(i + glm::vec2(1, 1));

    auto u_x = u.x;
    auto u_y = u.y;
    auto u_x_inv = 1.0f - u_x;
    auto u_y_inv = 1.0f - u_y;

    auto x1 = a * u_x_inv + b * u_x;
    auto x2 = c * u_x_inv + d * u_x;
    auto y1 = x1 * u_y_inv + x2 * u_y;

    // 使用余弦插值来代替线性插值，获得更平滑的过渡
    return glm::vec3(y1);
}

float fbm(glm::vec2 p)
{
    auto a = 0.f;
    auto fac = 0.45f;
    for (int i = 0; i < 4; ++i)
    {
        a += fac * noise(p).x;
        p = glm::vec2(1.5) * mat * p;
        fac *= 0.12f;  // 控制每一层噪声的衰减
    }
    return a;
}


Color render(int x, int y)
{
//    // Convert pixel coordinates to UV coordinates
//
//    // Camera position and ray direction
//    glm::vec3 ro = getCamera();
//    glm::vec3 rd = glm::normalize(glm::vec3(uv, 0) - ro);
    glm::vec2 uv = fixUV(x, y);
    // set camera
    auto target = glm::vec3(0, 0, 20);
    float camHight = -1.;
    float camRad = 1.5;
    auto camLoc = glm::vec3 (camRad, camHight, camRad);
    auto camMat = camera(target, camLoc, 0.);

    auto ray = camMat.getRay(uv);

    auto ro = ray.ro;
    auto rd = ray.rd;


    // auto color = sky_draw(rd);
    auto color = sky_draw_my(uv);

    // Ray marching for scene objects
    glm::vec3 light = glm::vec3(10, -15, -5); // Light position
    float t = rayMarch(ro, rd); // Ray march to find intersection

    if (t < tMax) // If the ray hits an object
    {
//        glm::vec3 p = ro + rd * t; // Intersection point
//        glm::vec3 n = calcNormal(p); // Surface normal
//        float diff = glm::dot(glm::normalize(light - p), n); // Diffuse lighting
//        diff *= softshadow(p, glm::normalize(light - p), 0.01, 100, 0.15); // Soft shadows
//        color = glm::vec3(1) * diff; // Apply lighting to object
//
//        // Ambient occlusion (AO)
//        float ao = calcAO(p, n); // Calculate AO
//        ao = calcAO2(n); // Optional: Additional AO calculation
//        color *= ao; // Apply AO to object color
        glm::vec3 p = ro + rd * t;
        glm::vec3 n = calcNormal(p);
        color = glm::vec3(0.67, 0.57, 0.44);

        float diff = glm::dot(
                glm::normalize(light - p),
                n);
        diff = glm::clamp(diff, 0.f, 1.f);
        auto lin = glm::vec3(0.);
        // 软阴影
        auto sh = softshadow(p, glm::normalize(light - p), 0.01, 100, 0.15); // Soft shadows

        lin += glm::vec3(diff * 1.3) * glm::vec3(sh, sh * sh * 0.5 + 0.5 * sh, sh * sh * 0.8 + 0.2 * sh);

        color *= lin;
        // 绘制雾气
//        color = fog_draw(color, t);
    }
    else {
        // clouds
        float cloudH = 100.;
//        auto step = glm::clamp(glm::vec3(std::abs((cloudH + ro.y) / rd.y)), 250.f, 450.f);
        auto step = glm::vec3(std::abs((cloudH + ro.y) / rd.y));
        glm::vec3 cloudUV = ro + step * rd;
        auto fb = fbm(glm::vec2(cloudUV.x, cloudUV.z) * glm::vec2(0.01));
//        std::cout << cloudUV.x << " " << cloudUV.y << " " << cloudUV.z << " " << fb << std::endl;
//        color = glm::mix(color, glm::vec3(1., 0.95, 1.), glm::smoothstep(0.2f, 0.8f, fb));
        color = glm::mix(color, glm::vec3(1., 0.95, 1.), fb);
    }

    // Convert final color to output format
    return fromVec(color);
}


#endif //RMRENDERER_RENDER2_H
