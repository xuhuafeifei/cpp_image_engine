//
// Created by 25080 on 2025/2/6.
//

/**
 * 绘制随机地面
 */

#ifndef RMRENDERER_RENDER3_H
#define RMRENDERER_RENDER3_H

#pragma once

#include "common.h"
#include "camera2.h"

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;

const glm::vec4 circle = glm::vec4({0, 0, 7, 2});
const float y_offset = 0.2;

glm::vec3 getCamera()
{
    return glm::vec3(0., 0 - y_offset, 0);
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

glm::vec3 noise(glm::vec2 pos)
{
    auto i = glm::floor(pos);
    auto f = glm::fract(pos);
    // auto u = f * f * (glm::vec2(3.0) - glm::vec2(2.0) * f);
    auto u = glm::smoothstep(0.f, 1.f, f);
    auto du = glm::vec2(6.) * u * (glm::vec2(1.) - u);

    auto a = random(i);
    auto b = random(i + glm::vec2(1, 0));
    auto c = random(i + glm::vec2(0, 1));
    auto d = random(i + glm::vec2(1, 1));

//    return glm::mix(a, b, u.x) + (c - a) * u.y * (1.f - u.x) + (d - b) * u.x * u.y;
    return glm::vec3(
            a + (b - a) * u.x * (1. - u.y) + (c - a) * (1. - u.x) * u.y + (d - a) * u.x * u.y,
            du * (glm::vec2(b - a, c - a) + (a - b - c + d) * glm::vec2(u.y, u.x))
    );
}

glm::mat2 mat = glm::mat2(.6, -0.8, 0.8, 0.6);

float ground(glm::vec3 x)
{
    auto p = glm::vec3(0.005) * x;
//    auto p = x;
    float a = 0.;
    float b = 1.;
    glm::vec2 d = glm::vec2(0);

    for (int i = 0; i < 8; ++i)
    {
        auto n = noise(p);
        d += glm::vec2(n.y, n.z);
        a += b * n.x / (1. + glm::dot(d, d));
        glm::vec2 p2(p.x, p.y);
        p2 = mat * p2 * glm::vec2(2.);
        p.x = p2.x;
        p.y = p2.y;
        b *= 0.5;
        /*
        // a += noise(p).x;
        // 将 glm::vec3 转换为 glm::vec2
        glm::vec2 p2(p.x, p.y);
        // 使用矩阵与向量的乘法
        p2 = mat * p2;
        // 将变换后的结果重新赋值给 p 的 x 和 y 分量
        p.x = p2.x;
        p.y = p2.y;
         */
    }

    return 120 * a;
}

float ground(glm::vec2 x)
{
    auto p = glm::vec2(0.008) * x;
    float a = 0.;
    float b = 1.;
    glm::vec2 d = glm::vec2(0);

    for (int i = 0; i < 8; ++i)
    {
        auto n = noise(x);
        d += glm::vec2(n.y, n.z);
        a += b * n.x / (1. + glm::dot(d, d));
        glm::vec2 p2(p.x, p.y);
        p2 = mat * p2 * glm::vec2(2.);
        p.x = p2.x;
        p.y = p2.y;
        b *= 0.5;
    }

    return 120 * a;
}

float groundH(glm::vec2 x)
{
    glm::vec2 p = glm::vec2(0.005) * x;
    auto a = 0.;
    auto b = 1.;
    glm::vec2 d = glm::vec2(0);

    for (int i = 0; i < 16; ++i)
    {
        glm::vec3 n = noise(p);
        d += glm::vec2(n.y, n.z);
        a += b * n.x / (1. + glm::dot(d, d));
        p = mat * p * glm::vec2(2.);
        b *= 0.5;
    }

    return 120 * a;
}

float groundH(glm::vec3 x)
{
    return groundH({x.x, x.y});
}

float map(glm::vec3 p)
{
    return groundH({p.x, p.z, 0}) - p.y;
}

float mapGround(glm::vec2 p)
{
    return groundH({p.x, p.y}) - p.y;
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

    for (int i = 0; i < maxMarchTime * 3 && t < tMax; i++)
    {
        glm::vec3 p = ro + rd * t;
        float d = map(p);
        if (d < delta * t)
            break;
        t += .2f * d;
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

glm::vec3 calcNormalGround(glm::vec3 p)
{
    glm::vec2 e = glm::vec2(1e-5, 0);
    return glm::normalize(
            glm::vec3(
                    mapGround({p.x + e.x, p.z + e.y}) - mapGround({p.x - e.x, p.z - e.y}),
                    -2.0 * e.x,
                    mapGround({p.x + e.y, p.z + e.x}) - mapGround({p.x - e.y, p.z - e.x})
                    )
            );
}

Color render(int x, int y)
{
    glm::vec2 uv = fixUV(x, y);

    // set camera
    auto target = glm::vec3(0, 0, 0);
    float camHight = -1.5;
    float camRad = 1.5;
    auto camLoc = glm::vec3 (camRad, camHight, camRad);
    auto camMat = camera(target, camLoc, 0.);

    glm::vec3 color(0);
    glm::vec3 light = glm::vec3(10, -15, 7);

    auto ray = camMat.getRay(uv);

    auto ro = ray.ro;
    auto rd = ray.rd;

    float t = rayMarch(ro, rd);

    if (t < tMax)
    {
        glm::vec3 p = ro + rd * t;
        glm::vec3 n = calcNormalGround(p);
        float diff = glm::dot(
                glm::normalize(light - p),
                n);
        color = glm::vec3(1) * diff;
    }
    return fromVec(color);
}

#endif //RMRENDERER_RENDER3_H
