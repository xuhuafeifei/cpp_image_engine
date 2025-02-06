//
// Created by 25080 on 2025/2/6.
//

/**
 * 绘制规则地面
 */

#ifndef RMRENDERER_RENDER3_H
#define RMRENDERER_RENDER3_H

#pragma once

#include "common.h"

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;

const glm::vec4 circle = glm::vec4({0, 0, 7, 2});

glm::vec3 getCamera()
{
    return glm::vec3(0., -5, -1.);
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

float random(glm::vec2 pos)
{

    auto f = 5820.34280 * sin(glm::dot(pos, glm::vec2(24.432, 409.43)));
    f = glm::fract(f);
    return glm::abs(f);
}

float noise(glm::vec2 pos)
{
    auto i = glm::floor(pos);
    auto f = glm::fract(pos);
    auto u = f * f * (glm::vec2(3.0) - glm::vec2(2.0) * f);

    auto a = random(i);
    auto b = random(i + glm::vec2(1, 0));
    auto c = random(i + glm::vec2(0, 1));
    auto d = random(i + glm::vec2(1, 1));

    return glm::mix(a, b, u.x) + (c - a) * u.y * (1.f - u.x) + (d - b) * u.x * u.y;
}

float ground(glm::vec3 p)
{
    return noise(p);
}

float map(glm::vec3 p)
{
    return 3. * ground({p.x, p.z, 0}) - p.y;
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
        if (d < delta)
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

Color render(int x, int y)
{
    glm::vec2 uv = fixUV(x, y);

    glm::vec3 color(0);
    glm::vec3 ro = getCamera();
    glm::vec3 rd = glm::normalize(glm::vec3(uv.x, uv.y - 5, 0) - ro);
    glm::vec3 light = glm::vec3(10, -15, 7);

    float t = rayMarch(ro, rd);

    if (t < tMax)
    {
        glm::vec3 p = ro + rd * t;
        glm::vec3 n = calcNormal(p);
        float diff = glm::dot(
                glm::normalize(light - p),
                n);
        color = glm::vec3(1) * diff;
    }
    return fromVec(color);
}

#endif //RMRENDERER_RENDER3_H
