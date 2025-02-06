/**
 * 绘制函数
 */
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
    return glm::vec3(0., 0, 0);
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

float ground(glm::vec3 p)
{
    return 2 * sin(p.x) * sin(p.y);
}

/**
 * 绘制直线
 */
float segment(glm::vec3 p, glm::vec3 a, glm::vec3 b, float w)
{
    auto ap = p - a;
    auto ab = b - a;
    auto cos_theta = glm::abs(glm::dot(ap, ab)) / glm::abs(glm::dot(ap, ap));
    auto aq_length = glm::length(ap) * cos_theta;
    auto aq = ab * glm::vec3(aq_length);
    // 点到直线的距离
    auto dist = glm::length(aq - ap);
    // 点在线上
    if (dist <= w) {
        return 1.;
    }
    return 0.;
}

float segment(glm::vec2 p, glm::vec2 a, glm::vec2 b, float w)
{
    auto ap = p - a;
    auto ab = b - a;
    auto cos_theta = glm::dot(ap, ab) / glm::dot(ap, ap);
    auto aq_length = glm::length(ab) * cos_theta;
    auto aq = ab * glm::vec2(aq_length);
    // 点到直线的距离
    auto dist = glm::length(aq - ap);
    // 点在线上
    if (dist < w) {
        return 1.;
    }
    return 0.;
}

float func(float x)
{
//    return 0.25f * sin((PI *  4) * x);
    return 0.75 * glm::smoothstep(0.f, 1.f, x) - 0.25;
}

float plotFunc1(glm::vec2 uv, int x)
{
    float f = 0.;
    auto fx = fixUV(x, 0.).x;
    auto nfx = fixUV(x + 1., 0.).x;
    f += segment(uv, glm::vec2(fx, func(fx)), glm::vec2(nfx, func(nfx)), 0.01);
    return f;
}

float plotFunc2(glm::vec2 uv, int x)
{
    float f = 0.;
    // 判断uv到整个范围内的所有线段的距离
    for (float x = 0.; x <= width; x += 1.)
    {
        auto fx = fixUV(x, 0.).x;
        auto nfx = fixUV(x + 1., 0.).x;
        f += segment(uv, glm::vec2(fx, func(fx)), glm::vec2(nfx, func(nfx)), 0.01);
    }
    // 累加, 放置出现超过1的情况
    return glm::clamp(f, 0.f, 1.f);
}

float funcPlot(glm::vec2 uv, int x)
{
    return plotFunc1(uv, x);
}

float map(glm::vec3 p)
{
    return func(p.x) - p.y;
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

Color render(int x, int y)
{
    glm::vec2 uv = fixUV(x, y);

    // auto color = glm::vec3(segment(uv, glm::vec2(0, 0), glm::vec2(2, 2), 0.5));
    auto color = glm::vec3(funcPlot(uv, x));

    return fromVec(color);
}