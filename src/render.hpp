#pragma once

#include "common.h"
#include "camera.h"

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;

const glm::vec4 circle = glm::vec4({0, -2, 7, 2});
const glm::vec3 cir = {circle.x, circle.y, circle.z};

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
    // {距离, tag}
    return glm::length(p - o) - r;
}

float sdfGround(glm::vec3 p)
{
    // 很操蛋, 我的这个y轴是向下为正
    // 所以这里应该是地面坐标 - p.y
    return - p.y;
}

/**
 * p 是交点
 * b 是矩阵的半长度. 也就是分别到x, y, z的半距离
 * center是矩阵的中心位置
 */
float sdBox(glm::vec3 p, glm::vec3 b, glm::vec3 center)
{
    auto q = glm::abs(p - center) - b;
    return mm::length(mm::max(q, glm::vec3(0.0) )) + mm::min(mm::max(q.x, mm::max(q.y,q.z)),0.0);
}

float sdBoxRound(glm::vec3 p, glm::vec3 b, glm::vec3 center, float rad) {
    auto q = glm::abs(p - center) - b;
    return mm::length(mm::max(q, glm::vec3(0.0) )) + mm::min(mm::max(q.x, mm::max(q.y,q.z)),0.0) - rad;
}

float sdBox( glm::vec3 p, glm::vec3 b )
{
    auto q = glm::abs(p) - b;
    return mm::length(mm::max(q, glm::vec3(0.0) )) + mm::min(mm::max(q.x, mm::max(q.y,q.z)),0.0);
}

float shape1(glm::vec3 p)
{
    // union
    auto d = sdfSphere(p, cir, circle.w);
    d = std::min(sdBox(p, {1,1,1}, {0, -4, 7}), d);
    return d;
}

float shape2(glm::vec3 p)
{
    // intersection
    auto d = sdfSphere(p, {5, -3, 6}, 1);
    d = std::max(sdBox(p, {.5,5,.5}, {5, -3.5, 6}), d);
    return d;
}

float shape3(glm::vec3 p)
{
    // difference
    auto d1 = sdfSphere(p, {-5, -3, 6}, 2);
    auto d2 = sdBox(p, {1,1,2}, {-5, -3, 6});
    auto d = std::max(d2 , -1 * d1);
    return d;
}

// quadratic polynomial
float smin( float a, float b, float k )
{
    float h = glm::clamp(0.5 + 0.5 * (b - a) / k, 0., 1.);
    return glm::mix(b, a, h) - k * h * (1. - h);
}

glm::vec2 vec2Min(glm::vec2 a, glm::vec2 b) {
    if (a.x <= b.x)
//    if (smin(a.x, b.x, 1.5))
    {
        return a;
    }
    else
    {
        return b;
    }
}

glm::vec2 map(glm::vec3 p)
{
    auto d = glm::vec2(shape1(p), 2);
    d = vec2Min(d, {shape2(p), 3});
    d = vec2Min(d, {shape3(p), 4});
    // d = vec2Min(d, {sdBox(p, {1,1,1}, {3, 0, 7}), 3});
    return d;
}

Color fromVec(glm::vec3 v)
{
    return Color{
        convert(v.x * 255),
        convert(v.y * 255),
        convert(v.z * 255),
        255};
}

// 初始softshadow, 但渲染时会出现阴影波纹, 因此被废弃
float softshadow_deprecated( glm::vec3 ro, glm::vec3 rd, float mint, float maxt, float k )
{
    float res = 1.0;
    float t = mint;
    for( int i=0; i<256 && t<maxt; i++ )
    {
        float h = map(ro + rd*t).x;
        if( h<0.001 )
            return 0.0;
        res = std::min( res, k*h/t );
        t += h;
    }
    return res;
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
        float h = map(k).x;
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

glm::vec2 rayMarch(glm::vec3 ro, glm::vec3 rd)
{
    float t = 0.1;
    float tmax = 40.;
    glm::vec2 res(-1.);
    // 如果光线方向为正(我的y轴正方向朝下的, md), 那么他一定会与地面相交(没hit other object)的情况下
    // 此时可以手动计算走过多少轮会与地面相交
    if (rd.y >= 0.) {
        float tp = std::abs(ro.y / rd.y);
        tmax = std::min(tp, tmax);
        // 地面flag为1
        // 必须设置为tp, 否则无法渲染阴影
        // todo: 没太懂为啥必须设置为tp
        res = {tp, 1};
    }

    for (int i = 0; i < maxMarchTime && t < tmax; i++)
    {
        glm::vec3 p = ro + rd * t;
        auto m = map(p);
        auto d = m.x;
        // hit object
        if (d < delta)
            return {t, m.y};
        t += d;
    }

    return res;
}

glm::vec3 calcNormal(glm::vec3 p)
{
    const glm::vec3 v1(1, -1, -1);
    const glm::vec3 v2(-1, -1, 1);
    const glm::vec3 v3(-1, 1, -1);
    const glm::vec3 v4(1);
    return glm::normalize(
        glm::vec3(
            map(p + v1 * delta).x * v1 +
            map(p + v2 * delta).x * v2 +
            map(p + v3 * delta).x * v3 +
            map(p + v4 * delta).x * v4));
}

float normalShadow(glm::vec3 p, glm::vec3 light)
{
    auto rm = rayMarch(p, glm::normalize(light - p));
    auto t = rm.x;
    if (t < tMax)
    {
        return 0.1;
    }
    return 1;
}

Color render(int x, int y)
{
    glm::vec2 uv = fixUV(x, y);

    glm::vec3 bg = glm::vec3(0.7, 0.7, 0.9);
    // 增加y轴上的渐变
    glm::vec3 color = bg - glm::normalize(uv).y * glm::vec3(0.13);
    camera cam;
    auto ray = cam.ray(uv);
    glm::vec3 light = glm::vec3(10, -15, -5);

    auto rm = rayMarch(ray.ro, glm::normalize(ray.rd));

    // rm.y = -1时, 表示没有hit object
    if (rm.y > 0)
    {
        auto t = rm.x;
        glm::vec3 p = ray.at(t);
        glm::vec3 n = (rm.y < 1.1) ? glm::vec3(0, -1, 0) : calcNormal(p);
        float diff = glm::dot(
            glm::normalize(light - p),
            n);
        // 绘制阴影
        diff *= softshadow(p, glm::normalize(light - p), 0.1, 100, 0.15);
        // 添加环境光
        float amb = 0.1;
        color = glm::vec3(1) * diff + amb;
        if (std::abs(rm.y - 1) < 0.1) {
            color += glm::vec3(0.23);
        }
        else if (std::abs(rm.y - 2) < 0.1) {
            color += glm::vec3(1, 0, 0);
        }
        else if (std::abs(rm.y - 3) < 0.1) {
            // 绿色
            color += glm::vec3(0.3, 0.3, 0.3);
        }
        else if (std::abs(rm.y - 4) < 0.1) {
            color += glm::vec3(.6, 1, .2);
        }
    }
    return fromVec(color);
}