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

#define GROUND_TIME 14

const int width = 800;
const int height = 600;
const float tMax = 200.0 / 2;
const float tMin = 0.1;
const int maxMarchTime = (128 / 2);
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
    auto u = glm::smoothstep(0.f, 1.f, f);
    auto du = glm::vec2(6.) * u * (glm::vec2(1.) - u);

    auto a = random(i);
    auto b = random(i + glm::vec2(1, 0));
    auto c = random(i + glm::vec2(0, 1));
    auto d = random(i + glm::vec2(1, 1));

    return glm::vec3(
            a + (b - a) * u.x * (1. - u.y) + (c - a) * (1. - u.x) * u.y + (d - a) * u.x * u.y,
            du * (glm::vec2(b - a, c - a) + (a - b - c + d) * glm::vec2(u.y, u.x))
    );
}

glm::mat2 mat = glm::mat2(.6, -0.8, 0.8, 0.6);

float ground(glm::vec3 x)
{
    auto p = glm::vec3(0.005) * x;
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
    }

    return 120 * a;
}

float ground(glm::vec2 x)
{
//    auto p = glm::vec2(0.008) * x;
    auto p = x;
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

//    return 80 * a;
    return a;
}

float groundH(glm::vec2 x)
{
    glm::vec2 p = glm::vec2(0.005) * x;
//    auto p = x;
    auto a = 0.;
    auto b = 1.;
    glm::vec2 d = glm::vec2(0);

    for (int i = 0; i < GROUND_TIME; ++i)
    {
        glm::vec3 n = noise(p);
        d += glm::vec2(n.y, n.z);
        a += b * n.x / (1. + glm::dot(d, d));
        p = mat * p * glm::vec2(2.);
        b *= 0.5;
    }

    return 80 * a;
//    return a;
}

float groundH(glm::vec3 x)
{
    return groundH({x.x, x.y});
}

float map(glm::vec3 p)
{
    return groundH({p.x, p.z, 0}) - 1 - p.y;
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

    for (int i = 0; i < maxMarchTime && t < tMax; i++)
    {
        glm::vec3 p = ro + rd * t;
        float d = map(p);
        if (d < delta * t)
            break;
        t += .4f * d;
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
                    mapGround({p.x - e.x, p.z - e.y}) - mapGround({p.x + e.x, p.z + e.y}),
                    -2.0 * e.x,
                    mapGround({p.x - e.y, p.z - e.x}) - mapGround({p.x + e.y, p.z + e.x})
                    )
            );
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
    glm::vec3 color = glm::vec3(163./255, 195./255, 247./255);
//    glm::vec3 color = glm::vec3(0.3, 0.5, 0.85) - glm::vec3(0.2);
    // 添加太阳(x, y)
    glm::vec2 sun = glm::vec2(-1, -1);
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

// 把x从0-正无穷放缩到0-1, 且单调递增
float sigmoid(float x)
{
//    return 1. - glm::exp(-pow(0.002 * x, 1.5));
    float alpha = 0.01; // 控制与环境融合的程度, alpha越小, 整体融合程度越小
    return 1 / (1 + glm::exp(-alpha * x));
}

glm::vec3 fog_draw(glm::vec3 color, float t)
{
    if (t <= (tMax / 5.0 * 2.5))
    {
        return color;
    }
    // 越远处, mix越大越接近灰色, 越近处, mix越小越接近本来的颜色
    return glm::mix(color, glm::vec3(0.5) * glm::vec3(0.5, 0.75, 1.), sigmoid(t));
}

float softshadow( glm::vec3 ro, glm::vec3 rd, float mint, float maxt, float w )
{
    float res = 1.0;
    float ph = 1e20;
    float t = mint;
    for( int i=0; i<(256) && t<maxt; i++ )
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

float fbm(glm::vec2 p)
{
    auto a = 0.f;
    auto fac = 0.45f;
    for (int i = 0; i < 4; ++i)
    {
        a += fac * noise(p).x;
        p = glm::vec2(1.5) * mat * p;
        fac *= 0.15;
    }
    return a;
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
    glm::vec3 light = glm::vec3(-10, -15, 20);

    auto ray = camMat.getRay(uv);

    auto ro = ray.ro;
    auto rd = ray.rd;

    float t = rayMarch(ro, rd);

    if (t < tMax)
    {
        glm::vec3 p = ro + rd * t;
        glm::vec3 n = calcNormalGround(p);
        color = glm::vec3(0.57, 0.47, 0.34);

        float diff = glm::dot(
                glm::normalize(light - p),
                n);
        diff = glm::clamp(diff, 0.f, 1.f);
        auto lin = glm::vec3(0.);
        lin += diff;
        // 软阴影
        auto sh = softshadow(p, glm::normalize(light - p), 0.01, 100, 2.); // Soft shadows

//        lin += glm::vec3(diff * 1.3) * glm::vec3(sh, sh * sh * 0.5 + 0.5 * sh, sh * sh * 0.8 + 0.2 * sh);
        lin += glm::vec3(diff * 0.75) * glm::vec3(sh);

        color *= lin;
        // 绘制雾气
        color = fog_draw(color, t);
    }
    else
    {
        // 绘制天空
        color = sky_draw_my(uv);
        // clouds
        float cloudH = 100.;
        // y轴朝下为正
        glm::vec3 cloudUV = ro + glm::vec3(std::abs((cloudH + ro.y) / rd.y)) * rd;
//        std::cout << cloudUV.x << " " << cloudUV.y << " " << cloudUV.z << std::endl;
//        color = glm::mix(color, glm::vec3(1., 0.95, 1.), glm::smoothstep(0.4f, 0.6f, fbm(glm::vec2(cloudUV.x, cloudUV.z) * glm::vec2(0.005))));
        color = glm::mix(color, glm::vec3(1., 0.95, 1.), glm::smoothstep(0.2f, 0.8f, fbm(glm::vec2(cloudUV.x, cloudUV.z) * glm::vec2(0.01))));
    }
    return fromVec(color);
}

#endif //RMRENDERER_RENDER3_H
