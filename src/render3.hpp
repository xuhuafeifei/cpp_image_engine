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

#define EPSILON 0.001f
#define MAX_DIST 200.0f
#define MAX_ITER 200

#define GROUND_TIME 8

const int width = 800;
const int height = 600;
const float tMax = 200.0 / 2;
const float tMin = 0.1;
const int maxMarchTime = (128 / 2);
const float delta = 0.001;

const glm::vec4 circle = glm::vec4({0, 0, 7, 2});
const float y_offset = 0.2;

//using namespace glm;

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

//glm::vec3 noise(glm::vec2 pos)
//{
//    auto i = glm::floor(pos);
//    auto f = glm::fract(pos);
//    auto u = glm::smoothstep(0.f, 1.f, f);
//    auto du = glm::vec2(6.) * u * (glm::vec2(1.) - u);
//
//    auto a = random(i);
//    auto b = random(i + glm::vec2(1, 0));
//    auto c = random(i + glm::vec2(0, 1));
//    auto d = random(i + glm::vec2(1, 1));
//
//    return glm::vec3(
//            a + (b - a) * u.x * (1. - u.y) + (c - a) * (1. - u.x) * u.y + (d - a) * u.x * u.y,
//            du * (glm::vec2(b - a, c - a) + (a - b - c + d) * glm::vec2(u.y, u.x))
//    );
//}
glm::vec3 noise(glm::vec2 pos) {
    glm::vec2 i = glm::floor(pos);
    glm::vec2 f = glm::fract(pos);
    glm::vec2 u = f * f * (3.0f - 2.0f * f);
    glm::vec2 du = 6.0f * u * (1.0f - u);

    float a = random(i);
    float b = random(i + glm::vec2(1.0f, 0.0f));
    float c = random(i + glm::vec2(0.0f, 1.0f));
    float d = random(i + glm::vec2(1.0f, 1.0f));

    return glm::vec3(a + (b - a) * u.x * (1.0f - u.y) +
                     (c - a) * (1.0f - u.x) * u.y +
                     (d - a) * u.x * u.y, du * (glm::vec2(b - a, c - a) +
                                                (a - b - c + d) * glm::vec2(u.y, u.x)));
}

glm::mat2 mat = glm::mat2(.6, -0.8, 0.8, 0.6);

glm::mat3 setCamera(glm::vec3 ro, glm::vec3 target, float cr) {
    glm::vec3 z = glm::normalize(target - ro);
    glm::vec3 up = glm::normalize(glm::vec3(glm::sin(cr), glm::cos(cr), 0.0f));
    glm::vec3 x = glm::cross(z, up);
    glm::vec3 y = glm::cross(x, z);
    return glm::mat3(x, y, z);
}

float ground(glm::vec2 x) {
    glm::vec2 p = 0.003f * x;
//    auto p = x;
    float a = 0.0f;
    float b = 1.0f;
    glm::vec2 d = glm::vec2(0);

    for(int i = 0; i < 8; i++) {
        glm::vec3 n = noise(p);
        d += glm::vec2(n.y, n.z);
        a += b * n.x / (1.0f + glm::dot(d, d));
        p = mat * p * 2.0f;
        b *= 0.5f;
    }

    return 120.0f * a;
//    return a;
//    return noise(x).x;
}


float ground(glm::vec3 x)
{
//    auto p = glm::vec3(0.005) * x;
//    float a = 0.;
//    float b = 1.;
//    glm::vec2 d = glm::vec2(0);
//
//    for (int i = 0; i < 8; ++i)
//    {
//        auto n = noise(p);
//        d += glm::vec2(n.y, n.z);
//        a += b * n.x / (1. + glm::dot(d, d));
//        glm::vec2 p2(p.x, p.y);
//        p2 = mat * p2 * glm::vec2(2.);
//        p.x = p2.x;
//        p.y = p2.y;
//        b *= 0.5;
//    }
//
//    return 120 * a;
    return ground({x.x, x.y});
}

//float ground(glm::vec2 x)
//{
////    auto p = glm::vec2(0.008) * x;
//    auto p = x;
//    float a = 0.;
//    float b = 1.;
//    glm::vec2 d = glm::vec2(0);
//
//    for (int i = 0; i < 8; ++i)
//    {
//        auto n = noise(x);
//        d += glm::vec2(n.y, n.z);
//        a += b * n.x / (1. + glm::dot(d, d));
//        glm::vec2 p2(p.x, p.y);
//        p2 = mat * p2 * glm::vec2(2.);
//        p.x = p2.x;
//        p.y = p2.y;
//        b *= 0.5;
//    }
//
////    return 80 * a;
//    return a;
//}


//float ground(glm::vec2 x)
//{
////    glm::vec2 p = glm::vec2(0.005) * x;
//    auto p = x;
//    auto a = 0.;
//    auto b = 1.;
//    glm::vec2 d = glm::vec2(0);
//
//    for (int i = 0; i < GROUND_TIME; ++i)
//    {
//        glm::vec3 n = noise(p);
//        d += glm::vec2(n.y, n.z);
//        a += b * n.x / (1. + glm::dot(d, d));
//        p = mat * p * glm::vec2(2.);
//        b *= 0.5;
//    }
//
//    return a;
////    return 80 * a;
//}

//float groundH(glm::vec2 x)
//{
//    glm::vec2 p = glm::vec2(0.005) * x;
////    auto p = x;
//    auto a = 0.;
//    auto b = 1.;
//    glm::vec2 d = glm::vec2(0);
//
//    for (int i = 0; i < GROUND_TIME; ++i)
//    {
//        glm::vec3 n = noise(p);
//        d += glm::vec2(n.y, n.z);
//        a += b * n.x / (1. + glm::dot(d, d));
//        p = mat * p * glm::vec2(2.);
//        b *= 0.5;
//    }
//
//    return 80 * a;
////    return a;
//}

float groundH(glm::vec2 x) {
    glm::vec2 p = 0.003f * x;
    float a = 0.0f;
    float b = 1.0f;
    glm::vec2 d = glm::vec2(0);

    for(int i = 0; i < 12; i++) {
        glm::vec3 n = noise(p);
        d += glm::vec2(n.y, n.z);
        a += b * n.x / (1.0f + glm::dot(d, d));
        p = mat * p * 2.0f;
        b *= 0.5f;
    }

    return 120.0f * a;
}

float groundL(const glm::vec2& x) {
    glm::vec2 p = 0.003f * x;
    float a = 0.0f;
    float b = 1.0f;
    glm::vec2 d(0.0f, 0.0f);  // 代替 glm::vec2(0)

    for (int i = 0; i < 3; i++) {
        glm::vec3 n = noise(p);  // 获取噪声值
        d += glm::vec2(n.y, n.z);  // 将噪声的 y 和 z 分量作为二维向量
        a += b * n.x / (1.0f + glm::dot(d, d));  // 计算 a

        // 使用 glm::mat2 对 p 进行变换
        p = mat * p * 2.0f;  // 对 p 进行二维矩阵变换

        b *= 0.5f;  // 更新 b
    }

    return 120.0f * a;  // 返回计算结果
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
//    float t = tMin;
//
//    for (int i = 0; i < maxMarchTime && t < tMax; i++)
//    {
//        glm::vec3 p = ro + rd * t;
//        float d = map(p);
//        if (d < delta * t)
//            break;
//        t += .4f * d;
//    }
//
//    return t;
    return 11111.;
}

float rayMarch(glm::vec3 ro, glm::vec3 rd, float tmin, float tmax) {
    float t = tmin;
    for(int i = 0; i < MAX_ITER && t < tmax; i++) {
        glm::vec3 p = ro + t * rd;
        float h = p.y - ground(glm::vec2(p.x, p.z));
        if(abs(h) < EPSILON * t)
            break;
        t += 0.3f * h;
    }
    return t;
}

//glm::vec3 calcNormal(glm::vec3 p)
//{
//    const glm::vec3 v1(1, -1, -1);
//    const glm::vec3 v2(-1, -1, 1);
//    const glm::vec3 v3(-1, 1, -1);
//    const glm::vec3 v4(1);
//    return glm::normalize(
//            glm::vec3(
//                    map(p + v1 * delta) * v1 +
//                    map(p + v2 * delta) * v2 +
//                    map(p + v3 * delta) * v3 +
//                    map(p + v4 * delta) * v4));
//}

glm::vec3 calcNorm(glm::vec3 p, float t) {
    glm::vec2 epsilon = glm::vec2(0.0027f * t, 0.0f);
    return glm::normalize(
            glm::vec3(
                        groundH(glm::vec2(p.x, p.z) - glm::vec2(epsilon.x, epsilon.y)) - groundH(glm::vec2(p.x, p.z) + glm::vec2(epsilon.x, epsilon.y)),
                        2.0f * epsilon.x,
                        groundH(glm::vec2(p.x, p.z) - glm::vec2(epsilon.y, epsilon.x)) - groundH(glm::vec2(p.x, p.z) + glm::vec2(epsilon.y, epsilon.x))
                    )
            );
}

//float softshadow( glm::vec3 ro, glm::vec3 rd, float mint, float maxt, float w )
//{
//    float res = 1.0;
//    float ph = 1e20;
//    float t = mint;
//    for( int i=0; i<(256) && t<maxt; i++ )
//    {
//        auto k = ro + rd*t;
//        float h = map(k);
//        if( h<0.001 )
//            return 0.0;
//        float y = h*h/(2.0*ph);
//        float d = std::sqrt(h*h-y*y);
//        res = std::min( res, d/(w*std::max(0.0f,t-y)) );
//        ph = h;
//        t += h;
//    }
//    return res;
//}

float softShadow(glm::vec3 ro, glm::vec3 rd, float dis) {
    float minStep = glm::clamp(0.01f * dis, 0.5f, 50.0f);
    float res = 1.0f;
    float t = 0.001f;
    for(int i = 0; i < 80; i++) {
        glm::vec3 p = ro + t * rd;
        float h = p.y - ground(glm::vec2(p.x, p.z));
        res = glm::min(res, 8.0f * h / t);
        t += glm::max(minStep, h);
        if(res < 0.001f || p.y > 200.0f)
            break;
    }
    return glm::clamp(res, 0.0f, 1.0f);
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

// new insert
float softshadow( glm::vec3 ro, glm::vec3 rd, float mint, float maxt)
{
    float t = mint;
    for( int i=0; i<256 && t<maxt; i++ )
    {
        auto v = ro + rd*t;
        float h = ground({v.x, v.y});
        if( h<0.001 )
            return 0.0;
        t += h;
    }
    return 1.0;
}

// new insert
glm::vec3 calcNormal(glm::vec3 p)
{
    const glm::vec3 v1(1, -1, -1);
    const glm::vec3 v2(-1, -1, 1);
    const glm::vec3 v3(-1, 1, -1);
    const glm::vec3 v4(1);
    return glm::normalize(
            glm::vec3(
                    ground(p + v1 * delta) * v1 +
                    ground(p + v2 * delta) * v2 +
                    ground(p + v3 * delta) * v3 +
                    ground(p + v4 * delta) * v4));
}

Color render(int x, int y, float iTime)
{
    auto uv = fixUV(x, y);
//    iTime = 1.5;

    glm::vec3 col = glm::vec3(0);

    float an = iTime * 0.04f;
    float r = 100.0f / 4;
    glm::vec2 pos2d = glm::vec2(r * glm::sin(an), r * glm::cos(an));
    float h = groundL(pos2d) + 10.0f;
    glm::vec3 ro = glm::vec3(pos2d.x, h, pos2d.y);
    glm::vec3 target = glm::vec3(r * glm::sin(an + 0.01f), h, r * glm::cos(an + 0.01f));
    glm::mat3 cam = setCamera(ro, target, 0.0f);

    float fl = 1.0f;
    glm::vec3 rd = glm::normalize(cam * glm::vec3(uv, fl));

    float tmin = 0.001f;
    float tmax = 1000.0f;

    float maxh = 300.0f;

    float tp = (maxh - ro.y) / rd.y;
    if(tp > 0.0f) {
        if(maxh > ro.y)
            tmax = glm::min(tmax, tp);
        else
            tmin = glm::max(tmin, tp);
    }
    float t = rayMarch(ro, rd, tmin, tmax);
    glm::vec3 sunlight = glm::normalize(glm::vec3(0.4f, 0.4f, -0.2f));
    float sundot = glm::clamp(glm::dot(rd, sunlight), 0.0f, 1.0f);


    if (t < tMax)
    {
        // new insert
        {
            glm::vec3 p = ro + t * rd;
            glm::vec3 n = calcNorm(p, t);
            float diff = glm::dot(
                    glm::normalize(sunlight - p),
                    n);
            // 绘制阴影
//         diff *= normalShadow(p, light);
            diff *= softshadow(p, glm::normalize(sunlight - p), 0.1, 100);
            float amb = 0.1;
            col = glm::vec3(1) * diff + amb;
        }


//        glm::vec3 p = ro + t * rd;
//        glm::vec3 n = calcNorm(p, t);
//        glm::vec3 difColor = glm::mix(glm::vec3(0.08f, 0.05f, 0.03f), glm::vec3(0.10f, 0.09f, 0.08f), noise(glm::vec2(p.x, p.z) * 0.02f).x);
//        float r = noise(glm::vec2(p.x, p.z) * 0.1f).x;
//
//        // rocks
//        col = (r * 0.55f + 0.75f) * 2.9f * difColor;
//        col = glm::mix(col, glm::vec3(0.09f, 0.06f, 0.03f) * (0.5f + 0.5f * r), n.y);
//        col = glm::mix(col, glm::vec3(0.085f, 0.045f, 0.015f) * (0.25f + 0.75f * r), glm::smoothstep(0.95f, 1.0f, n.y));
//        col *= 0.1f + 1.8f * glm::sqrt(fbm(glm::vec2(p.x, p.z) * 0.04f) * fbm(glm::vec2(p.x, p.z) * 0.01f));
//
//        // Snow
//        float h = glm::smoothstep(10.0f, 15.0f, p.y + 12.0f * fbm(0.01f * glm::vec2(p.x, p.z)));
////        float e = glm::smoothstep(1.0f - 0.5f * h, 1.0f - 0.1f * h, n.y);
////        float o = 0.15f + 0.3f * glm::smoothstep(0.0f, 0.02f, n.y + h * h);
//        float e = glm::smoothstep(1.0f - 0.5f * h, 1.0f - 0.1f * h, n.y);
//        float o = 0.25f + 0.008f * glm::smoothstep(0.0f, 0.1f, n.y + h * h);
//
//        float s = h * e * o;
//        col = glm::mix(col, 0.89f * glm::vec3(0.62f, 0.65f, 0.7f), s);
//
//        // Linear Lighting
//        glm::vec3 lin = glm::vec3(0.0f);
//
//        float dif = glm::clamp(glm::dot(sunlight, n), 0.0f, 1.0f);
//        float sh = softShadow(p + 0.01f * sunlight, sunlight, t);
//        float amb = glm::clamp(0.5f + 0.5f * n.y, 0.8f, 1.0f);
//        float bac = glm::clamp(0.2f + 0.8f * glm::dot(glm::vec3(-sunlight.x, 0.0f, sunlight.z), n), 0.0f, 1.0f);
//        lin += dif * glm::vec3(8.0f, 5.0f, 3.0f) * 1.8f * glm::vec3(sh, sh * sh * 0.5f + 0.5f * sh, sh * sh * 0.8f + 0.2f * sh);
//        lin += amb * glm::vec3(0.4f, 0.6f, 1.0f) * 2.f;
//        lin += bac * glm::vec3(0.4f, 0.5f, 0.6f);
//
//        col *= lin;
//
//        // half-angle
//        glm::vec3 hal = glm::normalize(sunlight - rd);
//
//        col += (0.7f + 0.3f * s) * (0.04f + 0.96f * glm::pow(glm::clamp(1.0f + glm::dot(hal, rd), 0.0f, 1.0f), 5.0f)) *
//               glm::vec3(7.0f, 5.0f, 3.0f) * sh * dif *
//               glm::pow(glm::clamp(glm::dot(n, hal), 0.0f, 1.0f), 16.0f);
//
//        col = glm::mix(col, 0.65f * glm::vec3(0.5f, 0.75f, 1.0f), 1.0f - glm::exp(-glm::pow(0.002f * t, 1.5f)));
    }
    else
    {
        // sky
//        col = glm::vec3(0.3f, 0.5f, 0.85f) - (-rd.y) * (-rd.y) * 0.5f;
//        col = glm::mix(col, 0.85f * glm::vec3(0.7f, 0.75f, 0.85f), glm::pow(1.0f - glm::max(rd.y, 0.0f), 4.0f));
////
//        // sun
//        col += 0.25f * glm::vec3(1.0f, 0.7f, 0.4f) * glm::pow(sundot, 5.0f);
//        col += 0.25f * glm::vec3(1.0f, 0.8f, 0.6f) * glm::pow(sundot, 64.0f);
//        col += 0.2f * glm::vec3(1.0f, 0.8f, 0.6f) * glm::pow(sundot, 512.0f);
//
//        // clouds
//        glm::vec2 skyPos = glm::vec2(ro.x, ro.z) + glm::vec2(rd.x, rd.z) * (200.0f - (-ro.y)) / (-rd.y) + iTime * 5.0f;
//        col = glm::mix(col, glm::vec3(1.0f, 0.95f, 1.0f), 0.5 * glm::smoothstep(0.1f, 0.5f, fbm(.005f * skyPos)));

    }

    // sun scatter
//    col += 0.15f * glm::vec3(1.0f, 0.7f, 0.3f) * glm::pow(sundot, 8.0f);
    col += 0.15;

    return fromVec(col);
}

#endif //RMRENDERER_RENDER3_H
