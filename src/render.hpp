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

glm::vec2 sdfSphere(glm::vec3 p, glm::vec3 o, float r)
{
    // {距离, tag}
    return {glm::length(p - o) - r, 1};
}

glm::vec2 sdfGround(glm::vec3 p)
{
    // 很操蛋, 我的这个y轴是向下为正
    // 所以这里应该是地面坐标 - p.y
    return {- p.y, 2};
}

glm::vec2 vec2Min(glm::vec2 a, glm::vec2 b) {
    if (a.x <= b.x) {
        return a;
    } else {
        return b;
    }
}

glm::vec2 map(glm::vec3 p)
{
    auto s = sdfSphere(p, {circle.x, circle.y, circle.z}, circle.w);
    return vec2Min(s, sdfGround(p));
}

Color fromVec(glm::vec3 v)
{
    return Color{
        convert(v.x * 255),
        convert(v.y * 255),
        convert(v.z * 255),
        255};
}

float softshadow( glm::vec3 ro, glm::vec3 rd, float mint, float maxt, float k )
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

glm::vec2 rayMarch(glm::vec3 ro, glm::vec3 rd)
{
    float t = 0;

    for (int i = 0; i < maxMarchTime && t < tMax; i++)
    {
        glm::vec3 p = ro + rd * t;
        auto m = map(p);
        auto d = m.x;
        // hit object
        if (d < delta)
            return {t, m.y};
        t += d;
    }

    return {t, -1};
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

    glm::vec3 color(0);
    camera cam;
    auto ray = cam.ray(uv);
    glm::vec3 light = glm::vec3(10, -15, 7);

    auto rm = rayMarch(ray.ro, ray.rd);

    // rm.y = -1时, 表示没有hit object
    if (rm.y > 0)
    {
        auto t = rm.x;
        glm::vec3 p = ray.at(t);
        glm::vec3 n = calcNormal(p);
        float diff = glm::dot(
            glm::normalize(light - p),
            n);
        // 绘制阴影
        // diff *= normalShadow(p, light);
        diff *= softshadow(p, glm::normalize(light - p), 0.1, 100, 5);
        // 增加环境光
        // float amb = 0.5 + 0.5 * dot(n, glm::vec3(0, 1, 0));
        // color = glm::vec3(1) * diff + amb * glm::vec3(0.5);
        float amb = 0.1;
        color = glm::vec3(1) * diff + amb;
        if (rm.y == 1) {
            color += glm::vec3(1, 0, 0);
        } else if (rm.y == 2) {
            color += glm::vec3(0, 0, 1);
        }
    }
    return fromVec(color);
}