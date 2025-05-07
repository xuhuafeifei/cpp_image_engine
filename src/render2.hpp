#pragma once

#include "common.h"
#include "camera.h"

#include <random>

// 设置每像素采样次数
const int SAMPLES_PER_PIXEL = 4;

// 简单的随机浮点数生成器
float randf() {
    static std::mt19937 gen(42); // 固定种子便于调试
    static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(gen);
}

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;

const glm::vec4 circle = glm::vec4({0, 0, 7, 2});

glm::vec3 getCamera()
{

    return glm::vec3(0., 0, -1.5);
}

glm::vec2 fixUV(int x, int y, float offsetX = 0.0f, float offsetY = 0.0f)
{
    float m = std::min(width, height);
    float u = (2.0 * (x + offsetX) - width) / m;
    float v = (2.0 * (y + offsetY) - height) / m;
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

float map(glm::vec3 p)
{
    return std::min(sdfSphere(p, {circle.x, circle.y, circle.z}, circle.w), sdfGround(p));
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
        float h = map(ro + rd*t);
        if( h<0.001 )
            return 0.0;
        res = std::min( res, k*h/t );
        t += h;
    }
    return res;
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

float normalShadow(glm::vec3 p, glm::vec3 light)
{
    float t = rayMarch(p, glm::normalize(light - p));
    if (t < tMax)
    {
        return 0.1;
    }
    return 1;
}

Color render(int x, int y)
{
    glm::vec3 colorAccum(0.0f);

    for (int s = 0; s < SAMPLES_PER_PIXEL; ++s)
    {
        // 随机偏移 [0, 1)
        float offsetX = randf();
        float offsetY = randf();

        glm::vec2 uv = fixUV(x, y, offsetX, offsetY);

        camera cam;
        auto ray = cam.ray(uv);
        glm::vec3 light = glm::vec3(10, -15, 7);

        float t = rayMarch(ray.ro, ray.rd);

        if (t < tMax)
        {
            glm::vec3 p = ray.at(t);
            glm::vec3 n = calcNormal(p);
            float diff = glm::dot(
                    glm::normalize(light - p),
                    n);

            // 使用软阴影（可选硬阴影替换）
            diff *= normalShadow(p, light);
            // diff *= softshadow(p, glm::normalize(light - p), 0.1, 100, 5);

            float amb = 0.1;
            colorAccum += glm::vec3(1) * diff + amb;
        }
        else
        {
            // 背景色
            colorAccum += glm::vec3(0.0f);
        }
    }

    // 取平均
    glm::vec3 finalColor = colorAccum / glm::vec3(SAMPLES_PER_PIXEL);

    return fromVec(finalColor);
}