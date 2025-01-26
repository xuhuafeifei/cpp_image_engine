#pragma once

#include "raylib.h"
#include "glm/glm.hpp"

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;

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

float sdf(glm::vec3 p)
{
    return sdfSphere(p, {0, 2, 7}, 2);
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
        float d = sdf(p);
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
            sdf(p + v1 * delta) * v1 +
            sdf(p + v2 * delta) * v2 +
            sdf(p + v3 * delta) * v3 +
            sdf(p + v4 * delta) * v4));
}

Color render(int x, int y)
{
    glm::vec2 uv = fixUV(x, y);

    glm::vec3 color(0);
    glm::vec3 ro(0, 0, -1.5);
    glm::vec3 rd = glm::normalize(glm::vec3(uv, 0) - ro);
    glm::vec3 light = glm::vec3(-10, 15, 7);

    float t = rayMarch(ro, rd);

    if (t < tMax)
    {
        glm::vec3 p = ro + rd * t;
        glm::vec3 n = calcNormal(p);
        float diff = glm::dot(
            glm::normalize(p - light),
            n);
        color = glm::vec3(1) * diff;
    }
    return fromVec(color);
}