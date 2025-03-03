/**
 * 绘制基本图形
 */
#pragma once

#include "common.h"
#include "camera.h"
#include "stb_image.h"

using namespace std;

const int width = 800;
const int height = 600;
const float tMax = 200.0;
const float tMin = 0.1;
const int maxMarchTime = 128;
const float delta = 0.001;
float myTime = 0.0f;

const glm::vec4 circle = glm::vec4({-0.10, 0.2, -0.3, 0.1});
const glm::vec3 cir = {circle.x, circle.y, circle.z};

glm::vec2 map(glm::vec3 p);

/*---------------------------------------------------*/
glm::vec2 smin(float a, float b, float k) {
    k *= 6.0f;
    float h = std::max(k - std::abs(a - b), 0.0f) / k;
    float m = h * h * h * 0.5f;
    float s = m * k * (1.0f / 3.0f);
    return (a < b) ? glm::vec2(a - s, m) : glm::vec2(b - s, 1.0f - m);
}

float hash1(float n) {
    return std::fmod(std::sin(n) * 43758.5453123f, 1.0f);
}

glm::vec3 forwardSF(float i, float n) {
    const float PHI = 1.6180339887498948482045868343656f;
    float phi = 2.0f * PI * std::fmod(i / PHI, 1.0f);
    float zi = 1.0f - (2.0f * i + 1.0f) / n;
    float sinTheta = std::sqrt(1.0f - zi * zi);
    return glm::vec3(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, zi);
}

glm::vec2 mapSpecial(const glm::vec3& q) {
    glm::vec2 res = glm::vec2(q.y, 2.0f);

    float d = glm::length(q - glm::vec3(0.0f, 0.1f + 0.05f * std::sin(myTime), 0.0f)) - 0.1f;

    return smin(res.x, d, 0.05f);
}

glm::vec2 intersect(const glm::vec3& ro, const glm::vec3& rd) {
    const float maxd = 10.0f;

    glm::vec2 res(0.0f);
    float t = 0.0f;
    for (int i = 0; i < 512; ++i) {
        glm::vec2 h = map(ro + rd * t);
        if (h.x < 0.0f || t > maxd) break;
        t += h.x;
        res = glm::vec2(t, h.y);
    }

    if (t > maxd) res = glm::vec2(-1.0f);
    return res;
}

float calcAO(const glm::vec3& pos, const glm::vec3& nor, float ran) {
    float ao = 0.0f;
    const int num = 32;
    for (int i = 0; i < num; ++i) {
        glm::vec3 ap = forwardSF(static_cast<float>(i) + ran, static_cast<float>(num));
        ap *= glm::sign(glm::dot(ap, nor)) * hash1(static_cast<float>(i));
        ao += glm::clamp(map(pos + nor * 0.01f + ap * 0.2f).x * 20.0f, 0.0f, 1.0f);
    }
    ao /= static_cast<float>(num);
    return glm::clamp(ao, 0.0f, 1.0f);
}

float calcAO(const glm::vec3& pos, const glm::vec3& nor) {
    float occ = 0.0f;
    float sca = 1.0f;

    for (int i = 0; i < 5; ++i) {
        float h = 0.01f + 0.12f * static_cast<float>(i) / 4.0f; // 计算采样距离
        float d = map(pos + h * nor).x;                        // 距离场值
        occ += (h - d) * sca;                                 // 累积遮挡值
        sca *= 0.95f;                                         // 缩减权重
        if (occ > 0.35f) break;                               // 提前退出
    }

    // 返回环境光遮蔽值
    return glm::clamp(1.0f - 3.0f * occ, 0.0f, 1.0f) * (0.5f + 0.5f * nor.y);
}

glm::vec4 generateRandom(float seed) {
    // 使用简单的哈希函数生成伪随机数
    auto hash = [](float n) -> float {
        return std::fmod(std::sin(n) * 43758.5453123f, 1.0f);
    };
    return glm::vec4(hash(seed + 1.0f), hash(seed + 2.0f), hash(seed + 3.0f), hash(seed + 4.0f));
}

// 加载图像数据
unsigned char* loadImage(const char* filePath, int& width, int& height, int& channels) {
    return stbi_load(filePath, &width, &height, &channels, 4); // 强制加载为 RGBA 格式
}

// 模拟 texelFetch 函数
glm::vec4 texelFetch(unsigned char* imageData, int width, int height, const glm::ivec2& texCoord) {
    // 确保坐标在有效范围内
    int x = texCoord.x % width;
    int y = texCoord.y % height;

    // 计算像素索引
    int index = (y * width + x) * 4; // 每个像素有 4 个通道 (RGBA)

    // 提取纹素值
    glm::vec4 texel(
            static_cast<float>(imageData[index + 0]) / 255.0f,     // R
            static_cast<float>(imageData[index + 1]) / 255.0f,     // G
            static_cast<float>(imageData[index + 2]) / 255.0f,     // B
            static_cast<float>(imageData[index + 3]) / 255.0f      // A
    );

    return texel;
}

float calcSoftshadow(const glm::vec3& ro, const glm::vec3& rd, float mint, float tmax) {
    // 边界体积检查
    float tp = (0.8f - ro.y) / rd.y;
    if (tp > 0.0f) {
        tmax = std::min(tmax, tp); // 更新最大距离
    }

    float res = 1.0f; // 初始遮挡值
    float t = mint;   // 当前光线步进距离

    for (int i = 0; i < 24; ++i) {
        float h = map(ro + rd * t).x; // 距离场值
        float s = glm::clamp(8.0f * h / t, 0.0f, 1.0f); // 计算遮挡因子
        res = std::min(res, s);                        // 更新最小遮挡值
        t += glm::clamp(h, 0.01f, 0.2f);               // 步进距离
        if (res < 0.004f || t > tmax) break;           // 提前退出
    }

    // 返回柔化阴影值
    res = glm::clamp(res, 0.0f, 1.0f);
    return res * res * (3.0f - 2.0f * res);
}

/*---------------------------------------------------*/

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

glm::mat2 mat = glm::mat2(.6, -0.8, 0.8, 0.6);

glm::mat3 setCamera(glm::vec3 ro, glm::vec3 target, float cr) {
    glm::vec3 z = glm::normalize(target - ro);
    glm::vec3 up = glm::normalize(glm::vec3(glm::sin(cr), glm::cos(cr), 0.0f));
    glm::vec3 x = glm::cross(z, up);
    glm::vec3 y = glm::cross(x, z);
    return glm::mat3(x, y, z);
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

float groundL(const glm::vec2& x) {
    return sdfGround(glm::vec3 (x.x, x.y, 0));
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
    d = std::min(sdBox(p, {1,1,-1}, {0, -4, -7}), d);
    return d;
}

float shape2(glm::vec3 p)
{
//    auto d = sdfSphere(p, {0.01, 0.4, -0.4}, 0.2);
//    d = std::max(sdBox(p, {.1,0.2,-.2}, {0.02, 0.25, -.38}), d);
    auto d = sdBox(p, {.05,0.4,.03}, {-0.1, 0.2, -.38});
    return d;
}

float shape3(glm::vec3 p)
{
//    auto d1 = sdfSphere(p, {0.02, 0.3, 0.5}, 0.1);
//    auto d2 = sdBox(p, {0.01,0.1,0.1}, {.1, .3, .3});
//    auto d = std::max(d2 , -1 * d1);
    auto d = sdBox(p, {0.01,0.1,0.1}, {.1, .3, .3});
    return d;
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
    auto d = glm::vec2(sdfSphere(p, cir, circle.w), 22);
    d = vec2Min(d, {shape2(p), 10});
    d = vec2Min(d, {shape3(p), 2});
    d = vec2Min(d, mapSpecial(p));
    return d;
//    return mapSpecial(p);
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

//Color render(int x, int y, float iTime)
//{
//
//    auto uv = fixUV(x, y);
////    iTime = 1.5;
//
//    glm::vec3 col = glm::vec3(0);
//
//    float an = iTime * 0.04f;
//    float r = 10.0f / 4;
//    float base_offset = 5.5f + 5.f + 3.f;
//    glm::vec2 pos2d = glm::vec2(r * glm::sin(an), r * glm::cos(an));
////    float h = groundL(pos2d);
//    float h = 1;
//    glm::vec3 ro = glm::vec3(pos2d.x, h - 4, pos2d.y + base_offset);
//    glm::vec3 target = glm::vec3(r * glm::sin(an + 0.01f), h, r * glm::cos(an + 0.01f));
////    cout << ro.x << " " << ro.y << " " << ro.z << endl;
////    cout << target.x << " " << target.y << " " << target.z << endl;
////    cout << " ---- " << endl;
//    glm::mat3 cam = setCamera(ro, target, 0.0f);
//
//    float fl = 1.0f;
//    glm::vec3 rd = glm::normalize(cam * glm::vec3(uv, fl));
//
//
////    glm::vec2 uv = fixUV(x, y);
//
//    glm::vec3 bg = glm::vec3(0.7, 0.7, 0.9);
//    // 增加y轴上的渐变
//    glm::vec3 color = bg - glm::normalize(uv).y * glm::vec3(0.13);
////    camera cam;
////    auto ray = cam.ray(uv);
//
//    glm::vec3 light = glm::vec3(10, -15, -5);
//
////    auto rm = rayMarch(ray.ro, glm::normalize(ray.rd));
//    auto rm = rayMarch(ro, rd);
//
//    // rm.y = -1时, 表示没有hit object
//    if (rm.y > 0)
//    {
//        auto t = rm.x;
//        glm::vec3 p = ro + t * rd;
//        glm::vec3 n = (rm.y < 1.1) ? glm::vec3(0, -1, 0) : calcNormal(p);
//        float diff = glm::dot(
//            glm::normalize(light - p),
//            n);
//        // 绘制阴影
//        diff *= softshadow(p, glm::normalize(light - p), 0.1, 100, 0.15);
//        // 添加环境光
//        float amb = 0.1;
//        color = glm::vec3(1) * diff + amb;
//        if (std::abs(rm.y - 1) < 0.1) {
//            color += glm::vec3(0.23);
//        }
//        else if (std::abs(rm.y - 2) < 0.1) {
//            color += glm::vec3(1, 0, 0);
//        }
//        else if (std::abs(rm.y - 3) < 0.1) {
//            // 绿色
//            color += glm::vec3(0.3, 0.3, 0.3);
//        }
//        else if (std::abs(rm.y - 4) < 0.1) {
//            color += glm::vec3(.6, 1, .2);
//        }
//    }
//    return fromVec(color);
//}

int w, h, channels;
unsigned char* imageData = loadImage("E:\\c++_code\\rmRenderer\\channel\\1.jpg", w, h, channels);

Color render(int x, int y, float iTime) {
    glm::vec2 uv = fixUV(x, y);
    myTime = iTime;

    // 加载图像
    glm::ivec2 texCoord(width, height); // 采样坐标
    glm::vec4 ran = texelFetch(imageData, width, height, texCoord);

    // glm::vec4 ran = generateRandom(static_cast<float>(x + y * 10000  + static_cast<int>(iTime * 10000)));

    float an = 1.f * iTime;
    glm::vec3 ro = glm::vec3(.4f * std::sin(an), 0.15f, .4f * std::cos(an));
    glm::vec3 ta = glm::vec3(0.0f, 0.05f, 0.0f);
    glm::vec3 ww = glm::normalize(ta - ro);
    glm::vec3 uu = glm::normalize(glm::cross(ww, glm::vec3(0.0f, 1.0f, 0.0f)));
    glm::vec3 vv = glm::normalize(glm::cross(uu, ww));
    glm::vec3 rd = glm::normalize(uv.x * uu + uv.y * vv + 1.7f * ww);

    // 拉远距离
    ro = ro * glm::vec3(2.);

//    cout << "ro: " << ro.x << " " << ro.y << " " << ro.z << endl;
//    cout << "rd: " << rd.x << " " << rd.y << " " << rd.z << endl;

    // Step 4: Ray marching
    glm::vec3 col(1.0f);
    glm::vec2 res = intersect(ro, rd);
    float t = res.x;
    if (t > 0.0f) {
        if (res.y < 1) {
            glm::vec3 pos = ro + t * rd;
            // debug
//            cout << "pos: " << pos.x << " " << pos.y << " " << pos.z << endl;

            glm::vec3 nor = calcNormal(pos);
            glm::vec3 ref = glm::reflect(rd, nor);
            float fre = glm::clamp(1.0f + glm::dot(nor, rd), 0.0f, 1.0f);
            float occ = calcAO(pos, nor, ran.y); occ = occ * occ;

            col = glm::mix(glm::vec3(0.0f, 0.05f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f), res.y);

            col = col * 0.72f + 0.2f * fre * glm::vec3(1.0f, 0.8f, 0.2f);
            glm::vec3 lin = 4.0f * glm::vec3(0.7f, 0.8f, 1.0f) * (0.5f + 0.5f * nor.y) * occ;
            lin += 0.8f * glm::vec3(1.0f, 1.0f, 1.0f) * fre * (0.6f + 0.4f * occ);
            col = col * lin;
            col += 2.0f * glm::vec3(0.8f, 0.9f, 1.0f) * glm::smoothstep(0.0f, 0.4f, ref.y) *
                   (0.06f + 0.94f * std::pow(fre, 5.0f)) * occ;
            col = glm::mix(col, glm::vec3(1.0f), 1.0f - std::exp2(-0.04f * t * t));

            col = glm::pow(col, glm::vec3(0.4545f));
            col *= 0.9f;
            col = glm::clamp(col, 0.0f, 1.0f);
            col = col * col * (3.0f - 2.0f * col);

            col += (ran.x - 0.5f) / 255.0f;
        } else {
            auto t = res.x;
            // glm::vec3 col = glm::vec3(0.7, 0.7, 0.9) - glm::max(rd.y,0.0f)*glm::vec3(0.3);
            col = glm::vec3(0.2) + glm::vec3(0.2) * sin(glm::vec3(res.y * 2.0) + glm::vec3(0.0, 1.0, 2.0));
            float ks = 1.0;

            glm::vec3 pos = ro + t*rd;
            glm::vec3 nor = calcNormal( pos );
            glm::vec3 ref = glm::reflect( rd, nor );

            // 环境光遮蔽
            float occ = calcAO(pos, nor);

            // 初始化光照贡献
            glm::vec3 lin(0.0f);

            // 太阳光
            {
                glm::vec3 lig = glm::normalize(glm::vec3(-4.5f, 1.4f, 5.6f)); // 光源方向
                glm::vec3 hal = glm::normalize(lig - rd);                     // 半角向量
                float dif = glm::clamp(glm::dot(nor, lig), 0.0f, 1.0f);       // 漫反射
                dif *= calcSoftshadow(pos, lig, 0.02f, 2.5f);                 // 柔化阴影
                float spe = std::pow(glm::clamp(glm::dot(nor, hal), 0.0f, 1.0f), 16.0f); // 高光
                spe *= dif; // 高光受漫反射影响
                spe *= 0.04f + 0.96f * std::pow(glm::clamp(1.0f - glm::dot(hal, lig), 0.0f, 1.0f), 5.0f); // 高光衰减
                lin += col * 2.20f * dif * glm::vec3(1.30f, 1.00f, 0.70f); // 漫反射贡献
                lin += 5.00f * spe * glm::vec3(1.30f, 1.00f, 0.70f) * ks;  // 高光贡献
            }

            // 天空光
            {
                float dif = std::sqrt(glm::clamp(0.5f + 0.5f * nor.y, 0.0f, 1.0f)); // 天空漫反射
                dif *= occ;                                                       // 环境光遮蔽
                float spe = glm::smoothstep(-0.2f, 0.2f, ref.y);                  // 反射高光
                spe *= dif;
                spe *= 0.04f + 0.96f * std::pow(glm::clamp(1.0f + glm::dot(nor, rd), 0.0f, 1.0f), 5.0f); // 高光衰减
                spe *= calcSoftshadow(pos, ref, 0.02f, 2.5f);                                           // 柔化反射阴影
                lin += col * 0.60f * dif * glm::vec3(0.40f, 0.60f, 1.15f); // 天空漫反射贡献
                lin += 2.00f * spe * glm::vec3(0.40f, 0.60f, 1.30f) * ks;  // 天空高光贡献
            }

            // 背光
            {
                float dif = glm::clamp(glm::dot(nor, glm::normalize(glm::vec3(0.5f, 0.0f, 0.6f))), 0.0f, 1.0f) *
                            glm::clamp(1.0f - pos.y, 0.0f, 1.0f); // 背光漫反射
                dif *= occ;                                      // 环境光遮蔽
                lin += col * 0.55f * dif * glm::vec3(0.25f, 0.25f, 0.25f); // 背光贡献
            }

            // 次表面散射 (SSS)
            {
                float dif = std::pow(glm::clamp(1.0f + glm::dot(nor, rd), 0.0f, 1.0f), 2.0f); // SSS 强度
                dif *= occ;                                                                 // 环境光遮蔽
                lin += col * 0.25f * dif * glm::vec3(1.00f, 1.00f, 1.00f);                   // SSS 贡献
            }

            col = lin;

            col = glm::mix( col, glm::vec3(0.7,0.7,0.9), 1.0-exp( -0.0001*t*t*t ) );
        }
    }


    return fromVec(col);
}