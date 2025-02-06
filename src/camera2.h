//
// Created by 25080 on 2025/2/6.
//

#ifndef RMRENDERER_CAMERA2_H
#define RMRENDERER_CAMERA2_H

#include "common.h"

class ray {
public:
    glm::vec3 rd;
    glm::vec3 ro;
    glm::vec3 at(float t) {
        return this->ro + t * this->rd;
    }
    ray(glm::vec3 ro, glm::vec3 rd) {
        this->ro = ro;
        this->rd = rd;
    }
};

class camera {
public:

    camera(glm::vec3 ta, float camRad, float camHeight, float cr) {
        this->camLoc = glm::vec3 (camRad, camHeight, camRad);
        auto ro = this->camLoc;
        this->z = glm::normalize(ta - ro);
        this->cp = glm::vec3(sin(cr), cos(cr), 0.);
        this->x = glm::normalize(glm::cross(this->z, this->cp));
        this->y = glm::cross(this->x, this->z);
        this->camMat = glm::mat3(this->x, this->y, this->z);
    }

    camera(glm::vec3 ta, glm::vec3 ro, float cr) {
        this->camLoc = ro;
        this->z = glm::normalize(ta - ro);
        this->cp = glm::vec3(sin(cr), cos(cr), 0.);
        this->x = glm::normalize(glm::cross(this->z, this->cp));
        this->y = glm::cross(this->x, this->z);
        this->camMat = glm::mat3(this->x, this->y, this->z);
    }

    ray getRay(glm::vec2 uv) {
        auto rd = glm::normalize(this->camMat * glm::vec3(uv, 1.));
        return ray(this->camLoc, rd);
    }
private:
    glm::vec3 x, y, z, cp, camLoc;
    glm::mat3 camMat;
};

#endif //RMRENDERER_CAMERA2_H
