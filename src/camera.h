//
// Created by 25080 on 2025/2/4.
//

#ifndef RMRENDERER_CAMERA_H
#define RMRENDERER_CAMERA_H

#include "common.h"

class ray {
public:
    ray(glm::vec3 ro, glm::vec3 rd) {
        this->ro = ro;
        this->rd = rd;
    }

    glm::vec3 at(float t) {
        return ro + t * rd;
    }
public:
    glm::vec3 ro;
    glm::vec3 rd;
};

class camera {
public:
    camera() {
        ro = glm::vec3(0., 0, -1.5);
    }

    ray ray(glm::vec2 uv) {
        glm::vec3 rd = glm::normalize(glm::vec3(uv, 0) - ro);
        return {ro, rd};
    }
private:
    glm::vec3 ro;
};

#endif //RMRENDERER_CAMERA_H
