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
        this->y_offset = -2.;
        ro = glm::vec3(0., 0. + y_offset, -1.5);
    }

    ray ray(glm::vec2 uv) {
        //
        glm::vec3 rd = glm::normalize(glm::vec3(uv + glm::vec2(0, 0. + this->y_offset), 0) - ro);
        return {ro, rd};
    }
private:
    glm::vec3 ro;
    // 控制镜头的垂直位置
    float y_offset;
};

#endif //RMRENDERER_CAMERA_H
