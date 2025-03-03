#include <iostream>
#include <vector>
#include <string>
#include <thread>

//#define MOUNTAIN_DRAW

# ifdef MOUNTAIN_DRAW
//    # define USE_IMAGE_EXPORT
    # define REVERSE_Y
    # include "render3.hpp"
# else
    # include "render.hpp"
# endif

void drawP(std::vector<Color> &pixels, int p, int tot, float iTime)
{
    int seg = height / tot;
    for (int y = seg * p; y < seg * (p + 1); y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y, iTime);
}

int main(int argc, char **argv) {
    std::vector<Color> pixels(width * height);

#ifdef USE_IMAGE_EXPORT
    std::string path = "E:\\c++_code\\rmRenderer\\picts\\";
    int count = 0;
    for (float iTime = 1.5; iTime <= 24.5 * 1000; iTime += 0.2) {
        std::thread a(drawP, std::ref(pixels), 0, 3, iTime);
        std::thread b(drawP, std::ref(pixels), 1, 3, iTime);
        std::thread c(drawP, std::ref(pixels), 2, 3, iTime);

        a.join();
        b.join();
        c.join();

        std::cerr << "iTime: " << iTime << std::endl;

        // 上下反转
        for (int y = 0; y < height / 2; y++)
            for (int x = 0; x < width; x++)
                std::swap(pixels[y * width + x], pixels[(height - y - 1) * width + x]);

        std::cerr << "swap: " << iTime << std::endl;

        InitWindow(width, height, "raylib - Save Image Example");

        Image image = {
                pixels.data(),
                width,
                height,
                TEXTURE_FILTER_BILINEAR,
                PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
        };

        Texture2D texture = LoadTextureFromImage(image);

        std::cerr << "Image_load.." << std::endl;

        auto fileName = path + "rm_" + std::to_string(count) + ".png";
        std::cerr << fileName << std::endl;
        ExportImage(image, fileName.c_str());

        std::cerr << "Exported image to: " << fileName << std::endl;

        CloseWindow();
        count += 1;
    }
#else
    std::thread a(drawP, std::ref(pixels), 0, 3, 1.5);
    std::thread b(drawP, std::ref(pixels), 1, 3, 1.5);
    std::thread c(drawP, std::ref(pixels), 2, 3, 1.5);

    a.join();
    b.join();
    c.join();

#ifdef REVERSE_Y
    // 上下反转
    for (int y = 0; y < height / 2; y++)
        for (int x = 0; x < width; x++)
            std::swap(pixels[y * width + x], pixels[(height - y - 1) * width + x]);
#endif

    InitWindow(width, height, "Render Result");

    Image image = {
        pixels.data(),
        width,
        height,
        TEXTURE_FILTER_BILINEAR,
        PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
    };

    Texture2D texture = LoadTextureFromImage(image);

    while (!WindowShouldClose()) {
        BeginDrawing();
        DrawTexture(texture, 0, 0, WHITE);
        EndDrawing();
    }

    CloseWindow();
#endif

    return 0;
}
