#include <iostream>
#include <vector>
#include <string>
//#include "render.hpp"
//#include "render2.hpp"
//#include "render4.hpp"
//#include "render5.hpp"
#include "render3.hpp"
#include <thread>

void drawP(std::vector<Color> &pixels, int p, int tot, float iTime)
{
    int seg = height / tot;
    for (int y = seg * p; y < seg * (p + 1); y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y, iTime);
}

//int main()
//{
//
//    std::vector <Color> pixels(width * height);
//    std::string path = "E:\\c++_code\\rmRenderer\\picts\\";
//    for (float iTime = 1.5; iTime <= 24.5 * 1000; iTime += 0.2) {
//        std::thread a(drawP, std::ref(pixels), 0, 3, iTime);
//        std::thread b(drawP, std::ref(pixels), 1, 3, iTime);
//        std::thread c(drawP, std::ref(pixels), 2, 3, iTime);
//
//        a.join();
//        b.join();
//        c.join();
//
//        std::cerr << "iTime: " << iTime << std::endl;
//
//        // 上下反转
//        for (int y = 0; y < height / 2; y++)
//            for (int x = 0; x < width; x++)
//                std::swap(pixels[y * width + x], pixels[(height - y - 1) * width + x]);
//
//        std::cerr << "swap: " << iTime << std::endl;
//
//        InitWindow(width, height, "raylib - Save Image Example");
//
//        Image image = {
//                pixels.data(),
//                width,
//                height,
//                TEXTURE_FILTER_BILINEAR,
//                PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
//        };
//
//        Texture2D texture = LoadTextureFromImage(image);
//
//        std::cerr << "Image_load.." << std::endl;
//
//        auto fileName = path + "rm_" + std::to_string(iTime) + ".png";
//        std::cerr << fileName << std::endl;
//        ExportImage(image, fileName.c_str());
//
//        std::cerr << "Exported image to: " << fileName << std::endl;
//
//        CloseWindow();
//    }
//}

int main(int argc, char **argv)
{
    std::vector<Color> pixels(width * height);

    std::thread a(drawP, std::ref(pixels), 0, 3, 1.5);
    std::thread b(drawP, std::ref(pixels), 1, 3, 1.5);
    std::thread c(drawP, std::ref(pixels), 2, 3, 1.5);

    a.join();
    b.join();
    c.join();

    // 上下反转
    for (int y = 0; y < height / 2; y++)
        for (int x = 0; x < width; x++)
            std::swap(pixels[y * width + x], pixels[(height - y - 1) * width + x]);

    InitWindow(width, height, "Render Result");

    Image image = {
        pixels.data(),
        width,
        height,
        TEXTURE_FILTER_BILINEAR,
        PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
    };

    Texture2D texture = LoadTextureFromImage(image);

    while (!WindowShouldClose())
    {
        BeginDrawing();
        DrawTexture(texture, 0, 0, WHITE);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
