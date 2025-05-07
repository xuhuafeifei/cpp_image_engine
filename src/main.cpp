#include <iostream>
#include <vector>
#include "render.hpp"

std::string path = "E:\\c++_code\\rmRenderer\\picts3\\";

int main(int argc, char **argv)
{
    std::vector<Color> pixels(width * height);
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y);

    // ✅ 保存为 PNG 图片
    std::string filename = path + "kang_ju_chi.png";

    InitWindow(width, height, "Render Result");


    Image image = {
        pixels.data(),
        width,
        height,
        TEXTURE_FILTER_BILINEAR,
        PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
    };

    Texture2D texture2 = LoadTextureFromImage(image);

    std::cout << filename << std::endl;

    ExportImage(image, filename.c_str());

    CloseWindow();

    return 0;

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
