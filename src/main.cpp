#include <iostream>
#include <vector>
//#include "render.hpp"
//#include "render2.hpp"
//#include "render4.hpp"
//#include "render5.hpp"
#include "render3.hpp"
#include <thread>

void drawHalf(std::vector<Color> &pixels)
{
    for (int y = 0; y < height / 2; y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y);
}

void drawHalf2(std::vector<Color> &pixels)
{
    for (int y = height / 2 + 1; y < height; y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y);
}

void drawP(std::vector<Color> &pixels, int p, int tot)
{
    int seg = height / tot;
    for (int y = seg * p; y < seg * (p + 1); y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y);
}

int main(int argc, char **argv)
{
    std::vector<Color> pixels(width * height);

    std::thread a(drawP, std::ref(pixels), 0, 3);
    std::thread b(drawP, std::ref(pixels), 1, 3);
    std::thread c(drawP, std::ref(pixels), 2, 3);

    a.join();
    b.join();
    c.join();

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
