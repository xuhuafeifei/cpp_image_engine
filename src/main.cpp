#include <iostream>
#include <vector>
//#include "render.hpp"
//#include "render2.hpp"
//#include "render4.hpp"
//#include "render5.hpp"
#include "render3.hpp"

int main(int argc, char **argv)
{
    std::vector<Color> pixels(width * height);
    for (int y = 0; y < height; y++)
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = render(x, y);

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
