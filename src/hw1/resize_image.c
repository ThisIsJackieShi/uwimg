#include <math.h>
#include "image.h"
// colaborated with Wang Yuan

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
    return get_pixel(im, round(x), round(y), c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    image result = make_image(w, h, im.c);
    float width_ratio = im.w / (float)w;
    float width_off = 0.5 * (width_ratio - 1);
    float height_ratio = im.h / (float)h;
    float height_off = 0.5 * (height_ratio - 1);
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < result.h; y++) {
            for (int x = 0; x < result.w; x++) {
                float nearest_pixel = 
                    nn_interpolate(im, x * width_ratio + width_off, 
                                    y * height_ratio + height_off, c);
                set_pixel(result, x, y, c, nearest_pixel);
            }
        }
    }
    return result;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    int cordix = x;
    int cordiy = y;
    float top_left = get_pixel(im, x, y, c);
    float top_right = get_pixel(im, x + 1, y, c);
    float bottom_left = get_pixel(im, x, y + 1, c);
    float bottom_right = get_pixel(im, x + 1, y + 1, c);
    
    float left = top_left * (cordiy + 1 - y) + bottom_left * (y - cordiy);
    float right = top_right * (cordiy + 1 - y) + bottom_right * (y - cordiy);
    float result = left * (cordix + 1 - x) + right * (x - cordix);

    return result;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image result = make_image(w, h, im.c);
    float width_ratio = im.w / (float)w;
    float width_off = 0.5 * (width_ratio - 1);
    float height_ratio = im.h / (float)h;
    float height_off = 0.5 * (height_ratio - 1);
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < result.h; y++) {
            for (int x = 0; x < result.w; x++) {
                float nearest_pixel = 
                    bilinear_interpolate(im, x * width_ratio + width_off, 
                                    y * height_ratio + height_off, c);
                set_pixel(result, x, y, c, nearest_pixel);
            }
        }
    }
    return result;
}

