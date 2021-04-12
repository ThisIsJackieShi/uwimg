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
    int x_min = (int) floorf(x);
    int x_max = x_min + 1;
    int y_min = (int) floorf(y);
    int y_max = y_min + 1;

    float q1 = (y - y_min) * get_pixel(im, x_min, y_min, c) + (y_max - y) * get_pixel(im, x_min, y_max, c);
    float q2 = (y - y_min) * get_pixel(im, x_max, y_min, c) + (y_max - y) * get_pixel(im, x_max, y_max, c);

    return (x - x_min) * q1 + (x_max - x) * q2;
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

