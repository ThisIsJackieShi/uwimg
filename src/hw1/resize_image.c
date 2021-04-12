#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    int xRound = (int) roundf(x);
    int yRound = (int) roundf(y);
    return get_pixel(im, x, y, c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
    return make_image(1,1,1);
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
    return make_image(1,1,1);
}

