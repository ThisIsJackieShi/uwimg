#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

typedef struct{
    float r, g, b;
} rgb_value;

typedef struct{
    float h, s, v;
} hsv_value;

// fit x into 0 to max - 1.
int fit(int x, int max);

// return the offset of the pixel of the coordinate in the image
int get_offset(image im, int x, int y, int c);

// get hsv values from rgb
hsv_value get_hsv_value(float r, float g, float b);

// get rgb values from hsv
rgb_value get_rgb_value(float h, float s, float v);

float get_pixel(image im, int x, int y, int c)
{
    x = fit(x, im.w);
    y = fit(y, im.h);
    c = fit(c, im.c);
    return im.data[get_offset(im, x, y, c)];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // invalid coordinate
    if (x < 0 || x >= im.w || y < 0 || y >= im.h || c < 0 || c >= im.c) {
        printf("invalid coordinate\n");
        return;
    }
    im.data[get_offset(im, x, y, c)] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));    
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // weights for computing gray scale, in RGB order
    float weight[3] = {0.299, 0.587, 0.114};
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < im.h; y++) {
            for (int x = 0; x < im.w; x++) {
                gray.data[get_offset(gray, x, y, 0)] += get_pixel(im, x, y, c) * weight[c];
            }
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    if (c < 0 || c >= im.c) {
        return;
    }
    for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
            float val = get_pixel(im, x, y, c) + v;
            set_pixel(im, x, y, c, val);
        }
    }
}

void clamp_image(image im)
{
    for (int i = 0; i < im.w * im.h * im.c; i++) {
        im.data[i] = MIN(MAX(im.data[i], 0), 1);
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
            hsv_value hsv = get_hsv_value(get_pixel(im, x, y, 0), get_pixel(im, x, y, 1), get_pixel(im, x, y, 2));
            set_pixel(im, x, y, 0, hsv.h);
            set_pixel(im, x, y, 1, hsv.s);
            set_pixel(im, x, y, 2, hsv.v);
        }
    }
}

void hsv_to_rgb(image im)
{
    for (int y = 0; y < im.h; y++) {
        for (int x = 0; x < im.w; x++) {
            rgb_value rgb = get_rgb_value(get_pixel(im, x, y, 0), get_pixel(im, x, y, 1), get_pixel(im, x, y, 2));
            set_pixel(im, x, y, 0, rgb.r);
            set_pixel(im, x, y, 1, rgb.g);
            set_pixel(im, x, y, 2, rgb.b);
        }
    }
}

int fit(int x, int max) 
{
    return MIN(MAX(0, x), max - 1);
}

int get_offset(image im, int x, int y, int c) 
{
    return c * (im.w * im.h) + y * (im.w) + x;
}

hsv_value get_hsv_value(float r, float g, float b)
{
    float v = three_way_max(r, g, b);
    float m = three_way_min(r, g, b);
    float C = v - m;
    float S = 0;
    if (v != 0) {
        S = C / v;
    }
    float h0 = 0;
    if (C != 0) {
        if (v == r) {
            h0 = (g - b) / C;
        } else if (v == g) {
            h0 = (b - r) / C + 2;
        } else {  // v == b
            h0 = (r - g) / C + 4;
        }
    }
    float H = h0 / 6.0;
    if (H < 0) {
        H += 1.0;
    }
    hsv_value res;
    res.h = H;
    res.s = S;
    res.v = v;
    return res;
}

// get rgb values from hsv
rgb_value get_rgb_value(float h, float s, float v) {
    // formula from https://www.rapidtables.com/convert/color/hsv-to-rgb.html
    rgb_value res;
    float C = v * s;
    float hue6 = h * 6.0;

    int int_hue = (int) hue6;
    assert(0 <= int_hue && int_hue < 6);
    while (hue6 > 2) {
        hue6 -= 2;
    }
    hue6 -= 1;
    if (hue6 < 0) {
        hue6 = -hue6;
    }
    float X = C * (1 - hue6);
    float m = v - C;
    float r0, g0, b0;
    switch (int_hue) {
        case 0:
            r0 = C;
            g0 = X;
            b0 = 0;
            break;
        case 1:
            r0 = X;
            g0 = C;
            b0 = 0;
            break;
        case 2:
            r0 = 0;
            g0 = C;
            b0 = X;
            break;
        case 3:
            r0 = 0;
            g0 = X;
            b0 = C;
            break;
        case 4:
            r0 = X;
            g0 = 0;
            b0 = C;
            break;
        case 5:
            r0 = C;
            g0 = 0;
            b0 = X;
            break;

        default: 
            printf("invalid hsv");
            return res;
    }

    res.r = r0 + m;
    res.g = g0 + m;
    res.b = b0 + m;

    return res;
}