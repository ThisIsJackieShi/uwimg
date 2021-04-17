#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853
#define PI 3.1415926

void l1_normalize(image im)
{
    int w = im.w;
    int h = im.h;
    float sum = 0;
    for(int c = 0; c < im.c; c++){
        for(int i = 0; i < h; i++){
            for(int j = 0; j < w; j++){
                sum += get_pixel(im, j, i, c);
            }
        }
    }
    if (sum == 0.0) {
        return;
    }
    for(int c = 0; c < im.c; c++){
        for(int i = 0; i < h; i++){
            for(int j = 0; j < w; j++){
                set_pixel(im, j, i, c, get_pixel(im, j, i, c) / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    image result = make_image(w, w, 1);
    for(int i = 0; i < w; i++){
        for(int j = 0; j < w; j++){
            set_pixel(result, j, i, 0, 1.0);
        }
    }
    l1_normalize(result);
    return result;
}

int filter_channel_fit(image filter, int channel){
    return filter.c == 1? 0:channel;
}

image convolve_image(image im, image filter, int preserve)
{
    image result = make_image(im.w, im.h, preserve == 1? im.c : 1);
    assert(filter.c == im.c || filter.c == 1);
    if(preserve == 1){
        for(int c = 0; c < result.c; c++){
            for(int i = 0; i < result.h; i++){
                for(int j = 0; j < result.w; j++){
                    float convolved = 0;
                    for(int m = -filter.w/2; m < filter.w/2 + 1; m++){
                        for(int n = -filter.h/2; n < filter.h/2 + 1; n++){
                            convolved += get_pixel(im, j+m, i+n, c) * get_pixel(filter, filter.w/2+m, filter.h/2+n, filter_channel_fit(filter, c));
                        }
                    }
                    set_pixel(result, j, i, c, convolved);
                }
            }
        }
    } else { // smash pic to one
        for(int i = 0; i < result.h; i++){
            for(int j = 0; j < result.w; j++){
                float convolved = 0;
                for(int c = 0; c < im.c; c++){
                    for(int m = -filter.w/2; m < filter.w/2 + 1; m++){
                        for(int n = -filter.h/2; n < filter.h/2 + 1; n++){
                            convolved += get_pixel(im, j+m, i+n, c) * get_pixel(filter, filter.w/2+m, filter.h/2+n, filter_channel_fit(filter, c));
                        }
                    }
                }
                set_pixel(result, j, i, 0, convolved);
            }
        }
    }
    return result;
}

image make_highpass_filter()
{
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, 0);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 0);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 1, 2, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 1, 1, 0, 4);
    return filter;
}

image make_sharpen_filter()
{
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, 0);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 0);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 1, 2, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 1, 1, 0, 5);
    return filter;
}

image make_emboss_filter()
{
    image filter = make_box_filter(3);
    set_pixel(filter, 0, 0, 0, -2);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 2);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 1, 2, 0, 1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, 1);
    set_pixel(filter, 1, 1, 0, 1);
    return filter;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: 
// Perserve: box-blur filter, sharpen filter, emboss filter; because we want to keep the colors of the graph
// Not Perserve: highpass filter. Because we only want to find the edges and their locations. 

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Highpass filter. Since it is smashed to one channel, we might need to process it to black-and-white in order to view it. 

float gaussian(float x, float y, float sigma) 
{
    float e = exp(- (pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2)));
    return e / (TWOPI * pow(sigma, 2));
}

image make_gaussian_filter(float sigma)
{
    int size = (int) ceil(sigma * 6.0);
    if (size % 2 == 0) {
        size++;
    }

    image filter = make_image(size, size, 1);
    for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
            set_pixel(filter, x, y, 0, gaussian(x - (size / 2), y - (size / 2), sigma));
        }
    }
    l1_normalize(filter);

    return filter;
}

image add_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image result = make_image(a.w, a.h, a.c);
    for (int c = 0; c < a.c; c++) {
        for (int h = 0; h < a.h; h++) {
            for (int w = 0; w < a.w; w++) {
                set_pixel(result, w, h, c, get_pixel(a, w, h, c) + get_pixel(b, w, h, c));
            }
        }
    }
    return result;
}

image sub_image(image a, image b)
{
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image result = make_image(a.w, a.h, a.c);
    for (int c = 0; c < a.c; c++) {
        for (int h = 0; h < a.h; h++) {
            for (int w = 0; w < a.w; w++) {
                set_pixel(result, w, h, c, get_pixel(a, w, h, c) - get_pixel(b, w, h, c));
            }
        }
    }
    return result;
}

image make_gx_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 0, 0, 0, -1);
    set_pixel(filter, 1, 0, 0, 0);
    set_pixel(filter, 2, 0, 0, 1);
    set_pixel(filter, 0, 1, 0, -2);
    set_pixel(filter, 1, 1, 0, 0);
    set_pixel(filter, 2, 1, 0, 2);
    set_pixel(filter, 0, 2, 0, -1);
    set_pixel(filter, 1, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 1);
    return filter;
}

image make_gy_filter()
{
    // TODO
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 0, 0, 0, -1);
    set_pixel(filter, 1, 0, 0, -2);
    set_pixel(filter, 2, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, 0);
    set_pixel(filter, 1, 1, 0, 0);
    set_pixel(filter, 2, 1, 0, 0);
    set_pixel(filter, 0, 2, 0, 1);
    set_pixel(filter, 1, 2, 0, 2);
    set_pixel(filter, 2, 2, 0, 1);
    return filter;
}

void feature_normalize(image im)
{
    float min_val = get_pixel(im, 0, 0, 0);
    float max_val = get_pixel(im, 0, 0, 0);
    for (int c = 0; c < im.c; c++) {
        for (int h = 0; h < im.h; h++) {
            for (int w = 0; w < im.w; w++) {
                min_val = min_val > get_pixel(im, w, h, c) ? get_pixel(im, w, h, c) : min_val;
                max_val = max_val < get_pixel(im, w, h, c) ? get_pixel(im, w, h, c) : max_val;
            }
        }
    }
    
    for (int c = 0; c < im.c; c++) {
        for (int h = 0; h < im.h; h++) {
            for (int w = 0; w < im.w; w++) {
                if (max_val == min_val) {
                    set_pixel(im, w, h, c, 0);
                } else {
                    float p = get_pixel(im, w, h, c);
                    p = (p - min_val) / (max_val - min_val);
                    set_pixel(im, w, h, c, p);
                }
            }
        }
    }
}

image *sobel_image(image im)
{
    image horizontal = convolve_image(im, make_gx_filter(), 0);
    image vertical = convolve_image(im, make_gy_filter(), 0);

    image* result = calloc(2, sizeof(image));
    result[0] = make_image(im.w, im.h, 1);
    result[1] = make_image(im.w, im.h, 1);

    for (int h = 0; h < im.h; h++) {
        for (int w = 0; w < im.w; w++) {
            float x = get_pixel(horizontal, w, h, 0);
            float y = get_pixel(vertical, w, h, 0);
            set_pixel(result[0], w, h, 0, sqrt(pow(x, 2) + pow(y, 2)));
            set_pixel(result[1], w, h, 0, atan2(y, x));
        }
    }

    free_image(horizontal);
    free_image(vertical);

    return result;
}

image colorize_sobel(image im)
{
    image result = make_image(im.w, im.h, 3);
    image* sobels = sobel_image(im);
    feature_normalize(sobels[0]);
    rgb_to_hsv(result);
    for (int h = 0; h < result.h; h++) {
        for (int w = 0; w < result.w; w++) {
            float angle = (get_pixel(sobels[1], w, h, 0) + (PI / 2)) / PI;
            float mag = get_pixel(sobels[0], w, h, 0);
            set_pixel(result, w, h, 0, angle);  // angel for hue
            set_pixel(result, w, h, 1, mag);    // mag for saturation
            set_pixel(result, w, h, 2, 1);      // 1 for value
        }
    }
    hsv_to_rgb(result);

    free_image(sobels[0]);
    free_image(sobels[1]);
    free(sobels);

    return result;
}
