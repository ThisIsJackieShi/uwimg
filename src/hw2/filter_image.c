#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
    int w = im.w;
    int h = im.h;
    for(int c = 0; c < im.c; c++){
        for(int i = 0; i < h; i++){
            for(int j = 0; j < w; j++){
                set_pixel(im, j, i, c, get_pixel(im, j, i, c)/(w * h));
            }
        }
    }
}

image make_box_filter(int w)
{
    // TODO
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
    // TODO
    image result = make_image(im.w, im.h, preserve == 1? im.c : 1);
    assert(filter.c == im.c || filter.c == 1);
    if(preserve == 1){
        for(int c = 0; c < result.c; c++){
            for(int i = 0; i < result.h; i++){
                for(int j = 0; j < result.w; j++){
                    float convolved = 0;
                    for(int m = -filter.w/2; m < filter.w/2 + 1; m++){
                        for(int n = -filter.w/2; n < filter.w/2 + 1; n++){
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
                for(int c = 0; c < result.c; c++){
                    for(int m = -filter.w/2; m < filter.w/2 + 1; m++){
                        for(int n = -filter.w/2; n < filter.w/2 + 1; n++){
                            convolved += get_pixel(im, j+m, i+n, c) * get_pixel(filter, filter.w/2+m, filter.h/2+n, filter_channel_fit(filter, c));
                        }
                    }
                }
                set_pixel(result, j, i, 1, convolved);
            }
        }
    }
    return result;
}

image make_highpass_filter()
{
    // TODO
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
    // TODO
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
    // TODO
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
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    return make_image(1,1,1);
}

image add_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO
    return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
}

image *sobel_image(image im)
{
    // TODO
    return calloc(2, sizeof(image));
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
