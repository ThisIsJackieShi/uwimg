#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>
#define TWOPI 6.2831853

// collaborated with Jiajie Shi

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    // TODO: optional, make separable 1d Gaussian.
    int size = (int) ceil(sigma * 6.0);
    if (size % 2 == 0) {
        size++;
    }
    image result = make_image(size, 1, 1);
    float denom = 1 / (sigma * sqrt(TWOPI * 2));
    for (int i = 0; i < size; i++) {
        float e =  exp(- (pow(i - size / 2, 2)) / (2 * pow(sigma, 2)));
        set_pixel(result, i, 0, 0, e * denom);
    }
    return result;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    if(0){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: optional, use two convolutions with 1d gaussian filter.
        // If you implement, disable the above if check.
        image g = make_1d_gaussian(sigma);
        image s = convolve_image(im, g, 1);
        image g2 = make_image(1, g.w, 1);
        for (int i = 0; i < g.w; i++) {
            set_pixel(g2, 0, i, 0, get_pixel(g, i, 0, 0));
        }
        image s2 = convolve_image(s, g2, 1);
        free_image(g);
        free_image(s);
        free_image(g2);
        return s2;
    }
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    // TODO: calculate structure matrix for im.
    image dx_filter = make_gx_filter();
    image dy_filter = make_gy_filter();
    image x_gradients = convolve_image(im, dx_filter, 0);
    image y_gradients = convolve_image(im, dy_filter, 0);
    for(int i = 0; i < im.h; i++){
        for(int j = 0; j < im.w; j++){
            set_pixel(S, j, i, 0, pow(get_pixel(x_gradients, j, i, 0), 2));
            set_pixel(S, j, i, 1, pow(get_pixel(y_gradients, j, i, 0), 2));
            set_pixel(S, j, i, 2, get_pixel(x_gradients, j, i, 0) * get_pixel(y_gradients, j, i, 0));
        }
    }
    return convolve_image(S, make_gaussian_filter(sigma), 1);
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    for(int i = 0; i < R.h; i++){
        for(int j = 0; j < R.w; j++){
            float x2 = get_pixel(S, j, i, 0);
            float y2 = get_pixel(S, j, i, 1);
            float xy = get_pixel(S, j, i, 2);
            set_pixel(R, j, i, 0, x2*y2 - xy*xy - 0.06 * pow(x2 + y2, 2));
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    for(int i = 0; i < r.h; i++){
        for(int j = 0; j < r.w; j++){
            float this_value = get_pixel(im, j, i, 0);
            int max = 1;
            for(int m = -w; m < w + 1; m++){
                for(int n = -w; n < w + 1; n++){
                    if(get_pixel(im, j+m, i+n, 0) > this_value){
                        max = 0;
                    }
                }
            }
            if (max == 0) {
                set_pixel(r, j, i, 0, -999999);
            }
        }
    }
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 0;
    for(int i = 0; i < Rnms.h; i++){
        for(int j = 0; j < Rnms.w; j++){
            float this_pixel = get_pixel(Rnms, j, i, 0);
            if(this_pixel > thresh){
                count++;
            }
        }
    }

    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.
    count = 0;
    for(int i = 0; i < Rnms.h; i++){
        for(int j = 0; j < Rnms.w; j++){
            float this_pixel = get_pixel(Rnms, j, i, 0);
            if(this_pixel > thresh){
                d[count] = describe_index(im, i * im.w + j);
                count++;
            }
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
