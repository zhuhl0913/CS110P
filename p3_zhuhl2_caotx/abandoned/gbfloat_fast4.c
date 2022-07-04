#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <algorithm>
#include <sys/time.h>
#include <time.h>
#include <immintrin.h>
//inplement dymanic

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define PI 3.14159

typedef struct FVec
{
    unsigned int length;
    unsigned int min_length;
    unsigned int min_deta;
    float* data;
    float* sum;
} FVec;

typedef struct Image
{
    unsigned int dimX, dimY, numChannels;
    float* data;
} Image;

void normalize_FVec(FVec v)
{
    // float sum = 0.0;
    unsigned int i,j;
    int ext = v.length / 2;
    v.sum[0] = v.data[ext];
    for (i = ext+1,j=1; i < v.length; i++,j++)
    {
        v.sum[j] = v.sum[j-1] + v.data[i]*2;
    }
    // for (i = 0; i <= ext; i++)
    // {
    //      v.data[i] /= v.sum[v.length - ext - 1 ] ;
    //      printf("%lf ",v.sum[i]);
    // }
}

float* get_pixel(Image img, int x, int y){
    if (x < 0)          {x = 0;}
    if (x >= img.dimX)  {x = img.dimX - 1;}
    if (y < 0)          {y = 0;}
    if (y >= img.dimY)  {y = img.dimY - 1;}
    return img.data + img.numChannels * (y * img.dimX + x);
}

float* get_pixel_array(float* data, int x, int y, unsigned X,unsigned int Y){
    if (x < 0)      {x = 0;}
    if (x >= X)     {x = X - 1;}
    if (y < 0)      {y = 0;}
    if (y >= Y)     {y = Y- 1;}
    return data + (y * X + x);
}

float gd(float a, float b, float x){
    float c = (x-b) / a;
    return exp((-.5) * c * c) / (a * sqrt(2 * PI));
}

FVec make_gv(float a, float x0, float x1, unsigned int length, unsigned int min_length){
    FVec v;
    v.length = length;
    v.min_length = min_length;
    if(v.min_length > v.length){
        v.min_deta = 0;
    }else{
        v.min_deta = ((v.length - v.min_length) / 2);
    }
    v.data = malloc(length * sizeof(float));
    v.sum = malloc((length / 2 + 1)* sizeof(float));
    float step = (x1 - x0) / ((float)length);
    int offset = length/2;

    # pragma omp parallel for
    for (int i = 0; i < length; i++)
    {
        v.data[i] = gd(a, 0.0f, (i-offset)*step);
    }
    normalize_FVec(v);
    return v;
}

Image img_sc(Image a)
{
    Image b = a;
    b.data = malloc(b.dimX * b.dimY * b.numChannels * sizeof(float));
    return b;
}

Image apply_gb(Image a, FVec gv)
{
    Image b = img_sc(a);

    int ext = gv.length / 2;
    int offset;
    unsigned int x, y, channel;
    float *pc;
    float sum;
    int i;

    float **data_a=(float**)malloc(a.numChannels*sizeof(float*));
    float **data_b_T=(float**)malloc(a.numChannels*sizeof(float*));

    #pragma omp parallel for
    for(i=0;i<a.numChannels;i++){
        data_a[i]=(float*)malloc(a.dimX*a.dimY*sizeof(float));
        data_b_T[i]=(float*)malloc(a.dimX*a.dimY*sizeof(float));
        for(int n=0;n<a.dimX*a.dimY;n++){
            data_a[i][n]=a.data[i+a.numChannels*n];
        }
    }

    unsigned int X = a.dimX;
    unsigned int Y = a.dimY;

    // data_a: data of original image but in CHW.
    // data_b_T: data of image after horizonal convolution, and transposed.


    // Do convolution on data_a and store the outcome in a transposed way in data_b_T
    
    # pragma omp parallel for collapse (2)
    for (channel = 0; channel < a.numChannels; channel++){
        for (y = 0; y < Y; y++){
            for (x = 0; x < X; x++){
                pc = get_pixel_array(data_b_T[channel],y,x,Y,X);
                unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
                deta = fmin(deta, gv.min_deta);
                sum = 0;
                //// ORIGINAL
                // for (i = deta; i < gv.length-deta; i++){
                //     offset = i - ext;
                //     sum += gv.data[i]/gv.sum[ext - deta] * get_pixel_array(data_a[channel], x + offset, y,X,Y)[0];
                // }
                // pc[0] = sum;

                //// LOOP UNROLLING
                for(i=deta;i<deta+(gv.length-2*deta)/8*8;i+=8){
                    offset=i-ext;
                    pc[0]+=gv.data[i]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset,y,a.dimX,a.dimY)[0]
                    +gv.data[i+1]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+1,y,a.dimX,a.dimY)[0]
                    +gv.data[i+2]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+2,y,a.dimX,a.dimY)[0]
                    +gv.data[i+3]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+3,y,a.dimX,a.dimY)[0]
                    +gv.data[i+4]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+4,y,a.dimX,a.dimY)[0]
                    +gv.data[i+5]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+5,y,a.dimX,a.dimY)[0]
                    +gv.data[i+6]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+6,y,a.dimX,a.dimY)[0]
                    +gv.data[i+7]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset+7,y,a.dimX,a.dimY)[0];
                }
                for(;i<gv.length-deta;i++){
                    offset=i-ext;
                    pc[0]+=gv.data[i]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset,y,a.dimX,a.dimY)[0];
                }
            }
        }
    }

    // Do convolution on data_b_T and store the outcome in (Image)b

    # pragma omp parallel for collapse (2)
    for (channel = 0; channel < a.numChannels; channel++){
        for (y = 0; y < Y; y++){
            for (x = 0; x < X; x++){
                pc = get_pixel(b,x,y);
                unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
                deta = fmin(deta, gv.min_deta);
                sum = 0;
                //// ORIGINAL
                // for (i = deta; i < gv.length-deta; i++){
                //     offset = i - ext;
                //     sum += gv.data[i]/gv.sum[ext - deta] * get_pixel_array(data_b_T[channel], y + offset, x,Y,X)[0];
                // }
                // pc[channel] = sum;

                //// LOOP UNROLLING
                for(i=deta;i<deta+(gv.length-2*deta)/8*8;i+=8){
                    offset=i-ext;
                    pc[channel]+=gv.data[i]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset,x,Y,X)[0]
                    +gv.data[i+1]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+1,x,Y,X)[0]
                    +gv.data[i+2]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+2,x,Y,X)[0]
                    +gv.data[i+3]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+3,x,Y,X)[0]
                    +gv.data[i+4]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+4,x,Y,X)[0]
                    +gv.data[i+5]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+5,x,Y,X)[0]
                    +gv.data[i+6]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+6,x,Y,X)[0]
                    +gv.data[i+7]/gv.sum[ext-deta]*get_pixel_array(data_b_T[channel],y+offset+7,x,Y,X)[0];
                }
                for(;i<gv.length-deta;i++){
                    offset=i-ext;
                    pc[channel]+=gv.data[i]/gv.sum[ext-deta]*get_pixel_array(data_a[channel],x+offset,y,a.dimX,a.dimY)[0];
                }

            }
        }
    }

    return b;
}

int main(int argc, char** argv)
{
    struct timeval start_time, stop_time, elapsed_time; 
    gettimeofday(&start_time,NULL);
    if (argc < 6)
    {
        printf("Usage: ./gb.exe <inputjpg> <outputname> <float: a> <float: x0> <float: x1> <unsigned int: dim>\n");
        exit(0);
    }

    float a, x0, x1;
    unsigned int dim, min_dim;

    sscanf(argv[3], "%f", &a);
    sscanf(argv[4], "%f", &x0);
    sscanf(argv[5], "%f", &x1);
    sscanf(argv[6], "%u", &dim);
    sscanf(argv[7], "%u", &min_dim);

    FVec v = make_gv(a, x0, x1, dim, min_dim);
    // print_fvec(v);
    Image img;
    img.data = stbi_loadf(argv[1], &(img.dimX), &(img.dimY), &(img.numChannels), 0);

    Image imgOut = apply_gb(img, v);
    stbi_write_jpg(argv[2], imgOut.dimX, imgOut.dimY, imgOut.numChannels, imgOut.data, 90);
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); 
    printf("%f \n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    free(imgOut.data);
    free(v.data);
    free(v.sum);
    return 0;
}