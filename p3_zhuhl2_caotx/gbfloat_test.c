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

void* Transpose(float* data, unsigned int X, unsigned int Y){
    float* temp_data = (float*)malloc((X * Y * 3) * sizeof(float));
    memcpy(temp_data, data, sizeof(float)* X * Y * 3);  // use temp_data to store the origin data, modify data directlY
    //  the x and y corresponds to the position in the old matrix
    for(int x = 0; x < X; x++){
        for(int y = 0; y < Y; y++){
            for(int c = 0; c < 3; c++){ // c denotes the colour channels that is using
                data[(y + x * Y) * 3 + c] = temp_data[(x + y * X) * 3 + c];
            } 
        } 
    }
    free(temp_data);
    temp_data = NULL;
}
void normalize_FVec_B(FVec v)
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
void normalize_FVec(FVec v)
{
    // float sum = 0.0;
    unsigned int i,j;
    int ext = v.length / 2;
    v.sum[0] = v.data[ext];
    for (i = ext+1,j=1; i < v.length; i++,j++)
    {
        v.sum[j] = v.sum[j-1] + v.data[3*i]*2;
    }
    // for (i = 0; i <= ext; i++)
    // {
    //      v.data[i] /= v.sum[v.length - ext - 1 ] ;
    //      printf("%lf ",v.sum[i]);
    // }
}

float* get_pixel(Image img, int x, int y)
{
    if (x < 0)
    {
        x = 0;
    }
    if (x >= img.dimX)
    {
        x = img.dimX - 1;
    }
    if (y < 0)
    {
        y = 0;
    }
    if (y >= img.dimY)
    {
        y = img.dimY - 1;
    }
    return img.data + img.numChannels * (y * img.dimX + x);
}

float* gget_pixel(float* data,int x,int y,unsigned int X,unsigned int Y){
    if (x < 0){x = 0;}
    if (x >= X){x = X - 1;}
    if (y < 0){y = 0;}
    if (y >= Y){y = Y - 1;}
    return data + 3 * (y * X + x);
}

float gd(float a, float b, float x)
{
    float c = (x-b) / a;
    return exp((-.5) * c * c) / (a * sqrt(2 * PI));
}
FVec make_gv_B(float a, float x0, float x1, unsigned int length, unsigned int min_length)
{
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

    for (int i = 0; i < length; i++)
    {
        v.data[i] = gd(a, 0.0f, (i-offset)*step);
    }
    normalize_FVec_B(v);
    return v;
}
FVec make_gv(float a, float x0, float x1, unsigned int length, unsigned int min_length)
{
    FVec v;
    v.length = length;
    v.min_length = min_length;
    if(v.min_length > v.length){
        v.min_deta = 0;
    }else{
        v.min_deta = ((v.length - v.min_length) / 2);
    }
    v.data = malloc(3 * length * sizeof(float));
    v.sum = malloc((length / 2 + 1)* sizeof(float));
    float step = (x1 - x0) / ((float)length);
    int offset = length/2;

    for (int i = 0; i < length; i++)
    {
        v.data[3*i] = gd(a, 0.0f, (i-offset)*step);
        v.data[3*i+1] = gd(a, 0.0f, (i-offset)*step);
        v.data[3*i+2] = gd(a, 0.0f, (i-offset)*step);
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

Image gb_h(Image a, FVec gv)
{
    Image b = img_sc(a);

    int ext = gv.length / 2;
    int offset;
    int x, y, channel;
    float *pc;
    float sum;
    int i;

    int X = b.dimX;
    int Y = b.dimY;

    float* data_origin = a.data;
    float* data_gb_h = (float*)malloc(X*Y*3*sizeof(float));
    float* data_out = b.data;



    for (y = 0; y < a.dimY; y++)
    {
        for (x = 0; x < a.dimX; x++)
        {
            pc = gget_pixel(data_gb_h, x, y, X,Y);
            int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
            deta = fmin(deta, gv.min_deta);
            float sum0 = 0, sum1 = 0, sum2 = 0;

            int aa,bb; // Left and right bounds.
            aa = deta-ext+x < 0    ?  0   :   x-(ext-deta)  ;
            bb = x+ext-deta > X-1  ?  X-1 :   x+(ext-deta)  ;
            
            int l = aa - x + ext;
            int r = bb - x + ext;

            // __m256 vec0, vec1, vec2, p0, p1, p2;
            // for (i = l; i < l+(r-l)/8*8; i+=8)
            // {
            //     offset = i - ext;
            //     //sum += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel(a, x + offset, y)[channel];

            //     // gv.data
            //     vec0 = _mm256_loadu_ps(gv.data+3*i);
            //     vec1 = _mm256_loadu_ps(gv.data+3*i+8);
            //     vec2 = _mm256_loadu_ps(gv.data+3*i+16);

            //     // pixels, data from original image
            //     p0 = _mm256_loadu_ps(gget_pixel(data_origin, x + offset, y, X, Y));
            //     p1 = _mm256_loadu_ps(gget_pixel(data_origin, x + offset, y, X, Y) + 8);
            //     p2 = _mm256_loadu_ps(gget_pixel(data_origin, x + offset, y, X, Y) + 16);

            //     // Multiplication
            //     vec0 = _mm256_mul_ps(vec0, p0);
            //     vec1 = _mm256_mul_ps(vec1, p1);
            //     vec2 = _mm256_mul_ps(vec2, p2);

            //     // Add to sum_vec
            //     sum0 += vec0[0]+vec0[3]+vec0[6]+vec1[1]+vec1[4]+vec1[7]+vec2[2]+vec2[5];
            //     sum1 += vec0[1]+vec0[4]+vec0[7]+vec1[2]+vec1[5]+vec2[0]+vec2[3]+vec2[6];
            //     sum2 += vec0[2]+vec0[5]+vec1[0]+vec1[3]+vec1[6]+vec2[1]+vec2[4]+vec2[7];
            // }


            // Tail cases
            for (i=l; i <= r; i++){
                offset = i - ext;
                sum0 += gv.data[i]*gget_pixel(data_origin,x+offset,y,X,Y)[0];
                sum1 += gv.data[i]*gget_pixel(data_origin,x+offset,y,X,Y)[1];
                sum2 += gv.data[i]*gget_pixel(data_origin,x+offset,y,X,Y)[2];
            }
            
            

            // Left and right edges
            if(aa==0){
                float temp=0;
                for(int ii=deta;ii<l;ii++){
                    temp+=gv.data[ii];
                }
                sum0 += temp*gget_pixel(data_origin,0,y, X, Y)[0];
                sum1 += temp*gget_pixel(data_origin,0,y, X, Y)[1];
                sum2 += temp*gget_pixel(data_origin,0,y, X, Y)[2];
                
            }
            if(bb==X-1){
                float temp=0;
                for(int ii=r+1;ii<gv.length-deta;ii++){
                    temp+=gv.data[ii];
                }
                sum0 += temp*gget_pixel(data_origin,X-1,y, X, Y)[0];
                sum1 += temp*gget_pixel(data_origin,X-1,y, X, Y)[1];
                sum2 += temp*gget_pixel(data_origin,X-1,y, X, Y)[2];
                
            }


            pc[0] = sum0 / gv.sum[ext - deta];
            pc[1] = sum1 / gv.sum[ext - deta];
            pc[2] = sum2 / gv.sum[ext - deta];
        }
    }

    b.data=data_gb_h;
    return b;
}



Image apply_gb(Image a, FVec gv)
{
    Image b = gb_h(a, gv);
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

    FVec v = make_gv_B(a, x0, x1, dim, min_dim);
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