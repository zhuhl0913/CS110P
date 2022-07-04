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

float* get_pixel_Transpose(Image img, int x, int y)
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
    return img.data + img.numChannels * (x * img.dimY + y);
}

float* get_pixel2(float* data,int x,int y,unsigned int dimX,unsigned int dimY){
    if (x < 0)
    {
        x = 0;
    }
    if (x >= dimX)
    {
        x = dimX - 1;
    }
    if (y < 0)
    {
        y = 0;
    }
    if (y >= dimY)
    {
        y = dimY - 1;
    }
    return data + (y * dimX + x);
}

float* get_pixel2_Transpose(float* data,int x,int y,unsigned int dimX,unsigned int dimY){
    if (x < 0)
    {
        x = 0;
    }
    if (x >= dimX)
    {
        x = dimX - 1;
    }
    if (y < 0)
    {
        y = 0;
    }
    if (y >= dimY)
    {
        y = dimY - 1;
    }
    return data + (x * dimY + y);
}


float gd(float a, float b, float x)
{
    float c = (x-b) / a;
    return exp((-.5) * c * c) / (a * sqrt(2 * PI));
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

void print_fvec(FVec v)
{
    unsigned int i;
    printf("\n");
    for (i = 0; i < v.length; i++)
    {
        printf("%f ", v.data[i]);
    }
    printf("\n");
}

Image img_sc(Image a)
{
    Image b = a;
    b.data = malloc(b.dimX * b.dimY * b.numChannels * sizeof(float));
    return b;
}

Image gb_h(Image a, FVec gv)
{
    float **data_a=(float**)malloc(a.numChannels*sizeof(float*));
  //  float **data_b=(float**)malloc(a.numChannels*sizeof(float*));//NCHW
    float **data_b_T=(float**)malloc(a.numChannels*sizeof(float*));
    Image b = img_sc(a);

    int ext = gv.length / 2;
    int offset;
    unsigned int x, y, channel;
    float *pc;
    float sum;
    int i;

    #pragma omp parallel for
    for(i=0;i<a.numChannels;i++){
        data_a[i]=(float*)malloc(a.dimX*a.dimY*sizeof(float));
        data_b_T[i]=(float*)malloc(a.dimX*a.dimY*sizeof(float));
       // data_b[i]=(float*)malloc(b.dimX*b.dimY*sizeof(float));
        for(int n=0;n<a.dimX*a.dimY;n++){
            data_a[i][n]=a.data[i+a.numChannels*n];
        }
    }
    /*
    for(i=0;i<a.dimX*a.dimY*a.numChannels;i++){
        data_a[i%a.numChannels][i/a.numChannels]=a.data[i];
       // data_b[i%3][i/3]=0;
    }
    */


    



    #pragma omp parallel for collapse(2)
        for(channel=0;channel<a.numChannels;channel++){
            for(y=0;y<a.dimY;y++){
                for(x=0;x<a.dimX;x++){
                    float* A = data_a[channel];
                    pc = get_pixel2_Transpose(data_b_T[channel],x,y,a.dimX,a.dimY);
                    unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
                    deta = fmin(deta, gv.min_deta);
                    sum = 0;
                    /* Original
                    for (i = deta; i < gv.length-deta; i++)
                    {
                        offset = i - ext;
                        sum+=gv.data[i]/gv.sum[ext-deta]*(float)(*get_pixel2(A,x+offset,y,a.dimX,a.dimY));
                    }
                    pc[channel] = sum;
                    */

                    // Unrolling
                    // for (i = deta; i < gv.length-deta-8; i+=8)
                    // {   
                    //     offset = i - ext;
                    //     sum+=gv.data[i]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+1]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+1,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+2]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+2,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+3]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+3,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+4]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+4,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+5]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+5,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+6]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+6,y,a.dimX,a.dimY));
                    //     sum+=gv.data[i+7]/gv.sum[ext-deta]*(float)*(get_pixel2(A,x+offset+7,y,a.dimX,a.dimY));  
                    // }
                    // for (;i < gv.length;i++){
                    //     offset = i - ext;
                    //     sum+=gv.data[i]/gv.sum[ext-deta]*(float)(*get_pixel2(A,x+offset,y,a.dimX,a.dimY));
                    // }
                    // pc[channel] = sum;


                    
                    
                    __m256 sum_vec = _mm256_setzero_ps();
                    for (i = deta; i < deta+(gv.length-2*deta)/8*8; i+=8){
                        offset=i-ext;
                        __m256 gv_data_vec = _mm256_loadu_ps(gv.data+i);
                        __m256 gv_sum_vec  = _mm256_set1_ps(gv.sum[ext-deta]);
                        __m256 data_a_vec;
                        
                        // {
                        //     data_a_vec[0]=(*get_pixel2(A,x+offset,y,a.dimX,a.dimY));
                        //     data_a_vec[1]=(*get_pixel2(A,x+offset+1,y,a.dimX,a.dimY));
                        //     data_a_vec[2]=(*get_pixel2(A,x+offset+2,y,a.dimX,a.dimY));
                        //     data_a_vec[3]=(*get_pixel2(A,x+offset+3,y,a.dimX,a.dimY));
                        //     data_a_vec[4]=(*get_pixel2(A,x+offset+4,y,a.dimX,a.dimY));
                        //     data_a_vec[5]=(*get_pixel2(A,x+offset+5,y,a.dimX,a.dimY));
                        //     data_a_vec[6]=(*get_pixel2(A,x+offset+6,y,a.dimX,a.dimY));
                        //     data_a_vec[7]=(*get_pixel2(A,x+offset+7,y,a.dimX,a.dimY));
                        // }


                        //printf("%f,%f,%f,%f\n",gv_sum_vec[0],gv_sum_vec[1],gv_sum_vec[2],gv_sum_vec[3]);
                        //sum_vec = _mm256_mul_ps(_mm256_div_ps(gv_data_vec,gv_sum_vec),data_a_vec);
                        sum_vec = _mm256_div_ps(gv_data_vec,gv_sum_vec);
                        *pc+=sum_vec[0]*(float)*get_pixel2(A,x+offset,y,a.dimX,a.dimY)
                                    +sum_vec[1]*(float)*get_pixel2(A,x+offset+1,y,a.dimX,a.dimY)
                                    +sum_vec[2]*(float)*get_pixel2(A,x+offset+2,y,a.dimX,a.dimY)
                                    +sum_vec[3]*(float)*get_pixel2(A,x+offset+3,y,a.dimX,a.dimY)
                                    +sum_vec[4]*(float)*get_pixel2(A,x+offset+4,y,a.dimX,a.dimY)
                                    +sum_vec[5]*(float)*get_pixel2(A,x+offset+5,y,a.dimX,a.dimY)
                                    +sum_vec[6]*(float)*get_pixel2(A,x+offset+6,y,a.dimX,a.dimY)
                                    +sum_vec[7]*(float)*get_pixel2(A,x+offset+7,y,a.dimX,a.dimY);
                    }
                    
                    
                    for (;i < gv.length-deta;i++){
                        offset = i - ext;
                        *pc+=gv.data[i]/gv.sum[ext-deta]*(float)(*get_pixel2(A,x+offset,y,a.dimX,a.dimY));
                    }
                    
                }
                
            }
        }
    /*
    #pragma omp parallel for collapse(2)
        for(channel=0;channel<a.numChannels;channel++){
            for(y=0;y<a.dimY;y++){
                for(x=0;x<a.dimX;x++){
                    float* A = data_b_T[channel];
                    pc = get_pixel(b,y,x);
                    unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
                    deta = fmin(deta, gv.min_deta);
                    sum = 0;
                    
                    __m256 sum_vec = _mm256_setzero_ps();
                    for (i = deta; i < deta+(gv.length-2*deta)/8*8; i+=8){
                        offset=i-ext;
                        __m256 gv_data_vec = _mm256_loadu_ps(gv.data+i);
                        __m256 gv_sum_vec  = _mm256_set1_ps(gv.sum[ext-deta]);
                        __m256 data_a_vec;
                        
                        sum_vec = _mm256_div_ps(gv_data_vec,gv_sum_vec);
                        pc[channel]+=sum_vec[0]*(float)*get_pixel2(A,x+offset,y,a.dimX,a.dimY)
                                    +sum_vec[1]*(float)*get_pixel2(A,x+offset+1,y,a.dimX,a.dimY)
                                    +sum_vec[2]*(float)*get_pixel2(A,x+offset+2,y,a.dimX,a.dimY)
                                    +sum_vec[3]*(float)*get_pixel2(A,x+offset+3,y,a.dimX,a.dimY)
                                    +sum_vec[4]*(float)*get_pixel2(A,x+offset+4,y,a.dimX,a.dimY)
                                    +sum_vec[5]*(float)*get_pixel2(A,x+offset+5,y,a.dimX,a.dimY)
                                    +sum_vec[6]*(float)*get_pixel2(A,x+offset+6,y,a.dimX,a.dimY)
                                    +sum_vec[7]*(float)*get_pixel2(A,x+offset+7,y,a.dimX,a.dimY);
                    }
                     
                    for (;i < gv.length-deta;i++){
                        offset = i - ext;
                        pc[channel]+=gv.data[i]/gv.sum[ext-deta]*(float)(*get_pixel2(A,x+offset,y,a.dimX,a.dimY));
                    }
                    
                }
                
            }
        }
    */


    /*
    for(i=0;i<a.numChannels;i++){
       // data_b[i]=(float*)malloc(b.dimX*b.dimY*sizeof(float));
        for(int n=0;n<a.dimX*a.dimY;n++){
            
            b.data[i+a.numChannels*n]=data_b_T[i][n];
        }
    }
    */


    # pragma omp parallel for collapse(3)
    for (y = 0; y < a.dimY; y++)
    {
        for (x = 0; x < a.dimX; x++)
        {
            for (channel = 0; channel < a.numChannels; channel++)
            {
                pc = get_pixel(b, x, y);
                unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
                deta = fmin(deta, gv.min_deta);
                sum = 0;
                for (i = deta; i < gv.length-deta; i++)
                {
                    offset = i - ext;
                    sum += gv.data[i]/gv.sum[ext - deta] * (float)get_pixel2(data_b_T[channel], y + offset, x, a.dimY, a.dimX)[channel];
                }
                pc[channel] = sum;
            }
        }
    }


    return b;
}

Image gb_v(Image a, FVec gv)
{
    Image b = img_sc(a);

    int ext = gv.length / 2;
    int offset;
    unsigned int x, y, channel;
    float* pc;
    float sum;
    int i;
    # pragma omp parallel for private (x,y,channel,pc,sum,offset,i)
    for (x = 0; x < a.dimX; x++)
    {
        for (y = 0; y < a.dimY; y++)
        {
            for (channel = 0; channel < a.numChannels; channel++)
            {
                pc = get_pixel(b, x, y);
                unsigned int deta = fmin(fmin(a.dimY-y-1, y),fmin(a.dimX-x-1, x));
                deta = fmin(deta, gv.min_deta);
                sum = 0;
                for (i = deta; i < gv.length-deta; i++)
                {
                    offset = i - ext;
                    sum += gv.data[i] /gv.sum[ext - deta] * (float)get_pixel(a, x, y + offset)[channel];
                }
                pc[channel] = sum;
            }
        }
    }
    return b;
}

Image apply_gb(Image a, FVec gv)
{
    Image b = gb_h(a, gv);
    //Image c = gb_v(b, gv);
    //free(b.data);
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