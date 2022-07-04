float* transpose(float* src, unsigned int N, unsigned int M){   // N is for dimX of the src matrix, M is for dimY of the src matrix
    float* dst=(float*)malloc(3*N*M*sizeof(float)); // so M is for dimX of the dst matrix, N is for dimY of the dst matrix
    #pragma omp parallel for
    for(int n = 0; n<N*M; n++) {    // if n denotes the index in dst, then x for dst = n % M, y for dst = n / M
        int i = n/M;    // i means y for dst, x for src
        int j = n%M;    // j means x for dst, y for src
        dst[3*n] = src[3*N*j + 3*i];
        dst[3*n+1] = src[3*N*j + 3*i + 1];
        dst[3*n+2] = src[3*N*j + 3*i + 2];
    }
    free(src);
    return dst;
}
