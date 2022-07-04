void* Transpose(float* data, unsigned int X, unsigned int Y){
    float* temp_data = (float*)malloc((X * Y * 3) * sizeof(float));
    memcpy(temp_data, data, sizeof(float)* X * Y * 3);  // use temp_data to store the origin data, modify data directlY
    //  the x and y corresponds to the position in the old matrix
    int flag, n,  block_count, block_remaining, x_start, y_start, inner_x, inner_y, x, y, blocksize;
    if( X >= Y) {
        n = Y;
        flag = 0; // that means the n*n matrix has the same n with Y
    }
    else{
        n = X;
        flag = 1;   // that means the n*n matrix has the same n with X
    }
    blocksize = n / 5; // we may change it whatever
    block_count = n / blocksize; 
    block_remaining = n % blocksize;
    // # pragma omp parallel for
    for(int x = 0; x < block_count; x++){
        for(int y = 0; y < block_count; y++){
            x_start = x * blocksize;
            y_start = y * blocksize;
            for(inner_x = x_start; inner_x < x_start + blocksize; inner_x ++){
                for(inner_y = y_start; inner_y < y_start + blocksize; inner_y++){
                    for(int c = 0; c < 3; c++){ // c denotes the colour channels that is using
                        data[(inner_x + inner_y * X) * 3 + c] = temp_data[(inner_y + inner_x * X) * 3 + c];
                    } 
                }
            }
        } 
    }
    // transpose the remaining elements if they exist
    if(block_remaining){
        x_start = 0 + block_count * blocksize;
        y_start = 0;
        
        for (inner_x = x_start; inner_x < n; inner_x++){
            for(inner_y = 0; inner_y < n; inner_y++){
                data[inner_x + inner_y * X] = temp_data[inner_x * X + inner_y];
            }
        }

        x_start = 0;
        y_start = 0 + block_count * blocksize;

        for (inner_x = x_start; inner_x < n - block_remaining; inner_x++){
            for(inner_y = y_start; inner_y < n; inner_y++){
                data[inner_x + inner_y * X] = temp_data[inner_x * X + inner_y];
            }
        }
    }

    // transpose the elements outside of the n*n matrix
    if(flag == 0){
        x_start = n;
        y_start = 0;
        for(inner_x = x_start; inner_x < X; inner_x ++){
            for(inner_y = y_start; inner_y < Y; inner_y ++){
                data[inner_x + inner_y * X] = temp_data[inner_x * X + inner_y];
            }
        }
    }
    if(flag == 1){
        x_start = 0;
        y_start = n;
        for(inner_x = x_start; inner_x < X; inner_x ++){
            for(inner_y = y_start; inner_y < Y; inner_y ++){
                data[inner_x + inner_y * X] = temp_data[inner_x * X + inner_y];
            }
        }
    }

    free(temp_data);
    temp_data = NULL;
}