/*  Project 1.1: RISC-V instructions to RISC-V compressed instructions in C89.
    The following is the starter code provided for you. To finish the task, you 
    should define and implement your own functions in translator.c, compression.c, 
    utils.c and their header files.
    Please read the problem description before you start.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "src/compression.h"
#include "src/utils.h"

#include "translator.h"

/*check if file can be correctly opened */
static int open_files(FILE** input, FILE** output, const char* input_name, const char* output_name){ 
    *input = fopen(input_name, "r");
    if (!*input){ /* open input file failed */
        printf("Error: unable to open input file: %s\n", input_name);
        return -1;
    }

    *output = fopen(output_name, "w");
    if (!*output){ /* open output file failed */
        printf("Error: unable to open output file: %s\n", output_name);
        fclose(*input);
        return -1;
    }
    return 0; /* no problem opening files */
}

static int close_files(FILE** input, FILE** output){
    fclose(*input);
    fclose(*output); /* close the files at the end */
    return 0;
}

static void print_usage_and_exit() {
    printf("Usage:\n");
    printf("Run program with translator <input file> <output file>\n"); /* print the correct usage of the program */
    exit(0);
}


/*Run the translator 
*/
int translate(const char*in, const char*out){
    FILE *input, *output;/* variable for input */
    int err = 0;
    char machine_code_32_input[32] = {0}; /* Temporary variable for input */
    unsigned int machine_code_32_int = 0;
    unsigned int *machine_code_32_int_array;
    int line_of_code = 0; /* number of lines of codes in input.s */
    int *classfication_array;

    /* other temporary variables */
    char c;
    int i = 0;


    if (in){    /* correct input file name */
        if(open_files(&input, &output, in, out) != 0)
            exit(1);

        /* Count lines of codes */
        while(!feof(input)){
            c = fgetc(input);
            if(c == '\n'){/* Count lines of codes */
                line_of_code++;
            }
        }/* Count lines of codes */
        
        fseek(input,-1,SEEK_END);
        c = fgetc(input);/* Count lines of codes */
        if(c!='\n'){
            line_of_code++;
        }/* Count lines of codes */

        rewind(input);/* re-scanf the input file */

        machine_code_32_int_array = (unsigned int*)malloc(line_of_code*sizeof(unsigned int));

        /* read the file line by line */
        while(!feof(input)){
            fscanf(input,"%s",machine_code_32_input);

            /* transform string into int */
            machine_code_32_int = strtol(machine_code_32_input,NULL,2);
            machine_code_32_int_array[i]=machine_code_32_int;
            i++;
        }

        /* now we have the number of lines and the (int)array saving those codes in binary form*/
        /* do the classfication, match each line of code to types defined in utils.h*/
        
        classfication_array = (int*)malloc(line_of_code * sizeof(int));

        /* maybe we can define functions in utils.c/.h*/
        /* the main function is classification_, which calls other judging functions */ 
        for (i = 0; i < line_of_code; i++){
            classfication_array[i] = machine_code_recognition(machine_code_32_int_array[i]);
        }   

        /* do the compression, update the offset according to the classification_array */

        /* ... */
        /* write correct result to the output file */

        for (i = 0; i < line_of_code; i++){
            int n;
            unsigned int code_output;
            code_output=code_compression(machine_code_32_int_array[i],/* do classification */
            classfication_array[i],i,classfication_array);
            /* do classification */
            if(classfication_array[i]==not_in_RVC||classfication_array[i]==sbtype_not_compressed||classfication_array[i]==ujtype_not_compressed){
                for(n=31;n>=0;n--){/* do classification */
                    fprintf(output,"%d",1&(code_output>>n));
                }    /* do classification */
            }/* ith code cannot be compressed, thus printing 32 b code. */
            else{/* write */
                for(n=15;n>=0;n--){/* write */
                    fprintf(output,"%d",1&(code_output>>n));
                }/* write */
            }/* ith code can be compressed, thus printing 16b code. */
            fprintf(output,"\n");
        }

        
        close_files(&input, &output);

        
        free(classfication_array); /* free allocated memory */
        free(machine_code_32_int_array); /* free allocated memory */

    }
    return err;
}

/* main func */
int main(int argc, char **argv){
    char* input_fname, *output_fname;
    int err;

    if (argc != 3) /* need correct arguments */
        print_usage_and_exit();

    input_fname = argv[1];
    output_fname = argv[2];

    err = translate(input_fname, output_fname); /* main translation process */
    if (err)
        printf("One or more errors encountered during translation operation.\n"); /* something wrong */
    else
        printf("Translation process completed successfully.\n"); /* correctly output */

    return 0;
}
