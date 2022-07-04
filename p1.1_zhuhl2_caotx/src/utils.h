#ifndef UTILS_H
#define UTILS_H

/* This macro aims to match the return value of function "machine_code_recognition" to those RVC types */

/* RVC instructions mapping to integer codings */
#define    add  1 /* add */
#define    mv   2 /* mv */
#define    jr   3 /* jr */
#define    jalr 4 /* jalr */
#define    li   5 /* li */
#define    lui  6 /* lui */
#define    addi 7 /* addi */
#define    slli 8 /* slli */
#define    lw   9 /* lw */
#define    sw   10 /* sw */
#define    and  11 /* and */
#define    or   12 /* or */
#define    xor  13 /* xor */
#define    sub  14 /* sub */
#define    beqz 15  /* offset required */
#define    bnez 16  /* offset required */
#define    srli 17 /* srli */
#define    srai 18 /* srai */
#define    andi 19 /* andi */
#define    j    20  /* offset required */
#define    jal  21  /* offset required */
#define not_in_RVC 0 /* making it more convenient to run if() */

#define sbtype_not_compressed 22
#define ujtype_not_compressed 23


/* This function reads the machine code, and return the corresponding type. */
int machine_code_recognition(unsigned int);

int is_add(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_mv(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_jr(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_jalr(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_li(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_lui(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_addi(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_slli(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_lw(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_sw(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_and(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_or(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_xor(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_sub(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_beqz(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_bnez(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_srli(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_srai(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_andi(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_j(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */
int is_jal(unsigned int);/* declaration of a function that determines whether a line of code is an compressable instruction. */

int is_sbtype(unsigned int);/* declaration of a function that determines whether a line of code belongs to sbtype instruction. */
int is_ujtype(unsigned int);/* declaration of a function that determines whether a line of code belongs to ujtype instruction. */

int is_valid_register(unsigned int);
/*  check whether the register referred in the code in valid(x8 - x15) */

unsigned int grab(unsigned int,int,int); /* takes the code, the left position, and the right position */
/* return the corresponding part of the code */
int get_signed_int(int d,int h); /* takes in an unsigned integer and its highest position.
return a interger that is stored as int_32 but is actually an int_h */
unsigned int reg_map(unsigned int);
/* map x8-x15 to 3b integer */
#endif
