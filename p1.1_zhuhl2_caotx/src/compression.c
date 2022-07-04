#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "compression.h"
#include "utils.h"

unsigned int code_compression(unsigned int code, int type, int count, int *c_array){

    if(type == not_in_RVC){ /* can't be compressed */ 
        return code;
    }
    else if(type == add){   /* it's add type */
        return CR_format(9, grab(code, 11, 7), grab(code, 24, 20), 2);
    }   
    else if(type == mv){ /*it's mv type*/
        return CR_format(8, grab(code, 11, 7), grab(code, 24, 20), 2);
    }
    else if(type == jr){    /*it's jr type*/
        return CR_format(8, grab(code, 19, 15), 0, 2);
    }
    else if(type == jalr){  /*it's jalr type*/
        return CR_format(9, grab(code, 19, 15), 0, 2);
    }
    else if(type == li){    /*it's li type*/
        return CI_format(2, grab(code, 25, 25), grab(code, 11, 7), grab(code, 24, 20), 1);
    }
    else if(type == lui){   /*it's lui type*/
        return CI_format(3, grab(code, 17, 17), grab(code, 11, 7), grab(code, 16, 12), 1);
    }
    else if(type == addi){  /*it's addi type*/
        return CI_format(0, grab(code, 25, 25), grab(code, 11, 7), grab(code, 24, 20), 1);
    }
    else if(type == slli){  /*it's slli type*/
        return CI_format(0, 0, grab(code, 11, 7), grab(code, 24, 20), 2);
    }
    else if(type == lw){    /*it's lw type*/
        return CL_format(2, grab(code, 25, 23), reg_map(grab(code, 19, 15)), (grab(code, 22, 22) * 2 + grab(code, 26, 26)), reg_map(grab(code, 11, 7)), 0);
    }
    else if(type == sw){    /*it's sw type*/
        return CS_format_1(6, (grab(code, 25, 25) << 2) + grab(code, 11, 10), reg_map(grab(code, 19, 15)), (grab(code, 9, 9) * 2 + grab(code, 26, 26)), reg_map(grab(code,24,20)), 0);
    }
    else if(type == and){   /*it's and type*/
        return CS_format_2(35, reg_map(grab(code,11,7)), 3, reg_map(grab(code,24,20)), 1);
    }
    else if(type == or){    /*it's or type*/
        return CS_format_2(35, reg_map(grab(code,11,7)), 2, reg_map(grab(code,24,20)), 1);
    }
    else if(type == xor){   /*it's xor type*/
        return CS_format_2(35, reg_map(grab(code,11,7)), 1, reg_map(grab(code,24,20)), 1);
    }
    else if(type == sub){   /*it's sub type*/
        return CS_format_2(35, reg_map(grab(code,11,7)), 0, reg_map(grab(code,24,20)), 1);
    }
    else if(type == beqz){  /*it's beqz type*/
        /* get the offset */
        int offset;
        int lines_jumped; /*save lines jumped*/
        offset= (int)((grab(code,11,8)<<1)+ /*offset's LSB is 2^0*/
        (grab(code,30,25)<<5) + 
        (grab(code,7,7)<<11) + 
        (grab(code,31,31)<<12));
        offset=get_signed_int(offset,12); /* indefinite whether offset is multiple of 4 */
        if (offset >= 0){
            lines_jumped = offset / 4;
            for (; lines_jumped > 0; lines_jumped --){
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset - 2); /*if lines passed can be compressed, subtract 2 from the offset*/
                count ++; /*adding count*/
            }/* offset upgraded */   
        }
        else{
            lines_jumped = - (offset / 4);
            for (; lines_jumped > 0; lines_jumped --){
                count --; /*subtracting count*/
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset + 2); /*if lines passed can be compressed, subtract 2 from the offset*/
            }/* offset upgraded */  
        }
        return CB_format_1(6, (grab(offset, 8, 8) << 2) + grab(offset, 4, 3), reg_map(grab(code,19,15)), (grab(offset, 7, 6) << 3) + (grab(offset, 2, 1) << 1) + (grab(offset, 5, 5)), 1); 
    }  
    else if(type == bnez){
        /* get the offset */
        int offset;
        int lines_jumped; /*save lines jumped*/
        offset= (int)((grab(code,11,8)<<1)+  /*offset's LSB is 2^0*/
        (grab(code,30,25)<<5) +  /*offset's LSB is 2^0*/
        (grab(code,7,7)<<11) +  /*offset's LSB is 2^0*/
        (grab(code,31,31)<<12)); /*offset's LSB is 2^0*/
        offset=get_signed_int(offset,12);/*offset's LSB is 2^0*/ 
        if (offset >= 0){ 
            lines_jumped = offset / 4; /*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset - 2); 
                count ++;/*adding count*/            
            }/* offset upgraded */   
        }
        else{
            lines_jumped = - (offset / 4);/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                count --;/*subtracting count*/
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset + 2); 
            }/* offset upgraded */  
        }
        return CB_format_1(7, (grab(offset, 8, 8) << 2) + grab(offset, 4, 3), reg_map(grab(code,19,15)), (grab(offset, 7, 6) << 3) + (grab(offset, 2, 1) << 1) + (grab(offset, 5, 5)), 1);
    }
    else if(type == srli){ /*it's srli type*/
        return CB_format_2(4, 0, 0, reg_map(grab(code, 11, 7)), grab(code, 24, 20), 1);
    }
    else if(type == srai){  /*it's srai type*/
        return CB_format_2(4, 0, 1, reg_map(grab(code, 11, 7)), grab(code, 24, 20), 1);
    }
    else if(type == andi){  /*it's andi type*/
        return CB_format_2(4, grab(code, 25, 25), 2, reg_map(grab(code, 11, 7)), grab(code, 24, 20), 1);
    }
    else if(type == j){  /*it's j type*/
        /* get the offset */
        int offset;
        int lines_jumped; /* save lines jumped*/
        offset= (int)((grab(code,30,21)<<1)+ /*offset's LSB is 2^0*/
        (grab(code,20,20)<<11)+
        (grab(code,19,12)<<12)+/*offset's LSB is 2^0*/
        (grab(code,31,31)<<20));
        offset=get_signed_int(offset,20);
        if (offset >= 0){
            lines_jumped = offset / 4;/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset - 2); 
                count ++; /*adding count*/
            }/* offset upgraded */   
        }
        else{
            lines_jumped = - (offset / 4);/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                count --;/*subtracting count*/
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset + 2); 
            }/* offset upgraded */  
        }
        return CJ_format(5, grab(offset,5,5) + (grab(offset, 3, 1) << 1) + (grab(offset, 7, 7) << 4) + (grab(offset, 6, 6) << 5) + (grab(offset, 10, 10) << 6) + (grab(offset, 9 , 8) << 7) + (grab(offset, 4, 4) << 9) + (grab(offset, 11, 11) << 10), 1);
    }
    else if(type == jal){ /* its jal type*/
        /* get the offset */
        int offset;
        int lines_jumped;  /* save lines jumped*/
        offset= (int)((grab(code,30,21)<<1)+ /*offset's LSB is 2^0*/
            (grab(code,20,20)<<11)+
            (grab(code,19,12)<<12)+/*offset's LSB is 2^0*/
            (grab(code,31,31)<<20));
        offset=get_signed_int(offset,20);
        if (offset >= 0){
            lines_jumped = offset / 4;/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset - 2); /*offset's LSB is 2^0*/
                count ++; /*adding count*/
            }/* offset upgraded */   
        }
        else{
            lines_jumped = - (offset / 4);/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                count --; /*subtracting count*/
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset + 2); /*offset's LSB is 2^0*/
            }/* offset upgraded */  
        }
        return CJ_format(1, grab(offset,5,5) + (grab(offset, 3, 1) << 1) + (grab(offset, 7, 7) << 4) + (grab(offset, 6, 6) << 5) + (grab(offset, 10, 10) << 6) + (grab(offset, 9 , 8) << 7) + (grab(offset, 4, 4) << 9) + (grab(offset, 11, 11) << 10), 1);
    }
    else if(type == sbtype_not_compressed){
        /* get the offset */
        int offset;
        int lines_jumped; /* save lines jumped*/
        offset= (int)((grab(code,11,8)<<1)+  /*offset's LSB is 2^0*/
        (grab(code,30,25)<<5) + 
        (grab(code,7,7)<<11) + /*offset's LSB is 2^0*/
        (grab(code,31,31)<<12));
        offset=get_signed_int(offset,12);
        if(offset>0){
            lines_jumped=offset/4;/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset - 2); /*offset's LSB is 2^0*/
                count ++;/*adding count*/
            }
        }
        else{
            lines_jumped = - (offset / 4);/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                count --;/*subtracting count*/
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset + 2); /*offset's LSB is 2^0*/
            }/* offset upgraded */  
        }
        code = code & (33550463);/* 0b00000001111111111111000001111111 */
        code += grab(offset,4,1) << 8;
        code += grab(offset,11,11) << 7; /*code's offset modified*/
        code += grab(offset,10,5) << 25;/*code's offset modified*/
        code += (unsigned int)(grab(offset,12,12) << 31);/*code's offset modified*/
        return code;
    }
    else if(type == ujtype_not_compressed){
        /* get the offset */
        int offset;
        int lines_jumped; /* save lines jumped*/
        offset= (int)((grab(code,30,21)<<1)+ /*offset's LSB is 2^0*/
            (grab(code,20,20)<<11)+/*code's offset modified*/
            (grab(code,19,12)<<12)+/*code's offset modified*/
            (grab(code,31,31)<<20));/*code's offset modified*/
        offset=get_signed_int(offset,20);/*code's offset modified*/
        if(offset>0){
            lines_jumped=offset/4; /*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset - 2); /*offset modified*/
                count ++; /*adding count*/
            }
        }
        else{
            lines_jumped = - (offset / 4);/*lines jumped*/
            for (; lines_jumped > 0; lines_jumped --){
                count --; /*subtracting count*/
                offset = (c_array[count] == 0 ||c_array[count] == ujtype_not_compressed || c_array[count] == sbtype_not_compressed)? (offset) : (offset + 2); /*offset modified*/
            }/* offset upgraded */  
        }
        code = code & 4095;/*code's offset modified*/
        code += grab(offset,19,12) << 12;
        code += grab(offset,11,11) << 20;/*code's offset modified*/
        code += grab(offset,10,1)  << 21;/*code's offset modified*/
        code += (unsigned int)(grab(offset,20,20) << 31);/*code's offset modified*/
        return code;
    }


    return code;
}

unsigned int CR_format(unsigned int funct4, unsigned int rd, unsigned int rs2, unsigned int opcode){ /* written by zhuhl2 */
    return (opcode << 0) + (rs2 << 2) + (rd << 7) + (funct4 << 12);
}
unsigned int CI_format(unsigned int funct3, unsigned int imm_1, unsigned int rd, unsigned int imm_5, unsigned int opcode){  /* written by zhuhl2 */
    return (opcode << 0) + (imm_5 << 2) + (rd << 7) + (imm_1 << 12) + (funct3 << 13); 
}
unsigned int CL_format(unsigned int funct3, unsigned int imm_3, unsigned int rs1, unsigned int imm_2, unsigned int rd, unsigned int opcode){ /* written by zhuhl2 */
    return (opcode << 0) + (rd << 2) + (imm_2 << 5) + (rs1 << 7) + (imm_3 << 10) + (funct3 << 13);
}
unsigned int CS_format_1(unsigned int funct3, unsigned int imm_3, unsigned int rs1, unsigned int imm_2, unsigned int rs2, unsigned int opcode){ /* written by zhuhl2 */
    return (opcode << 0) + (rs2 << 2) + (imm_2 << 5) + (rs1 << 7) + (imm_3 << 10) + (funct3 << 13);
}
unsigned int CS_format_2(unsigned int funct6, unsigned int rd, unsigned int funct2, unsigned int rs2, unsigned int opcode){ /* written by zhuhl2 */
    return (opcode << 0) + (rs2 << 2) + (funct2 << 5) + (rd << 7) + (funct6 << 10);
}
unsigned int CB_format_1(unsigned int funct3, unsigned int imm_3, unsigned int rd, unsigned int imm_5, unsigned int opcode){ /* written by zhuhl2 */
    return (opcode << 0) + (imm_5 << 2) + (rd << 7) + (imm_3 << 10) + (funct3 << 13);
}
unsigned int CB_format_2(unsigned int funct3, unsigned int imm_1, unsigned int funct2, unsigned int rd, unsigned int imm_5, unsigned int opcode){ /* written by zhuhl2 */
    return (opcode << 0) + (imm_5 << 2) + (rd << 7) + (funct2 << 10) + (imm_1 << 12) + (funct3 << 13);
}
unsigned int CJ_format(unsigned int funct3, unsigned int jump_target, unsigned int opcode){/* written by zhuhl2 */ 
    return (opcode << 0) + (jump_target << 2) + (funct3 << 13);
}

