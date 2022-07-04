#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

/* Your code here... */

 
int machine_code_recognition(unsigned int origin){ /* used to classify the code */


    /* c.add */
    if(is_add(origin)){
        return add;
    }
       /* c.mv */
    else if(is_mv(origin)){
        return mv;
    } /* c.jr */
    else if(is_jr(origin)){
        return jr;
    } /* c.mv */
    else if(is_jalr(origin)){
        return jalr;
    } /* branch */
    else if(is_li(origin)){
        return li;
    } /* branch */
    else if(is_lui(origin)){
        return lui;
    } /* branch */
    else if(is_addi(origin)){
        return addi;
    } /* branch */
    else if(is_slli(origin)){
        return slli;
    } /* branch */
    else if(is_lw(origin)){
        return lw;
    } /* branch */
    else if(is_sw(origin)){
        return sw;
    } /* branch */
    else if(is_and(origin)){
        return and;
    } /* branch */
    else if(is_or(origin)){
        return or;
    } /* branch */
    else if(is_xor(origin)){
        return xor;
    } /* branch */
    else if(is_sub(origin)){
        return sub;
    } /* branch */
    else if(is_beqz(origin)){
        return beqz;
    } /* branch */
    else if(is_bnez(origin)){
        return bnez;
    } /* branch */
    else if(is_srli(origin)){
        return srli;
    } /* branch */
    else if(is_srai(origin)){
        return srai;
    } /* branch */
    else if(is_andi(origin)){
        return andi;
    } /* branch */
    else if(is_j(origin)){
        return j;
    } /* branch */
    else if(is_jal(origin)){
        return jal;
    } /* branch */

    /* branch */
    if(is_sbtype(origin)){return sbtype_not_compressed;}
    /* branch */
    else if(is_ujtype(origin)){return ujtype_not_compressed;}

    /* not_in_RVC */
    return not_in_RVC;
}

int is_add(unsigned int code){  /* written by caotx */
    unsigned int opcode,funct7,funct3,rd,rs1,rs2; /* declare parts of the code */

    /* use grab to select parts of the code */
    opcode = grab(code,6,0);
    rd = grab(code,11,7);
    funct3 = grab(code,14,12);
    /* use grab to select parts of the code */
    rs1 = grab(code,19,15);
    rs2 = grab(code,24,20);
    funct7 = grab(code,31,25);

    if(opcode != 51 /* 0b0110011 */){    
        return not_in_RVC;
    }/* judge the opcode */
    if(rd == 0 || rs2 == 0){return not_in_RVC;}
    if(rd != rs1){return not_in_RVC;}
    /* judge funct3  */
    if(funct3 != 0 /* 0b000 */){return not_in_RVC;}
    if(funct7 != 0 /* 0b0000000 */){    /* judge funct7 */
        return not_in_RVC;}
    return add;/* is an add */
}

int is_mv(unsigned int code){   /* written by caotx */
    unsigned int opcode,funct7,funct3,rd,rs1,rs2; /* declare parts of the code */

    /* use grab to select parts of the code */
    opcode = grab(code,6,0);
    rd = grab(code,11,7);
    funct3 = grab(code,14,12);
    rs1 = grab(code,19,15);/* use grab to select parts of the code */
    rs2 = grab(code,24,20);
    funct7 = grab(code,31,25);

    if(opcode != 51 /* 0b0110011 */){return not_in_RVC;}/* check opcode */
    if(rd == 0 || rs2 == 0){return not_in_RVC;}
    if(rs1 != 0){return not_in_RVC;}/* check funct3  */
    if(funct3 != 0 /* 0b000 */){return not_in_RVC;}
    if(funct7 != 0 /* 0b0000000 */){return not_in_RVC;} /* check funct7 */
    return mv;
}

int is_jr(unsigned int code){   /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1; /* declare parts of the code */
    int offset;
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/* use grab to select parts of the code */
    funct3 = grab(code, 14, 12);/* use grab to select parts of the code */
    rs1 = grab(code, 19, 15);
    offset = grab(code, 31, 20);/* use grab to select parts of the code */
    /* opcode should be 0b1100111 */
    if(opcode != 103 /* 0b1100111 */){return not_in_RVC;}
    /* rs1 !=0 and rd == 0 */
    if(rs1 == 0){return not_in_RVC;}
    if(rd != 0){return not_in_RVC;}
    /* if(!is_valid_register(rs1)){return not_in_RVC;} */
    /* funct3 should be 0b000 */
    if(funct3 != 0 /* 0b000 */){return not_in_RVC;}
    /* offset should be zero */
    if(offset != 0 /* 0b000000000000 */){return not_in_RVC;}
    /* return jr if it satisfy all the conditions */
    return jr;  
}

int is_jalr(unsigned int code){ /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1; /* declare parts of the code */
    int offset;
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0); /*grab*/
    rd = grab(code, 11, 7);/*grab*/
    funct3 = grab(code, 14, 12);/*grab*/
    rs1 = grab(code, 19, 15);/*grab*/
    offset = grab(code, 31, 20);/*grab*/
    /* judge the opcode */
    if(opcode != 103 /* 0b1100111 */){return not_in_RVC;}
    /* judge the registers */
    if(rs1 == 0){return not_in_RVC;}
    if(rd != 1){return not_in_RVC;}
    /* if(!is_valid_register(rs1)){return not_in_RVC;} */
    /* judge the funct3 */
    if(funct3 != 0 /* 0b000 */){return not_in_RVC;}
    /*judge the offset */
    if(offset != 0 /* 0b000000000000 */){return not_in_RVC;}
    /* return jalr if it satisfy all the conditions */
    return jalr;  
}

int is_li(unsigned int code){   /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1; /* declare parts of the code */
    int imm;
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);/*grab*/
    rd = grab(code, 11, 7);/*grab*/
    funct3 = grab(code, 14, 12);/*grab*/
    rs1 = grab(code, 19, 15);/*grab*/
    imm = grab(code, 31, 20);/*grab*/
    /* judge the opcode */
    if(opcode != 19 /* 0b0010011 */){return not_in_RVC;}
    /* judge the registers */
    if(rs1 != 0){return not_in_RVC;}
    if(rd == 0){return not_in_RVC;}
    /* if(!is_valid_register(rd)){return not_in_RVC;} */
    /* judge the funct3 */
    if(funct3 != 0 /* 0b000 */){return not_in_RVC;}
    imm = get_signed_int(imm, 11);
    /* imm should be zero except for imm[5:0] */
    if(imm < -32 /* 0b1000000 */){return not_in_RVC;}
    if(imm > 31){return not_in_RVC;}
    /* return li if it satisfy all the conditions */
    return li;  
}

int is_lui(unsigned int code){  /* written by zhuhl2 */
    int imm_20;
    unsigned int opcode, rd; /* declare parts of the code */
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/*grab*/
    imm_20 = grab(code, 31, 12);
    /* judge the opcode */
    if(opcode != 55 /* 0b0110111 */){return not_in_RVC;}
    /* judge the register rd */
    if(rd == 0){return not_in_RVC;}
    if(rd == 2){return not_in_RVC;}
    /* if(!is_valid_register(rd)){return not_in_RVC;} */
    /* imm_20 should be nonzero and limited to imm_20[5:0] */
    if(imm_20 == 0){return not_in_RVC;}/*imm*/
    imm_20 = get_signed_int(imm_20, 19);

    if(imm_20 < -32){return not_in_RVC;}
    if(imm_20 > 31){return not_in_RVC;}
    /* return lui if it satisfy all the conditions */
    return lui;  
}

int is_addi(unsigned int code){ /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1; /* declare parts of the code */
    int imm;
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);/*grab*/
    rd = grab(code, 11, 7);/*grab*/
    funct3 = grab(code, 14, 12);/*grab*/
    rs1 = grab(code, 19, 15);/*grab*/
    imm = grab(code, 31, 20);/*grab*/
    /* judge the opcode */
    if(opcode != 19 /* 0b0010011 */){return not_in_RVC;}
    /* judge the registers */
    if(rs1 != rd){return not_in_RVC;}
    if(rd == 0){return not_in_RVC;}
    /* if(!is_valid_register(rd)){return not_in_RVC;} */
    /* judge the funct3 */
    if(funct3 !=  0 /* 0b000 */){return not_in_RVC;}
    /* imm should be nonzero and limited to the first 6 bits */
    imm = get_signed_int(imm, 11);
    if(imm == 0 /* 0b000000000000 */){return not_in_RVC;}
    if(imm < -32){return not_in_RVC;}
    if(imm > 31){return not_in_RVC;}
    /* return addi if it satisfy all the conditions */
    return addi;  
}

int is_slli(unsigned int code){ /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1, highest, shamt; /* declare parts of the code */
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);/*grab*/
    rd = grab(code, 11, 7);
    funct3 = grab(code, 14, 12); /*grab*/
    rs1 = grab(code, 19, 15); /*rs1 = grab(code, 19, 15); */
    shamt = grab(code, 24, 20);
    highest = grab(code, 31, 25);
    /* judge the opcode */
    if(opcode != 19 /* 0b0010011 */){return not_in_RVC;}
    /* judge the registers */
    if(rs1 != rd){return not_in_RVC;}
    if(rd == 0){return not_in_RVC;}
    /* if(!is_valid_register(rs1)){return not_in_RVC;} */
    /* judge the funct3 */
    if(funct3 != 1 /* 0b001 */){return not_in_RVC;}
    /* shamt should be limited to shamt[4:0] */
    if(shamt > 31){return not_in_RVC;}
    /* the highest should be 0 */
    if(highest != 0 /* 0b0000000 */){return not_in_RVC;}
    /* return slli if it satisfy all the conditions */
    return slli;  
}

int is_lw(unsigned int code){   /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1; /* declare parts of the code */
    int offset;
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/*grab*/
    funct3 = grab(code, 14, 12);
    /* use grab to select parts of the code */
    rs1 = grab(code, 19, 15);
    /* use grab to select parts of the code */
    offset = grab(code, 31, 20);
    /* judge the opcode */
    if(opcode != 3 /* 0b0000011 */){return not_in_RVC;}
    /* judge the registers */
    if(!is_valid_register(rd)){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;}
    /* judge the funct3 */
    if(funct3 != 2 /* 0b010 */){return not_in_RVC;}
    /* the last two bits of offset should be zero and the offset should be limited to offset[6:0] */
    offset = get_signed_int(offset, 11);
    if((offset & 3 /*0b000000000011*/) != 0 /*0b000000000000*/ ){return not_in_RVC;}
    if(offset < 0){return not_in_RVC;}
    if(offset > get_signed_int(127,7)){return not_in_RVC;}
    /* return lw if it satisfy all the conditions */
    return lw; 
}

int is_sw(unsigned int code){   /* written by caotx */
    unsigned int funct3,opcode, imm_5, imm_7;/* variable declaration */
    int offset;
    unsigned int rs1,rs2;
    /* use grab to select parts of the code */
    funct3=grab(code,14,12);
    imm_5 = grab(code, 11, 7);
    /* use grab to select parts of the code */
    imm_7 = grab(code, 31, 25);
    opcode=grab(code,6,0);
    /* use grab to select parts of the code */
    rs2=grab(code,24,20);
    rs1=grab(code,19,15);
    /* Check opcode */
    if(opcode!=35/*0b0100011*/){return not_in_RVC;}
    /* Check funct3 */
    if(funct3!=2/*0b010*/){return not_in_RVC;}
    /* register should be from x8 to x15 */
    if(!is_valid_register(rs1)){return not_in_RVC;}
    if(!is_valid_register(rs2)){return not_in_RVC;}
    /* the last two bits of offset should be zero and the offset should be limited to offset[6:0] */
    offset = (imm_7 << 5) + imm_5;
    offset = get_signed_int(offset, 11);
    if((offset & 3 /*0b000000000011*/) != 0 /*0b000000000000*/ ){return not_in_RVC;}
    if(offset < 0){return not_in_RVC;}
    if(offset > get_signed_int(127,7)){return not_in_RVC;}
    /* return sw if it satisfy all the conditions */
    return sw;
}

int is_and(unsigned int code){  /* written by caotx */
    unsigned int opcode,funct3,rd;/* use grab to select parts of the code */
    unsigned int funct7,rs1,rs2;
    funct3=grab(code,14,12);/* use grab to select parts of the code */
    opcode=grab(code,6,0);
    rs1 = grab(code,19,15);/* use grab to select parts of the code */
    rs2 = grab(code,24,20);
    funct7 = grab(code,31,25);/*funct7*/
    rd = grab(code,11,7);
    if(opcode!=51/*0b0110011*/){return not_in_RVC;}
    if(funct7!=0/*0b0000000*/){return not_in_RVC;}
    if(funct3!=7/*0b111*/){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;}/* exceeded */
    if(!is_valid_register(rs2)){return not_in_RVC;}/* exceeded */
    if(!is_valid_register(rd)){return not_in_RVC;}/* exceeded */
    if(rs1!=rd){return not_in_RVC;}/*invalid*/
    return and;
}

int is_or(unsigned int code){   /* written by caotx */
    unsigned int opcode,funct3,rd;/* use grab to select parts of the code */
    unsigned int funct7,rs1,rs2;
    funct3=grab(code,14,12);
    opcode=grab(code,6,0);/* use grab to select parts of the code */
    rs1 = grab(code,19,15);
    rs2 = grab(code,24,20);
    funct7 = grab(code,31,25);/* use grab to select parts of the code */
    rd = grab(code,11,7);
    if(opcode!=51/*0b0110011*/){return not_in_RVC;}
    if(funct7!=0/*0b0000000*/){return not_in_RVC;}
    if(funct3!=6/*0b110*/){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;}/*invalid*/
    if(!is_valid_register(rs2)){return not_in_RVC;}/*invalid*/
    if(!is_valid_register(rd)){return not_in_RVC;}/*invalid*/
    if(rs1!=rd){return not_in_RVC;}/*invalid*/
    return or;
}

int is_xor(unsigned int code){  /* written by caotx */
    unsigned int opcode,funct3,rd;/* use grab to select parts of the code */
    unsigned int funct7,rs1,rs2;
    funct3=grab(code,14,12);/* use grab to select parts of the code */
    opcode=grab(code,6,0);
    rs1 = grab(code,19,15);/* use grab to select parts of the code */
    rs2 = grab(code,24,20);
    funct7 = grab(code,31,25);/* use grab to select parts of the code */
    rd = grab(code,11,7);
    if(opcode!=51/*0b0110011*/){return not_in_RVC;}
    if(funct7!=0/*0b0000000*/){return not_in_RVC;}
    if(funct3!=4/*0b100*/){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;}
    if(!is_valid_register(rs2)){return not_in_RVC;}/*invalid*/
    if(!is_valid_register(rd)){return not_in_RVC;}
    if(rs1!=rd){return not_in_RVC;}/*invalid*/
    return xor;
}

int is_sub(unsigned int code){  /* written by caotx */
    unsigned int opcode,funct3,rd;
    unsigned int funct7,rs1,rs2;/* use grab to select parts of the code */
    funct3=grab(code,14,12);
    opcode=grab(code,6,0);/* use grab to select parts of the code */
    rs1 = grab(code,19,15);
    rs2 = grab(code,24,20);/* use grab to select parts of the code */
    funct7 = grab(code,31,25);
    rd = grab(code,11,7);
    if(opcode!=51/*0b0110011*/){return not_in_RVC;}/*invalid*/
    if(funct7!=32/*0b0100000*/){return not_in_RVC;}
    if(funct3!=0/*0b000*/){return not_in_RVC;}/*invalid*/
    if(!is_valid_register(rs1)){return not_in_RVC;}
    if(!is_valid_register(rs2)){return not_in_RVC;}/*invalid*/
    if(!is_valid_register(rd)){return not_in_RVC;}
    if(rs1!=rd){return not_in_RVC;}/*invalid*/
    return sub;
}

int is_beqz(unsigned int code){ /* written by caotx */
    /* declaration */
    unsigned int opcode,funct3;
    unsigned int rs1,rs2;
    int offset;
    /* use grab to select parts of the code */
    funct3=grab(code,14,12);
    opcode=grab(code,6,0);/* use grab to select parts of the code */
    rs1 = grab(code,19,15);/* use grab to select parts of the code */
    rs2 = grab(code,24,20);/* obtain offset from 32b code */
    offset= (int)((grab(code,11,8)<<1)+ (grab(code,30,25)<<5) + (grab(code,7,7)<<11) + (grab(code,31,31)<<12));
    offset=get_signed_int(offset,12);
    

    if(funct3!=0){return not_in_RVC;}/*invalid*/
    if(opcode!=99/*0b1100011*/){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;}/*invalid*/
    if(rs2!=0){return not_in_RVC;}
   /* if(offset<-256)*/
    /*{return not_in_RVC;}*/
   /* if(offset>255)*/
  /*  {return not_in_RVC;} */


    return beqz;
}

int is_bnez(unsigned int code){ /* written by caotx */
    unsigned int opcode,funct3;/* declaration */
    unsigned int rs1,rs2;
    int offset;/* declaration */
    funct3=grab(code,14,12);/* use grab to select parts of the code */
    opcode=grab(code,6,0);
    rs1 = grab(code,19,15);/* use grab to select parts of the code */
    rs2 = grab(code,24,20); /* obtain offset from 32b code */
    offset= (int)((grab(code,11,8)<<1)+ (grab(code,30,25)<<5) + (grab(code,7,7)<<11) + (grab(code,31,31)<<12));
    offset=get_signed_int(offset,12);
    

    if(funct3!=1/*0b001*/){return not_in_RVC;}
    if(opcode!=99/*0b1100011*/){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;}
    if(rs2!=0){return not_in_RVC;}/*invalid*/
   /* if(offset<-256)*/
    /*{return not_in_RVC;}*/
   /* if(offset>255)*/
  /*  {return not_in_RVC;} */

    return bnez;
}

int is_srli(unsigned int code){ /* written by zhuhl2 */
    unsigned int opcode, rd, funct3, rs1, highest, shamt; /* declare parts of the code */
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/* use grab to select parts of the code */
    funct3 = grab(code, 14, 12);
    rs1 = grab(code, 19, 15);/* use grab to select parts of the code */
    shamt = grab(code, 24, 20);
    highest = grab(code, 31, 25);/* use grab to select parts of the code */
    /* judge the opcode */
    if(opcode != 19/*0b0010011*/){return not_in_RVC;}
    /* judge the registers */
    if(rs1 != rd){return not_in_RVC;}
    if(rd == 0){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;} 
    /* judge the funct3 */
    if(funct3 != 5/*0b101*/){return not_in_RVC;}
    /* shamt should be limited to shamt[4:0] */
    if(shamt > 31){return not_in_RVC;}
    /*judge the highest */
    if(highest != 0/*0b0000000*/){return not_in_RVC;}
    /* return srli if it satisfy all the conditions */
    return srli;   
}

int is_srai(unsigned int code){ /* written by zhuhl2 */ /* this funct hasn't been tested */
    unsigned int opcode, rd, funct3, rs1, highest,shamt; /* declare parts of the code */
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/* use grab to select parts of the code */
    funct3 = grab(code, 14, 12);
    rs1 = grab(code, 19, 15);/* use grab to select parts of the code */
    shamt = grab(code, 24, 20);
    highest = grab(code, 31, 25);/* use grab to select parts of the code */
    /* opcode should be 0b0010011 */
    if(opcode != 19/*0b0010011*/){return not_in_RVC;}
    /* rd == rs1, and should be valid */
    if(rs1 != rd){return not_in_RVC;}
    if(rd == 0){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;} 
    /* funct3 should be 0b101 */
    if(funct3 != 5/*0b101*/){return not_in_RVC;}
    /* shamt should be limited to shamt[4:0] */
    if(shamt > 31) {return not_in_RVC;}
    /*highest should be 0b0100000 */
    if(highest != 32/*0b0100000*/){return not_in_RVC;}
    /* return srai if it satisfy all the conditions */
    return srai;
}

int is_andi(unsigned int code){ /* written by zhuhl2 */  /* this funct hasn't been tested */
    unsigned int opcode, rd, funct3, rs1; /* declare parts of the code */
    int imm;
    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/* use grab to select parts of the code */
    funct3 = grab(code, 14, 12);
    rs1 = grab(code, 19, 15);/* use grab to select parts of the code */
    imm = grab(code, 31, 20);
    /* opcode should be 0b0010011 */
    if(opcode != 19/*0b0010011*/){return not_in_RVC;}
    /* rd == rs1, and should be valid */
    if(rs1 != rd){return not_in_RVC;}
    if(rd == 0){return not_in_RVC;}
    if(!is_valid_register(rs1)){return not_in_RVC;} 
    /* funct3 should be 0b111 */
    if(funct3 != 7/*0b111*/){return not_in_RVC;}
    imm = get_signed_int(imm, 11);
    /* imm should be zero except for imm[5:0] */
    if(imm < -32 /* 0b1000000 */){return not_in_RVC;}
    if(imm > 31){return not_in_RVC;}
    /* return andi if it satisfy all the conditions */
    return andi;   
}

int is_j(unsigned int code){    /* written by caotx */
    int offset;
    unsigned int rd, opcode;
/* use grab to select parts of the code */

    offset= (int)((grab(code,30,21)<<1)+(grab(code,20,20)<<11)+(grab(code,19,12)<<12)+(grab(code,31,31)<<20));
    offset=get_signed_int(offset,20);

/* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);/* use grab to select parts of the code */
    if(opcode!=111/*0b1101111*/){return not_in_RVC;}
    if(rd!=0){return not_in_RVC;}
 /*   if(offset<get_signed_int(2048,11))*/
 /*   {return not_in_RVC;}*/
 /*   if(offset>get_signed_int(2047,11))*/
  /*  {return not_in_RVC;}*/


    return j;
}

int is_jal(unsigned int code){  /* written by caotx */
    int offset;
    unsigned int rd, opcode;

/* use grab to select parts of the code */
    offset= (int)((grab(code,30,21)<<1)+(grab(code,20,20)<<11)+(grab(code,19,12)<<12)+(grab(code,31,31)<<20));
    offset=get_signed_int(offset,20);

    /* use grab to select parts of the code */
    opcode = grab(code, 6, 0);
    rd = grab(code, 11, 7);
    if(opcode!=111/*0b1101111*/){return not_in_RVC;}/*invalid*/
    if(rd!=1){return not_in_RVC;}/*invalid*/
 /*   if(offset<get_signed_int(2048,11))*/
 /*   {return not_in_RVC;}*/
 /*   if(offset>get_signed_int(2047,11))*/
  /*  {return not_in_RVC;}*/


    return jal;
}

int is_sbtype(unsigned int code){
    unsigned int opcode;
    opcode = grab(code, 6, 0);/* use grab to select parts of the code */
    if(opcode!=99){return not_in_RVC;}/*invalid*/
    return sbtype_not_compressed;
}
int is_ujtype(unsigned int code){
    unsigned int opcode;
    opcode = grab(code, 6, 0);/* use grab to select parts of the code */
    if(opcode!=111){return not_in_RVC;}/*invalid*/
    return sbtype_not_compressed;
}

/*is_valid_register*/
int is_valid_register(unsigned int reg){
    return (8 <= reg) && (reg <= 15);
}
/*grab*/
unsigned int grab(unsigned int code, int l,int r){
    unsigned int result;
    result = (code >> r) & ((1 << (l - r + 1)) - 1);/*grab*/
    return result;/*grab*/
}
/*get_signed_int*/
int get_signed_int(int d,int highest_pos){
    if((d>>highest_pos)==1){/*get_signed_int*/
        return d|(((0xFFFFFFFF)>>(highest_pos+1))<<(highest_pos+1));
    }/*get_signed_int*/
    else return d;
}
/* reg_map */
unsigned int reg_map(unsigned int reg){
    return reg-8;
}
