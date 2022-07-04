#   Project 1.2 RVC instructions to RISC-V instructions in RISC-V

#   You can suppose you will get the number of machine codes you should translate from 
#   "lines_of_rvc_codes" and get the 16 or 32 bits-machine codes from "rvc_codes". 
#   Both of them are stored in file "input.S"

#   When we test your project, we will change the input files content.

# DO NOT REMOVE THIS LINE WHEN SUBMITTING
.import input.S

# Your code to convert the RVC codes to RVI codes.

#                                               OP   FUNCT3 FUNCT7
# 1	    c.add	CR	add rd, rd, rs2	        R	0110011	000	0000000
# 2	    c.mv	CR	add rd, x0, rs2	        R	0110011	000	0000000
# 3	    c.jr	CR	jalr x0, 0 (rs1)	    I	1100111	000	 
# 4	    c.jalr	CR	jalr x1, 0 (rs1)	    I	1100111	000	 
# 5	    c.li	CI	addi rd, x0, imm	    I	0010011	000	 
# 6	    c.lui	CI	lui rd, nzimm	        U	0110111	 	 
# 7	    c.addi	CI	addi rd, rd, nzimm	    I	0010011	000	 
# 8	    c.slli	CI	slli rd, rd, shamt	    I	0010011	001	0000000
# 9	    c.lw	CL	lw rd', offset (rs1')	I	0000011	010	 
# 10 	c.sw	CS	sw rs2', offset (rs1')	S	0100011	010	 
# 11	c.and	CA	and rd', rd', rs2'	    R	0110011	111	0000000
# 12	c.or	CA	or rd', rd', rs2'	    R	0110011	110	0000000
# 13	c.xor	CA	xor rd', rd', rs2'	    R	0110011	100	0000000
# 14	c.sub	CA	sub rd', rd', rs2'	    R	0110011	000	0100000
# 15	c.beqz	CB	beq rs1', x0, offset	SB	1100011	000	 
# 16	c.bnez	CB	bne rs1', x0, offset	SB	1100011	001	 
# 17	c.srli	CB	srli rd', rd', shamt	I	0010011	101	0000000
# 18	c.srai	CB	srai rd', rd', shamt	I	0010011	101	0100000
# 19	c.andi	CB	andi rd', rd', imm	    I	0010011	111	 
# 20	c.j	CJ	jal x0, offset	            UJ	1101111	 	 
# 21	c.jal	CJ	jal x1, offset	        UJ	1101111	 	 
# 22	32bits SB-type
# 23	32bits UJ-type
# 24	32bits instruction with no adjustments

# s1 = lines_of_rvc_codes
# s2 is the beginning address of rvc_codes
# s3: type_array (the address of the classification array, each slot is of 1 byte.) 
# s4: normalized_array (each slot is of 4 bytes.)

# t0 = i
# use t1 to adjust the pointer in s2
# use t2 to adjust the pointer in s3
# use t3 to adjust the pointer in s4
.text
main:
    # load lines_of_rvc_codes
    la  a0, lines_of_rvc_codes
    lw  s1, 0(a0) # s1 = lines_of_rvc_codes

    # load rvc_codes
    la  s2, rvc_codes # s2 is the beginning address of rvc_codes
    
    # create two arrays (their sizes are lines_of_rvc_codes) on the heap and save the types and code respectively
    li  a0, 9
    mv  a1, s1
    ecall
    mv  s3, a0 # s3: type_array (the address of the classification array, each slot is of 1 byte.) 
    
    li  a0, 9
    li  t1, 4
    mul a1, a1, t1 # a1 = 4 * a1
    ecall
    mv  s4, a0 # s4: normalized_code (each slot is of 4 bytes.)

    # use for loop to judge the type of the rvc_codes( whether 16 or 32, which type)
    li  t0, 0 # i = 0
    mv  t1, s2 # use t1 to adjust the pointer in s2
    mv  t2, s3 # use t2 to adjust the pointer in s3
    mv  t3, s4 # use t3 to adjust the pointer in s4

    judge_loop: # for(int i = 0; i < lines_of_rvc_codes; i++)
        beq t0, s1, out_judge_loop
        lh t4, 0(t1) # use t4 to save the 16bits 
        andi t5, t4, 3 # use t5 to see whether the last two bits of t2 is 0b11
        li  t6, 3 # let t6 = 3( 0b11)
        bne t5, t6, code_is_16b

        # so the code is 32bits, save the code into normalized_array
        lw  t4, 0(t1)
        sw  t4, 0(t3)
        addi    t3, t3, 4 # move t3 to next word

        li  a0, 0
        # now judge whether it's SB/UJ type (we may discuss later to pass which register as return value)
        jal SB_32b
        jal UJ_32b
        jal OTHER_32b

        judge_end_32b:

        # suppose the return value is in a0, save the type in the type_array
        mv  t4, a0 # since a0 <= 32, its value can be expressed with less than 1 bit
        sb  t4, 0(t2)
        addi    t2, t2, 1 # move t2 to next byte

        # go to the next loop
        addi    t1, t1, 4 # move t1 to next 16 bits
        addi    t0, t0, 1 # i++
        j judge_loop

        # so the code is 16bits, save the code into normalized_array
        code_is_16b: 
            li  t4, 0
            lhu  t4, 0(t1)
            sw  t4, 0(t3)
            addi    t3, t3, 4 # move t3 to next word

            li  a0, 0
            # now judge which type it is (we may discuss later to pass which register as return value)
            jal add_16b
            jal mv_16b
            jal jr_16b
            jal jalr_16b
            jal li_16b
            jal lui_16b
            jal addi_16b
            jal slli_16b
            jal lw_16b
            jal sw_16b
            jal and_16b
            jal or_16b
            jal xor_16b
            jal sub_16b
            jal beqz_16b
            jal bnez_16b
            jal srli_16b
            jal srai_16b
            jal andi_16b
            jal j_16b
            jal jal_16b
            
            judge_end:

            # suppose the return value is in a0, save the type in the type_array
            mv  t4, a0 # since a0 <= 32, its value can be expressed with less than 1 bit
            sb  t4, 0(t2)
            addi    t2, t2, 1 # move t2 to next byte

            # go to the next loop
            addi    t1, t1, 2 # move t1 to next 16 bits
            addi    t0, t0, 1 # i++
            j judge_loop

    out_judge_loop:



    # expansion - remember to update the jump amount
    li  t0, 0
    mv  t2, s3 # use t2 to adjust the pointer in s3
    mv  t3, s4 # use t3 to adjust the pointer in s4
    expansion_loop:
        beq t0, s1, out_expansion_loop
        lb  t4, 0(t2) # type value
        lw  t5, 0(t3) # code

        # t4, t5 --> a7 (temp: s7~s11)
        li  t1, 1
        bne t1, t4, expand_add_exit
        jal expand_add
        expand_add_exit:

        li  t1, 2
        bne t1, t4, expand_mv_exit
        jal expand_mv
        expand_mv_exit:

        li  t1, 3
        bne t1, t4, expand_jr_exit
        jal expand_jr
        expand_jr_exit:

        li  t1, 4
        bne t1, t4, expand_jalr_exit
        jal expand_jalr
        expand_jalr_exit:

        li  t1, 5
        bne t1, t4, expand_li_exit
        jal expand_li
        expand_li_exit:

        li  t1, 6
        bne t1, t4, expand_lui_exit
        jal expand_lui
        expand_lui_exit:

        li  t1, 7
        bne t1, t4, expand_addi_exit
        jal expand_addi
        expand_addi_exit:

        li  t1, 8
        bne t1, t4, expand_slli_exit
        jal expand_slli
        expand_slli_exit:

        li  t1, 9
        bne t1, t4, expand_lw_exit
        jal expand_lw
        expand_lw_exit:

        li  t1, 10
        bne t1, t4, expand_sw_exit
        jal expand_sw
        expand_sw_exit:

        li  t1, 11
        bne t1, t4, expand_and_exit
        jal expand_and
        expand_and_exit:

        li  t1, 12
        bne t1, t4, expand_or_exit
        jal expand_or
        expand_or_exit:

        li  t1, 13
        bne t1, t4, expand_xor_exit
        jal expand_xor
        expand_xor_exit:

        li  t1, 14
        bne t1, t4, expand_sub_exit
        jal expand_sub
        expand_sub_exit:

        li  t1, 15
        bne t1, t4, expand_beqz_exit
        jal expand_beqz
        expand_beqz_exit:
        
        li  t1, 16
        bne t1, t4, expand_bnez_exit
        jal expand_bnez
        expand_bnez_exit:
        
        li  t1, 17
        bne t1, t4, expand_srli_exit
        jal expand_srli
        expand_srli_exit:
        
        li  t1, 18
        bne t1, t4, expand_srai_exit
        jal expand_srai
        expand_srai_exit:
        
        li  t1, 19
        bne t1, t4, expand_andi_exit
        jal expand_andi
        expand_andi_exit:
        
        li  t1, 20
        bne t1, t4, expand_j_exit
        jal expand_j
        expand_j_exit:
        
        li  t1, 21
        bne t1, t4, expand_jal_exit
        jal expand_jal
        expand_jal_exit:
        
        li  t1, 22
        bne t1, t4, modify_SBtype_exit
        jal modify_SBtype
        modify_SBtype_exit:
        
        li  t1, 23
        bne t1, t4, modify_UJtype_exit
        jal modify_UJtype
        modify_UJtype_exit:
        
        li  t1, 24
        bne t1, t4, no_change_exit
        jal no_change
        no_change_exit:
        


        

        # print modified codes to binary, parameter is a7
        jal print_bin
        addi t2, t2, 1 # move the pointer of type array
        addi t3, t3, 4 # move the pointer of code array
        addi t0, t0, 1 # i++

        j expansion_loop

    out_expansion_loop:


    # Exit with error code 0
    add     a1, x0, x0
    addi    a0, x0, 17
    ecall


print_bin:
    mv  s7, a7
    li  s8, 31
    li  s9, -1

    addi sp, sp, -8
    sw  a0, 0(sp)
    sw  a1, 4(sp)
    li  a0, 1
    print_bin_loop:
        beq     s8, s9, print_bin_loop_end
        srl    a1, s7, s8
        andi    a1, a1, 1
        ecall # print
        
        addi    s8, s8, -1
        j print_bin_loop
    print_bin_loop_end:

    li  a1, '\n'
    li 	a0, 11
    ecall

    lw  a0, 0(sp)
    lw  a1, 4(sp)
    addi sp, sp, 4
    
    jr ra






# expansion functions
expand_add:
    # s11 <-- RS1/RD
    li  s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7 

    # s10 <-- RS2
    li  s7, 0b0000000001111100
    and s10, t5, s7
    srli s10, s10, 2

    # s9 = temporary reg storing the result
    li  s9, 0b00000000000000000000000000110011
    # add RD
    slli s11, s11, 7
    add s9, s9, s11
    srli s11, s11, 7
    # add RS1
    slli s11, s11, 15
    add s9, s9, s11
    srli s11, s11, 15
    # add RS2
    slli s10, s10, 20
    add s9, s9, s10
    srli s10, s10, 20

    # all steps are complete, assign s9 to return reg
    mv  a7, s9
    jr ra

expand_mv:
    # s11 <-- RD
    li  s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7 

    # s10 <-- RS2
    li  s7, 0b0000000001111100
    and s10, t5, s7
    srli s10, s10, 2

    # s9 = temporary reg storing the result
    li  s9, 0b00000000000000000000000000110011
    # add RD
    slli s11, s11, 7
    add s9, s9, s11
    srli s11, s11, 7
    # add RS2
    slli s10, s10, 20
    add s9, s9, s10
    srli s10, s10, 20

    # all steps are complete, assign s9 to return reg
    mv  a7, s9
    jr ra

expand_jr: # jalr x0, 0 (rs1)
    # s11 <-- RS1
    li  s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7

    # s9
    li  s9, 0b00000000000000000000000001100111
    # add RS1
    slli s11, s11, 15
    add s9, s9, s11
    srli s11, s11, 15

    mv  a7, s9
    jr ra

expand_jalr: # jalr x1, 0 (rs1)
    # s11 <-- RS1
    li  s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7

    # s9
    li  s9, 0b00000000000000000000000011100111 
    # add RS1
    slli s11, s11, 15
    add s9, s9, s11
    srli s11, s11, 15

    mv  a7, s9
    jr  ra

expand_li: # addi rd, x0, imm
    # s11 <-- rd
    li s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7
    # s10 <-- imm
    li s10, 0
    li s7, 0b0000000001111100
    and s8, t5, s7
    srli s8, s8, 2
    add s10, s10, s8 # add imm[4:0]

    li s7, 0b0001000000000000
    and s8, t5, s7
    srli s8, s8, 7
    add s10, s10, s8 # add imm[5]

    # highest bit complement
    slli    s10, s10, 26 # move to leftmost
    srai    s10, s10, 26 # arithmetic right shift

    # s9
    li  s9, 0b00000000000000000000000000010011
    # add RD
    slli s11, s11, 7
    add s9, s9, s11
    srli s11, s11, 7
    # add imm
    slli s10, s10, 20
    add s9, s9, s10
    srli s10, s10, 20

    mv  a7, s9
    jr  ra

expand_lui: # lui rd, nzimm
    # s10 <-- imm
    li  s7, 0b0000000001111100 # imm[16:12]
    and s10, t5, s7
    srli s10, s10, 2

    li s7, 0b0001000000000000
    and s8, t5, s7
    srli s8, s8, 7
    add s10, s10, s8 # add imm[17]

    # imm highest bit complement
    slli    s10, s10, 26 # move to leftmost
    srai    s10, s10, 26 # arithmetic right shift


    # s11 <-- rd
    li  s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7

    li  s9, 0b00000000000000000000000000110111
    # add imm
    slli s10, s10, 12
    add s9, s9, s10
    srli s10, s10, 12

    # add rd
    slli s11, s11, 7
    add s9, s9, s11
    srli s11, s11, 7

    mv  a7, s9
    jr  ra

expand_addi: # addi rd, rd, nzimm
    # s11 <-- rd
    li s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7
    # s10 <-- imm
    li s10, 0
    li s7, 0b0000000001111100
    and s8, t5, s7
    srli s8, s8, 2
    add s10, s10, s8 # add imm[4:0]

    li s7, 0b0001000000000000
    and s8, t5, s7
    srli s8, s8, 7
    add s10, s10, s8 # add imm[5]

    # highest bit complement
    slli    s10, s10, 26 # move to leftmost
    srai    s10, s10, 26 # arithmetic right shift

    # s9
    li  s9, 0b00000000000000000000000000010011
    # add RD
    slli s11, s11, 7
    add s9, s9, s11
    srli s11, s11, 7
    slli s11, s11, 15
    add s9, s9, s11
    srli s11, s11, 15

    # add imm
    slli s10, s10, 20
    add s9, s9, s10
    srli s10, s10, 20

    mv  a7, s9
    jr  ra

expand_slli: # slli rd, rd, shamt
    # s11 <-- rd
    li s7, 0b0000111110000000
    and s11, t5, s7
    srli s11, s11, 7

    # s10 <-- shamt
    li s10, 0
    li s7, 0b0000000001111100
    and s8, t5, s7
    srli s8, s8, 2
    add s10, s10, s8 # get shamt[4:0]

    li s7, 0b0001000000000000
    and s8, t5, s7
    srli s8, s8, 12
    add s10, s10, s8 # get shamt[5]

    # highest bit complement
    slli    s10, s10, 26 # move to leftmost
    srai    s10, s10, 26 # arithmetic right shift

    li  s9, 0b00000000000000000001000000010011
    # add RS1/RD
    slli s11, s11, 15
    add s9, s9, s11
    srli s11, s11, 15
    slli s11, s11, 7
    add s9, s9, s11
    srli s11, s11, 7

    # add shamt
    slli s10, s10, 20
    add s9, s9, s10
    srli s10, s10, 20

    mv  a7, s9
    jr  ra

expand_lw:
    # s11 <-- offset
    li  s11, 0
    li  s7, 0b0000000001000000
    and s8, t5, s7
    srli s11, s8, 6

    li  s7, 0b0001110000000000
    and s8, t5, s7
    srli s8, s8, 9
    add s11, s11, s8

    li  s7, 0b0000000000100000
    and s8, t5, s7
    srli s8, s8, 1
    add s11, s11, s8

    # s10 <-- rd
    li  s7, 0b0000000000011100
    and s10, t5, s7
    srli s10, s10, 2
    # do mapping on register
    addi s10, s10, 8

    # s9 <-- rs1
    li  s7, 0b0000001110000000
    and s9, t5, s7
    srli s9, s9, 7
    addi s9, s9, 8

    # s8 <-- result
    li  s8, 0b00000000000000000010000000000011
    # add imm
    slli    s11, s11, 22
    add s8, s8, s11
    srli    s11, s11, 22
    # add rs1
    slli    s9, s9, 15
    add s8, s8, s9
    srli    s9, s9, 15
    # add rd
    slli    s10, s10, 7
    add s8, s8, s10
    srli    s10, s10, 7

    mv  a7, s8
    jr  ra


expand_sw:
    # s11 <-- offset
    li  s11, 0
    li  s7, 0b0000000001000000
    and s11, t5, s7
    srli s11, s11, 6

    li  s7, 0b0001110000000000
    and s8, t5, s7
    srli s8, s8, 9
    add s11, s11, s8

    li  s7, 0b0000000000100000
    and s8, t5, s7
    srli s8, s8, 1
    add s11, s11, s8

    # s10 <-- rs2
    li  s7, 0b0000000000011100
    and s10, t5, s7
    srli s10, s10, 2
    # do mapping on register
    addi s10, s10, 8

    # s9 <-- rs1
    li  s7, 0b0000001110000000
    and s9, t5, s7
    srli s9, s9, 7
    addi s9, s9, 8

    # s8 <-- result
    li  s8, 0b00000000000000000010000000100011
    # add rs1
    slli    s9, s9, 15
    add s8, s8, s9
    srli    s9, s9, 15
    # add rs2
    slli    s10, s10, 20
    add s8, s8, s10
    srli    s10, s10, 20
    # add imm[11:5]
    srli    s7, s11, 3
    slli    s7, s7, 25
    add s8, s8, s7

    # add imm[4:0]
    andi    s7, s11, 0b111 # get imm[4:2]
    slli    s7, s7, 9
    add s8, s8, s7

    mv  a7, s8
    jr  ra


# t4 is type value
# t5 is machine code
# a7 is our return value
# t4, t5 --> a7 (temp: s7~s11)

expand_and: # type is and_16b
    li  a7, 0 # a7 = 0

    li  s7, 0b0110011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b11100 
    and s8, t5, s7 # get [4:2] -> rs2
    srli    s8, s8, 2 # get rs2 [2:0]
    addi    s8, s8, 8 # get the real rs2
    slli    s8, s8, 20 # real rs2 is in [24:20]
    add a7, a7, s8 # add rs2

    li  s7, 0b111000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3

    li  s7, 0b1110000000
    and s8, t5, s7 # get [9:7] -> rd
    srli    s8, s8, 7 # get rd [2:0]
    addi    s8, s8, 8 # get the real rs1/rd
    mv s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1
    
    li  s7, 0b00000000000000000000000000000000 # FUNCT7
    add a7, a7, s7 # add FUNCT7

    jr  ra # return 


expand_or:
    li  a7, 0 # a7 = 0

    li  s7, 0b0110011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b11100 
    and s8, t5, s7 # get [4:2] -> rs2
    srli    s8, s8, 2 # get rs2 [2:0]
    addi    s8, s8, 8 # get the real rs2
    slli    s8, s8, 20 # real rs2 is in [24:20]
    add a7, a7, s8 # add rs2

    li  s7, 0b110000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3

    li  s7, 0b1110000000
    and s8, t5, s7 # get [9:7] -> rd
    srli    s8, s8, 7 # get rd [2:0]
    addi    s8, s8, 8 # get the real rs1/rd
    mv s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1
    
    li  s7, 0b00000000000000000000000000000000 # FUNCT7
    add a7, a7, s7 # add FUNCT7

    jr  ra # return 


expand_xor:
    li  a7, 0 # a7 = 0

    li  s7, 0b0110011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b11100 
    and s8, t5, s7 # get [4:2] -> rs2
    srli    s8, s8, 2 # get rs2 [2:0]
    addi    s8, s8, 8 # get the real rs2
    slli    s8, s8, 20 # real rs2 is in [24:20]
    add a7, a7, s8 # add rs2

    li  s7, 0b100000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3

    li  s7, 0b1110000000
    and s8, t5, s7 # get [9:7] -> rd
    srli    s8, s8, 7 # get rd [2:0]
    addi    s8, s8, 8 # get the real rs1/rd
    mv s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1
    
    li  s7, 0b00000000000000000000000000000000 # FUNCT7
    add a7, a7, s7 # add FUNCT7

    jr  ra # return 


expand_sub:
    li  a7, 0 # a7 = 0

    li  s7, 0b0110011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b11100 
    and s8, t5, s7 # get [4:2] -> rs2
    srli    s8, s8, 2 # get rs2 [2:0]
    addi    s8, s8, 8 # get the real rs2
    slli    s8, s8, 20 # real rs2 is in [24:20]
    add a7, a7, s8 # add rs2

    li  s7, 0b000000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3

    li  s7, 0b1110000000
    and s8, t5, s7 # get [9:7] -> rd
    srli    s8, s8, 7 # get rd [2:0]
    addi    s8, s8, 8 # get the real rs1/rd
    mv s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1
    
    li  s7, 0b01000000000000000000000000000000 # FUNCT7
    add a7, a7, s7 # add FUNCT7

    jr  ra # return 


expand_beqz: # use a5-a6, s7-s11
    li  a7, 0 # a7 = 0

    li  s7, 0b1100011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b1110000000
    and s8, t5, s7 # get [9:7] -> rs1
    srli    s8, s8, 7 # get rs1 [2:0]
    addi    s8, s8, 8 # get the real rs1
    slli    s8, s8, 15 # the real rs1 is in [19:15]
    add a7, a7, s8 # add rs1

    li  s7, 0b000000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3

    # rs2 of beqz is 0, don't need to consider it

    # update the offset
    # first load the old offset
    # use s11 to save the offset in 16b_code
    li  s11, 0
    li  s7, 0b11000
    and s8, t5, s7 # get [4:3] -> imm[2:1]
    srli    s8, s8, 2 # real imm[2:1] is in [2:1]
    add s11, s11, s8 # add imm[2:1]

    li  s7, 0b110000000000
    and s8, t5, s7 # get [11:10] -> imm[4:3]
    srli    s8, s8, 7 # real imm[4:3] is in [4:3]
    add s11, s11, s8 # add imm[4:3]

    li  s7, 0b100
    and s8, t5, s7 # get [2] -> imm[5]
    slli    s8, s8, 3 # real imm[5] is in [5]
    add s11, s11, s8 # add imm[5]

    li  s7, 0b1100000
    and s8, t5, s7 # get [6:5] -> imm[7:6]
    slli    s8, s8, 1 # real imm[7:6] is in [7:6]
    add s11, s11, s8 # add imm[7:6]

    li  s7, 0b1000000000000
    and s8, t5, s7 # get [12] -> imm[8]
    srli    s8, s8, 4 # real imm[8] is in [8]
    add s11, s11, s8 # add imm[8]

    # now update the offset
    li  s7, 0b100000000 # negative
    # judge whether imm[8] (s8) is 0b100000000(negative) or 0b000000000(positive)
    beq s8, s7, beqz_negative   # s11 -> offset, a7 -> output, t2 -> points to the type array, s10 -> new offset
    # the offset is zero or positive
    # so the higher bits of s11 are zeros, s11 remains the same
    # since the jump_line is not certain, we can use the offset as stopping condition in while loop
    li  s7, 0 # use s7 to record the offset used
    mv  s8, t2 # use s8 to move the pointer in type array
    mv  s10, s11 # now the new offset = old offset
    positive_while_beqz: # while (s7 != s11)
        beq s7, s11, beqz_16b_offset_update
        li  a6, 22 
        lb  a5, 0(s8) # a5 save the type of the current code
        bge a5, a6, beqz_32b_pos_if # if the type >= 22 then it's 32b_code
        # so it's 16b_code, also need to update s10
        addi    s10, s10, 2
        addi    s8, s8, 1
        addi    s7, s7, 2
        j positive_while_beqz # while loop
        beqz_32b_pos_if: # only need to update s8 and s7, jump to the next loop
            addi    s8, s8, 1
            addi    s7, s7, 4
            j positive_while_beqz # while loop
    beqz_negative: # the offset is negative
        ori s11, s11, 0b11111111111111111111111000000000# so the higher bits of s11 are ones
        li  s7, 0 # use s7 to record the offset used
        mv  s8, t2 # use s8 to move the pointer in type array
        mv  s10, s11 # now the new offset = old offset
        negative_while_beqz: # while (s7 != s11)
            beq s7, s11, beqz_16b_offset_update
            li  a6, 22
            addi, s8, s8, -1 
            lb  a5, 0(s8) # a5 save the type of the former code
            bge a5, a6, beqz_32b_neg_if # if the type >= 22 then it's 32b_code
            # so it's 16b_code, also need to update s10
            addi    s10, s10, -2
            addi    s7, s7, -2
            j negative_while_beqz
            beqz_32b_neg_if: # only need to update s7, jump to the next loop
                addi    s7, s7, -4
                j negative_while_beqz
    beqz_16b_offset_update:
        # update the imm as the format shows
        li  s7, 0b1000000000000
        and s8, s10, s7 # imm[12]
        slli    s8, s8, 19 # 12 -> 31
        add a7, a7, s8 # add imm[12]
        
        li  s7, 0b11111100000
        and s8, s10, s7 # imm[10:5]
        slli    s8, s8, 20 # [10:5] -> [30:25]
        add a7, a7, s8 # add imm[10:5]
        
        li  s7, 0b11110
        and s8, s10, s7 # imm[4:1]
        slli    s8, s8, 7 # [4:1] -> [11:8]
        add a7, a7, s8 # add imm[4:1]

        li  s7, 0b100000000000
        and s8, s10, s7 # imm[11]
        srli    s8, s8, 4 # 11 -> 7
        add a7, a7, s8 # add imm[11]
        jr  ra # return     


expand_bnez: # use a5-a6, s7-s11
    li  a7, 0 # a7 = 0

    li  s7, 0b1100011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b1110000000
    and s8, t5, s7 # get [9:7] -> rs1
    srli    s8, s8, 7 # get rs1 [2:0]
    addi    s8, s8, 8 # get the real rs1
    slli    s8, s8, 15 # the real rs1 is in [19:15]
    add a7, a7, s8 # add rs1

    li  s7, 0b001000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3

    # rs2 of bnez is 0, don't need to consider it

    # update the offset
    # first load the old offset
    # use s11 to save the offset in 16b_code
    li  s11, 0
    li  s7, 0b11000
    and s8, t5, s7 # get [4:3] -> imm[2:1]
    srli    s8, s8, 2 # real imm[2:1] is in [2:1]
    add s11, s11, s8 # add imm[2:1]

    li  s7, 0b110000000000
    and s8, t5, s7 # get [11:10] -> imm[4:3]
    srli    s8, s8, 7 # real imm[4:3] is in [4:3]
    add s11, s11, s8 # add imm[4:3]

    li  s7, 0b100
    and s8, t5, s7 # get [2] -> imm[5]
    slli    s8, s8, 3 # real imm[5] is in [5]
    add s11, s11, s8 # add imm[5]

    li  s7, 0b1100000
    and s8, t5, s7 # get [6:5] -> imm[7:6]
    slli    s8, s8, 1 # real imm[7:6] is in [7:6]
    add s11, s11, s8 # add imm[7:6]

    li  s7, 0b1000000000000
    and s8, t5, s7 # get [12] -> imm[8]
    srli    s8, s8, 4 # real imm[8] is in [8]
    add s11, s11, s8 # add imm[8]

    # now update the offset
    li  s7, 0b100000000 # negative
    # judge whether imm[8] (s8) is 0b100000000(negative) or 0b000000000(positive)
    beq s8, s7, bnez_negative   # s11 -> offset, a7 -> output, t2 -> points to the type array, s10 -> new offset
    # the offset is zero or positive
    # so the higher bits of s11 are zeros, s11 remains the same
    # since the jump_line is not certain, we can use the offset as stopping condition in while loop
    li  s7, 0 # use s7 to record the offset used
    mv  s8, t2 # use s8 to move the pointer in type array
    mv  s10, s11 # now the new offset = old offset
    positive_while_bnez: # while (s7 != s11)
        beq s7, s11, bnez_16b_offset_update
        li  a6, 22 
        lb  a5, 0(s8) # a5 save the type of the current code
        bge a5, a6, bnez_32b_pos_if # if the type >= 22 then it's 32b_code
        # so it's 16b_code, also need to update s10
        addi    s10, s10, 2
        addi    s8, s8, 1
        addi    s7, s7, 2
        j positive_while_bnez # while loop
        bnez_32b_pos_if: # only need to update s8 and s7, jump to the next loop
            addi    s8, s8, 1
            addi    s7, s7, 4
            j positive_while_bnez # while loop
    bnez_negative: # the offset is negative
        ori s11, s11, 0b11111111111111111111111000000000# so the higher bits of s11 are ones
        li  s7, 0 # use s7 to record the offset used
        mv  s8, t2 # use s8 to move the pointer in type array
        mv  s10, s11 # now the new offset = old offset
        negative_while_bnez: # while (s7 != s11)
            beq s7, s11, bnez_16b_offset_update
            li  a6, 22
            addi, s8, s8, -1 
            lb  a5, 0(s8) # a5 save the type of the former code
            bge a5, a6, bnez_32b_neg_if # if the type >= 22 then it's 32b_code
            # so it's 16b_code, also need to update s10
            addi    s10, s10, -2
            addi    s7, s7, -2
            j negative_while_bnez
            bnez_32b_neg_if: # only need to update s7, jump to the next loop
                addi    s7, s7, -4
                j negative_while_bnez
    bnez_16b_offset_update:
        # update the imm as the format shows
        li  s7, 0b1000000000000
        and s8, s10, s7 # imm[12]
        slli    s8, s8, 19 # 12 -> 31
        add a7, a7, s8 # add imm[12]
        
        li  s7, 0b11111100000
        and s8, s10, s7 # imm[10:5]
        slli    s8, s8, 20 # [10:5] -> [30:25]
        add a7, a7, s8 # add imm[10:5]
        
        li  s7, 0b11110
        and s8, s10, s7 # imm[4:1]
        slli    s8, s8, 7 # [4:1] -> [11:8]
        add a7, a7, s8 # add imm[4:1]

        li  s7, 0b100000000000
        and s8, s10, s7 # imm[11]
        srli    s8, s8, 4 # 11 -> 7
        add a7, a7, s8 # add imm[11]
        jr  ra # return         


expand_srli: 
    li  a7, 0 # a7 = 0

    li  s7, 0b0010011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b1110000000 
    and s8, t5, s7 # get [9:7] -> rd/rs1
    srli    s8, s8, 7 # get rd/rs1 [2:0]
    addi    s8, s8, 8 # get the real rd/rs1
    mv  s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1

    li  s7, 0b101000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3
    
    li  s7, 0b1111100
    and s8, t5, s7 # get [6:2] -> shamt[4:0]
    # shamt[5] = 0, so don't need to consider it
    slli    s8, s8, 18 # real shamt is in [24:20]
    add a7, a7, s8 # add shamt

    li  s7, 0b00000000000000000000000000000000 # FUNCT7
    add a7, a7, s7 # add FUNCT7

    jr  ra # return 


expand_srai:
    li  a7, 0 # a7 = 0

    li  s7, 0b0010011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b1110000000 
    and s8, t5, s7 # get [9:7] -> rd/rs1
    srli    s8, s8, 7 # get rd/rs1 [2:0]
    addi    s8, s8, 8 # get the real rd/rs1
    mv  s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1

    li  s7, 0b101000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3
    
    li  s7, 0b1111100
    and s8, t5, s7 # get [6:2] -> shamt[4:0]
    # shamt[5] = 0, so don't need to consider it
    slli    s8, s8, 18 # real shamt is in [24:20]
    add a7, a7, s8 # add shamt

    li  s7, 0b01000000000000000000000000000000 # FUNCT7
    add a7, a7, s7 # add FUNCT7

    jr  ra # return 


expand_andi:
    li  a7, 0 # a7 = 0

    li  s7, 0b0010011 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b1110000000 
    and s8, t5, s7 # get [9:7] -> rd/rs1
    srli    s8, s8, 7 # get rd/rs1 [2:0]
    addi    s8, s8, 8 # get the real rd/rs1
    mv  s9, s8 # s9 = s8
    slli    s8, s8, 7 # real rd is in [11:7]
    slli    s9, s9, 15 # real rs1 is in [19:15]
    add a7, a7, s8 # add rd
    add a7, a7, s9 # add rs1

    li  s7, 0b111000000000000 # FUNCT3
    add a7, a7, s7 # add FUNCT3
    
    li  s7, 0b1111100
    and s8, t5, s7 # get[6:2] -> imm[4:0]
    srli    s8, s8, 2 # get imm[4:0] [4:0]
    li  s7, 0b1000000000000
    and s9, t5, s7 # get [12] -> imm[5]
    srli    s9, s9, 7 # get imm[5] [5]
    # check the highest bit
    li  s7, 0b100000
    beq s9, s7, andi_negative 
    # so the imm isn't negative
    add s9, s9, s8 # s9 = imm[5:0] in [5:0]
    slli    s9, s9, 20 # real imm[11:0] is in [31:20]
    add a7, a7, s9 # add imm[11:0]
    jr  ra # return 

    andi_negative: # the imm is negative
        add s9, s9, s8 # s9 = imm[5:0] in [5:0]
        li  a4, 0b111111000000
        or  s9, s9, a4 # imm now signed extended [11:0]
        slli    s9, s9, 20 # real imm[11:0] is in [31:20]
        add a7, a7, s9 # add imm[11:0]
        jr  ra # return 


expand_j: # use a5-a6, s7-s11
    li  a7, 0 # a7 = 0

    li  s7, 0b1101111 # opcode
    add a7, a7, s7 # add opcode

    # rd of j is 0, don't need to consider it

    # update the offset
    # first load the old offset
    # use s11 to save the offset in 16b_code
    li  s11, 0
    li  s7, 0b111000
    and s8, t5, s7 # get [5:3] -> imm[3:1]
    srli    s8, s8, 2 # real imm[3:1] is in [3:1]
    add s11, s11, s8 # add imm[3:1]

    li  s7, 0b100000000000
    and s8, t5, s7 # get [11] -> imm[4]
    srli    s8, s8, 7 # real imm[4] is in [4]
    add s11, s11, s8 # add imm[4]

    li  s7, 0b100
    and s8, t5, s7 # get [2] -> imm[5]
    slli    s8, s8, 3 # real imm[5] is in [5]
    add s11, s11, s8 # add imm[5]

    li  s7, 0b10000000
    and s8, t5, s7 # get [7] -> imm[6]
    srli    s8, s8, 1 # real imm[6] is in [6]
    add s11, s11, s8 # add imm[6]

    li  s7, 0b1000000
    and s8, t5, s7 # get [6] -> imm[7]
    slli    s8, s8, 1 # real imm[7] is in [7]
    add s11, s11, s8 # add imm[7]

    li  s7, 0b11000000000
    and s8, t5, s7 # get [10:9] -> imm[9:8]
    srli    s8, s8, 1 # real imm[9:8] is in [9:8]
    add s11, s11, s8 # add imm[9:8]

    li  s7, 0b100000000
    and s8, t5, s7 # get [8] -> imm[10]
    slli    s8, s8, 2 # real imm[10] is in [10]
    add s11, s11, s8 # add imm[10]

    li  s7, 0b1000000000000
    and s8, t5, s7 # get [12] -> imm[11]
    srli    s8, s8, 1 # real imm[11] is in [11]
    add s11, s11, s8 # add imm[11]

    # now update the offset
    li  s7, 0b100000000000 # negative
    # judge whether imm[8] (s8) is 0b100000000(negative) or 0b000000000(positive)
    beq s8, s7, j_negative   # s11 -> offset, a7 -> output, t2 -> points to the type array, s10 -> new offset
    # the offset is zero or positive
    # so the higher bits of s11 are zeros, s11 remains the same
    # since the jump_line is not certain, we can use the offset as stopping condition in while loop
    li  s7, 0 # use s7 to record the offset used
    mv  s8, t2 # use s8 to move the pointer in type array
    mv  s10, s11 # now the new offset = old offset
    positive_while_j: # while (s7 != s11)
        beq s7, s11, j_16b_offset_update
        li  a6, 22 
        lb  a5, 0(s8) # a5 save the type of the current code
        bge a5, a6, j_32b_pos_if # if the type >= 22 then it's 32b_code
        # so it's 16b_code, also need to update s10
        addi    s10, s10, 2
        addi    s8, s8, 1
        addi    s7, s7, 2
        j positive_while_j # while loop
        j_32b_pos_if: # only need to update s8 and s7, jump to the next loop
            addi    s8, s8, 1
            addi    s7, s7, 4
            j positive_while_j # while loop
    j_negative: # the offset is negative
        li  a4, 0b11111111111111111111000000000000
        or s11, s11, a4 # so the higher bits of s11 are ones
        li  s7, 0 # use s7 to record the offset used
        mv  s8, t2 # use s8 to move the pointer in type array
        mv  s10, s11 # now the new offset = old offset
        negative_while_j: # while (s7 != s11)
            beq s7, s11, j_16b_offset_update
            li  a6, 22
            addi, s8, s8, -1 
            lb  a5, 0(s8) # a5 save the type of the former code
            bge a5, a6, j_32b_neg_if # if the type >= 22 then it's 32b_code
            # so it's 16b_code, also need to update s10
            addi    s10, s10, -2
            addi    s7, s7, -2
            j negative_while_j
            j_32b_neg_if: # only need to update s7, jump to the next loop
                addi    s7, s7, -4
                j negative_while_j
    j_16b_offset_update:
        # update the imm as the format shows
        li  s7, 0b100000000000000000000
        and s8, s10, s7 # imm[20]
        slli    s8, s8, 11 # 20 -> 31
        add a7, a7, s8 # add imm[20]

        li  s7, 0b11111111110
        and s8, s10, s7 # imm[10:1]
        slli    s8, s8, 20 # [10:1] -> [30:21]
        add a7, a7, s8 # add imm[10:1]

        li  s7, 0b100000000000
        and s8, s10, s7 # imm[11]
        slli    s8, s8, 9 # 11 -> 20
        add a7, a7, s8 # add imm[11]

        li  s7, 0b11111111000000000000
        and s8, s10, s7 # imm[19:12]
        # don't need to shift
        add a7, a7, s8 # add imm[19:12]
        jr  ra # return            


expand_jal: # use a5-a6, s7-s11
    li  a7, 0 # a7 = 0

    li  s7, 0b1101111 # opcode
    add a7, a7, s7 # add opcode

    li  s7, 0b10000000 # rd is x1
    add a7, a7, s7 # add rd
    # rd of j is 0, don't need to consider it

    # update the offset
    # first load the old offset
    # use s11 to save the offset in 16b_code
    li  s11, 0
    li  s7, 0b111000
    and s8, t5, s7 # get [5:3] -> imm[3:1]
    srli    s8, s8, 2 # real imm[3:1] is in [3:1]
    add s11, s11, s8 # add imm[3:1]

    li  s7, 0b100000000000
    and s8, t5, s7 # get [11] -> imm[4]
    srli    s8, s8, 7 # real imm[4] is in [4]
    add s11, s11, s8 # add imm[4]

    li  s7, 0b100
    and s8, t5, s7 # get [2] -> imm[5]
    slli    s8, s8, 3 # real imm[5] is in [5]
    add s11, s11, s8 # add imm[5]

    li  s7, 0b10000000
    and s8, t5, s7 # get [7] -> imm[6]
    srli    s8, s8, 1 # real imm[6] is in [6]
    add s11, s11, s8 # add imm[6]

    li  s7, 0b1000000
    and s8, t5, s7 # get [6] -> imm[7]
    slli    s8, s8, 1 # real imm[7] is in [7]
    add s11, s11, s8 # add imm[7]

    li  s7, 0b11000000000
    and s8, t5, s7 # get [10:9] -> imm[9:8]
    srli    s8, s8, 1 # real imm[9:8] is in [9:8]
    add s11, s11, s8 # add imm[9:8]

    li  s7, 0b100000000
    and s8, t5, s7 # get [8] -> imm[10]
    slli    s8, s8, 2 # real imm[10] is in [10]
    add s11, s11, s8 # add imm[10]

    li  s7, 0b1000000000000
    and s8, t5, s7 # get [12] -> imm[11]
    srli    s8, s8, 1 # real imm[11] is in [11]
    add s11, s11, s8 # add imm[11]

    # now update the offset
    li  s7, 0b100000000000 # negative
    # judge whether imm[8] (s8) is 0b100000000(negative) or 0b000000000(positive)
    beq s8, s7, jal_negative   # s11 -> offset, a7 -> output, t2 -> points to the type array, s10 -> new offset
    # the offset is zero or positive
    # so the higher bits of s11 are zeros, s11 remains the same
    # since the jump_line is not certain, we can use the offset as stopping condition in while loop
    li  s7, 0 # use s7 to record the offset used
    mv  s8, t2 # use s8 to move the pointer in type array
    mv  s10, s11 # now the new offset = old offset
    positive_while_jal: # while (s7 != s11)
        beq s7, s11, jal_16b_offset_update
        li  a6, 22 
        lb  a5, 0(s8) # a5 save the type of the current code
        bge a5, a6, jal_32b_pos_if # if the type >= 22 then it's 32b_code
        # so it's 16b_code, also need to update s10
        addi    s10, s10, 2
        addi    s8, s8, 1
        addi    s7, s7, 2
        j positive_while_jal # while loop
        jal_32b_pos_if: # only need to update s8 and s7, jump to the next loop
            addi    s8, s8, 1
            addi    s7, s7, 4
            j positive_while_jal # while loop
    jal_negative: # the offset is negative
        li  a4, 0b11111111111111111111000000000000
        or s11, s11, a4# so the higher bits of s11 are ones
        li  s7, 0 # use s7 to record the offset used
        mv  s8, t2 # use s8 to move the pointer in type array
        mv  s10, s11 # now the new offset = old o0b11111111111111111111111000000000he former code
        negative_while_jal: # while (s7 != s11)
            beq s7, s11, jal_16b_offset_update
            li  a6, 22 
            lb  a5, 0(s8) # a5 save the type of the former code
            bge a5, a6, jal_32b_neg_if # if the type >= 22 then it's 32b_code
            # so it's 16b_code, also need to update s10
            addi    s10, s10, -2
            addi    s7, s7, -2
            j negative_while_jal
            jal_32b_neg_if: # only need to update s7, jump to the next loop
                addi    s7, s7, -4
                j negative_while_jal
    jal_16b_offset_update:
        # update the imm as the format shows
        li  s7, 0b100000000000000000000
        and s8, s10, s7 # imm[20]
        slli    s8, s8, 11 # 20 -> 31
        add a7, a7, s8 # add imm[20]

        li  s7, 0b11111111110
        and s8, s10, s7 # imm[10:1]
        slli    s8, s8, 20 # [10:1] -> [30:21]
        add a7, a7, s8 # add imm[10:1]

        li  s7, 0b100000000000
        and s8, s10, s7 # imm[11]
        slli    s8, s8, 9 # 11 -> 20
        add a7, a7, s8 # add imm[11]

        li  s7, 0b11111111000000000000
        and s8, s10, s7 # imm[19:12]
        # don't need to shift
        add a7, a7, s8 # add imm[19:12]
        jr  ra # return    


modify_SBtype: # use a5-a6, s7-s11
    # update the offset
    # first load the old offset
    # use s11 to save the offset in 32b_code
    li  s11, 0
    li  s7, 0b10000000
    and s8, t5, s7 # get [7] -> imm[11]
    slli    s8, s8, 4 # real imm[11] is in [11]
    add s11, s11, s8 # add imm[11]

    li  s7, 0b111100000000
    and s8, t5, s7 # get [11:8] -> imm[4:1]
    srli    s8, s8, 7 # real imm[4:1] is in [4:1]
    add s11, s11, s8 # add imm[4:1]

    li  s7, 0b1111110000000000000000000000000
    and s8, t5, s7 # get [30:25] -> imm[10:5]
    srli    s8, s8, 20 # real imm[10:5] is in [10:5]
    add s11, s11, s8 # add imm[10:5]

    li  s7, 0b10000000000000000000000000000000
    and s8, t5, s7 # get [31] -> imm[12]
    srli    s8, s8, 19 # real imm[12] is in [12]
    add s11, s11, s8 # add imm[12]

    # now update the offset
    li  s7, 0b1000000000000 # negative
    # judge whether imm[8] (s8) is 0b100000000(negative) or 0b000000000(positive)
    beq s8, s7, SB_negative   # s11 -> offset, a7 -> output, t2 -> points to the type array, s10 -> new offset
    # the offset is zero or positive
    # so the higher bits of s11 are zeros, s11 remains the same
    # since the jump_line is not certain, we can use the offset as stopping condition in while loop
    li  s7, 0 # use s7 to record the offset used
    mv  s8, t2 # use s8 to move the pointer in type array
    mv  s10, s11 # now the new offset = old offset
    positive_while_SB: # while (s7 != s11)
        beq s7, s11, SB_16b_offset_update
        li  a6, 22 
        lb  a5, 0(s8) # a5 save the type of the current code
        bge a5, a6, SB_32b_pos_if # if the type >= 22 then it's 32b_code
        # so it's 16b_code, also need to update s10
        addi    s10, s10, 2
        addi    s8, s8, 1
        addi    s7, s7, 2
        j positive_while_SB # while loop
        SB_32b_pos_if: # only need to update s8 and s7, jump to the next loop
            addi    s8, s8, 1
            addi    s7, s7, 4
            j positive_while_SB # while loop
    SB_negative: # the offset is negative
        li  a4, 0b11111111111111111110000000000000
        or  s11, s11, a4 # so the higher bits of s11 are ones
        li  s7, 0 # use s7 to record the offset used
        mv  s8, t2 # use s8 to move the pointer in type array
        mv  s10, s11 # now the new offset = old offset
        negative_while_SB: # while (s7 != s11)
            beq s7, s11, SB_16b_offset_update
            li  a6, 22
            addi, s8, s8, -1 
            lb  a5, 0(s8) # a5 save the type of the former code
            bge a5, a6, SB_32b_neg_if # if the type >= 22 then it's 32b_code
            # so it's 16b_code, also need to update s10
            addi    s10, s10, -2
            addi    s7, s7, -2
            j negative_while_SB
            SB_32b_neg_if: # only need to update s7, jump to the next loop
                addi    s7, s7, -4
                j negative_while_SB
    SB_16b_offset_update:
        mv  a7, t5
        li  s7, 0b00000001111111111111000001111111 # ignore the imm in a7
        and a7, a7, s7 # ignore
        # update the imm as the format shows
        li  s7, 0b1000000000000
        and s8, s10, s7 # imm[12]
        slli    s8, s8, 19 # 12 -> 31
        add a7, a7, s8 # add imm[12]
        
        li  s7, 0b11111100000
        and s8, s10, s7 # imm[10:5]
        slli    s8, s8, 20 # [10:5] -> [30:25]
        add a7, a7, s8 # add imm[10:5]
        
        li  s7, 0b11110
        and s8, s10, s7 # imm[4:1]
        slli    s8, s8, 7 # [4:1] -> [11:8]
        add a7, a7, s8 # add imm[4:1]

        li  s7, 0b100000000000
        and s8, s10, s7 # imm[11]
        srli    s8, s8, 4 # 11 -> 7
        add a7, a7, s8 # add imm[11]
        jr  ra # return


modify_UJtype: # use a5-a6, s7-s11
    # update the offset
    # first load the old offset
    # use s11 to save the offset in 16b_code
    li  s11, 0
    li  s7, 0b10000000000000000000000000000000
    and s8, t5, s7 # get [31] -> imm[20]
    srli    s8, s8, 11 # real imm[20] is in [20]
    add s11, s11, s8 # add imm[20]

    li  s7, 0b1111111111000000000000000000000
    and s8, t5, s7 # get [30:21] -> imm[10:1]
    srli    s8, s8, 20 # real imm[10:1] is in [10:1]
    add s11, s11, s8 # add imm[10:1]

    li  s7, 0b100000000000000000000
    and s8, t5, s7 # get [20] -> imm[11]
    srli    s8, s8, 9 # real imm[11] is in [11]
    add s11, s11, s8 # add imm[11]

    li  s7, 0b11111111000000000000
    and s8, t5, s7 # get [19：12] -> imm[19：12]
    add s11, s11, s8 # add imm[19：12]
    
    # now update the offset
    li  s7, 0b1000000000000 # negative
    # judge whether imm[8] (s8) is 0b100000000(negative) or 0b000000000(positive)
    beq s8, s7, UJ_negative   # s11 -> offset, a7 -> output, t2 -> points to the type array, s10 -> new offset
    # the offset is zero or positive
    # so the higher bits of s11 are zeros, s11 remains the same
    # since the jump_line is not certain, we can use the offset as stopping condition in while loop
    li  s7, 0 # use s7 to record the offset used
    mv  s8, t2 # use s8 to move the pointer in type array
    mv  s10, s11 # now the new offset = old offset
    positive_while_UJ: # while (s7 != s11)
        beq s7, s11, UJ_16b_offset_update
        li  a6, 22 
        lb  a5, 0(s8) # a5 save the type of the current code
        bge a5, a6, UJ_32b_pos_if # if the type >= 22 then it's 32b_code
        # so it's 16b_code, also need to update s10
        addi    s10, s10, 2
        addi    s8, s8, 1
        addi    s7, s7, 2
        j positive_while_UJ # while loop
        UJ_32b_pos_if: # only need to update s8 and s7, jump to the next loop
            addi    s8, s8, 1
            addi    s7, s7, 4
            j positive_while_UJ # while loop
    UJ_negative: # the offset is negative
        li  a4, 0b11111111111111111110000000000000
        or s11, s11, a4 # so the higher bits of s11 are ones
        li  s7, 0 # use s7 to record the offset used
        mv  s8, t2 # use s8 to move the pointer in type array
        mv  s10, s11 # now the new offset = old o0b11111111111111111111111000000000he former code
        negative_while_UJ:
        beq s7, s11, UJ_16b_offset_update
            li  a6, 22
            addi, s8, s8, -1
            lb  a5, 0(s8) # a5 saves the type of the former code
            bge a5, a6, UJ_32b_neg_if # if the type >= 22 then it's 32b_code
            # so it's 16b_code, also need to update s10
            addi    s10, s10, -2
            addi    s7, s7, -2
            j negative_while_UJ
            UJ_32b_neg_if: # only need to update s7, jump to the next loop
                addi    s7, s7, -4
                j negative_while_UJ
    UJ_16b_offset_update:
        mv  a7, t5 # ignore the imm
        li  s7, 0b00000000000000000000111111111111
        and a7, a7, s7 # ignore the imm
        # update the imm as the format shows
        li  s7, 0b100000000000000000000
        and s8, s10, s7 # imm[20]
        slli    s8, s8, 11 # 20 -> 31
        add a7, a7, s8 # add imm[20]

        li  s7, 0b11111111110
        and s8, s10, s7 # imm[10:1]
        slli    s8, s8, 20 # [10:1] -> [30:21]
        add a7, a7, s8 # add imm[10:1]

        li  s7, 0b100000000000
        and s8, s10, s7 # imm[11]
        slli    s8, s8, 9 # 11 -> 20
        add a7, a7, s8 # add imm[11]

        li  s7, 0b11111111000000000000
        and s8, s10, s7 # imm[19:12]
        # don't need to shift
        add a7, a7, s8 # add imm[19:12]
        jr  ra # return  


no_change: 
    mv  a7, t5
    jr  ra # return     






# classification functions


add_16b: # s7 - s11
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b10
    bne s11, s7, add_16b_fail # check opcode 10
    li  s10, 0b1111000000000000 # 1111 0000 0000 0000
    and s9, t4, s10 # get [15:12]
    li  s7, 0b1001000000000000 # 1001 0000 0000 0000
    bne s9, s7, add_16b_fail # check FUNCT4 1001

    li  s10, 0b0000000001111100
    and s9, t4, s10 # get [6:2]
    li  s7, 0b0000000000000000
    beq s9, s7, add_16b_fail # check RS2 != x0

    li a0, 1
    j judge_end # classification succeeds

    add_16b_fail:
        jr ra # return


mv_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b10
    bne s11, s7, mv_16b_fail # check opcode 10
    li  s10, 0b1111000000000000 # 1111 0000 0000 0000
    and s9, t4, s10 # get [15:12]
    li  s7, 0b1000000000000000
    bne s9, s7, mv_16b_fail # check FUNCT4 1000
    li  s8, 0b0000000001111100
    and s9, s8, t4
    srli    s9, s9, 2
    beq s9, x0, mv_16b_fail

    li  s10, 0b0000000001111100
    and s9, t4, s10 # get [6:2]
    li  s7, 0b0000000000000000
    beq s9, s7, mv_16b_fail # check RS2 != x0


    li a0, 2
    j judge_end # classification succeeds

    mv_16b_fail:
        jr ra # return


jr_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b10
    bne s11, s7, jr_16b_fail # check opcode 10
    li  s10, 0b1111000000000000 # 1111 0000 0000 0000
    and s9, t4, s10 # get [15:12]
    li  s7, 0b1000000000000000
    bne s9, s7, jr_16b_fail # check FUNCT4 1000
    
    li  s10, 0b0000000001111100
    and s9, t4, s10 # get [6:2]
    li  s7, 0b0000000000000000
    bne s9, s7, jr_16b_fail # check RS2 0

    li a0, 3
    j judge_end # classification succeeds

    jr_16b_fail:
        jr ra # return

jalr_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b10
    bne s11, s7, jalr_16b_fail # check opcode 10
    li  s10, 0b1111000000000000 # 1111 0000 0000 0000
    and s9, t4, s10 # get [15:12]
    li  s7, 0b1001000000000000
    bne s9, s7, jalr_16b_fail # check FUNCT4 1001
    
    li  s10, 0b0000000001111100
    and s9, t4, s10 # get [6:2]
    li  s7, 0b0000000000000000
    bne s9, s7, jalr_16b_fail # check RS2 x1

    li a0, 4
    j judge_end # classification succeeds

    jalr_16b_fail:
        jr ra # return


li_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, li_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b0100000000000000
    bne s9, s7, li_16b_fail # check FUNCT3 010

    li a0, 5
    j judge_end # judge succeeds

    li_16b_fail:
        jr ra # return

lui_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, lui_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b0110000000000000
    bne s9, s7, lui_16b_fail # check FUNCT3 011

    li a0, 6
    j judge_end # judge succeeds

    lui_16b_fail:
        jr ra # return


addi_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, addi_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b0000000000000000
    bne s9, s7, addi_16b_fail # check FUNCT3 000

    li a0, 7
    j judge_end # judge succeeds

    addi_16b_fail:
        jr ra # return


slli_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b10
    bne s11, s7, slli_16b_fail # check OPCODE 10

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b0000000000000000
    bne s9, s7, slli_16b_fail # check FUNCT3 000

    li a0, 8
    j judge_end # judge succeeds

    slli_16b_fail:
        jr ra # return

lw_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b00
    bne s11, s7, lw_16b_fail # check OPCODE 00

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b0100000000000000
    bne s9, s7, lw_16b_fail # check FUNCT3 010

    li a0, 9
    j judge_end # judge succeeds

    lw_16b_fail:
        jr ra # return

sw_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b00
    bne s11, s7, sw_16b_fail # check OPCODE 00

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1100000000000000
    bne s9, s7, sw_16b_fail # check FUNCT3 110

    li a0, 10
    j judge_end # judge succeeds

    sw_16b_fail:
        jr ra # return

and_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, and_16b_fail # check OPCODE 01

    li  s10, 0b1111110000000000
    and s9, t4, s10 # get [15:10]
    li  s7, 0b1000110000000000
    bne s9, s7, and_16b_fail # check FUNCT6 100011

    li  s10, 0b1100000
    and s9, t4, s10 # get [6:5]
    li  s7, 0b1100000
    bne s9, s7, and_16b_fail # check FUNCT2 11

    li a0, 11
    j judge_end # judge succeeds

    and_16b_fail:
        jr ra # return


or_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, or_16b_fail # check OPCODE 01

    li  s10, 0b1111110000000000
    and s9, t4, s10 # get [15:10]
    li  s7, 0b1000110000000000
    bne s9, s7, or_16b_fail # check FUNCT6 100011

    li  s10, 0b1100000
    and s9, t4, s10 # get [6:5]
    li  s7, 0b1000000
    bne s9, s7, or_16b_fail # check FUNCT2 10
    
    li a0, 12
    j judge_end # judge succeeds

    or_16b_fail:
        jr ra # return



xor_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, xor_16b_fail # check OPCODE 01

    li  s10, 0b1111110000000000
    and s9, t4, s10 # get [15:10]
    li  s7, 0b1000110000000000
    bne s9, s7, xor_16b_fail # check FUNCT6 100011

    li  s10, 0b1100000
    and s9, t4, s10 # get [6:5]
    li  s7, 0b0100000
    bne s9, s7, xor_16b_fail # check FUNCT2 01
    
    li a0, 13
    j judge_end # judge succeeds

    xor_16b_fail:
        jr ra # return



sub_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, sub_16b_fail # check OPCODE 01

    li  s10, 0b1111110000000000
    and s9, t4, s10 # get [15:10]
    li  s7, 0b1000110000000000
    bne s9, s7, sub_16b_fail # check FUNCT6 100011

    li  s10, 0b1100000
    and s9, t4, s10 # get [6:5]
    li  s7, 0b0000000
    bne s9, s7, sub_16b_fail # check FUNCT2 00
    
    li a0, 14
    j judge_end # judge succeeds

    sub_16b_fail:
        jr ra # return



beqz_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, beqz_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1100000000000000
    bne s9, s7, beqz_16b_fail # check FUNCT3 110

    li a0, 15
    j judge_end # judge succeeds

    beqz_16b_fail:
        jr ra # return



bnez_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, bnez_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1110000000000000
    bne s9, s7, bnez_16b_fail # check FUNCT3 111

    li a0, 16
    j judge_end # judge succeeds

    bnez_16b_fail:
        jr ra # return



srli_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, srli_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1000000000000000
    bne s9, s7, srli_16b_fail # check FUNCT3 100

    li  s10, 0b110000000000
    and s9, t4, s10 # get [11:10]
    li  s7, 0b000000000000
    bne s9, s7, srli_16b_fail # check FUNCT2 00
    
    li  a0, 17
    j judge_end # judge succeeds

    srli_16b_fail:
        jr ra # return



srai_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, srai_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1000000000000000
    bne s9, s7, srai_16b_fail # check FUNCT3 100

    li  s10, 0b110000000000
    and s9, t4, s10 # get [11:10]
    li  s7, 0b010000000000
    bne s9, s7, srai_16b_fail # check FUNCT2 01
    
    li  a0, 18
    j judge_end # judge succeeds

    srai_16b_fail:
        jr ra # return



andi_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, andi_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1000000000000000
    bne s9, s7, andi_16b_fail # check FUNCT3 100

    li  s10, 0b110000000000
    and s9, t4, s10 # get [11:10]
    li  s7, 0b100000000000
    bne s9, s7, andi_16b_fail # check FUNCT2 10
    
    li  a0, 19
    j judge_end # judge succeeds

    andi_16b_fail:
        jr ra # return



j_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, j_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b1010000000000000
    bne s9, s7, j_16b_fail # check FUNCT3 101

    li  a0, 20
    j judge_end # judge succeeds

    j_16b_fail:
        jr ra # return



jal_16b:
    andi s11, t4, 0b11 # get [1:0]
    li  s7, 0b01
    bne s11, s7, jal_16b_fail # check OPCODE 01

    li  s10, 0b1110000000000000
    and s9, t4, s10 # get [15:13]
    li  s7, 0b0010000000000000
    bne s9, s7, jal_16b_fail # check FUNCT3 001

    li  a0, 21
    j judge_end # judge succeeds

    jal_16b_fail:
        jr ra # return

SB_32b:
    li  s7, 0b1111111
    and s11, t4, s7 # get [6:0]
    li  s7, 0b1100011
    bne s11, s7, SB_32b_fail # check opcode 1100011
    
    li  a0, 22
    j judge_end_32b # judge succeeds

    SB_32b_fail:
    jr ra # return

UJ_32b:
    li  s7, 0b1111111
    and s11, t4, s7 # get [6:0]
    li  s7, 0b1101111
    bne s11, s7, UJ_32b_fail # check opcode 1101111
        
    li  a0, 23
    j judge_end_32b # judge succeeds

    UJ_32b_fail:
    jr ra # return

OTHER_32b:
    li  a0, 24
    jr ra # return


