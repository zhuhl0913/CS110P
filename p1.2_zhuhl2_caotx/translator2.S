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
        ori s9, s9, 0b111111000000 # imm now signed extended [11:0]
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
        ori s11, s11, 0b11111111111111111111000000000000# so the higher bits of s11 are ones
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
        ori s11, s11, 0b11111111111111111111000000000000# so the higher bits of s11 are ones
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
        ori s11, s11, 0b11111111111111111110000000000000# so the higher bits of s11 are ones
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

    li  s7, 0b1111110000000000000000000000000
    and s8, t5, s7 # get [30:25] -> imm[10:5]
    srli    s8, s8, 20 # real imm[10:5] is in [10:5]
    add s11, s11, s8 # add imm[10:5]

    li  s7, 0b111100000000
    and s8, t5, s7 # get [11:8] -> imm[4:1]
    srli    s8, s8, 7 # real imm[4:1] is in [4:1]
    add s11, s11, s8 # add imm[4:1]

    li  s7, 0b10000000
    and s8, t5, s7 # get [7] -> imm[11]
    slli    s8, s8, 4 # real imm[11] is in [11]
    add s11, s11, s8 # add imm[11]
    
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
        ori s11, s11, 0b11111111111111111110000000000000# so the higher bits of s11 are ones
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