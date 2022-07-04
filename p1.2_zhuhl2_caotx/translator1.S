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
    srli s8, s8, 12
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
    srli s10, s10, 12

    li s7, 0b0001000000000000
    and s8, t5, s7
    srli s8, s8, 12
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
    slli s10, s10, 10
    add s9, s9, s10
    srli s10, s10, 10

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
    srli s8, s8, 12
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
    add s9, s9, s11
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
    srli s8, s8, 10
    addi s11, s11, s8

    li  s7, 0b0000000000100000
    and s8, t5, s7
    srli s8, s8, 1
    addi s11, s11, s7

    # s10 <-- rd
    li  s7, 0b0000000000011100
    and s8, s8, s7
    srli s8, s8, 2
    # do mapping on register
    addi s8, s8, 8

    # s9 <-- rs1
    li  s7, 0b0000001110000000
    and s8, s8, s7
    srli s8, s8, 7
    addi s8, s8, 8

    # s8 <-- result
    li  s8, 0b00000000000000000010000000000011
    # add imm
    slli    s11, s11, 20
    add s8, s8, s11
    srli    s11, s11, 20
    # add rs1
    slli    s9, s9, 15
    add s8, s8, s9
    srli    s9, s9, 15
    # add rd
    slli    s10, s10, 7
    add s8, s8, s10
    srli    s10, s10, 7

    mv  a7, s10
    jr  ra


expand_sw:
    # s11 <-- offset
    li  s11, 0
    li  s7, 0b0000000001000000
    and s8, t5, s7
    srli s11, s8, 6

    li  s7, 0b0001110000000000
    and s8, t5, s7
    srli s8, s8, 10
    addi s11, s11, s8

    li  s7, 0b0000000000100000
    and s8, t5, s7
    srli s8, s8, 1
    addi s11, s11, s7

    # s10 <-- rs2
    li  s7, 0b0000000000011100
    and s8, s8, s7
    srli s8, s8, 2
    # do mapping on register
    addi s8, s8, 8

    # s9 <-- rs1
    li  s7, 0b0000001110000000
    and s8, s8, s7
    srli s8, s8, 7
    addi s8, s8, 8

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