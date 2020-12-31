/*
 * TMS320C28X disassembler
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "qemu/osdep.h"
#include "disas/dis-asm.h"
#include "qemu/bitops.h"
#include "cpu.h"

typedef disassemble_info DisasContext;

static void get_loc_string(char* str, uint32_t loc, uint32_t loc_type) {
    switch((loc & 0b11000000)>>6) {
        case 0b00: //b00 III III @6bit
            sprintf(str, "@0x%x", loc & 0x3f);
            break;
        case 0b01: //b01 III III *-SP[6bit]
            sprintf(str, "*-SP[0x%x]", loc & 0x3f);
            break;
        case 0b10: //b10 ... ...
            switch((loc & 0b00111000)>>3) {
                case 0b000: //b10 000 AAA *XARn++
                    sprintf(str, "*XAR%d++", loc & 0b00000111);
                    break;
                case 0b001: //b10 001 AAA *--XARn
                    sprintf(str, "*--XAR%d", loc & 0b00000111);
                    break;
                case 0b010: //b10 010 AAA *+XARn[AR0]
                    sprintf(str, "*+XAR%d[AR0]", loc & 0b00000111);
                    break;
                case 0b011: //b10 011 AAA *+XARn[AR1]
                    sprintf(str, "*+XAR%d[AR1]", loc & 0b00000111);
                    break;
                case 0b100: //@XARn @ARn
                    if (loc_type == LOC16) {
                        sprintf(str, "@AR%d", loc & 0b00000111);
                    }
                    else {
                        sprintf(str, "@XAR%d", loc & 0b00000111);
                    }
                    break;
                case 0b101: //b10 101 ...
                    switch(loc & 0b00000111) {
                        case 0b000: //@AH
                            sprintf(str, "@AH");
                            break;
                        case 0b001: //@ACC @AL
                            if (loc_type == LOC16) {
                                sprintf(str, "@AL");
                            }
                            else {
                                sprintf(str, "@ACC");
                            }
                            break;
                        case 0b010:
                            sprintf(str, "@PH");
                            break;
                        case 0b011:
                            if (loc_type == LOC16) {
                                sprintf(str, "@PL");
                            }
                            else {
                                sprintf(str, "@P");
                            }
                            break;
                        case 0b100:
                            if (loc_type == LOC16) {
                                sprintf(str, "@TH");
                            }
                            else {
                                sprintf(str, "@XT");
                            }
                            break;
                        case 0b101: //@SP
                            sprintf(str, "@SP");
                            break;
                        case 0b110: //b10 101 110 *BR0++
                            sprintf(str, "*BR0++");
                            break;
                        case 0b111: //b10 101 111 *BR0--
                            sprintf(str, "*BR0--");
                            break;
                    }
                    break;
                case 0b110: //b10 110 RRR *,ARPn
                    sprintf(str, "*,ARP%d",loc & 0b00000111);
                    break;
                case 0b111: //b10 111 ...
                    switch(loc & 0b00000111) {
                        case 0b000: //b10 111 000 *
                            sprintf(str, "*");
                            break;
                        case 0b001: //b10 111 001 *++
                            sprintf(str, "*++");
                            break;
                        case 0b010: //b10 111 010 *--
                            sprintf(str, "*--");
                            break;
                        case 0b011: //b10 111 011 *0++
                            sprintf(str, "*0++");
                            break;
                        case 0b100: //b10 111 100 *0--
                            sprintf(str, "*0--");
                            break;
                        case 0b101: //b10 111 101 *SP++
                            sprintf(str, "*SP++");
                            break;
                        case 0b110: //b10 111 110 *--SP
                            sprintf(str, "*--SP");
                            break;
                        case 0b111: //b10 111 111 *AR6%++
                            sprintf(str, "*AR6%%++");
                            break;
                    }
                    break;
            }
            break;
        case 0b11: //b11 III AAA
            sprintf(str, "*+XAR%d[0x%x]", loc & 0b00000111, (loc >> 3) & 0b00000111);
            break;
    }
}

static void get_cond_string(char *str, uint32_t cond) {
    cond = cond & 0xf;
    switch(cond) {
        case 0: //NEQ
            sprintf(str, "NEQ");
            break;
        case 1: //EQ
            sprintf(str, "EQ");
            break;
        case 2: //GT
            sprintf(str, "GT");
            break;
        case 3: //GEQ
            sprintf(str, "GEQ");
            break;
        case 4: //LT
            sprintf(str, "LT");
            break;
        case 5: //LEQ
            sprintf(str, "LEQ");
            break;
        case 6: //HI
            sprintf(str, "HI");
            break;
        case 7: //HIS,C
            sprintf(str, "HIS(C)");
            break;
        case 8: //LO,NC
            sprintf(str, "LO(NC)");
            break;
        case 9: //LOS
            sprintf(str, "LOS");
            break;
        case 10: //NOV
            sprintf(str, "NOV");
            break;
        case 11: //OV
            sprintf(str, "OV");
            break;
        case 12: //NTC
            sprintf(str, "NTC");
            break;
        case 13: //TC
            sprintf(str, "TC");
            break;
        case 14: //NBIO
            sprintf(str, "NBIO");
            break;
        case 15: //UNC
            sprintf(str, "UNC");
            break;
        default://can not reach
            break;
    }
}

static void get_stf_bit_string(char *str, uint32_t mode) {
    mode = mode & 0xff;
    uint32_t i = 0;
    if (((mode>>0) & 1) == 1) {
        sprintf(str+i, "LVF,");
        i += 4;
    }
    if (((mode>>1) & 1) == 1) {
        sprintf(str+i, "LUF,");
        i += 4;
    }
    if (((mode>>2) & 1) == 1) {
        sprintf(str+i, "NF,");
        i += 3;
    }
    if (((mode>>3) & 1) == 1) {
        sprintf(str+i, "NI,");
        i += 3;
    }
    if (((mode>>4) & 1) == 1) {
        sprintf(str+i, "ZF,");
        i += 3;
    }
    if (((mode>>5) & 1) == 1) {
        sprintf(str+i, "ZI,");
        i += 3;
    }
    if (((mode>>6) & 1) == 1) {
        sprintf(str+i, "CI,");
        i += 3;
    }
    if (((mode>>7) & 1) == 1) {
        sprintf(str+i, "TF,");
        i += 3;
    }
    str[i-1] = 0;
}

static void get_status_bit_string(char *str, uint32_t mode) {
    mode = mode & 0xff;
    uint32_t i = 0;
    if (((mode>>0) & 1) == 1) {
        sprintf(str+i, "SXM,");
        i += 4;
    }
    if (((mode>>1) & 1) == 1) {
        sprintf(str+i, "OVM,");
        i += 4;
    }
    if (((mode>>2) & 1) == 1) {
        sprintf(str+i, "TC,");
        i += 3;
    }
    if (((mode>>3) & 1) == 1) {
        sprintf(str+i, "C,");
        i += 2;
    }
    if (((mode>>4) & 1) == 1) {
        sprintf(str+i, "INTM,");
        i += 5;
    }
    if (((mode>>5) & 1) == 1) {
        sprintf(str+i, "DBGM,");
        i += 5;
    }
    if (((mode>>6) & 1) == 1) {
        sprintf(str+i, "PAGE0,");
        i += 6;
    }
    if (((mode>>7) & 1) == 1) {
        sprintf(str+i, "VMAP,");
        i += 5;
    }
    str[i-1] = 0;
}

static bool get_fpu_reg_name(char *str, uint32_t addr)
{
    switch (addr)
    {
        case 0xf00:
        {
            sprintf(str, "RB");
            return true;
        }
        case 0xf02:
        {
            sprintf(str, "STF");
            return true;
        }
        case 0xf12:
        {
            sprintf(str, "R0H");
            return true;
        }
        case 0xf16:
        {
            sprintf(str, "R1H");
            return true;
        }
        case 0xf1A:
        {
            sprintf(str, "R2H");
            return true;
        }
        case 0xf1E:
        {
            sprintf(str, "R3H");
            return true;
        }
        case 0xf22:
        {
            sprintf(str, "R4H");
            return true;
        }
        case 0xf26:
        {
            sprintf(str, "R5H");
            return true;
        }
        case 0xf2a:
        {
            sprintf(str, "R6H");
            return true;
        }
        case 0xf2e:
        {
            sprintf(str, "R7H");
            return true;
        }
        default:
        {
            return false;
        }
    }
}

static void get_condf_string(char *str, uint32_t cond) {
    cond = cond & 0xf;
    switch(cond) {
        case 0: //NEQ
            sprintf(str, "NEQ");
            break;
        case 1: //EQ
            sprintf(str, "EQ");
            break;
        case 2: //GT
            sprintf(str, "GT");
            break;
        case 3: //GEQ
            sprintf(str, "GEQ");
            break;
        case 4: //LT
            sprintf(str, "LT");
            break;
        case 5: //LEQ
            sprintf(str, "LEQ");
            break;
        case 10: //TF
            sprintf(str, "TF");
            break;
        case 11: //NTF
            sprintf(str, "NTF");
            break;
        case 12: //LU
            sprintf(str, "LU");
            break;
        case 13: //LV
            sprintf(str, "LV");
            break;
        case 14: //UNC
            sprintf(str, "UNC");
            break;
        case 15: //UNCF
            sprintf(str, "UNCF");
            break;
        default://can not reach
            break;
    }
}

int print_insn_tms320c28x(bfd_vma addr, disassemble_info *info)
{
    fprintf_function fprintf_func = info->fprintf_func;
    void *stream = info->stream;
    bfd_byte buffer[4];
    uint32_t insn,insn32;
    int status;
    int length = 2;
    status = info->read_memory_func(addr, buffer, 4, info);
    if (status != 0) {
        info->memory_error_func(status, addr, info);
        return -1;
    }
    char str[40];
    char str2[40];

    unsigned long v;
    v = (unsigned long) buffer[2];
    v |= (unsigned long) buffer[3] << 8;
    v |= (unsigned long) buffer[0] << 16;
    v |= (unsigned long) buffer[1] << 24;
    insn32 = v;
    insn = insn32 >> 16;

    // fprintf_func (stream,  ".long, 0x%04x\t", insn);

    switch ((insn & 0xf000) >> 12) {
        case 0b0000:
            switch ((insn & 0x0f00) >> 8) {
                case 0b0000: //0000 0000 .... ....
                    switch ((insn & 0x00c0) >> 6) {
                        case 0b00:
                            if(((insn >> 5) & 1) == 0) { //0000 0000 000. ....
                                if(((insn >> 4) & 1) == 0) { //0000 0000 0000 ....
                                    switch(insn & 0x000f) {
                                        case 0: //0000 0000 0000 0000
                                            fprintf_func(stream, "0x%04x;     ITRAP0", insn);
                                            break;
                                        case 1: //0000 0000 0000 0001, ABORTI P124
                                            fprintf_func(stream, "0x%04x;     ABORTI", insn);
                                            break;
                                        case 2: //0000 0000 0000 0010, POP IFR
                                            fprintf_func(stream, "0x%04x;     POP IFR", insn);
                                            break;
                                        case 3: //0000 0000 0000 0011 POP AR1H:AR0H
                                            fprintf_func(stream, "0x%04x;     POP AR1H:AR0H", insn);
                                            break;
                                        case 4: //0000 0000 0000 0100 PUSH RPC
                                            fprintf_func(stream, "0x%04x;     PUSH RPC", insn);
                                            break;
                                        case 5: //0000 0000 0000 0101 PUSH AR1H:AR0H
                                            fprintf_func(stream, "0x%04x;     PUSH AR1H:AR0H", insn);
                                            break;
                                        case 6: //0000 0000 0000 0110 LRETR
                                            fprintf_func(stream, "0x%04x;     LRETR", insn);
                                            break;
                                        case 7: //0000 0000 0000 0111 POP RPC
                                            fprintf_func(stream, "0x%04x;     POP RPC", insn);
                                            break;
                                        default: //0000 0000 0000 1nnn CCCC CCCC CCCC CCCC BANZ 16bitOffset,ARn--
                                        {
                                            uint32_t n = insn & 0b111;
                                            int16_t offset = (uint16_t)(insn32&0xffff);
                                            length = 4;
                                            fprintf_func(stream, "0x%08x; BANZ #%d,AR%d--", insn32, offset, n);
                                            break;
                                        }
                                    }
                                }
                                else { //0000 0000 0001 CCCC INTR INTx
                                    uint32_t n = insn & 0xf;
                                    fprintf_func(stream, "0x%04x;     INTR %s", insn, INTERRUPT_NAME[n]);
                                }
                            }
                            else {// 0000 0000 001C CCCC, TRAP #VectorNumber
                                uint32_t n = insn & 0b11111;
                                fprintf_func(stream, "0x%04x;     TRAP %s", insn, INTERRUPT_NAME[n]);
                            }
                            break;
                        case 0b01: /*0000 0000 01.. .... LB 22bit */
                        {
                            uint32_t dest = insn32 & 0x3fffff;
                            fprintf_func(stream, "0x%08x; LB @0x%x", insn32, dest);
                            length = 4;
                            break;
                        }
                        case 0b10: // 0000 0000 10CC CCCC CCCC CCCC CCCC CCCC, LC 22bit
                        {
                            uint32_t imm = insn32 & 0x3fffff;
                            fprintf_func(stream, "0x%08x; LC @0x%x", insn32, imm);
                            length = 4;
                            break;
                        }
                        case 0b11: // 0000 0000 11CC CCCC CCCC CCCC CCCC CCCC, FFC XAR7,22bit
                        {
                            uint32_t imm = insn32 & 0x3fffff;
                            fprintf_func(stream, "0x%08x; FFC XAR7,@0x%x", insn32, imm);
                            length = 4;
                            break;
                        }
                    }
                    break;
                case 0b0001: //0000 0001 LLLL LLLL SUBU ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     SUBU ACC,%s", insn, str);
                    break;
                }
                case 0b0010: //0000 0010 CCCC CCCC MOVB ACC,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     MOVB ACC,#%d", insn, imm);
                    break;
                }
                case 0b0011: //0000 0011 LLLL LLLL SUBL ACC,loc32
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC32);
                    fprintf_func(stream, "0x%04x;     SUBL ACC,%s", insn, str);
                    break;
                }
                case 0b0100:
                    break;
                case 0b0101: //0000 0101 LLLL LLLL ADD ACC,loc16<<#16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     ADD ACC,%s<<#16", insn, str);
                    break;
                }
                case 0b0110: //0000 0110 .... .... MOVL ACC,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0xbe)
                    {
                        fprintf_func(stream, "0x%04x;     POP ACC", insn);
                    }
                    else
                    {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL ACC,%s", insn, str);
                    }
                    break;
                }
                case 0b0111: //0000 0111 LLLL LLLL ADDL ACC,loc32
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC32);
                    fprintf_func(stream, "0x%04x;     ADDL ACC,%s", insn, str);
                    break;
                }
                case 0b1000: //0000 1000 LLLL LLLL 32bit ADD loc16,#16bitSigned
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    uint32_t imm = insn32 & 0xffff;
                    if ((imm & 0x8000) >> 15) { //neg value
                        imm = imm | 0xffff0000;
                    }
                    length = 4;
                    fprintf_func(stream, "0x%04x;     ADD %s, #%d", insn, str, imm);
                    break;
                }
                case 0b1001: //0000 1001 CCCC CCCC ADDB ACC,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     ADDB ACC,#%d", insn, imm);
                    break;
                }
                case 0b1010: //0000 1010 LLLL LLLL INC loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     INC %s", insn, str);
                    break;
                }
                case 0b1011: //0000 1011 LLLL LLLL DEC loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     DEC %s", insn, str);
                    break;
                }
                case 0b1100: //0000 1100 LLLL LLLL ADDCU ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     ADDCU ACC,%s", insn, str);
                    break;
                }
                case 0b1101: //0000 1101 LLLL LLLL ADDU ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     ADDU ACC,%s", insn, str);
                    break;
                }
                case 0b1110://0000 1110 LLLL LLLL MOVU ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVU ACC,%s", insn, str);
                    break;
                }
                case 0b1111://0000 1111 LLLL LLLL CMPL ACC,loc32
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC32);
                    fprintf_func(stream, "0x%04x;     CMPL ACC,%s", insn, str);
                    break;
                }
            }
            break;
        case 0b0001:
            switch ((insn & 0x0f00) >> 8) {
                case 0b0000: //0001 0000 LLLL LLLL MOVA T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10101100) {
                        fprintf_func(stream, "0x%04x;     ADDL ACC,P<<PM", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC16);
                        fprintf_func(stream, "0x%04x;     MOVA T,%s", insn, str);
                    }
                    break;
                }
                case 0b0001: //0001 0001 LLLL LLLL MOVS T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10101100)
                    {
                        fprintf_func(stream, "0x%04x;     SUBL ACC,P<<PM", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC16);
                        fprintf_func(stream, "0x%04x;     MOVS T,%s", insn, str);
                    }
                    break;
                }
                case 0b0010: //0001 0010 LLLL LLLL MPY ACC,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPY ACC,T,%s", insn, str);
                    break;
                }
                case 0b0011: //0001 0011 LLLL LLLL MPYS P,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPYS P,T,%s", insn, str);
                    break;
                }
                case 0b0100: //0001 0100 LLLL LLLL CCCC CCCC CCCC CCCC MAC P,loc16,0:pma
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t addr = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; MAC P,%s,0:#0x%04x", insn, str, addr);
                    length = 4;
                    break;
                }
                case 0b0101: //0001 0101 LLLL LLLL CCCC CCCC CCCC CCCC MPYA P,loc16,#16bit
                {
                    uint32_t mode = insn & 0xff;
                    int32_t imm = (int16_t)(insn32 & 0xffff);
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; MPYA P,%s,#%d", insn, str, imm);
                    length = 4;
                    break;
                }
                case 0b0110: //0001 0110 LLLL LLLL MOVP T,loc16
                {

                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10101100)
                    {
                        fprintf_func(stream, "0x%04x;     MOVL ACC,P<<PM", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC16);
                        fprintf_func(stream, "0x%04x;     MOVP T,%s", insn, str);
                    }
                    break;
                }
                case 0b0111: //0001 0111 LLLL LLLL MPYA P,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPYA P,T,%s", insn, str);
                    break;
                }
                case 0b1000: //0001 1000 LLLL LLLL CCCC CCCC CCCC CCCC AND loc16,#16bitSigned
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t imm = insn32 & 0xffff;
                    length = 4;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; AND %s,#0x%x", insn, str, imm);
                    break;
                }
                case 0b1001: //0001 1001 CCCC CCCC SUBB ACC,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     SUBB ACC,#0x%x", insn, imm);
                    break;
                }
                case 0b1010: //0001 1010 LLLL LLLL CCCC CCCC CCCC CCCC OR loc16,#16bit
                {
                    length = 4;
                    int32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; OR %s,#0x%x", insn32, str, imm);
                    break;
                }
                case 0b1011: //0001 1011 LLLL LLLL CCCC CCCC CCCC CCCC CMP loc16,#16bitSigned
                {
                    length = 4;
                    int32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; CMP %s,#0x%x", insn32, str, imm);
                    break;
                }
                case 0b1100: //0001 1100 LLLL LLLL CCCC CCCC CCCC CCCC XOR loc16,#16bit
                {
                    length = 4;
                    int32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; XOR %s,#0x%x", insn32, str, imm);
                    break;
                }
                case 0b1101: //0001 1101 LLLL LLLL SBBU ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     SBBU ACC,%s", insn, str);
                    break;
                }
                case 0b1110: //0001 1110 .... .... MOVL loc32,ACC
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH ACC", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOV %s,ACC", insn, str);
                    }
                    break;
                }
                case 0b1111: //0001 1111 LLLL LLLL SUBCU ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     SUBCU ACC,%s", insn, str);
                    break;
                }
            }
            break;
        case 0b0010:
            switch ((insn & 0x0f00) >> 8) {
                case 0b0000: //0010 0000 LLLL LLLL MOV loc16,IER
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV %s,IER", insn, str);
                    break;
                }
                case 0b0001: //0010 0001 LLLL LLLL MOV loc16,T
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV %s,P", insn, str);
                    break;
                }
                case 0b0010: //0010 0010 LLLL LLLL PUSH loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     PUSH %s", insn, str);
                    break;
                }
                case 0b0011: //0010 0011 LLLL LLLL MOV IER,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV IER,%s", insn, str);
                    break;
                }
                case 0b0100: //0010 0100 LLLL LLLL PREAD loc16,*XAR7
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     PREAD %s,*XAR7", insn, str);
                    break;
                }
                case 0b0101: /* 0010 0101 .... .... MOV ACC, loc16<<#16 */
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV ACC,%s<<#16", insn, str);
                }
                    break;
                case 0b0110: //0010 0110 LLLL LLLL PWRITE *XAR7,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     PWRITE *XAR7,%s", insn, str);
                    break;
                }
                case 0b0111: //0010 1111 LLLL LLLL MOV PL,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV PL,%s", insn, str);
                    break;
                }
                case 0b1000: /* MOV loc16, #16bit p260 */
                {
                    uint32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; MOV %s,#0x%x", insn32, str,imm);
                    length = 4;
                    break;
                }
                case 0b1001: //0010 1001 CCCC CCCC CLRC mode
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10000)
                    {
                        fprintf_func(stream, "0x%04x;     EINT", insn);
                    }
                    else
                    {
                        get_status_bit_string(str, mode);
                        fprintf_func(stream, "0x%04x;     CLRC %s", insn, str);
                    }
                    break;
                }
                case 0b1010: //0010 1010 LLLL LLLL POP loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     POP %s", insn, str);
                    break;
                }
                case 0b1011: //0010 1011 LLLL LLLL MOV loc16,#0
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV %s,#0", insn, str);
                    break;
                }
                case 0b1100: //0010 1100 LLLL LLLL CCCC CCCC CCCC CCCC LOOPZ loc16,#16bit
                {
                    length = 4;
                    uint32_t mask = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; LOOPZ %s,#0x%d", insn, str, mask);
                    break;
                }
                case 0b1101: //0010 1101 LLLL LLLL MOV T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV T,%s", insn, str);
                    break;
                }
                case 0b1110: //0010 1110 LLLL LLLL CCCC CCCC CCCC CCCC LOOPNZ loc16,#16bit
                {
                    length = 4;
                    uint32_t mask = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; LOOPNZ %s,#0x%d", insn, str, mask);
                    break;
                }
                case 0b1111: //0010 1111 LLLL LLLL MOV PH,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV PH,%s", insn, str);
                    break;
                }
            }
            break;
        case 0b0011:
            switch ((insn & 0xf00) >> 8) {
                case 0b0000: //0011 0000 LLLL LLLL MPYXU ACC,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPYXU ACC,T,%s", insn, str);
                    break;
                }
                case 0b0001: //0011 0001 CCCC CCCC MPYB P,T,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     MPYB P,T,#%d", insn, imm);
                    break;
                }
                case 0b0010: //0011 0010 LLLL LLLL MPYXU P,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPYXU P,T,%s", insn, str);
                    break;
                }
                case 0b0011: //0011 0011 LLLL LLLL MPY P,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPY P,T,%s", insn, str);
                    break;
                }
                case 0b0100: //0011 0100 LLLL LLLL CCCC CCCC CCCC CCCC MPY ACC,loc16,#16bit
                {
                    uint32_t mode = insn & 0xff;
                    int32_t imm = (int16_t)(insn32 & 0xffff);
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; MPY ACC,%s,#%d", insn, str, imm);
                    length = 4;
                    break;
                }
                case 0b0101: //0011 0101 CCCC CCCC MPYB ACC,T,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     MPYB ACC,T,#%d", insn, imm);
                    break;
                }
                case 0b0110: //0011 0111 LLLL LLLL MPYU ACC,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPYU ACC,T,%s", insn, str);
                    break;
                }
                case 0b0111: //0011 0111 LLLL LLLL MPYU P,T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MPYU P,T,%s", insn, str);
                    break;
                }
                case 0b1000: //0011 1000 LLLL LLLL MOVB AL.MSB,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB AL.MSB,%s", insn, str);
                    break;
                }
                case 0b1001: //0011 1001 LLLL LLLL MOVB AH.MSB,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB AL.MSB,%s", insn, str);
                    break;
                }
                case 0b1010: //0011 1010 LLLL LLLL MOVL loc32, XAR0
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR0", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR0", insn, str);
                    }
                    break;
                }
                case 0b1011: //0011 1011 CCCC CCCC SETC mode
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b00010000) {
                        fprintf_func(stream, "0x%04x;     DINT", insn);
                    }
                    else {
                        get_status_bit_string(str, mode);
                        fprintf_func(stream, "0x%04x;     SETC %s", insn, str);
                    }
                    break;
                }
                case 0b1100://0011 110a LLLL LLLL MOVB loc16,AL.LSB
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB %s,AL.LSB", insn, str);
                    break;
                }
                case 0b1101://0011 110a LLLL LLLL MOVB loc16,AH.LSB
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB %s,AH.LSB", insn, str);
                    break;
                }
                case 0b1110: //0011 1110 .... ....
                {
                    switch ((insn & 0x00f0) >> 4) {
                        case 0b0000: //0011 1110 0000 SHFT CCCC CCCC CCCC CCCC AND ACC,#16bit<<0...15
                        {
                            uint32_t imm = insn32 & 0xffff;
                            uint32_t shift = insn & 0xf;
                            fprintf_func(stream, "0x%08x; AND ACC,#%d<<#%d", insn32, imm, shift);
                            length = 4;
                            break;
                        }
                        case 0b0001: //0011 1110 0001 SHFT CCCC CCCC CCCC CCCC OR ACC,#16bit<<#0...15
                        {
                            uint32_t imm = insn32 & 0xffff;
                            uint32_t shift = insn & 0xf;
                            fprintf_func(stream, "0x%08x; OR ACC,#%d<<#%d", insn32, imm, shift);
                            length = 4;
                            break;
                        }
                        case 0b0010: //0011 1110 0010 SHFT CCCC CCCC CCCC CCCC XOR ACC,#16bit<<#0...15
                        {
                            length = 4;
                            uint32_t shift = insn & 0xf;
                            uint32_t imm  = insn32 & 0xffff;
                            fprintf_func(stream, "0x%08x; XOR ACC,#%d<<#%d", insn32, imm, shift);
                            break;
                        }
                        case 0b0011: //0011 1110 0011 ....
                        {
                            uint32_t addr = insn32 & 0xffff;
                            uint32_t n = insn & 0b111;
                            if (((insn >> 3) & 1) == 0) //0011 1110 0011 0nnn CCCC CCCC CCCC CCCC XBANZ pma,*,ARPn
                            {
                                fprintf_func(stream, "0x%08x; XBANZ 0x%x,*,ARP%d", insn32, addr,n);
                            }
                            else // XBANZ pma,*++,ARPn
                            {
                                fprintf_func(stream, "0x%08x; XBANZ 0x%x,*++,ARP%d", insn32, addr,n);
                            }
                            length = 4;
                            break;
                        }
                        case 0b0100: //0011 1110 0100 ....
                        {
                            uint32_t addr = insn32 & 0xffff;
                            uint32_t n = insn & 0b111;
                            if (((insn >> 3) & 1) == 0) //0011 1110 0100 0nnn CCCC CCCC CCCC CCCC XBANZ pma,*--,ARPn
                            {
                                fprintf_func(stream, "0x%08x; XBANZ 0x%x,*--,ARP%d", insn32, addr,n);
                            }
                            else // XBANZ pma,*0++,ARPn
                            {
                                fprintf_func(stream, "0x%08x; XBANZ 0x%x,*0++,ARP%d", insn32, addr,n);
                            }
                            length = 4;
                            break;
                        }
                        case 0b0101://0011 1110 0101 ....
                        {
                            if (((insn >> 3) & 1) == 1) {//0011 1110 0101 1nnn MOV XARn,PC
                                uint32_t n = insn & 0b111;
                                fprintf_func(stream, "0x%04x;     MOV XAR%d,PC", insn, n);
                            }
                            else {//0011 1110 0101 0nnn CCCC CCCC CCCC CCCC XBANZ pma,*0--,ARPn
                                length = 4;
                                uint32_t addr = insn32 & 0xffff;
                                uint32_t n = insn & 0b111;
                                fprintf_func(stream, "0x%08x; XBANZ 0x%x,*0--,ARP%d", insn32, addr,n);
                            }
                            break;
                        }
                        case 0b0110://0011 1110 0110 ....
                        {
                            if (((insn >> 3) & 1) == 0) //0011 1110 0110 0RRR LCR *XARn
                            {
                                uint32_t n = insn & 0b111;
                                fprintf_func(stream, "0x%04x;     LCR *XAR%d", insn, n);
                            }
                            else //0011 1110 0110 1nnn CCCC CCCC CCCC CCCC XCALL pma,*,ARPn
                            {
                                length = 4;
                                uint32_t addr = insn32 & 0xffff;
                                uint32_t n = insn & 0b111;
                                fprintf_func(stream, "0x%08x; XCALL 0x%x,*,ARP%d", insn32, addr,n);
                            }
                            break;
                        }
                        case 0b0111: //0011 1110 0111 ....
                        {
                            if (((insn >> 3) & 1) == 0) //0011 1110 0111 0nnn CCCC CCCC CCCC CCCC XB pma,*,ARPn
                            {
                                uint32_t n = insn & 0b111;
                                uint32_t addr = insn32 & 0xffff;
                                fprintf_func(stream, "0x%08x; LCR 0x%x,*,ARP%d", insn, addr, n);
                                length = 4;
                            }
                            else
                            {

                            }
                            break;
                        }
                    }
                    break;
                }
                case 0b1111: //0011 1111 LLLL LLLL MOV loc16,P
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOV %s,P", insn, str);
                    break;
                }
            }
            break;
        case 0b0100: //0100 BBBB LLLL LLLL TBIT loc16,#bit
        {
            uint32_t mode = insn & 0xff;
            uint32_t bit = (insn >> 8) & 0xf;
            get_loc_string(str, mode, LOC16);
            fprintf_func(stream, "0x%04x;     TBIT %s,#%d", insn, str, bit);
            break;
        }
        case 0b0101:
            switch ((insn & 0x0f00) >> 8) {
                case 0b0000://0101 0000 CCCC CCCC ORB AL,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     ORB AL,#0x%x", insn, imm);
                    break;
                }
                case 0b0001://0101 0001 CCCC CCCC ORB AH,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     ORB AH,#0x%x", insn, imm);
                    break;
                }
                case 0b0010:
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     CMP AL,#0x%x", insn, imm);
                    break;
                }
                case 0b0011://0101 001A CCCC CCCC CMPB AX,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     CMP AH,#0x%x", insn, imm);
                    break;
                }
                case 0b0100:
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     CMP AL,%s", insn, str);
                }
                case 0b0101://0101 010A LLLL LLLL CMP AX,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     CMP AH,%s", insn, str);
                }
                case 0b0110://0101 0110 .... ....
                    switch ((insn & 0x00f0) >> 4) {
                        case 0b0000: //0101 0110 0000 ....
                        {
                            switch (insn & 0x000f) {
                                case 0b0000: //0101 0110 0000 0000 0000 SHFT LLLL LLLL SUB ACC,loc16<<#1...15
                                {
                                    if (((insn32 & 0xf000) >> 12) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        uint32_t shift = (insn32 >> 8) & 0xf;
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; SUB ACC,%s<<#%d", insn32, str, shift);
                                    }
                                    break;
                                }
                                case 0b0001: //0101 0110 0000 0001 0000 0000 LLLL LLLL ADDL loc32,ACC
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str,mode,LOC32);
                                        fprintf_func(stream, "0x%08x; ADDL %s,ACC", insn32, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0010: //0101 0110 0000 0010 0000 0000 LLLL LLLL MOV OVC,loc16
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; MOV OVC,%s", insn32, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0011: /* 0101 0110 0000 0011  MOV ACC, loc16<<#1...15 */
                                {
                                    if ((insn32 & 0xf000) == 0) { //this 4bit should be 0
                                        uint32_t shift = (insn32 & 0x0f00) >> 8;
                                        uint32_t mode = (insn32 & 0xff);
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; MOV ACC,%s<<#%d", insn32, str,shift);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0100: // 0101 0110 0000 0100 ADD ACC,loc16<<#1...15
                                {
                                    if ((insn32 & 0xf000) == 0) { //this 4bit should be 0
                                        uint32_t shift = (insn32 & 0x0f00) >> 8;
                                        uint32_t mode = (insn32 & 0xff);
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; ADD ACC,%s<<#%d", insn32, str,shift);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0101: //0101 0110 0000 0101 0000 0000 LLLL LLLL IMPYL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str,mode,LOC32);
                                        fprintf_func(stream, "0x%08x; IMPYL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0110: //0101 0110 0000 0110 0000 0000 LLLL LLLL MOV ACC,loc16<<T
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; MOV ACC,%s<<T", insn32, str);
                                    }
                                    length = 4;
                                    break;
                                }
                                case 0b0111: //0101 0110 0000 0111 .... //MAC P,loc16,*XAR7/++
                                {
                                    uint32_t mode2 = (insn32 >> 8) & 0xff;
                                    uint32_t mode1 = insn32 & 0xff;
                                    if (mode2 == 0b11000111)
                                    {
                                        length = 4;
                                        get_loc_string(str, mode1, LOC16);
                                        fprintf_func(stream, "0x%08x; MAC P,%s,*XAR7", insn32, str);
                                    }
                                    else if (mode2 == 0b10000111)
                                    {
                                        length = 4;
                                        get_loc_string(str, mode1, LOC16);
                                        fprintf_func(stream, "0x%08x; MAC P,%s,*XAR7++", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1000: //0101 0110 0000 1000 CCCC CCCC CCCC CCCC AND ACC,#16bit<<#16
                                {
                                    uint32_t imm = insn32 & 0xff;
                                    fprintf_func(stream, "0x%08x; AND ACC,#%d<<#16", insn32, imm);
                                    length = 4;
                                    break;
                                }
                                case 0b1001: //0101 0110 0000 1001 0000 BBBB LLLL LLLL TCLR loc16,#bit
                                {
                                    if (((insn32 >> 12) & 0xf) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        uint32_t bit = (insn32 >> 8) & 0xf;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; TCLR %s,#%d", insn32, str, bit);
                                    }
                                    break;
                                }
                                case 0b1010: //0101 0110 0000 1010 CCCC CCCC CCCC CCCC XBANZ pma,*++
                                {
                                    length = 4;
                                    uint32_t addr = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; XBANZ 0x%x,*++", insn32, addr);
                                    break;
                                }
                                case 0b1011: //0101 0110 0000 1010 CCCC CCCC CCCC CCCC XBANZ pma,*--
                                {
                                    length = 4;
                                    uint32_t addr = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; XBANZ 0x%x,*--", insn32, addr);
                                    break;
                                }
                                case 0b1100: //0101 0110 0000 1100 CCCC CCCC CCCC CCCC XBANZ pma,*
                                {
                                    length = 4;
                                    uint32_t addr = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; XBANZ 0x%x,*", insn32, addr);
                                    break;
                                }
                                case 0b1101: //0101 0110 0000 1101 0000 BBBB LLLL LLLL TSET loc16,#bit
                                {
                                    if (((insn32 >> 12) & 0xf) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        uint32_t bit = (insn32 >> 8) & 0xf;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; TSET %s,#%d", insn32, str, bit);
                                    }
                                    break;
                                }
                                case 0b1110: //0101 0110 0000 1110 CCCC CCCC CCCC CCCC XBANZ pma,*0++
                                {
                                    length = 4;
                                    uint32_t addr = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; XBANZ 0x%x,*0++", insn32, addr);
                                    break;
                                }
                                case 0b1111: //0101 0110 0000 1111 CCCC CCCC CCCC CCCC XBANZ pma,*0--
                                {
                                    length = 4;
                                    uint32_t addr = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; XBANZ 0x%x,*0--", insn32, addr);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0001: //0101 0110 0001 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //0101 0110 0001 0000 ASRL ACC,T
                                {
                                    fprintf_func(stream, "0x%04x;     ASRL ACC:P,T", insn);
                                    break;
                                }
                                case 0b0001: //0101 0110 0001 0001 xxxx xxxx LLLL LLLL SQRS loc16
                                {
                                    uint32_t mode = insn32 & 0xff;
                                    length = 4;
                                    get_loc_string(str, mode, LOC16);
                                    fprintf_func(stream, "0x%08x; SQRS %s", insn32,str);
                                    break;
                                }
                                case 0b0011: //0101 0110 0001 0011 0000 0000 LLLL LLLL ZALR ACC,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; ZALR ACC,%s", insn32,str);
                                    }
                                    break;
                                }
                                case 0b0100: //0101 0110 0001 0100 XB *AL
                                {
                                    fprintf_func(stream, "0x%04x;     XB *AL", insn);
                                    break;
                                }
                                case 0b0101: //0101 0110 0001 0101 0000 0000 LLLL LLLL SQRA loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; SQRA %s", insn32,str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0110: //0101 0110 0001 0110 CLRC AMODE
                                {
                                    fprintf_func(stream, "0x%04x;     CLRC AMODE", insn);
                                    break;
                                }
                                case 0b0111: //0101 0110 0001 0111 0000 0000 LLLL LLLL SUBCUL ACC,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; SUBCUL ACC,%s", insn32,str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b1000: //0101 0110 0001 1000 CMPR 2
                                {
                                    fprintf_func(stream, "0x%04x;     CMPR 2", insn);
                                    break;
                                }
                                case 0b1001: //0101 0110 0001 1001 CMPR 1
                                {
                                    fprintf_func(stream, "0x%04x;     CMPR 1", insn);
                                    break;
                                }
                                case 0b1010: //0101 0110 0001 1010 SETC M0M1MAP
                                {
                                    fprintf_func(stream, "0x%04x;     SETC M0M1MAP", insn);
                                    break;
                                }
                                case 0b1011: //0101 0110 0001 1011 CLRC XF
                                {
                                    fprintf_func(stream, "0x%04x;     CLRC XF", insn);
                                    break;
                                }
                                case 0b1100: //0101 0110 0001 1100 CMPR 3
                                {
                                    fprintf_func(stream, "0x%04x;     CMPR 3", insn);
                                    break;
                                }
                                case 0b1101: //0101 0110 0001 1101 CMPR 0
                                {
                                    fprintf_func(stream, "0x%04x;     CMPR 0", insn);
                                    break;
                                }
                                case 0b1110: //0101 0110 0001 1110 LPADDR
                                {
                                    fprintf_func(stream, "0x%04x;     LPADDR", insn);
                                    break;
                                }
                                case 0b1111: //0101 0110 0001 1111 SETC Objmode
                                {
                                    fprintf_func(stream, "0x%04x;     SETC Objmode", insn);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0010: //0101 0110 0010 ....
                            switch (insn & 0xf) {
                                case 0b0000: //0101 0110 0010 0000 NORM ACC,*--
                                {
                                    fprintf_func(stream, "0x%04x;     NORM ACC,*--", insn);
                                    break;
                                }
                                case 0b0001: //0101 0110 0010 0001 xxxx xxxx LLLL LLLL MOVX TL,loc16
                                {
                                    length = 4;
                                    uint32_t mode = insn32 & 0xff;
                                    get_loc_string(str,mode,LOC16);
                                    fprintf_func(stream, "0x%08x; MOVX TL,%s", insn32, str);
                                    break;
                                }
                                case 0b0010: //0101 0110 0010 0010 LSRL ACC,T
                                {
                                    fprintf_func(stream, "0x%04x;     LSRL ACC,T", insn);
                                    break;
                                }
                                case 0b0011: //0101 0110 0010 0011 32bit ADD ACC,loc16<<T
                                {
                                    length = 4;
                                    if ((insn32 & 0xff00) == 0) { //this 8bit should be 0
                                        uint32_t loc16 = (insn32 & 0xff);
                                        get_loc_string(str,loc16,LOC16);
                                        fprintf_func(stream, "0x%08x; ADD ACC,%s<<T", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0100: //0101 0110 0010 0100 NORM ACC,*
                                {
                                    fprintf_func(stream, "0x%04x;     NORM ACC,*", insn);
                                    break;
                                }
                                case 0b0101: //0101 0110 0010 0101 0000 0000 LLLL LLLL TBIT loc16,T
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; TBIT %s,T", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0110: //0101 0110 0010 0110 SETC XF
                                {
                                    fprintf_func(stream, "0x%04x;     SETC XF", insn);
                                    break;
                                }
                                case 0b0111: //0101 0110 0010 0111 0000 0000 LLLL LLLL SUB ACC,loc16<<T
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; SUB ACC,%s<<T", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1000: //0101 0110 0010 1000 0000 0000 LLLL LLLL MOVU loc16,OVC
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MOVU %s,OVC", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1001: //0101 0110 0010 1001 0000 0000 LLLL LLLL MOV loc16,OVC
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MOV %s,OVC", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1010:
                                case 0b1011: //0101 0110 0010 101A 0000 COND LLLL LLLL MOV loc16,AX,COND
                                {
                                    uint32_t mode = insn32 & 0xff;
                                    uint32_t cond = (insn32 >> 8) & 0xf;
                                    if (((insn32 >> 12) & 0xf) == 0) {
                                        get_loc_string(str, mode, LOC16);
                                        get_cond_string(str2, cond);
                                        uint32_t is_AH = insn & 1;
                                        length = 4;
                                        if (is_AH) {
                                            fprintf_func(stream, "0x%08x; MOV %s,AH,%s", insn32, str, str2);
                                        }
                                        else {
                                            fprintf_func(stream, "0x%08x; MOV %s,AL,%s", insn32, str, str2);
                                        }
                                    }
                                    break;
                                }
                                case 0b1100: //0101 0110 0010 1100 ASR64 ACC:P,T
                                {
                                    fprintf_func(stream, "0x%04x;     ASR64 ACC:P,T", insn);
                                    break;
                                }
                                case 0b1101: //0101 0110 0010 1101 0000 0SHF LLLL LLLL MOV loc16,ACC<<2...8
                                {
                                    if (((insn32 & 0xffff) >> 11) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        uint32_t shift = (insn32 >> 8) & 0b111;
                                        shift += 1;
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; MOV %s,ACC<<%d", insn32, str, shift);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b1111://0101 0110 0010 1111 0000 0SHF LLLL LLLL MOVH loc16,ACC<<2...8
                                {
                                    if (((insn32 & 0xffff) >> 11) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        uint32_t shift = (insn32 >> 8) & 0b111;
                                        shift += 1;
                                        get_loc_string(str,mode,LOC16);
                                        fprintf_func(stream, "0x%08x; MOVH %s,ACC<<%d", insn32, str, shift);
                                        length = 4;
                                    }
                                    break;
                                }
                            }
                            break;
                        case 0b0011: //0101 0110 0011 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //0101 0110 0011 0000 NORM ACC,*0--
                                {
                                    fprintf_func(stream, "0x%04x;     NORM ACC,*0--", insn);
                                    break;
                                }
                                case 0b0010: //0101 0110 0011 0010 NEGTC ACC
                                {
                                    fprintf_func(stream, "0x%04x;     NEGTC ACC", insn);
                                    break;
                                }
                                case 0b0011: //0101 0110 0011 0011 ZAPA
                                {
                                    fprintf_func(stream, "0x%04x;     ZAPA", insn);
                                    break;
                                }
                                case 0b0101: //0101 0110 0011 0101 CSB ACC
                                {
                                    fprintf_func(stream, "0x%04x;     CSB ACC", insn);
                                    break;
                                }
                                case 0b0100: //0101 0110 0011 0100 XCALL *AL
                                {
                                    fprintf_func(stream, "0x%04x;     XCALL *AL", insn);
                                    break;
                                }
                                case 0b0110: //0101 0110 0011 0110 CLRC Objmode
                                {
                                    fprintf_func(stream, "0x%04x;     CLRC Objmode", insn);
                                    break;
                                }   
                                case 0b1000: //0101 0110 0011 100A MOV PM,AX
                                case 0b1001: //0101 0110 0011 100A MOV PM,AX
                                {
                                    uint32_t is_AH = insn & 1;
                                    if (is_AH) {
                                        fprintf_func(stream, "0x%04x;     MOV PM,AH", insn);
                                    }
                                    else {
                                        fprintf_func(stream, "0x%04x;     MOV PM,AL", insn);
                                    }
                                    break;
                                }
                                case 0b1011: //0101 0110 0011 1011 LSLL ACC,T
                                {
                                    fprintf_func(stream, "0x%04x;     LSLL ACC,T", insn);
                                    break;
                                }
                                case 0b1100: //0101 0110 0011 1100 0000 0000 LLLL LLLL XPREAD loc16,*AL
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; XPREAD %s,*AL", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1101: //0101 0110 0011 1101 0000 0000 LLLL LLLL XPWRITE *AL,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; XPWRITE *AL,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1110: //0101 0110 0011 1110 SAT64 ACC:P
                                {
                                    fprintf_func(stream, "0x%04x;     SAT64 ACC:P", insn);
                                    break;
                                }
                                case 0b1111: //0101 0110 0011 1111 CLRC M0M1MAP
                                {
                                    fprintf_func(stream, "0x%04x;     CLRC M0M1MAP", insn);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0100: //0101 0110 0100 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //0101 0110 0100 0000 ADDCL ACC,loc32
                                {
                                    uint32_t mode = insn32 & 0xff;
                                    length = 4;
                                    get_loc_string(str,mode,LOC32);
                                    fprintf_func(stream, "0x%08x; ADDCL ACC,%s", insn32, str);
                                    break;
                                }
                                case 0b0001: //0101 0110 0100 0001 0000 0000 LLLL LLLL SUBL loc32,ACC
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; SUBL %s,ACC", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0010: //0101 0110 0100 0010 0000 0000 LLLL LLLL QMPYXUL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; QMPYXUL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0011: //0101 0110 0100 0011 0000 0000 LLLL LLLL IMPYSL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; IMPYSL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0100: //0101 0110 0100 0100 0000 0000 LLLL LLLL IMPYL ACC,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; IMPYL ACC,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0101: //0101 0110 0100 0101 0000 0000 LLLL LLLL QMPYSL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; QMPYSL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0110: //0101 0110 0100 0110 0000 0000 LLLL LLLL QMPYAL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; QMPYAL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0111: //0101 0110 0100 0111 0000 0000 LLLL LLLL QMPYUL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; QMPYUL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1000: //0101 0110 0100 1000 0000 COND LLLL LLLL MOVL loc32,ACC,COND
                                {
                                    uint32_t mode = insn32 & 0xff;
                                    uint32_t cond = (insn32 >> 8) & 0xf;
                                    if (((insn32 >> 12) & 0xf) == 0) {
                                        get_loc_string(str, mode, LOC32);
                                        get_cond_string(str2, cond);
                                        fprintf_func(stream, "0x%08x; MOVL %s,ACC,%s", insn32, str, str2);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b1001: //0101 0110 0100 1001 0000 0000 LLLL LLLL SUBRL loc32,ACC
                                {
                                    uint32_t mode = insn32 & 0xff;
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; SUBRL %s,ACC", insn32, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b1010: //0101 0110 0100 1010 CCCC CCCC CCCC CCCC OR ACC,#16bit<<#16
                                {
                                    uint32_t imm = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; OR ACC,#%d<<16", insn32, imm);
                                    length = 4;
                                    break;
                                }
                                case 0b1011: //0101 0110 0100 1011 .... DMAC ACC:P,loc32,*XAR7/++
                                {
                                    uint32_t mode1 = insn32 & 0xff;
                                    uint32_t mode2 = (insn32 >> 8) & 0xff;
                                    get_loc_string(str, mode1, LOC32);
                                    if (mode2 == 0b11000111)
                                    {
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; DMAC ACC:P,%s,*XAR7", insn32, str);
                                    }
                                    else if(mode2 == 0b10000111)
                                    {
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; DMAC ACC:P,%s,*XAR7++", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1100: //0101 0110 0100 1100 0000 0000 LLLL LLLL IMPYAL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; IMPYAL P,XT,%s", insn32, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b1101: //0101 0110 0100 1101 .... IMACL P,loc32,*XAR7/++
                                {
                                    uint32_t mode1 = insn32 & 0xff;
                                    uint32_t mode2 = (insn32 >> 8) & 0xff;
                                    get_loc_string(str, mode1, LOC32);
                                    if (mode2 == 0b11000111)
                                    {
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; IMACL P,%s,*XAR7", insn32, str);
                                    }
                                    else if(mode2 == 0b10000111)
                                    {
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; IMACL P,%s,*XAR7++", insn32, str);
                                    }
                                    break;
                                }
                                case 0b1110: //0101 0110 0100 1110 CCCC CCCC CCCC CCCC XOR ACC,#16bit<<#16
                                {
                                    length = 4;
                                    uint32_t imm = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; XOR ACC,#%d<<#16", insn32, imm);
                                    break;
                                }
                                case 0b1111: //0101 0110 0100 1111 .... QMACL P,loc32,*XAR7/++
                                {
                                    uint32_t mode1 = insn32 & 0xff;
                                    uint32_t mode2 = (insn32 >> 8) & 0xff;
                                    get_loc_string(str, mode1, LOC32);
                                    if (mode2 == 0b11000111)
                                    {
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; QMACL P,%s,*XAR7", insn32, str);
                                    }
                                    else if(mode2 == 0b10000111)
                                    {
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; QMACL P,%s,*XAR7++", insn32, str);
                                    }
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0101: //0101 0110 0101 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //0101 0110 0101 0000 0000 0000 LLLL LLLL MINL ACC,loc32
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; MINL ACC,%s", insn, str);
                                    }
                                    break;
                                }
                                case 0b0001: //0101 0110 0101 0001 0000 0000 LLLL LLLL MAXCUL P,loc32
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; MAXCUL P,%s", insn, str);
                                    }
                                    break;
                                }
                                case 0b0010: //0101 0110 0101 0010 LSL64 ACC:P,T
                                {
                                    fprintf_func(stream, "0x%04x;     LSL64 ACC:P,T", insn);
                                    break;
                                }
                                case 0b0011: //0101 0110 0101 0011 xxxx xxxx LLLL LLLL ADDUL ACC,loc32
                                {
                                    uint32_t mode = insn32 & 0xff;
                                    get_loc_string(str, mode, LOC32);
                                    fprintf_func(stream, "0x%08x; ADDUL ACC,%s", insn, str);
                                    length = 4;
                                    break;
                                }
                                case 0b0100: //0101 0110 0101 0100 0000 0000 LLLL LLLL SUBBL ACC,loc32
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; SUBBL ACC,%s", insn, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0101: //0101 0110 0101 0101 0000 0000 LLLL LLLL SUBUL ACC,loc32
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; SUBUL ACC,%s", insn, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b0110: //0101 0110 0101 0110 MOV TL,#0
                                {
                                    fprintf_func(stream, "0x%04x;     MOV TL,#0", insn);
                                    break;
                                }
                                case 0b0111: //0101 0110 0101 0111 0000 0000 LLLL LLLL ADDUL P,loc32
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; ADDUL P,%s", insn, str);
                                    }
                                    length = 4;
                                    break;
                                }
                                case 0b1000: //0101 0110 0101 1000 NEG64 ACC:P
                                {
                                    fprintf_func(stream, "0x%04x;     NEG64 ACC:P", insn);
                                    break;
                                }
                                case 0b1001: //0101 0110 0101 1001 XXXX XXXX LLLL LLLL MINCUL P,loc32
                                {
                                    length = 4;
                                    uint32_t mode = insn32 & 0xff;
                                    get_loc_string(str, mode, LOC32);
                                    fprintf_func(stream, "0x%08x; MINCUL P,%s", insn, str);
                                    break;
                                }
                                case 0b1010: //0101 0110 0101 1010 NORM ACC,*++
                                {
                                    fprintf_func(stream, "0x%04x;     NORM ACC,*++", insn);
                                    break;
                                }
                                case 0b1011: //0101 0110 0101 1011 LSR64 ACC:P,T
                                {
                                    fprintf_func(stream, "0x%04x;     LSR64 ACC:P,T", insn);
                                    break;
                                }
                                case 0b1100: //0101 0110 0101 1100 CLRC OVC
                                {
                                    fprintf_func(stream, "0x%04x;     CLRC OVC", insn);
                                    break;
                                }
                                case 0b1101: //0101 0110 0101 1101 0000 0000 LLLL LLLL SUBUL P,loc32
                                {
                                    if (((insn32 & 0xff00) >> 8) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; SUBUL P,%s", insn, str);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b1110: //0101 0110 0101 1110 CMP64 ACC:P
                                {
                                    fprintf_func(stream, "0x%04x;     CMP64 ACC:P", insn);
                                    break;
                                }
                                case 0b1111: //0101 0110 0101 1111 ABSTC  ACC
                                {
                                    fprintf_func(stream, "0x%04x;     ABCTC ACC", insn);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0110: //0101 0110 0110 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0001: //0101 0110 0110 0001 0000 0000 LLLL LLLL MAXL ACC,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        length = 4;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; MAXL ACC,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0010: //0101 0110 0110 0010 0000 0000 LLLL LLLL MOVU OVC,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        length = 4;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MOVU OVC,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0011: //0101 0110 0110 0011 0000 0000 LLLL LLLL QMPYL ACC,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; QMPYL ACC,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0101: //0101 0110 0110 0101 0000 0000 LLLL LLLL IMPYXUL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0) {
                                        uint32_t mode = insn32 & 0xff;
                                        length = 4;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; IMPYXUL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                                case 0b0111: //0101 0110 0110 0111 0000 0000 LLLL LLLL QMPYL P,XT,loc32
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC32);
                                        fprintf_func(stream, "0x%08x; QMPYL P,XT,%s", insn32, str);
                                    }
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0111: //0101 0110 0111 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //0101 0110 0111 0000 FLIP AL
                                {
                                    fprintf_func(stream, "0x%04x;     FLIP AL", insn);
                                    break;
                                }
                                case 0b0001: //0101 0110 0111 0001 FLIP AH
                                {
                                    fprintf_func(stream, "0x%04x;     FLIP AH", insn);
                                    break;
                                }
                                case 0b0010://0101 0110 0111 0010 0000 0000 LLLL LLLL MAX AL,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MAX AL,%s", insn, str);
                                    }
                                    break;
                                }
                                case 0b0011://0101 0110 0111 0011 0000 0000 LLLL LLLL MAX AH,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MAX AH,%s", insn, str);
                                    }
                                    break;
                                }
                                case 0b0100://0101 0110 0111 0100 0000 0000 LLLL LLLL MIN AL,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MIN AL,%s", insn, str);
                                    }
                                    break;
                                }
                                case 0b0101://0101 0110 0111 0101 0000 0000 LLLL LLLL MIN AH,loc16
                                {
                                    if (((insn32 >> 8) & 0xff) == 0)
                                    {
                                        length = 4;
                                        uint32_t mode = insn32 & 0xff;
                                        get_loc_string(str, mode, LOC16);
                                        fprintf_func(stream, "0x%08x; MIN AH,%s", insn, str);
                                    }
                                    break;
                                }
                                case 0b0111://0101 0110 0111 0111 NORM ACC,*0++
                                {
                                    fprintf_func(stream, "0x%04x;     NORM ACC,*0++", insn);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b1000: //0101 0110 1000 SHFT ASR64 ACC:P,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%04x;     ASR64 ACC:P,#%d", insn, shift);
                            break;
                        }
                        case 0b1001: //0101 0110 1001 SHFT
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%04x;     LSR64 ACC:P,#%d", insn, shift);
                            break;
                        }
                        case 0b1010: //0101 0110 1010 SHFT LSL64 ACC:P,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%04x;     LSL64 ACC:P,#%d", insn, shift);
                            break;
                        }
                        case 0b1011: //0101 0110 1011 COND CCCC CCCC LLLL LLLL MOVB loc16,#8bit,COND
                        {
                            uint32_t cond = insn & 0xf;
                            uint32_t imm = (insn32 >> 8) & 0xff;
                            uint32_t mode = insn32 & 0xff;
                            get_loc_string(str, mode, LOC16);
                            get_cond_string(str2, cond);
                            fprintf_func(stream, "0x%08x; MOVB %s,#%d,%s", insn32, str, imm, str2);
                            length = 4;
                            break;
                        }
                        case 0b1100: //0101 0110 1100 COND  BF 16bitOffset,COND
                        {
                            int16_t imm = (uint16_t)(insn32&0xffff);
                            int16_t cond = (insn&0xf);
                            length = 4;
                            get_cond_string(str, cond);
                            fprintf_func(stream, "0x%08x; BF #%d,%s", insn32, imm, str);
                            break;
                        }
                        case 0b1101: //0101 0110 1101 COND CCCC CCCC CCCC CCCC XB pma,COND
                        {
                            uint32_t addr = insn32 & 0xffff;
                            uint32_t cond = insn & 0xf;
                            get_cond_string(str, cond);
                            fprintf_func(stream, "0x%08x; XB 0x%x,%s", insn32, addr, str);
                            length = 4;
                            break;
                        }
                        case 0b1110: //0101 0110 1110 COND CCCC CCCC CCCC CCCC XCALL pma,COND
                        {
                            uint32_t addr = insn32 & 0xffff;
                            uint32_t cond = insn & 0xf;
                            get_cond_string(str, cond);
                            fprintf_func(stream, "0x%08x; XCALL 0x%x,%s", insn32, addr, str);
                            length = 4;
                            break;
                        }
                        case 0b1111: //0101 0110 1111 COND XRETC COND
                        {
                            uint32_t cond = insn & 0xf;
                            if (cond == 0xff)
                            {
                                fprintf_func(stream, "0x%04x;     XRET", insn);
                            }
                            else
                            {
                                get_cond_string(str, cond);
                                fprintf_func(stream, "0x%04x;     XRETC %s", insn, str);
                            }
                            break;
                        }
                    }
                    break;
                case 0b0111://0101 0111 LLLL LLLL MOVH loc16,P
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVH %s,P", insn, str);
                    break;
                }
                case 0b1000:
                case 0b1001:
                case 0b1010:
                case 0b1011:
                case 0b1100:
                case 0b1101://0101 1nnn LLLL LLLL MOVZ ARn,loc16
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t n = (insn >> 8) & 0b111;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVZ AR%d,%s", insn, n, str);
                    break;
                }
                case 0b1110://0101 1110 LLLL LLLL MOV AR6,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOV AR6,%s", insn, str);
                    break;
                }
                case 0b1111://0101 1111 LLLL LLLL MOV AR7,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOV AR7,%s", insn, str);
                    break;
                }
            }
            break;
        case 0b0110: //0110 COND CCCC CCCC SB 8bitOffset,COND
        {
            int16_t offset = insn & 0xff;
            if ((offset >> 7) == 1) {
                offset = offset | 0xff00;
            }
            int16_t cond = (insn>>8)&0xf;
            get_cond_string(str, cond);
            fprintf_func(stream, "0x%04x;     SB #%d,%s", insn, offset, str);
            break;
        }
        case 0b0111:
            switch((insn & 0x0f00) >> 8) {
                case 0b0000: //0111 0000 LLLL LLLL XOR AL,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     XOR AL,%s", insn, str);
                    break;
                }
                case 0b0001: //0111 0001 LLLL LLLL XOR AH,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     XOR AH,%s", insn, str);
                    break;
                }
                case 0b0010: //0111 0010 LLLL LLLL ADD loc16,AL
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     ADD %s,AL", insn, str);
                    break;
                }
                case 0b0011: //0111 0011 LLLL LLLL ADD loc16,AH
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     ADD %s,AH", insn, str);
                    break;
                }
                case 0b0100://0111 0100 LLLL LLLL SUB loc16,AL
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     SUB %s,AL", insn, str);
                    break;
                }
                case 0b0101://0111 0101 LLLL LLLL SUB loc16,AH
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     SUB %s,AH", insn, str);
                    break;
                }
                case 0b0110: //0111 0110 .... ....
                    if (((insn >> 7) & 0x1) == 1) { //0111 0110 1... .... MOVL XARn,#22bit 
                        uint32_t n = ((insn >> 6) & 0b1) + 6;//XAR6,XAR7
                        uint32_t imm = insn32 & 0x3fffff;
                        fprintf_func(stream, "0x%08x; MOVL XAR%d, #0x%x", insn32, n,imm);
                        length = 4;
                    }
                    else { // 0111 0110 0... ....
                        if (((insn >> 6) & 0x1) == 1) { //0111 0110 01.. .... LCR #22bit
                            length = 4;
                            uint32_t imm = insn32 & 0x3fffff;
                            fprintf_func(stream, "0x%08x; LCR @0x%x", insn32, imm);
                        }
                        else { //0111 0110 00.. ....
                            switch(insn & 0x3f) {
                                case 0b000000: //0111 0110 0000 0000 POP ST1
                                {
                                    fprintf_func(stream, "0x%04x;     POP ST1", insn);
                                    break;
                                }
                                case 0b000001: //0111 0110 0000 0001 POP DP:ST1
                                {
                                    fprintf_func(stream, "0x%04x;     POP DP:ST1", insn);
                                    break;
                                }
                                case 0b000010: //0111 0110 0000 0010 IRET
                                {
                                    fprintf_func(stream, "0x%04x;     IRET", insn);
                                    break;
                                }
                                case 0b000011: //0111 0110 0000 0011 POP DP
                                {
                                    fprintf_func(stream, "0x%04x;     POP DP", insn);
                                    break;
                                }
                                case 0b000100: //0111 0110 0000 0100 LC *XAR7
                                {
                                    fprintf_func(stream, "0x%04x;     LC *XAR7", insn);
                                    break;
                                }
                                case 0b000101: //0111 0110 0000 0111 POP AR3:AR2
                                {
                                    fprintf_func(stream, "0x%04x;     POP AR3:AR2", insn);
                                    break;
                                }
                                case 0b000110: //0111 0110 0000 0111 POP AR5:AR4
                                {
                                    fprintf_func(stream, "0x%04x;     POP AR5:AR4", insn);
                                    break;
                                }
                                case 0b000111: //0111 0110 0000 0111 POP AR1:AR0
                                {
                                    fprintf_func(stream, "0x%04x;     POP AR1:AR0", insn);
                                    break;
                                }
                                case 0b001000: //0111 0110 0000 1000 PUSH ST1
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH ST1", insn);
                                    break;
                                }
                                case 0b001001: //0111 0110 0000 1001 PUSH DP:ST1
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH DP:ST1", insn);
                                    break;
                                }
                                case 0b001010: //0111 0110 0000 1010 PUSH IFR
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH IFR", insn);
                                    break;
                                }
                                case 0b001011: //0111 0110 0000 1011 PUSH DP
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH DP", insn);
                                    break;
                                }
                                case 0b001100: //0111 0110 0000 1100 PUSH AR5:AR4
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH AR5:AR4", insn);
                                    break;
                                }
                                case 0b001101: //0111 0110 0000 1101 PUSH AR1:AR0
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH AR1:AR0", insn);
                                    break;
                                }
                                case 0b001110: //0111 0110 0000 1110 PUSH DBGIER
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH DBGIER", insn);
                                    break;
                                }
                                case 0b001111: //0111 0110 0000 1111 PUSH AR3:AR2
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH AR3:AR2", insn);
                                    break;
                                }
                                case 0b010000: //0111 0110 0001 0000 LRETE
                                {
                                    fprintf_func(stream, "0x%04x;     LRETE", insn);
                                    break;
                                }
                                case 0b010001: //0111 0110 0001 0001 POP P
                                {
                                    fprintf_func(stream, "0x%04x;     POP P", insn);
                                    break;
                                }
                                case 0b010010: //0111 0110 0001 0010 POP DBGIER
                                {
                                    fprintf_func(stream, "0x%04x;     POP DBGIER", insn);
                                    break;
                                }
                                case 0b010011: //0111 0110 0001 0011 POP ST0
                                {
                                    fprintf_func(stream, "0x%04x;     POP ST0", insn);
                                    break;
                                }
                                case 0b010100: //0111 0110 0001 0100 LRET
                                {
                                    fprintf_func(stream, "0x%04x;     LRET", insn);
                                    break;
                                }
                                case 0b010101: //0111 0110 0001 0101 POP T:ST0
                                {
                                    fprintf_func(stream, "0x%04x;     POP T:ST0", insn);
                                    break;
                                }
                                case 0b010111: //0111 0110 0001 0111 NASP
                                {
                                    fprintf_func(stream, "0x%04x;     NASP", insn);
                                    break;
                                }
                                case 0b011000: //0111 0110 0001 1000 PUSH ST0
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH ST0", insn);
                                    break;
                                }
                                case 0b011001: //0111 0110 0001 1001 PUSH T:ST0
                                {
                                    fprintf_func(stream, "0x%04x;     PUSH T:ST0", insn);
                                    break;
                                }
                                case 0b011010: //0111 0110 0001 1010 EDIS
                                {
                                    fprintf_func(stream, "0x%04x;     EDIS", insn);
                                    break;
                                }
                                case 0b011011: //0111 0110 0001 1011 ASP
                                {
                                    fprintf_func(stream, "0x%04x;     ASP", insn);
                                    break;
                                }
                                case 0b011111: //0111 0110 0001 1111 MOVW DP,#16bit
                                {
                                    uint32_t imm = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; MOVW DP, #0x%04x", insn32, imm);
                                    length = 4;
                                    break;
                                }
                                case 0b100000: //0111 0110 0010 0000 LB *XAR7
                                {
                                    fprintf_func(stream, "0x%04x;     LB *XAR7", insn);
                                    break;
                                }
                                case 0b100001: //0111 0110 0010 0001 IDLE
                                {
                                    fprintf_func(stream, "0x%04x;     IDLE", insn);
                                    break;
                                }
                                case 0b100010: //0111 0110 0010 0010 EALLOW
                                {
                                    fprintf_func(stream, "0x%04x;     EALLOW", insn);
                                    break;
                                }
                                case 0b100011: //0111 0110 0010 0011 CCCC CCCC CCCC CCCC OR IER,#16bit
                                {
                                    uint32_t imm = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; OR IER,#0x%04x", insn32, imm);
                                    length = 4;
                                    break;
                                }
                                case 0b100100: //0111 0110 0010 0100 ESTOP1
                                {
                                    fprintf_func(stream, "0x%04x;     ESTOP1", insn);
                                    break;
                                }
                                case 0b100101: //0111 0110 0010 0101 ESTOP0
                                {
                                    fprintf_func(stream, "0x%04x;     ESTOP0", insn);
                                    break;
                                }
                                case 0b100110: //0111 0110 0010 0110 CCCC CCCC CCCC CCCC AND IER,#16bit
                                {
                                    uint32_t imm = insn32 & 0xffff;
                                    length = 4;
                                    fprintf_func(stream, "0x%08x; AND IER,#0x%04x", insn32, imm);
                                    break;
                                }
                                case 0b100111: //0111 0110 0010 0111 CCCC CCCC CCCC CCCC OR IFR,#16bit
                                {
                                    uint32_t imm = insn32 & 0xffff;
                                    length = 4;
                                    fprintf_func(stream, "0x%08x; OR IFR,#0x%04x", insn32, imm);
                                    break;
                                }
                                case 0b101111: //0111 0110 0010 1111 CCCC CCCC CCCC CCCC AND IFR,#16bit
                                {
                                    uint32_t imm = insn32 & 0xffff;
                                    length = 4;
                                    fprintf_func(stream, "0x%08x; AND IFR,#0x%04x", insn32, imm);
                                    break;
                                }
                                case 0b111111://0111 0110 0011 1111 CCCC CCCC CCCC CCCC IACK #16bit
                                {
                                    length = 4;
                                    uint32_t imm = insn32 & 0xffff;
                                    fprintf_func(stream, "0x%08x; IACK #0x%04x", insn32, imm);
                                    break;
                                }
                            }
                        }
                    }
                    break;
                case 0b0111: //0111 0111 LLLL LLLL NOP {*ind}{ARPn}
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0)
                    {
                        fprintf_func(stream, "0x%04x;     NOP", insn);
                    }
                    else
                    {
                        get_loc_string(str,mode,LOC16);
                        fprintf_func(stream, "0x%04x;     NOP %s", insn, str);
                    }
                    break;
                }
                case 0b1000: 
                case 0b1001: 
                case 0b1010: 
                case 0b1011: 
                case 0b1100: 
                case 0b1101: 
                case 0b1110: 
                case 0b1111: //0111 1nnn LLLL LLLL MOV loc16,ARn
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t n = (insn >> 8) & 0b111;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV %s,AR%d", insn, str, n);
                }
            }
            break;
        case 0b1000:
            switch ((insn & 0x0f00) >> 8) {
                case 0b0000: //1000 0000 LLLL LLLL MOVZ AR7,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVZ AR7,%s", insn, str);
                    break;
                }
                case 0b0001: //1000 0001 LLLL LLLL ADD ACC,loc16<<#0
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     ADD ACC,%s<<#0", insn, str);
                    break;
                }
                case 0b0010: //1000 0010 LLLL LLLL MOVL XAR3,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110)
                    {
                        fprintf_func(stream, "0x%04x;     POP XAR3", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR3,%s", insn, str);
                    }
                    break;
                }
                case 0b0011: //1000 0011 LLLL LLLL MOVL XAR5,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110)
                    {
                        fprintf_func(stream, "0x%04x;     POP XAR5", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR5,%s", insn, str);
                    }
                    break;
                }
                case 0b0101: /*1000 0101 .... .... MOV ACC, loc16<<#0 */
                {                  
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV ACC,%s<<#0", insn, str);
                    break;
                }
                case 0b0100: //1000 0100 LLLL LLLL CCCC CCCC CCCC CCCC XMAC P,loc16,*(pma)
                {
                    length = 4;
                    uint32_t mode = insn & 0xff;
                    uint32_t addr = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x;     XMAC P,%s,*(0x%x)", insn32, str, addr);
                    break;
                }
                case 0b0110: //1000 0110 LLLL LLLL MOVL XAR2,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110)
                    {
                        fprintf_func(stream, "0x%04x;     POP XAR2", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR2,%s", insn, str);
                    }
                    break;
                }
                case 0b0111: //1000 0111 LLLL LLLL MOVL XT,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110) {
                        fprintf_func(stream, "0x%04x;     POP XT", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XT,%s", insn, str);
                    }
                    break;
                }
                case 0b1000: //1000 1000 LLLL LLLL MOVZ AR6,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVZ AR6,%s", insn, str);
                    break;
                }
                case 0b1001: //1000 1001 LLLL LLLL AND ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     AND ACC,%s", insn, str);
                    break;
                }
                case 0b1010: //1000 1010 LLLL LLLL MOVL XAR4,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110)
                    {
                        fprintf_func(stream, "0x%04x;     POP XAR4", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR4,%s", insn, str);
                    }
                    break;
                }
                case 0b1011: //1000 1011 LLLL LLLL MOVL XAR1,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110)
                    {
                        fprintf_func(stream, "0x%04x;     POP XAR1", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR1,%s", insn, str);
                    }
                    break;
                }
                case 0b1100: //1000 1100 LLLL LLLL CCCC CCCC CCCC CCCC MPY P,loc16,#16bit
                {
                    uint32_t mode = insn & 0xff;
                    int32_t imm = (int16_t)(insn32 & 0xffff);
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; MPY P,%s,#%d", insn, str, imm);
                    length = 4;
                    break;
                }
                case 0b1101: /*1000 1101 .... .... MOVL XARn,#22bit */
                {
                    uint32_t n = (insn >> 6) & 0b11;
                    uint32_t imm = insn32 & 0x3fffff;
                    fprintf_func(stream, "0x%08x; MOVL XAR%d, #0x%x", insn32, n, imm);
                    length = 4;
                    break;
                }
                case 0b1110: //1000 1110 LLLL LLLL MOVL XAR0,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110)
                    {
                        fprintf_func(stream, "0x%04x;     POP XAR0", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR0,%s", insn, str);
                    }
                    break;
                }
                case 0b1111:
                    switch((insn >> 6) & 0b11) {
                        case 0b00: //1000 1111 00cc cccc cccc cccc cccc cccc MOVL XAR4,#22bit
                        {
                            uint32_t imm = insn32 & 0x3fffff;
                            fprintf_func(stream, "0x%08x; MOVL XAR4, #0x%x", insn32, imm);
                            length = 4;
                            break;
                        }
                        case 0b01: //1000 1111 01cc cccc cccc cccc cccc cccc MOVL XAR5,#22bit
                        {
                            uint32_t imm = insn32 & 0x3fffff;
                            fprintf_func(stream, "0x%08x; MOVL XAR5, #0x%x", insn32, imm);
                            length = 4;
                            break;
                        }
                        case 0b10: //1000 1111 10nn nmmm CCCC CCCC CCCC CCCC BAR 16bitOffset,ARn,ARm,EQ
                        {
                            uint32_t n = (insn >> 3) & 0b111;
                            uint32_t m = insn & 0b111;
                            int16_t imm = (uint16_t)(insn32&0xffff);
                            fprintf_func(stream, "0x%08x; BAR #%d,AR%d,AR%d,EQ", insn32, imm, n, m);
                            length = 4;
                            break;
                        }
                        case 0b11: //1000 1111 11nn nmmm CCCC CCCC CCCC CCCC BAR 16bitOffset,ARn,ARm,NEQ
                        {
                            uint32_t n = (insn >> 3) & 0b111;
                            uint32_t m = insn & 0b111;
                            int16_t imm = (uint16_t)(insn32&0xffff);
                            fprintf_func(stream, "0x%08x; BAR #%d,AR%d,AR%d,NEQ", insn32, imm, n, m);
                            length = 4;
                            break;
                        }
                    }
                    break;
            }
            break;
        case 0b1001:
            switch ((insn >> 9) & 0b111) {
                case 0b000: //1001 000a CCCC CCCC ANDB AX,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    uint32_t is_AH = ((insn >> 8) & 1);
                    if (is_AH)
                        fprintf_func(stream, "0x%04x;     ANDB AH,#%d", insn, imm);
                    else 
                        fprintf_func(stream, "0x%04x;     ANDB AL,#%d", insn, imm);
                    break;
                }
                case 0b001: //1001 001. .... ....  MOV AX,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    if (((insn >> 8) & 1) == 1){ //1001 0011 .... .... ah
                        fprintf_func(stream, "0x%04x;     MOV AH,%s", insn, str);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     MOV AL,%s", insn, str);
                    }
                    break;
                }
                case 0b010: //1001 010. .... .... ADD AX,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    if (((insn >> 8) & 1) == 1){ //1001 0011 .... .... ah
                        fprintf_func(stream, "0x%04x;     ADD AH,%s", insn, str);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     ADD AL,%s", insn, str);
                    }
                    break;
                }
                case 0b011: //1001 011. .... .... MOV loc16,AX
                {                    
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    if (((insn >> 8) & 1) == 1){ //1001 0011 .... .... ah
                        fprintf_func(stream, "0x%04x;     MOV %s,AH", insn, str);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     MOV %s,AL", insn, str);
                    }
                    break;
                }
                case 0b100: //1001 100A LLLL LLLL OR loc16,AX
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    if (((insn >> 8) & 1) == 1){ //1001 0011 .... .... ah
                        fprintf_func(stream, "0x%04x;     OR %s,AH", insn, str);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     OR %s,AL", insn, str);
                    }
                    break;
                }
                case 0b101: //1001 101A CCCC CCCC MOVB AX,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    if (((insn >> 8) & 1) == 1){ //1001 0011 .... .... ah
                        fprintf_func(stream, "0x%04x;     MOVB AH,#%d", insn, imm);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     MOVB AL,#%d", insn, imm);
                    }
                    break;
                }
                case 0b110: //1001 110A CCCC CCCC ADDB AX,#8bitSigned
                {
                    int32_t imm = insn & 0xff;
                    if ((imm >> 7) == 1) {
                        imm = imm | 0xffffffffff00;
                    }
                    if (((insn >> 8) & 1) == 1){ //1001 0011 .... .... ah
                        fprintf_func(stream, "0x%04x;     ADDB AH,#%d", insn, imm);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     ADDB AL,#%d", insn, imm);
                    }
                    break;
                }
                case 0b111: //1001 111A LLLL LLLL SUB AX,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    if (((insn >> 8) & 1) == 1){ //ah
                        fprintf_func(stream, "0x%04x;     SUB AH,%s", insn, str);
                    }
                    else { // al
                        fprintf_func(stream, "0x%04x;     SUB AL,%s", insn, str);
                    }
                    break;
                }
            }
            break;
        case 0b1010:
            switch ((insn & 0xf00) >> 8) {
                case 0b0000: //1010 0000 LLLL LLLL MOVL loc32, XAR5
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR5", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR5", insn, str);
                    }
                    break;
                }
                case 0b0010: //1010 0010 LLLL LLLL MOVL loc32, XAR3
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR3", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR3", insn, str);
                    }
                    break;
                }
                case 0b0011: //1010 0011 LLLL LLLL MOVL P,loc32
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC32);
                    fprintf_func(stream, "0x%04x;     MOVL P,%s", insn, str);
                    break;
                }
                case 0b0101: //1010 0101 LLLL LLLL DMOV loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     DMOV %s", insn, str);
                    break;
                }
                case 0b0100: //1010 0100 LLLL LLLL CCCC CCCC CCCC CCCC XMACD P,loc16,*(pma)
                {
                    length = 4;
                    uint32_t mode = insn & 0xff;
                    uint32_t addr = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; XMACD P,%s,*(0x%x)", insn32, str, addr);
                    break;
                }
                case 0b0110: //1010 0110 LLLL LLLL MOVDL XT,loc32
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC32);
                    fprintf_func(stream, "0x%04x;     MOVDL XT,%s", insn, str);
                    break;
                }
                case 0b0111: //1010 0111 LLLL LLLL MOVAD T,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVAD T,%s", insn, str);
                    break;
                }
                case 0b1000: //1010 1000 LLLL LLLL MOVL loc32, XAR4
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR4", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR4", insn, str);
                    }
                    break;
                }
                case 0b1001: //1010 1001 LLLL LLLL MOVL loc32,P
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC32);
                    fprintf_func(stream, "0x%04x;     MOVL %s,P", insn, str);
                    break;
                }
                case 0b1010: //1010 1010 LLLL LLLL MOVL loc32, XAR2
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR2", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR2", insn, str);
                    }
                    break;
                }
                case 0b1011: //1010 1011 LLLL LLLL MOVL loc32,XT
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XT", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XT", insn, str);
                    }
                    break;
                }
                case 0b1100: //1010 1100 MMMM MMMM LLLL LLLL LLLL LLLL XPREAD loc16,*(pma)
                {
                    length = 4;
                    uint32_t mode = insn & 0xff;
                    uint32_t addr = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; XPREAD %s,*(#0x%x)", insn32, str, addr);
                    break;
                }
                case 0b1101: //1010 1101 FFFF FFFF MOVST0 FLAG
                {
                    uint32_t flag = insn & 0xff;
                    get_stf_bit_string(str, flag);
                    fprintf_func(stream, "0x%04x;     MOVST0 %s", insn, str);
                    break;
                }
                case 0b1110: //1010 1110 LLLL LLLL SUB ACC,loc16<<#0
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     SUB ACC,%s<<0", insn, str);
                    break;
                }
                case 0b1111: //1010 1111 LLLL LLLL OR ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     OR ACC,%s", insn, str);
                    break;
                }
            }
            break;
        case 0b1011:
            switch ((insn & 0xf00) >> 8) {
                case 0b0000: //1011 0000 LLLL LLLL CCCC CCCC CCCC CCCC UOUT *(PA),loc16
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t addr  = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; UOUT *(0x%x),%s", insn32, addr, str);
                    length = 4;
                    break;
                }
                case 0b0001: //1011 0001 LLLL LLLL MOV loc16,ACC<<1
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOV %s,ACC<<1", insn, str);
                    break;
                }
                case 0b0010: //1011 0010 LLLL LLLL MOVL loc32, XAR1
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR1", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR1", insn, str);
                    }
                    break;
                }
                case 0b0011: //1011 0011 LLLL LLLL MOVH loc16,ACC<<1
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVH %s,ACC<<1", insn, str);
                    break;
                }
                case 0b0100: //1011 0100 LLLL LLLL CCCC CCCC CCCC CCCC IN loc16,*(PA)
                {
                    length = 4;
                    uint32_t mode = insn & 0xff;
                    uint32_t pa = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; IN %s,*(0x%x)", insn32, str, pa);
                    break;
                }
                case 0b0101: //1011 0101 .... ....
                {
                    if (((insn >>7) & 1) == 0)
                    {//1011 0101 0bbb bbbb 0000 0000 loc16 RPTB label,loc16
                        length = 4;
                        uint32_t size = insn & 0x7f;
                        uint32_t mode = insn32 & 0xff;
                        get_loc_string(str,mode,LOC16);
                        fprintf_func(stream, "0x%08x; RPTB #%d,@%s", insn32, size, str);
                    }
                    else
                    {//1011 0101 1bbb bbbb cccc cccc cccc cccc RPTB label,#RC
                        length = 4;
                        uint32_t size = insn & 0x7f;
                        uint32_t repeat_count = insn32 & 0xffff;
                        fprintf_func(stream, "0x%08x; RPTB #%d,#0x%x", insn32, size, repeat_count);
                    }
                    break;
                }
                case 0b0110: //1011 0110 CCCC CCCC MOVB XAR7,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     MOVB XAR7,#%d", insn, imm);
                    break;
                }
                case 0b0111: //1011 0111 LLLL LLLL XOR ACC,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     XOR ACC,%s", insn, str);
                    break;
                }
                case 0b1000:
                case 0b1001:
                case 0b1010:
                case 0b1011://1011 10CC CCCC CCCC MOVZ DP,#10bit
                {
                    uint32_t imm = insn & 0x3ff;
                    fprintf_func(stream, "0x%04x;     MOVZ DP,#%d", insn, imm);
                    break;
                }
                case 0b1100://1011 1100 LLLL LLLL CCCC CCCC CCCC CCCC OUT *(PA),loc16
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t pa = insn32 & 0xffff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%08x; OUT *(0x%x),%s", insn32, pa, str);
                    length = 4;
                    break;
                }
                case 0b1101://1011 1101 loc32 iiii iiii iiii iiii MOV32 RaH, ACC
                {
                    length = 4;
                    uint32_t loc32 = insn & 0xff;
                    uint32_t addr = insn32 & 0xffff;
                    get_loc_string(str,loc32,LOC32);
                    if (get_fpu_reg_name(str2, addr))
                    {
                        fprintf_func(stream, "0x%08x; MOV32 %s,%s", insn32, str2, str);
                    }
                    else
                    {
                        fprintf_func(stream, "0x%08x; MOV32 *(0x%x),%s", insn32, addr, str);
                    }
                    break;
                }
                case 0b1110: //1011 1110 CCCC CCCC MOVB XAR6,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     MOVB XAR6,#%d", insn, imm);
                    break;
                }
                case 0b1111: //1011 1111 loc32 iiii iiii iiii iiii MOV32 loc32,*(0:16bitAddr)/ MOV32 ACC,RaH
                {
                    length = 4;
                    uint32_t loc32 = insn & 0xff;
                    uint32_t addr = insn32 & 0xffff;
                    get_loc_string(str,loc32,LOC32);
                    if (get_fpu_reg_name(str2, addr))
                    {
                        fprintf_func(stream, "0x%08x; MOV32 %s,%s", insn32, str, str2);
                    }
                    else
                    {
                        fprintf_func(stream, "0x%08x; MOV32 %s,*(0x%x)", insn32, str, addr);
                    }
                    break;
                }
            }
            break;
        case 0b1100:
            switch ((insn & 0xf00) >> 8) {
                case 0b0000:
                case 0b0001: //1100 000A LLLL LLLL AND loc16,AX
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    uint32_t is_AH = (insn >> 8) & 1;
                    if (is_AH)
                        fprintf_func(stream, "0x%04x;     AND %s,AH", insn, str);
                    else 
                        fprintf_func(stream, "0x%04x;     AND %s,AL", insn, str);
                    break;
                }
                case 0b0010: //1100 0010 LLLL LLLL MOVL loc32, XAR6
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR6", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR6", insn, str);
                    }
                    break;
                }
                case 0b0011: //1100 0011 LLLL LLLL MOVL loc32, XAR7
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111101) {
                        fprintf_func(stream, "0x%04x;     PUSH XAR7", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL %s,XAR7", insn, str);
                    }
                    break;
                }
                case 0b0100: //1100 0100 LLLL LLLL MOVL XAR6,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110) {
                        fprintf_func(stream, "0x%04x;     POP XAR6", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR6,%s", insn, str);
                    }
                    break;
                }
                case 0b0101: //1100 0101 LLLL LLLL MOVL XAR7,loc32
                {
                    uint32_t mode = insn & 0xff;
                    if (mode == 0b10111110) {
                        fprintf_func(stream, "0x%04x;     POP XAR7", insn);
                    }
                    else {
                        get_loc_string(str,mode,LOC32);
                        fprintf_func(stream, "0x%04x;     MOVL XAR7,%s", insn, str);
                    }
                    break;
                }
                case 0b0110: //1100 0110 LLLL LLLL MOVB AL.LSB,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB AL.LSB,%s", insn, str);
                    break;
                }
                case 0b0111: //1100 0111 LLLL LLLL MOVB AH.LSB,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB AH.LSB,%s", insn, str);
                    break;
                }
                case 0b1000://1100 100a LLLL LLLL MOVB loc16,AL.MSB
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB %s,AL.MSB", insn, str);
                    break;
                }
                case 0b1001://1100 100a LLLL LLLL MOVB loc16,AH.MSB
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     MOVB %s,AH.MSB", insn, str);
                    break;
                }
                case 0b1010: //1100 1010 LLLL LLLL OR AL,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     OR AL,%s", insn, str);
                    break;
                }
                case 0b1011: //1100 1011 LLLL LLLL OR AH,loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     OR AH,%s", insn, str);
                    break;
                }
                case 0b1100:
                case 0b1101: //1100 110a LLLL LLLL CCCC CCCC CCCC CCCC AND AX,loc16,#16bit
                {
                    length = 4;
                    uint32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    uint32_t is_AH = (insn >> 8) & 1;
                    get_loc_string(str, mode, LOC16);
                    if (is_AH)
                        fprintf_func(stream, "0x%08x; AND AH,%s,#0x%04x", insn32, str, imm);
                    else
                        fprintf_func(stream, "0x%08x; AND AL,%s,#0x%04x", insn32, str, imm);
                    break;
                }
                case 0b1110:
                case 0b1111: //1100 111a LLLL LLLL AND AX,loc16
                {
                    uint32_t mode = insn & 0xff;
                    uint32_t is_AH = (insn >> 8) & 1;
                    get_loc_string(str, mode, LOC16);
                    if (is_AH)
                        fprintf_func(stream, "0x%04x;     AND AH,%s", insn, str);
                    else 
                        fprintf_func(stream, "0x%04x;     AND AL,%s", insn, str);
                    break;
                }

            }
            break;
        case 0b1101:
            if (((insn >> 11) & 1) == 0) {//1101 0... .... ....
                uint32_t imm = insn & 0xff;
                uint32_t n = (insn >> 8) & 0b111;
                if (n < 6) {//1101 0nnn CCCC CCCC MOVB XAR1...5 ,#8bit
                    fprintf_func(stream, "0x%04x;     MOVB XAR[%d],#%d", insn, n, imm);
                }
                else {//1101 0110/1 CCCC CCCC MOVB AR6/7, #8bit
                    fprintf_func(stream, "0x%04x;     MOVB AR[%d],#%d", insn, n, imm);
                }
            }
            else {//1101 1... .... ....
                uint32_t imm = insn & 0x7f;
                uint32_t n = (insn >> 8) & 0b111;
                if (((insn & 0xff) >> 7) == 0) { //ADDB XARn,#7bit
                    fprintf_func(stream, "0x%04x;     ADDB XAR[%d],#%d", insn, n, imm);
                }
                else { //SUBB XARn,#7bit
                    fprintf_func(stream, "0x%04x;     SUBB XAR[%d],#%d", insn, n, imm);
                }
            }
            break;
        case 0b1110:
        {
            switch((insn & 0x0f00) >> 8) {
                case 0b0000: //1110 0000 .... ....
                {
                    switch((insn & 0x00f0) >> 4) {
                        case 0b0001: //1110 0000 0001 fffe eedd daaa mem32 ADDF32 RdH, ReH, RfH || MOV32 mem32, RaH
                        {
                            length = 4;
                            uint32_t mem32 = insn32 & 0xff;
                            uint32_t a = (insn32 >> 8) & 0b111;
                            uint32_t d = (insn32 >> 11) & 0b111;
                            uint32_t e = (insn32 >> 14) & 0b111;
                            uint32_t f = (insn >> 1) & 0b111;
                            get_loc_string(str, mem32, LOC32);
                            fprintf_func(stream, "0x%08x; ADDF32 R%dH,R%dH,R%dH ||MOV32 %s,R%dH", insn32, d, e, f, str, a);
                            break;
                        }
                        case 0b0010: //1110 0000 0010 fffe eedd daaa mem32 SUBF32 RdH,ReH,RfH || MOV32 mem32,RaH
                        {
                            length = 4;
                            uint32_t mem32 = insn32 & 0xff;
                            uint32_t a = (insn32 >> 8) & 0b111;
                            uint32_t d = (insn32 >> 11) & 0b111;
                            uint32_t e = (insn32 >> 14) & 0b111;
                            uint32_t f = (insn >> 1) & 0b111;
                            get_loc_string(str, mem32, LOC32);
                            fprintf_func(stream, "0x%08x; SUBF32 R%dH,R%dH,R%dH ||MOV32 %s,R%dH", insn32, d, e, f, str, a);
                            break;
                        }

                    }
                    break;
                }
                case 0b0010: //1110 0010 .... ....
                {
                    switch((insn & 0x00f0) >> 4) {
                        case 0b0000:
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //1110 0010 0000 0000 0000 0000 mem32 MOV32 mem32,STF
                                {
                                    if ((insn32 >> 8) & 0xff)
                                    {
                                        length = 4;
                                        uint32_t mem32 = insn32 & 0xff;
                                        get_loc_string(str, mem32, LOC32);
                                        fprintf_func(stream, "0x%08x; MOV32 %s,STF", insn32, str);
                                        break;
                                    }
                                }
                                case 0b0011: //1110 0010 0000 0011 0000 0aaa mem32 MOV32 mem32,RaH
                                {
                                    length = 4;
                                    uint32_t a = (insn32 >> 8) & 0b111;
                                    uint32_t mem32 = insn32 & 0xff;
                                    get_loc_string(str, mem32, LOC32);
                                    fprintf_func(stream, "0x%08x; MOV32 %s,R%dH", insn32, str, a);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0001: //1110 0010 0001 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0011: //1110 0010 0001 0011 0000 0aaa mem16 MOV16 mem16,RaH
                                {
                                    length = 4;
                                    uint32_t a = (insn32 >> 8) & 0b111;
                                    uint32_t mem16 = insn32 & 0xff;
                                    get_loc_string(str, mem16, LOC16);
                                    fprintf_func(stream, "0x%08x; MOV16 %s,R%dH", insn32, str, a);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b0010: //1110 0010 0010 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0011: //1110 0010 0010 0011 0000 0aaa mem32 MOVD32 RaH, mem32
                                {
                                    length = 4;
                                    uint32_t a = (insn32 >> 8) & 0b111;
                                    uint32_t mem32 = insn32 & 0xff;
                                    get_loc_string(str, mem32, LOC32);
                                    fprintf_func(stream, "0x%08x; MOVD32 R%dH,%s", insn32, a, str);
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b1000: //1110 0010 1000 ....
                        {
                            switch (insn & 0xf) {
                                case 0b0000: //1110 0010 1000 0000 0000 0000 mem32 MOV32 STF, mem32
                                {
                                    length = 4;
                                    uint32_t mem32 = insn32 & 0xff;
                                    get_loc_string(str, mem32, LOC32);
                                    fprintf_func(stream, "0x%08x; MOV32 STF,%s", insn32, str);
                                    break;
                                }
                                case 0b0100: //1110 0010 1000 0100 0000 0aaa mem32 UI32TOF32 RaH,mem32
                                {
                                    if (((insn32 & 0xffff) >> 11) == 0)
                                    {
                                        length = 4;
                                        uint32_t mem32 = insn32 & 0xff;
                                        uint32_t a = (insn32 >> 8) & 0b111;
                                        get_loc_string(str, mem32, LOC32);
                                        fprintf_func(stream, "0x%08x; UI32TOF32 R%dH,%s", insn32, a, str);
                                    }
                                    break;
                                }
                                case 0b1000: //1110 0010 1000 1000 0000 0aaa mem32 I32TOF32 RaH,mem32
                                {
                                    if (((insn32 & 0xffff) >> 11) == 0)
                                    {
                                        length = 4;
                                        uint32_t mem32 = insn32 & 0xff;
                                        uint32_t a = (insn32 >> 8) & 0b111;
                                        get_loc_string(str, mem32, LOC32);
                                        fprintf_func(stream, "0x%08x; I32TOF32 R%dH,%s", insn32, a, str);
                                    }
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b1010: //1110 0010 1010 CNDF 0000 0aaa mem32 MOV32 RaH, mem32 {, CNDF} 
                        {
                            length = 4;
                            uint32_t condf = insn & 0xf;
                            uint32_t a = (insn32 >> 8) & 0b111;
                            uint32_t mem32 = insn32 & 0xff;
                            //todo
                            get_loc_string(str, mem32, LOC32);
                            get_condf_string(str2, condf);
                            fprintf_func(stream, "0x%08x; MOV32 R%dH,%s,%s", insn32, a, str, str2);
                            break;
                        }
                        case 0b1100: //1110 0010 1100 ....
                        {
                            if ((insn & 0xf) == 0b1000)//1110 0010 1100 1000 0000 0aaa mem16 I16TOF32 RaH,mem16
                            {
                                if (((insn32 & 0xffff) >> 11) == 0)
                                {
                                    length = 4;
                                    uint32_t mem16 = insn32 & 0xff;
                                    uint32_t a = (insn32 >> 8) & 0b111;
                                    get_loc_string(str, mem16, LOC16);
                                    fprintf_func(stream, "0x%08x; I32TOF32 R%dH,%s", insn32, a, str);
                                }
                            }
                            if ((insn & 0xf) == 0b0100)//1110 0010 1100 0100 0000 0aaa mem16 UI16TOF32 RaH,mem16
                            {
                                if (((insn32 & 0xffff) >> 11) == 0)
                                {
                                    length = 4;
                                    uint32_t mem16 = insn32 & 0xff;
                                    uint32_t a = (insn32 >> 8) & 0b111;
                                    get_loc_string(str, mem16, LOC16);
                                    fprintf_func(stream, "0x%08x; UI32TOF32 R%dH,%s", insn32, a, str);
                                }
                            }
                            break;
                        }
                    }
                    break;
                }
                case 0b0011: //1110 0011 .... ....
                {
                    switch((insn & 0x00f0) >> 4) {
                        case 0b0001: //1110 0011 0001 fffe eedd daaa mem32 ADDF32 RdH, ReH, RfH || MOV32 RaH, mem32
                        {
                            length = 4;
                            uint32_t mem32 = insn32 & 0xff;
                            uint32_t a = (insn32 >> 8) & 0b111;
                            uint32_t d = (insn32 >> 11) & 0b111;
                            uint32_t e = (insn32 >> 14) & 0b111;
                            uint32_t f = (insn >> 1) & 0b111;
                            get_loc_string(str, mem32, LOC32);
                            fprintf_func(stream, "0x%08x; ADDF32 R%dH,R%dH,R%dH ||MOV32 R%dH,%s", insn32, d, e, f, a, str);
                            break;
                        }
                        case 0b0010: //1110 0011 0010 fffe eedd daaa mem32 SUBF32 RdH,ReH,RfH||MOV32 RaH,mem32
                        {
                            length = 4;
                            uint32_t mem32 = insn32 & 0xff;
                            uint32_t a = (insn32 >> 8) & 0b111;
                            uint32_t d = (insn32 >> 11) & 0b111;
                            uint32_t e = (insn32 >> 14) & 0b111;
                            uint32_t f = (insn >> 1) & 0b111;
                            get_loc_string(str, mem32, LOC32);
                            fprintf_func(stream, "0x%08x; SUBF32 R%dH,R%dH,R%dH ||MOV32 R%dH,%s", insn32, d, e, f, a, str);
                            break;
                        }
                    }
                    break;
                }
                case 0b0101: //1110 0101 .... ....
                {
                    switch((insn & 0x00f0) >> 4) {
                        case 0b1010://1110 0101 1010 0aaa CMPF32 RaH,#0.0
                        {
                            if (((insn>>3) & 1) == 0)
                            {
                                uint32_t a = insn & 0b111;
                                fprintf_func(stream, "0x%04x;     CMPF32 R%dH,#0.0", insn, a);
                            }
                            break;
                        }
                    }
                    break;
                }
                case 0b0110: //1110 0110 .... .... 
                {
                    switch ((insn & 0x00c0) >> 6) {
                        case 0b00:
                        {
                            break;
                        }
                        case 0b01://1110 0110 01.. ....
                        {
                            break;
                        }
                        case 0b10://1110 0110 10.. ....
                        {
                            switch (insn & 0x003f) {
                                case 0b001000:
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1000 0000 0000 00bb baaa F32TOI32 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; F32TOI32 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001001: 
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1001 0000 0000 00bb baaa I32TOF32 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; I32OTF32 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001010:
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1010 0000 0000 00bb baaa F32TOUI32 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; F32TOUI32 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001011: 
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1011 0000 0000 00bb baaa UI32TOF32 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; UI32OTF32 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001100: 
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1100 0000 0000 00bb baaa F32TOI16 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; F32TOI16 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    else if(((insn32 & 0xffff) >> 6) == 0b1000000000)//1110 0110 1000 1100 1000 0000 00bb baaa F32TOI16R RaH,RbH
                                    {   
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; F32TOI16R R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001101:
                                {          
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1101 0000 0000 00bb baaa I16TOF32 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; I16TOF32 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001110:
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1110 0000 0000 00bb baaa F32TOUI16 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; F32TOUI16 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    else if(((insn32 & 0xffff) >> 6) == 0b1000000000)//1110 0110 1000 1110 1000 0000 00bb baaa F32TOUI16R RaH,RbH
                                    {   
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; F32TOUI16R R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b001111:
                                {          
                                    if (((insn32 & 0xffff) >> 6) == 0)//1110 0110 1000 1111 0000 0000 00bb baaa UI16TOF32 RaH,RbH
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; UI16TOF32 R%dH,R%dH", insn32, a, b);
                                        length = 4;
                                    }
                                    break;
                                }
                                case 0b010010: //1110 0110 1001 0010 0000 0000 00bb baaa EISQRTF32 RaH,RbH
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)
                                    {
                                        length = 4;
                                        uint32_t a = insn32 & 0b111;
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        fprintf_func(stream, "0x%08x; EISQRTF32 R%dH,R%dH", insn32, a, b);
                                    }
                                    break;
                                }
                                case 0b010011: //1110 0110 1001 0011 0000 0000 00bb baaa EINVF32 RaH,RbH
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)
                                    {
                                        length = 4;
                                        uint32_t a = insn32 & 0b111;
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        fprintf_func(stream, "0x%08x; EINVF32 R%dH,R%dH", insn32, a, b);
                                    }
                                    break;
                                }
                                case 0b010101: //1110 0110 1001 0101 0000 0000 00bb baaa ABSF32 RaH,RbH
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; ABSF32 R%dH,R%dH", insn32, a, b);
                                    }
                                    break;
                                }
                                case 0b010100: //1110 0110 1001 0100 0000 0000 00bb baaa CMPF32 RaH, RbH
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)
                                    {
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        length = 4;
                                        fprintf_func(stream, "0x%08x; CMPF32 R%dH,R%dH", insn32, a, b);
                                    }
                                    break;
                                }
                            }
                            break;
                        }
                        case 0b11://1110 0110 11.. ....
                        {
                            switch(((insn & 0x0030) >> 4))
                            {
                                case 0b00: //1110 0110 1100 CNDF 0000 0000 00bb baaa MOV32 RaH,RbH{,CNDF}
                                {
                                    if (((insn32 & 0xffff) >> 6) == 0)
                                    {
                                        length = 4;
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        uint32_t cndf = insn & 0xf;
                                        get_condf_string(str, cndf);
                                        fprintf_func(stream, "0x%08x; MOV32 R%dH,R%dH,%s", insn32, a, b, str);
                                    }
                                    break;
                                }
                                case 0b01:
                                {
                                    break;
                                }
                                case 0b10:
                                {
                                    break;
                                }
                                case 0b11://1110 0110 1111 ....
                                {
                                    if (((insn&0xf) == 1) && (((insn32 & 0xffff) >> 6) == 0))//1110 0110 1111 0001 0000 0000 00bb baaa
                                    {
                                        length = 4;
                                        uint32_t b = (insn32 >> 3) & 0b111;
                                        uint32_t a = insn32 & 0b111;
                                        fprintf_func(stream, "0x%08x; FRACF32 R%dH,R%dH", insn32, a, b);
                                    }
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    break;
                }
                case 0b0111: //1110 0111 .... ....
                {
                    switch((insn & 0x00f0) >> 4) {
                        case 0b0000: //1110 0111 0000 0000 0000 000c ccbb baaa MPYF32 RaH, RbH, RcH
                        {
                            length = 4;
                            uint32_t a = insn32 & 0b111;
                            uint32_t b = (insn32 >> 3) & 0b111;
                            uint32_t c = (insn32 >> 6) & 0b111;
                            fprintf_func(stream, "0x%08x; MPYF32 R%dH,R%dH,R%dH", insn32, a, b, c);
                            break;
                        }
                        case 0b0001: //1110 0111 0001 0000 0000 000c ccbb baaa ADDF32 RaH, RbH, RcH
                        {
                            length = 4;
                            uint32_t a = insn32 & 0b111;
                            uint32_t b = (insn32 >> 3) & 0b111;
                            uint32_t c = (insn32 >> 6) & 0b111;
                            fprintf_func(stream, "0x%08x; ADDF32 R%dH,R%dH,R%dH", insn32, a, b, c);
                            break;
                        }
                        case 0b0100: //1110 0111 0100 ....
                        {
                            if (((insn32>>18) & 0b11) == 0)//1110 0111 0100 00ff feee dddc ccbb baaa MPYF32 RaH, RbH, RcH || ADDF32 RdH, ReH, RfH 
                            {
                                length = 4;
                                uint32_t a = insn32 & 0b111;
                                uint32_t b = (insn32 >> 3) & 0b111;
                                uint32_t c = (insn32 >> 6) & 0b111;
                                uint32_t d = (insn32 >> 9) & 0b111;
                                uint32_t e = (insn32 >> 12) & 0b111;
                                uint32_t f = (insn32 >> 15) & 0b111;
                                fprintf_func(stream, "0x%08x; MPYF32 R%dH,R%dH,R%dH || ADDF32 R%dH,R%dH,R%dH", insn32, a, b, c, d, e, f);
                            }
                            break;
                        }
                        case 0b0010: //1110 0111 0010 ....
                        {
                            if (((insn32 >> 9) & 0b11111111111) == 0)//1110 0111 0010 0000 0000 000c ccbb baaa SUBF32 RaH,RbH,RcH
                            {
                                length = 4;
                                uint32_t c = (insn32 >> 6) & 0b111;
                                uint32_t b = (insn32 >> 3) & 0b111;
                                uint32_t a = insn32 & 0b111;
                                fprintf_func(stream, "0x%08x; SUBF32 R%dH,R%dH,R%dH", insn32, a, b, c);
                            }
                            break;
                        }
                    }
                    break;
                }
                case 0b1000: //1110 1000 .... .....
                {
                    switch((insn & 0x00f0) >> 4) {
                        case 0b0000: //1110 1000 0000 xIII IIII IIII IIII Iaaa MOVIZ RaH, #16FHiHex/ MOVXI RaH, #16FLoHex
                        {
                            length = 4;
                            uint32_t a = insn32 & 0b111;
                            uint32_t imm = (insn32 >> 3) & 0xffff;
                            if (((insn >> 3) & 1) == 0)
                            {
                                fprintf_func(stream, "0x%08x; MOVIZ R%dH,#0x%4x", insn32, a, imm);
                            }
                            else
                            {
                                fprintf_func(stream, "0x%08x; MOVXI R%dH,#0x%4x", insn32, a, imm);
                            }
                            break;
                        }
                        case 0b0001: //1110 1000 0001 0III IIII IIII IIII Iaaa CMPF32 RaH,#16FHi
                        {
                            if (((insn >> 3) & 1) == 0)
                            {
                                length = 4;
                                uint32_t imm = (insn32 >> 3) & 0xffff;
                                uint32_t a = insn32 & 0b111;
                                fprintf_func(stream, "0x%08x; CMPF32 R%dH,#0x%4x", insn32, a, imm);
                            }
                            break;
                        }
                        case 0b0100:
                        case 0b0101:
                        case 0b0110:
                        case 0b0111: //1110 1000 10II IIII IIII IIII IIbb baaa MPYF32 RaH, #16FHi, RbH
                        {
                            length = 4;
                            uint32_t hi = (insn32>>6) & 0xffff;
                            uint32_t a = insn32 & 0b111;
                            uint32_t b = (insn32 >> 3) & 0b111;
                            fprintf_func(stream, "0x%08x; MPYF32 R%dH,#0x%4x,R%dH", insn32, a, hi, b);
                            break;
                        }
                        case 0b1000:
                        case 0b1001:
                        case 0b1010:
                        case 0b1011: //1110 1000 10II IIII IIII IIII IIbb baaa ADDF32 RaH, #16FHi, RbH
                        {
                            length = 4;
                            uint32_t hi = (insn32>>6) & 0xffff;
                            uint32_t a = insn32 & 0b111;
                            uint32_t b = (insn32 >> 3) & 0b111;
                            fprintf_func(stream, "0x%08x; ADDF32 R%dH,#0x%4x,R%dH", insn32, a, hi, b);
                            break;
                        }
                        case 0b1100:
                        case 0b1101:
                        case 0b1110:
                        case 0b1111://1110 1000 11II IIII IIII IIII IIbb baaa SUBF32 RaH,#16FHi,RbH
                        {
                            length = 4;
                            uint32_t hi = (insn32>>6) & 0xffff;
                            uint32_t a = insn32 & 0b111;
                            uint32_t b = (insn32 >> 3) & 0b111;
                            fprintf_func(stream, "0x%08x; SUBF32 R%dH,#0x%4x,R%dH", insn32, a, hi, b);
                            break;
                        }
                    }
                    break;
                }
                case 0b1010://1110 1010 LLLL LLLL SUBR loc16,AL
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     SUBR %s,AL", insn, str);
                    break;
                }
                case 0b1011://1110 1011 LLLL LLLL SUBR loc16,AH
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     SUBR %s,AH", insn, str);
                    break;
                }
                case 0b1100://1110 1100 CCCC CCCC SBF 8bitoffset,EQ
                {
                    int8_t offset = (insn & 0xff);
                    fprintf_func(stream, "0x%04x;     SBF %d,EQ", insn, offset);
                    break;
                }
                case 0b1101://1110 1101 CCCC CCCC SBF 8bitoffset,NEQ
                {
                    int8_t offset = (insn & 0xff);
                    fprintf_func(stream, "0x%04x;     SBF %d,NEQ", insn, offset);
                    break;
                }
                case 0b1110://1110 1110 CCCC CCCC SBF 8bitoffset,TC
                {
                    int8_t offset = (insn & 0xff);
                    fprintf_func(stream, "0x%04x;     SBF %d,TC", insn, offset);
                    break;
                }
                case 0b1111://1110 1111 CCCC CCCC SBF 8bitoffset,NTC
                {
                    int8_t offset = (insn & 0xff);
                    fprintf_func(stream, "0x%04x;     SBF %d,NTC", insn, offset);
                    break;
                }
            }
            break;
        }
        case 0b1111:
            switch ((insn & 0x0f00) >> 8) {
                case 0b0000: //1111 0000 CCCC CCCC XORB AL,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     XORB AL,#0x%x", insn, imm);
                    break;
                }
                case 0b0001: //1111 0001 CCCC CCCC XORB AH,#8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     XORB AH,#0x%x", insn, imm);
                    break;
                }
                case 0b0010: //1111 0010 LLLL LLLL XOR loc16,AL
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     XOR %s,AL", insn, str);
                    break;
                }
                case 0b0011: //1111 0011 LLLL LLLL XOR loc16,AH
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%04x;     XOR %s,AH", insn, str);
                    break;
                }
                case 0b0100: //1111 0100 LLLL LLLL 32bit MOV *(0:16bit),loc16
                {
                    uint32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%08x; MOV *(0:#0x%04x),%s", insn32, imm, str);
                    length = 4;
                    break;
                }
                case 0b0101: //1111 0101 LLLL LLLL 32bit MOV loc16,*(0:16bit)
                {
                    uint32_t imm = insn32 & 0xffff;
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str, mode, LOC16);
                    fprintf_func(stream, "0x%08x; MOV %s,*(0:#0x%04x)", insn32, str, imm);
                    length = 4;
                    break;
                }
                case 0b0110: //1111 0110 CCCC CCCC RPT #8bit
                {
                    uint32_t value = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     RPT #%d", insn,value);
                    break;
                }
                case 0b0111: //1111 0111 CCCC CCCC RPT loc16
                {
                    uint32_t mode = insn & 0xff;
                    get_loc_string(str,mode,LOC16);
                    fprintf_func(stream, "0x%04x;     RPT %s", insn, str);
                    break;
                }
                case 0b1000: //1111 10CC CCCC CCCC MOV DP,#10bit
                case 0b1001: //1111 10CC CCCC CCCC MOV DP,#10bit
                case 0b1010: //1111 10CC CCCC CCCC MOV DP,#10bit
                case 0b1011: //1111 10CC CCCC CCCC MOV DP,#10bit
                {
                    uint32_t imm = insn & 0x3ff;// 10bit
                    fprintf_func(stream, "0x%04x;     MOV DP,#0x%x", insn, imm);
                    break;
                }
                case 0b1100: //1111 1100 IIII IIII ADRK #8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     ADRK #%d", insn, imm);
                    break;
                }
                case 0b1101: //1111 1101 CCCC CCCC SBRK #8bit
                {
                    uint32_t imm = insn & 0xff;
                    fprintf_func(stream, "0x%04x;     SBRK #%d", insn, imm);
                    break;
                }
                case 0b1110: //1111 1110 .... ....
                {
                    uint32_t imm = insn & 0x7f;
                    if (((insn & 0xff) >> 7) == 0) { //ADDB SP,#7bit
                        fprintf_func(stream, "0x%04x;     ADDB SP,#%d", insn, imm);
                    }
                    else { //SUBB SP,#7bit
                        fprintf_func(stream, "0x%04x;     SUBB SP,#%d", insn, imm);
                    }
                    break;
                }
                case 0b1111://1111 1111 .... ....
                    switch ((insn & 0x00f0) >> 4) {
                        case 0b0000: //1111 1111 0000 SHFT CCCC CCCC CCCC CCCC SUB ACC,#16bit<<#0...15
                        {
                            int32_t imm = insn32 & 0xffff;
                            uint32_t shift = insn & 0x000f;
                            fprintf_func(stream, "0x%08x; SUB ACC, #0x%x<<#%d", insn32, imm, shift);
                            length = 4;
                            break;
                        }
                        case 0b0001: //1111 1111 0001 SHFT 32bit ADD ACC, #16bit<#0...15
                        {
                            int32_t imm = insn32 & 0xffff;
                            uint32_t shift = insn & 0x000f;
                            fprintf_func(stream, "0x%08x; ADD ACC, #0x%x<<#%d", insn32, imm, shift);
                            length = 4;
                            break;
                        }
                        case 0b0010: /*1111 1111 0010 .... MOV ACC, #16bit<#0...15 */
                        {
                            uint32_t imm = insn32 & 0xffff;
                            uint32_t shift = insn & 0x000f;
                            fprintf_func(stream, "0x%08x; MOV ACC, #0x%x<<#%d", insn32, imm, shift);
                            length = 4;
                            break;
                        }
                        case 0b0011: //1111 1111 0011 SHFT LSL ACC,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%08x; LSL ACC,#%d", insn, shift);
                            break;
                        }
                        case 0b0100: //1111 1111 0100 SHFT SFR ACC,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%08x; SFR ACC,#%d", insn, shift);
                            break;
                        }
                        case 0b0101: //1111 1111 0101 ....
                            switch (insn & 0xf) {
                                case 0b0000: //1111 1111 0101 0000 LSL ACC,T
                                {
                                    fprintf_func(stream, "0x%04x;     LSL ACC,T", insn);
                                    break;
                                }
                                case 0b0001: //1111 1111 0101 0001 SFR ACC,T
                                {
                                    fprintf_func(stream, "0x%04x;     SFR ACC,T", insn);
                                    break;
                                }
                                case 0b0010: //1111 1111 0101 0011 ROR ACC
                                {
                                    fprintf_func(stream, "0x%04x;     ROR ACC", insn);
                                    break;
                                }
                                case 0b0011: //1111 1111 0101 0011 ROL ACC
                                {
                                    fprintf_func(stream, "0x%04x;     ROL ACC", insn);
                                    break;
                                }
                                case 0b0100: //1111 1111 0101 0100 NEG ACC
                                {
                                    fprintf_func(stream, "0x%04x;     NEG ACC", insn);
                                    break;
                                }
                                case 0b0101: //1111 1111 0101 0101 NOT ACC
                                {
                                    fprintf_func(stream, "0x%04x;     NOT ACC", insn);
                                    break;
                                }
                                case 0b0110: //1111 1111 0101 0110 ABS ACC
                                {
                                    fprintf_func(stream, "0x%04x;     ABS ACC", insn);
                                    break;
                                }
                                case 0b0111: //1111 1111 0101 0111 SAT ACC
                                {
                                    fprintf_func(stream, "0x%04x;     SAT ACC", insn);
                                    break;
                                }
                                case 0b1000: //1111 1111 0101 1000 TEST ACC
                                {
                                    fprintf_func(stream, "0x%04x;     TEST ACC", insn);
                                    break;
                                }
                                case 0b1001: //1111 1111 0101 1001 CMPL ACC,P<<PM
                                {
                                    fprintf_func(stream, "0x%04x;     CMPL ACC,P<<PM", insn);
                                    break;
                                }
                                case 0b1010: //1111 1111 0101 1010 MOVL P,ACC
                                {
                                    fprintf_func(stream, "0x%04x;     MOVL P,ACC", insn);
                                    break;
                                }
                                case 0b1100: //1111 1111 0101 1100 NEG AL
                                {
                                    fprintf_func(stream, "0x%04x;     NEG AL", insn);
                                    break;
                                }
                                case 0b1101: //1111 1111 0101 1101 NEG AH
                                {
                                    fprintf_func(stream, "0x%04x;     NEG AH", insn);
                                    break;
                                }
                                case 0b1110: //1111 1111 0101 1110 NOT AL
                                {
                                    fprintf_func(stream, "0x%04x;     NOT AL", insn);
                                    break;
                                }
                                case 0b1111: //1111 1111 0101 1111 NOT AH
                                {
                                    fprintf_func(stream, "0x%04x;     NOT AH", insn);
                                    break;
                                }
                            }
                            break;
                        case 0b0110: //1111 1111 0110 ....
                        {
                            if (((insn & 0xf) >> 3) == 1) { //1111 1111 0110 1... SPM shift
                                int32_t pm = insn & 0b111;
                                pm = 1 - pm;
                                if (pm > 0) {
                                    fprintf_func(stream, "0x%04x;     SPF +%d", insn, pm);
                                }
                                else {
                                    fprintf_func(stream, "0x%04x;     SPF %d", insn, pm);
                                }
                            }
                            else {//1111 1111 0110 0...
                                switch(insn & 0b111) {
                                    case 0b010://1111 1111 0110 0010 LSR AL,T
                                    {
                                        fprintf_func(stream, "0x%04x;     LSR AL,T", insn);
                                        break;
                                    }
                                    case 0b011://1111 1111 0110 0010 LSR AH,T
                                    {
                                        fprintf_func(stream, "0x%04x;     LSR AH,T", insn);
                                        break;
                                    }
                                    case 0b100:
                                    {
                                        fprintf_func(stream, "0x%04x;     ASR AL,T", insn);
                                        break;
                                    }
                                    case 0b101: //1111 1111 0110 010a ASR AX,T
                                    {
                                        fprintf_func(stream, "0x%04x;     ASR AH,T", insn);
                                        break;
                                    }
                                    case 0b110:
                                    {
                                        fprintf_func(stream, "0x%04x;     LSL AL,T", insn);
                                        break;
                                    }
                                    case 0b111: //1111 1111 0110 011a LSL AX,T
                                    {
                                        fprintf_func(stream, "0x%04x;     LSL AH,T", insn);
                                        break;
                                    }
                                }
                            }
                            break;
                        }
                        case 0b0111: //1111 1111 0111 1/0nnn NORM ACC,XARn++/--
                        {
                            uint32_t type = insn & 0xf;
                            if (type < 8)
                            {
                                fprintf_func(stream, "0x%04x;     NORM ACC,XAR%d--", insn, type);
                            }
                            else
                            {
                                fprintf_func(stream, "0x%04x;     NORM ACC,XAR%d++", insn, type-8);
                            }
                            break;
                        }
                        case 0b1000:
                        case 0b1001://1111 1111 100A SHFT LSL AX,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            uint32_t is_AH = (insn >> 4) & 1;
                            if (is_AH) 
                                fprintf_func(stream, "0x%04x;     LSL AH,#%d", insn, shift);
                            else
                                fprintf_func(stream, "0x%04x;     LSL AL,#%d", insn, shift);
                            break;
                        }
                        case 0b1010:
                        case 0b1011://1111 1111 101A SHFT ASR AX,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            uint32_t is_AH = (insn >> 4) & 1;
                            if (is_AH) 
                                fprintf_func(stream, "0x%04x;     ASR AH,#%d", insn, shift);
                            else
                                fprintf_func(stream, "0x%04x;     ASR AL,#%d", insn, shift);
                            break;
                        }
                        case 0b1100://1111 1111 1100 SHFT LSR AL,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%04x;     LSR AL,#%d", insn, shift);
                            break;
                        }
                        case 0b1101://1111 1111 1101 SHFT LSR AH,#1...16
                        {
                            uint32_t shift = (insn & 0xf) + 1;
                            fprintf_func(stream, "0x%04x;     LSR AH,#%d", insn, shift);
                            break;
                        }
                        case 0b1110://1111 1111 1110 COND CCCC CCCC CCCC CCCC B 16bitOffset,COND
                        {
                            int16_t imm = (uint16_t)(insn32&0xffff);
                            int16_t cond = (insn&0xf);
                            length = 4;
                            get_cond_string(str, cond);
                            fprintf_func(stream, "0x%08x; B #%d,%s", insn32, imm, str);
                            break;
                        }
                    }
                    break;
            }
            break;
    }



    return length;
}