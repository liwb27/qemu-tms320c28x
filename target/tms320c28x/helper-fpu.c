#include "qemu/osdep.h"
#include "cpu.h"
#include "exec/helper-proto.h"
#include "exception.h"
#include "fpu/softfloat.h"

uint32_t HELPER(fpu_absf)(CPUTms320c28xState *env, uint32_t value)
{
    float32 ret = float32_abs(make_float32(value));
    //NF=0
    cpu_set_stf(env, 0, NF_BIT, NF_MASK);
    // ZF = 0;
    // if ( RaH[30:23] == 0) ZF = 1;
    if((ret & 0x7f800000) == 0)
    {
        cpu_set_stf(env, 1, ZF_BIT, ZF_MASK);
    }
    else
    {
        cpu_set_stf(env, 0, ZF_BIT, ZF_MASK);
    }
    return ret;
}

uint32_t HELPER(fpu_addf)(CPUTms320c28xState *env, uint32_t a, uint32_t b)
{
    //set round mode
    if (cpu_get_stf(env, RND32_BIT, RND32_MASK) == 0) //truncate
    {
        env->fp_status.float_rounding_mode = float_round_to_zero;
    }
    else
    {
        env->fp_status.float_rounding_mode = float_round_nearest_even;
    }
    float32 ret = float32_add(a, b, &env->fp_status);
    //LVF
    if (env->fp_status.float_exception_flags & float_flag_overflow)
    {
        cpu_set_stf(env, 1, LVF_BIT, LVF_MASK);
    }
    //LUF
    if (env->fp_status.float_exception_flags & float_flag_underflow)
    {
        cpu_set_stf(env, 1, LUF_BIT, LUF_MASK);
    }
    return ret;
}

void HELPER(fpu_cmpf)(CPUTms320c28xState *env, uint32_t a, uint32_t b)
{
    //Negative zero will be treated as positive zero. ---------checked
    //A denormalized value will be treated as positive zero.
    env->fp_status.flush_inputs_to_zero = 1;
    //Not-a-Number (NaN) will be treated as infinity.
    if (a == 0x7FBFFFFF)
        a = 0x7F800000;
    if (b == 0x7FBFFFFF)
        b = 0x7F800000;
    //compare
    int flag = float32_compare(a, b, &env->fp_status);
    switch(flag)
    {
        case float_relation_less://If(RaH < RbH) {ZF=0, NF=1}
            cpu_set_stf(env, 0, ZF_BIT, ZF_MASK);
            cpu_set_stf(env, 1, NF_BIT, NF_MASK);
            break;
        case float_relation_equal://If(RaH == RbH) {ZF=1, NF=0}
            cpu_set_stf(env, 1, ZF_BIT, ZF_MASK);
            cpu_set_stf(env, 0, NF_BIT, NF_MASK);
            break;
        case float_relation_greater://If(RaH > RbH) {ZF=0, NF=0}
            cpu_set_stf(env, 0, ZF_BIT, ZF_MASK);
            cpu_set_stf(env, 0, NF_BIT, NF_MASK);
            break;
        case float_relation_unordered:
        default:
            g_assert_not_reached();
    }
    //restore status
    env->fp_status.flush_inputs_to_zero = 0;
}

uint32_t HELPER(fpu_einvf)(CPUTms320c28xState *env, uint32_t value)
{
    uint32_t luf, lvf;
    luf = 0;
    if(((value >> 25) & 0b111111) == 0b111111) //[30:25] = 6'b1
    {
        if (((value >> 23) & 0b11) == 0b11) //[24:23] = 2'b1
        {
            luf = 1;
        }
        else if(((value >> 23) & 0b11) == 0b10) //[24:23] = b10
        {
            luf = 1;
        }
        else if ((((value >> 23) & 0b11) == 0b01) && (((value >> 16) & 0b1111111) != 0))//[24:23] = b01 & [22:16]!=0
        {
            luf = 1;
        }
    }
    lvf = 0;
    if (((value >> 23) & 0xff) == 0)
    {
        lvf = 1;
    }
    cpu_set_stf(env, luf, LUF_BIT, LUF_MASK);
    cpu_set_stf(env, lvf, LVF_BIT, LVF_MASK);
    uint32_t D;
    switch((value >> 16) & 0b1111111)
    {
        case 0x00: D = 0b0000000; break;
        case 0x01: D = 0b1111110; break;
        case 0x02: D = 0b1111100; break;
        case 0x03: D = 0b1111010; break;
        case 0x04: D = 0b1111000; break;
        case 0x05: D = 0b1110110; break;
        case 0x06: D = 0b1110101; break;
        case 0x07: D = 0b1110011; break;
        case 0x08: D = 0b1110001; break;
        case 0x09: D = 0b1101111; break;
        case 0x0a: D = 0b1101101; break;
        case 0x0b: D = 0b1101100; break;
        case 0x0c: D = 0b1101010; break;
        case 0x0d: D = 0b1101000; break;
        case 0x0e: D = 0b1100111; break;
        case 0x0f: D = 0b1100101; break;
        case 0x10: D = 0b1100100; break;
        case 0x11: D = 0b1100010; break;
        case 0x12: D = 0b1100000; break;
        case 0x13: D = 0b1011111; break;
        case 0x14: D = 0b1011101; break;
        case 0x15: D = 0b1011100; break;
        case 0x16: D = 0b1011010; break;
        case 0x17: D = 0b1011001; break;
        case 0x18: D = 0b1011000; break;
        case 0x19: D = 0b1010110; break;
        case 0x1a: D = 0b1010101; break;
        case 0x1b: D = 0b1010011; break;
        case 0x1c: D = 0b1010010; break;
        case 0x1d: D = 0b1010001; break;
        case 0x1e: D = 0b1001111; break;
        case 0x1f: D = 0b1001110; break;
        case 0x20: D = 0b1001101; break;
        case 0x21: D = 0b1001100; break;
        case 0x22: D = 0b1001010; break;
        case 0x23: D = 0b1001001; break;
        case 0x24: D = 0b1001000; break;
        case 0x25: D = 0b1000111; break;
        case 0x26: D = 0b1000101; break;
        case 0x27: D = 0b1000100; break;
        case 0x28: D = 0b1000011; break;
        case 0x29: D = 0b1000010; break;
        case 0x2a: D = 0b1000001; break;
        case 0x2b: D = 0b1000000; break;
        case 0x2c: D = 0b0111111; break;
        case 0x2d: D = 0b0111101; break;
        case 0x2e: D = 0b0111100; break;
        case 0x2f: D = 0b0111011; break;
        case 0x30: D = 0b0111010; break;
        case 0x31: D = 0b0111001; break;
        case 0x32: D = 0b0111000; break;
        case 0x33: D = 0b0110111; break;
        case 0x34: D = 0b0110110; break;
        case 0x35: D = 0b0110101; break;
        case 0x36: D = 0b0110100; break;
        case 0x37: D = 0b0110011; break;
        case 0x38: D = 0b0110010; break;
        case 0x39: D = 0b0110001; break;
        case 0x3a: D = 0b0110000; break;
        case 0x3b: D = 0b0101111; break;
        case 0x3c: D = 0b0101110; break;
        case 0x3d: D = 0b0101101; break;
        case 0x3e: D = 0b0101100; break;
        case 0x3f: D = 0b0101100; break;
        case 0x40: D = 0b0101011; break;
        case 0x41: D = 0b0101010; break;
        case 0x42: D = 0b0101001; break;
        case 0x43: D = 0b0101000; break;
        case 0x44: D = 0b0100111; break;
        case 0x45: D = 0b0100110; break;
        case 0x46: D = 0b0100101; break;
        case 0x47: D = 0b0100101; break;
        case 0x48: D = 0b0100100; break;
        case 0x49: D = 0b0100011; break;
        case 0x4a: D = 0b0100010; break;
        case 0x4b: D = 0b0100001; break;
        case 0x4c: D = 0b0100001; break;
        case 0x4d: D = 0b0100000; break;
        case 0x4e: D = 0b0011111; break;
        case 0x4f: D = 0b0011110; break;
        case 0x50: D = 0b0011110; break;
        case 0x51: D = 0b0011101; break;
        case 0x52: D = 0b0011100; break;
        case 0x53: D = 0b0011011; break;
        case 0x54: D = 0b0011011; break;
        case 0x55: D = 0b0011010; break;
        case 0x56: D = 0b0011001; break;
        case 0x57: D = 0b0011000; break;
        case 0x58: D = 0b0011000; break;
        case 0x59: D = 0b0010111; break;
        case 0x5a: D = 0b0010110; break;
        case 0x5b: D = 0b0010110; break;
        case 0x5c: D = 0b0010101; break;
        case 0x5d: D = 0b0010100; break;
        case 0x5e: D = 0b0010100; break;
        case 0x5f: D = 0b0010011; break;
        case 0x60: D = 0b0010010; break;
        case 0x61: D = 0b0010010; break;
        case 0x62: D = 0b0010001; break;
        case 0x63: D = 0b0010000; break;
        case 0x64: D = 0b0010000; break;
        case 0x65: D = 0b0001111; break;
        case 0x66: D = 0b0001110; break;
        case 0x67: D = 0b0001110; break;
        case 0x68: D = 0b0001101; break;
        case 0x69: D = 0b0001101; break;
        case 0x6a: D = 0b0001100; break;
        case 0x6b: D = 0b0001011; break;
        case 0x6c: D = 0b0001011; break;
        case 0x6d: D = 0b0001010; break;
        case 0x6e: D = 0b0001010; break;
        case 0x6f: D = 0b0001001; break;
        case 0x70: D = 0b0001001; break;
        case 0x71: D = 0b0001000; break;
        case 0x72: D = 0b0000111; break;
        case 0x73: D = 0b0000111; break;
        case 0x74: D = 0b0000110; break;
        case 0x75: D = 0b0000110; break;
        case 0x76: D = 0b0000101; break;
        case 0x77: D = 0b0000101; break;
        case 0x78: D = 0b0000100; break;
        case 0x79: D = 0b0000100; break;
        case 0x7a: D = 0b0000011; break;
        case 0x7b: D = 0b0000011; break;
        case 0x7c: D = 0b0000010; break;
        case 0x7d: D = 0b0000010; break;
        case 0x7e: D = 0b0000001; break;
        case 0x7f: D = 0b0000001; break;
    }

    uint32_t result = 0;
    //   Result[31:23] = 9'h0;
    //   Result[22:0] = {(LUF || LVF)? 7'b0: OP_EISQRTF32? ISQRT_SEED: OP_EINVF32? INV_SEED: 7'b0, 16'b0};
    if (luf == 1 || lvf == 1)
    {
        //Result[22:0]=7'b0
    }
    else
    {
        result = (result & 0xff80ffff) | (D<<16);
    }
    
    //     if(isunderflow) begin
    if (luf == 1)
    {
    //       Result[31:23] = 9'h0;
    }
    else
    {
    //       if(Src_DenormalizedorZero) begin
        if (float32_is_zero_or_denormal(value))
        {
    //         Result[31:23] = 9'h0ff;
            result = (result & 0x7fffff) | 0x7f800000;
        }
        else//Src_Normalized && !Src_Zero && !isunderflow
        {
    //         Result[31] = Src[31];
            result = (result & 0x7fffffff) | (value & 0x80000000);
    //         Result[30:23] = {Adder_12_Result_INV[7:0]}; 
    //         {11'h7ff,!(|Src[22:16])} + {4'b0,~Src[30:23]}
            uint32_t tmp = (((value >> 16) & 0b1111111) == 0) ? 0xfff : 0xffe;
            uint32_t tmp2 = (~((value >> 23) & 0xff)) & 0xfff;
            result = (result & 0x807fffff) | ((tmp + tmp2) & 0xff) << 23;
        }
    }
    return result;
}

uint32_t HELPER(fpu_eisqrtf)(CPUTms320c28xState *env, uint32_t value)
{
    uint32_t luf,lvf;
    luf = 0;
    if(((value >> 23) & 0xff) == 0xff) //[30:23] = 8'b1
    {
        if ((value >> 31) == 0) //[31] == 0
        {
            luf = 1;
        }
    }
    lvf = 0;
    if (float32_is_zero_or_denormal(value) || ((value >> 31) == 1)) //(Src_DenormalizedorZero || Src[31])
    {
        lvf = 1;
    }
    cpu_set_stf(env, luf, LUF_BIT, LUF_MASK);
    cpu_set_stf(env, lvf, LVF_BIT, LVF_MASK);
    uint32_t D;
    switch ((value >> 16) & 0xff)
    {
        case 0x00: D = 0b0110101; break;
        case 0x01: D = 0b0110100; break;
        case 0x02: D = 0b0110100; break;
        case 0x03: D = 0b0110011; break;
        case 0x04: D = 0b0110010; break;
        case 0x05: D = 0b0110010; break;
        case 0x06: D = 0b0110001; break;
        case 0x07: D = 0b0110000; break;
        case 0x08: D = 0b0110000; break;
        case 0x09: D = 0b0101111; break;
        case 0x0a: D = 0b0101110; break;
        case 0x0b: D = 0b0101110; break;
        case 0x0c: D = 0b0101101; break;
        case 0x0d: D = 0b0101100; break;
        case 0x0e: D = 0b0101100; break;
        case 0x0f: D = 0b0101011; break;
        case 0x10: D = 0b0101011; break;
        case 0x11: D = 0b0101010; break;
        case 0x12: D = 0b0101001; break;
        case 0x13: D = 0b0101001; break;
        case 0x14: D = 0b0101000; break;
        case 0x15: D = 0b0101000; break;
        case 0x16: D = 0b0100111; break;
        case 0x17: D = 0b0100111; break;
        case 0x18: D = 0b0100110; break;
        case 0x19: D = 0b0100110; break;
        case 0x1a: D = 0b0100101; break;
        case 0x1b: D = 0b0100100; break;
        case 0x1c: D = 0b0100100; break;
        case 0x1d: D = 0b0100011; break;
        case 0x1e: D = 0b0100011; break;
        case 0x1f: D = 0b0100010; break;
        case 0x20: D = 0b0100010; break;
        case 0x21: D = 0b0100001; break;
        case 0x22: D = 0b0100001; break;
        case 0x23: D = 0b0100000; break;
        case 0x24: D = 0b0100000; break;
        case 0x25: D = 0b0011111; break;
        case 0x26: D = 0b0011111; break;
        case 0x27: D = 0b0011110; break;
        case 0x28: D = 0b0011110; break;
        case 0x29: D = 0b0011110; break;
        case 0x2a: D = 0b0011101; break;
        case 0x2b: D = 0b0011101; break;
        case 0x2c: D = 0b0011100; break;
        case 0x2d: D = 0b0011100; break;
        case 0x2e: D = 0b0011011; break;
        case 0x2f: D = 0b0011011; break;
        case 0x30: D = 0b0011010; break;
        case 0x31: D = 0b0011010; break;
        case 0x32: D = 0b0011010; break;
        case 0x33: D = 0b0011001; break;
        case 0x34: D = 0b0011001; break;
        case 0x35: D = 0b0011000; break;
        case 0x36: D = 0b0011000; break;
        case 0x37: D = 0b0010111; break;
        case 0x38: D = 0b0010111; break;
        case 0x39: D = 0b0010111; break;
        case 0x3a: D = 0b0010110; break;
        case 0x3b: D = 0b0010110; break;
        case 0x3c: D = 0b0010101; break;
        case 0x3d: D = 0b0010101; break;
        case 0x3e: D = 0b0010101; break;
        case 0x3f: D = 0b0010100; break;
        case 0x40: D = 0b0010100; break;
        case 0x41: D = 0b0010011; break;
        case 0x42: D = 0b0010011; break;
        case 0x43: D = 0b0010011; break;
        case 0x44: D = 0b0010010; break;
        case 0x45: D = 0b0010010; break;
        case 0x46: D = 0b0010010; break;
        case 0x47: D = 0b0010001; break;
        case 0x48: D = 0b0010001; break;
        case 0x49: D = 0b0010000; break;
        case 0x4a: D = 0b0010000; break;
        case 0x4b: D = 0b0010000; break;
        case 0x4c: D = 0b0001111; break;
        case 0x4d: D = 0b0001111; break;
        case 0x4e: D = 0b0001111; break;
        case 0x4f: D = 0b0001110; break;
        case 0x50: D = 0b0001110; break;
        case 0x51: D = 0b0001110; break;
        case 0x52: D = 0b0001101; break;
        case 0x53: D = 0b0001101; break;
        case 0x54: D = 0b0001101; break;
        case 0x55: D = 0b0001100; break;
        case 0x56: D = 0b0001100; break;
        case 0x57: D = 0b0001100; break;
        case 0x58: D = 0b0001011; break;
        case 0x59: D = 0b0001011; break;
        case 0x5a: D = 0b0001011; break;
        case 0x5b: D = 0b0001010; break;
        case 0x5c: D = 0b0001010; break;
        case 0x5d: D = 0b0001010; break;
        case 0x5e: D = 0b0001001; break;
        case 0x5f: D = 0b0001001; break;
        case 0x60: D = 0b0001001; break;
        case 0x61: D = 0b0001001; break;
        case 0x62: D = 0b0001000; break;
        case 0x63: D = 0b0001000; break;
        case 0x64: D = 0b0001000; break;
        case 0x65: D = 0b0000111; break;
        case 0x66: D = 0b0000111; break;
        case 0x67: D = 0b0000111; break;
        case 0x68: D = 0b0000110; break;
        case 0x69: D = 0b0000110; break;
        case 0x6a: D = 0b0000110; break;
        case 0x6b: D = 0b0000110; break;
        case 0x6c: D = 0b0000101; break;
        case 0x6d: D = 0b0000101; break;
        case 0x6e: D = 0b0000101; break;
        case 0x6f: D = 0b0000100; break;
        case 0x70: D = 0b0000100; break;
        case 0x71: D = 0b0000100; break;
        case 0x72: D = 0b0000100; break;
        case 0x73: D = 0b0000011; break;
        case 0x74: D = 0b0000011; break;
        case 0x75: D = 0b0000011; break;
        case 0x76: D = 0b0000011; break;
        case 0x77: D = 0b0000010; break;
        case 0x78: D = 0b0000010; break;
        case 0x79: D = 0b0000010; break;
        case 0x7a: D = 0b0000010; break;
        case 0x7b: D = 0b0000001; break;
        case 0x7c: D = 0b0000001; break;
        case 0x7d: D = 0b0000001; break;
        case 0x7e: D = 0b0000001; break;
        case 0x7f: D = 0b0000000; break;
        case 0x80: D = 0b0000000; break;
        case 0x81: D = 0b1111111; break;
        case 0x82: D = 0b1111110; break;
        case 0x83: D = 0b1111101; break;
        case 0x84: D = 0b1111100; break;
        case 0x85: D = 0b1111011; break;
        case 0x86: D = 0b1111010; break;
        case 0x87: D = 0b1111001; break;
        case 0x88: D = 0b1111000; break;
        case 0x89: D = 0b1110111; break;
        case 0x8a: D = 0b1110111; break;
        case 0x8b: D = 0b1110110; break;
        case 0x8c: D = 0b1110101; break;
        case 0x8d: D = 0b1110100; break;
        case 0x8e: D = 0b1110011; break;
        case 0x8f: D = 0b1110010; break;
        case 0x90: D = 0b1110001; break;
        case 0x91: D = 0b1110001; break;
        case 0x92: D = 0b1110000; break;
        case 0x93: D = 0b1101111; break;
        case 0x94: D = 0b1101110; break;
        case 0x95: D = 0b1101101; break;
        case 0x96: D = 0b1101100; break;
        case 0x97: D = 0b1101100; break;
        case 0x98: D = 0b1101011; break;
        case 0x99: D = 0b1101010; break;
        case 0x9a: D = 0b1101001; break;
        case 0x9b: D = 0b1101001; break;
        case 0x9c: D = 0b1101000; break;
        case 0x9d: D = 0b1100111; break;
        case 0x9e: D = 0b1100110; break;
        case 0x9f: D = 0b1100110; break;
        case 0xa0: D = 0b1100101; break;
        case 0xa1: D = 0b1100100; break;
        case 0xa2: D = 0b1100100; break;
        case 0xa3: D = 0b1100011; break;
        case 0xa4: D = 0b1100010; break;
        case 0xa5: D = 0b1100001; break;
        case 0xa6: D = 0b1100001; break;
        case 0xa7: D = 0b1100000; break;
        case 0xa8: D = 0b1011111; break;
        case 0xa9: D = 0b1011111; break;
        case 0xaa: D = 0b1011110; break;
        case 0xab: D = 0b1011101; break;
        case 0xac: D = 0b1011101; break;
        case 0xad: D = 0b1011100; break;
        case 0xae: D = 0b1011100; break;
        case 0xaf: D = 0b1011011; break;
        case 0xb0: D = 0b1011010; break;
        case 0xb1: D = 0b1011010; break;
        case 0xb2: D = 0b1011001; break;
        case 0xb3: D = 0b1011000; break;
        case 0xb4: D = 0b1011000; break;
        case 0xb5: D = 0b1010111; break;
        case 0xb6: D = 0b1010111; break;
        case 0xb7: D = 0b1010110; break;
        case 0xb8: D = 0b1010110; break;
        case 0xb9: D = 0b1010101; break;
        case 0xba: D = 0b1010100; break;
        case 0xbb: D = 0b1010100; break;
        case 0xbc: D = 0b1010011; break;
        case 0xbd: D = 0b1010011; break;
        case 0xbe: D = 0b1010010; break;
        case 0xbf: D = 0b1010010; break;
        case 0xc0: D = 0b1010001; break;
        case 0xc1: D = 0b1010000; break;
        case 0xc2: D = 0b1010000; break;
        case 0xc3: D = 0b1001111; break;
        case 0xc4: D = 0b1001111; break;
        case 0xc5: D = 0b1001110; break;
        case 0xc6: D = 0b1001110; break;
        case 0xc7: D = 0b1001101; break;
        case 0xc8: D = 0b1001101; break;
        case 0xc9: D = 0b1001100; break;
        case 0xca: D = 0b1001100; break;
        case 0xcb: D = 0b1001011; break;
        case 0xcc: D = 0b1001011; break;
        case 0xcd: D = 0b1001010; break;
        case 0xce: D = 0b1001010; break;
        case 0xcf: D = 0b1001001; break;
        case 0xd0: D = 0b1001001; break;
        case 0xd1: D = 0b1001000; break;
        case 0xd2: D = 0b1001000; break;
        case 0xd3: D = 0b1000111; break;
        case 0xd4: D = 0b1000111; break;
        case 0xd5: D = 0b1000110; break;
        case 0xd6: D = 0b1000110; break;
        case 0xd7: D = 0b1000110; break;
        case 0xd8: D = 0b1000101; break;
        case 0xd9: D = 0b1000101; break;
        case 0xda: D = 0b1000100; break;
        case 0xdb: D = 0b1000100; break;
        case 0xdc: D = 0b1000011; break;
        case 0xdd: D = 0b1000011; break;
        case 0xde: D = 0b1000010; break;
        case 0xdf: D = 0b1000010; break;
        case 0xe0: D = 0b1000010; break;
        case 0xe1: D = 0b1000001; break;
        case 0xe2: D = 0b1000001; break;
        case 0xe3: D = 0b1000000; break;
        case 0xe4: D = 0b1000000; break;
        case 0xe5: D = 0b0111111; break;
        case 0xe6: D = 0b0111111; break;
        case 0xe7: D = 0b0111111; break;
        case 0xe8: D = 0b0111110; break;
        case 0xe9: D = 0b0111110; break;
        case 0xea: D = 0b0111101; break;
        case 0xeb: D = 0b0111101; break;
        case 0xec: D = 0b0111101; break;
        case 0xed: D = 0b0111100; break;
        case 0xee: D = 0b0111100; break;
        case 0xef: D = 0b0111011; break;
        case 0xf0: D = 0b0111011; break;
        case 0xf1: D = 0b0111011; break;
        case 0xf2: D = 0b0111010; break;
        case 0xf3: D = 0b0111010; break;
        case 0xf4: D = 0b0111001; break;
        case 0xf5: D = 0b0111001; break;
        case 0xf6: D = 0b0111001; break;
        case 0xf7: D = 0b0111000; break;
        case 0xf8: D = 0b0111000; break;
        case 0xf9: D = 0b0111000; break;
        case 0xfa: D = 0b0110111; break;
        case 0xfb: D = 0b0110111; break;
        case 0xfc: D = 0b0110110; break;
        case 0xfd: D = 0b0110110; break;
        case 0xfe: D = 0b0110110; break;
        case 0xff: D = 0b0110101; break; 
    }
    uint32_t result = 0;
    if (luf == 1 || lvf == 1)
    {
        //Result[22:0]=7'b0
    }
    else
    {
        result = (result & 0xff80ffff) | (D<<16);
    }
    if (float32_is_zero_or_denormal(value))
    {
        //Result[31:23] = 9'h0ff;
        result = (result & 0x7fffff) | 0x7f800000;
    }
    else
    {
        if (float32_is_infinity(value) || float32_is_any_nan(value) || float32_is_neg(value)) //Src_Infinity || Src_NaN || negative
        {
            // Result[31:23] = 9'h0;
            result = (result & 0x7fffff);
        }
        else //positive normalized
        {
            // Result[31] = 1'b0;
            result = result & 0x7fffffff;
            //when the lsb of exponent is 1,and the 8 msb of the fraction is full zero,
            //the new exponent=127-(e-127)/2;
            //else if the the lsb of the exponent is 1, then the new exponent=126-(e-127)/2;
            //else the new exponent=126-(e-128)/2;
            // wire [11:0] A_ISQRT= {3'h0,6'b101111,!(|({!Src[23],Src[22:16]}))|!Src[23]?2'b10:2'b01,Src[23]};                            
            uint32_t tmp = 0b000101111000;
            tmp = tmp | ((value >> 23) & 1);
            if (((~(value >> 23) & 1)<<8 | ((value >> 16) & 0b1111111)) == 0 || (((value >> 23) & 1) == 0))
            {
                tmp = tmp | 0b100;
            }
            else
            {
                tmp = tmp | 0b010;
            }
            // wire [11:0] B_ISQRT= {3'b0,OP_EISQRTF32,~Src[30:23]};
            uint32_t tmp2 = (~(value >> 23) & 0xff) | 0b0001<<8;
            // Result[30:23] = {Adder_12_Result_ISQRT[8:1]};
            uint32_t tmp3 = ((tmp + tmp2 + 1) >> 1) & 0xff;
            result = (result & 0x807fffff) | (tmp3 << 23);
        }
    }
    return result;
}

uint32_t HELPER(fpu_fracf32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    float frac = *(float *)&value - (int)*(float *)&value;
    return *(uint32_t *)&frac;
}

uint32_t HELPER(fpu_subf)(CPUTms320c28xState *env, uint32_t a, uint32_t b)
{
    //set round mode
    if (cpu_get_stf(env, RND32_BIT, RND32_MASK) == 0) //truncate
    {
        env->fp_status.float_rounding_mode = float_round_to_zero;
    }
    else
    {
        env->fp_status.float_rounding_mode = float_round_nearest_even;
    }
    float32 ret = float32_sub(a, b, &env->fp_status);
    //LVF
    if (env->fp_status.float_exception_flags & float_flag_overflow)
    {
        cpu_set_stf(env, 1, LVF_BIT, LVF_MASK);
    }
    //LUF
    if (env->fp_status.float_exception_flags & float_flag_underflow)
    {
        cpu_set_stf(env, 1, LUF_BIT, LUF_MASK);
    }
    return ret;
}

uint32_t HELPER(fpu_f32toi16)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    int ret = float32_to_int16(value, &env->fp_status);
    return ret;
}

uint32_t HELPER(fpu_f32toi16r)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_nearest_even;
    int ret = float32_to_int16(value, &env->fp_status);
    return ret;
}

uint32_t HELPER(fpu_f32toi32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    int ret = float32_to_int32(value, &env->fp_status);
    return ret;
}

uint32_t HELPER(fpu_f32toui16)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    int ret = float32_to_uint16(value, &env->fp_status);
    return ret;
}

uint32_t HELPER(fpu_f32toui16r)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_nearest_even;
    int ret = float32_to_uint16(value, &env->fp_status);
    return ret;
}

uint32_t HELPER(fpu_f32toui32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    int ret = float32_to_uint32(value, &env->fp_status);
    return ret;
}

uint32_t HELPER(fpu_i16tof32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    return int16_to_float32(value, &env->fp_status);
}

uint32_t HELPER(fpu_ui16tof32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    return uint16_to_float32(value, &env->fp_status);
}

uint32_t HELPER(fpu_i32tof32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    return int32_to_float32(value, &env->fp_status);
}

uint32_t HELPER(fpu_ui32tof32)(CPUTms320c28xState *env, uint32_t value)
{
    env->fp_status.float_rounding_mode = float_round_to_zero;
    return uint32_to_float32(value, &env->fp_status);
}