/*
 * TMS320C28x helper defines
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, see <http://www.gnu.org/licenses/>.
 */

//interrupt
DEF_HELPER_1(aborti, void, env)
// DEF_HELPER_FLAGS_3(exception, TCG_CALL_NO_WG, void, env, i32)
DEF_HELPER_3(exception, void, env, i32, i32)

// addressing mode
DEF_HELPER_3(addressing_mode, i32, env, i32, i32)
DEF_HELPER_2(ld_loc16_byte_addressing, i32, env, i32) //value = low 8bit
DEF_HELPER_2(ld_loc16, i32, env, i32)//value = low 16bit
DEF_HELPER_2(ld_loc32, i32, env, i32)
DEF_HELPER_3(st_loc16, void, env, i32, i32)
DEF_HELPER_3(st_loc16_byte_addressing, void, env, i32, i32) //value = low 8bit
DEF_HELPER_3(st_loc32, void, env, i32, i32)

DEF_HELPER_2(test_cond,i32,env,i32) //affect V bit
DEF_HELPER_2(test_condf,i32,env,i32) //not affect stf

// test bit
DEF_HELPER_2(test_N_Z_16, void, env, i32)
DEF_HELPER_2(test_N_Z_32, void, env, i32)
DEF_HELPER_3(test_N_Z_64, void, env, i32, i32)

DEF_HELPER_4(test_C_V_16, void, env, i32, i32, i32)
DEF_HELPER_4(test_sub_C_V_16, void, env, i32, i32, i32)
DEF_HELPER_4(test_C_32, void, env, i32, i32, i32)
DEF_HELPER_4(test_C_32_shift16, void, env, i32, i32, i32)
DEF_HELPER_4(test_sub_C_32, void, env, i32, i32, i32)
DEF_HELPER_4(test_sub_C_32_shift16, void, env, i32, i32, i32)
DEF_HELPER_4(test_V_32, void, env, i32, i32, i32)
DEF_HELPER_4(test_sub_V_32, void, env, i32, i32, i32)
DEF_HELPER_4(test_C_V_32, void, env, i32, i32, i32)
DEF_HELPER_4(test_sub_C_V_32, void, env, i32, i32, i32)
DEF_HELPER_4(test_OVC_OVM_32, void, env, i32, i32, i32) // affect acc value
DEF_HELPER_4(test_sub_OVC_OVM_32, void, env, i32, i32, i32) // affect acc value
DEF_HELPER_5(test2_C_V_OVC_OVM_32, void, env, i32, i32, i32, i32) // affect acc value, used for 3 op add
DEF_HELPER_5(test2_sub_C_V_OVC_OVM_32, void, env, i32, i32, i32, i32) // affect acc value, used for 3 op add
DEF_HELPER_4(test_OVCU_32, void, env, i32, i32, i32) // affect acc value
DEF_HELPER_4(test_sub_OVCU_32, void, env, i32, i32, i32) // affect acc value

//load
DEF_HELPER_2(extend_low_sxm, i32, env, i32) //signed extend value with sxm
DEF_HELPER_1(ld_xar_arp, i32, env) //load XAR[arp]

//store
DEF_HELPER_2(st_xar_arp, void, env, i32) //store XAR[arp]

//print
DEF_HELPER_1(print, void, i32)
DEF_HELPER_1(print_env, void, env)

//math
DEF_HELPER_1(abs_acc, void, env)
DEF_HELPER_1(abstc_acc, void, env)
DEF_HELPER_2(subcu, void, env, i32)
DEF_HELPER_2(subcul, void, env, i32)

//mov
DEF_HELPER_4(mov_16bit_loc16, void, env, i32, i32, i32)// MOV *(0:16bit),loc16
DEF_HELPER_4(mov_loc16_16bit, void, env, i32, i32, i32)// MOV loc16,*(0:16bit)

//cmp
DEF_HELPER_3(cmp16_N_Z_C, void, env, i32, i32)// cmp, test N,Z,C
DEF_HELPER_3(cmp32_N_Z_C, void, env, i32, i32)// cmpl, test N,Z,C
DEF_HELPER_1(cmp64_acc_p, void, env)// cmp64 ACC:P, test N,Z,V
DEF_HELPER_2(maxcul_p_loc32, void, env, i32)// maxcul p,loc32
DEF_HELPER_2(mincul_p_loc32, void, env, i32)// mincul p,loc32

//fpu
DEF_HELPER_2(fpu_absf, i32, env, i32)// ret = absf(value), modify NF,ZF bit
DEF_HELPER_3(fpu_addf, i32, env, i32, i32)// ret = addf(a,b), modify LUF,LVF bit
DEF_HELPER_3(fpu_cmpf, void, env, i32, i32)// cmp(a,b), modify NF,ZF bit
DEF_HELPER_2(fpu_einvf, i32, env, i32)// einvf32, modify LUF,LVF bit
DEF_HELPER_2(fpu_eisqrtf, i32, env, i32)// einvf32, modify LUF,LVF bit
DEF_HELPER_2(fpu_fracf32, i32, env, i32)//RaH = fractional portion of Rbh
DEF_HELPER_3(fpu_subf, i32, env, i32, i32)// ret = subf(a,b), modify LUF,LVF bit
DEF_HELPER_3(fpu_mpyf, i32, env, i32, i32)// ret = mpyf(a,b), modify LUF,LVF bit

//fpu max/min
DEF_HELPER_3(fpu_maxf, i32, env, i32, i32)// a = max(a,b), modify ZF,NF bit
DEF_HELPER_3(fpu_minf, i32, env, i32, i32)// a = max(a,b), modify ZF,NF bit

//fpu float to integer
DEF_HELPER_2(fpu_f32toi16, i32, env, i32)//ret = f32toi16(a), round to zero
DEF_HELPER_2(fpu_f32toi16r, i32, env, i32)//ret = f32toi16(a), round to nearest even
DEF_HELPER_2(fpu_f32toi32, i32, env, i32)//ret = f32toi32(a), round to zero
DEF_HELPER_2(fpu_f32toui16, i32, env, i32)//ret = f32toui16(a), round to zero
DEF_HELPER_2(fpu_f32toui16r, i32, env, i32)//ret = f32toui16(a), round to nearest even
DEF_HELPER_2(fpu_f32toui32, i32, env, i32)//ret = f32toui32(a), round to zero
//fpu integer to float
DEF_HELPER_2(fpu_i16tof32, i32, env, i32)
DEF_HELPER_2(fpu_ui16tof32, i32, env, i32)
DEF_HELPER_2(fpu_i32tof32, i32, env, i32)//ret = i32tof32(a)
DEF_HELPER_2(fpu_ui32tof32, i32, env, i32)

