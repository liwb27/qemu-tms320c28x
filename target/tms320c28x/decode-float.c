// ABSF32 RaH,RbH
static void gen_absf32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_absf(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

// ADDF32 RaH, #16FHi, RbH
static void gen_addf32_rah_16fhi_rbh(DisasContext *ctx, uint32_t a, uint32_t hi, uint32_t b)
{
    TCGv tmp = cpu_tmp[0];
    hi = hi << 16;
    tcg_gen_movi_i32(tmp, hi);
    gen_helper_fpu_addf(cpu_rh[a], cpu_env, tmp, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

// ADDF32 RaH, RbH, RcH
static void gen_addf32_rah_rbh_rch(DisasContext *ctx, uint32_t a, uint32_t b, uint32_t c)
{
    gen_helper_fpu_addf(cpu_rh[a], cpu_env, cpu_rh[b], cpu_rh[c]);
    gen_sync_fpu_mem(a);
}

//ADDF32 RdH, ReH, RfH || MOV32 mem32, RaH
static void gen_addf32_rdh_reh_rfh_mov32_mem32_rah(DisasContext *ctx, uint32_t d, uint32_t e, uint32_t f, uint32_t mem32, uint32_t a)
{
    if(is_reg_addressing_mode(mem32, LOC32))
    {
        return;
    }
    //save to mem32 first
    gen_st_loc32(mem32, cpu_rh[a]);
    //add
    gen_helper_fpu_addf(cpu_rh[d], cpu_env, cpu_rh[e], cpu_rh[f]);
    gen_sync_fpu_mem(d);
}

//ADDF32 RdH, ReH, RfH || MOV32  RaH,mem32
static void gen_addf32_rdh_reh_rfh_mov32_rah_mem32(DisasContext *ctx, uint32_t d, uint32_t e, uint32_t f, uint32_t a, uint32_t mem32)
{
    if(is_reg_addressing_mode(mem32, LOC32))
    {
        return;
    }
    //save to rah
    gen_ld_loc32(cpu_rh[a], mem32);
    gen_test_nf_ni_zf_zi(cpu_rh[a]);
    gen_sync_fpu_mem(a);
    //add
    gen_helper_fpu_addf(cpu_rh[d], cpu_env, cpu_rh[e], cpu_rh[f]);
    gen_sync_fpu_mem(d);
}

//CMPF32 RaH, RbH
static void gen_cmpf32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_cmpf(cpu_env, cpu_rh[a], cpu_rh[b]);
}

//CMPF32 RaH, #16FHi
static void gen_cmpf32_rah_16fhi(DisasContext *ctx, uint32_t a, uint32_t imm)
{
    TCGv tmp = cpu_tmp[0];
    tcg_gen_movi_i32(tmp, imm << 16);
    gen_helper_fpu_cmpf(cpu_env, cpu_rh[a], tmp);
}

//CMPF32 RaH, #0.0
static void gen_cmpf32_rah_0(DisasContext *ctx, uint32_t a)
{
    TCGv tmp = cpu_tmp[0];
    tcg_gen_movi_i32(tmp, 0);
    gen_helper_fpu_cmpf(cpu_env, cpu_rh[a], tmp);
}

//EINVF32 RaH,RbH
static void gen_einvf32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_einvf(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//EISQRTF32 RaH,RbH
static void gen_eisqrtf32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_eisqrtf(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//F32TOI16 RaH,RbH
static void gen_f32toi16_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_f32toi16(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//F32TOI16R RaH,RbH
static void gen_f32toi16r_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_f32toi16r(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//F32TOI32 RaH,RbH
static void gen_f32toi32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_f32toi32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//F32TOUI16 RaH,RbH
static void gen_f32toui16_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_f32toui16(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//F32TOUI16R RaH,RbH
static void gen_f32toui16r_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_f32toui16r(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//F32TOUI32 RaH,RbH
static void gen_f32toui32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_f32toui32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//FRACF32 RaH,RbH
static void gen_fracf32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_fracf32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//I16TOF32 RaH,RbH
static void gen_i16tof32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_i16tof32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//I16TOF32 RaH,mem16
static void gen_i16tof32_rah_mem16(DisasContext *ctx, uint32_t a, uint32_t mem16)
{
    TCGv tmp = cpu_tmp[0];
    gen_ld_loc16(tmp, mem16);
    gen_helper_fpu_i16tof32(cpu_rh[a], cpu_env, tmp);
    gen_sync_fpu_mem(a);
}

//I32TOF32 RaH,mem32
static void gen_i32tof32_rah_mem32(DisasContext *ctx, uint32_t a, uint32_t mem32)
{
    TCGv tmp = cpu_tmp[0];
    gen_ld_loc32(tmp, mem32);
    gen_helper_fpu_i32tof32(cpu_rh[a], cpu_env, tmp);
    gen_sync_fpu_mem(a);
}

//I32TOF32 RaH,RbH
static void gen_i32tof32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_i32tof32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//UI32TOF32 RaH,RbH
static void gen_ui32tof32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_ui32tof32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

// MOV16 mem16, RaH
static void gen_mov16_mem16_rah(DisasContext *ctx, uint32_t mem16, uint32_t a)
{
    // TCGv tmp = cpu_tmp[0];
    // tcg_gen_andi_i32(tmp, cpu_rh[a], 0xffff);

    if(!is_reg_addressing_mode(mem16, LOC16))
    {
        // get low 16bit in helper func
        gen_st_loc16(mem16, cpu_rh[a]);
    }
}

// MOV32 *(0:16bitAddr), loc32
// MOV32 RaH, ACC
// MOV32 RaH, XARn
// MOV32 RaH, XT
static void gen_mov32_addr16_loc32(DisasContext *ctx, uint32_t addr, uint32_t mode)
{
    TCGv tmp = cpu_tmp[0];
    TCGv addr_param = cpu_tmp[1];
    tcg_gen_movi_i32(addr_param, addr);
    //load loc32
    gen_ld_loc32(tmp, mode);
    //mask stf
    gen_mask_st32_stf(addr, tmp);
    //store
    gen_st32u_swap(tmp, addr_param);
    //sync mem with fpu register
    gen_sync_fpu_reg(addr, tmp);
}

// MOV32 loc32, *(0:16bitAddr)
// MOV32 ACC, RaH
// MOV32 XARn, RaH
// MOV32 XT, RaH
static void gen_mov32_loc32_addr16(DisasContext *ctx, uint32_t addr, uint32_t loc32)
{
    TCGv tmp = cpu_tmp[0];
    TCGv addr_param = cpu_tmp[0];
    tcg_gen_movi_i32(addr_param, addr);
    gen_ld32u_swap(tmp, addr_param);
    gen_st_loc32(loc32, tmp);
    gen_test_acc_N_Z(loc32);
}

// MOV32 mem32, RaH
static void gen_mov32_mem32_rah(DisasContext *ctx, uint32_t mem32, uint32_t a)
{
    if(!is_reg_addressing_mode(mem32, LOC32))
    {
        gen_st_loc32(mem32, cpu_rh[a]);
    }
}

// MOV32 mem32, STF
static void gen_mov32_mem32_stf(DisasContext *ctx, uint32_t mem32)
{
    if(!is_reg_addressing_mode(mem32, LOC32))
    {
        gen_st_loc32(mem32, cpu_stf);
    }
}

// MOV32 RaH, mem32 {, CNDF} 
static void gen_mov32_rah_mem32_cndf(DisasContext *ctx, uint32_t a, uint32_t mem32, uint32_t cndf)
{
    if(is_reg_addressing_mode(mem32, LOC32))
    {
        return;
    }
    TCGLabel *done = gen_new_label();
    TCGv test = cpu_tmp[0];
    TCGv condf_tcg = cpu_tmp[1];
    TCGv tmp = cpu_tmp[2];
    tcg_gen_movi_i32(condf_tcg, cndf);
    //if (CNDF == TRUE)
    gen_helper_test_condf(test, cpu_env, condf_tcg);
    gen_ld_loc32(tmp, mem32);//if test false, still do the load
    tcg_gen_brcondi_i32(TCG_COND_EQ, test, 0, done);
    //if (CNDF == TRUE) RaH = [mem32]
    tcg_gen_mov_i32(cpu_rh[a], tmp);//save to reg
    gen_sync_fpu_mem(a);//save to mem
    //modify stf
    if (cndf == 15) //UNCF
    {
        gen_test_nf_ni_zf_zi(cpu_rh[a]);
    }
    gen_set_label(done);
}

// MOV32 RaH, RbH {, CNDF}
static void gen_mov32_rah_rbh_cndf(DisasContext *ctx, uint32_t a, uint32_t b, uint32_t cndf)
{
    TCGLabel *done = gen_new_label();
    TCGv test = cpu_tmp[0];
    TCGv condf_tcg = cpu_tmp[1];
    tcg_gen_movi_i32(condf_tcg, cndf);
    //if (CNDF == TRUE)
    gen_helper_test_condf(test, cpu_env, condf_tcg);
    tcg_gen_brcondi_i32(TCG_COND_EQ, test, 0, done);
    //if (CNDF == TRUE) RaH = RbH
    tcg_gen_mov_i32(cpu_rh[a], cpu_rh[b]);//save to reg
    //modify stf
    if (cndf == 15) //UNCF
    {
        gen_test_nf_ni_zf_zi(cpu_rh[a]);
    }
    gen_set_label(done);
}

// MOV32 STF, mem32
static void gen_mov32_stf_mem32(DisasContext *ctx, uint32_t mem32)
{
    if(is_reg_addressing_mode(mem32, LOC32))
    {
        return;
    }
    TCGv tmp = cpu_tmp[2];
    gen_ld_loc32(tmp, mem32);
    //SHDWS bit is not affected by loading the status register either from memory or from the shadow values.
    tcg_gen_andi_i32(tmp, tmp, 0x0000027f);//mask stf
    tcg_gen_mov_i32(cpu_stf, tmp);//save to reg
    gen_sync_fpu_mem(-1);//save to mem
}

// MOVD32 RaH, mem32
static void gen_movd32_rah_mem32(DisasContext *ctx, uint32_t a, uint32_t mem32)
{
    if(is_reg_addressing_mode(mem32, LOC32))
    {
        return;
    }
    TCGv addr = cpu_tmp[0];
    gen_get_loc_addr(addr, mem32, LOC32); //get mem32
    // RaH = [mem32]
    gen_ld32u_swap(cpu_rh[a], addr);
    gen_sync_fpu_mem(a);
    // [mem32+2] = [mem32]
    tcg_gen_addi_i32(addr, addr, 2);
    gen_st32u_swap(cpu_rh[a], addr);
    //modify stf
    gen_test_nf_ni_zf_zi(cpu_rh[a]);
}


// MOVIZ RaH, #16FHiHex
static void gen_moviz_rah_16fhihex(DisasContext *ctx, uint32_t a, uint32_t hi)
{
    //RaH[31:16] = #16FHiHex
    //RaH[15:0] = 0
    hi = hi << 16;
    tcg_gen_movi_i32(cpu_rh[a], hi);
    gen_sync_fpu_mem(a);
}

// MOVST0 FLAG
static void gen_movst0_flag(DisasContext *ctx, uint32_t flag)
{
    TCGv tmp = cpu_tmp[0];
    TCGv tmp2 = cpu_tmp[1];
    flag = flag & 0xff;
    uint32_t lvf,luf,nf,ni,zf,zi;
    lvf = (flag>>0) & 1;
    luf = (flag>>1) & 1;
    if (lvf == 1)
    {
        if (luf == 1)//luf,lvf
        {
            gen_get_bit(tmp, cpu_stf, LUF_BIT, LUF_MASK);
            gen_get_bit(tmp2, cpu_stf, LVF_BIT, LVF_MASK);
            tcg_gen_or_i32(tmp, tmp, tmp2);//luf=1 or lvf=1, set st0_v=1
            gen_set_bit(cpu_st0, V_BIT, V_MASK, tmp);
            gen_seti_bit(cpu_stf, LUF_BIT, LUF_MASK, 0);//clear luf,lvf
            gen_seti_bit(cpu_stf, LVF_BIT, LVF_MASK, 0);
        }
        else//lvf
        {
            gen_get_bit(tmp, cpu_stf, LVF_BIT, LVF_MASK);
            gen_set_bit(cpu_st0, V_BIT, V_MASK, tmp);
            gen_seti_bit(cpu_stf, LVF_BIT, LVF_MASK, 0);//clear lvf
        }
    }
    else 
    {
        if (luf == 1)//luf
        {
            gen_get_bit(tmp, cpu_stf, LUF_BIT, LUF_MASK);
            gen_set_bit(cpu_st0, V_BIT, V_MASK, tmp);
            gen_seti_bit(cpu_stf, LUF_BIT, LUF_MASK, 0);//clear luf
        }
    }
    nf = (flag>>2) & 1;
    ni = (flag>>3) & 1;
    if (nf == 1)
    {
        if (ni == 1)//NF,NI
        {
            gen_get_bit(tmp, cpu_stf, NF_BIT, NF_MASK);
            gen_get_bit(tmp2, cpu_stf, NI_BIT, NI_MASK);
            tcg_gen_or_i32(tmp, tmp, tmp2);//NF=1 or NI=1, set st0_N=1
            gen_set_bit(cpu_st0, N_BIT, N_MASK, tmp);
        }
        else//NF
        {
            gen_get_bit(tmp, cpu_stf, NF_BIT, NF_MASK);
            gen_set_bit(cpu_st0, N_BIT, N_MASK, tmp);
        }
    }
    else 
    {
        if (ni == 1)//NI
        {
            gen_get_bit(tmp, cpu_stf, NI_BIT, NI_MASK);
            gen_set_bit(cpu_st0, N_BIT, N_MASK, tmp);
        }
    }
    zf = (flag>>4) & 1;
    zi = (flag>>5) & 1;
    if (zf == 1)
    {
        if (zi == 1)//ZF,ZI
        {
            gen_get_bit(tmp, cpu_stf, ZF_BIT, ZF_MASK);
            gen_get_bit(tmp2, cpu_stf, ZI_BIT, ZI_MASK);
            tcg_gen_or_i32(tmp, tmp, tmp2);//ZF=1 or ZI=1, set st0_Z=1
            gen_set_bit(cpu_st0, Z_BIT, Z_MASK, tmp);
        }
        else//ZF
        {
            gen_get_bit(tmp, cpu_stf, ZF_BIT, ZF_MASK);
            gen_set_bit(cpu_st0, Z_BIT, Z_MASK, tmp);
        }
    }
    else 
    {
        if (zi == 1)//ZI
        {
            gen_get_bit(tmp, cpu_stf, ZI_BIT, ZI_MASK);
            gen_set_bit(cpu_st0, Z_BIT, Z_MASK, tmp);
        }
    }
    if (((flag>>6) & 1) == 1) { //CI
        gen_get_bit(tmp, cpu_stf, TF_BIT, TF_MASK);
        gen_set_bit(cpu_st0, C_BIT, C_MASK, tmp);
    }
    if (((flag>>7) & 1) == 1) { //TF
        gen_get_bit(tmp, cpu_stf, TF_BIT, TF_MASK);
        gen_set_bit(cpu_st0, TC_BIT, TC_MASK, tmp);
    }
}

// MOVXI RaH, #16FHiHex
static void gen_movxi_rah_16flohex(DisasContext *ctx, uint32_t a, uint32_t lo)
{
    // RaH[15:0] = #16FLoHex
    // RaH[31:16] = Unchanged
    lo = lo & 0xffff;
    tcg_gen_andi_i32(cpu_rh[a], cpu_rh[a], 0xffff0000);
    tcg_gen_ori_i32(cpu_rh[a], cpu_rh[a], lo);
    gen_sync_fpu_mem(a);
}

// SUBF32 RaH, RbH, RcH
static void gen_subf32_rah_rbh_rch(DisasContext *ctx, uint32_t a, uint32_t b, uint32_t c)
{
    gen_helper_fpu_subf(cpu_rh[a], cpu_env, cpu_rh[b], cpu_rh[c]);
    gen_sync_fpu_mem(a);
}

// SUBF32 RaH, #16FHi, RbH
static void gen_subf32_rah_16fhi_rbh(DisasContext *ctx, uint32_t a, uint32_t hi, uint32_t b)
{
    TCGv tmp = cpu_tmp[0];
    hi = hi << 16;
    tcg_gen_movi_i32(tmp, hi);
    gen_helper_fpu_subf(cpu_rh[a], cpu_env, tmp, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//SUBF32 RdH, ReH, RfH || MOV32  RaH,mem32
static void gen_subf32_rdh_reh_rfh_mov32_rah_mem32(DisasContext *ctx, uint32_t d, uint32_t e, uint32_t f, uint32_t a, uint32_t mem32)
{
    if(is_reg_addressing_mode(mem32, LOC32))
    {
        return;
    }
    //save to rah
    gen_ld_loc32(cpu_rh[a], mem32);
    gen_test_nf_ni_zf_zi(cpu_rh[a]);
    gen_sync_fpu_mem(a);
    //add
    gen_helper_fpu_addf(cpu_rh[d], cpu_env, cpu_rh[e], cpu_rh[f]);
    gen_sync_fpu_mem(d);
}

//UI16TOF32 RaH,mem16
static void gen_ui16tof32_rah_mem16(DisasContext *ctx, uint32_t a, uint32_t mem16)
{
    TCGv tmp = cpu_tmp[0];
    gen_ld_loc16(tmp, mem16);
    gen_helper_fpu_ui16tof32(cpu_rh[a], cpu_env, tmp);
    gen_sync_fpu_mem(a);
}

//UI16TOF32 RaH,RbH
static void gen_ui16tof32_rah_rbh(DisasContext *ctx, uint32_t a, uint32_t b)
{
    gen_helper_fpu_ui16tof32(cpu_rh[a], cpu_env, cpu_rh[b]);
    gen_sync_fpu_mem(a);
}

//UI32TOF32 RaH,mem32
static void gen_ui32tof32_rah_mem32(DisasContext *ctx, uint32_t a, uint32_t mem32)
{
    TCGv tmp = cpu_tmp[0];
    gen_ld_loc32(tmp, mem32);
    gen_helper_fpu_ui32tof32(cpu_rh[a], cpu_env, tmp);
    gen_sync_fpu_mem(a);
}