
static void gen_intr(DisasContext *ctx, uint32_t int_n) {
    gen_reset_rptc(ctx);

    if (int_n > 0 && int_n < 15) { // n == 1... 14
        gen_exception(ctx, int_n + 100);
        // TCGv tmp = tcg_const_i32(int_n);
        // TCGv pc = tcg_const_i32(ctx->base.pc_next >> 1);// current pc
        // gen_helper_intr(cpu_env, tmp, pc);
        // tcg_temp_free(tmp);
        // tcg_temp_free(pc);

        ctx->base.is_jmp = DISAS_JUMP;
    }
    else
    {
        //exception
    }
    
}

static void gen_iret(DisasContext *ctx)
{
    TCGv tmp = tcg_temp_new_i32();
    TCGv tmp2 = tcg_temp_new_i32();
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //PC = [SP]
    gen_ld32u_swap(cpu_pc, cpu_sp);
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //DBGSTAT:IER = [SP]
    gen_ld32u_swap(tmp, cpu_sp);
    tcg_gen_shri_i32(cpu_dbgier, tmp, 16);
    tcg_gen_andi_i32(cpu_ier, tmp, 0xffff);
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //DP:ST1 = [SP]
    gen_ld32u_swap(tmp, cpu_sp);
    tcg_gen_shri_i32(cpu_dp, tmp, 16);
    tcg_gen_andi_i32(cpu_st1, tmp, 0xffff);
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //AR1:AR0 = [SP]
    gen_ld32u_swap(tmp, cpu_sp);
    tcg_gen_shri_i32(tmp2, tmp, 16);
    gen_st_reg_low_half(cpu_xar[0], tmp);
    gen_st_reg_low_half(cpu_xar[1], tmp2);
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //PH:PL = [SP]
    gen_ld32u_swap(cpu_p, cpu_sp);
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //AH:AL = [SP]
    gen_ld32u_swap(cpu_acc, cpu_sp);
    //SP = SP - 2
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 2);
    //T:ST0 = [SP]
    gen_ld32u_swap(tmp, cpu_sp);
    tcg_gen_andi_i32(cpu_st0, tmp, 0xffff);
    tcg_gen_shri_i32(tmp, tmp, 16);
    gen_st_reg_high_half(cpu_xt, tmp);
    //SP = SP - 1
    tcg_gen_subi_i32(cpu_sp, cpu_sp, 1);

    // when interrupt is exits, restore RAS to RA bit
    gen_get_bit(tmp, cpu_rb, RAS_BIT, RAS_MASK);
    gen_set_bit(cpu_rb, RA_BIT, RAS_MASK, tmp);
    //SAVE RNDF32=1 (in TMS320C28x FPU Primer (Rev. A))
    gen_seti_bit(cpu_stf, RND32_BIT, RND32_MASK, 1);

    tcg_temp_free(tmp);
    tcg_temp_free(tmp2);
    ctx->base.is_jmp = DISAS_JUMP;
}

static void gen_trap(DisasContext *ctx, uint32_t int_n)
{
    gen_reset_rptc(ctx);

    if (int_n <= 31) { //
        gen_exception(ctx, int_n);
        ctx->base.is_jmp = DISAS_JUMP;
    }
    else
    {
        //exception
    }
}