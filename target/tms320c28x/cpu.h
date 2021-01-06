/*
 * QEMU TMS320C28X CPU
 *
 * Copyright (c) 2019 Liwb <bolia@yeah.net>
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

#ifndef TMS320C28X_CPU_H
#define TMS320C28X_CPU_H

#include "exec/cpu-defs.h"
#include "hw/core/cpu.h"

struct Tms320c28xCPU;

#define TYPE_TMS320C28X_CPU "tms320c28x-cpu"

//用来类型转换的宏
#define TMS320C28X_CPU_CLASS(klass) \
    OBJECT_CLASS_CHECK(Tms320c28xCPUClass, (klass), TYPE_TMS320C28X_CPU)
#define TMS320C28X_CPU(obj) \
    OBJECT_CHECK(Tms320c28xCPU, (obj), TYPE_TMS320C28X_CPU)
#define TMS320C28X_CPU_GET_CLASS(obj) \
    OBJECT_GET_CLASS(Tms320c28xCPUClass, (obj), TYPE_TMS320C28X_CPU)

/**
 * Tms320c28xCPUClass:
 * @parent_realize: The parent class' realize handler.
 * @parent_reset: The parent class' reset handler.
 *
 * A Tms320c28x CPU model.
 */
typedef struct Tms320c28xCPUClass {
    /*< private >*/
    CPUClass parent_class;
    /*< public >*/

    DeviceRealize parent_realize;
    void (*parent_reset)(CPUState *cpu);
} Tms320c28xCPUClass;

// used in tcg.h
#define TARGET_INSN_START_EXTRA_WORDS 1

// mmu mode
enum {
    MMU_NOMMU_IDX = 0,
    MMU_SUPERVISOR_IDX = 1,
    MMU_USER_IDX = 2,
};

/* Interrupt */
#define NR_IRQS  32 //max irq number

/* st0 */
enum ST0_BIT {
    SXM_BIT = 0,
    OVM_BIT = 1,
    TC_BIT = 2,
    C_BIT = 3,
    Z_BIT = 4,
    N_BIT = 5,
    V_BIT = 6,
    PM_BIT = 7,
    OVC_BIT = 10,
};

enum ST0_MASK {
    SXM_MASK = 1 << SXM_BIT,
    OVM_MASK = 1 << OVM_BIT,
    TC_MASK = 1 << TC_BIT,
    C_MASK = 1 << C_BIT,
    Z_MASK = 1 << Z_BIT,
    N_MASK = 1 << N_BIT,
    V_MASK = 1 << V_BIT,
    PM_MASK = 0b111 << PM_BIT,
    OVC_MASK = 0b111111 << OVC_BIT,
};

/* st1 */
enum ST1_BIT{
    INTM_BIT = 0,
    DBGM_BIT = 1,
    PAGE0_BIT = 2,
    VMAP_BIT = 3,
    SPA_BIT = 4,
    LOOP_BIT = 5,
    EALLOW_BIT = 6,
    IDLESTAT_BIT = 7,
    AMODE_BIT = 8,
    OBJMODE_BIT = 9,
    Reserved_BIT = 10,
    M0M1MAP_BIT = 11,
    XF_BIT = 12,
    ARP_BIT = 13,
};

enum ST1_MASK{
    INTM_MASK = 1 << 0,
    DBGM_MASK = 1 << 1,
    PAGE0_MASK = 1 << 2,
    VMAP_MASK = 1 << 3,
    SPA_MASK = 1 << 4,
    LOOP_MASK = 1 << 5,
    EALLOW_MASK = 1 << 6,
    IDLESTAT_MASK = 1 << 7,
    AMODE_MASK = 1 << 8,
    OBJMODE_MASK = 1 << 9,
    Reserved_MASK = 1 << 10,
    M0M1MAP_MASK = 1 << 11,
    XF_MASK = 1 << 12,
    ARP_MASK = 0b111 << 13,
};

/* stf */
enum STF_BIT{
    LVF_BIT = 0,
    LUF_BIT = 1,
    NF_BIT = 2,
    ZF_BIT = 3,
    NI_BIT = 4,
    ZI_BIT = 5,
    TF_BIT = 6,
    RND32_BIT = 9,
    SHDWS_BIT = 31,
};

enum STF_MASK{
    LVF_MASK = 1 << 0,
    LUF_MASK = 1 << 1,
    NF_MASK = 1 << 2,
    ZF_MASK = 1 << 3,
    NI_MASK = 1 << 4,
    ZI_MASK = 1 << 5,
    TF_MASK = 1 << 6,
    RND32_MASK = 1 << 9,
    SHDWS_MASK = 1 << 31,
};

/* RB */
enum RB_BIT {
    RC_BIT = 0,
    RE_BIT = 16,
    RSIZE_BIT = 23,
    RA_BIT = 30,
    RAS_BIT = 31,
};

enum RB_MASK {
    RC_MASK = 0xffff << 0, //16bit
    RE_MASK = 0b1111111 << 16, //7bit
    RSIZE_MASK = 0b1111111 << 23, //7bit
    RA_MASK = 1 << 30,
    RAS_MASK = 1 << 31,
};

/* loc */
enum {
    LOC16 = 1,
    LOC32 = 2,
};

typedef struct CPUTms320c28xState {
    target_ulong acc; /* Accumulator todo:ah,al */
    target_ulong xar[8]; /* Auxiliary todo:ar[8] */
    target_ulong dp; /* Data-page pointer */
    target_ulong ifr; /* Interrupt flag register */
    target_ulong ier; /* Interrupt enable register */
    target_ulong dbgier; /* Debug interrupt enable register */
    target_ulong p; /* Product register todo:ph,pl*/
    
    target_ulong pc;          /* Program counter 22bit, reset to 0x3f ffc0*/
    target_ulong rpc;         /* Return program counte 22bit*/

    target_ulong sp; /* Stack pointer, reset to 0x0400 */
    target_ulong st0; /* Status register 0 */
    target_ulong st1; /* Status register 1, reset to 0x080b */

    target_ulong xt;         /* Multiplicand register, todo:t,tl*/

    target_ulong rh[8]; /* Floating-point result register */
    target_ulong rh_shadow[8]; /* Floating-point result register shadow*/
    target_ulong stf; /* Floating-point status register */
    target_ulong stf_shadow; /* Floating-point status register shadow*/
    target_ulong rb; /* repeat block */
    float_status fp_status;

    target_ulong rptc;
    // target_ulong shadow[8];
    target_ulong tmp[8];
    /* Fields up to this point are cleared by a CPU reset */
    struct {} end_reset_fields;

    void *irq[NR_IRQS];          /* Interrupt irq input */
    target_ulong insn_length;//used in irq, to determine length of current instruction
} CPUTms320c28xState;

/**
 * Tms320c28xCPU:
 * @env: #CPUTms320c28xState
 *
 * A Tms320c28x CPU.
 */
typedef struct Tms320c28xCPU {
    /*< private >*/
    CPUState parent_obj;
    /*< public >*/

    CPUNegativeOffsetState neg;
    CPUTms320c28xState env;
} Tms320c28xCPU;


void cpu_tms320c28x_list(void);
void tms320c28x_cpu_do_interrupt(CPUState *cpu); //callback when handle interrupt, used in tms320c28x_cpu_class_init
bool tms320c28x_cpu_exec_interrupt(CPUState *cpu, int int_req); // callback in cpu_exec, used in tms320c28x_cpu_class_init
void tms320c28x_cpu_dump_state(CPUState *cpu, FILE *f, int flags); //used in tms320c28x_cpu_class_init
hwaddr tms320c28x_cpu_get_phys_page_debug(CPUState *cpu, vaddr addr);  //used in tms320c28x_cpu_class_init
int tms320c28x_cpu_gdb_read_register(CPUState *cpu, uint8_t *buf, int reg);//used in tms320c28x_cpu_class_init
int tms320c28x_cpu_gdb_write_register(CPUState *cpu, uint8_t *buf, int reg);//used in tms320c28x_cpu_class_init
void tms320c28x_translate_init(void); ////used in tms320c28x_cpu_class_init
bool tms320c28x_cpu_tlb_fill(CPUState *cs, vaddr address, int size,
                           MMUAccessType access_type, int mmu_idx,
                           bool probe, uintptr_t retaddr);  //used in tms320c28x_cpu_class_init
int cpu_tms320c28x_signal_handler(int host_signum, void *pinfo, void *puc); //todo ??
int print_insn_tms320c28x(bfd_vma addr, disassemble_info *info); //use in tms320c28x_disas_set_info



#define cpu_list cpu_tms320c28x_list
#define cpu_signal_handler cpu_tms320c28x_signal_handler

#ifndef CONFIG_USER_ONLY
extern const VMStateDescription vmstate_tms320c28x_cpu;

/* hw/tms320c28x_pic.c */
void cpu_tms320c28x_pic_init(Tms320c28xCPU *cpu);

// /* hw/tms320c28x_timer.c */
// void cpu_tms320c28x_clock_init(Tms320c28xCPU *cpu); //todo
// uint32_t cpu_tms320c28x_count_get(Tms320c28xCPU *cpu); //todo
// void cpu_tms320c28x_count_set(Tms320c28xCPU *cpu, uint32_t val); //todo
// void cpu_tms320c28x_count_update(Tms320c28xCPU *cpu); //todo
// void cpu_tms320c28x_timer_update(Tms320c28xCPU *cpu); //todo
// void cpu_tms320c28x_count_start(Tms320c28xCPU *cpu); //todo
// void cpu_tms320c28x_count_stop(Tms320c28xCPU *cpu); //todo
#endif

#define TMS320C28X_CPU_TYPE_SUFFIX "-" TYPE_TMS320C28X_CPU
#define TMS320C28X_CPU_TYPE_NAME(model) model TMS320C28X_CPU_TYPE_SUFFIX
#define CPU_RESOLVING_TYPE TYPE_TMS320C28X_CPU

typedef CPUTms320c28xState CPUArchState;
typedef Tms320c28xCPU ArchCPU;

#include "exec/cpu-all.h"

#define CPU_GET_STATUS(reg, bit) cpu_get_##reg(env, bit##_BIT, bit##_MASK)
#define CPU_SET_STATUS(reg, bit, value) cpu_set_##reg(env, value, bit##_BIT, bit##_MASK)

static inline uint32_t cpu_get_st0(const CPUTms320c28xState *env, uint32_t bit, uint32_t mask)
{
    return (env->st0 & mask) >> bit;
}

static inline uint32_t cpu_get_st1(const CPUTms320c28xState *env, uint32_t bit, uint32_t mask)
{
    return (env->st1 & mask) >> bit;
}

static inline uint32_t cpu_get_stf(const CPUTms320c28xState *env, uint32_t bit, uint32_t mask)
{
    return (env->stf & mask) >> bit;
}

static inline uint32_t cpu_get_rb(const CPUTms320c28xState *env, uint32_t bit, uint32_t mask)
{
    return (env->rb & mask) >> bit;
}

static inline void cpu_set_st0(CPUTms320c28xState *env,uint32_t value, uint32_t bit, uint32_t mask)
{
    value = (value << bit) & mask;
    env->st0 = (env->st0 & ~mask) | value;
}

static inline void cpu_set_st1(CPUTms320c28xState *env,uint32_t value, uint32_t bit, uint32_t mask)
{
    value = (value << bit) & mask;
    env->st1 = (env->st1 & ~mask) | value;
}

static inline void cpu_set_stf(CPUTms320c28xState *env,uint32_t value, uint32_t bit, uint32_t mask)
{
    value = (value << bit) & mask;
    env->stf = (env->stf & ~mask) | value;
}

static inline void cpu_set_rb(CPUTms320c28xState *env,uint32_t value, uint32_t bit, uint32_t mask)
{
    value = (value << bit) & mask;
    env->rb = (env->rb & ~mask) | value;
}

// 获取cpu状态，在查找tb时会比较pc,cs_base,flags
static inline void cpu_get_tb_cpu_state(CPUTms320c28xState *env,
                                        target_ulong *pc,
                                        target_ulong *cs_base, uint32_t *flags)
{
    *pc = env->pc * 2;
    // todo: cs_base,flags
    *cs_base = 0;
    *flags = 0;
}

//获取mmu状态，返回无mmu
static inline int cpu_mmu_index(CPUTms320c28xState *env, bool ifetch)
{
    int ret = MMU_NOMMU_IDX;  /* mmu is disabled */
    return ret;
}

#define CPU_INTERRUPT_INT   CPU_INTERRUPT_TGT_EXT_0

#define EXCP_INTERRUPT_RESET 0
#define EXCP_INTERRUPT_INT1 1
#define EXCP_INTERRUPT_INT2 2
#define EXCP_INTERRUPT_INT3 3
#define EXCP_INTERRUPT_INT4 4
#define EXCP_INTERRUPT_INT5 5
#define EXCP_INTERRUPT_INT6 6
#define EXCP_INTERRUPT_INT7 7
#define EXCP_INTERRUPT_INT8 8
#define EXCP_INTERRUPT_INT9 9
#define EXCP_INTERRUPT_INT10 10
#define EXCP_INTERRUPT_INT11 11
#define EXCP_INTERRUPT_INT12 12
#define EXCP_INTERRUPT_INT13 13
#define EXCP_INTERRUPT_INT14 14
#define EXCP_INTERRUPT_DLOGINT 15
#define EXCP_INTERRUPT_RTOSINT 16
#define EXCP_INTERRUPT_RESERVED 17
#define EXCP_INTERRUPT_NMI 18
#define EXCP_INTERRUPT_ILLEGAL 19
#define EXCP_INTERRUPT_USER1 20
#define EXCP_INTERRUPT_USER2 21
#define EXCP_INTERRUPT_USER3 22
#define EXCP_INTERRUPT_USER4 23
#define EXCP_INTERRUPT_USER5 24
#define EXCP_INTERRUPT_USER6 25
#define EXCP_INTERRUPT_USER7 26
#define EXCP_INTERRUPT_USER8 27
#define EXCP_INTERRUPT_USER9 28
#define EXCP_INTERRUPT_USER10 29
#define EXCP_INTERRUPT_USER11 30
#define EXCP_INTERRUPT_USER12 31

static const char * const INTERRUPT_NAME[EXCP_INTERRUPT_USER12+1] = {
    "RESET","INT1","INT2","INT3","INT4","INT5","INT6","INT7","INT8","INT9","INT10",
    "INT11","INT12","INT13","INT14","DLOGINT","RTOSINT","Reserved","NMI",
    "ILLEGAL","USER1","USER2","USER3","USER4","USER5","USER6","USER7","USER8","USER9",
    "USER10","USER11","USER12"
};

#endif /* TMS320C28X_CPU_H */
