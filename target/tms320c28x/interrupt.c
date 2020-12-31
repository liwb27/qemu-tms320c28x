/*
 * tms320c28x interrupt.
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

#include "qemu/osdep.h"
#include "cpu.h"
#include "exec/exec-all.h"
#include "exec/gdbstub.h"
#include "qemu/host-utils.h"
#include "ldst.h"

// handle exception/interrupt, cs->exception_index
void tms320c28x_cpu_do_interrupt(CPUState *cs)
{
#ifndef CONFIG_USER_ONLY
    Tms320c28xCPU *cpu = TMS320C28X_CPU(cs);
    CPUTms320c28xState *env = &cpu->env;
    int exception = cs->exception_index;

    if (exception >= EXCP_INTERRUPT_RESET && exception <= EXCP_INTERRUPT_RTOSINT + 100) {
        if (exception <= EXCP_INTERRUPT_USER12)
        {
            qemu_log_mask(CPU_LOG_INT, "INTERRUPT: %s\n", INTERRUPT_NAME[exception]);
        }
        if (exception > EXCP_INTERRUPT_USER12) {//use +100 to intr and hw int
            exception = exception - 100;
            qemu_log_mask(CPU_LOG_INT, "INTERRUPT: %s\n", INTERRUPT_NAME[exception]);
            // clear corresponding IFR bit
            uint32_t mask = 1 << (exception - 1);
            env->ifr = env->ifr & (~mask);
            // clear corresponding IER bit
            env->ier = env->ier & (~mask);
        }

        // when interrupt is taken, save RA to RAS bit
        int ra = cpu_get_rb(env, RA_BIT, RA_MASK);
        cpu_set_rb(env, ra, RAS_BIT, RAS_MASK);
        
        // temp = pc + 1 or + 2
        int temp = env->pc + env->insn_length; 
        // sp = sp + 1, for odd address
        if ((env->sp & 1) == 1) {
            env->sp += 1;
        }
        // [SP] = T:ST0
        st32_swap(env, env->sp, (env->xt & 0xffff0000) | (env->st0 & 0xffff));
        env->sp += 2;
        // [SP] = AH:AL
        st32_swap(env, env->sp, env->acc);
        env->sp += 2;
        // [SP] = PH:PL
        st32_swap(env, env->sp, env->p);
        env->sp += 2;
        // [SP] = AR1:AR0
        uint32_t ar1 = env->xar[1] & 0xffff;
        uint32_t ar0 = env->xar[0] & 0xffff;
        st32_swap(env, env->sp, (ar1<<16) | ar0);
        env->sp += 2;
        // [SP] = DP:ST1
        st32_swap(env, env->sp, (env->dp<<16) | (env->st1 & 0xffff));
        env->sp += 2;
        // [SP] = DBGSTAT:IER
        st32_swap(env, env->sp, (env->dbgier<<16) | (env->ier & 0xffff));
        env->sp += 2;
        // [SP] = temp
        st32_swap(env, env->sp, temp);
        env->sp += 2;
        // INTM = 0
        cpu_set_st1(env, 0, INTM_BIT, INTM_MASK);
        // DBGM = 1
        cpu_set_st1(env, 0, DBGM_BIT, DBGM_MASK);
        // EALLOW = 0
        cpu_set_st1(env, 0, EALLOW_BIT, EALLOW_MASK);
        // LOOP = 0
        cpu_set_st1(env, 0, LOOP_BIT, LOOP_MASK);
        // IDLESTAT = 0
        cpu_set_st1(env, 0, IDLESTAT_BIT, IDLESTAT_MASK);
        // PC = fetched vector
        uint32_t vector_base = 0x3fffc0;
        uint32_t addr = vector_base + exception * 2;
        env->pc = ld32_swap(env, addr);
        
    } else {
        cpu_abort(cs, "Unhandled exception 0x%x\n", exception);
    }
#endif

    cs->exception_index = -1;
}


// handle int from hw
bool tms320c28x_cpu_exec_interrupt(CPUState *cs, int interrupt_request)
{
    Tms320c28xCPU *cpu = TMS320C28X_CPU(cs);
    CPUTms320c28xState *env = &cpu->env;
    int exception = cs->exception_index - 100;
    if (interrupt_request == CPU_INTERRUPT_INT) {
        // set the corresponding IFR bit
        env->ifr = env->ifr | (1 << (exception - 1));
        // is interrupt enabled?
        int ier_bit = (env->ier >> (exception - 1)) & 1;
        int intm_bit = cpu_get_st1(env, INTM_BIT, INTM_MASK);
        if (ier_bit == 1 && intm_bit == 0)//interrupt enabled
        {
            tms320c28x_cpu_do_interrupt(cs);
            return true;
        }
    }
    return false;
}
