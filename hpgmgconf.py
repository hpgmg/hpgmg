def mkdir_p(path):
    import os
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == os.errno.EEXIST:
            pass
        else: raise

def main():
    import argparse, os
    parser = argparse.ArgumentParser(description='Configure High-performance Geometric Multigrid (HPGMG)')
    parser.add_argument('--arch', help='Name of this configuration', default=None)
    parser.add_argument('--petsc-dir', help='PETSC_DIR', default=os.environ.get('PETSC_DIR',''))
    parser.add_argument('--petsc-arch', help='PETSC_ARCH', default=os.environ.get('PETSC_ARCH',''))
    cf = parser.add_argument_group('Compilers and flags')
    cf.add_argument('--CC', help='Path to C compiler', default='mpicc')
    cf.add_argument('--CFLAGS', help='Flags for C compiler', default='-Wall')
    cf.add_argument('--LDFLAGS', help='Flags to pass to linker', default='')
    cf.add_argument('--LDLIBS', help='Libraries to pass to linker', default='')
    fv = parser.add_argument_group('Finite Volume options')
    fv.add_argument('--no-fv', action='store_false', dest='fv', help='Do not build the Finite-Volume solver')
    fv.add_argument('--fv-mpi', action='store_true', help='Use MPI', default=True)
    fv.add_argument('--fv-cycle', help='Use U-, V-, or F-cycles', default='F', type=str)
    fv.add_argument('--fv-bicgstab', action='store_true', help='Use BiCGStab as a bottom (coarse grid) solver', default=False)
    fv.add_argument('--fv-subcomm', action='store_true', help='build a subcommunicator for each level in the MG v-cycle to minimize the scope of MPI_AllReduce()', default=True)
    fv.add_argument('--fv-smoother', help='Use Cheby, GSRB, Jacobi, or L1Jacobi for smoother', default='Cheby', type=str)
    fv.add_argument('--fv-timer', help='Use "x86" or "bgq" timers', default='x86', type=str);
    args = parser.parse_args()
    if args.arch is None:
        args.arch = args.petsc_arch
    if not args.arch:
        args.arch = 'build'
    mkdir_p(args.arch)
    configure(args)

def configure(args):
    import sys,os
    mkdir_p(os.path.join(args.arch, 'include'))
    with open(os.path.join(args.arch, 'include', 'hpgmgconf.h'), 'w') as f:
        f.write(hpgmgconf_h(args))
        open(os.path.join(args.arch,'Makefile'), 'w').write('\n'.join([
            'HPGMG_ARCH = %s' % args.arch,
            'HPGMG_CC = %s' % args.CC,
            'HPGMG_CFLAGS = %s' % args.CFLAGS,
            'HPGMG_LDFLAGS = %s' % args.LDFLAGS,
            'HPGMG_LDLIBS = %s' % args.LDLIBS,
            'PETSC_DIR = %s' % args.petsc_dir,
            'PETSC_ARCH = %s' % args.petsc_arch,
            'PYTHON = %s' % sys.executable,
            'SRCDIR = %s' % os.path.abspath(os.path.dirname(__name__)),
            'CONFIG_FV = %s' % ('y' if args.fv else ''),
            'CONFIG_X86 = %s' % ('y' if args.fv_timer == 'x86' else ''),
            'CONFIG_BGQ = %s' % ('y' if args.fv_timer == 'bgq' else ''),
            'include $(PETSC_DIR)/conf/variables',
            'include $(SRCDIR)/base.mk',
        ]))
    reconfname = os.path.join(args.arch,'reconfigure-%s.py' % args.arch)
    open(reconfname, 'w').write('\n'.join([
                '#!'+sys.executable,
                'import os,sys',
                'from argparse import Namespace',
                "sys.path.insert(0, os.path.abspath('.'))",
                'import hpgmgconf',
                'hpgmgconf.configure(%r)' % args,
                ]))
    os.chmod(reconfname,0o755)
    print('Configuration complete in: %s' % os.path.realpath(args.arch))
    print('To build: make -j3 -C %s' % args.arch)

def hpgmgconf_h(args):
    defines = []
    if args.fv_mpi:
        defines.append('MPI')
    if args.fv_bicgstab:
        defines.append('USE_BICGSTAB')
    if args.fv_subcomm:
        defines.append('USE_SUBCOMM')
    defines.append('USE_%sCYCLES' % args.fv_cycle.upper())
    defines.append('USE_%s' % args.fv_smoother.upper())
    return """#ifndef _hpgmgconf_h
#define _hpgmgconf_h
%s
#endif
""" % '\n'.join('#define __%s 1'%d for d in defines)
