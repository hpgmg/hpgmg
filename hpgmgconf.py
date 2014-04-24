import sys,os

def mkdir_p(path):
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
    parser.add_argument('--with-hpm', help='Use libHPM profiling library on Blue Gene')
    cf = parser.add_argument_group('Compilers and flags')
    cf.add_argument('--CC', help='Path to C compiler', default=os.environ.get('CC',''))
    cf.add_argument('--CFLAGS', help='Flags for C compiler', default=os.environ.get('CFLAGS',''))
    cf.add_argument('--LDFLAGS', help='Flags to pass to linker', default=os.environ.get('LDFLAGS',''))
    cf.add_argument('--LDLIBS', help='Libraries to pass to linker', default=os.environ.get('LDLIBS',''))
    fv = parser.add_argument_group('Finite Volume options')
    fv.add_argument('--no-fv', action='store_false', dest='fv', help='Do not build the Finite-Volume solver')
    fv.add_argument('--no-fv-mpi', action='store_false', dest='fv_mpi', help='Use MPI')
    fv.add_argument('--fv-cycle', help='Multigrid cycle type', choices=['V','F','W'], default='F')
    fv.add_argument('--no-fv-subcomm', action='store_false', dest='fv_subcomm', help='Build a subcommunicator for each level in the MG v-cycle to minimize the scope of MPI_AllReduce()')
    fv.add_argument('--fv-coarse-solver', help='Use BiCGStab as a bottom (coarse grid) solver', choices=['bicgstab','cabicgstab','cg','cacg'], default='bicgstab')
    fv.add_argument('--fv-smoother', help='Multigrid smoother', choices=['cheby','gsrb','jacobi','l1jacobi'], default='cheby')
    fv.add_argument('--fv-timer', help='Timer implementation', choices=['omp','mpi'], default='mpi')
    args = parser.parse_args()
    if args.arch is None:
        args.arch = args.petsc_arch
    if not args.arch:
        args.arch = 'build'
    mkdir_p(args.arch)
    configure(args)

def configure(args):
    open(os.path.join(args.arch,'Makefile'), 'w').write(makefile(args))
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

def makefile(args):
    if args.petsc_dir and not args.CC:
        CC = '$(PCC)'
    else:
        CC = args.CC
    m = ['HPGMG_ARCH = %s' % args.arch,
        'HPGMG_CC = %s' % CC,
        'HPGMG_CFLAGS = %s' % args.CFLAGS,
        'HPGMG_LDFLAGS = %s' % args.LDFLAGS,
        'HPGMG_LDLIBS = %s' % args.LDLIBS,
        'PETSC_DIR = %s' % args.petsc_dir,
        'PETSC_ARCH = %s' % args.petsc_arch,
        'PYTHON = %s' % sys.executable,
        'SRCDIR = %s' % os.path.abspath(os.path.dirname(__name__)),]
    if args.petsc_dir:
        m.append('HPGMG_CFLAGS += $(PCC_FLAGS) $(CCPPFLAGS)')
    if args.with_hpm:
        m.append('CONFIG_HPM = y')
        hpm_lib = args.with_hpm
        if not isinstance(hpm_lib,str): # ALCF location
            hpm_lib = '/soft/perftools/hpctw/lib/libmpihpm.a /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a'
        for p in hpm_lib.split():
            assert os.path.exists(p), "HPM path '%s' not found" % p
        m.append('HPGMG_LDLIBS += ' + hpm_lib)
        m.append('HPGMG_CPPFLAGS += -DCONFIG_HPM')
    if args.fv:
        m.append('CONFIG_FV = y')
    if args.petsc_dir:
        m.append('CONFIG_FE = y')
    m.append('CONFIG_TIMER_%s = y' % args.fv_timer.upper())
    m.append('CONFIG_FV_CFLAGS = ' + hpgmg_fv_cflags(args))
    if args.petsc_dir:
        m.append('include $(PETSC_DIR)/conf/variables')
    m.append('include $(SRCDIR)/base.mk\n')
    return '\n'.join(m)

def hpgmg_fv_cflags(args):
    defines = []
    if args.fv_mpi:
        defines.append('USE_MPI')
    defines.append('USE_%s' % args.fv_coarse_solver.upper())
    if args.fv_subcomm:
        defines.append('USE_SUBCOMM')
    defines.append('USE_%sCYCLES' % args.fv_cycle.upper())
    defines.append('USE_%s' % args.fv_smoother.upper())
    defines.append('STENCIL_FUSE_DINV')
    defines.append('STENCIL_FUSE_BC')
    return ' '.join('-D%s=1'%d for d in defines)
