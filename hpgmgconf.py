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
    parser.add_argument('--CC', help='Path to C compiler', default='mpicc')
    parser.add_argument('--CFLAGS', help='Flags for C compiler', default='-Wall -Werror')
    parser.add_argument('--LDFLAGS', help='Flags to pass to linker', default='')
    parser.add_argument('--LDLIBS', help='Libraries to pass to linker', default='')
    parser.add_argument('--petsc-dir', help='PETSC_DIR', default=os.environ.get('PETSC_DIR',''))
    parser.add_argument('--petsc-arch', help='PETSC_ARCH', default=os.environ.get('PETSC_ARCH',''))
    args = parser.parse_args()
    if args.arch is None:
        args.arch = args.petsc_arch
    if not args.arch:
        args.arch = 'build'
    mkdir_p(args.arch)
    configure(args)

def configure(args):
    import sys,os
    open(os.path.join(args.arch,'Makefile'), 'w').write('\n'.join([
                'HPGMG_ARCH = %s' % args.arch,
                'PETSC_DIR = %s' % args.petsc_dir,
                'PETSC_ARCH = %s' % args.petsc_arch,
                'PYTHON = %s' % sys.executable,
                'SRCDIR = %s' % os.path.abspath(os.path.dirname(__name__)),
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
