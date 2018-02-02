#!/usr/bin/env python

import re

CREATE_RE = re.compile(r'PetscErrorCode OpCreate_([^_]\w*) ?\(')

def mangle(name):
    return name.lower().replace('_','-'), name

def build_ops(files):
    ops = []
    for src in files:
        with open(src) as f:
            for line in f:
                m = CREATE_RE.match(line)
                if m:
                    ops.append(mangle(m.groups()[0]))
    return ops

def genregister(outname, files):
    ops = build_ops(files)
    with open(outname, 'w') as out:
        out.write("""#include <op/fefas-op.h>

%(opdecl)s

PetscErrorCode OpRegisterAll_Generated()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  %(opreg)s
  PetscFunctionReturn(0);
}
""" % dict(opdecl='\n'.join('PetscErrorCode OpCreate_%s(Op);'%(o[1],) for o in ops),
           opreg='\n  '.join(['ierr = OpRegister("%s",OpCreate_%s);CHKERRQ(ierr);'%o for o in ops])))

if __name__ == '__main__':
    import sys
    genregister(sys.argv[1], sys.argv[2:])
