def parse_logfile(fname):
    import re
    FP = r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)'
    PERFLINE = []
    PERFLINE.append(re.compile(r'Q2 G''\[([\d ]{4})([\d ]{4})([\d ]{4})\] P\[ *(\d+) +(\d+) +(\d+)\]  '+FP+r' s +'+FP+r' GF +'+FP+r' MEq/s'))
    PERFLINE.append(re.compile(r'Q2 G''\[([\d ]{5})([\d ]{5})([\d ]{5})\] P\[ *(\d+) +(\d+) +(\d+)\]  '+FP+r' s +'+FP+r' GF +'+FP+r' MEq/s'))
    HOSTLINE = re.compile(r'.*on a ([a-z\-_0-9]+) named [^ ]+ with (\d+) processors')
    Dofs = []
    GFlops = []
    MEqs = []
    Procs = None
    HostName = 'unknown'
    with open(fname) as f:
        while 'Starting performance sampling' not in next(f):
            pass
        while True:
            line = next(f)
            for perfline in PERFLINE:
                m = re.match(perfline,line)
                if m: break
            if not m: break
            g0,g1,g2, p0,p1,p2, time, gflops, meqs = m.groups()
            g = (float(g0)*2+1)*(float(g1)*2+1)*(float(g2)*2+1)
            p = int(p0)*int(p1)*int(p2)
            Dofs.append(g)
            GFlops.append(float(gflops))
            MEqs.append(float(meqs))
            if Procs is None:
                Procs = p
            elif p != Procs:
                raise RuntimeError('Procs varies within file "%s"' % (fname,))
        while True:
            line = next(f)
            m = re.match(HOSTLINE,line)
            if m:
                HostName, p = m.groups()
                assert int(p) == Procs
                break

    return Dofs, GFlops, MEqs, HostName, Procs

def plot(args):
    symbols = iter(['ro', 'bv', 'ks', 'g^', 'bx'])
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    plt.title('HPGMG-FE Performance')
    if args.perprocess:
        plt.xlabel('Number of equations/process')
    else:
        plt.xlabel('Global number of equations')
    ax2 = ax1.twinx()
    #ax1.set_autoscaley_on(False)
    ax1.set_ylabel('MEquations/second')
    all_dofs = []
    all_gflops = []
    all_meqs = []
    max_meqs = 0
    for f in args.logfiles:
        dofs, gflops, meqs, hostname, procs = parse_logfile(f)
        if args.perprocess:
            dofs = [d/procs for d in dofs]
        all_dofs += dofs
        all_gflops += gflops
        all_meqs += meqs
        if args.loglog:
            ax1.loglog(dofs, meqs, next(symbols), label='%s np=%d'%(hostname, procs))
        else:
            ax1.semilogx(dofs, meqs, next(symbols), label='%s np=%d'%(hostname, procs))
    flops_per_meqn = all_gflops[-1] / all_meqs[-1]
    ax1.set_xlim(0.9*min(all_dofs),1.05*max(all_dofs))
    ax2.set_xlim(0.9*min(all_dofs),1.05*max(all_dofs))
    ax2.set_autoscaley_on(False)
    if args.loglog:
        ax2.set_yscale('log')
        ax1.legend(loc='lower right')
    else:
        ax1.legend(loc='upper left')
    ax1.set_ylim(0.9*min(all_meqs),1.1*max(all_meqs))
    ax2.set_ylim(0.9*min(all_meqs)*flops_per_meqn,1.1*max(all_meqs)*flops_per_meqn)
    ax2.set_ylabel('GFlop/s')
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser('FE-FAS Performance Analyzer')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('--loglog', action='store_true', help='Use logarithmic y axis (x is always logarithmic)')
    parser.add_argument('--perprocess', action='store_true', help='Use problem size per process for x axis')
    parser.add_argument('logfiles', nargs='+', type=str, help='List of files to process, usually including -log_summary')
    args = parser.parse_args()
    plot(args)
