def parse_logfile(fname):
    import re
    FP = r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)'
    PERFLINE = re.compile(r'Q2 G''\[([\d ]{4})([\d ]{4})([\d ]{4})\] P\[ *(\d+) +(\d+) +(\d+)\]  '+FP+r' s +'+FP+r' GF +'+FP+r' MEq/s')
    HOSTLINE = re.compile(r'.*on a ([a-z\-_0-9]+) named \w+ with (\d+) processors')
    Sizes = []
    GFlops = []
    MEqs = []
    Procs = None
    HostName = 'unknown'
    with open(fname) as f:
        while 'Starting performance sampling' not in next(f):
            pass
        while True:
            line = next(f)
            m = re.match(PERFLINE,line)
            if not m:
                break
            g0,g1,g2, p0,p1,p2, time, gflops, meqs = m.groups()
            g = float(g0)*float(g1)*float(g2)
            p = int(p0)*int(p1)*int(p2)
            Sizes.append(g)
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

    return Sizes, GFlops, MEqs, HostName, Procs

def plot(logfiles, outputfile):
    symbols = iter(['ro', 'bv', 'ks', 'g^'])
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots()
    plt.title('FE-FAS Performance')
    plt.xlabel('Global number of elements')
    ax2 = ax1.twinx()
    #ax1.set_autoscaley_on(False)
    ax1.set_ylabel('MEquations/second')
    all_sizes = []
    all_gflops = []
    all_meqs = []
    max_meqs = 0
    for f in logfiles:
        sizes, gflops, meqs, hostname, procs = parse_logfile(f)
        all_sizes += sizes
        all_gflops += gflops
        all_meqs += meqs
        ax1.semilogx(sizes, meqs, next(symbols), label='%s np=%d'%(hostname, procs))
    flops_per_meqn = all_gflops[-1] / all_meqs[-1]
    ax1.set_xlim(0,max(all_sizes))
    ax2.set_xlim(0,max(all_sizes))
    ax2.set_autoscaley_on(False)
    ax1.set_ylim(0,1.1*max(all_meqs))
    ax2.set_ylim(0,1.1*max(all_meqs)*flops_per_meqn)
    ax2.set_ylabel('GFlop/s')
    ax1.legend(loc='upper left')
    if outputfile:
        plt.savefig(outputfile)
    else:
        plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser('FE-FAS Performance Analyzer')
    parser.add_argument('-o', '--output', type=str, help='Output file')
    parser.add_argument('logfiles', nargs='+', type=str, help='List of files to process, usually including -log_summary')
    args = parser.parse_args()
    plot(args.logfiles, args.output)
