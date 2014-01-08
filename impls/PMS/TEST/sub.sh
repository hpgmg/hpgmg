#PBS -q regular
#PBS -l mppwidth=49152
##PBS -l mppwidth=768
#PBS -l walltime=0:30:00
#PBS -N hpgmg
#PBS -V 
#PBS -j eo
#PBS -A m1792

cd $PBS_O_WORKDIR

# -S 4 for Hopper, -S 8 for Edison

#
# 1 solve, N=256
#
aprun -n 32768 -S 4 ../pms.ex -verbose 1 -nxloc 256 -num_solves 1 >& out_Fcycle_032768_nx256_1solves
/bin/mv Convergence.history Convergence.Fcycle.03276cores8.256N.1solve.history
/bin/mv Run.history Run.Fcycle.03276cores8.256N.1solve.history

aprun -n 4096 -S 4 ../pms.ex -verbose 1 -nxloc 256 -num_solves 1 >& out_Fcycle_004096_nx256_1solves
/bin/mv Convergence.history Convergence.Fcycle.004096cores.256N.1solve.history
/bin/mv Run.history Run.Fcycle.004096cores.256N.1solve.history

aprun -n 512 -S 4 ../pms.ex -verbose 1 -nxloc 256 -num_solves 1 >& out_Fcycle_000512_nx256_1solves
/bin/mv Convergence.history Convergence.Fcycle.000512cores.256N.1solve.history
/bin/mv Run.history Run.Fcycle.000512cores.256N.1solve.history

aprun -n 64 -S 4 ../pms.ex -verbose 1 -nxloc 256 -num_solves 1 >& out_Fcycle_000064_nx256_1solves
/bin/mv Convergence.history Convergence.Fcycle.000064cores.256N.1solve.history
/bin/mv Run.history Run.Fcycle.000064cores.256N.1solve.history

aprun -n 8 -S 4 ../pms.ex -verbose 1 -nxloc 256 -num_solves 1 >& out_Fcycle_000008_nx256_1solves
/bin/mv Convergence.history Convergence.Fcycle.000008cores.256N.1solve.history
/bin/mv Run.history Run.Fcycle.000008cores.256N.1solve.history
#
# 512 solves, N=32
#
aprun -n 32768 -S 4 ../pms.ex -verbose 1 -nxloc 32 -num_solves 512 >& out_Fcycle_032768_nx256_512solves
/bin/mv Convergence.history Convergence.Fcycle.032768cores.032N.512solve.history
/bin/mv Run.history Run.Fcycle.032768cores.032N.512solve.history

aprun -n 4096 -S 4 ../pms.ex -verbose 1 -nxloc 32 -num_solves 512 >& out_Fcycle_004096_nx032_512solves
/bin/mv Convergence.history Convergence.Fcycle.004096cores.032N.512solve.history
/bin/mv Run.history Run.Fcycle.004096cores.032N.512solve.history

aprun -n 512 -S 4 ../pms.ex -verbose 1 -nxloc 32 -num_solves 512 >& out_Fcycle_000512_nx032_512solves
/bin/mv Convergence.history Convergence.Fcycle.000512cores.032N.512solve.history
/bin/mv Run.history Run.Fcycle.000512cores.032N.512solve.history

aprun -n 64 -S 4 ../pms.ex -verbose 1 -nxloc 32 -num_solves 512 >& out_Fcycle_000064_nx032_512solves
/bin/mv Convergence.history Convergence.Fcycle.000064cores.032N.512solve.history
/bin/mv Run.history Run.Fcycle.000064cores.032N.512solve.history

aprun -n 8 -S 4 ../pms.ex -verbose 1 -nxloc 32 -num_solves 512 >& out_Fcycle_000008_nx032_512solves
/bin/mv Convergence.history Convergence.Fcycle.000008cores.032N.512solve.history
/bin/mv Run.history Run.Fcycle.000008cores.032N.512solve.history
#
# V-cycles, 1 solve, N=256
#
aprun -n 32768 -S 4 ../pms.ex -verbose 1 -nxloc 256 -nvcycles 100 -nfcycles 0 -rtol 1.e-6 -num_solves 1 >& out_Vcycles_032768_nx256_1solves
/bin/mv Convergence.history Convergence.Vcycles.03276cores8.256N.1solve.history
/bin/mv Run.history Run.Vcycles.03276cores8.256N.1solve.history

aprun -n 4096 -S 4 ../pms.ex -verbose 1 -nxloc 256 -nvcycles 100 -nfcycles 0 -rtol 1.e-6 -num_solves 1 >& out_Vcycles_004096_nx256_1solves
/bin/mv Convergence.history Convergence.Vcycles.004096cores.256N.1solve.history
/bin/mv Run.history Run.Vcycles.004096cores.256N.1solve.history

aprun -n 512 -S 4 ../pms.ex -verbose 1 -nxloc 256 -nvcycles 100 -nfcycles 0 -rtol 1.e-6 -num_solves 1 >& out_Vcycles_000512_nx256_1solves
/bin/mv Convergence.history Convergence.Vcycles.000512cores.256N.1solve.history
/bin/mv Run.history Run.Vcycles.000512cores.256N.1solve.history

aprun -n 64 -S 4 ../pms.ex -verbose 1 -nxloc 256 -nvcycles 100 -nfcycles 0 -rtol 1.e-6 -num_solves 1 >& out_Vcycles_000064_nx256_1solves
/bin/mv Convergence.history Convergence.Vcycles.000064cores.256N.1solve.history
/bin/mv Run.history Run.Vcycles.000064cores.256N.1solve.history

aprun -n 8 -S 4 ../pms.ex -verbose 1 -nxloc 256 -nvcycles 100 -nfcycles 0 -rtol 1.e-6 -num_solves 1 >& out_Vcycles_000008_nx256_1solves
/bin/mv Convergence.history Convergence.Vcycles.000008cores.256N.1solve.history
/bin/mv Run.history Run.Vcycles.000008cores.256N.1solve.history
