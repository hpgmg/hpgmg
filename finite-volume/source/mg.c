//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "timers.h"
#include "defines.h"
#include "level.h"
#include "operators.h"
#include "solvers.h"
#include "mg.h"
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int sendRank;
  int sendBoxID;
  int sendBox;
  int recvRank;
  int recvBoxID;
  int recvBox;
  int i,j,k; // offsets used to index into the coarse box
} RP_type;

int qsortRP(const void *a, const void*b){
  RP_type *rpa = (RP_type*)a;
  RP_type *rpb = (RP_type*)b;
  // sort first by sendRank
  if(rpa->sendRank < rpb->sendRank)return(-1);
  if(rpa->sendRank > rpb->sendRank)return( 1);
  // then by sendBoxID
  if(rpa->sendBoxID < rpb->sendBoxID)return(-1);
  if(rpa->sendBoxID > rpb->sendBoxID)return( 1);
  return(0);
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
void MGPrintTiming(mg_type *all_grids){
  int level,num_levels = all_grids->num_levels;
  uint64_t _timeStart=CycleTime();sleep(1);uint64_t _timeEnd=CycleTime();
  double SecondsPerCycle = (double)1.0/(double)(_timeEnd-_timeStart);
  double scale = SecondsPerCycle/(double)all_grids->MGSolves_performed; // prints average performance per MGSolve

   
  #if 0 // find the maximum time taken by any process in any function... 
  #warning Will print the max time spent by any process for each operation/level
  for(level=0;level<num_levels;level++)max_level_timers(all_grids->levels[level]);
  #endif

  if(all_grids->my_rank!=0)return;
  double time,total;
          printf("\n\n");
          printf("                          ");for(level=0;level<(num_levels  );level++){printf("%12d ",level);}printf("\n");
        //printf("v-cycles initiated        ");for(level=0;level<(num_levels  );level++){printf("%12d ",all_grids->levels[level]->vcycles_from_this_level/all_grids->MGSolves_performed);}printf("\n");
          printf("box dimension             ");for(level=0;level<(num_levels  );level++){printf("%10d^3 ",all_grids->levels[level]->box_dim);}printf("       total\n");
  total=0;printf("------------------        ");for(level=0;level<(num_levels+1);level++){printf("------------ ");}printf("\n");
  total=0;printf("smooth                    ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.smooth;               total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("residual                  ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.residual;             total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("applyOp                   ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.apply_op;             total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("BLAS1                     ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.blas1;                total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("BLAS3                     ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.blas3;                total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("Boundary Conditions       ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.boundary_conditions;  total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("Restriction               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_total;    total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  local restriction       ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_local;    total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #ifdef USE_MPI
  total=0;printf("  pack MPI buffers        ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_pack;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers      ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_unpack;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_send;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_recv;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall             ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.restriction_wait;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  total=0;printf("Interpolation             ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_total;  total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  local interpolation     ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_local;  total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #ifdef USE_MPI
  total=0;printf("  pack MPI buffers        ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_pack;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers      ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_unpack; total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_send;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_recv;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall             ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.interpolation_wait;   total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  total=0;printf("Ghost Zone Exchange       ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_total;      total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  local exchange          ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_local;      total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #ifdef USE_MPI
  total=0;printf("  pack MPI buffers        ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_pack;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers      ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_unpack;     total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_send;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv               ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_recv;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall             ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.ghostZone_wait;       total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  #ifdef USE_MPI
  total=0;printf("MPI_collectives           ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.collectives;          total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);
  #endif
  total=0;printf("------------------        ");for(level=0;level<(num_levels+1);level++){printf("------------ ");}printf("\n");
  total=0;printf("Total by level            ");for(level=0;level<(num_levels  );level++){time=scale*(double)all_grids->levels[level]->cycles.Total;                total+=time;printf("%12.6f ",time);}printf("%12.6f\n",total);

  printf("\n");
  printf( "   Total time in MGBuild  %12.6f seconds\n",SecondsPerCycle*(double)all_grids->cycles.MGBuild);
  printf( "   Total time in MGSolve  %12.6f seconds\n",scale*(double)all_grids->cycles.MGSolve);
  printf( "      number of v-cycles  %12d\n"  ,all_grids->levels[0]->vcycles_from_this_level/all_grids->MGSolves_performed);
  printf( "Bottom solver iterations  %12d\n"  ,all_grids->levels[num_levels-1]->Krylov_iterations/all_grids->MGSolves_performed);
  #if defined(USE_CABICGSTAB) || defined(USE_CACG)
  printf( "     formations of G[][]  %12d\n"  ,all_grids->levels[num_levels-1]->CAKrylov_formations_of_G/all_grids->MGSolves_performed);
  #endif
  printf("\n");
  double numDOF = (double)all_grids->levels[0]->dim.i*(double)all_grids->levels[0]->dim.j*(double)all_grids->levels[0]->dim.k;
  printf( "            Performance   %12.3e DOF/s\n",numDOF/(scale*(double)all_grids->cycles.MGSolve));
  printf("\n\n");fflush(stdout);
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
void MGResetTimers(mg_type *all_grids){
  int level;
  for(level=0;level<all_grids->num_levels;level++)reset_level_timers(all_grids->levels[level]);
//all_grids->cycles.MGBuild     = 0;
  all_grids->cycles.MGSolve     = 0;
  all_grids->MGSolves_performed = 0;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
void build_interpolation(mg_type *all_grids){
  int level;
  for(level=0;level<all_grids->num_levels;level++){
  all_grids->levels[level]->interpolation.num_recvs           = 0;
  all_grids->levels[level]->interpolation.num_sends           = 0;
  all_grids->levels[level]->interpolation.blocks[0]           = NULL;
  all_grids->levels[level]->interpolation.blocks[1]           = NULL;
  all_grids->levels[level]->interpolation.blocks[2]           = NULL;
  all_grids->levels[level]->interpolation.num_blocks[0]       = 0;
  all_grids->levels[level]->interpolation.num_blocks[1]       = 0;
  all_grids->levels[level]->interpolation.num_blocks[2]       = 0;
  all_grids->levels[level]->interpolation.allocated_blocks[0] = 0;
  all_grids->levels[level]->interpolation.allocated_blocks[1] = 0;
  all_grids->levels[level]->interpolation.allocated_blocks[2] = 0;


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct pack, send(to level-1), and local...
  if( (level>0) && (all_grids->levels[level]->num_my_boxes>0) ){ // not top  *and*  I have boxes to send
    // construct a list of fine boxes to be coarsened and sent to me...
    int numFineBoxes = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*
                       (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*
                       (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*
                                                               all_grids->levels[level]->num_my_boxes;
        int *fineRanks = (    int*)malloc(numFineBoxes*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *fineBoxes = (RP_type*)malloc(numFineBoxes*sizeof(RP_type)); 
        numFineBoxes       = 0;
    int numFineBoxesLocal  = 0;
    int numFineBoxesRemote = 0;
    int coarseBox;
    for(coarseBox=0;coarseBox<all_grids->levels[level]->num_my_boxes;coarseBox++){
      int bi,bj,bk;
      int   coarseBoxID = all_grids->levels[level]->my_boxes[coarseBox].global_box_id;
      int   coarseBox_i = all_grids->levels[level]->my_boxes[coarseBox].low.i / all_grids->levels[level]->box_dim;
      int   coarseBox_j = all_grids->levels[level]->my_boxes[coarseBox].low.j / all_grids->levels[level]->box_dim;
      int   coarseBox_k = all_grids->levels[level]->my_boxes[coarseBox].low.k / all_grids->levels[level]->box_dim;
      for(bk=0;bk<all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;bk++){
      for(bj=0;bj<all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;bj++){
      for(bi=0;bi<all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;bi++){
        int fineBox_i = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*coarseBox_i + bi;
        int fineBox_j = (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*coarseBox_j + bj;
        int fineBox_k = (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*coarseBox_k + bk;
        int fineBoxID =  fineBox_i + fineBox_j*all_grids->levels[level-1]->boxes_in.i + fineBox_k*all_grids->levels[level-1]->boxes_in.i*all_grids->levels[level-1]->boxes_in.j;
        int fineBox   = -1;int f;for(f=0;f<all_grids->levels[level-1]->num_my_boxes;f++)if( all_grids->levels[level-1]->my_boxes[f].global_box_id == fineBoxID )fineBox=f; // try and find the index of a fineBox global_box_id == fineBoxID
        fineBoxes[numFineBoxes].sendRank  = all_grids->levels[level  ]->rank_of_box[coarseBoxID];
        fineBoxes[numFineBoxes].sendBoxID = coarseBoxID;
        fineBoxes[numFineBoxes].sendBox   = coarseBox;
        fineBoxes[numFineBoxes].recvRank  = all_grids->levels[level-1]->rank_of_box[  fineBoxID];
        fineBoxes[numFineBoxes].recvBoxID = fineBoxID;
        fineBoxes[numFineBoxes].recvBox   = fineBox;
        fineBoxes[numFineBoxes].i         = bi*all_grids->levels[level-1]->box_dim/2;
        fineBoxes[numFineBoxes].j         = bj*all_grids->levels[level-1]->box_dim/2;
        fineBoxes[numFineBoxes].k         = bk*all_grids->levels[level-1]->box_dim/2;
                  numFineBoxes++;
        if(all_grids->levels[level-1]->rank_of_box[fineBoxID] != all_grids->levels[level]->my_rank){
          fineRanks[numFineBoxesRemote++] = all_grids->levels[level-1]->rank_of_box[fineBoxID];
        }else{numFineBoxesLocal++;}
      }}}
    } // my (coarse) boxes
    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(fineBoxes,numFineBoxes      ,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(fineRanks,numFineBoxesRemote,sizeof(    int),qsortInt);
    int numFineRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numFineBoxesRemote;neighbor++)if(fineRanks[neighbor] != _rank){_rank=fineRanks[neighbor];fineRanks[numFineRanks++]=fineRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->interpolation.num_sends     =                         numFineRanks;
    all_grids->levels[level]->interpolation.send_ranks    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->interpolation.send_sizes    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->interpolation.send_buffers  =        (double**)malloc(numFineRanks*sizeof(double*));
    if(numFineRanks>0){
    if(all_grids->levels[level]->interpolation.send_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.send_ranks\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.send_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.send_sizes\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.send_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.send_buffers\n",level);exit(0);}
    }

    int elementSize = all_grids->levels[level-1]->box_dim*all_grids->levels[level-1]->box_dim*all_grids->levels[level-1]->box_dim;
    double * all_send_buffers = (double*)malloc(numFineBoxesRemote*elementSize*sizeof(double));
          if(numFineBoxesRemote*elementSize>0)
          if(all_send_buffers==NULL){fprintf(stderr,"malloc failed - interpolation/all_send_buffers\n");exit(0);}
                      memset(all_send_buffers,0,numFineBoxesRemote*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve
    //printf("level=%d, rank=%2d, send_buffers=%6d\n",level,all_grids->my_rank,numFineBoxesRemote*elementSize*sizeof(double));

    // for each neighbor, construct the pack list and allocate the MPI send buffer... 
    for(neighbor=0;neighbor<numFineRanks;neighbor++){
      int fineBox;
      int offset = 0;
      all_grids->levels[level]->interpolation.send_buffers[neighbor] = all_send_buffers;
      for(fineBox=0;fineBox<numFineBoxes;fineBox++)if(fineBoxes[fineBox].recvRank==fineRanks[neighbor]){
        // pack the MPI send buffer...
        append_block_to_list(&(all_grids->levels[level]->interpolation.blocks[0]),&(all_grids->levels[level]->interpolation.allocated_blocks[0]),&(all_grids->levels[level]->interpolation.num_blocks[0]),
          /* dim.i         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.j         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.k         = */ all_grids->levels[level-1]->box_dim/2,
          /* read.box      = */ fineBoxes[fineBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ fineBoxes[fineBox].i,
          /* read.j        = */ fineBoxes[fineBox].j,
          /* read.k        = */ fineBoxes[fineBox].k,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].kStride,
          /* read.scale    = */ 1,
          /* write.box     = */ -1,
          /* write.ptr     = */ all_grids->levels[level]->interpolation.send_buffers[neighbor],
          /* write.i       = */ offset,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ all_grids->levels[level-1]->box_dim,
          /* write.kStride = */ all_grids->levels[level-1]->box_dim*all_grids->levels[level-1]->box_dim,
          /* write.scale   = */ 2,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->interpolation.send_ranks[neighbor] = fineRanks[neighbor];
      all_grids->levels[level]->interpolation.send_sizes[neighbor] = offset;
      all_send_buffers+=offset;
    } // neighbor
    {
      int fineBox;
      for(fineBox=0;fineBox<numFineBoxes;fineBox++)if(fineBoxes[fineBox].recvRank==all_grids->my_rank){
        // local interpolations...
        append_block_to_list(&(all_grids->levels[level]->interpolation.blocks[1]),&(all_grids->levels[level]->interpolation.allocated_blocks[1]),&(all_grids->levels[level]->interpolation.num_blocks[1]),
          /* dim.i         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.j         = */ all_grids->levels[level-1]->box_dim/2,
          /* dim.k         = */ all_grids->levels[level-1]->box_dim/2,
          /* read.box      = */ fineBoxes[fineBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ fineBoxes[fineBox].i,
          /* read.j        = */ fineBoxes[fineBox].j,
          /* read.k        = */ fineBoxes[fineBox].k,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].sendBox].kStride,
          /* read.scale    = */ 1,
          /* write.box     = */ fineBoxes[fineBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ 0,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ all_grids->levels[level-1]->my_boxes[fineBoxes[fineBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level-1]->my_boxes[fineBoxes[fineBox].recvBox].kStride,
          /* write.scale   = */ 2,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
      }
    } // local to local interpolation

    // free temporary storage...
    free(fineBoxes);
    free(fineRanks);
  } // pack/send/local


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct recv(from level+1) and unpack...
  if( (level<all_grids->num_levels-1) && (all_grids->levels[level]->num_my_boxes>0) ){ // not bottom  *and*  I have boxes to receive

    // construct the list of coarsened boxes and neighboring ranks that will be interpolated and sent to me...
    int numCoarseBoxes = all_grids->levels[level]->num_my_boxes; // I may receive a block for each of my boxes
        int *coarseRanks = (    int*)malloc(numCoarseBoxes*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *coarseBoxes = (RP_type*)malloc(numCoarseBoxes*sizeof(RP_type)); 
        numCoarseBoxes       = 0;
    int fineBox;
    for(fineBox=0;fineBox<all_grids->levels[level]->num_my_boxes;fineBox++){
      int   fineBoxID = all_grids->levels[level]->my_boxes[fineBox].global_box_id;
      int   fineBox_i = all_grids->levels[level]->my_boxes[fineBox].low.i / all_grids->levels[level]->box_dim;
      int   fineBox_j = all_grids->levels[level]->my_boxes[fineBox].low.j / all_grids->levels[level]->box_dim;
      int   fineBox_k = all_grids->levels[level]->my_boxes[fineBox].low.k / all_grids->levels[level]->box_dim;
      int coarseBox_i = fineBox_i*all_grids->levels[level+1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;
      int coarseBox_j = fineBox_j*all_grids->levels[level+1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;
      int coarseBox_k = fineBox_k*all_grids->levels[level+1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;
      int coarseBoxID =  coarseBox_i + coarseBox_j*all_grids->levels[level+1]->boxes_in.i + coarseBox_k*all_grids->levels[level+1]->boxes_in.i*all_grids->levels[level+1]->boxes_in.j;
      if(all_grids->levels[level]->my_rank != all_grids->levels[level+1]->rank_of_box[coarseBoxID]){
        coarseBoxes[numCoarseBoxes].sendRank  = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
        coarseBoxes[numCoarseBoxes].sendBoxID = coarseBoxID;
        coarseBoxes[numCoarseBoxes].sendBox   = -1; 
        coarseBoxes[numCoarseBoxes].recvRank  = all_grids->levels[level  ]->rank_of_box[  fineBoxID];
        coarseBoxes[numCoarseBoxes].recvBoxID = fineBoxID;
        coarseBoxes[numCoarseBoxes].recvBox   = fineBox;
        coarseRanks[numCoarseBoxes] = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
                    numCoarseBoxes++;
      }
    } // my (fine) boxes

    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(coarseBoxes,numCoarseBoxes,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(coarseRanks,numCoarseBoxes,sizeof(    int),qsortInt);
    int numCoarseRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numCoarseBoxes;neighbor++)if(coarseRanks[neighbor] != _rank){_rank=coarseRanks[neighbor];coarseRanks[numCoarseRanks++]=coarseRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->interpolation.num_recvs     =                         numCoarseRanks;
    all_grids->levels[level]->interpolation.recv_ranks    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->interpolation.recv_sizes    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->interpolation.recv_buffers  =        (double**)malloc(numCoarseRanks*sizeof(double*));
    if(numCoarseRanks>0){
    if(all_grids->levels[level]->interpolation.recv_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.recv_ranks\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.recv_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.recv_sizes\n",level);exit(0);}
    if(all_grids->levels[level]->interpolation.recv_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->interpolation.recv_buffers\n",level);exit(0);}
    }

    int elementSize = all_grids->levels[level]->box_dim*all_grids->levels[level]->box_dim*all_grids->levels[level]->box_dim;
    double * all_recv_buffers = (double*)malloc(numCoarseBoxes*elementSize*sizeof(double)); 
          if(numCoarseBoxes*elementSize>0)
          if(all_recv_buffers==NULL){fprintf(stderr,"malloc failed - interpolation/all_recv_buffers\n");exit(0);}
                      memset(all_recv_buffers,0,numCoarseBoxes*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve
    //printf("level=%d, rank=%2d, recv_buffers=%6d\n",level,all_grids->my_rank,numCoarseBoxes*elementSize*sizeof(double));

    // for each neighbor, construct the unpack list and allocate the MPI recv buffer... 
    for(neighbor=0;neighbor<numCoarseRanks;neighbor++){
      int coarseBox;
      int offset = 0;
      all_grids->levels[level]->interpolation.recv_buffers[neighbor] = all_recv_buffers;
      for(coarseBox=0;coarseBox<numCoarseBoxes;coarseBox++)if(coarseBoxes[coarseBox].sendRank==coarseRanks[neighbor]){
        // unpack MPI recv buffer...
        append_block_to_list(&(all_grids->levels[level]->interpolation.blocks[2]),&(all_grids->levels[level]->interpolation.allocated_blocks[2]),&(all_grids->levels[level]->interpolation.num_blocks[2]),
          /* dim.i         = */ all_grids->levels[level]->box_dim,
          /* dim.j         = */ all_grids->levels[level]->box_dim,
          /* dim.k         = */ all_grids->levels[level]->box_dim,
          /* read.box      = */ -1,
          /* read.ptr      = */ all_grids->levels[level]->interpolation.recv_buffers[neighbor],
          /* read.i        = */ offset,
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ all_grids->levels[level]->box_dim,
          /* read.kStride  = */ all_grids->levels[level]->box_dim*all_grids->levels[level]->box_dim,
          /* read.scale    = */ 1,
          /* write.box     = */ coarseBoxes[coarseBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ 0,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].recvBox].kStride,
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->interpolation.recv_ranks[neighbor] = coarseRanks[neighbor];
      all_grids->levels[level]->interpolation.recv_sizes[neighbor] = offset;
      all_recv_buffers+=offset;
    } // neighbor

    // free temporary storage...
    free(coarseBoxes);
    free(coarseRanks);
  } // recv/unpack


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //MPI_Barrier(MPI_COMM_WORLD);
  //if(all_grids->my_rank==0){printf("================================================================================\n");}
  //print_communicator(0x7,all_grids->my_rank,level,&all_grids->levels[level]->interpolation);
  //if((all_grids->my_rank==0)&&(level==all_grids->num_levels-1))print_communicator(2,all_grids->my_rank,level,&all_grids->levels[level]->interpolation);
  //MPI_Barrier(MPI_COMM_WORLD);
  } // all levels


  #ifdef USE_MPI
  for(level=0;level<all_grids->num_levels;level++){
    all_grids->levels[level]->interpolation.requests = NULL;
    all_grids->levels[level]->interpolation.status   = NULL;
    if(level<all_grids->num_levels-1){  // i.e. bottom never calls interpolation()
    // by convention, level_f allocates a combined array of requests for both level_f recvs and level_c sends...
    int nMessages = all_grids->levels[level+1]->interpolation.num_sends + all_grids->levels[level]->interpolation.num_recvs;
    all_grids->levels[level]->interpolation.requests = (MPI_Request*)malloc(nMessages*sizeof(MPI_Request));
    all_grids->levels[level]->interpolation.status   = (MPI_Status *)malloc(nMessages*sizeof(MPI_Status ));
    }
  }
  #endif
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
void build_restriction(mg_type *all_grids, int restrictionType){
  int level;
  for(level=0;level<all_grids->num_levels;level++){
  all_grids->levels[level]->restriction[restrictionType].num_recvs           = 0;
  all_grids->levels[level]->restriction[restrictionType].num_sends           = 0;
  all_grids->levels[level]->restriction[restrictionType].blocks[0]           = NULL;
  all_grids->levels[level]->restriction[restrictionType].blocks[1]           = NULL;
  all_grids->levels[level]->restriction[restrictionType].blocks[2]           = NULL;
  all_grids->levels[level]->restriction[restrictionType].allocated_blocks[0] = 0;
  all_grids->levels[level]->restriction[restrictionType].allocated_blocks[1] = 0;
  all_grids->levels[level]->restriction[restrictionType].allocated_blocks[2] = 0;
  all_grids->levels[level]->restriction[restrictionType].num_blocks[0]       = 0; // number of unpack/insert operations  = number of boxes on level+1 that I don't own and restrict to 
  all_grids->levels[level]->restriction[restrictionType].num_blocks[1]       = 0; // number of unpack/insert operations  = number of boxes on level+1 that I own and restrict to
  all_grids->levels[level]->restriction[restrictionType].num_blocks[2]       = 0; // number of unpack/insert operations  = number of boxes on level-1 that I don't own that restrict to me


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct pack, send, and local...
  if( (level<all_grids->num_levels-1) && (all_grids->levels[level]->num_my_boxes>0) ){ // not bottom  *and*  I have boxes to send

    // construct the list of coarsened boxes and neighboring ranks...
    int numCoarseBoxes = (all_grids->levels[level]->boxes_in.i/all_grids->levels[level+1]->boxes_in.i)*
                         (all_grids->levels[level]->boxes_in.j/all_grids->levels[level+1]->boxes_in.j)*
                         (all_grids->levels[level]->boxes_in.k/all_grids->levels[level+1]->boxes_in.k)*
                          all_grids->levels[level]->num_my_boxes;
        int *coarseRanks = (    int*)malloc(numCoarseBoxes*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *coarseBoxes = (RP_type*)malloc(numCoarseBoxes*sizeof(RP_type)); 
        numCoarseBoxes       = 0;
    int numCoarseBoxesLocal  = 0;
    int numCoarseBoxesRemote = 0;
    int fineBox;
    for(fineBox=0;fineBox<all_grids->levels[level]->num_my_boxes;fineBox++){
      int   fineBoxID = all_grids->levels[level]->my_boxes[fineBox].global_box_id;
      int   fineBox_i = all_grids->levels[level]->my_boxes[fineBox].low.i / all_grids->levels[level]->box_dim;
      int   fineBox_j = all_grids->levels[level]->my_boxes[fineBox].low.j / all_grids->levels[level]->box_dim;
      int   fineBox_k = all_grids->levels[level]->my_boxes[fineBox].low.k / all_grids->levels[level]->box_dim;
      int coarseBox_i = fineBox_i*all_grids->levels[level+1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;
      int coarseBox_j = fineBox_j*all_grids->levels[level+1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;
      int coarseBox_k = fineBox_k*all_grids->levels[level+1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;
      int coarseBoxID =  coarseBox_i + coarseBox_j*all_grids->levels[level+1]->boxes_in.i + coarseBox_k*all_grids->levels[level+1]->boxes_in.i*all_grids->levels[level+1]->boxes_in.j;
      int coarseBox   = -1;int c;for(c=0;c<all_grids->levels[level+1]->num_my_boxes;c++)if( all_grids->levels[level+1]->my_boxes[c].global_box_id == coarseBoxID )coarseBox=c; // try and find the coarseBox index of a box with global_box_id == coaseBoxID
      coarseBoxes[numCoarseBoxes].sendRank  = all_grids->levels[level  ]->rank_of_box[  fineBoxID];
      coarseBoxes[numCoarseBoxes].sendBoxID = fineBoxID;
      coarseBoxes[numCoarseBoxes].sendBox   = fineBox;
      coarseBoxes[numCoarseBoxes].recvRank  = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
      coarseBoxes[numCoarseBoxes].recvBoxID = coarseBoxID;
      coarseBoxes[numCoarseBoxes].recvBox   = coarseBox;  // -1 if off-node
      coarseBoxes[numCoarseBoxes].i         = (all_grids->levels[level]->box_dim/2)*( fineBox_i % (all_grids->levels[level]->boxes_in.i/all_grids->levels[level+1]->boxes_in.i) );
      coarseBoxes[numCoarseBoxes].j         = (all_grids->levels[level]->box_dim/2)*( fineBox_j % (all_grids->levels[level]->boxes_in.j/all_grids->levels[level+1]->boxes_in.j) );
      coarseBoxes[numCoarseBoxes].k         = (all_grids->levels[level]->box_dim/2)*( fineBox_k % (all_grids->levels[level]->boxes_in.k/all_grids->levels[level+1]->boxes_in.k) );
                  numCoarseBoxes++;
      if(all_grids->levels[level]->my_rank != all_grids->levels[level+1]->rank_of_box[coarseBoxID]){
        coarseRanks[numCoarseBoxesRemote++] = all_grids->levels[level+1]->rank_of_box[coarseBoxID];
      }else{numCoarseBoxesLocal++;}
    } // my (fine) boxes

    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(coarseBoxes,numCoarseBoxes      ,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(coarseRanks,numCoarseBoxesRemote,sizeof(    int),qsortInt);
    int numCoarseRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numCoarseBoxesRemote;neighbor++)if(coarseRanks[neighbor] != _rank){_rank=coarseRanks[neighbor];coarseRanks[numCoarseRanks++]=coarseRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->restriction[restrictionType].num_sends     =                         numCoarseRanks;
    all_grids->levels[level]->restriction[restrictionType].send_ranks    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].send_sizes    =            (int*)malloc(numCoarseRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].send_buffers  =        (double**)malloc(numCoarseRanks*sizeof(double*));
    if(numCoarseRanks>0){
    if(all_grids->levels[level]->restriction[restrictionType].send_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].send_ranks\n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].send_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].send_sizes\n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].send_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].send_buffers\n",level);exit(0);}
    }

    int elementSize;
    int restrict_dim_i=-1;
    int restrict_dim_j=-1;
    int restrict_dim_k=-1;
    switch(restrictionType){
      case RESTRICT_CELL   : restrict_dim_i = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level]->box_dim/2);break;
      case RESTRICT_FACE_I : restrict_dim_i = (1+all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level]->box_dim/2);break;
      case RESTRICT_FACE_J : restrict_dim_i = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (1+all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level]->box_dim/2);break;
      case RESTRICT_FACE_K : restrict_dim_i = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level]->box_dim/2);
                             restrict_dim_k = (1+all_grids->levels[level]->box_dim/2);break;
    }
    elementSize = restrict_dim_i*restrict_dim_j*restrict_dim_k;
   
    double * all_send_buffers = (double*)malloc(numCoarseBoxes*elementSize*sizeof(double));
          if(numCoarseBoxes*elementSize>0)
          if(all_send_buffers==NULL){fprintf(stderr,"malloc failed - restriction/all_send_buffers\n");exit(0);}
                      memset(all_send_buffers,0,numCoarseBoxes*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve

    // for each neighbor, construct the pack list and allocate the MPI send buffer... 
    for(neighbor=0;neighbor<numCoarseRanks;neighbor++){
      int coarseBox;
      int offset = 0;
      all_grids->levels[level]->restriction[restrictionType].send_buffers[neighbor] = all_send_buffers;
      for(coarseBox=0;coarseBox<numCoarseBoxes;coarseBox++)if(coarseBoxes[coarseBox].recvRank==coarseRanks[neighbor]){
        // restrict to MPI send buffer...
        append_block_to_list( &(all_grids->levels[level]->restriction[restrictionType].blocks[0]),
                              &(all_grids->levels[level]->restriction[restrictionType].allocated_blocks[0]),
                              &(all_grids->levels[level]->restriction[restrictionType].num_blocks[0]),
          /* dim.i         = */ restrict_dim_i, 
          /* dim.j         = */ restrict_dim_j, 
          /* dim.k         = */ restrict_dim_k, 
          /* read.box      = */ coarseBoxes[coarseBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ 0,
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].kStride,
          /* read.scale    = */ 2,
          /* write.box     = */ -1,
          /* write.ptr     = */ all_grids->levels[level]->restriction[restrictionType].send_buffers[neighbor],
          /* write.i       = */ offset,
          /* write.j       = */ 0,
          /* write.k       = */ 0,
          /* write.jStride = */ restrict_dim_i,
          /* write.kStride = */ restrict_dim_i*restrict_dim_j, 
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->restriction[restrictionType].send_ranks[neighbor] = coarseRanks[neighbor];
      all_grids->levels[level]->restriction[restrictionType].send_sizes[neighbor] = offset;
      all_send_buffers+=offset;
    }
    // for construct the local restriction list... 
    {
      int coarseBox;
      for(coarseBox=0;coarseBox<numCoarseBoxes;coarseBox++)if(coarseBoxes[coarseBox].recvRank==all_grids->levels[level+1]->my_rank){
        // restrict to local...
        append_block_to_list( &(all_grids->levels[level]->restriction[restrictionType].blocks[1]),
                              &(all_grids->levels[level]->restriction[restrictionType].allocated_blocks[1]),
                              &(all_grids->levels[level]->restriction[restrictionType].num_blocks[1]),
          /* dim.i         = */ restrict_dim_i, 
          /* dim.j         = */ restrict_dim_j, 
          /* dim.k         = */ restrict_dim_k, 
          /* read.box      = */ coarseBoxes[coarseBox].sendBox,
          /* read.ptr      = */ NULL,
          /* read.i        = */ 0, 
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].jStride,
          /* read.kStride  = */ all_grids->levels[level]->my_boxes[coarseBoxes[coarseBox].sendBox].kStride,
          /* read.scale    = */ 2,
          /* write.box     = */ coarseBoxes[coarseBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ coarseBoxes[coarseBox].i,
          /* write.j       = */ coarseBoxes[coarseBox].j,
          /* write.k       = */ coarseBoxes[coarseBox].k,
          /* write.jStride = */ all_grids->levels[level+1]->my_boxes[coarseBoxes[coarseBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level+1]->my_boxes[coarseBoxes[coarseBox].recvBox].kStride,
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
      }
    } // local to local

    // free temporary storage...
    free(coarseBoxes);
    free(coarseRanks);
  } // send/pack/local




  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // construct recv and unpack...
  if( (level>0) && (all_grids->levels[level]->num_my_boxes>0) ){ // not top  *and*  I have boxes to receive
    // construct a list of fine boxes to be coarsened and sent to me...
    int numFineBoxesMax = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*
                          (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*
                          (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*
                                                                  all_grids->levels[level]->num_my_boxes;
        int *fineRanks = (    int*)malloc(numFineBoxesMax*sizeof(    int)); // high water mark (assumes every neighboring box is a different process)
    RP_type *fineBoxes = (RP_type*)malloc(numFineBoxesMax*sizeof(RP_type)); 
    int numFineBoxesRemote = 0;
    int coarseBox;
    for(coarseBox=0;coarseBox<all_grids->levels[level]->num_my_boxes;coarseBox++){
      int bi,bj,bk;
      int   coarseBoxID = all_grids->levels[level]->my_boxes[coarseBox].global_box_id;
      int   coarseBox_i = all_grids->levels[level]->my_boxes[coarseBox].low.i / all_grids->levels[level]->box_dim;
      int   coarseBox_j = all_grids->levels[level]->my_boxes[coarseBox].low.j / all_grids->levels[level]->box_dim;
      int   coarseBox_k = all_grids->levels[level]->my_boxes[coarseBox].low.k / all_grids->levels[level]->box_dim;
      for(bk=0;bk<all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k;bk++){
      for(bj=0;bj<all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j;bj++){
      for(bi=0;bi<all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i;bi++){
        int fineBox_i = (all_grids->levels[level-1]->boxes_in.i/all_grids->levels[level]->boxes_in.i)*coarseBox_i + bi;
        int fineBox_j = (all_grids->levels[level-1]->boxes_in.j/all_grids->levels[level]->boxes_in.j)*coarseBox_j + bj;
        int fineBox_k = (all_grids->levels[level-1]->boxes_in.k/all_grids->levels[level]->boxes_in.k)*coarseBox_k + bk;
        int fineBoxID =  fineBox_i + fineBox_j*all_grids->levels[level-1]->boxes_in.i + fineBox_k*all_grids->levels[level-1]->boxes_in.i*all_grids->levels[level-1]->boxes_in.j;
        if(all_grids->levels[level-1]->rank_of_box[fineBoxID] != all_grids->levels[level]->my_rank){
          fineBoxes[numFineBoxesRemote].sendRank  = all_grids->levels[level-1]->rank_of_box[  fineBoxID];
          fineBoxes[numFineBoxesRemote].sendBoxID = fineBoxID;
          fineBoxes[numFineBoxesRemote].sendBox   = -1; // I don't know the off-node box index
          fineBoxes[numFineBoxesRemote].recvRank  = all_grids->levels[level  ]->rank_of_box[coarseBoxID];
          fineBoxes[numFineBoxesRemote].recvBoxID = coarseBoxID;
          fineBoxes[numFineBoxesRemote].recvBox   = coarseBox;
          fineBoxes[numFineBoxesRemote].i         = bi*all_grids->levels[level-1]->box_dim/2;
          fineBoxes[numFineBoxesRemote].j         = bj*all_grids->levels[level-1]->box_dim/2;
          fineBoxes[numFineBoxesRemote].k         = bk*all_grids->levels[level-1]->box_dim/2;
          fineRanks[numFineBoxesRemote] = all_grids->levels[level-1]->rank_of_box[fineBoxID];
                    numFineBoxesRemote++;
        }
      }}}
    } // my (coarse) boxes
    // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
    qsort(fineBoxes,numFineBoxesRemote,sizeof(RP_type),qsortRP );
    // sort the lists of neighboring ranks and remove duplicates...
    qsort(fineRanks,numFineBoxesRemote,sizeof(    int),qsortInt);
    int numFineRanks=0;
    int _rank=-1;int neighbor=0;
    for(neighbor=0;neighbor<numFineBoxesRemote;neighbor++)if(fineRanks[neighbor] != _rank){_rank=fineRanks[neighbor];fineRanks[numFineRanks++]=fineRanks[neighbor];}

    // allocate structures...
    all_grids->levels[level]->restriction[restrictionType].num_recvs     =                         numFineRanks;
    all_grids->levels[level]->restriction[restrictionType].recv_ranks    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].recv_sizes    =            (int*)malloc(numFineRanks*sizeof(int));
    all_grids->levels[level]->restriction[restrictionType].recv_buffers  =        (double**)malloc(numFineRanks*sizeof(double*));
    if(numFineRanks>0){
    if(all_grids->levels[level]->restriction[restrictionType].recv_ranks  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].recv_ranks  \n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].recv_sizes  ==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].recv_sizes  \n",level);exit(0);}
    if(all_grids->levels[level]->restriction[restrictionType].recv_buffers==NULL){fprintf(stderr,"malloc failed - all_grids->levels[%d]->restriction[restrictionType].recv_buffers\n",level);exit(0);}
    }

    int elementSize;
    int restrict_dim_i=-1;
    int restrict_dim_j=-1;
    int restrict_dim_k=-1;
    switch(restrictionType){
      case RESTRICT_CELL   : restrict_dim_i = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level-1]->box_dim/2);break;
      case RESTRICT_FACE_I : restrict_dim_i = (1+all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level-1]->box_dim/2);break;
      case RESTRICT_FACE_J : restrict_dim_i = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (1+all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (  all_grids->levels[level-1]->box_dim/2);break;
      case RESTRICT_FACE_K : restrict_dim_i = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_j = (  all_grids->levels[level-1]->box_dim/2);
                             restrict_dim_k = (1+all_grids->levels[level-1]->box_dim/2);break;
    }
    elementSize = restrict_dim_i*restrict_dim_j*restrict_dim_k;

    double * all_recv_buffers = (double*)malloc(numFineBoxesRemote*elementSize*sizeof(double));
          if(numFineBoxesRemote*elementSize>0)
          if(all_recv_buffers==NULL){fprintf(stderr,"malloc failed - restriction/all_recv_buffers\n");exit(0);}
                      memset(all_recv_buffers,0,numFineBoxesRemote*elementSize*sizeof(double)); // DO NOT DELETE... you must initialize to 0 to avoid getting something like 0.0*NaN and corrupting the solve
    //printf("level=%d, rank=%2d, recv_buffers=%6d\n",level,all_grids->my_rank,numFineBoxesRemote*elementSize*sizeof(double));

    // for each neighbor, construct the unpack list and allocate the MPI recv buffer... 
    for(neighbor=0;neighbor<numFineRanks;neighbor++){
      int fineBox;
      int offset = 0;
      all_grids->levels[level]->restriction[restrictionType].recv_buffers[neighbor] = all_recv_buffers;
      for(fineBox=0;fineBox<numFineBoxesRemote;fineBox++)if(fineBoxes[fineBox].sendRank==fineRanks[neighbor]){
        // unpack MPI recv buffer...
        append_block_to_list( &(all_grids->levels[level]->restriction[restrictionType].blocks[2]),
                              &(all_grids->levels[level]->restriction[restrictionType].allocated_blocks[2]),
                              &(all_grids->levels[level]->restriction[restrictionType].num_blocks[2]),
          /* dim.i         = */ restrict_dim_i, 
          /* dim.j         = */ restrict_dim_j, 
          /* dim.k         = */ restrict_dim_k, 
          /* read.box      = */ -1,
          /* read.ptr      = */ all_grids->levels[level]->restriction[restrictionType].recv_buffers[neighbor],
          /* read.i        = */ offset,
          /* read.j        = */ 0,
          /* read.k        = */ 0,
          /* read.jStride  = */ restrict_dim_i,
          /* read.kStride  = */ restrict_dim_i*restrict_dim_j, 
          /* read.scale    = */ 1,
          /* write.box     = */ fineBoxes[fineBox].recvBox,
          /* write.ptr     = */ NULL,
          /* write.i       = */ fineBoxes[fineBox].i,
          /* write.j       = */ fineBoxes[fineBox].j,
          /* write.k       = */ fineBoxes[fineBox].k,
          /* write.jStride = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].recvBox].jStride,
          /* write.kStride = */ all_grids->levels[level]->my_boxes[fineBoxes[fineBox].recvBox].kStride,
          /* write.scale   = */ 1,
          /* blockcopy_i   = */ BLOCKCOPY_TILE_I, // default
          /* blockcopy_j   = */ BLOCKCOPY_TILE_J, // default
          /* blockcopy_k   = */ BLOCKCOPY_TILE_K, // default
          /* subtype       = */ 0
        );
        offset+=elementSize;
      }
      all_grids->levels[level]->restriction[restrictionType].recv_ranks[neighbor] = fineRanks[neighbor];
      all_grids->levels[level]->restriction[restrictionType].recv_sizes[neighbor] = offset;
      all_recv_buffers+=offset;
    } // neighbor

    // free temporary storage...
    free(fineBoxes);
    free(fineRanks);
  } // recv/unpack



  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //MPI_Barrier(MPI_COMM_WORLD);
  //if(all_grids->my_rank==0){printf("================================================================================\n");}
  //if(                         (level==all_grids->num_levels-2))print_communicator(2,all_grids->my_rank,level,&all_grids->levels[level]->restriction);
  //if((all_grids->my_rank==0)&&(level==all_grids->num_levels-1))print_communicator(1,all_grids->my_rank,level,&all_grids->levels[level]->restriction);
  //MPI_Barrier(MPI_COMM_WORLD);
  } // level loop


  #ifdef USE_MPI
  for(level=0;level<all_grids->num_levels;level++){
    all_grids->levels[level]->restriction[restrictionType].requests = NULL;
    all_grids->levels[level]->restriction[restrictionType].status   = NULL;
    if(level<all_grids->num_levels-1){ // bottom never calls restriction()
    // by convention, level_f allocates a combined array of requests for both level_f sends and level_c recvs...
    int nMessages = all_grids->levels[level+1]->restriction[restrictionType].num_recvs + all_grids->levels[level]->restriction[restrictionType].num_sends;
    all_grids->levels[level]->restriction[restrictionType].requests = (MPI_Request*)malloc(nMessages*sizeof(MPI_Request));
    all_grids->levels[level]->restriction[restrictionType].status   = (MPI_Status *)malloc(nMessages*sizeof(MPI_Status ));
    }
  }
  #endif
}


//------------------------------------------------------------------------------------------------------------------------------
void MGBuild(mg_type *all_grids, level_type *fine_grid, double a, double b, int minCoarseGridDim){
  int  maxLevels=100;
  int     nProcs[100];
  int      dim_i[100];
  int boxes_in_i[100];
  int    box_dim[100];
  int box_ghosts[100];
  all_grids->my_rank = fine_grid->my_rank;
  all_grids->cycles.MGBuild = 0;
  uint64_t _timeStartMGBuild = CycleTime();

  // calculate how deep we can make the v-cycle...
  int level=1;
                             int coarse_dim = fine_grid->dim.i;
//if(fine_grid->dim.j<coarse_dim)coarse_dim = fine_grid->dim.j;
//if(fine_grid->dim.k<coarse_dim)coarse_dim = fine_grid->dim.k;
  while( (coarse_dim>=2*minCoarseGridDim) && ((coarse_dim&0x1)==0) ){ // grid dimension is even and big enough...
    level++;
    coarse_dim = coarse_dim / 2;
  }if(level<maxLevels)maxLevels=level;

      nProcs[0] = fine_grid->num_ranks;
       dim_i[0] = fine_grid->dim.i;
  boxes_in_i[0] = fine_grid->boxes_in.i;
     box_dim[0] = fine_grid->box_dim;
  box_ghosts[0] = fine_grid->box_ghosts;

  // build the list of levels...
  all_grids->levels = (level_type**)malloc(maxLevels*sizeof(level_type*));
  if(all_grids->levels == NULL){fprintf(stderr,"malloc failed - MGBuild/all_grids->levels\n");exit(0);}
  all_grids->num_levels=1;
  all_grids->levels[0] = fine_grid;


  // build a table to guide the construction of the v-cycle...
  int doRestrict=1;if(maxLevels<2)doRestrict=0; // i.e. can't restrict if there is only one level !!!
  #ifdef USE_UCYCLES
  while(doRestrict){
    level = all_grids->num_levels;
    doRestrict=0;
    if( (box_dim[level-1] % 2 == 0) ){
          nProcs[level] =     nProcs[level-1];
           dim_i[level] =      dim_i[level-1]/2;
         box_dim[level] =    box_dim[level-1]/2;
      boxes_in_i[level] = boxes_in_i[level-1];
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }
    if(box_dim[level] < box_ghosts[level])doRestrict=0;
    if(doRestrict)all_grids->num_levels++;
  }
  #else // TRUE V-Cycle...
  while(doRestrict){
    level = all_grids->num_levels;
    doRestrict=0;
    int fine_box_dim    =    box_dim[level-1];
    int fine_nProcs     =     nProcs[level-1];
    int fine_dim_i      =      dim_i[level-1];
    int fine_boxes_in_i = boxes_in_i[level-1];
    if( (fine_box_dim % 2 == 0) && (fine_box_dim > MG_AGGLOMERATION_START) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // Boxes are too big to agglomerate
          nProcs[level] = fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (fine_boxes_in_i % 2 == 0) && ((fine_box_dim)>=stencil_get_radius()) ){ // 8:1 box agglomeration
          nProcs[level] = fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim;
      boxes_in_i[level] = fine_boxes_in_i/2;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (coarse_dim != 1) && (fine_dim_i == 2*coarse_dim) && ((fine_dim_i/2)>=stencil_get_radius()) ){ // agglomerate everything
          nProcs[level] = 1;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_dim_i/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = 1;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (coarse_dim != 1) && (fine_dim_i == 4*coarse_dim) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // restrict box dimension, and run on fewer ranks
          nProcs[level] = coarse_dim<fine_nProcs ? coarse_dim : fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (coarse_dim != 1) && (fine_dim_i == 8*coarse_dim) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // restrict box dimension, and run on fewer ranks
          nProcs[level] = coarse_dim*coarse_dim<fine_nProcs ? coarse_dim*coarse_dim : fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }else
    if( (fine_box_dim % 2 == 0) && ((fine_box_dim/2)>=stencil_get_radius()) ){ // restrict box dimension, and run on the same number of ranks
          nProcs[level] = fine_nProcs;
           dim_i[level] = fine_dim_i/2;
         box_dim[level] = fine_box_dim/2; // FIX, verify its not less than the stencil radius
      boxes_in_i[level] = fine_boxes_in_i;
      box_ghosts[level] = box_ghosts[level-1];
             doRestrict = 1;
    }
    if(dim_i[level]<minCoarseGridDim)doRestrict=0;
    if(doRestrict)all_grids->num_levels++;
  }
  #endif
  //if(all_grids->my_rank==0){for(level=0;level<all_grids->num_levels;level++){
  //  printf("level %2d: %4d^3 using %4d^3.%4d^3 spread over %4d processes\n",level,dim_i[level],boxes_in_i[level],box_dim[level],nProcs[level]);
  //}}


  // now build all the coarsened levels...
  for(level=1;level<all_grids->num_levels;level++){
    all_grids->levels[level] = (level_type*)malloc(sizeof(level_type));
    if(all_grids->levels[level] == NULL){fprintf(stderr,"malloc failed - MGBuild/doRestrict\n");exit(0);}
    create_level(all_grids->levels[level],boxes_in_i[level],box_dim[level],box_ghosts[level],all_grids->levels[level-1]->numVectors,all_grids->levels[level-1]->boundary_condition.type,all_grids->levels[level-1]->my_rank,nProcs[level]);
    all_grids->levels[level]->h = 2.0*all_grids->levels[level-1]->h;
  }


  // bottom solver (level = all_grids->num_levels-1) gets extra grids...
  create_vectors(all_grids->levels[all_grids->num_levels-1],all_grids->levels[all_grids->num_levels-1]->numVectors + IterativeSolver_NumVectors() );

  // build the restriction and interpolation communicators...
  build_restriction(all_grids,RESTRICT_CELL  ); // cell-centered
  build_restriction(all_grids,RESTRICT_FACE_I); // face-centered, normal to i
  build_restriction(all_grids,RESTRICT_FACE_J); // face-centered, normal to j
  build_restriction(all_grids,RESTRICT_FACE_K); // face-centered, normal to k
  build_interpolation(all_grids);


  // build subcommunicators...
  #ifdef USE_MPI
  #ifdef USE_SUBCOMM
  if(all_grids->my_rank==0){fprintf(stdout,"\n");}
  for(level=1;level<all_grids->num_levels;level++){
    double comm_split_start = MPI_Wtime();
    if(all_grids->my_rank==0){fprintf(stdout,"  Building MPI subcommunicator for level %d...",level);fflush(stdout);}
    all_grids->levels[level]->active=0;
    int ll;for(ll=level;ll<all_grids->num_levels;ll++)if(all_grids->levels[ll]->num_my_boxes>0)all_grids->levels[level]->active=1;
    MPI_Comm_split(MPI_COMM_WORLD,all_grids->levels[level]->active,all_grids->levels[level]->my_rank,&all_grids->levels[level]->MPI_COMM_ALLREDUCE);
    double comm_split_end = MPI_Wtime();
    double comm_split_time_send = comm_split_end-comm_split_start;
    double comm_split_time = 0;
    MPI_Allreduce(&comm_split_time_send,&comm_split_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(all_grids->my_rank==0){fprintf(stdout,"done (%0.6f seconds)\n",comm_split_time);fflush(stdout);}
  }
  #endif
  #endif


  // rebuild various coefficients for the operator... must occur after build_restriction !!!
  if(all_grids->my_rank==0){fprintf(stdout,"\n");}
  for(level=1;level<all_grids->num_levels;level++){
    rebuild_operator(all_grids->levels[level],(level>0)?all_grids->levels[level-1]:NULL,a,b);
  }
  if(all_grids->my_rank==0){fprintf(stdout,"\n");}


  // used for quick test for poisson
  for(level=0;level<all_grids->num_levels;level++){
    all_grids->levels[level]->alpha_is_zero = (dot(all_grids->levels[level],VECTOR_ALPHA,VECTOR_ALPHA) == 0.0);
  }

  
  all_grids->cycles.MGBuild += (uint64_t)(CycleTime()-_timeStartMGBuild);
}


//------------------------------------------------------------------------------------------------------------------------------
void MGVCycle(mg_type *all_grids, int e_id, int R_id, double a, double b, int level){
  if(!all_grids->levels[level]->active)return;
  uint64_t _LevelStart;

  // bottom solve...
  if(level==all_grids->num_levels-1){
    uint64_t _timeBottomStart = CycleTime();
    IterativeSolver(all_grids->levels[level],e_id,R_id,a,b,MG_DEFAULT_BOTTOM_NORM);
    all_grids->levels[level]->cycles.Total += (uint64_t)(CycleTime()-_timeBottomStart);
    return;
  }

  // down...
  _LevelStart = CycleTime();
       smooth(all_grids->levels[level  ],e_id,R_id,a,b);
     residual(all_grids->levels[level  ],VECTOR_TEMP,e_id,R_id,a,b);
  restriction(all_grids->levels[level+1],R_id,all_grids->levels[level],VECTOR_TEMP,RESTRICT_CELL);
    zero_vector(all_grids->levels[level+1],e_id);
  all_grids->levels[level]->cycles.Total += (uint64_t)(CycleTime()-_LevelStart);

  // recursion...
  MGVCycle(all_grids,e_id,R_id,a,b,level+1);

  // up...
  _LevelStart = CycleTime();
  interpolation_vcycle(all_grids->levels[level  ],e_id,1.0,all_grids->levels[level+1],e_id);
         smooth(all_grids->levels[level  ],e_id,R_id,a,b);
  all_grids->levels[level]->cycles.Total += (uint64_t)(CycleTime()-_LevelStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void MGSolve(mg_type *all_grids, int u_id, int F_id, double a, double b, double dtol, double rtol){
  all_grids->MGSolves_performed++;
  if(!all_grids->levels[0]->active)return;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int e_id = u_id; // __u FIX
  int R_id = VECTOR_F_MINUS_AV;
  int v;
  int maxVCycles = 20;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_MPI
  double MG_Start_Time = MPI_Wtime();
  #endif
  if(all_grids->levels[0]->my_rank==0){fprintf(stdout,"MGSolve... ");}
  uint64_t _timeStartMGSolve = CycleTime();

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  double norm_of_F     = 1.0;
  double norm_of_DinvF = 1.0;
  if(dtol>0){
    mul_vectors(all_grids->levels[0],VECTOR_TEMP,1.0,F_id,VECTOR_DINV); // D^{-1}F
    norm_of_DinvF = norm(all_grids->levels[0],VECTOR_TEMP);		// ||D^{-1}F||
  }
  if(rtol>0)norm_of_F = norm(all_grids->levels[0],F_id);		// ||F||

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // make initial guess for e (=0) and setup the RHS
   zero_vector(all_grids->levels[0],e_id);                  // ee = 0
  scale_vector(all_grids->levels[0],R_id,1.0,F_id);         // R_id = F_id

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // now do v-cycles to calculate the correction...
  for(v=0;v<maxVCycles;v++){   
    int level = 0;
    all_grids->levels[level]->vcycles_from_this_level++;

    // do the v-cycle...
    MGVCycle(all_grids,e_id,R_id,a,b,level);

    // now calculate the norm of the residual...
    uint64_t _timeStart = CycleTime();
    if( (all_grids->levels[level]->boundary_condition.type==BC_PERIODIC) && ((a==0) || (all_grids->levels[level]->alpha_is_zero==1)) ){
      // Poisson with Periodic Boundary Conditions... by convention, we assume the solution sums to zero... so eliminate any constants from the solution...
      double average_value_of_e = mean(all_grids->levels[level],e_id);
      shift_vector(all_grids->levels[level],e_id,e_id,-average_value_of_e);
    }
    residual(all_grids->levels[level],VECTOR_TEMP,e_id,F_id,a,b);
    if(dtol>0)mul_vectors(all_grids->levels[level],VECTOR_TEMP,1.0,VECTOR_TEMP,VECTOR_DINV); //  Using ||D^{-1}(b-Ax)||_{inf} as convergence criteria...
    double norm_of_residual = norm(all_grids->levels[level],VECTOR_TEMP);
    uint64_t _timeNorm = CycleTime();
    all_grids->levels[level]->cycles.Total += (uint64_t)(_timeNorm-_timeStart);
    if(all_grids->levels[level]->my_rank==0){
      if(   v>0){fprintf(stdout,"\n           ");}
      if(rtol>0){fprintf(stdout,"v-cycle=%2d  norm=%1.15e  rel=%1.15e  ",v+1,norm_of_residual,norm_of_residual/norm_of_F    );}
            else{fprintf(stdout,"v-cycle=%2d  norm=%1.15e  rel=%1.15e  ",v+1,norm_of_residual,norm_of_residual/norm_of_DinvF);}
    }
    if(norm_of_residual/norm_of_F < rtol)break;
    if(norm_of_residual           < dtol)break;
  } // maxVCycles
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  all_grids->cycles.MGSolve += (uint64_t)(CycleTime()-_timeStartMGSolve);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_MPI
  if(all_grids->levels[0]->my_rank==0){fprintf(stdout,"done (%f seconds)\n",MPI_Wtime()-MG_Start_Time);} // used to monitor variability in individual solve times
  #else
  if(all_grids->levels[0]->my_rank==0){fprintf(stdout,"done\n");}
  #endif
}


//------------------------------------------------------------------------------------------------------------------------------
void FMGSolve(mg_type *all_grids, int u_id, int F_id, double a, double b, double dtol, double rtol){
  all_grids->MGSolves_performed++;
  if(!all_grids->levels[0]->active)return;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int maxVCycles=0;
  int v;
  int level;
  int e_id = u_id;
  int R_id = VECTOR_F_MINUS_AV;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_MPI
  double FMG_Start_Time = MPI_Wtime();
  #endif
  if(all_grids->levels[0]->my_rank==0){fprintf(stdout,"FMGSolve... ");}
  uint64_t _timeStartMGSolve = CycleTime();

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  double norm_of_F     = 1.0;
  double norm_of_DinvF = 1.0;
  if(dtol>0){
    mul_vectors(all_grids->levels[0],VECTOR_TEMP,1.0,F_id,VECTOR_DINV); // D^{-1}F
    norm_of_DinvF = norm(all_grids->levels[0],VECTOR_TEMP);		// ||D^{-1}F||
  }
  if(rtol>0)norm_of_F = norm(all_grids->levels[0],F_id);		// ||F||

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // initialize the RHS for the f-cycle to f...
  uint64_t _LevelStart = CycleTime();
  level=0;
  scale_vector(all_grids->levels[0],R_id,1.0,F_id);              // R_id = F_id
  all_grids->levels[0]->cycles.Total += (uint64_t)(CycleTime()-_LevelStart);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


  // restrict RHS to bottom (coarsest grids)
  for(level=0;level<(all_grids->num_levels-1);level++){
    uint64_t _LevelStart = CycleTime();
    restriction(all_grids->levels[level+1],R_id,all_grids->levels[level],R_id,RESTRICT_CELL);
    all_grids->levels[level]->cycles.Total += (uint64_t)(CycleTime()-_LevelStart);
  }


  // solve coarsest grid...
    uint64_t _timeBottomStart = CycleTime();
    level = all_grids->num_levels-1;
    if(level>0)zero_vector(all_grids->levels[level],e_id);//else use whatever was the initial guess
    IterativeSolver(all_grids->levels[level],e_id,R_id,a,b,MG_DEFAULT_BOTTOM_NORM);  // -1 == exact solution
    all_grids->levels[level]->cycles.Total += (uint64_t)(CycleTime()-_timeBottomStart);


  // now do the F-cycle proper...
  for(level=all_grids->num_levels-2;level>=0;level--){
    // high-order interpolation
    _LevelStart = CycleTime();
    interpolation_fcycle(all_grids->levels[level],e_id,0.0,all_grids->levels[level+1],e_id);
    all_grids->levels[level]->cycles.Total += (uint64_t)(CycleTime()-_LevelStart);

    // v-cycle
    all_grids->levels[level]->vcycles_from_this_level++;
    MGVCycle(all_grids,e_id,R_id,a,b,level);
  }


  // now do the post-F V-cycles
  for(v=-1;v<maxVCycles;v++){
    int level = 0;

    // do the v-cycle...
    if(v>=0){
    all_grids->levels[level]->vcycles_from_this_level++;
    MGVCycle(all_grids,e_id,R_id,a,b,level);
    }

    // now calculate the norm of the residual...
    uint64_t _timeStart = CycleTime();
    if( (all_grids->levels[level]->boundary_condition.type==BC_PERIODIC) && ((a==0) || (all_grids->levels[level]->alpha_is_zero==1)) ){
      // Poisson with Periodic Boundary Conditions... by convention, we assume the solution sums to zero... so eliminate any constants from the solution...
      double average_value_of_e = mean(all_grids->levels[level],e_id);
      shift_vector(all_grids->levels[level],e_id,e_id,-average_value_of_e);
    }
    residual(all_grids->levels[level],VECTOR_TEMP,e_id,F_id,a,b);
    if(dtol>0)mul_vectors(all_grids->levels[level],VECTOR_TEMP,1.0,VECTOR_TEMP,VECTOR_DINV); //  Using ||D^{-1}(b-Ax)||_{inf} as convergence criteria...
    double norm_of_residual = norm(all_grids->levels[level],VECTOR_TEMP);
    uint64_t _timeNorm = CycleTime();
    all_grids->levels[level]->cycles.Total += (uint64_t)(_timeNorm-_timeStart);
    if(all_grids->levels[level]->my_rank==0){
      if(  v>=0){fprintf(stdout,"\n            v-cycle=%2d  ",v+1);}
            else{fprintf(stdout,              "f-cycle     ");}
      if(rtol>0){fprintf(stdout,"norm=%1.15e  rel=%1.15e  ",norm_of_residual,norm_of_residual/norm_of_F    );}
            else{fprintf(stdout,"norm=%1.15e  rel=%1.15e  ",norm_of_residual,norm_of_residual/norm_of_DinvF);}
    }
    if(norm_of_residual/norm_of_F < rtol)break;
    if(norm_of_residual           < dtol)break;
  }


  
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  all_grids->cycles.MGSolve += (uint64_t)(CycleTime()-_timeStartMGSolve);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_MPI
  if(all_grids->levels[0]->my_rank==0){fprintf(stdout,"done (%f seconds)\n",MPI_Wtime()-FMG_Start_Time);} // used to monitor variability in individual solve times
  #else
  if(all_grids->levels[0]->my_rank==0){fprintf(stdout,"done\n");}
  #endif
}
//------------------------------------------------------------------------------------------------------------------------------
