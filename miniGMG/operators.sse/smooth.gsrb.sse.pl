eval 'exec perl $0 $*'
        if 0;
$NumArgs = $#ARGV+1;

$ARGV0 = $ARGV[0];
$ARGV0 =~ tr/A-Z/a-z/;


#==========================================================================================================================================
open(F,">./gsrb.c");
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
print F "// Samuel Williams\n";
print F "// SWWilliams\@lbl.gov\n";
print F "// Lawrence Berkeley National Lab\n";
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
print F "#include <stdint.h>\n";
print F "#include \"../timer.h\"\n";
print F "#define MAX(a, b)  (((a) > (b)) ? (a) : (b))\n";
print F "#define MIN(a, b)  (((a) < (b)) ? (a) : (b))\n";
print F "#define MAX_THREADS 256\n";
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
$SpinBarrier = 1;
$GlobalU = 10;
&smooth_gsrb_wavefront($GlobalU,0);
&smooth_gsrb_wavefront($GlobalU,1);
print F "//==================================================================================================\n";
print F "void smooth(domain_type * domain, int level, int phi_id, int rhs_id, double a, double b){\n";
print F "  int CollaborativeThreadingBoxSize = 100000; // i.e. never\n";
print F "  #ifdef __COLLABORATIVE_THREADING\n";
print F "    CollaborativeThreadingBoxSize = 1 << __COLLABORATIVE_THREADING;\n";
print F "  #endif\n";

print F "  uint64_t _timeStart = CycleTime();\n";
print F "  int box,s;\n";
print F "  int ghosts = domain->ghosts;\n";
print F "  // if communication-avoiding, need RHS for stencils in ghost zones\n";
print F "  if(ghosts>1)exchange_boundary(domain,level,rhs_id,1,1,1);\n";

print F "  for(s=0;s<numSmooths;s+=ghosts){\n";
print F "    exchange_boundary(domain,level,phi_id,1,ghosts>1,ghosts>1);  // corners/edges if doing communication-avoiding...\n";
print F "    if(domain->subdomains[0].levels[level].dim.i >= CollaborativeThreadingBoxSize){\n";
print F "      uint64_t _timeStart = CycleTime();\n";
print F "      for(box=0;box<domain->subdomains_per_rank;box++){__box_smooth_GSRB_multiple_threaded(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,s);}\n";
print F "      domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);\n";
print F "    }else{\n";
print F "      // now do ghosts communication-avoiding smooths on each box...\n";
print F "      uint64_t _timeStart = CycleTime();\n";
print F "      #pragma omp parallel for private(box)\n";
print F "      for(box=0;box<domain->subdomains_per_rank;box++){__box_smooth_GSRB_multiple(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,s);}\n";
print F "      domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);\n";
print F "    }\n";
print F "  }\n";

print F "}\n";
print F "//==================================================================================================\n";
close(F);


#==========================================================================================================================================
sub smooth_gsrb_wavefront{
  local($U,$THREADED)=@_;
  $Uminus1 = $U-1;
  $Alignment = 2; # i.e. align to multiples of four...
  $AlignmentMinus1 = $Alignment-1;
  if(($U % $Alignment) != 0){printf("Warning, SSE code's unrolling must be a multiple of 4\n");return;}
  if($THREADED){$FunctionName = "__box_smooth_GSRB_multiple_threaded";}
           else{$FunctionName = "__box_smooth_GSRB_multiple";}
  print F "void $FunctionName(box_type *box, int phi_id, int rhs_id, double a, double b, int sweep){\n";

  if($THREADED && $SpinBarrier){
  print F "  volatile int64_t KPlaneFinishedByThread[MAX_THREADS];\n";
  }

  if($THREADED){
  if($SpinBarrier){print F "  #pragma omp parallel shared(KPlaneFinishedByThread)\n";}
              else{print F "  #pragma omp parallel\n";}
  print F "  {\n";
  print F "    int pencil = box->pencil;\n";
  print F "    int plane = box->plane;\n";
  print F "    int ghosts = box->ghosts;\n";
  print F "    int DimI = box->dim.i;\n";
  print F "    int DimJ = box->dim.j;\n";
  print F "    int DimK = box->dim.k;\n";
  print F "    double h2inv = 1.0/(box->h*box->h);\n";
  print F "    double * __restrict__ phi    = box->grids[  phi_id] + ghosts*plane;\n";
  print F "    double * __restrict__ rhs    = box->grids[  rhs_id] + ghosts*plane;\n";
  print F "    double * __restrict__ alpha  = box->grids[__alpha ] + ghosts*plane;\n";
  print F "    double * __restrict__ beta_i = box->grids[__beta_i] + ghosts*plane;\n";
  print F "    double * __restrict__ beta_j = box->grids[__beta_j] + ghosts*plane;\n";
  print F "    double * __restrict__ beta_k = box->grids[__beta_k] + ghosts*plane;\n";
  print F "    double * __restrict__ lambda = box->grids[__lambda] + ghosts*plane;\n";
  print F "    const __m128d       a_splat2 =            _mm_loaddup_pd(&a);\n";
  print F "    const __m128d b_h2inv_splat2 = _mm_mul_pd(_mm_loaddup_pd(&b),_mm_loaddup_pd(&h2inv));\n";
  print F "    int id      = omp_get_thread_num();\n";
  print F "    int threads = omp_get_num_threads();\n";

  if($SpinBarrier){
  print F "    // only works if (ij_end-ij_start)>=pencil;\n";
  print F "    int  left = MAX(        0,id-1);\n";
  print F "    int right = MIN(threads-1,id+1);\n";
  print F "    if(ghosts==1){right=id;left=id;}\n";
  print F "    if(ghosts>1){\n";
  print F "      KPlaneFinishedByThread[id]=-100;\n";
  print F "      #pragma omp barrier\n";
  print F "    }\n";
  }

  print F "    int global_ij_start[8];\n";
  print F "    int global_ij_end[8];\n";
  print F "    int ij_start[8];\n";
  print F "    int ij_end[8];\n";
  print F "    int planeInWavefront=0;\n";
  print F "      global_ij_start[planeInWavefront] = (                   (1)*pencil)&~$AlignmentMinus1;\n";
  print F "      global_ij_end[  planeInWavefront] = ((ghosts+DimJ+ghosts-1)*pencil);\n";
  print F "      int TotalUnrollings = ((global_ij_end[planeInWavefront]-global_ij_start[planeInWavefront]+$U-1)/$U);\n";
  print F "      ij_start[planeInWavefront] = global_ij_start[planeInWavefront] + $U*( (id          )*(TotalUnrollings)/(threads));\n";
  print F "      ij_end[  planeInWavefront] = global_ij_start[planeInWavefront] + $U*( (id+1        )*(TotalUnrollings)/(threads));\n";
  print F "      if(ij_end[planeInWavefront]>global_ij_end[planeInWavefront])ij_end[planeInWavefront]=global_ij_end[planeInWavefront];\n";
  print F "    for(planeInWavefront=1;planeInWavefront<ghosts;planeInWavefront++){\n";
  print F "      ij_start[planeInWavefront] = ij_start[0];\n";
  print F "      ij_end[  planeInWavefront] = ij_end[0];\n";
  print F "    }\n";

  }else{ # sequential version...
 #print F "  {\n";
  print F "    int pencil = box->pencil;\n";
  print F "    int plane = box->plane;\n";
  print F "    int ghosts = box->ghosts;\n";
  print F "    int DimI = box->dim.i;\n";
  print F "    int DimJ = box->dim.j;\n";
  print F "    int DimK = box->dim.k;\n";
  print F "    double h2inv = 1.0/(box->h*box->h);\n";
  print F "    double  * __restrict__ phi          = box->grids[  phi_id] + ghosts*plane;\n";
  print F "    double  * __restrict__ rhs          = box->grids[  rhs_id] + ghosts*plane;\n";
  print F "    double  * __restrict__ alpha        = box->grids[__alpha ] + ghosts*plane;\n";
  print F "    double  * __restrict__ beta_i       = box->grids[__beta_i] + ghosts*plane;\n";
  print F "    double  * __restrict__ beta_j       = box->grids[__beta_j] + ghosts*plane;\n";
  print F "    double  * __restrict__ beta_k       = box->grids[__beta_k] + ghosts*plane;\n";
  print F "    double  * __restrict__ lambda       = box->grids[__lambda] + ghosts*plane;\n";
  print F "    const __m128d       a_splat2 =            _mm_loaddup_pd(&a);\n";
  print F "    const __m128d b_h2inv_splat2 = _mm_mul_pd(_mm_loaddup_pd(&b),_mm_loaddup_pd(&h2inv));\n";
  print F "    int global_ij_start[8];\n";
  print F "    int global_ij_end[8];\n";
  print F "    int ij_start[8];\n";
  print F "    int ij_end[8];\n";
  print F "    int planeInWavefront;for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){\n";
  print F "      global_ij_start[planeInWavefront] = (                   (1+planeInWavefront)*pencil)&~$AlignmentMinus1;\n";
  print F "      global_ij_end[planeInWavefront]   = ((ghosts+DimJ+ghosts-1-planeInWavefront)*pencil);\n";
  print F "      ij_start[planeInWavefront] = global_ij_start[planeInWavefront];\n";
  print F "      ij_end[planeInWavefront]   = global_ij_end[planeInWavefront];\n";
  print F "    }\n";
  }


  $PF_Streams=0;
  print F "    #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)\n";
  print F "    double * __restrict__ DRAM_PREFETCH_POINTERS[20];\n";
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] =    phi+plane-pencil;\n";$PF_Streams++;
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] = beta_k+plane       ;\n";$PF_Streams++;
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] = beta_j             ;\n";$PF_Streams++;
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] = beta_i             ;\n";$PF_Streams++;
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] =  alpha             ;\n";$PF_Streams++;
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] =    rhs             ;\n";$PF_Streams++;
  print F "    DRAM_PREFETCH_POINTERS[$PF_Streams] = lambda             ;\n";$PF_Streams++;
  print F "    #endif\n";


  print F "    int leadingK;\n";
  print F "    int kLow  =     -(ghosts-1);\n";
  print F "    int kHigh = DimK+(ghosts-1);\n";
  print F "    for(leadingK=kLow;leadingK<kHigh;leadingK++){\n";
  print F "      #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)\n";
  print F "      int DRAM_prefetch_stream=0;\n";
  print F "      if(leadingK>=(kHigh-1))DRAM_prefetch_stream=$PF_Streams; // don't prefetch next plane when on last plane\n";
  print F "      int DRAM_prefetch_ijk_start = ij_start[0] + (leadingK+1)*plane;\n";
  print F "      int DRAM_prefetch_ijk_end   = ij_end[0]   + (leadingK+1)*plane;\n";
  print F "      int DRAM_prefetch_ijk       = DRAM_prefetch_ijk_start;\n";
  print F "      #endif\n";
  print F "      for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){\n";
  print F "        int k=(leadingK-planeInWavefront);\n";
  print F "        if(k>=kLow){\n";
                   &smooth_gsrb_VL_sse($GlobalU);
  print F "        } // active plane\n";
  print F "      }\n";
  if($THREADED){
  print F "      if(ghosts>1){\n";
  if($SpinBarrier){
  print F "        KPlaneFinishedByThread[id]=leadingK;\n";
  print F "        while( (KPlaneFinishedByThread[left ]<leadingK) || (KPlaneFinishedByThread[right]<leadingK) ){_mm_pause();}; // pause() in case HT is in use...\n";
  }else{
  print F "        #pragma omp barrier\n";
  }
  print F "      }\n";
  }

  print F "    } // leadingK\n";
  if($THREADED){
  print F "  } // omp parallel region\n";
  }
  print F "}\n\n\n";
}


#==========================================================================================================================================
sub smooth_gsrb_VL_sse{
  local($U)=@_;

                        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        print  F "        uint64_t invertMask = 0-((k^planeInWavefront^sweep^1)&0x1);\n";
                        print  F "        const __m128d invertMask2 = _mm_loaddup_pd((double*)&invertMask);\n";
                        print  F "        double * __restrict__ RedBlackFP = box->RedBlack_FP[(k^planeInWavefront^sweep^1)&0x1];\n";
                        print  F "        int kplane=k*plane;\n";
                        print  F "        int ij           = ij_start[planeInWavefront];\n";
                        print  F "        int _ij_end      = ij_end[  planeInWavefront];\n";
                        print  F "        int ijk=ij+kplane;\n";
                        print  F "        while(ij<_ij_end){ // smooth a vector...\n";
                        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        $PF_number    = int(($PF_Streams*$U+31)/32); # $PF_Streams (7) planes prefetched for every 4 planes computed
                        $PF_increment = int(($PF_Streams*$U   )/4);  # exact increment... e.g. increment by 14 == redundant prefetches
                       #$PF_increment = 8*$PF_number;                # increment matches number of prefetches == too many prefetches
                             print  F "          #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)\n";
                             print  F "          #warning will attempt to prefetch the next plane from DRAM one component at a time\n";
                             printf(F "          if(DRAM_prefetch_stream<$PF_Streams){\n");
                             printf(F "            double * _base = DRAM_PREFETCH_POINTERS[DRAM_prefetch_stream] + DRAM_prefetch_ijk;\n");
for($x=0;$x<$PF_number;$x++){printf(F "            _mm_prefetch((const char*)(_base+%2d),_MM_HINT_T1);\n",8*$x);}
                             printf(F "            DRAM_prefetch_ijk+=$PF_increment;\n");
                             printf(F "            if(DRAM_prefetch_ijk>DRAM_prefetch_ijk_end){DRAM_prefetch_stream++;DRAM_prefetch_ijk=DRAM_prefetch_ijk_start;}\n");
                             printf(F "          }\n");
                             print  F "          #endif\n";
                        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        printf(F "          //careful... assumes the compiler maps _mm128_load_pd to unaligned vmovupd and not the aligned version (should be faster when pencil is a multiple of 2 doubles (16 bytes)\n");
                        printf(F "          #if 1\n");
  for($x=0;$x<$U;$x+=2){printf(F "                __m128d helmholtz_%02d;\n",$x);}
                        printf(F "                                         //_mm_prefetch((const char*)(   phi+ijk+1+$U),_MM_HINT_T0);\n");
for($x=-2;$x<$U+4;$x+=2){ # pipelined loads/shuffles of phi
 if(($x< 0)           ){printf(F "          const __m128d         phi_m2 = _mm_load_pd(phi+ijk+%3d);\n",$x);}
 if(($x>=0)&&($x<$U+2)){printf(F "          const __m128d         phi_%02d = _mm_load_pd(phi+ijk+%3d);\n",$x,$x);}
 if(($x==2)           ){printf(F "          const __m128d         phi_m1 = _mm_shuffle_pd(phi_m2,phi_00,1);\n");}
 if(($x>=4)           ){printf(F "          const __m128d         phi_%02d = _mm_shuffle_pd(phi_%02d,phi_%02d,1);\n",$x-3,$x-4,$x-2);}}
                        printf(F "                                         //_mm_prefetch((const char*)(beta_i+ijk+1+$U),_MM_HINT_T0);\n");
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d =                           _mm_mul_pd(_mm_sub_pd(phi_%02d,phi_%02d),_mm_loadu_pd(beta_i+ijk+       %2d)); \n",$x,$x+1,$x,$x+1);}
                   $x=0;printf(F "                        helmholtz_%02d = _mm_sub_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(phi_%02d,phi_%02s),_mm_load_pd( beta_i+ijk+       %2d)));\n",$x,$x,$x,"m1",$x);
  for($x=2;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_sub_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(phi_%02d,phi_%02d),_mm_load_pd( beta_i+ijk+       %2d)));\n",$x,$x,$x,$x-1,$x);}
                        printf(F "          #else\n");
  for($x=0;$x<$U;$x+=2){printf(F "                __m128d helmholtz_%02d;\n",$x);}
                        printf(F "                                         //_mm_prefetch((const char*)(   phi+ijk+1+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(beta_i+ijk+1+$U),_MM_HINT_T0);\n");
                        printf(F "          // this version performs unalligned accesses for phi+/-1, betai+1 and phi+/-pencil\n",$x);
  for($x=0;$x<$U;$x+=2){printf(F "          const __m128d       phi_%02d = _mm_load_pd(phi+ijk+%3d);\n",$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+       %2d),             phi_%02d           ),_mm_loadu_pd(beta_i+ijk+       %2d)); \n",$x,$x+1,$x,$x+1);}
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_sub_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(             phi_%02d           ,_mm_loadu_pd(phi+ijk+       %2d)),_mm_load_pd( beta_i+ijk+       %2d)));\n",$x,$x,$x,$x-1,$x);}
                        printf(F "          #endif\n");
                        printf(F "                                         //_mm_prefetch((const char*)(   phi+ijk-pencil+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(   phi+ijk+pencil+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(beta_j+ijk       +$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(beta_j+ijk+pencil+$U),_MM_HINT_T0);\n");
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_add_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+%2d),             phi_%02d           ),_mm_load_pd( beta_j+ijk+pencil+%2d)));\n",$x,$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_sub_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(             phi_%02d           ,_mm_load_pd( phi+ijk-pencil+%2d)),_mm_load_pd( beta_j+ijk+       %2d)));\n",$x,$x,$x,$x,$x);}
                        printf(F "                                         //_mm_prefetch((const char*)(   phi+ijk-plane+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(   phi+ijk+plane+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(beta_k+ijk      +$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(beta_k+ijk+plane+$U),_MM_HINT_T0);\n");
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_add_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+%2d),             phi_%02d           ),_mm_load_pd( beta_k+ijk+ plane+%2d)));\n",$x,$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_sub_pd(helmholtz_%02d,_mm_mul_pd(_mm_sub_pd(             phi_%02d           ,_mm_load_pd( phi+ijk- plane+%2d)),_mm_load_pd( beta_k+ijk       +%2d)));\n",$x,$x,$x,$x,$x);}
                        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        print  F "          #ifdef __GSRB_FP\n";
                        print  F "          #warning GSRB using precomputed 64b FP array for Red-Black\n";
                        printf(F "                                         //_mm_prefetch((const char*)(      alpha+ijk+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(        rhs+ijk+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(     lambda+ijk+$U),_MM_HINT_T0);\n");
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_mul_pd(helmholtz_%02d,b_h2inv_splat2);\n",$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+%2d)),phi_%02d),helmholtz_%02d);\n",$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                __m128d       new_%02d = _mm_mul_pd(_mm_load_pd(lambda+ijk+%2d),_mm_sub_pd(helmholtz_%02d,_mm_load_pd(rhs+ijk+%2d)));\n",$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "          const __m128d  RedBlack_%02d = _mm_load_pd(RedBlackFP+ij+%2d);\n",$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                              new_%02d = _mm_sub_pd(phi_%02d,_mm_mul_pd(RedBlack_%02d,new_%02d));\n",$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                                              _mm_store_pd(phi+ijk+%2d,new_%02d);\n",$x,$x);}
                        print  F "          ij+=$U;\n";
                        print  F "          ijk+=$U;\n";
                        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        print  F "          #else\n";
                        print  F "          #warning GSRB using precomputed 64b integer mask array for Red-Black\n";
                        printf(F "                                         //_mm_prefetch((const char*)(      alpha+ijk+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(        rhs+ijk+$U),_MM_HINT_T0);\n");
                        printf(F "                                         //_mm_prefetch((const char*)(     lambda+ijk+$U),_MM_HINT_T0);\n");
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_mul_pd(helmholtz_%02d,b_h2inv_splat2);\n",$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                        helmholtz_%02d = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+%2d)),phi_%02d),helmholtz_%02d);\n",$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                __m128d       new_%02d = _mm_mul_pd(_mm_load_pd(lambda+ijk+%2d),_mm_sub_pd(helmholtz_%02d,_mm_load_pd(rhs+ijk+%2d)));\n",$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                              new_%02d = _mm_sub_pd(phi_%02d,new_%02d);\n",$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "          const __m128d RedBlack_%02d = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+%2d) ));\n",$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                        new_%02d  = _mm_or_pd(_mm_and_pd(RedBlack_%02d,new_%02d),_mm_andnot_pd(RedBlack_%02d,phi_%02d));\n",$x,$x,$x,$x,$x);}
  for($x=0;$x<$U;$x+=2){printf(F "                                              _mm_store_pd(phi+ijk+%2d,new_%02d);\n",$x,$x);}
                        print  F "          ij+=$U;\n";
                        print  F "          ijk+=$U;\n";
                        print  F "          #endif\n";
                        #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        print  F "        }\n";

return;



}
#==========================================================================================================================================
