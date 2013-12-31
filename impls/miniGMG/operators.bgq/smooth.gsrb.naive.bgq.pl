eval 'exec perl $0 $*'
        if 0;
$NumArgs = $#ARGV+1;

$ARGV0 = $ARGV[0];
$ARGV0 =~ tr/A-Z/a-z/;


#==========================================================================================================================================
$GlobalU = 4;
$INTERLEAVED_SMT     = 1;
$SMT_TOTAL_THREADS   = 4; # BGQ=4 
$SMT_STENCIL_THREADS = 3; 
$SMT_SCOUT_THREADS   = 0; # 0=none, 1=prefetch wavefront from L2 into L1
if($SMT_SCOUT_THREADS > 1){
  print "ERROR, SMT_SCOUT_THREADS($SMT_SCOUT_THREADS) exceeds the maximum number of scout threads per core (1)\n";
  exit(0);
}
if( ($SMT_STENCIL_THREADS + $SMT_SCOUT_THREADS) > $SMT_TOTAL_THREADS){
  print "ERROR, SMT_STENCIL_THREADS($SMT_STENCIL_THREADS) + SMT_SCOUT_THREADS($SMT_SCOUT_THREADS) exceeds the maximum number of threads per core ($SMT_TOTAL_THREADS)\n";
  exit(0);
}
#==========================================================================================================================================
open(F,">./smooth.c");
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
print F "// Samuel Williams\n";
print F "// SWWilliams\@lbl.gov\n";
print F "// Lawrence Berkeley National Lab\n";
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
print F "#include <stdint.h>\n";
print F "#include \"../timer.h\"\n";
print F "#include <spi/include/l1p/sprefetch.h>\n";
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
&smooth_gsrb($GlobalU,0);
&smooth_gsrb($GlobalU,1);
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
print F "void smooth(domain_type * domain, int level, int phi_id, int rhs_id, double a, double b, int sweep){\n";
print F "  uint64_t _timeStart = CycleTime();\n";
print F "  int CollaborativeThreadingBoxSize = 100000; // i.e. never\n";
print F "  #ifdef __COLLABORATIVE_THREADING\n";
print F "    CollaborativeThreadingBoxSize = 1 << __COLLABORATIVE_THREADING;\n";
print F "  #endif\n";
print F "  int box;\n";
print F "  if(domain->subdomains[0].levels[level].dim.i >= CollaborativeThreadingBoxSize){\n";
print F "    for(box=0;box<domain->numsubdomains;box++){__box_smooth_GSRB_multiple_threaded(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,hLevel,sweep);}\n";
print F "  }else{\n";
print F "    #pragma omp parallel for private(box)\n";
print F "    for(box=0;box<domain->numsubdomains;box++){__box_smooth_GSRB_multiple(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,hLevel,sweep);}\n";
print F "  }\n";
print F "  domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);\n";
print F "}\n";
print F "//------------------------------------------------------------------------------------------------------------------------------\n";
close(F);
 

#==========================================================================================================================================
sub smooth_gsrb{
  local($U,$THREADED)=@_;
  $Alignment = 4; # i.e. align to multiples of four...
  $AlignmentMinus1 = $Alignment-1;
  $Uminus1 = $U-1;
 #if(($U % $Alignment) != 0){printf("Warning, BGQ code's unrolling must be a multiple of 4\n");return;}
  if($THREADED){$FunctionName = "__box_smooth_GSRB_multiple_threaded";}
           else{$FunctionName = "__box_smooth_GSRB_multiple";}
  print F "void $FunctionName(box_type *box, int phi_id, int rhs_id, double a, double b, int sweep){\n";

if($THREADED){
  print F "  #pragma omp parallel\n";
  print F "  {\n";
  print F "    int s,color;\n";
  print F "    int pencil = box->pencil;\n";
  print F "    int plane = box->plane;\n";
  print F "    int plane8   = plane<<3;\n";
  print F "    int pencil8  = pencil<<3;\n";
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
  print F "    uint64_t* __restrict__ RedBlackMask = box->RedBlack_64bMask;\n";
  print F "        vector4double       a_splat4 = vec_splats(a);\n";
  print F "        vector4double b_h2inv_splat4 = vec_mul(vec_splats(b),vec_splats(h2inv));\n";
  print F "    int id      = omp_get_thread_num();\n";
  print F "    int threads = omp_get_num_threads();\n";
  print F "    int smt_id      = id &  0x3;\n";
  print F "    int smt_leader  = id & ~0x3;\n";

  print F "    int ij_start[8];\n";
  print F "    int ij_end[8];\n";
  print F "    int k_first[8];\n";
  print F "    int k_last[8];\n";
  print F "    for(s=0;s<ghosts;s++){\n";
  print F "      int TotalPencilsPerPlane = (DimJ+2*(ghosts-s-1));\n";
  print F "      int TotalPencilsPerSweep = (DimK+2*(ghosts-s-1))*TotalPencilsPerPlane;\n";
  if($INTERLEAVED_SMT){
  print F "      int pencil_start = ( (smt_leader  )*(TotalPencilsPerSweep)/(threads));\n";
  print F "      int pencil_end   = ( (smt_leader+4)*(TotalPencilsPerSweep)/(threads));\n";
  }else{
  print F "      int pencil_start = ( (id          )*(TotalPencilsPerSweep)/(threads));\n";
  print F "      int pencil_end   = ( (id+1        )*(TotalPencilsPerSweep)/(threads));\n";
  }
  print F "      ij_start[s] = (pencil*(((pencil_start) % TotalPencilsPerPlane)+1+s)+1+s)&~$AlignmentMinus1;\n";
  print F "      k_first[s]  =          ((pencil_start) / TotalPencilsPerPlane)-ghosts+1+s;\n";
  print F "      ij_end[s]   = (pencil*(((pencil_end-1) % TotalPencilsPerPlane)+2+s)+1+s)&~$AlignmentMinus1;\n";
  print F "      if( ((pencil_end-1) % TotalPencilsPerPlane) == (TotalPencilsPerPlane-1) ){\n";
  print F "      ij_end[s]   = (pencil*(ghosts+DimJ+ghosts-s-1)-s-1);\n";
  print F "      }\n";
  print F "      k_last[s]   =          ((pencil_end-1) / TotalPencilsPerPlane)-ghosts+1+s;\n";
  print F "    }\n";


  print F "    uint32_t old_stream_depth;L1P_GetStreamDepth(&old_stream_depth);\n";
  print F "    if(smt_id==0){\n";
 #print F "      L1P_SetStreamAdaptiveMode(0);\n";
 #print F "      L1P_SetStreamPolicy(L1P_stream_disable);\n";
 #print F "      L1P_SetStreamPolicy(L1P_stream_optimistic);\n";
 #print F "      L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);\n";
 #print F "      L1P_SetStreamPolicy(L1P_stream_confirmed);\n";
 #print F "      L1P_SetStreamDepth(2);\n";
  if($INTERLEAVED_SMT){
  print F "      L1P_SetStreamDepth(2);\n";
  }
  print F "    }\n";



}else{
  print F "    int s,color;\n";
  print F "    int pencil = box->pencil;\n";
  print F "    int plane = box->plane;\n";
  print F "    int plane8   = plane<<3;\n";
  print F "    int pencil8  = pencil<<3;\n";
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
  print F "    uint64_t* __restrict__ RedBlackMask = box->RedBlack_64bMask;\n";
  print F "        vector4double       a_splat4 = vec_splats(a);\n";
  print F "        vector4double b_h2inv_splat4 = vec_mul(vec_splats(b),vec_splats(h2inv));\n";
  print F "    int smt_id      = 0;\n";
  print F "    int smt_leader  = 0;\n";

  print F "    int ij_start[8];\n";
  print F "    int ij_end[8];\n";
  print F "    int k_first[8];\n";
  print F "    int k_last[8];\n";
  print F "    for(s=0;s<ghosts;s++){\n";
  print F "      ij_start[s] = (pencil*(                   1+s)+s+1)&~$AlignmentMinus1;\n";
  print F "      ij_end[s]   = (pencil*(ghosts+DimJ+ghosts-1-s)-s-1)&~$AlignmentMinus1;\n";
  print F "      k_first[s]  =            -ghosts+1+s;\n";
  print F "      k_last[s]   =        DimK+ghosts-1-s;\n";
  print F "    }\n";
}



  print F "    for(s=0,color=sweep;s<ghosts;s++,color++){\n";
if($THREADED){
  #print F "      printf(\"%3d:  sweep0:%4d,%3d...%4d,%3d\\n\",id,ij_start[0],k_first[0],ij_end[0],k_last[0]);\n";
  #print F "      printf(\"%3d:  sweep0:%4d,%3d...%4d,%3d  sweep1:%4d,%3d...%4d,%3d  sweep2:%4d,%3d...%4d,%3d  sweep3:%4d,%3d...%4d,%3d\\n\",id,ij_start[0],k_first[0],ij_end[0],k_last[0],ij_start[1],k_first[1],ij_end[1],k_last[1],ij_start[2],k_first[2],ij_end[2],k_last[2],ij_start[3],k_first[3],ij_end[3],k_last[3]);\n";
  print F "      if(ghosts>1){\n";
  print F "      #pragma omp barrier\n";
  print F "      }\n";
  #print F "      exit(0);\n";
}


  if($THREADED && $INTERLEAVED_SMT && $SMT_SCOUT_THREADS){
  print F "      if(smt_id == 3){\n";
  print F "      }else\n";
  }
  if($THREADED && $INTERLEAVED_SMT){
  print F "      if(smt_id<$SMT_STENCIL_THREADS){\n";
  }else{
  print F "      {\n";
  }
  print F "        int k;\n";
  print F "        for(k=k_first[s];k<=k_last[s];k++){\n";
  $GlobalS=$GlobalU;
  if($THREADED && $INTERLEAVED_SMT){$GlobalS=$SMT_STENCIL_THREADS*$GlobalU;}
  print F "          int ij,_ij_end;\n";
  print F "          if(k==k_first[s]){ij      = ij_start[s];}else{ij      = (pencil*(                   s+1)+s+1)&~$AlignmentMinus1;}\n";
  print F "          if(k==k_last[s] ){_ij_end = ij_end[s];}  else{_ij_end = (pencil*(ghosts+DimJ+ghosts-s-1)-s-1);}\n";
  print F "          ij = ij+$GlobalU*smt_id;\n";
  print F "          int ijk = ij+k*plane;\n";
                     &smooth_gsrb_VL_bgq($GlobalU,$GlobalS);
  print F "        }\n";
  print F "      } // stencils\n";



  print F "      } // sweep\n";
if($THREADED){
  print F "    if(smt_id==0)L1P_SetStreamDepth(old_stream_depth);\n";
  print F "    }\n";
}
  print F "  }\n\n";
  return;

}


#==========================================================================================================================================
sub smooth_gsrb_VL_bgq{
  local($U,$S)=@_;

                         print  F "        uint64_t invertMask = 0-((1^color^k)&0x1);\n";
                         print  F "        vector4double invertMask4 = vec_splats(*((double*)&invertMask));\n";

                         print  F "        while(ij<_ij_end){ // smooth an interleaved vector...\n";


                         printf(F "                                                         __dcbtst((void*)(   phi+ijk+$S));\n");
  # load the requisite data for the laplacian...
                   $x=-4;printf(F "                vector4double                 phi_m4   = vec_lda(%3d,(double*)(   phi+ijk));\n",8*$x);
  for($x=0;$x<=$U;$x+=4){printf(F "                vector4double                 phi_%02d = vec_lda(%3d,(double*)(   phi+ijk));\n",$x,8*$x);}
  for($x=0;$x<=$U;$x+=4){printf(F "                vector4double              beta_i_%02d = vec_lda(%3d,(double*)(beta_i+ijk));\n",$x,8*$x);}
                    $x=0;printf(F "                vector4double                 phi_m1   = vec_sldw(   phi_m4  ,   phi_%02d,3);\n",$x);
  for($x=4;$x< $U;$x+=4){printf(F "                vector4double                 phi_%02d = vec_sldw(   phi_%02d,   phi_%02d,3);\n",$x-1,$x-4,$x);} # phi - 1
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double                 phi_%02d = vec_sldw(   phi_%02d,   phi_%02d,1);\n",$x+1,$x,$x+4);} # phi + 1
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double              beta_i_%02d = vec_sldw(beta_i_%02d,beta_i_%02d,1);\n",$x+1,$x,$x+4);} # beta_i + 1

  for($x=0;$x<=$U;$x+=4){printf(F "                vector4double        phi_jm1_temp_%02d = vec_ld( %3d-pencil8,(double*)(   phi+ijk));\n",$x,8*$x);} # round down to 32-byte boundary
  for($x=0;$x<=$U;$x+=4){printf(F "                vector4double        phi_jp1_temp_%02d = vec_ld( %3d+pencil8,(double*)(   phi+ijk));\n",$x,8*$x);} # round down to 32-byte boundary
  for($x=0;$x<=$U;$x+=4){printf(F "                vector4double     beta_j_jp1_temp_%02d = vec_ld( %3d+pencil8,(double*)(beta_j+ijk));\n",$x,8*$x);} # round down to 32-byte boundary
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double              beta_j_%02d = vec_lda(%3d        ,(double*)(beta_j+ijk));\n",$x,8*$x);} # should be aligned

                         printf(F "                vector4double          phi_jm1_permute = vec_lvsl( 0-pencil8,(double*)(   phi+ijk));\n");
                         printf(F "                vector4double          phi_jp1_permute = vec_lvsl( 0+pencil8,(double*)(   phi+ijk));\n");
                         printf(F "                vector4double       beta_j_jp1_permute = vec_lvsl( 0+pencil8,(double*)(beta_j+ijk));\n");
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double             phi_jm1_%02d = vec_perm(   phi_jm1_temp_%02d,   phi_jm1_temp_%02d,   phi_jm1_permute);\n",$x,$x,$x+4);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double             phi_jp1_%02d = vec_perm(   phi_jp1_temp_%02d,   phi_jp1_temp_%02d,   phi_jp1_permute);\n",$x,$x,$x+4);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double          beta_j_jp1_%02d = vec_perm(beta_j_jp1_temp_%02d,beta_j_jp1_temp_%02d,beta_j_jp1_permute);\n",$x,$x,$x+4);}

  for($x=0;$x< $U;$x+=4){printf(F "                vector4double              beta_k_%02d = vec_lda(%3d       ,(double*)(beta_k+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double          beta_k_kp1_%02d = vec_lda(%3d+plane8,(double*)(beta_k+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double             phi_km1_%02d = vec_lda(%3d-plane8,(double*)(   phi+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double             phi_kp1_%02d = vec_lda(%3d+plane8,(double*)(   phi+ijk));\n",$x,8*$x);}

 

  # form laplacian...
  # derivative in i-dimension... need unaligned accesses.  Therefore, round down, load extra word, permute
  # derivative in j-dimension... need unaligned accesses as pencil is not a multiple of 32 bytes.  Therefore, round down, load extra word, permute
  # derivative in k-direction... plane is a multiple of 32 bytes.  Therefore, +/-plane are alligned accesses
  for($x=0;$x<$U;$x+=4){ printf(F "                vector4double           laplacian_%02d;\n",$x);}
  for($x=0;$x<$U;$x+=4){ printf(F "                vector4double         laplacian_i_%02d;\n",$x);}
  for($x=0;$x<$U;$x+=4){ printf(F "                vector4double         laplacian_k_%02d;\n",$x);}
  for($x=0;$x<$U;$x+=4){ printf(F "                vector4double         laplacian_j_%02d;\n",$x);}

                    $x=0;printf(F "                                      laplacian_i_%02d = vec_mul(  beta_i_%02d    ,phi_m1                       );\n",$x,$x,$x,$x);       #  +beta_i[ijk       ]*phi[ijk-1     ]
  for($x=4;$x< $U;$x+=4){printf(F "                                      laplacian_i_%02d = vec_mul(  beta_i_%02d    ,phi_%02d                     );\n",$x,$x,$x-1,$x);} 
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_i_%02d = vec_nmsub(beta_i_%02d    ,phi_%02d    ,laplacian_i_%02d);\n",$x,$x,$x,$x);}      #  -beta_i[ijk       ]*phi[ijk       ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_i_%02d = vec_nmsub(beta_i_%02d    ,phi_%02d    ,laplacian_i_%02d);\n",$x,$x+1,$x,$x);}    #  -beta_i[ijk+1     ]*phi[ijk       ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_i_%02d = vec_madd( beta_i_%02d    ,phi_%02d    ,laplacian_i_%02d);\n",$x,$x+1,$x+1,$x);}  #  +beta_i[ijk+1     ]*phi[ijk+1     ]

  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_j_%02d = vec_mul(  beta_j_%02d    ,phi_jm1_%02d                 );\n",$x,$x,$x,$x);}      #  +beta_j[ijk       ]*phi[ijk-pencil]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_j_%02d = vec_nmsub(beta_j_%02d    ,phi_%02d    ,laplacian_j_%02d);\n",$x,$x,$x,$x);}      #  -beta_j[ijk       ]*phi[ijk       ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_j_%02d = vec_nmsub(beta_j_jp1_%02d,phi_%02d    ,laplacian_j_%02d);\n",$x,$x,$x,$x);}      #  -beta_j[ijk+pencil]*phi[ijk       ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_j_%02d = vec_madd( beta_j_jp1_%02d,phi_jp1_%02d,laplacian_j_%02d);\n",$x,$x,$x,$x);}      #  +beta_j[ijk+pencil]*phi[ijk+pencil]

  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_k_%02d = vec_mul(  beta_k_%02d    ,phi_km1_%02d                 );\n",$x,$x,$x,$x);}      #  +beta_k[ijk       ]*phi[ijk-plane ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_k_%02d = vec_nmsub(beta_k_%02d    ,phi_%02d    ,laplacian_k_%02d);\n",$x,$x,$x,$x);}      #  -beta_k[ijk       ]*phi[ijk       ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_k_%02d = vec_nmsub(beta_k_kp1_%02d,phi_%02d    ,laplacian_k_%02d);\n",$x,$x,$x,$x);}      #  -beta_k[ijk+plane ]*phi[ijk       ]
  for($x=0;$x< $U;$x+=4){printf(F "                                      laplacian_k_%02d = vec_madd( beta_k_kp1_%02d,phi_kp1_%02d,laplacian_k_%02d);\n",$x,$x,$x,$x);}      #  +beta_k[ijk+plane ]*phi[ijk+plane ]


  # form the helmholtz
 #for($x=0;$x< $U;$x+=4){printf(F "                                              phi_%02d = vec_lda(%3d,(double*)(   phi+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double               alpha_%02d = vec_lda(%3d,(double*)( alpha+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double         a_alpha_phi_%02d = vec_mul(vec_mul(a_splat4,phi_%02d),alpha_%02d);\n",$x,$x,$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double           helmholtz_%02d = vec_nmsub(b_h2inv_splat4,vec_add(laplacian_i_%02d,laplacian_j_%02d),vec_nmsub(b_h2inv_splat4,laplacian_k_%02d,a_alpha_phi_%02d));\n",$x,$x,$x,$x,$x,$x);}


  # GSRB relaxation part...
 #for($x=0;$x< $U;$x+=4){printf(F "                                              phi_%02d = vec_lda(%3d,(double*)(   phi+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double                 rhs_%02d = vec_lda(%3d,(double*)(   rhs+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double              lambda_%02d = vec_lda(%3d,(double*)(lambda+ijk));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double            RedBlack_%02d = vec_xor(invertMask4,vec_lda(%2d,(double*)(RedBlackMask+ij) ));\n",$x,8*$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double phi_plus_lambda_rhs_%02d = vec_madd(lambda_%02d,rhs_%02d,phi_%02d);\n",$x,$x,$x,$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                vector4double                 new_%02d = vec_nmsub(lambda_%02d,helmholtz_%02d,phi_plus_lambda_rhs_%02d);\n",$x,$x,$x,$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                                              new_%02d = vec_sel(new_%02d,phi_%02d,RedBlack_%02d);\n",$x,$x,$x,$x);}
  for($x=0;$x< $U;$x+=4){printf(F "                                                         vec_sta(new_%02d,%2d,(double*)(phi+ijk));\n",$x,8*$x);}


                         print  F "          ij+=$S;ijk+=$S;\n";
                         print  F "        }\n";

return;
}
