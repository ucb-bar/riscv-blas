	# Custom assembler routines here

	.text
	.align 10
	.globl vsetvlen
vsetvlen:
	vpset vp0
	vstop

        .align 12
        .globl gemm
	.globl sgemm_pre
	.globl sgemm_loop
        .globl sgemm_postalpha
	.globl sgemm_postbeta
	.globl sgemm_post
        .globl sgemm_preedge
	.globl sgemm_loopedge
        .globl sgemm_postalphaedge
	.globl sgemm_postbetaedge
	.globl sgemm_postedge
        .globl dgemm_pre
	.globl dgemm_loop
        .globl dgemm_postalpha
	.globl dgemm_postbeta
	.globl dgemm_post
        .globl dgemm_preedge
	.globl dgemm_loopedge
        .globl dgemm_postalphaedge
	.globl dgemm_postbetaedge
	.globl dgemm_postedge
gemm:
sgemm_pre:
        vadd vv0, vs0, vs0
        vadd vv1, vs0, vs0
        vadd vv2, vs0, vs0
        vadd vv3, vs0, vs0
        vstop
sgemm_loop:
        vlw vv4, va0
        vfmadd.s vv0, vv4, vs3, vv0
        vfmadd.s vv1, vv4, vs4, vv1
        vfmadd.s vv2, vv4, vs5, vv2
        vfmadd.s vv3, vv4, vs6, vv3
        vstop
sgemm_postalpha:
        vfmul.s vv0, vv0, vs1
        vfmul.s vv1, vv1, vs1
        vfmul.s vv2, vv2, vs1
        vfmul.s vv3, vv3, vs1
        vstop
sgemm_postbeta:
        vlw vv4, va0
        vfmadd.s vv0, vv4, vs2, vv0
        vlw vv4, va1
        vfmadd.s vv1, vv4, vs2, vv1
        vlw vv4, va2
        vfmadd.s vv2, vv4, vs2, vv2
        vlw vv4, va3
        vfmadd.s vv3, vv4, vs2, vv3
        vstop
sgemm_post:
        vsw vv0, va0
        vsw vv1, va1
        vsw vv2, va2
        vsw vv3, va3
        vstop
sgemm_preedge:
        vadd vv0, vs0, vs0
        vstop
sgemm_loopedge:
        vlw vv4, va0
        vfmadd.s vv0, vv4, vs3, vv0
        vstop
sgemm_postalphaedge:
        vfmul.s vv0, vv0, vs1
        vstop
sgemm_postbetaedge:
        vlw vv4, va0
        vfmadd.s vv0, vv4, vs2, vv0
        vstop
sgemm_postedge:
        vsw vv0, va0
        vstop
dgemm_pre:
        vadd vv0, vs0, vs0
        vadd vv1, vs0, vs0
        vadd vv2, vs0, vs0
        vadd vv3, vs0, vs0
        vstop
dgemm_loop:
        vld vv4, va0
        vfmadd.d vv0, vv4, vs3, vv0
        vfmadd.d vv1, vv4, vs4, vv1
        vfmadd.d vv2, vv4, vs5, vv2
        vfmadd.d vv3, vv4, vs6, vv3
        vstop
dgemm_postalpha:
        vfmul.d vv0, vv0, vs1
        vfmul.d vv1, vv1, vs1
        vfmul.d vv2, vv2, vs1
        vfmul.d vv3, vv3, vs1
        vstop
dgemm_postbeta:
        vld vv4, va0
        vfmadd.d vv0, vv4, vs2, vv0
        vld vv4, va1
        vfmadd.d vv1, vv4, vs2, vv1
        vld vv4, va2
        vfmadd.d vv2, vv4, vs2, vv2
        vld vv4, va3
        vfmadd.d vv3, vv4, vs2, vv3
        vstop
dgemm_post:
        vsd vv0, va0
        vsd vv1, va1
        vsd vv2, va2
        vsd vv3, va3
        vstop
dgemm_preedge:
        vadd vv0, vs0, vs0
        vstop
dgemm_loopedge:
        vld vv4, va0
        vfmadd.d vv0, vv4, vs3, vv0
        vstop
dgemm_postalphaedge:
        vfmul.d vv0, vv0, vs1
        vstop
dgemm_postbetaedge:
        vld vv4, va0
        vfmadd.d vv0, vv4, vs2, vv0
        vstop
dgemm_postedge:
        vsd vv0, va0
        vstop





        .align 12
        .globl blas1
        .globl scopy_unit
        .globl scopy_stride
        .globl sscal_unit
        .globl sscal_stride
        .globl dcopy_unit
        .globl dcopy_stride
        .globl dscal_unit
        .globl dscal_stride
        .globl saxpy_loop
        .globl saxpy_stride_loop
        .globl daxpy_loop
        .globl daxpy_stride_loop
        .globl sdot_pre
        .globl sdot_loop
        .globl sdot_post
        .globl sdot_reduce_loop
        .globl sdot_stride_loop
blas1:
scopy_unit:
        vlw vv0, va0
        vsw vv0, va1
        vstop
scopy_stride:
        vlstw vv0, va0, va2
        vsstw vv0, va1, va3
        vstop
sscal_unit:
        vlw vv0, va0
        vfmul.s vv0, vv0, vs1
        vsw vv0, va0
        vstop
sscal_stride:
        vlstw vv0, va0, va1
        vfmul.s vv0, vv0, vs1
        vsstw vv0, va0, va1
        vstop
dcopy_unit:
        vld vv0, va0
        vsd vv0, va1
        vstop
dcopy_stride:
        vlstd vv0, va0, va2
        vsstd vv0, va1, va3
        vstop
dscal_unit:
        vld vv0, va0
        vfmul.d vv0, vv0, vs1
        vsd vv0, va0
        vstop
dscal_stride:
        vlstd vv0, va0, va1
        vfmul.d vv0, vv0, vs1
        vsstd vv0, va0, va1
        vstop
saxpy_loop:
        vlw vv0, va0
        vlw vv1, va1
        vfmadd.s vv1, vs1, vv0, vv1
        vsw vv1, va1
        vstop
saxpy_stride_loop:
        vlstw vv0, va0, va2
        vlstw vv1, va1, va3
        vfmadd.s vv1, vs1, vv0, vv1
        vsstw vv1, va1, va3
        vstop
daxpy_loop:
        vld vv0, va0
        vld vv1, va1
        vfmadd.d vv1, vs1, vv0, vv1
        vsd vv1, va1
        vstop  
daxpy_stride_loop:
        vlstd vv0, va0, va2
        vlstd vv1, va1, va3
        vfmadd.d vv1, vs1, vv0, vv1
        vsstd vv1, va1, va3
        vstop
sdot_pre:
        vfsub.s vv2, vv2, vv2
        vstop
sdot_loop:
        vlw vv0, va0
        vlw vv1, va1
        vfmadd.s vv2, vv0, vv1, vv2
        vstop
sdot_post:
        vsw vv2, va2
        vstop 
sdot_reduce_loop:
        vlw vv1, va1
        vfadd.s vv2, vv2, vv1
        vsw vv2, va2
        vstop
sdot_stride_loop:
        vlstw vv0, va0, va3
        vlstw vv1, va1, va4
        vfmadd.s vv2, vv0, vv1, vv2
        vstop





        .align 12
        .globl blas2
        .globl sgemv_zero_loop
        .globl sgemv_beta_loop
        .globl sgemv_stride_zero_loop
        .globl sgemv_stride_beta_loop
blas2:
sgemv_zero_loop:
        vlw vv0, va0
        vfsub.s vv0, vv0, vv0
        vsw vv0, va0
        vstop
sgemv_beta_loop:
        vlw vv0, va0
        vfmul.s vv0, vv0, vs1
        vsw vv0, va0
        vstop
sgemv_stride_zero_loop:
        vlstw vv0, va0, va1
        vfsub.s vv0, vv0, vv0
        vsstw vv0, va0, va1
        vstop
sgemv_stride_beta_loop:
        vlstw vv0, va0, va1
        vfmul.s vv0, vv0, vs1
        vsstw vv0, va0, va1
        vstop
