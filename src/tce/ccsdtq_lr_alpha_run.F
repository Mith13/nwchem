      subroutine ccsdtq_lr_alpha_run(rtdb,d_d0,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_t3,d_t4,
     1           d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     1           d_tr1,d_tr2,d_tr3,d_tr4,
     2           k_d0_offset,k_a0_offset,k_f1_offset,k_v2_offset,
     3           k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5           k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset,
     6           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     7           size_tr1,size_tr2,size_tr3,size_tr4,size_a0,l_pole)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer rtdb
      integer i,j,dummy,axis
      integer omegacount,omegasign,dynaxis,dynfreq
      integer irrep_g
      parameter (irrep_g=0)
      integer d_d0,d_a0,d_f1,d_v2,d_d1(10),size_a0
      integer d_t1,d_t2,d_t3,d_t4
      integer d_lambda1,d_lambda2,d_lambda3,d_lambda4
      integer d_tr1(38),d_tr2(38),d_tr3(38),d_tr4(38)
      integer k_d0_offset,k_a0_offset
      integer k_f1_offset,k_v2_offset,k_d1_offset(10)
      integer k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset
      integer k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset
      integer k_tr1_offset(19),k_tr2_offset(19),k_tr3_offset(19),
     1        k_tr4_offset(19)
      integer size_tr1(19),size_tr2(19),size_tr3(19),size_tr4(19)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision cpu, wall
      double precision omega
      double precision debye
      parameter (debye=2.54176568d0)
      logical nodezero
      character*4 irrepname
      character*3 axisname(19)  ! Axis
      data axisname/'X','Y','Z','XX','XY','XZ','YY','YZ','ZZ',
     1              'XXX','XXY','XXZ','XYY','XYZ',
     2              'XZZ','YYY','YYZ','YZZ','ZZZ'/
      character*255 filename

      integer l_pole,axis_start,axis_max
      if (l_pole.eq.1) then
        axis_start = 1
        axis_max = 3
      else if (l_pole.eq.2) then
        axis_start = 4
        axis_max = 9
      else if (l_pole.eq.3) then
        axis_start = 10
        axis_max = 19
      else
        call errquit('ccsdtq_lr_alpha: Wrong multipole',CALC_ERR)
      endif

c
      nodezero=(ga_nodeid().eq.0)
c
!          if (nodezero) write(LuOut,*) "lineresp is go"
      do omegacount=1,anumfreq
!            if (nodezero) write(LuOut,*) "omegacount is go"
        omega = afreq(omegacount)
        if (omega.eq.(0.0d0)) then
          dynfreq = 1
        else
          dynfreq = 2
        endif
!
        do axis = axis_start, axis_max
!            if (nodezero) write(LuOut,*) "axis is go"
         if (respaxis(axis)) then
!            if (nodezero) write(LuOut,*) "respaxis(axis) is go"
          irrep_d=sym_abelian_axis(geom,axis)
          call sym_irrepname(geom,irrep_d+1,irrepname)
          if (nodezero.and.util_print('mod1',print_default)) then
            write(LuOut,*)
            write(LuOut,9440) axisname(axis),irrepname
          endif
          irrep_o=irrep_d
          irrep_x=irrep_d
!
          do omegasign=1,dynfreq
!            if (nodezero) write(LuOut,*) "omegasign is go"
            if (omegasign.eq.1) then
              omega = (1.0d0)*afreq(omegacount)
              dynaxis = 0
              if (l_pole.eq.2) then
                  dynaxis = 3 
              else if (l_pole.eq.3) then 
                  dynaxis = 9 
              endif
            elseif (omegasign.eq.2) then
              omega = (-1.0d0)*afreq(omegacount)
              dynaxis = 3
              if (l_pole.eq.2) then
                 dynaxis = 9
              else if (l_pole.eq.3) then
                 dynaxis = 19
              endif
            endif
            if (nodezero) write(LuOut,9431) omega
!            if (nodezero) write(LuOut,*) "omega = ",omega
!            if (nodezero) write(LuOut,*) "dynaxis = ",dynaxis
!
!
            call ccsdtq_lr_iter(axis,dynaxis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_t3,d_t4,
     1           d_tr1,d_tr2,d_tr3,d_tr4,
     2           k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     4           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     5           size_tr1,size_tr2,size_tr3,size_tr4,axis_start-1)
          enddo ! omegasign loop
         endif ! respaxis(axis)
       enddo ! axis loop
       if (nopolar) then
         if (nodezero) write(LuOut,9701)
       else if (l_pole.eq.1) then
         call  ccsdtq_lr_alpha_eval(rtdb,omega,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_t3,d_t4,
     1           d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     2           d_tr1,d_tr2,d_tr3,d_tr4,
     3           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5           k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset,
     6           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     7           size_a0)
       else 
         call  ccsdtq_lr_alpha_eval_l(rtdb,omega,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_t3,d_t4,
     1           d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     2           d_tr1,d_tr2,d_tr3,d_tr4,
     3           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5           k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset,
     6           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     7           size_a0,l_pole)
       endif
      enddo !omegacount
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
 9701 format(1x,'Skipping evaluation of polarizability by request')
      return
      end

      subroutine ccsdtq_lr_iter(axis,dynaxis,omega,
     1           d_f1,d_v2,d_d1,d_t1,d_t2,d_t3,d_t4,
     1           d_tr1,d_tr2,d_tr3,d_tr4,
     2           k_f1_offset,k_v2_offset,k_d1_offset,
     3           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     4           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     5           size_tr1,size_tr2,size_tr3,size_tr4,axis_start)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
c#include "tce_diis.fh"
#include "tce_prop.fh"
#include "tce_restart.fh"
c
      integer i,j,dummy,axis,iter
      integer omegacount,omegasign,dynaxis
      integer irrep_g
      parameter (irrep_g=0)
      integer d_f1,d_v2,d_d1(10),d_t1,d_t2,d_t3,d_t4
      integer d_tr1(38),d_tr2(38),d_tr3(38),d_tr4(38)
      integer d_rr1(19),d_rr2(19),d_rr3(19),d_rr4(19)
      integer k_f1_offset,k_v2_offset,k_d1_offset(10)
      integer k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset
      integer k_tr1_offset(19),k_tr2_offset(19),k_tr3_offset(19),
     1        k_tr4_offset(19)
      integer size_tr1(19),size_tr2(19),size_tr3(19),size_tr4(19)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision residual,rr1,rr2,rr3,rr4
      double precision ddotfile
      external ddotfile
      character*255 filename
      character*6 rr1filename(19) ! File name stub
      data rr1filename/'rr1x  ','rr1y  ','rr1z  ',
     1                 'rr1xx ','rr1xy ','rr1xz ',
     2                 'rr1yy ','rr1yz ','rr1zz ',
     3                 'rr1xxx','rr1xxy','rr1xxz','rr1xyy','rr1xyz',
     4                 'rr1xzz','rr1yyy','rr1yyz','rr1yzz','rr1zzz'/

      character*6 rr2filename(19) ! File name stub
      data rr2filename/'rr2x  ','rr2y  ','rr2z  ',
     1                 'rr2xx ','rr2xy ','rr2xz ',
     2                 'rr2yy ','rr2yz ','rr2zz ',
     3                 'rr2xxx','rr2xxy','rr2xxz','rr2xyy','rr2xyz',
     4                 'rr2xzz','rr2yyy','rr2yyz','rr2yzz','rr2zzz'/
      character*6 rr3filename(19) ! File name stub
      data rr3filename/'rr3x  ','rr3y  ','rr3z  ',
     1                 'rr3xx ','rr3xy ','rr3xz ',
     2                 'rr3yy ','rr3yz ','rr3zz ',
     3                 'rr3xxx','rr3xxy','rr3xxz','rr3xyy','rr3xyz',
     4                 'rr3xzz','rr3yyy','rr3yyz','rr3yzz','rr3zzz'/
      character*6 rr4filename(19) ! File name stub
      data rr4filename/'rr4x  ','rr4y  ','rr4z  ',
     1                 'rr4xx ','rr4xy ','rr4xz ',
     2                 'rr4yy ','rr4yz ','rr4zz ',
     3                 'rr4xxx','rr4xxy','rr4xxz','rr4xyy','rr4xyz',
     4                 'rr4xzz','rr4yyy','rr4yyz','rr4yzz','rr4zzz'/
      integer axis_start
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c

        call tce_diis_init()
        do iter=1,maxiter
          cpu=-util_cpusec()
          wall=-util_wallsec()
          if (nodezero.and.(iter.eq.1))
     &          write(LuOut,9400) "CCSDTQ-LR"
!
          call tce_filename(rr1filename(axis),filename)
          call createfile(filename,d_rr1(axis),size_tr1(axis))
          call tce_zero(d_rr1(axis),size_tr1(axis))
!        write(LuOut,*) "daxpyfile"
          call daxpyfile(1,(-1.0d0)*omega,d_tr1(axis+dynaxis),
     1             d_rr1(axis),size_tr1(axis))
!        write(LuOut,*) "eomccsdtq_x1"
          call eomccsdtq_x1(d_f1,d_rr1(axis),d_t1,d_t2,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),k_f1_offset,
     3             k_tr1_offset(axis),
     4             k_t1_offset,k_t2_offset,k_v2_offset,
     5             k_tr1_offset(axis),k_tr2_offset(axis),
     6             k_tr3_offset(axis))
!        write(LuOut,*) "ccsdtq1_o1"
          call ccsdtq_o1(d_rr1(axis),d_d1(axis-axis_start),
     &             d_t1,d_t2,
     &             k_tr1_offset(axis),k_d1_offset(axis-axis_start),
     &             k_t1_offset,k_t2_offset)
!
          call tce_filename(rr2filename(axis),filename)
          call createfile(filename,d_rr2(axis),size_tr2(axis))
          call tce_zero(d_rr2(axis),size_tr2(axis))
!        write(LuOut,*) "daxpyfile"
          call daxpyfile(1,(-1.0d0)*omega,d_tr2(axis+dynaxis),
     1             d_rr2(axis),size_tr2(axis))
!        write(LuOut,*) "eomccsdtq_x2"
          call eomccsdtq_x2(d_f1,d_rr2(axis),d_t1,d_t2,d_t3,d_v2,
     1             d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr2_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_v2_offset,
     5             k_tr1_offset(axis),k_tr2_offset(axis),
     6             k_tr3_offset(axis),k_tr4_offset(axis))
!        write(LuOut,*) "ccsdtq_o2"
          call ccsdtq_o2(d_rr2(axis),d_d1(axis-axis_start),
     &             d_t1,d_t2,d_t3,
     &             k_tr2_offset(axis),k_d1_offset(axis-axis_start),
     &             k_t1_offset,k_t2_offset,k_t3_offset)
!
          call tce_filename(rr3filename(axis),filename)
          call createfile(filename,d_rr3(axis),size_tr3(axis))
          call tce_zero(d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "daxpyfile"
          call daxpyfile(1,(-1.0d0)*omega,d_tr3(axis+dynaxis),
     1             d_rr3(axis),size_tr3(axis))
!        write(LuOut,*) "eomccsdtq_x3"
          call eomccsdtq_x3(d_f1,d_rr3(axis),d_t1,d_t2,d_t3,d_t4,
     1             d_v2,d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr3_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5             k_v2_offset,
     6             k_tr1_offset(axis),k_tr2_offset(axis),
     7             k_tr3_offset(axis),k_tr4_offset(axis))
!        write(LuOut,*) "ccsdtq_o3"
          call ccsdtq_o3(d_rr3(axis),d_d1(axis-axis_start),
     &             d_t1,d_t2,d_t3,d_t4,
     &             k_tr3_offset(axis),k_d1_offset(axis-axis_start),
     &             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset)
!
          call tce_filename(rr4filename(axis),filename)
          call createfile(filename,d_rr4(axis),size_tr4(axis))
          call tce_zero(d_rr4(axis),size_tr4(axis))
!        write(LuOut,*) "daxpyfile"
          call daxpyfile(1,(-1.0d0)*omega,d_tr4(axis+dynaxis),
     1             d_rr4(axis),size_tr4(axis))
!        write(LuOut,*) "eomccsdtq_x4"
          call eomccsdtq_x4(d_f1,d_rr4(axis),d_t1,d_t2,d_t3,d_t4,
     1             d_v2,d_tr1(axis+dynaxis),d_tr2(axis+dynaxis),
     2             d_tr3(axis+dynaxis),d_tr4(axis+dynaxis),
     3             k_f1_offset,k_tr4_offset(axis),
     4             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5             k_v2_offset,
     6             k_tr1_offset(axis),k_tr2_offset(axis),
     7             k_tr3_offset(axis),k_tr4_offset(axis))
!        write(LuOut,*) "ccsdtq_o4"
          call ccsdtq_o4(d_rr4(axis),d_d1(axis-axis_start),
     &             d_t1,d_t2,d_t3,d_t4,
     &             k_tr4_offset(axis),k_d1_offset(axis-axis_start),
     &             k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset)
!
!        write(LuOut,*) "reconcilefile"
          call reconcilefile(d_rr1(axis),size_tr1(axis))
          call reconcilefile(d_rr2(axis),size_tr2(axis))
          call reconcilefile(d_rr3(axis),size_tr3(axis))
          call reconcilefile(d_rr4(axis),size_tr4(axis))
!        write(LuOut,*) "tce_residual_tr1"
          call tce_residual_tr1(d_rr1(axis),k_tr1_offset(axis),rr1)
!        write(LuOut,*) "tce_residual_tr2"
          call tce_residual_tr2(d_rr2(axis),k_tr2_offset(axis),rr2)
!        write(LuOut,*) "tce_residual_tr3"
          call tce_residual_tr3(d_rr3(axis),k_tr3_offset(axis),rr3)
!        write(LuOut,*) "tce_residual_tr4"
          call tce_residual_tr4(d_rr4(axis),k_tr4_offset(axis),rr4)
          residual = max(rr1,rr2,rr3,rr4)
          cpu=cpu+util_cpusec()
          wall=wall+util_wallsec()
          if (nodezero) write(LuOut,9420) iter,residual,cpu,wall
          if (residual .lt. thresh) then
            if (nodezero) then
              write(LuOut,9410)
            endif
            if (ampnorms) then
              call tce_residual_tr1(d_tr1(axis+dynaxis),
     1                                  k_tr1_offset(axis),rr1)
              call tce_residual_tr2(d_tr2(axis+dynaxis),
     1                                  k_tr2_offset(axis),rr2)
              call tce_residual_tr3(d_tr3(axis+dynaxis),
     1                                  k_tr3_offset(axis),rr3)
              call tce_residual_tr4(d_tr4(axis+dynaxis),
     1                                  k_tr4_offset(axis),rr4)
              if (nodezero) then
                write(LuOut,9082) "T(1) singles",rr1
                write(LuOut,9082) "T(1) doubles",rr2
                write(LuOut,9082) "T(1) triples",rr3
                write(LuOut,9082) "T(1) quadruples",rr4
              endif
            endif

            call deletefile(d_rr4(axis))
            call deletefile(d_rr3(axis))
            call deletefile(d_rr2(axis))
            call deletefile(d_rr1(axis))
            call tce_diis_tidy()
            if (save_tr(1)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T1(1) now...'
              endif
              call tr1_restart_save(d_tr1(axis+dynaxis),
     1                 k_tr1_offset(axis),size_tr1(axis),
     2                 axis+dynaxis,handle_tr1(axis),irrep_x)
            endif
            if (save_tr(2)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T2(1) now...'
              endif
              call tr2_restart_save(d_tr2(axis+dynaxis),
     1                 k_tr2_offset(axis),size_tr2(axis),
     2                 axis+dynaxis,handle_tr2(axis),irrep_x)
            endif
            if (save_tr(3)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T3(1) now...'
              endif
              call tr3_restart_save(d_tr3(axis),k_tr3_offset(axis),
     1                 size_tr3(axis),axis,handle_tr3(axis),irrep_x)
            endif
            if (save_tr(4)) then
              if(nodezero) then
                write(LuOut,*) 'Saving T4(1) now...'
              endif
              call tr4_restart_save(d_tr4(axis),k_tr4_offset(axis),
     1                 size_tr4(axis),axis,handle_tr4(axis),irrep_x)
            endif
            return 
          endif
          if (save_tr(1).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T1(1) now...'
            endif
            call tr1_restart_save(d_tr1(axis+dynaxis),
     1               k_tr1_offset(axis),size_tr1(axis),
     2               axis+dynaxis,handle_tr1(axis),irrep_x)
          endif
          if (save_tr(2).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T2(1) now...'
            endif
            call tr2_restart_save(d_tr2(axis+dynaxis),
     1               k_tr2_offset(axis),size_tr2(axis),
     2               axis+dynaxis,handle_tr2(axis),irrep_x)
          endif
          if (save_tr(3).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T3(1) now...'
            endif
            call tr3_restart_save(d_tr3(axis),k_tr3_offset(axis),
     1               size_tr3(axis),axis,handle_tr3(axis),irrep_x)
          endif
          if (save_tr(4).and.(mod(iter,save_interval).eq.0)) then
            if(nodezero) then
              write(LuOut,*) 'Saving T4(1) now...'
            endif
            call tr4_restart_save(d_tr4(axis),k_tr4_offset(axis),
     1               size_tr4(axis),axis,handle_tr4(axis),irrep_x)
          endif
          call tce_diis2(.false.,iter,.true.,.true.,.true.,.true.,
     1             d_rr1(axis),d_tr1(axis+dynaxis),k_tr1_offset(axis),
     2             size_tr1(axis),
     3             d_rr2(axis),d_tr2(axis+dynaxis),k_tr2_offset(axis),
     4             size_tr2(axis),
     5             d_rr3(axis),d_tr3(axis+dynaxis),k_tr3_offset(axis),
     6             size_tr3(axis),
     7             d_rr4(axis),d_tr4(axis+dynaxis),k_tr4_offset(axis),
     8             size_tr4(axis),omega)
          call deletefile(d_rr4(axis))
          call deletefile(d_rr3(axis))
          call deletefile(d_rr2(axis))
          call deletefile(d_rr1(axis))
          if (nodezero) call util_flush(LuOut)
        enddo ! iter loop
      call errquit('tce_energy: maxiter exceeded',iter,CALC_ERR)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9082 format(1x,'amplitude norm of ',A9,' = ',f25.15)
 9100 format(1x,i4,2f18.13,2f8.1)
 9120 format(1x,A)
 9121 format(/,1x,A)
 9122 format(1x,A,i4)
 9400 format(/,1x,A,' iterations',/,
     1  1x,'---------------------------------------------',/
     2  1x,'Iter          Residuum            Cpu    Wall',/
     3  1x,'---------------------------------------------')
 9410 format(
     1  1x,'---------------------------------------------',/
     2  1x,'Iterations converged')
 9420 format(1x,i4,f25.13,2f8.1)
 9431 format(/,1x,'Frequency = ',f15.7,' / au')
 9440 format(1x,A3,' axis ( ',A4,'symmetry)')
      return
      end

      subroutine ccsdtq_lr_alpha_eval(rtdb,omega,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_t3,d_t4,
     1           d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     2           d_tr1,d_tr2,d_tr3,d_tr4,
     3           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5           k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset,
     6           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     7           size_a0)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_prop.fh"
c
      integer rtdb
      integer i,j,dummy,axis
      integer dynfreq,dynaxis,axisA,axisB,omegasign,omegacount
      integer d_a0,d_f1,d_v2,d_d1(3),size_a0
      integer d_t1,d_t2,d_t3,d_t4
      integer d_lambda1,d_lambda2,d_lambda3,d_lambda4
      integer d_tr1(6),d_tr2(6),d_tr3(6),d_tr4(6)
      integer k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset(3)
      integer k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset
      integer k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset
      integer k_tr1_offset(3),k_tr2_offset(3),k_tr3_offset(3),
     1        k_tr4_offset(3)
      integer irrep_g
      parameter (irrep_g=0)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision alpha1,alpha2
      double precision alpha(3,3)       ! Dipole polarizability tensor
      double precision alphacopy(3,3)   ! Dipole polarizability tensor copy
      double precision alphaiso         ! Isotropic dipole polarizability
      double precision alphaani         ! Anisotropic dipole polarizability
      double precision alphaevr(3)      ! Dipole polarizability tensor eigenvalues (real)
      double precision alphaevi(3)      ! Dipole polarizability tensor eigenvalues (imag)
      double precision aeigvec(3,3)     ! dummy dgeev eigenvectors
      double precision auXnm            ! Conversion factor from a.u. (frequency) to nm (wavelength)
      parameter (auXnm=45.563353d0)
      double precision au2ang           ! Conversion factor from bohr to Angstrom
      parameter (au2ang=5.29177249d-01)
      double precision au2ang3          ! Conversion factor from bohr^3 to Angstrom^3
      parameter (au2ang3=au2ang*au2ang*au2ang)
      double precision tol_omega
      parameter (tol_omega=1.0d-50)
      double precision ddotfile
      external ddotfile
      character*3 axisname(3)  ! Axis
      data axisname/'X','Y','Z'/
      character*255 filename
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
       
      do axisA = 1, 3
      do axisB = 1, axisA
      alpha(axisA,axisB)=0.0d0
      if (respaxis(axisA).and.respaxis(axisB)) then
        irrep_a=sym_abelian_axis(geom,axisA)
        irrep_b=sym_abelian_axis(geom,axisB)
        irrep_y=irrep_g
!
!            if (nodezero) write(LuOut,*) "axisA = ",axisA
!            if (nodezero) write(LuOut,*) "axisB = ",axisB
!
        call tce_filename('a0',filename)
        call createfile(filename,d_a0,size_a0)
!
          alpha1=0.0d0
          alpha2=0.0d0
!
          if (omega.ne.0.0d0) then
            dynaxis = 3
          else
            dynaxis = 0
          endif
!
          irrep_c=irrep_a
          irrep_oa=irrep_a
          irrep_ob=irrep_b
          irrep_tra=irrep_a
          irrep_trb=irrep_b
!              if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha1"
          call ccsdtq_lr_alpha(d_f1,d_a0,d_d1(axisA),d_d1(axisB),
     1         d_t1,d_t2,d_t3,d_t4,
     2         d_tr1(axisA+dynaxis),d_tr2(axisA+dynaxis),
     3         d_tr3(axisA+dynaxis),d_tr4(axisA+dynaxis),
     4         d_tr1(axisB),d_tr2(axisB),d_tr3(axisB),d_tr4(axisB),
     5         d_v2,d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     6         k_f1_offset,k_a0_offset,k_d1_offset(axisA),
     7         k_d1_offset(axisB),k_t1_offset,k_t2_offset,
     8         k_t3_offset,k_t4_offset,k_tr1_offset(axisA),
     9         k_tr2_offset(axisA),k_tr3_offset(axisA),
     1         k_tr4_offset(axisA),k_tr1_offset(axisB),
     2         k_tr2_offset(axisB),k_tr3_offset(axisB),
     3         k_tr4_offset(axisB),k_v2_offset,k_l1_offset,
     4         k_l2_offset,k_l3_offset,k_l4_offset)
!
          call reconcilefile(d_a0,size_a0)
          call get_block(d_a0,alpha1,1,0)
          call tce_zero(d_a0,size_a0)
!              write(LuOut,*) "alpha1 = ",alpha1
!
          if (omega.eq.0.0d0) then
            alpha2=alpha1
          else
!                if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha2"
            call ccsdtq_lr_alpha(d_f1,d_a0,d_d1(axisA),d_d1(axisB),
     1         d_t1,d_t2,d_t3,d_t4,
     2         d_tr1(axisA),d_tr2(axisA),d_tr3(axisA),d_tr4(axisA),
     3         d_tr1(axisB+dynaxis),d_tr2(axisB+dynaxis),
     4         d_tr3(axisB+dynaxis),d_tr4(axisB+dynaxis),
     5         d_v2,d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     6         k_f1_offset,k_a0_offset,k_d1_offset(axisA),
     7         k_d1_offset(axisB),k_t1_offset,k_t2_offset,
     8         k_t3_offset,k_t4_offset,k_tr1_offset(axisA),
     9         k_tr2_offset(axisA),k_tr3_offset(axisA),
     1         k_tr4_offset(axisA),k_tr1_offset(axisB),
     2         k_tr2_offset(axisB),k_tr3_offset(axisB),
     3         k_tr4_offset(axisB),k_v2_offset,k_l1_offset,
     4         k_l2_offset,k_l3_offset,k_l4_offset)
!
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha2,1,0)
            call tce_zero(d_a0,size_a0)
          endif ! axisA.eq.axisB
!              write(LuOut,*) "alpha2 = ",alpha2
!
          alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2)
!
        call deletefile(d_a0)
!
!            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
!            if (nodezero) write(LuOut,9020) cpu, wall
!
      endif ! respaxis(axis)
      enddo ! axisB loop
      enddo ! axisA loop
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
!
      do i = 1, 3
        do j = 1, i
          alphacopy(i,j)=alpha(i,j)
          alphacopy(j,i)=alpha(i,j)
        enddo
      enddo
!
      call hnd_diag(alphacopy,alphaevr,3,.false.,.false.)
!
      alphaiso = (alphaevr(1)+alphaevr(2)+alphaevr(3))/3.0d0
      alphaani = (alphaevr(1)-alphaevr(2))*(alphaevr(1)-alphaevr(2))
     1             + (alphaevr(1)-alphaevr(3))*(alphaevr(1)-alphaevr(3))
     2             + (alphaevr(2)-alphaevr(3))*(alphaevr(2)-alphaevr(3))
      alphaani = dsqrt(0.5d0*alphaani)
!
      if ((.not.(respaxis(1).and.respaxis(2).and.respaxis(3)))
     1        .and.nodezero) write(LuOut,9911)
!
      if (nodezero) then
         if (dabs(omega).gt.tol_omega) then
          write(LuOut,9434) "CCSDTQ Linear Response",
     1      afreq(omegacount),auXnm/afreq(omegacount),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,
     2      alphaani,au2ang3*alphaani
         else
          write(LuOut,9436) "CCSDTQ Linear Response",
     1      dabs(omega),
     2      alpha(1,1),alpha(2,1),alpha(3,1),
     3      au2ang3*alpha(1,1),au2ang3*alpha(2,1),au2ang3*alpha(3,1),
     4      alpha(2,1),alpha(2,2),alpha(3,2),
     5      au2ang3*alpha(2,1),au2ang3*alpha(2,2),au2ang3*alpha(3,2),
     6      alpha(3,1),alpha(3,2),alpha(3,3),
     7      au2ang3*alpha(3,1),au2ang3*alpha(3,2),au2ang3*alpha(3,3),
     8      alphaevr(1),alphaevr(2),alphaevr(3),
     9      au2ang3*alphaevr(1),au2ang3*alphaevr(2),au2ang3*alphaevr(3),
     1      alphaiso,au2ang3*alphaiso,
     2      alphaani,au2ang3*alphaani
         endif
      endif
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)
 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9434 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength = ',f15.7,' / nm',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     2  1x,'                    atomic units (bohr^3)       '
     2    ,'       |                   angstroms^3           ',/
     2  1x,'                 X              Y              Z',
     2  1x,'      |             X              Y              Z',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     4  1x,'X      ',3f15.7,3x,'|',3f15.7,/
     5  1x,'Y      ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Z      ',3f15.7,3x,'|',3f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     6  1x,'Eigs = ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Isotropic   = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     6  1x,'Anisotropic = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------')
 9436 format(/,1x,A,' polarizability / au ',/
     1  1x,'Frequency  = ',f15.7,' / au',/
     1  1x,'Wavelength =        Infinity / nm',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     2  1x,'                    atomic units (bohr^3)       '
     2    ,'       |                   angstroms^3           ',/
     2  1x,'                 X              Y              Z',
     2  1x,'      |             X              Y              Z',/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     4  1x,'X      ',3f15.7,3x,'|',3f15.7,/
     5  1x,'Y      ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Z      ',3f15.7,3x,'|',3f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------',/
     6  1x,'Eigs = ',3f15.7,3x,'|',3f15.7,/
     6  1x,'Isotropic   = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     6  1x,'Anisotropic = ',8x,1f15.7,3x,15x,'|',15x,1f15.7,/
     3  1x,'-----------------------------------------------'
     3    ,'--------|-----------------------------------------------')
 9911 format(/,1x,'---> Warning: you have not solved ',
     1            'the response equations for all axes.  ',
     2            'Please analyze the results carefully as ',
     3            'the average and anisotropic polarizabilities ',
     4            'are surely wrong.',/)
      return
      end
      
      subroutine ccsdtq_lr_alpha_eval_l(rtdb,omega,d_a0,d_f1,d_v2,d_d1,
     1           d_t1,d_t2,d_t3,d_t4,
     1           d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     2           d_tr1,d_tr2,d_tr3,d_tr4,
     3           k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset,
     4           k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset,
     5           k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset,
     6           k_tr1_offset,k_tr2_offset,k_tr3_offset,k_tr4_offset,
     7           size_a0,l_pole)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_prop.fh"
c
      integer rtdb
      integer i,j,dummy,axis
      integer dynfreq,dynaxis,axisA,axisB,omegasign,omegacount
      integer d_a0,d_f1,d_v2,d_d1(10),size_a0
      integer d_t1,d_t2,d_t3,d_t4
      integer d_lambda1,d_lambda2,d_lambda3,d_lambda4
      integer d_tr1(38),d_tr2(38),d_tr3(38),d_tr4(38)
      integer k_a0_offset,k_f1_offset,k_v2_offset,k_d1_offset(10)
      integer k_t1_offset,k_t2_offset,k_t3_offset,k_t4_offset
      integer k_l1_offset,k_l2_offset,k_l3_offset,k_l4_offset
      integer k_tr1_offset(19),k_tr2_offset(19),k_tr3_offset(19),
     1        k_tr4_offset(19)
      integer irrep_g
      parameter (irrep_g=0)
      integer sym_abelian_axis
      external sym_abelian_axis
      double precision omega
      double precision cpu, wall
      double precision alpha(19,19), alpha1,alpha2
      double precision auXnm            ! Conversion factor from a.u. (frequency) to nm (wavelength)
      parameter (auXnm=45.563353d0)
      double precision au2ang           ! Conversion factor from bohr to Angstrom
      parameter (au2ang=5.29177249d-01)
      double precision au2ang3          ! Conversion factor from bohr^3 to Angstrom^3
      parameter (au2ang3=au2ang*au2ang*au2ang)
      double precision tol_omega
      parameter (tol_omega=1.0d-50)
      double precision ddotfile
      external ddotfile
      character*3 axisname(19)  ! Axis
      data axisname/'X','Y','Z','XX','XY','XZ','YY','YZ','ZZ',
     1              'XXX','XXY','XXZ','XYY','XYZ',
     2              'XZZ','YYY','YYZ','YZZ','ZZZ'/

      character*255 filename
      logical nodezero
      integer l_pole,axis_max,axis_start,pos_offset

      
      pos_offset = 0
      if (l_pole.eq.2) then
        pos_offset = 3
        axis_start = 4
        axis_max = 9
      else if (l_pole.eq.3) then
        pos_offset = 9
        axis_start = 10
        axis_max = 19
      else
        call errquit('ccsdtq_lr_alpha: Wrong multipole',l_pole,CALC_ERR)
      endif

c
      nodezero=(ga_nodeid().eq.0)
c

      do axisA = axis_start,axis_max
        do axisB = axis_start, axisA
        alpha(axisA,axisB)=0.0d0
         if (respaxis(axisA).and.respaxis(axisB)) then
          irrep_a=sym_abelian_axis(geom,axisA)
          irrep_b=sym_abelian_axis(geom,axisB)
          irrep_y=irrep_g
!
!              if (nodezero) write(LuOut,*) "axisA = ",axisA
!              if (nodezero) write(LuOut,*) "axisB = ",axisB
!
          call tce_filename('a0',filename)
          call createfile(filename,d_a0,size_a0)
!
            alpha1=0.0d0
            alpha2=0.0d0
!
          if (omega.ne.0.0d0) then
            if (l_pole.eq.3) then
                dynaxis = 19
            else if (l_pole.eq.2) then
                dynaxis = 9
            else
              call errquit('ccsdtq_lr_alpha: Wrong mom',l_pole,CALC_ERR)
            endif
          else
            dynaxis = 0
            if (l_pole.eq.2) then
                dynaxis = 3 
            else if (l_pole.eq.3) then 
                dynaxis = 9 
            else 
                call errquit('ccsdtq_lr_alpha: Wrong mom',
     1                 l_pole,CALC_ERR)
            endif
          endif
!
          irrep_c=irrep_a
          irrep_oa=irrep_a
          irrep_ob=irrep_b
          irrep_tra=irrep_a
          irrep_trb=irrep_b
!              if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
          call ccsdtq_lr_alpha(d_f1,d_a0,
     1         d_d1(axisA-(axis_start-1)),d_d1(axisB-(axis_start-1)),
     1         d_t1,d_t2,d_t3,d_t4,
     2         d_tr1(axisA+dynaxis),d_tr2(axisA+dynaxis),
     3         d_tr3(axisA+dynaxis),d_tr4(axisA+dynaxis),
     4         d_tr1(axisB+pos_offset),d_tr2(axisB+pos_offset),
     4         d_tr3(axisB+pos_offset),d_tr4(axisB+pos_offset),
     5         d_v2,d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     6         k_f1_offset,k_a0_offset,
     6         k_d1_offset(axisA-(axis_start-1)),
     7         k_d1_offset(axisB-(axis_start-1)),
     7         k_t1_offset,k_t2_offset,
     8         k_t3_offset,k_t4_offset,k_tr1_offset(axisA),
     9         k_tr2_offset(axisA),k_tr3_offset(axisA),
     1         k_tr4_offset(axisA),k_tr1_offset(axisB),
     2         k_tr2_offset(axisB),k_tr3_offset(axisB),
     3         k_tr4_offset(axisB),k_v2_offset,k_l1_offset,
     4         k_l2_offset,k_l3_offset,k_l4_offset)
!
          call reconcilefile(d_a0,size_a0)
          call get_block(d_a0,alpha1,1,0)
          call tce_zero(d_a0,size_a0)
!              write(LuOut,*) "alpha1 = ",alpha1
!
          if (omega.eq.0.0d0) then
            alpha2=alpha1
          else
!                if (nodezero) write(LuOut,*) "ccsdtq_lr_alpha"
            call ccsdtq_lr_alpha(d_f1,d_a0,
     1         d_d1(axisA-(axis_start-1)),d_d1(axisB-(axis_start-1)),
     1         d_t1,d_t2,d_t3,d_t4,
     2         d_tr1(axisA+pos_offset),d_tr2(axisA+pos_offset),
     2         d_tr3(axisA+pos_offset),d_tr4(axisA+pos_offset),
     3         d_tr1(axisB+dynaxis),d_tr2(axisB+dynaxis),
     4         d_tr3(axisB+dynaxis),d_tr4(axisB+dynaxis),
     5         d_v2,d_lambda1,d_lambda2,d_lambda3,d_lambda4,
     6         k_f1_offset,k_a0_offset,
     6         k_d1_offset(axisA-(axis_start-1)),
     7         k_d1_offset(axisB-(axis_start-1)),
     7         k_t1_offset,k_t2_offset,
     8         k_t3_offset,k_t4_offset,k_tr1_offset(axisA),
     9         k_tr2_offset(axisA),k_tr3_offset(axisA),
     1         k_tr4_offset(axisA),k_tr1_offset(axisB),
     2         k_tr2_offset(axisB),k_tr3_offset(axisB),
     3         k_tr4_offset(axisB),k_v2_offset,k_l1_offset,
     4         k_l2_offset,k_l3_offset,k_l4_offset)
!
            call reconcilefile(d_a0,size_a0)
            call get_block(d_a0,alpha2,1,0)
            call tce_zero(d_a0,size_a0)
          endif ! axisA.eq.axisB
!              write(LuOut,*) "alpha2 = ",alpha2
!
          alpha(axisA,axisB)=(-0.5d0)*(alpha1+alpha2)
!
        call deletefile(d_a0)
!
!            write(LuOut,*) "alpha(axisA,axisB) = ",alpha(axisA,axisB)
!            if (nodezero) write(LuOut,9020) cpu, wall
!
      endif ! respaxis(axis)
      enddo ! axisB loop
      enddo ! axisA loop
      cpu=cpu+util_cpusec()
      wall=wall+util_wallsec()
!
!
      if (nodezero) then
          if (l_pole.eq.3) then
          write(LuOut,9343) "CCSDTQ Linear Response due to octupole",
     1                       dabs(omega)
          else if (l_pole.eq.2) then
          write(LuOut,9343) "CCSDTQ Linear Response due to quadrupole",
     1                       dabs(omega)
          endif
          write(LuOut,9344)
          do axisA = axis_start,axis_max
            do axisB = axis_start, axisA
            if (respaxis(axisA).and.respaxis(axisB)) then
               write(LuOut,9444) trim(axisname(axisA)),axisname(axisB),
     1                alpha(axisA,axisB),alpha(axisA,axisB)*au2ang3
            endif
            enddo
          enddo
          write(LuOut,9344)
          write(LuOUt,*) ""
      endif
      
      if (nodezero) write(LuOut,9020) cpu, wall
      call util_flush(LuOut)

 9020 format(1x,'Cpu & wall time / sec',2f15.1)
 9343 format(/,1x,A,/,/,1x,'Frequency =',f15.7,' / au',/
     1  1x,'---------|---------------------|---------------------',/
     2  1x,'  Axes   |        bohr^3       |       angstroms^3   ')
 9344 format(1x,'---------|---------------------|---------------------')
 9444 format(1x,' ',A3,'/',A3,' | ',f16.7,3x,' | ',f16.7,3x)
      return
      end
