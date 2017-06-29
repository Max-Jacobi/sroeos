!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SRO_EOS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE Solve_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_constants_Mod
  USE Phase_Space_Input_Mod
  USE ALLOCATABLE_Tables_Mod
  USE READ_NUCLEAR_DATA_TABLE_MOD, ONLY : nuclei
  USE Make_Tables_Mod, ONLY : output_directory
  USE NSE_EOS_MOD
  USE Isotopes_Mod
  USE share
  USE omp_lib

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SOLVE

    IMPLICIT NONE

!   some array sizes
    INTEGER(I4B) :: nn, nt, ny, nn3, nt3, ny3

!   allocatable arrays to store solutions
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: logn, logt, ye, logn3, logt3, ye3
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: p3,s3,e3,muh3,mun3,mup3
    LOGICAL,  DIMENSION(:,:,:), ALLOCATABLE :: fix

!   name of output files
    CHARACTER(LEN=64) :: filename0, filename1, filename2
    CHARACTER(LEN=64) :: command, string_ye

!   loop variables and auxiliary integers
    INTEGER(I4B)  :: i,j,k,i0,j0,k0,i1,j1,k1,i3,j3,k3,l,thread

!
    REAL(DP), ALLOCATABLE :: xmass(:)

!   spacing used to compute numerical derivatives in log10n, log10T, Yp
    REAL(DP), PARAMETER   :: dlnn = 1.0D-3, dlnt = 1.0D-3, dy = 1.0D-3

    REAL(DP)   :: xy3, xn3, xt3
    REAL(DP)   :: xn, xt, xy
    REAL(DP)   :: press, eps, entropy, mu_hat
    REAL(DP)   :: mu_n, mu_p

    REAL(DP)   :: xxn, xxp, xxa, xxh, xxl
    REAL(DP)   :: xahbar, xzhbar, xalbar, xzlbar

    REAL(DP)   :: mu_p1, mu_p2, mu_n1, mu_n2
    REAL(DP)   :: dmup_dt, dmun_dt
    REAL(DP)   :: xt_1, xt_2, delta_t

!
    INTEGER(I4B) :: modi, modj, modk, modi0, modj0, modk0
    LOGICAL :: loop

    ALLOCATE(xmass(num_isotopes))
    xmass(:) = 0.0d0

!   Set size of ALLOCATABLE arrays
!   Some arrays are 3 times their final size
!   to store values for the numerical derivatives.
    nn3 = 3*(pointsrho  - 1) + 1 ; nn = (nn3 + 2) / 3
    nt3 = 3*(pointstemp - 1) + 1 ; nt = (nt3 + 2) / 3
    ny3 = 3*(pointsye   - 1) + 1 ; ny = (ny3 + 2) / 3

    WRITE (*,*)
    WRITE (*,*) "Points in density:        ", nn
    WRITE (*,*) "Points in tempearture:    ", nt
    WRITE (*,*) "Points in protonfraction: ", ny
    WRITE (*,*)

    WRITE (*,*)
    WRITE (*,*) "Each thread is dynamically assigned 3 proton fractions:"
    WRITE (*,*) " The ones set in the input/space.in file as well as"
    WRITE (*,*) " one slightly below and slightly above the ones set"
    WRITE (*,*) " that are used to compute the numerical derivatives."
    WRITE (*,*)

    ALLOCATE(fix(nn,nt,ny))

!   Open ascii output files to output error messages
    command = 'mkdir -p '//adjustl(trim(output_directory)) 
    CALL SYSTEM(command)
    IF (write_solutions_to_file) THEN
      command = 'mkdir -p '//adjustl(trim(output_directory))//'/NO_SOL'
      CALL SYSTEM(command)
      command = 'mkdir -p '//adjustl(trim(output_directory))//'/ERROR'
      CALL SYSTEM(command)
      command = 'mkdir -p '//adjustl(trim(output_directory))//'/SOL'
      CALL SYSTEM(command)
    ENDIF
!   Iterate over (Yp,T,n) grid and save output arrays
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,3) DEFAULT(SHARED)&
!$OMP PRIVATE(xy,xn,xt,loop,xmass,press,eps,entropy,mu_n,mu_p,mu_hat) &
!$OMP PRIVATE(i,j,i0,j0,k0,i1,j1,k1,i3,j3,k3,l,thread) &
!$OMP PRIVATE(modi,modj,modk,modi0,modj0,modk0) &
!$OMP PRIVATE(xxn,xxp,xxa,xxl,xxh,xahbar,xzhbar,xalbar,xzlbar) &
!$OMP PRIVATE(filename0,filename1,filename2,string_ye) &
!$OMP PRIVATE(xy3,xn3,xt3,p3,s3,e3,muh3,mup3,mun3) &
!$OMP PRIVATE(xt_1,xt_2,delta_t) &
!$OMP PRIVATE(mu_p1,mu_p2,dmup_dt) &
!$OMP PRIVATE(mu_n1,mu_n2,dmun_dt) &
!$OMP SHARED(p_tab,s_tab,e_tab) &
!$OMP SHARED(muh_tab,mun_tab,mup_tab) &
!$OMP SHARED(dpdn_tab,dsdn_tab,dmudn_tab) &
!$OMP SHARED(dpdt_tab,dsdt_tab,dmudt_tab) &
!$OMP SHARED(dpdy_tab,dsdy_tab,dmudy_tab) &
!$OMP SHARED(ah_tab,zh_tab,al_tab,zl_tab) &
!$OMP SHARED(xn_tab,xp_tab,xa_tab,xh_tab,xl_tab,fix)
  DO k=ny3+1,0,-1 !0,ny3+1,1
    ! ALLOCATE auxiliary arrays and set them to zero
    ALLOCATE(  p3(0:nn3+1,0:nt3+1,0:2),   s3(0:nn3+1,0:nt3+1,0:2),   e3(0:nn3+1,0:nt3+1,0:2))
    ALLOCATE(muh3(0:nn3+1,0:nt3+1,0:2), mun3(0:nn3+1,0:nt3+1,0:2), mup3(0:nn3+1,0:nt3+1,0:2))
     p3 = zero ; s3 = zero ; e3 = zero ; muh3 = zero ; mun3 = zero ; mup3 = zero
    ! get thread number (one for each proton fraction)
    thread = omp_get_thread_num()
    ! compute proton fraction for given point
    modk = mod(k,3) ; k0 = (k+2) ; k3 = k0/3 ; modk0 = mod(k0,3)
    IF (modk==0) xy = Yp_min + (Yp_max-Yp_min)*REAL(k-0)/REAL(ny3-1) - dy/2.d0
    IF (modk==1) xy = Yp_min + (Yp_max-Yp_min)*REAL(k-1)/REAL(ny3-1)
    IF (modk==2) xy = Yp_min + (Yp_max-Yp_min)*REAL(k-2)/REAL(ny3-1) + dy/2.d0

!   get filenames for solution and error outputs
    IF (write_solutions_to_file) THEN
      WRITE (string_ye,"(1F8.6)") xy

      filename0 =  TRIM(ADJUSTL(output_directory))//"/NO_SOL/"//TRIM(ADJUSTL(string_ye))
      filename1 = TRIM(ADJUSTL(output_directory))//"/ERROR/"//TRIM(ADJUSTL(string_ye))
      filename2 = TRIM(ADJUSTL(output_directory))//"/SOL/"//TRIM(ADJUSTL(string_ye))

      OPEN(1000*thread+10,FILE=filename0,STATUS='replace')
      OPEN(1000*thread+11,FILE=filename1,STATUS='replace')
      OPEN(1000*thread+12,FILE=filename2,STATUS='replace')
    ENDIF

    WRITE (*,"('Thread',1i2,' computing loop ',1i4,' of ',1i4, ' for y = ',1F8.6)") &
          thread, k, ny3+1, xy

    DO i=nn3+1,0,-1
      modi = mod(i,3) ; i0 = (i+2) ; i3 = i0/3 ; modi0 = mod(i0,3)

      IF (modi==0) xn = 10.D0**(Log10n_min + (Log10n_max-Log10n_min)*REAL(i-0)/REAL(nn3-1) - dlnn/2.d0)
      IF (modi==1) xn = 10.D0**(Log10n_min + (Log10n_max-Log10n_min)*REAL(i-1)/REAL(nn3-1))
      IF (modi==2) xn = 10.D0**(Log10n_min + (Log10n_max-Log10n_min)*REAL(i-2)/REAL(nn3-1) + dlnn/2.d0)

      xn3 = 10.D0**(Log10n_min + (Log10n_max - Log10n_min)*REAL(i-1)/REAL(nn3-1))

      ! dummies used to compute temperature derivatives of chemical potentials
      ! used to speed up computation of solutions at lower T
      mu_p1 = 0.d0 ; mu_p2 = 0.d0
      mu_n1 = 0.d0 ; mu_n2 = 0.d0
      xt_1  = 0.d0 ; xt_2  = 0.d0

      DO j=nt3+1,0,-1
        ! set counters for number of reduction and reorder operations
        nr_calls = 0
        no_calls = 0

        modj = mod(j,3) ; j0 = j+2 ; j3 = j0/3 ; modj0= mod(j0,3)

        IF (modj==0) xt = 10.D0**(Log10T_min + (Log10T_max - Log10T_min)*REAL(j-0)/REAL(nt3-1) - dlnt/2.d0)
        IF (modj==1) xt = 10.D0**(Log10T_min + (Log10T_max - Log10T_min)*REAL(j-1)/REAL(nt3-1))
        IF (modj==2) xt = 10.D0**(Log10T_min + (Log10T_max - Log10T_min)*REAL(j-2)/REAL(nt3-1) + dlnt/2.d0)

        xt3 = 10.D0**(Log10T_min + (Log10T_max - Log10T_min)*REAL(j-1)/REAL(nt3-1))

        loop = .false.

        IF ( (modi==1) .and. (modj==1) ) loop = .true.
        IF ( (modj==1) .and. (modk==1) ) loop = .true.
        IF ( (modk==1) .and. (modi==1) ) loop = .true.

        IF (.not.loop) CYCLE

        delta_t = xt - xt_1
        mu_p = mu_p + dmup_dt*delta_t
        mu_n = mu_n + dmun_dt*delta_t

        CALL CALC_IONEOS(xn,xt,xy,xmass,press,eps,entropy,mu_n,mu_p,mu_hat)

        p3(i,j,modk)  = press
        s3(i,j,modk)  = entropy
        e3(i,j,modk)  = eps

        muh3(i,j,modk)  = mu_hat
        mun3(i,j,modk)  = mu_n
        mup3(i,j,modk)  = mu_p

!       calculate light nuclei fractions
!       and their average charge and mass numbers
        xxl = 0.d0 ; xalbar = 0.d0 ; xzlbar = 0.d0

        DO l = llmin, llmax
          ! if proton pr neutron do not count
          IF (nuclei(l)%A==1.d0) THEN
            ! obtain neutron fraction
            IF (nuclei(l)%Z==0.d0) xxn = xmass(l)
            ! obtain proton fraction
            IF (nuclei(l)%Z==1.d0) xxp = xmass(l)
            CYCLE
          ENDIF
          ! obtain alpha particle fraction
          IF (nuclei(l)%Z==2.d0.and.nuclei(l)%A==4.d0) THEN
            xxa = xmass(l)
            CYCLE
          ENDIF
          xalbar = xalbar + nuclei(l)%A*xmass(l)
          xzlbar = xzlbar + nuclei(l)%Z*xmass(l)
          xxl    = xxl    + xmass(l)
        ENDDO

        IF (xxl>0.d0) THEN
          xalbar = xalbar/xxl ; xzlbar = xzlbar/xxl
        ELSE
          xalbar = 0.d0 ; xzlbar = 0.d0
        ENDIF

!       calculate heavy nuclei fractions
!       and their average charge and mass numbers
        xxh = 0.d0 ; xahbar = 0.d0 ; xzhbar = 0.d0

        DO l = llmax+1, num_isotopes
          xahbar = xahbar + nuclei(l)%A*xmass(l)
          xzhbar = xzhbar + nuclei(l)%Z*xmass(l)
          xxh   = xxh   + xmass(l)
        ENDDO

        IF (xxh>0.d0) THEN
          xahbar = xahbar/xxh ; xzhbar = xzhbar/xxh
        ELSE
          xahbar = 0.d0 ; xzhbar = 0.d0
        ENDIF

!       check if solution is acceptable
!       if not record array values to fix later
!       by interpolating using nearby points.
        IF (ABS(1.d0-(xxp+xxn+xxh+xxa+xxl))>1.d-6) THEN
          IF (write_solutions_to_file) THEN
            IF (write_solutions_to_file) &
              WRITE (1000*thread+10,"(11ES13.5,2I6)") &
              xn, xt, xy, mu_n, mu_p, ABS(1.d0-xxp-xxn-xxh-xxa-xxl), &
              MAX(xxn,1.0d-99), MAX(xxp,1.0d-99), &
              MAX(xxa,1.0d-99), MAX(xxl,1.0d-99), &
              MAX(xxh,1.0d-99), nr_calls, no_calls
            FLUSH(1000*thread+10)
          ENDIF
          IF (modi==0) i1 = i
          IF (modi==1) i1 = (i-1)
          IF (modi==2) i1 = (i-2)
          i1 = i1/3+1
          IF (modj==0) j1 = j
          IF (modj==1) j1 = (j-1)
          IF (modj==2) j1 = (j-2)
          j1 = j1/3+1
          IF (modk==0) k1 = k
          IF (modk==1) k1 = (k-1)
          IF (modk==2) k1 = (k-2)
          k1 = k1/3+1
          fix(i1,j1,k1) = .true.
        ENDIF

        IF (modk0==0 .and. modj0==0 .and. modi0==0) THEN
          p_tab(i3,j3,k3)  = p3(i,j,1)
          s_tab(i3,j3,k3)  = s3(i,j,1)
          e_tab(i3,j3,k3)  = e3(i,j,1)

          muh_tab(i3,j3,k3) = muh3(i,j,1)
          mun_tab(i3,j3,k3) = mun3(i,j,1)
          mup_tab(i3,j3,k3) = mup3(i,j,1)

          xn_tab(i3,j3,k3) = xxn
          xp_tab(i3,j3,k3) = xxp
          xa_tab(i3,j3,k3) = xxa
          xh_tab(i3,j3,k3) = xxh
          xl_tab(i3,j3,k3) = xxl
          zh_tab(i3,j3,k3) = xzhbar
          ah_tab(i3,j3,k3) = xahbar
          zl_tab(i3,j3,k3) = xzlbar
          al_tab(i3,j3,k3) = xalbar
        ENDIF

        IF (modk==2) THEN
          dpdn_tab (i3,j3,k3) = (p3(i+1,j,1)-p3(i-1,j,1))/dlnn/xn3
          dsdn_tab (i3,j3,k3) = (s3(i+1,j,1)-s3(i-1,j,1))/dlnn/xn3
          dmudn_tab(i3,j3,k3) = (muh3(i+1,j,1)-muh3(i-1,j,1))/dlnn/xn3

          dpdt_tab (i3,j3,k3) = (p3(i,j+1,1)-p3(i,j-1,1))/dlnt/xt3
          dsdt_tab (i3,j3,k3) = (s3(i,j+1,1)-s3(i,j-1,1))/dlnt/xt3
          dmudt_tab(i3,j3,k3) = (muh3(i,j+1,1)-muh3(i,j-1,1))/dlnt/xt3

          dpdy_tab (i3,j3,k3) = (p3(i,j,2)-p3(i,j,0))/dy
          dsdy_tab (i3,j3,k3) = (s3(i,j,2)-s3(i,j,0))/dy
          dmudy_tab(i3,j3,k3) = (muh3(i,j,2)-muh3(i,j,0))/dy
        ENDIF

        IF (write_solutions_to_file) &
          WRITE (1000*thread+12,"(11ES13.5,2I6)") &
          xn, xt, xy, mu_n, mu_p, ABS(1.d0-xxp-xxn-xxh-xxa-xxl), &
          MAX(xxn,1.0d-99), MAX(xxp,1.0d-99), &
          MAX(xxa,1.0d-99), MAX(xxl,1.0d-99), &
          MAX(xxh,1.0d-99), nr_calls, no_calls

!       compute temperature derivatives to use
!       as guess for next temperature iteration
!       should improve seeking solution at
!       low Temperatures significantly
        xt_2  = xt_1
        mu_p2 = mu_p1
        mu_n2 = mu_n1

        xt_1  = xt
        mu_p1 = mu_p
        mu_n1 = mu_n

        dmup_dt = (mu_p1 - mu_p2) / (xt_1 - xt_2)
        dmun_dt = (mu_n1 - mu_n2) / (xt_1 - xt_2)
      ENDDO
      IF (write_solutions_to_file) THEN
        FLUSH(1000*thread+10) ; FLUSH(1000*thread+11) ; FLUSH(1000*thread+12)
      ENDIF
    ENDDO
    DEALLOCATE(p3,s3,e3,muh3,mun3,mup3)
    IF (write_solutions_to_file) THEN
      CLOSE(1000*thread+10) ; CLOSE(1000*thread+11) ; CLOSE(1000*thread+12)
    ENDIF
  ENDDO
!$OMP END PARALLEL DO


! set Ye array
  write (*,*) 'write ye array'
  DO i = 1, ny
    yp_tab(i) = Yp_min + (Yp_max-Yp_min)*REAL(i-1)/REAL(ny-1)
  ENDDO

! set log10(n) array
  write (*,*) 'write log10(n) array'
  DO i = 1, nn
    log10n_tab(i) = Log10n_min + (Log10n_max-Log10n_min)*REAL(i-1)/REAL(nn-1)
  ENDDO

! set log10(T) array
  write (*,*) 'write log10t array'
  DO i = 1, nt
    log10t_tab(i) = Log10T_min + (Log10T_max-Log10T_min)*REAL(i-1)/REAL(nt-1)
  ENDDO

! TODO: put this into a subroutine
! if no solution found for some point use interpolation in density
!  to have an average value for that point
  WRITE (*,*) "Adjusting points with problems."
  WRITE (*,*) "See file 'points_tried_to_fix.dat' for details. "
  OPEN(10,FILE=TRIM(ADJUSTL(output_directory))//"/points_tried_to_fix.dat",status='replace')
  WRITE (10,"(A50)") 'Trying to fix values at (log10(n),log10(T),Ye):'
  DO k = 2, ny-1
    WRITE (*,*) 'Iteration', k-1, ' of ', ny-2
    DO j = 2, nt-1
      DO i = 2, nn-1
        IF (fix(i,j,k)) WRITE (10,"(3ES14.6)") 10.d0**log10n_tab(i), 10.d0*log10t_tab(j), yp_tab(k)
        IF (fix(i,j,k) .and. .not.fix(i-1,j,k) .and. .not.fix(i+1,j,k)) THEN

         fix (i,j,k) = .false.

         p_tab(i,j,k) = (p_tab(i+1,j,k)+p_tab(i-1,j,k))/2.d0
         s_tab(i,j,k) = (s_tab(i+1,j,k)+s_tab(i-1,j,k))/2.d0
         e_tab(i,j,k) = (e_tab(i+1,j,k)+e_tab(i-1,j,k))/2.d0

         mun_tab(i,j,k) = (mun_tab(i+1,j,k)+mun_tab(i-1,j,k))/2.d0
         mup_tab(i,j,k) = (mup_tab(i+1,j,k)+mup_tab(i-1,j,k))/2.d0
         muh_tab(i,j,k) = (muh_tab(i+1,j,k)+muh_tab(i-1,j,k))/2.d0

         xn_tab(i,j,k) = (xn_tab(i+1,j,k)+xn_tab(i-1,j,k))/2.d0
         xp_tab(i,j,k) = (xp_tab(i+1,j,k)+xp_tab(i-1,j,k))/2.d0
         xa_tab(i,j,k) = (xa_tab(i+1,j,k)+xa_tab(i-1,j,k))/2.d0
         xh_tab(i,j,k) = (xh_tab(i+1,j,k)+xh_tab(i-1,j,k))/2.d0
         xl_tab(i,j,k) = (xl_tab(i+1,j,k)+xl_tab(i-1,j,k))/2.d0
         zh_tab(i,j,k) = (zh_tab(i+1,j,k)+zh_tab(i-1,j,k))/2.d0
         ah_tab(i,j,k) = (ah_tab(i+1,j,k)+ah_tab(i-1,j,k))/2.d0
         zl_tab(i,j,k) = (zl_tab(i+1,j,k)+zl_tab(i-1,j,k))/2.d0
         al_tab(i,j,k) = (al_tab(i+1,j,k)+al_tab(i-1,j,k))/2.d0

         dpdn_tab(i,j,k) = (dpdn_tab(i-1,j,k)+dpdn_tab(i+1,j,k))/2.d0
         dpdt_tab(i,j,k) = (dpdt_tab(i-1,j,k)+dpdt_tab(i+1,j,k))/2.d0
         dpdy_tab(i,j,k) = (dpdy_tab(i-1,j,k)+dpdy_tab(i+1,j,k))/2.d0

         dsdn_tab(i,j,k) = (dsdn_tab(i-1,j,k)+dsdn_tab(i+1,j,k))/2.d0
         dsdt_tab(i,j,k) = (dsdt_tab(i-1,j,k)+dsdt_tab(i+1,j,k))/2.d0
         dsdy_tab(i,j,k) = (dsdy_tab(i-1,j,k)+dsdy_tab(i+1,j,k))/2.d0

         dmudn_tab(i,j,k) = (dmudn_tab(i-1,j,k)+dmudn_tab(i+1,j,k))/2.d0
         dmudt_tab(i,j,k) = (dmudt_tab(i-1,j,k)+dmudt_tab(i+1,j,k))/2.d0
         dmudy_tab(i,j,k) = (dmudy_tab(i-1,j,k)+dmudy_tab(i+1,j,k))/2.d0

        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CLOSE(10)

! check for points still not fixed
  OPEN(10,FILE=TRIM(ADJUSTL(output_directory))//"/not_fixed.dat")
  WRITE (10,"(A50)") 'Problems at (log10(n),log10(T),Ye):'
  do i = 2, nn-1
    do j = 2, nt-1
      do k = 2, ny-1
        IF (fix(i,j,k)) WRITE (*,"(3ES14.6)") 10.d0**log10n_tab(i), 10.d0**log10t_tab(j), yp_tab(k)
      ENDDO
    ENDDO
  ENDDO
  CLOSE(10)

END SUBROUTINE SOLVE

END MODULE SOLVE_MOD
