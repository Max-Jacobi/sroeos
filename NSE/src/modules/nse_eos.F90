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
MODULE NSE_EOS_MOD

  USE Kind_types_Mod, ONLY : I4B, DP
  USE FUNCTIONS_MOD
  USE READ_NUCLEAR_DATA_TABLE_MOD, ONLY : nuclei
  USE SHARE
  USE Make_Tables_Mod, ONLY : write_solutions_to_file

  IMPLICIT NONE

  INTEGER(I4B), PRIVATE :: num_isotopes = 0
!$OMP THREADPRIVATE(num_isotopes)
  REAL(DP), PRIVATE, DIMENSION(:), allocatable :: aion_nse, zion_nse
!$OMP THREADPRIVATE(aion_nse, zion_nse)
  REAL(DP), PRIVATE, DIMENSION(:), allocatable :: wion_nse, DwionDT_nse
!$OMP THREADPRIVATE(wion_nse,DwionDT_nse)
  REAL(DP), PRIVATE, DIMENSION(:), allocatable :: mion_nse, beion_nse
!$OMP THREADPRIVATE(mion_nse, beion_nse)

CONTAINS

  SUBROUTINE local_nse(xn, xtemp, xye, xmass, mu_n, mu_p)

    USE Physical_Constants_Mod
    USE Fermi_Integrals_Mod
#ifdef _OPENMP
    USE omp_lib
#endif
    
    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: xn, xtemp, xye
    REAL(DP), INTENT(INOUT), DIMENSION(SIZE(nuclei))  :: xmass
    REAL(DP), INTENT(INOUT) :: mu_n, mu_p
    REAL(DP), DIMENSION(SIZE(nuclei)) :: xmasstemp

    INTEGER(I4B) :: i, j, k, l
    INTEGER(I4B) :: err, errj, errk
    INTEGER(I4B) :: thread
    INTEGER(I4B), SAVE :: ineut, iprot

    REAL(DP) :: mui, ni, nn, np, Tc, Tstart, Ec
    REAL(DP) :: Tt, munt, mupt, sumy
    REAL(DP) :: aion, zion, nion, mion, wion, massk
    REAL(DP) :: f(2), jac(2,2), det

    LOGICAL :: logy

    REAL(DP),  SAVE :: muno, mupo, To
    REAL(DP), DIMENSION(3335), SAVE :: ymass

    INTEGER(I4B), DIMENSION(3335), SAVE :: order, order_old
    INTEGER(I4B), SAVE :: kmax, kmax_old
!$OMP THREADPRIVATE(kmax,ymass,order,order_old,muno,mupo,To)

#ifdef _OPENMP    
    thread = omp_get_thread_num()
#else
    thread = 1
#endif
    
!   find proton and neutron locations
!   TODO: move this nse_read_nuclear_data module so it is done only once
!         add some way to make sure that ineut = 1 and iprot = 2
!         as this is assumed below, but not checked if true.
    DO i=1,num_isotopes
      IF (aion_nse(i)==1) THEN
        IF (zion_nse(i)==1) THEN
          iprot = i
        ELSE
          ineut = i
        ENDIF
      ENDIF
    ENDDO

!   get order of most abundant isotopes from previous step
    order_old = order
    kmax_old  = kmax

!   if T > 200MeV only use lightest 10 isotopes
    IF (xtemp>2.d2) THEN
      kmax = num_isotopes
      kmax = min(10,kmax)
      DO j = 1, 10 ; order(j) = j ; ENDDO
    ENDIF

!   set solver to use first 10 isotopes in list
!   plus any others with number fraction larger than 1.d-16

!   set error = -1. err < 0 Implies no solution found.
!                   err = 0 Implies a solution was found
    err = -1

!   set Tc to input temperature xtemp
    Tc = xtemp

!   Compute solution for current T
!   using as start value solution from last call
    err = nr_iterate()

!   If no solution found use as initial guess
!    matter containing only protons and neutrons
    IF (err .ne. 0) THEN
      mu_n = calc_mu((1.d0-xye)*xn,mion_nse(ineut),wion_nse(ineut),xtemp,0.d0)
      mu_p = calc_mu(      xye *xn,mion_nse(iprot),wion_nse(iprot),xtemp,0.d0)
      err = nr_iterate()
    ENDIF

!   If no solution found use as initial guess
!    matter containing only protons and neutrons
    IF (err .ne. 0) THEN
      ! compute total number fraction and check if within 1e-6 of 1.
      sumy = sum(ymass(1:kmax))
      IF (abs(1.d0-sumy)<1.d-6) logy = .true.
      do l=3,kmax
        k = order(l)
        if (k==0) exit
        massk = xmass(k)
        if (massk > 1.d0 .or. massk < 1.d-20) exit
        aion = aion_nse(k)
        zion = zion_nse(k)
        nion = aion - zion
        mion = mion_nse(k)
        wion = wion_nse(k)

        nn = (1.d0-xye-xye*nion/zion)*xn
        np = (1.d0 - zion/nion*(1.d0 + xye))*xn
        if (nn>0.d0) then
          ni = xn*xye/zion
          mu_n = calc_mu(nn,mion_nse(ineut),wion_nse(ineut),xtemp,0.d0)
          Ec = calc_Ec(xn, xtemp, xye, aion, zion, 0.16d0)
          mui  = calc_mu(ni,mion,wion,xtemp,Ec)
          mu_p = (mui - mu_n*nion)/zion
        else
          ni = xn*(1.d0-xye)/nion
          mu_p = calc_mu(np,mion_nse(iprot),wion_nse(iprot),xtemp,0.d0)
          Ec = calc_Ec(xn, xtemp, xye, aion, zion, 0.16d0)
          mui  = calc_mu(ni,mion,wion,xtemp,Ec)
          mu_n = (mui - mu_p*zion)/nion
        endif
        err = nr_iterate()
        if (err==0) exit
      enddo
    endif

!   If no solution found using methods above
!   use previous higher T solution and reduce temperature
!   in 10 steps to reach desired temperature
101 continue
    if (err .ne. 0) then
      order = order_old
      kmax  = kmax_old
      Tstart = To
      mu_n = muno !calc_mu((1.d0-xye)*xn,mion_nse(ineut),wion_nse(ineut),Tstart,0.d0)
      mu_p = mupo !calc_mu(      xye*xn,mion_nse(iprot),wion_nse(iprot),Tstart,0.d0)
      do k=11,1,-1
        Tc = 10.d0**(LOG10(xtemp)+(DBLE(k-1)/10.d0)*(LOG10(Tstart)-LOG10(xtemp)))
        err = nr_iterate()
        if (err .ne. 0) then
          mu_n = munt
          mu_p = mupt
          Tc = Tt
          exit
        endif
        ! SAVE last one that suceeded
        munt = mu_n
        mupt = mu_p
        Tt = Tc
      enddo
    endif

!   If no solution found using methods above
!   start from a high temperature and reduce it
!   in small steps to reach desired temperature
!   (Sometimes this succedes even if loop above does not.)
102 continue
    if (err .ne. 0) then
      Tstart = 1.D1**2.4D0
      kmax = num_isotopes
      kmax = min(10,kmax)
      do j = 1, kmax ; order(j) = j ; enddo
      mu_n = calc_mu((1.d0-xye)*xn,mion_nse(ineut),wion_nse(ineut),Tstart,0.d0)
      mu_p = calc_mu(      xye *xn,mion_nse(iprot),wion_nse(iprot),Tstart,0.d0)
      do k=100,0,-1
        Tc = 2.4D0 + (LOG10(XTEMP)-2.4D0)/REAL(100)*REAL(100-k)
        Tc = 1.d1**Tc
        err = nr_iterate()
        if (err .ne. 0) then
          mu_n = munt
          mu_p = mupt
          Tc = Tt
          exit
        endif
        ! SAVE last one that suceeded
        munt = mu_n
        mupt = mu_p
        Tt = Tc
      enddo
    endif

    IF (err .ne. 0 .AND. write_solutions_to_file) THEN
      WRITE (1000*thread+11,"(1F7.4,2ES15.7,I2)") xye, xn, xtemp, err
    ELSE
      muno = mu_n
      mupo = mu_p
      To   = xtemp
    ENDIF

!   calculate mass fractions
    errk = mass_frac()

!   order isotopes
    errj = reorder()

    RETURN

CONTAINS

    FUNCTION nr_iterate() RESULT (ierr)
      integer :: m, n, ierr, imax, iini, iter
      REAL(DP) :: tol, tol1, tol2, tol1old, tol2old
      REAL(DP) :: aion, zion, nion, mion, wion
      REAL(DP) :: calc_n,calc_Ec

      imax = 100
      tol1old = 1.d50 ;  tol2old = 1.d50
      ! set tolerance to be smaller for lower temperatures
      ! this seems to work fine for now
      tol = 1.d-9
      if (Tc<100.d0) tol = 1.d0*10.d0**(- 8.d0 - log10(Tc)/2.d0)
      ! NR iteration
      iini = 1
      iter = 0
102   ierr = 0
      do i=iini,imax
        nr_calls = nr_calls+1
        f = 0.d0
        jac = 0.d0

        errk = mass_frac()

        jac = jac/Tc
        det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)

        mu_n = mu_n - (f(1)*jac(2,2) - f(2)*jac(1,2))/det
        mu_p = mu_p - (f(2)*jac(1,1) - f(1)*jac(2,1))/det

        tol1 = abs(f(1))/xn ; tol2 = abs(f(2))/(xn*xye)
        if (tol1*tol1+tol2*tol2<tol*tol) exit
        if ((mu_n .ne. mu_n) .or. (mu_p .ne. mu_p)) then
          ierr = 2
          exit
        endif
      enddo

!   Calculate mass fractions
!     used to recalculate species that contribute the most
      if (ierr == 0) then
      f = 0.d0
!     calculate mass fractions
      errk = mass_frac()

!     reorder nuclei
      errj = reorder()

      tol1 = abs(f(1))/xn ; tol2 = abs(f(2))/(xn*xye)
!     check if number of iterations larger than max of it tolerance has been achieved
      if (i>imax .or. (tol1*tol1+tol2*tol2)>tol*tol) then
        ierr = 1
        iini = i+1
        ! iter = iter + 1
        ! try another 50 steps with updated list of isotopes
        if (imax>599) then
          !make tolerance larger for low T.
          if (Tc<1.d-2 .and. (tol1*tol1+tol2*tol2)<4.d0*tol*tol) ierr = 0
          return
        endif
        imax = iini+100
        ! do at most 10 retrials
        if (iter<10) then
          iter = iter + 1
          goto 102
        endif
      endif
      endif
    end FUNCTION nr_iterate

    FUNCTION reorder() RESULT(jerr)

      integer :: ij, ik, jerr

!     get list of isotopes with number fraction larger than 1.d-16
      xmasstemp = xmass
      ymass(:) = 0.d0
      order(:) = 0
      no_calls = no_calls+1
      do ij = 1, num_isotopes
        if (ij<=10) then
          ik = ij
        else
          ik = maxloc(xmasstemp(11:num_isotopes),1)+10
          xmasstemp(ik) = -dble(ij)
          kmax = ij
        endif
        order(ij) = ik
        ymass(ij) = xmass(ik)
        if (ij>10.and.xmass(ik)<1.d-16) exit
      enddo

      jerr = 0

    END FUNCTION reorder

    FUNCTION mass_frac() RESULT(kerr)

      IMPLICIT NONE

      INTEGER(I4B) :: jk, kerr

    ! Calculate mass fractions
      f = 0.d0
      ymass(:) = 0.d0

      do jk=1,num_isotopes
        aion = aion_nse(jk)
        zion = zion_nse(jk)
        nion = aion - zion
        mion = mion_nse(jk)
        wion = wion_nse(jk)
        mui = nion*mu_n + zion*mu_p
        Ec = calc_Ec(xn, Tc, xye, aion, zion, 0.16d0)
        ni = calc_n(mui,mion,wion,Tc,Ec)
        f(1) = f(1) + aion*ni
        f(2) = f(2) + zion*ni
        xmass(jk) = ni / xn * aion
        jac(1,1) = jac(1,1) + aion*ni*nion
        jac(2,1) = jac(2,1) + zion*ni*nion
        jac(1,2) = jac(1,2) + aion*ni*zion
        jac(2,2) = jac(2,2) + zion*ni*zion
      enddo

      f(1) = f(1) - xn
      f(2) = f(2) - xn*xye

    END FUNCTION mass_frac

  END SUBROUTINE local_nse

  SUBROUTINE calc_ioneos(xn, xtemp, xye, xmass, press, eps, &
       entropy, mu_n, mu_p, mu_hat)

    !use electron_eos_mod
    use Physical_Constants_Mod
    use READ_NUCLEAR_DATA_TABLE_MOD
    use Fermi_Integrals_Mod
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none

    REAL(DP), INTENT(IN)    :: xn, xtemp, xye
    REAL(DP), INTENT(INOUT) :: xmass(SIZE(nuclei))
    REAL(DP), INTENT(out)   :: press, eps, entropy, mu_hat
    REAL(DP), INTENT(INOUT) :: mu_n, mu_p

    integer :: iiso, thread
    logical, SAVE :: already_found_n_p_indices = .false.
    integer, SAVE :: neutron_index, proton_index
!$OMP THREADPRIVATE(already_found_n_p_indices)
!$OMP THREADPRIVATE(neutron_index,proton_index)
    REAL(DP)  :: xni, xpressi, xepsi, xentropyi
    REAL(DP)  :: pressc,dpressc

    REAL(DP) rhocgs,TK,mu_i,mu_i_eq, Ec, dEcdln, d2Ecdln2
    REAL(DP) :: zfermi12, ifermi12, zfermi32, etat, tau, dz

#ifdef _OPENMP    
    thread = omp_get_thread_num()
#else
    thread = 1
#endif
    
    ! Intialize internal arrays, this is basically vestigial
    if (num_isotopes .ne. SIZE(nuclei) .or. (.not. allocated(aion_nse))) then

      if (allocated(aion_nse)) deallocate(aion_nse, zion_nse, wion_nse, &
          DwiondT_nse, mion_nse, beion_nse)

      num_isotopes = SIZE(nuclei)
      allocate(aion_nse(num_isotopes), zion_nse(num_isotopes), &
          wion_nse(num_isotopes), DwiondT_nse(num_isotopes), &
          mion_nse(num_isotopes), beion_nse(num_isotopes))

      do iiso=1,num_isotopes
        aion_nse(iiso) = nuclei(iiso)%A
        zion_nse(iiso) = nuclei(iiso)%Z
        mion_nse(iiso) = nuclei(iiso)%mass
        beion_nse(iiso) = nuclei(iiso)%be
      enddo

    endif

    ! get partition FUNCTIONs
    wion_nse(:) = 1.0d0
    DwiondT_nse(:) = 0.0d0
    call getpart(num_isotopes,xtemp,nuclei,wion_nse,DwiondT_nse)

    if (.not. already_found_n_p_indices) then
      ! find neutron and proton indices
      neutron_index = -1
      proton_index = -1

      do iiso=1,num_isotopes
        if (aion_nse(iiso) .eq. 1) then
          if (zion_nse(iiso) .eq. 0) then
            neutron_index = iiso
          elseif (zion_nse(iiso) .eq. 1) then
            proton_index = iiso
          endif
        endif
      enddo

      if (neutron_index .eq. -1) then
        write(6,*) "Could not find neutron index in ioneos"
        stop
      endif

      if (proton_index .eq. -1) then
        write(6,*) "Could not find proton index in ioneos"
        stop
      endif

      already_found_n_p_indices = .true.
    endif

    dpressc = 0.0d0
    pressc  = 0.0d0
    press   = 0.0d0
    eps     = 0.0d0
    entropy = 0.0d0

    call local_nse(xn,xtemp,xye,xmass,mu_n,mu_p)
    !if (ABS(SUM(xmass)-1.d0)>1.d-4) then
      !print *,'bad abundances passed to calc_ioneos', SUM(xmass)
    !endif

    do iiso=1,num_isotopes
      ! n_i = n X_i / A_i
      xni = xn * xmass(iiso) / aion_nse(iiso)

      Ec = calc_Ec(xn, xtemp, xye, aion_nse(iiso), zion_nse(iiso), 0.16d0)
      dEcdln = calc_dEcdln(xn, xtemp, xye, aion_nse(iiso), zion_nse(iiso), 0.16d0, d2Ecdln2)
      ! partial pressure : P_i = n_i * T
      xpressi   = xni * xtemp

      if (xmass(iiso)/aion_nse(iiso) > 1.d-20) then
      ! internal energy in MeV / fm^3 with neutron-proton correction
      ! and binding energy
        xepsi     = 1.5d0 * xtemp  + Ec &
              - zion_nse(iiso) * (neutron_mass_MeV - proton_mass_MeV) &
              - beion_nse(iiso) + DwionDT_nse(iiso)/wion_nse(iiso)*xtemp**2
      ! entropy per particle of mass mion_nse(iiso)
        xentropyi = 2.5d0 + log(wion_nse(iiso)/xni * &
            (mion_nse(iiso) * xtemp / sac_const)**(3.0/2.0)) &
            + DwionDT_nse(iiso)/wion_nse(iiso)*xtemp
      else ! If the abundance is very small, the species won't matter and we want to avoid NaNs
        xepsi     = 0.d0
        xentropyi = 0.d0
      endif

      ! Use RESULTs for arbitrarily degenerate nucleons
      if (aion_nse(iiso) == 1 .and. xni > 1.d-40) then
        etat = (calc_mu(xni, mion_nse(iiso), wion_nse(iiso), xtemp, 0.d0) &
            - mion_nse(iiso))/xtemp
        tau = (2.d0*mion_nse(iiso)*xtemp/HBC**2.d0)**2.5d0 &
            *wion_nse(iiso)/(4.d0*PI**2)*fermi_three_halves(etat)
        xepsi = HBC**2.d0*tau/(2.d0*mion_nse(iiso))/xni &
            - zion_nse(iiso)*(neutron_mass_MeV - proton_mass_MeV)
        xpressi = HBC**2*tau/(3.d0*mion_nse(iiso))
        xentropyi = 5.d0/(6.d0*mion_nse(iiso)*xtemp)*HBC**2.d0*tau/xni - etat
      endif

      dpressc = dpressc + xni*d2Ecdln2
      pressc  = pressc  + xni*dEcdln
      press   = press   + xpressi
      eps     = eps     + xepsi * xmass(iiso) / aion_nse(iiso)

      ! want: entropy per baryon
      entropy = entropy + xentropyi * xmass(iiso)  / aion_nse(iiso)

    enddo

    press = press + pressc

    mu_hat = mu_n - mu_p ! Not actually muhat, don't want to change names in code everywhere

  end SUBROUTINE calc_ioneos

  SUBROUTINE SWAP_R(A, B)
    IMPLICIT NONE
    REAL, INTENT (INOUT)      :: A, B
    REAL                      :: TEMP
    TEMP = A ; A = B ; B = TEMP
  END SUBROUTINE SWAP_R

end module NSE_EOS_MOD
