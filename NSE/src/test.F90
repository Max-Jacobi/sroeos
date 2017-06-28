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
program nse_test

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_constants_Mod
  USE Read_Input_Mod
  USE Test_Input_Mod
  USE READ_NUCLEAR_DATA_TABLE_MOD
  USE NSE_EOS_MOD
  USE SHARE
  USE Print_test_Mod

  implicit none

  ! this is a test program that will return
  ! internal energy, pressure, and entropy
  ! for the ion component in a given
  ! for a given rho, T, Ye composition in NSE.
  !
  ! To get this right, one needs to set
  ! in ioneos.F90: do_rad = .false. and
  ! do_ele = .false.
  !
  real(dp)   :: xtemp, xtot
  real(dp)   :: xye, xyeplus, xyeminus
  real(dp)   :: xn, xnplus, xnminus
  real(dp)   :: xt, xtplus, xtminus
  real(dp)   :: finish, start
  real(dp), allocatable :: xmass(:), xmasstemp(:)
  integer(i4b), allocatable :: order(:)
  real(dp)   :: press,eps,entropy, ep, em, dsn
  ! helpers/dummies
  real(dp)   :: mu_n, mu_p, mu_n0, mu_p0, mu_hat
  real(dp)   :: pplus, pminus, dpdn, dpdt, dpdy
  real(dp)   :: splus, sminus, dsdn, dsdt, dsdy
  real(dp)   :: eplus, eminus, dedn, dedt, dedy
  real(dp)   :: muhplus, muhminus, dmuhdn, dmuhdt, dmuhdy
  real(dp)   :: dmup_dt, dmun_dt, mu_p1, mu_p2, mu_n1, mu_n2
  real(dp)   :: xt_1, xt_2, delta_t
  logical    :: first_init_nsenet = .true.
  integer(i4b)  :: i,j,k,imax
  integer(i4b)  :: fail_nse_local, num_isotopes
  character(80) :: isolist_filename != "tables/list_206.set"
  character(80) :: part_filename = "tables/partition_MESA"
  character(80) :: mass_filename = "tables/mass_table_MESA"

  real(dp) :: x_n, x_p, x_a, x_l, x_h, A_l, Z_l, A_h, Z_h

  fail_nse_local = 0

  IS_TEST = .TRUE.

  CALL READ_TEST_INPUT

  xye = proton_fraction
  xtemp = temperature
  xn = density

  call read_nuclear_data()

  write(*,*)
  write(*,*) "*** successfully read nuclear data"
  write(*,*)
  num_isotopes = SIZE(nuclei)
  allocate(xmass(num_isotopes))
  xmass(:) = 0.0d0
  call cpu_time(start)
  imax = 100
  IF (LOG10(xtemp) >= 1.0D0) imax = 10

  WRITE (*,*) 'Computing solution from 10**2.4 MeV to ', xtemp, ' MeV in ', imax, 'log10(T) steps.'  

  mu_p1 = 0.d0 ; mu_p2 = 0.d0
  mu_n1 = 0.d0 ; mu_n2 = 0.d0
  xt_1  = 0.d0 ; xt_2  = 0.d0

  WRITE (*,*)
  WRITE (*,*) 'Writing solutions for higher Temperatures in TEST_SOLUTIONS.DAT file.'
  WRITE (*,*)

  OPEN (667,FILE='TEST_SOLUTIONS.DAT')
  WRITE (667,*) &
   '      T (MeV)              mu_n                mu_p         reduction calls   reorder calls'
  do i = imax, 0, -1
    nr_calls = 0
    no_calls = 0
    IF (imax>0) THEN
      xt = 2.4D0 + (LOG10(XTEMP)-2.4D0)/REAL(imax)*REAL(imax-i)
    ELSE
      xt = xtemp
    ENDIF
    xt = 1.d1**xt
    delta_t = xt - xt_1
    mu_p = mu_p + dmup_dt*delta_t
    mu_n = mu_n + dmun_dt*delta_t
    call calc_ioneos(xn,xt,xye,xmass,press,eps,entropy,mu_n,mu_p,mu_hat)
    WRITE (667,'(3ES20.12,2I16)') xt, mu_n, mu_p, nr_calls, no_calls
    xt_2  = xt_1
    mu_p2 = mu_p1
    mu_n2 = mu_n1
    xt_1  = xt
    mu_p1 = mu_p
    mu_n1 = mu_n
    dmup_dt = (mu_p1 - mu_p2) / (xt_1 - xt_2)
    dmun_dt = (mu_n1 - mu_n2) / (xt_1 - xt_2)
  enddo
  CLOSE(667)

  call cpu_time(finish)
  write(*,*)
  print '("Run time = ",f18.9," seconds.")',finish-start
  write(*,*)

  WRITE(*,*)
  WRITE(*,*) 'Computing composition.'
  WRITE(*,*)

  allocate(xmasstemp(num_isotopes))
  allocate(order(num_isotopes))
  xtot = 0.d0
  xmasstemp = xmass
  order = 0

  do j = 1, num_isotopes
    k = maxloc(xmasstemp(1:num_isotopes),1)
    xmasstemp(k) = -j
    imax = j
    order(j) = k
  enddo

  OPEN (667,FILE='TEST_COMPOSITION.DAT')
  WRITE (667,*) &
   '  order    isotope   mass fraction   cummulative mass fraction '

  WRITE(*,*)
  WRITE(*,*) 'Writing composition for isotopes with number fraction'
  WRITE(*,*) 'larger than 1.d-10 in TEST_COMPOSITION.DAT file.'
  WRITE(*,*)

  xtot = 0.d0

  x_l = 0.d0
  x_h = 0.d0

  A_l = 0.d0
  Z_l = 0.d0
  A_h = 0.d0
  Z_h = 0.d0

  do i=1,num_isotopes
     j = order(i)
     xtot = xtot + xmass(j)
     if (nuclei(j)%name == 'n') x_n = xmass(j)
     if (nuclei(j)%name == 'p') x_p = xmass(j)
     if (nuclei(j)%name == 'he4') x_a = xmass(j)
     if (nuclei(j)%z <= 6.d0 .and. nuclei(j)%name .ne. 'n' .and. &
         nuclei(j)%name .ne. 'p' .and. nuclei(j)%name .ne. 'he4') then
       x_l = x_l + xmass(j)
       A_l = A_l + xmass(j)*(nuclei(j)%n+nuclei(j)%z)
       Z_l = Z_l + xmass(j)*(nuclei(j)%z)
     endif
     if (nuclei(j)%z >  6.d0) then
       x_h = x_h + xmass(j)
       A_h = A_h + xmass(j)*(nuclei(j)%n+nuclei(j)%z)
       Z_h = Z_h + xmass(j)*(nuclei(j)%z)
     endif
     if(xmass(j).gt.1.0d-10.or.i<20) then
        write(667,"(I8,A9,1P10E20.10)") j,nuclei(j)%name,xmass(j),xtot
     endif
  enddo

  CLOSE(667)

  WRITE(*,*)
  WRITE(*,*) 'Computing temperature derivatives'
  WRITE(*,*)

  xtplus = 1.01d0*xt
  call calc_ioneos(xn,xtplus,xye,xmass,pplus,eplus,splus,mu_n0,mu_p0,muhplus)

  xtminus = 0.99d0*xt
  call calc_ioneos(xn,xtminus,xye,xmass,pminus,eminus,sminus,mu_n0,mu_p0,muhminus)

  dpdt = (pplus-pminus)/(xtplus-xtminus)
  dsdt = (splus-sminus)/(xtplus-xtminus)
  dedt = (eplus-eminus)/(xtplus-xtminus)
  dmuhdt = (muhplus-muhminus)/(xtplus-xtminus)

  WRITE(*,*)
  WRITE(*,*) 'Computing density derivatives'
  WRITE(*,*)

  xnplus = 1.01d0*xn
  call calc_ioneos(xnplus,xt,xye,xmass,pplus,eplus,splus,mu_n0,mu_p0,muhplus)

  xnminus = 0.99d0*xn
  call calc_ioneos(xnminus,xt,xye,xmass,pminus,eminus,sminus,mu_n0,mu_p0,muhminus)

  dpdn = (pplus-pminus)/(xnplus-xnminus)
  dsdn = (splus-sminus)/(xnplus-xnminus)
  dedn = (eplus-eminus)/(xnplus-xnminus)
  dmuhdn = (muhplus-muhminus)/(xnplus-xnminus)

  WRITE(*,*)
  WRITE(*,*) 'Computing proton fraction derivatives'
  WRITE(*,*)

  xyeplus = 1.01d0*xye
  call calc_ioneos(xn,xt,xyeplus,xmass,pplus,eplus,splus,mu_n0,mu_p0,muhplus)

  xyeminus = 0.99d0*xye
  call calc_ioneos(xn,xt,xyeminus,xmass,pminus,eminus,sminus,mu_n0,mu_p0,muhminus)

  dpdy = (pplus-pminus)/(xyeplus-xyeminus)
  dsdy = (splus-sminus)/(xyeplus-xyeminus)
  dedy = (eplus-eminus)/(xyeplus-xyeminus)
  dmuhdy = (muhplus-muhminus)/(xyeplus-xyeminus)

  WRITE (*,*) 'RESULTS:'

  CALL PRINT_TEST( xn, xtemp, xye, x_n, x_p, x_a, x_l, x_h, &
                   A_h, Z_h, A_l, Z_l, press, entropy, eps, &
                   mu_n, mu_p, mu_hat, dmuhdn, dmuhdt, dmuhdy, &
                   dpdn, dpdt, dpdy, dsdn, dsdt, dsdy, dedn, dedt, dedy)

end program nse_test
