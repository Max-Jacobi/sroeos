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
MODULE READ_NUCLEAR_DATA_TABLE_MOD

  USE Kind_types_Mod, ONLY : I4B, DP
  USE Input_Files_Mod

  IMPLICIT NONE

! define nucleus type
  TYPE nucleus
    CHARACTER(len=5)    :: name
    REAL(DP)              :: A
    REAL(DP)              :: Z
    REAL(DP)              :: N
    REAL(DP)              :: spin    ! ground spin: 2*J + 1
    REAL(DP)              :: mass    ! mass of nuclei in MeV
    REAL(DP)              :: be      ! total binding energy
    REAL(DP), ALLOCATABLE :: PartFun(:)
    LOGICAL               :: have_part_fun
  END TYPE nucleus

  TYPE(nucleus), DIMENSION(:), ALLOCATABLE, SAVE :: nuclei
  REAL(DP), DIMENSION(:), ALLOCATABLE :: PartT9
  INTEGER(I4B), PARAMETER :: NPartT9 = 24

CONTAINS

  SUBROUTINE READ_NUCLEAR_DATA ()

    USE Physical_Constants_Mod, ONLY: proton_mass_MeV, neutron_mass_MeV
    USE Isotopes_Mod
    USE Test_Input_Mod, ONLY : IS_TEST
    USE Make_Tables_Mod, ONLY : output_directory

    IMPLICIT NONE

! internal vars
    INTEGER(I4B) :: i,j,k,l,niso
    INTEGER(I4B), PARAMETER :: maxnuc = 7753
    TYPE(nucleus), ALLOCATABLE :: tmpnuc(:)
    CHARACTER(len=256) :: sbuf, command
    CHARACTER(len=5),ALLOCATABLE :: isoname(:)
    REAL(DP), ALLOCATABLE :: aion(:),zion(:),nion(:),dm(:),be(:)
    REAL(DP), ALLOCATABLE :: jion(:),tmppart(:,:)

! first READ isotope list and check that we have the
! right number of isotopes
! TODO: CHANGE HOW LIST IS READ SO CAN USE MAXNUC == NUMBER OF NUCLEONS

    ALLOCATE(tmpnuc(maxnuc))

    OPEN(666,FILE=TRIM(ADJUSTL(input_iso_list)),STATUS='old')
    IF (IS_TEST) THEN
      OPEN(667,FILE='TEST_ISOTOPES.DAT')
    ELSE
      WRITE(command,*) 'cp ', TRIM(ADJUSTL(input_iso_list)), '  ', &
                              TRIM(ADJUSTL(output_directory))
    ENDIF

    i=1

    WRITE (*,*)
    WRITE (*,*) 'Reading isotopes in the input list'
    WRITE (*,*)

    DO
      READ(666,*,END=10) tmpnuc(i)%name,tmpnuc(i)%A,tmpnuc(i)%Z,tmpnuc(i)%N
      IF (IS_TEST) WRITE(667,"(1A6,3F7.2)") &
                         tmpnuc(i)%name,tmpnuc(i)%A,tmpnuc(i)%Z,tmpnuc(i)%N
      ! check if N + Z = A
      IF (tmpnuc(i)%A-tmpnuc(i)%Z-tmpnuc(i)%N>1.d-10) THEN
        WRITE (*,*) 'Error in list isotope', tmpnuc(i)%name
        STOP
      ENDIF
      i = i+1
    ENDDO

10  CLOSE(666)
    CLOSE(667)

    IF (IS_TEST) THEN
      WRITE(*,*)
      WRITE(*,*) 'Isotope list written to TEST_ISOTOPES.DAT file.'
      WRITE(*,*)
    ELSE
      WRITE(*,*)
      WRITE(*,*) 'Isotope list list.in copied to folder ', &
                  TRIM(ADJUSTL(output_directory))
      WRITE(*,*)
      WRITE(command,*) 'cp ', TRIM(ADJUSTL(input_iso_list)), '  ', &
                              TRIM(ADJUSTL(output_directory))
    ENDIF

    WRITE (*,*)
    WRITE (*,*) 'Isotope list contains ', i-1, ' isotopes.'
    WRITE (*,*)

    niso = i-1
    num_isotopes = niso
    ! check below should never be an issue!!
    IF (num_isotopes.ne.niso) THEN
      WRITE(6,*) "READ list problem: num_isotopes != number of isotopes"
      STOP
    ELSE
      IF (ALLOCATED(nuclei)) DEALLOCATE (nuclei)
      ALLOCATE(nuclei(num_isotopes))
      nuclei(1:num_isotopes) = tmpnuc(1:num_isotopes)
    ENDIF
    DEALLOCATE(tmpnuc)
!
!   check for light isotopes (Z<=6)
    llmin = 0 ; llmax = 0
    DO l=1,num_isotopes
      IF (nuclei(l)%Z<=6.d0 .and. llmin==0) llmin = l
      IF (nuclei(l)%Z<=6.d0 .and. llmin/=0) llmax = l
    ENDDO

    WRITE (*,*) 'Light Isotopes considered'
    WRITE (*,*) 'Name  Z  A  ', llmin, llmax
    DO l = llmin, llmax
      WRITE (*,*) nuclei(l)%name,nuclei(l)%Z,nuclei(l)%A
    ENDDO
    WRITE (*,*)
    WRITE (*,*) 'proton  number fraction stored in xp parameter'
    WRITE (*,*) 'neutron number fraction stored in xn parameter'
    WRITE (*,*) 'alpha particle number fraction stored in xa parameter'
    WRITE (*,*) 'other light isotopes number fraction stored in xl parameter'
    WRITE (*,*) 'heavy nuclei number fraction stored in xh parameter'
    WRITE (*,*)

! now READ masses of isotopes in the list
! and figure out how many entries in mass table
    OPEN(666,file=TRIM(ADJUSTL(input_iso_properties)),status='old')
    i=1
    DO
      READ(666,*,END=20) sbuf
      i=i+1
    ENDDO

20  CLOSE(666)
    niso = i-1

!   allocate tables to store isotopes properties
    ALLOCATE(aion(niso))
    ALLOCATE(zion(niso))
    ALLOCATE(nion(niso))
    ALLOCATE(dm(niso))
    ALLOCATE(be(niso))
    ALLOCATE(isoname(niso))

!   READ isotope masses
    OPEN(666,file=TRIM(ADJUSTL(input_iso_properties)),status='old')
    DO i=1,niso
       READ(666,*) isoname(i),aion(i),zion(i),nion(i),dm(i),be(i)
    ENDDO

!03  format(a5,3F5.0,F16.5,F14.5)
    ! now add to nuclei array
    k = 0

    WRITE (*,*) 'ISOTOPES USED'
    WRITE (*,*) ' name    mass (MeV) '
    WRITE (*,*) 
    TARGET: DO i=1,niso
      LIST: DO j=1,num_isotopes
        IF(TRIM(ADJUSTL(nuclei(j)%name)) .eq. TRIM(ADJUSTL(isoname(i)))) THEN
          nuclei(j)%be = be(i)
          nuclei(j)%mass = (nuclei(j)%A-nuclei(j)%Z) * neutron_mass_MeV &
                          + nuclei(j)%Z * proton_mass_MeV - nuclei(j)%be
          write (*,*) nuclei(j)%name, nuclei(j)%mass
          k = k +1
          CYCLE TARGET
        ENDIF
      ENDDO LIST
    ENDDO TARGET

    DEALLOCATE(aion)
    DEALLOCATE(zion)
    DEALLOCATE(nion)
    DEALLOCATE(dm)
    DEALLOCATE(be)
    DEALLOCATE(isoname)

    IF(k.ne.num_isotopes) THEN
       WRITE(6,*) "Fatal Error!"
       WRITE(6,*) k,num_isotopes
       WRITE(6,*) "Did not find nuclear data for all requested isotopes! :-("
       STOP
    ENDIF

!   now READ partition function information
!   for now use hard coded number of T9 entries
!   and ALLOCATE the arrays inside the nuclei structure
    DO i=1,num_isotopes
       ALLOCATE(nuclei(i)%PartFun(NPartT9))
       nuclei(i)%have_part_fun = .false.
    ENDDO

!   hard coded temperatures in partition table
    ALLOCATE(PartT9(NPartT9))
    PartT9 = (/0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7, &
         &      0.8,0.9,1.0,1.5,2.0,2.5,3.0,3.5, &
         &      4.0,4.5,5.0,6.0,7.0,8.0,9.0,10.0/)

!   figure out how many entries in partition function table
    OPEN(666,file=TRIM(ADJUSTL(input_partition)),status='old')

    i=1

    DO
      ! assume a particular format for the partition table
      ! -- generally not a good idea :-(
      READ(666,*,END=40) sbuf
      READ(666,*,END=40) sbuf
      READ(666,*,END=40) sbuf
      READ(666,*,END=40) sbuf
      READ(666,*,END=40) sbuf
      i=i+1
    ENDDO

40  CLOSE(666)

    niso = i-1

    ALLOCATE(aion(niso))
    ALLOCATE(zion(niso))
    ALLOCATE(jion(niso))
    ALLOCATE(tmppart(niso,NPartT9))
    ALLOCATE(isoname(niso))
    OPEN(666,file=TRIM(ADJUSTL(input_partition)),status='old')

    DO i=1,niso
       READ(666,'(a5)') isoname(i)
       READ(666,*) zion(i),aion(i),jion(i)
       READ(666,*) (tmppart(i,j),j=1,NPartT9)
    ENDDO
    CLOSE(666)

    ! now sort this in
    k=0

    TARGET2: DO i=1,niso
      LIST2: DO j=1,num_isotopes
        IF (nuclei(j)%Z .eq. zion(i) .AND. nuclei(j)%A .eq. aion(i) ) THEN
          nuclei(j)%spin = 2.0d0*jion(i) + 1.0d0
          nuclei(j)%PartFun(:) = tmppart(i,:)
          nuclei(j)%have_part_fun = .true.
          k = k+1
          CYCLE TARGET2
        ENDIF
      ENDDO LIST2
    ENDDO TARGET2

    DEALLOCATE(aion)
    DEALLOCATE(zion)
    DEALLOCATE(jion)
    DEALLOCATE(tmppart)
    DEALLOCATE(isoname)

    IF (k.ne.num_isotopes) THEN
      WRITE(6,*) "Fatal Error!"
      WRITE(6,*) k,num_isotopes
      WRITE(6,*) "Did not find partition function data for all requested isotopes! :-("
      DO i=1,num_isotopes
        IF (.NOT.nuclei(i)%have_part_fun) then
          WRITE(6,*) nuclei(i)%name
          nuclei(i)%spin = 1.0d0
          nuclei(i)%PartFun(:) = 1.0d0
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE READ_NUCLEAR_DATA

  SUBROUTINE GETPART (num_isotopes,xtemp,nuclei,wion,DwionDT)

    USE Physical_Constants_Mod, ONLY: temp_MeV_to_Kelvin

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN)    :: num_isotopes
    REAL(DP), INTENT(IN)        :: xtemp
    TYPE(nucleus), INTENT(IN)   :: nuclei(num_isotopes)
    REAL(DP), INTENT(out)       :: wion(num_isotopes)
    REAL(DP), INTENT(out)       :: DwionDT(num_isotopes)
    ! local vars
    REAL(DP) :: T9,T9max
    REAL(DP) :: x0,x1,p0,p1,m0,m1,pl,pu,xl,xu

    REAL(DP) :: dmid, dup, dlo, xx
    REAL(DP) :: h00, h10, h01, h11
    REAL(DP) ::  dh00, dh10, dh01, dh11
    REAL(DP) ::  ddh00, ddh10, ddh01, ddh11

    INTEGER(I4B) :: i,ill,ilow,iup,iuu,s9

!   convert MeV to 10^9 Kelvin
    T9 = xtemp * temp_MeV_to_Kelvin * 1.0d-9
!   Set a maximum value to compute partition functions to avoid over/underflow
!    Partition functions at T>T9Max have the same values as they do at T9max
    T9max = 10.d0**2.d0
    s9 = size(PartT9)
    IF (T9 > T9max) T9 = T9max
    IF (T9 >= PartT9(s9)) THEN
      ! Perform a linear extrapolation in log space
      x0 = partt9(s9-1)
      x1 = partt9(s9)
      DO i=1, num_isotopes
        p0 = LOG(nuclei(i)%partfun(s9-1))
        p1 = LOG(nuclei(i)%partfun(s9))
        m0 = (p1-p0)/(x1-x0)
        wion(i)    = p1 + m0*(T9-x1)
        DwionDT(i) = m0
        wion(i) =  nuclei(i)%spin * EXP(wion(i))
        DwionDT(i) = wion(i) * DwionDT(i) * temp_MeV_to_Kelvin * 1.0d-9
        if(T9>=T9max) DwionDT(i) = 0.d0
      ENDDO
      RETURN
!   If Temperature is smaller than minimum in partition function table
!    set Temperature to its value at Tmin = PartT9(1)
    ELSEIF (T9 <= PartT9(1)) THEN
      DO i=1, num_isotopes
        wion(i) = nuclei(i)%spin * nuclei(i)%partfun(1)
        DwionDt(i) = 0.d0
      ENDDO
      RETURN
    ENDIF

!   If T9min < T < T9max then use interpolation scheme
!    to find partition function at T

    ilow = findlower(PartT9,T9)

    ilow = MIN(SIZE(PartT9)-1, ilow)
    iup = ilow+1
    iuu = MIN(iup+1, SIZE(PartT9))
    ill = MAX(ilow-1, 1)

!    if (ilow == 1) ill = iup
! Make derivative at the lower bound of the table zero

    x0 = PartT9(ilow)
    x1 = PartT9(iup)
    xl = PartT9(ill)
    xu = PartT9(iuu)
    xx = (T9 - x0)/(x1 - x0)
    dmid = x1 - x0

    h00 = 2.d0*xx**3 - 3.d0*xx**2 + 1.d0
    h10 = xx**3 - 2.d0*xx**2 + xx
    h01 = -2.d0*xx**3 + 3.d0*xx**2
    h11 = xx**3 - xx**2

    dh00 = 6.d0*xx**2 - 6.d0*xx
    dh10 = 3.d0*xx**2 - 4.d0*xx + 1.d0
    dh01 =-6.d0*xx**2 + 6.d0*xx
    dh11 = 3.d0*xx**2 - 2.d0*xx

    ddh00 = 12.d0*xx - 6.d0
    ddh10 = 6.d0*xx - 4.d0
    ddh01 =-12.d0*xx + 6.d0
    ddh11 = 6.d0*xx - 2.d0

    DO i=1, num_isotopes
      p0 = LOG(nuclei(i)%partfun(ilow))
      p1 = LOG(nuclei(i)%partfun(iup))
      pl = LOG(nuclei(i)%partfun(ill))
      pu = LOG(nuclei(i)%partfun(iuu))
      m0 = (p1-pl)/(x1-xl) !0.5d0*((p1-p0)/(x1-x0)+(p0-pl)/(x0-xl))
      m1 = (pu-p0)/(xu-x0) !0.5d0*((pu-p1)/(xu-x1)+(p1-p0)/(x1-x0))
      IF (x0==xl) m0 = 0.d0
      IF (x1==xu) m1 = (p1-p0)/(x1-x0)
      wion(i)    = p0* h00 + m0*dmid* h10 + p1* h01 + m1*dmid* h11
      DwionDT(i) = p0*dh00 + m0*dmid*dh10 + p1*dh01 + m1*dmid*dh11
      wion(i) = nuclei(i)%spin * exp(wion(i))
      DwionDt(i) = wion(i) * DwionDt(i) / dmid * temp_MeV_to_Kelvin * 1.0d-9
    ENDDO

CONTAINS

    FUNCTION findlower(xx,x) RESULT (lower)
      ! find the position of x in xx

      IMPLICIT NONE

      INTEGER(I4B) :: lower,n,jl,jn,ju
      REAL(DP), INTENT(IN) :: x
      REAL(DP), dimension(:), INTENT(IN) :: xx
      LOGICAL :: ascnd

      n = SIZE(xx)
      ascnd = (xx(n) >= xx(1))
      jl = 0
      ju = n+1
      DO
         IF(ju-jl <= 1) EXIT
         jn = (ju+jl)/2 !.0d0
         IF (ascnd .eqv. (x >= xx(jn))) THEN
            jl = jn
         ELSE
            ju = jn
         ENDIF
      ENDDO
      if (x == xx(1)) then
         lower = 1
      else if (x == xx(n)) then
         lower = n-1
      else
         lower = jl
      ENDIF
      lower = MIN(SIZE(xx), MAX(1, lower))
    END FUNCTION findlower

  END SUBROUTINE GETPART

END MODULE READ_NUCLEAR_DATA_TABLE_MOD
