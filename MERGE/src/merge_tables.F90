!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SRO_EOS is distributed in the hope that it will be USEful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE MERGE_TABLES_MOD

  USE SNA_TABLE_MOD
  USE NSE_TABLE_MOD
  USE MERGE_TABLE_MOD
  USE Tables_Input_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE MERGE

    IMPLICIT NONE

    IF (.NOT.only_NSE .OR. only_SNA) THEN
      CALL SNA_EOS_TABLE
    ENDIF

    IF (.NOT.only_SNA .OR. only_NSE) THEN
      CALL NSE_EOS_TABLE
    ENDIF

!   obtain table for density larger than nt_max
    IF (.NOT.only_SNA .OR. .NOT. only_NSE) THEN
      CALL MERGE_EOS_TABLE
    ENDIF

  END SUBROUTINE MERGE

END MODULE MERGE_TABLES_MOD
