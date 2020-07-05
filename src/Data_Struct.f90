!
! =============================================================================
!
! Data_Struct
! Computational Systems Design Laboratory(CSDL)
! by Hyungmin Jun(hjun@jbnu.ac.kr)
!
! =============================================================================
!
! Copyright 2020 CSDL. All rights reserved.
!
! License - GPL version 3
! This program is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
module Data_Struct

    integer, parameter :: nNPE = 4      ! # of nodes per element
    integer, parameter :: nDPN = 2      ! # of DOFs per node
    integer, parameter :: nDIM = 2      ! Problem dimension
    integer, parameter :: nGP  = 2      ! # of Gauss points

    ! NodeType structure
    type :: NodeType
        double precision :: x(nDIM)     ! Nodal position (x, y)
        double precision :: pm(nDPN)    ! Nodal force (Px, Py)
        integer          :: bc(nDPN)    ! Displacement BC (u, v) (1=fixed, 0=free)
        integer          :: eq_n(nDPN)  ! Equation number (u, v)
    end type NodeType

    ! ElemType structure
    type :: ElemType
        integer          :: cn(nNPE)    ! Connectivity
        double precision :: q(nDIM)     ! Distributed load in x- and y-directions
    end type ElemType

    ! PropType structure
    type :: PropType
        double precision :: young       ! Young's modulus
        double precision :: poisson     ! Poison's ratio
        double precision :: thick       ! Thickness
    end type PropType

! ---------------------------------------------------------------------------------------

end module Data_Struct