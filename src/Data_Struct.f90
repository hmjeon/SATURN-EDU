
! --------------------------------------------------------------------------------
!            Module for defining user-defined types
!            - created by Phill-Seung Lee, 19/Feb/1998 
!            - FORTRAN 95
! --------------------------------------------------------------------------------

module Data_Struct

    integer, parameter :: nNPE = 4      ! # of nodes per element
    integer, parameter :: nDPN = 2      ! # of DOFs per node
    integer, parameter :: nDIM = 2      ! Problem dimension
    integer, parameter :: nGP  = 2      ! # of Gauss points


    ! define a type for node
    type :: NodeType
        double precision :: x(2)       ! nodal position (x, y), nDIM
        double precision :: pm(2)      ! nodal force (Px, Py), nDPN
        integer          :: bc(2)      ! displacement BC (u, v) (1=fixed, 0=free), nDPN
        integer          :: eq_n(2)    ! equation number (u, v), nDPN
    end type NodeType

    ! define a type for element
    type :: ElementType
        integer          :: cn(4)      ! connectivity, nNPE
        double precision :: thickness  ! thickness
        double precision :: q(2)       ! distributed load in x- and y-directions, nDIM
    end type ElementType

    ! define a type for material
    type :: MaterialType
        double precision :: Young      ! Young's modulus
        double precision :: Poisson    ! Poison's ratio
    end type MaterialType

end module Data_Struct