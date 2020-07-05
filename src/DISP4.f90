!
! =============================================================================
!
! DISP4
! Computational Systems Design Laboratory(CSDL)
! by Hyungmin Jun(hjun@jbnu.ac.kr)
!
!               s
!               |
!           2***|***1
!           *   |   *
!           *   ---------->r
!           *       *
!           3*******4
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

module DISP4

    use Data_Struct

    implicit none

    public Plane_Stiffness
    public Plane_Load
    public Plane_Stress

    private Material_Law
    private Gauss22_Point
    private Strain_Displacement
    private dHxy_Matrix
    private Shape_Function
    private dHrs_Matrix
    private Det

contains

! --------------------------------------------------------------------------------


! Stiffness matrix (u1, u2, u3, u4, v1, v2, v3, v4)
function Plane_Stiffness(Young, Poisson, thick, node) result(Ke)
    double precision, intent(in) :: Young, Poisson, thick
    double precision, intent(in) :: node(nNPE,nDIM)

    double precision :: Ke(nNPE*nDPN, nNPE*nDPN)
    double precision :: r, s, weight           ! Integration point, weight factor
    double precision :: det_j                  ! Determinant of Jacobian
    double precision :: B(3, nNPE*nDPN)        ! B-matrix (strain-displacement matrix)
    double precision :: C(3,3)                 ! Material matrix (material law)
    integer :: i, j

    ! Material matrix
    call Material_Law(Young, Poisson, C)

    Ke(:,:) = 0.0d0

    ! Numerical integration
    do i = 1, nGP
        do j = 1, nGP

            ! Gauss points and weight factor
            call Gauss22_Point(i, j, r, s, weight)

            ! Strain displacement matrix, B matrix
            call Strain_Displacement(r, s, node, det_j, B)

            ! Gauss integration
            Ke = Ke + weight * thick * matmul(matmul(transpose(B), C), B) * det_j
        end do
    end do
end function Plane_Stiffness

! --------------------------------------------------------------------------------

! Equivalent nodal loads
function Plane_Load(node, q) result(nodal_load)
    double precision, intent(in) :: node(nNPE, nDIM)    ! Node position
    double precision, intent(in) :: q(nDIM)             ! Body forces in x- and y-directions

    double precision :: nodal_load(nDPN*nNPE)   ! Element load vector 
    double precision :: H(nNPE)                 ! Shape functions
    double precision :: r, s, weight            ! Gauss point, weight factor
    double precision :: Jacob(2,2)              ! Jacobian matrix
    integer :: i, j, k

    nodal_load(:) = 0.0d0

    ! Numerical integration
    do i = 1, nGP
        do j = 1, nGP
            ! Gauss points
            call Gauss22_Point(i, j, r, s, weight)

            ! Shape function and Jacobian matrix
            H = Shape_Function(r, s)
            Jacob = Jacobian(r, s, node)

            ! Equivalent nodal load vector
            do k=1, nNPE
                nodal_load(k)      = nodal_load(k) + weight * dabs(Det(Jacob)) * H(k) * q(1)
                nodal_load(k+nNPE) = nodal_load(k+nNPE) + weight * dabs(Det(Jacob)) * H(k) * q(2)
            end do
        end do
    end do
end function Plane_Load

! --------------------------------------------------------------------------------

! Element stress (Sx, Sy, Sxy)
function Plane_Stress(young, poisson, node, disp) result(stress)
    double precision, intent(in) :: young, poisson      ! Young's modulus, Possion's ratio
    double precision, intent(in) :: node(nNPE,nDIM)     ! nodal position of element
    double precision, intent(in) :: disp(nDPN*nNPE)     ! displacement vector

    double precision :: stress(nNPE,3)      ! Stress (Sx, Sy, Sxy)
    double precision :: C(3,3), B(3, nDPN*nNPE)
    double precision :: buf_det
    double precision :: r, s, weight
    integer :: i, j

    ! Material matrix
    call Material_Law(Young, Poisson, C)

    do i = 1, nGP
        do j = 1, nGP
            ! Gauss points where stresses are out
            call Gauss22_Point(i, j, r, s, weight)

            ! B-matrix
            call Strain_Displacement(r, s, node, buf_det, B)

            ! Stresses
            stress(i*2+j-2, :) = matmul(matmul(C,B), disp)
        end do
    end do
end function Plane_Stress

! --------------------------------------------------------------------------------

! Material law for plane stress condition
subroutine Material_Law(Young, Poisson, C)
    double precision, intent(in)  :: Young, Poisson     ! Material constants
    double precision, intent(out) :: C(3,3)             ! Material matrix

    C(:,:) = 0.0d0
    C(1,1) = Young / (1.0d0 - Poisson**2)
    C(2,2) = C(1,1)
    C(1,2) = Poisson * C(1,1)
    C(2,1) = C(1,2)
    C(3,3) = 0.5d0 * Young / (1.0d0 + Poisson)
end subroutine Material_Law

! --------------------------------------------------------------------------------

! Gauss integration point 2 * 2 and weight factor
subroutine Gauss22_Point(i, j, r, s, weight)
    integer, intent(in) :: i, j
    double precision, intent(out) :: r, s, weight

    double precision :: GaussPoint(2), w(2)
    data GaussPoint / -0.577350269d0, 0.577350269d0 /
    data w          / 1.000000000d0, 1.000000000d0 /

    weight = w(i) * w(j)

    r = GaussPoint(i)
    s = GaussPoint(j)
end subroutine Gauss22_Point

! --------------------------------------------------------------------------------

! B-matrix (strain-displacement matrix)
subroutine Strain_Displacement(r, s, node, det_j, B)
    double precision, intent(in)  :: r, s, node(nNPE,nDIM)
    double precision, intent(out) :: det_j, B(3, nDPN*nNPE)

    integer:: i
    double precision :: dHxy(nDIM,nNPE)

    ! calculate dHxy. dHxy(1,:)=dH/dx, dH(2,:)=dHxy/dy
    call dHxy_Matrix(r, s, node, det_j, dHxy)

    ! B-matrix
    B(:,:) = 0.0d0

    do i=1, nNPE
        B(1, i)        = dHxy(1, i)
        B(2, i + nNPE) = dHxy(2, i)
        B(3, i)        = dHxy(2, i)
        B(3, i + NNPE) = dHxy(1, i)
    end do
end subroutine Strain_Displacement

! --------------------------------------------------------------------------------

! dHxy matirx. dHxy(1,:)=dH/dx, dHxy(2,:)=dH/dy
subroutine dHxy_Matrix(r, s, node, det_j, dHxy)
    double precision, intent(in)  :: r, s, node(nNPE,nDIM)
    double precision, intent(out) :: det_j, dHxy(nDIM,nNPE) 

    double precision :: buf, dHrs(nDIM,nNPE), Jacob(2,2)

    ! Derivative of shape function, dH/dr and dH/ds
    dHrs = dHrs_Matrix(r, s)

    ! Jacob = Jacobian Matrix
    Jacob = Jacobian(r, s, node)

    ! Jacob => inverse of jacobian Matrix
    det_j      = Jacob(1,1) * Jacob(2,2) - Jacob(1,2) * Jacob(2,1)
    Jacob(1,2) = -Jacob(1,2)
    Jacob(2,1) = -Jacob(2,1)
    buf        = Jacob(1,1)
    Jacob(1,1) = Jacob(2,2)
    Jacob(2,2) = buf
    Jacob      = Jacob / det_j
 
    ! dHxy(1,:)=dH/dx, dHxy(2,:)=dH/dy
    dHxy = matmul(Jacob, dHrs)
end subroutine dHxy_Matrix

! --------------------------------------------------------------------------------

! H matrix, determinant of Jacobian
function Jacobian(r, s, node) result(Jacob)
    double precision, intent(in) :: r, s
    double precision, intent(in) :: node(nNPE,nDIM)

    double precision :: dHrs(nDIM,nNPE), Jacob(2,2)

    ! Derivative of shape functions, dH/dr and dH/ds
    dHrs = dHrs_Matrix(r, s)

    ! Jacobian matrix
    Jacob = matmul(dHrs, node)
end function Jacobian

! --------------------------------------------------------------------------------

! Derivatives of shape functions, dHrs(1,:)=dH/dr and dHrs(2,:)=dH/ds
function dHrs_Matrix(r, s) result (dHrs)
    double precision, intent(in) :: r, s
    double precision :: dHrs(nDIM,nNPE)

    dHrs(1,1) =  0.25d0 * (1.0d0 + s)
    dHrs(1,2) = -0.25d0 * (1.0d0 + s)
    dHrs(1,3) = -0.25d0 * (1.0d0 - s)
    dHrs(1,4) =  0.25d0 * (1.0d0 - s)

    dHrs(2,1) =  0.25d0 * (1.0d0 + r)
    dHrs(2,2) =  0.25d0 * (1.0d0 - r)
    dHrs(2,3) = -0.25d0 * (1.0d0 - r)
    dHrs(2,4) = -0.25d0 * (1.0d0 + r)
end function dHrs_Matrix

! --------------------------------------------------------------------------------

! Shape functions
function Shape_Function(r, s) result(H)
    double precision, intent(in) :: r, s   ! Natural coordinate

    double precision :: H(nNPE)

    H(1) = 0.25d0 * (1.0d0 + r) * (1.0d0 + s)
    H(2) = 0.25d0 * (1.0d0 - r) * (1.0d0 + s)
    H(3) = 0.25d0 * (1.0d0 - r) * (1.0d0 - s)
    H(4) = 0.25d0 * (1.0d0 + r) * (1.0d0 - s)
end function Shape_Function

! --------------------------------------------------------------------------------

! Determinant of 2 by 2 matrix
function Det(Jacob) result(det_j)
    double precision, intent(in) :: Jacob(2,2)

    double precision :: det_j

    ! Determinant
    det_j = Jacob(1,1) * Jacob(2,2) - Jacob(1,2) * Jacob(2,1)
end function Det

! --------------------------------------------------------------------------------

end Module DISP4