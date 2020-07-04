!
! =============================================================================
!
! DISP4
! by Hyungmin Jun (hjun@jbnu.ac.kr)
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
! Copyright 2019 Hyungmin Jun. All rights reserved.
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

    public Plane_Stiffness, Plane_Stress, Plane_Load

    private Gauss22_Point, Shape_Function ,dHrs_Matrix
    private Strain_displacement, dHxy_Matrix, Material_law

contains

! --------------------------------------------------------------------------------

! Stiffness matrix (u1, u2, u3, u4, v1, v2, v3, v4)
subroutine Plane_Stiffness(Young, Poisson, thick, node, Ke)
    double precision, intent(in)  :: Young, Poisson, thick
    double precision, intent(in)  :: node(nNPE,nDIM)
    double precision, intent(out) :: Ke(nNPE*nDPN, nNPE*nDPN)
  
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
            call Strain_displacement(r, s, node, det_j, B)

            ! Gauss integration
            Ke = Ke + weight * thick * matmul(matmul(transpose(B), C), B) * det_j
        end do
    end do
end subroutine Plane_Stiffness

! --------------------------------------------------------------------------------

! calculate stress (Sx, Sy, Sxy)
subroutine Plane_Stress(Young, Poisson, Node, Displace, Stress)
    double precision, intent(in)  :: Young, Poisson         ! Young's modulus, Possion's ratio
    double precision, intent(in)  :: node(nNPE,nDIM)        ! nodal position of element
    double precision, intent(in)  :: displace(nDPN * nNPE)  ! displacement vector
    double precision, intent(out) :: Stress(nNPE,3)         ! stress (Sx, Sy, Sxy)
    double precision :: C(3,3), B(3,nDPN * nNPE)
    double precision :: buf_det
    double precision :: r, s, weight
    integer :: i, j
 
    ! calculate a material matrix
    call Material_law(Young, Poisson, C)
   
    do i=1, 2
        do j=1, 2
            ! set postions where stresses are out
            call Gauss22_Point(i, j, r, s, weight)              ! at Gauss points 
         
            ! calculate B-matrix
            call strain_displacement(r, s, node, buf_det, B)
         
            ! calculate stresses
            Stress(i * 2 + j - 2, :) = matmul( matmul(C,B), displace )
        end do
    end do
end subroutine Plane_Stress

! --------------------------------------------------------------------------------

! equivalent nodal loads
subroutine Plane_Load(node, q, nodal_load)
    double precision, intent(in)  :: node(nNPE,nDIM)          ! node position
    double precision, intent(in)  :: q(nDIM)                  ! body forces in x- and y-directions
    double precision, intent(out) :: nodal_load(nDPN * nNPE)  ! element load vector
 
    double precision :: H(nNPE)                     ! shape functions
    double precision :: r, s, weight                ! Gauss point, weight factor
    
    double precision :: Jacob(2,2)                       ! determinant of Jacobian
    integer           :: i, j, k

    ! numerical integration
    nodal_load(:) = 0.0d0

    do i=1, 2
        do j=1, 2
            ! Gauss points
            call Gauss22_Point(i, j, r, s, weight)      
         
            ! determinant of Jacobian and shape function (H)
            H = Shape_Function(r, s)
            Jacob = Jacobian(r, s, node)
         
            ! equivalent nodal load vector
            do k=1, nNPE
                nodal_load(k)      = nodal_load(k) + weight * dabs(Det(Jacob)) * H(k) * q(1)
                nodal_load(k+nNPE) = nodal_load(k+nNPE) + weight * dabs(Det(Jacob)) * H(k) * q(2)
            end do
        end do
    end do
end subroutine Plane_Load

! --------------------------------------------------------------------------------

! Gauss integration point 2 * 2 and weight factor
subroutine Gauss22_Point(i, j, r, s, weight)
    integer, intent(in)             :: i, j
    double precision, intent(out)  :: r, s, weight

    double precision :: GaussPoint(2), w(2)
    data GaussPoint / -0.577350269d0, 0.577350269d0 /
    data w          / 1.000000000d0, 1.000000000d0 /

    weight = w(i) * w(j)

    r = GaussPoint(i)
    s = GaussPoint(j)
end subroutine Gauss22_Point

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

! Derivatives of shape functions, dHrs(1,:)=dH/dr and dHrs(2,:)=dH/ds
function dHrs_Matrix(r, s) result (dHrs)
    double precision, intent(in)  :: r, s
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

! H matrix, determinant of Jacobian
function Jacobian(r, s, node) result(Jacob)
    double precision, intent(in)  :: r, s
    double precision, intent(in)  :: node(nNPE,nDIM)

    double precision :: dHrs(nDIM,nNPE), Jacob(2,2)

    ! Derivative of shape functions, dH/dr and dH/ds
    dHrs = dHrs_Matrix(r, s)

    ! Jacobian matrix
    Jacob = matmul(dHrs, node)
end function Jacobian

! --------------------------------------------------------------------------------

! H matrix, determinant of Jacobian
function Det(Jacob) result(det_j)
    double precision, intent(in)  :: Jacob(2,2)

    double precision :: det_j

    ! Determinant of jacobian matrix
    det_j = Jacob(1,1) * Jacob(2,2) - Jacob(1,2) * Jacob(2,1)
end function Det

! --------------------------------------------------------------------------------

! dHxy matirx. dHxy(1,:)=dH/dx, dHxy(2,:)=dH/dy
subroutine dHxy_Matrix(r, s, node, det_j, dHxy)
    double precision, intent(in)  :: r, s, Node(nNPE,nDIM)
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

! B-matrix (strain-displacement matrix)
subroutine Strain_displacement(r, s, node, det_j, B)
    double precision, intent(in)  :: r, s, node(nNPE,nDIM)
    double precision, intent(out) :: det_j, B(3,nDPN * nNPE)
    integer:: i
    double precision :: dHxy(nDIM,nNPE)

    ! calculate dHxy. dHxy(1,:)=dH/dx, dH(2,:)=dHxy/dy
    call dHxy_Matrix(r, s, node, det_j, dHxy)

    ! B-matrix
    B(:,:)   = 0.0d0

    do i=1, nNPE
        B(1,i)          = dHxy(1,i)
        B(2,i + nNPE)   = dHxy(2,i)
        B(3,i)          = dHxy(2,i)
        B(3,i + NNPE)   = dHxy(1,i)
    end do

end subroutine Strain_displacement

! --------------------------------------------------------------------------------

! material law for plane stress condition
subroutine Material_law(Young, Poisson, C)
    double precision, intent(in)  :: Young, Poisson  ! material constants
    double precision, intent(out) :: C(3,3)          ! material matrix

    C(:,:) = 0.0d0
    C(1,1) = Young / (1.0d0 - Poisson**2)
    C(2,2) = C(1,1)
    C(1,2) = Poisson * C(1,1)
    C(2,1) = C(1,2)
    C(3,3) = 0.5d0 * Young / (1.0d0 + Poisson)
end subroutine Material_law

! --------------------------------------------------------------------------------

end Module DISP4  ! --------------- End of Module --------------------------
