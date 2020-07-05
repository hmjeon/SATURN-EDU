!
! =============================================================================
!
! Math_Kernel
! Last Updated : 07/22/2019, by Hyungmin Jun (hyungminjun@outlook.com)
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

module Math_Kernel

    use Data_Struct

    implicit none

    public Int2Str
    public Dble2Str
    public ES2Str
    public E2Str
    public Rad2Deg
    public Deg2Rad
    public Init_Node
    public Init_Elem

    public Material_Law
    public Gauss3_Point
    public Gauss3_Point_xy
    public Gauss7_Point
    public Gauss7_Point_xy
    public Director
    public Cross_Product
    public Normalize
    public Size_Vector
    public Determinant_33
    public Determinant
    public Det
    public Inverse_22
    public Inverse_33
    public Inverse_LU
    public InvMat
    public Find_Diff
    public Equal_Vector
    public Find_Minimum
    public Swap
    public Sort
    public Sort2
    public Cross
    public IsPlanar
    public GaussianMask
    public AxialRotate
    public Rotate_Vector
    public LeastSquare

    private Triangularization

    double precision, parameter :: pi  = 3.141592653589793d0
    double precision, parameter :: eps = 0.0000001d0

contains

! ---------------------------------------------------------------------------------------

! Integer to string
function Int2Str(int) result(str)
    integer, intent(in) :: int

    character(100) :: str

    write(str, "(i)"), int
end function Int2Str

! ---------------------------------------------------------------------------------------

! double precision to string
function Dble2Str(dbl, opt1, opt2) result(str)
    double precision, intent(in) :: dbl
    integer, intent(in) :: opt1, opt2

    character(len=100) :: str, form

    form = "(f"//trim(adjustl(Int2Str(opt1)))//"."//trim(adjustl(Int2Str(opt2)))//")"
    !write(str, "(f<opt1>.<opt2>)"), dbl
    write(str, trim(form)), dbl
end function Dble2Str

! ---------------------------------------------------------------------------------------

! Scientific format to string
function ES2Str(dbl, opt1, opt2) result(str)
    double precision, intent(in) :: dbl
    integer, intent(in) :: opt1, opt2

    character(len=100) :: str, form

    form = "(es"//trim(adjustl(Int2Str(opt1)))//"."//trim(adjustl(Int2Str(opt2)))//")"
    !write(str, "(es<opt1>.<opt2>)"), dbl
    write(str, trim(form)), dbl
end function ES2Str

! ---------------------------------------------------------------------------------------

! Scientific format to string
function E2Str(dbl, opt1, opt2) result(str)
    double precision, intent(in) :: dbl
    integer, intent(in) :: opt1, opt2

    character(len=100) :: str, form

    form = "(e"//trim(adjustl(Int2Str(opt1)))//"."//trim(adjustl(Int2Str(opt2)))//")"
    !write(str, "(e<opt1>.<opt2>)"), dbl
    write(str, trim(form)), dbl
end function E2Str

! ---------------------------------------------------------------------------------------

! Radian to degree
! Last updated on Friday 29 Apr 2016 by Hyungmin
function Rad2Deg(rad) result(deg)
    double precision, intent(in) :: rad
    
    double precision :: deg
    deg = rad * 180.0d0 / pi
end function Rad2Deg

! ---------------------------------------------------------------------------------------

! Degree to radian
! Last updated on Monday 2 May 2016 by Hyungmin
function Deg2Rad(deg) result(rad)
    double precision, intent(in) :: deg

    double precision :: rad
    rad = deg * pi / 180.0d0
end function Deg2Rad

! ---------------------------------------------------------------------------------------

! Initialize NodeType structure
subroutine Init_Node(node)
    type (NodeType), intent(inout) :: node(:)

    integer :: i

    do i = 1, size(node)
        node(i)%x(:)    = 0.0d0
        node(i)%pm(:)   = 0.0d0
        node(i)%bc(:)   = 0
        node(i)%eq_n(:) = 0
    end do
end subroutine Init_Node

! ---------------------------------------------------------------------------------------

! Initialize ElemType structure
subroutine Init_Elem(elem)
    type(ElemType), intent(inout) :: elem(:)

    integer :: i

    do i = 1, size(elem)
        elem(i)%q(:)  = 0.0d0
        elem(i)%cn(:) = 0
    end do
end subroutine Init_Elem

! ---------------------------------------------------------------------------------------

! Material law of local coordinate
subroutine Material_Law(Young, Poisson, C)
    double precision, intent(in)  :: Young, Poisson
    double precision, intent(out) :: C(6,6)

    ! Shear correction factor
    double precision :: k

    k      = 1.0d0
    C(:,:) = 0.0d0
    C(1,1) = Young / (1.0d0 - Poisson**2.0d0)
    C(2,2) = C(1,1)
    C(1,2) = Poisson * C(1,1)
    C(2,1) = C(1,2)
    C(4,4) = Young / (1.0d0 + Poisson) / 2.0d0
    C(5,5) = k * C(4,4)
    C(6,6) = C(5,5)
end subroutine Material_Law

! ---------------------------------------------------------------------------------------

! Gauss integration point 3 and weight factor
subroutine Gauss3_Point(i, k, r, s, t, w)
    integer, intent(in) :: i, k
    double precision, intent(out) :: r, s, t, w

    double precision :: GaussPoint2(2)
    data GaussPoint2 / -0.577350269189626d0, 0.577350269189626d0 /
    
    select case (i)
        case (1); r = 0.166666666666667d0; s = 0.166666666666667d0
        case (2); r = 0.666666666666667d0; s = 0.166666666666667d0
        case (3); r = 0.166666666666667d0; s = 0.666666666666667d0
    end select

    t = GaussPoint2(k)
    w = 0.333333333333333d0
    w = 0.5d0 * w
end subroutine Gauss3_Point

! ---------------------------------------------------------------------------------------

subroutine Gauss3_Point_xy(i, x, y, w)
    integer, intent(in) :: i
    double precision, intent(out) :: x, y, w

    w  = 0.333333333333333d0
    select case (i)
        case (1) 
            x = 0.166666666666667d0
            y = 0.166666666666667d0
        case (2)
            x = 0.666666666666667d0
            y = 0.166666666666667d0
        case (3)
            x = 0.166666666666667d0
            y = 0.666666666666667d0
    end select

    w = 0.5d0*w
end subroutine Gauss3_Point_xy

! ---------------------------------------------------------------------------------------

! Gauss integration point 7 and weight factor
subroutine Gauss7_Point(i, k, r, s, t, w)
    integer, intent(in) :: i, k
    double precision, intent(out) :: r, s, t, w

    double precision :: GaussPoint2(2)
    data GaussPoint2 / -0.577350269189626d0, 0.577350269189626d0 /

    select case (i)
        case (1); r = 0.1012865073235d0; s = 0.1012865073235d0; w = 0.1259391805448d0
        case (2); r = 0.7974269853531d0; s = 0.1012865073235d0; w = 0.1259391805448d0
        case (3); r = 0.1012865073235d0; s = 0.7974269853531d0; w = 0.1259391805448d0
        case (4); r = 0.4701420641051d0; s = 0.0597158717898d0; w = 0.1323941527885d0
        case (5); r = 0.4701420641051d0; s = 0.4701420641051d0; w = 0.1323941527885d0
        case (6); r = 0.0597158717898d0; s = 0.4701420641051d0; w = 0.1323941527885d0
        case (7); r = 0.3333333333333d0; s = 0.3333333333333d0; w = 0.2250000000000d0
    end select

    t = GaussPoint2(k)
    w = 0.5d0 * w
end subroutine Gauss7_Point

! ---------------------------------------------------------------------------------------

! Gauss integration point 7 and weight factor
subroutine Gauss7_Point_xy(i, r, s, w)
    integer, intent(in) :: i
    double precision, intent(out) :: r, s, w

    select case (i)
        case (1); r = 0.1012865073235d0; s = 0.1012865073235d0; w = 0.1259391805448d0
        case (2); r = 0.7974269853531d0; s = 0.1012865073235d0; w = 0.1259391805448d0
        case (3); r = 0.1012865073235d0; s = 0.7974269853531d0; w = 0.1259391805448d0
        case (4); r = 0.4701420641051d0; s = 0.0597158717898d0; w = 0.1323941527885d0
        case (5); r = 0.4701420641051d0; s = 0.4701420641051d0; w = 0.1323941527885d0
        case (6); r = 0.0597158717898d0; s = 0.4701420641051d0; w = 0.1323941527885d0
        case (7); r = 0.3333333333333d0; s = 0.3333333333333d0; w = 0.2250000000000d0
    end select

    w = 0.5d0 * w
end subroutine Gauss7_Point_xy

! ---------------------------------------------------------------------------------------

! Director vector
subroutine Director(Vn, V1, V2)
    double precision, intent(in)  :: Vn(3)  ! Director vector
    double precision, intent(out) :: V1(3)  ! V1
    double precision, intent(out) :: V2(3)  ! V2

    double precision e2(3) / 0.0d0, 1.0d0, 0.0d0 /

    ! V1, V2 from Vn
    if(Equal_Vector(Vn, e2)) then

        ! For the special case Vn = e2
        V1(1) = 0.0d0; V1(2) = 0.0d0; V1(3) = 1.0d0
    else

        ! Other cases
        call Cross_Product(e2(:), Vn(:), V1(:))
        call Normalize(V1(:))
    end if

    call Cross_Product(Vn(:), V1(:), V2(:))
end subroutine Director

! ---------------------------------------------------------------------------------------

! cross product
subroutine Cross_Product(x, y, z)
    double precision, intent(in)  :: x(3), y(3)
    double precision, intent(out) :: z(3)
    
    if(x(1)==y(1) .and. x(2)==y(2) .and. x(3)==y(3)) then

        z(1) = 0.0d0; z(2) = 0.0d0; z(3) = 1.0d0
    else

        z(1) = x(2)*y(3) - x(3)*y(2)
        z(2) = x(3)*y(1) - x(1)*y(3)
        z(3) = x(1)*y(2) - x(2)*y(1)
    end if
end subroutine Cross_Product

! ---------------------------------------------------------------------------------------

! Normalization of vector
subroutine Normalize(v)
    double precision, intent(inout) :: v(3)
    double precision :: length

    length = sqrt(v(1)**2.0d0 + v(2)**2.0d0 + v(3)**2.0d0)
    v(:) = v(:) / length
end subroutine Normalize

! ---------------------------------------------------------------------------------------

! Size vector
function Size_Vector(vec) result(size)
    double precision, intent(in) :: vec(3)

    double precision :: size

    size = sqrt(vec(1)**2.0d0 + vec(2)**2.0d0 + vec(3)**2.0d0)
end function Size_Vector

! ---------------------------------------------------------------------------------------

! Determinant
function Determinant_33(mat) result(deter)
    double precision, intent(in)  :: mat(3,3)

    double precision :: deter

    deter &
        = mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
        - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
        + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
end function Determinant_33

! ---------------------------------------------------------------------------------------

! Function to find the determinant of a square matrix
double precision function Det(matrix)
    double precision :: matrix(:,:)

    double precision :: m, temp
    integer :: i, j, k, l, n
    logical :: detexists = .true.

    n = ubound(matrix, 1)
    select case (n)
    case (2)
        det = matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
    case (3)
        det = matrix(1,1)*matrix(2,2)*matrix(3,3)+ &
              matrix(1,2)*matrix(2,3)*matrix(3,1)+ &
              matrix(1,3)*matrix(2,1)*matrix(3,2)- &
              matrix(1,1)*matrix(2,3)*matrix(3,2)- &
              matrix(1,2)*matrix(2,1)*matrix(3,3)- &
              matrix(1,3)*matrix(2,2)*matrix(3,1)
    case (4:)
        l = 1

        ! Convert to upper triangular form
        do k = 1, n-1
            if (matrix(k,k) == 0) then
                detexists = .false.
                do i = k+1, n
                    if (matrix(i,k) /= 0) then
                        do j = 1, n
                            temp = matrix(i,j)
                            matrix(i,j)= matrix(k,j)
                            matrix(k,j) = temp
                        end do
                        detexists = .true.
                        l=-l
                        exit
                    endif
                end do
                if (detexists .eqv. .false.) then
                    det = 0
                    return
                end if
            endif
            do j = k+1, n
                m = matrix(j,k)/matrix(k,k)
                do i = k+1, n
                    matrix(j,i) = matrix(j,i) - m*matrix(k,i)
                end do
            end do
        end do

        ! Calculate determinant by finding product of diagonal elements
        det = l
        do i = 1, n
            det = det * matrix(i,i)
        end do
    end select
end function det

! ---------------------------------------------------------------------------------------

! The determinant of a double precision square matrix mat by Gauss method with full pivoting
function Determinant(mat, eps) result(det)
    double precision, intent(in) :: mat(:,:)
    double precision, intent(in) :: eps

    double precision, allocatable :: diag(:,:)
    integer, allocatable :: kp(:)
    integer, allocatable :: lp(:)
    double precision :: det
    integer :: i, count, dim, flag

    dim = ubound(mat, 1)

    ! Allocate local matrix diag and vectors kp, lp
    allocate(diag(dim, dim))
    allocate(Kp(dim))
    allocate(Lp(dim))

    ! Triangularization subroutine
    flag = Triangularization(mat, eps, diag, kp, lp)

    if(flag == 0) then
        ! If the matrix singular, det = 0
        det = 0.0d0
    else
        ! If the matrix regular
        det = 1.0d0
        do i = 1, dim
            det = det * diag(i,i)
        end do

        count = 0
        do i = 1, dim - 1
            if(lp(i) /= i) count = count + 1
            if(kp(i) /= i) count = count + 1
        end do

        ! If count is odd
        if(mod(count, 2) /= 0) det = -det
    end if

    ! Dellocate memory
    deallocate(diag)
    deallocate(Kp)
    deallocate(Lp)
end function Determinant

! ---------------------------------------------------------------------------------------

! Inverse matrix
subroutine Inverse_22(a, ainv, ok_flag)
    double precision, dimension(2,2), intent(in)  :: a
    double precision, dimension(2,2), intent(out) :: ainv
    logical, intent(out) :: ok_flag

    double precision, parameter :: eps = 1.0d-10
    double precision, dimension(2,2) :: cofactor
    double precision :: det

    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    if (abs(det) <= eps) then
        ainv = 0.0d0
        ok_flag = .false.
        return
    end if

    cofactor(1,1) = +a(2,2)
    cofactor(1,2) = -a(2,1)
    cofactor(2,1) = -a(1,2)
    cofactor(2,2) = +a(1,1)

    ainv = transpose(cofactor) / det

    ok_flag = .true.
end subroutine Inverse_22

! ---------------------------------------------------------------------------------------

! Inverse matrix
subroutine Inverse_33(mat, inv_mat)
    double precision, intent(in)  :: mat(3,3)
    double precision, intent(out) :: inv_mat(3,3)

    double precision :: det
    
    inv_mat(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
    inv_mat(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
    inv_mat(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
    inv_mat(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
    inv_mat(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
    inv_mat(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
    inv_mat(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
    inv_mat(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
    inv_mat(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

    det = mat(1,1) * (mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) &
        - mat(1,2) * (mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) &
        + mat(1,3) * (mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))

    inv_mat(:,:) = inv_mat(:,:) / det
end subroutine Inverse_33

! ---------------------------------------------------------------------------------------

! Inverse matrix based on Doolittle LU factorization for Ax=b
subroutine Inverse_LU(mat)
    double precision, intent(inout) :: mat(:,:)

    double precision, allocatable :: L(:,:), U(:,:), a(:,:), c(:,:)
    double precision, allocatable :: b(:), d(:), x(:)
    double precision :: coeff
    integer :: i, j, k, n

    n = ubound(mat, 1)

    allocate(L(n,n))
    allocate(U(n,n))
    allocate(a(n,n))
    allocate(c(n,n))
    allocate(b(n))
    allocate(d(n))
    allocate(x(n))

    a = mat

    ! Step 0: Initialization for matrices L and U and b
    L = 0.0d0
    U = 0.0d0
    b = 0.0d0

    ! Step 1: Forward elimination
    do k = 1, n - 1
        do i = k + 1, n
            coeff  = a(i,k) / a(k,k)
            L(i,k) = coeff
            do j = k + 1, n
                a(i,j) = a(i,j) - coeff * a(k,j)
            end do
        end do
    end do

    ! Step 2: Prepare L and U matrices

    ! L matrix is a matrix of the elimination coefficient
    do i = 1, n
        L(i,i) = 1.0d0
    end do

    ! U matrix is the upper triangular part of A
    do j = 1, n
        do i = 1, j
            U(i,j) = a(i,j)
        end do
    end do

    ! Step 3: Compute columns of the inverse matrix C
    do k = 1, n
        b(k) = 1.0d0
        d(1) = b(1)

        ! Step 3a: Solve Ld=b using the forward substitution
        do i = 2, n
            d(i) = b(i)
            do j = 1, i - 1
                d(i) = d(i) - L(i,j) * d(j)
            end do
        end do

        ! Step 3b: Solve Ux=d using the back substitution
        x(n) = d(n) / U(n,n)
        do i = n - 1, 1, -1
            x(i) = d(i)
            do j = n, i + 1, -1
                x(i) = x(i) - U(i,j) * x(j)
            end do
            x(i) = x(i) / u(i,i)
        end do

        ! Step 3c: Fill the solutions x(n) into column k of C
        do i = 1, n
            c(i,k) = x(i)
        end do
        b(k) = 0.0d0
    end do

    mat = c
end subroutine Inverse_LU

! ---------------------------------------------------------------------------------------

! Inverse matrix based on Doolittle LU factorization for Ax=b
subroutine InvMat(MatIn,MatOut,n)
    integer :: n
    double precision :: MatIn(n,n), MatOut(n,n)
    double precision :: L(n,n), U(n,n), b(n), d(n), x(n)
    double precision :: coeff
    integer :: i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L = 0.0d0
    U = 0.0d0
    b = 0.0d0

    ! step 1: forward elimination
    do k = 1, n - 1
        do i = k + 1, n
            coeff  = MatIn(i,k) / MatIn(k,k)
            L(i,k) = coeff
            do j = k + 1, n
                MatIn(i,j) = MatIn(i,j) - coeff * MatIn(k,j)
            end do
        end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i = 1, n
        L(i,i) = 1.0d0
    end do

    ! U matrix is the upper triangular part of A
    do j = 1, n
        do i = 1, j
            U(i,j) = MatIn(i,j)
        end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k = 1, n
        b(k) = 1.0d0
        d(1) = b(1)

        ! Step 3a: Solve Ld=b using the forward substitution
        do i = 2, n
            d(i) = b(i)
            do j = 1, i - 1
                d(i) = d(i) - L(i,j)*d(j)
            end do
        end do

        ! Step 3b: Solve Ux=d using the back substitution
        x(n) = d(n) / U(n,n)
        do i = n - 1, 1, -1
            x(i) = d(i)
            do j = n,i + 1, -1
                x(i) = x(i) - U(i,j) * x(j)
            end do
            x(i) = x(i) / u(i,i)
        end do

        ! Step 3c: fill the solutions x(n) into column k of C
        do i = 1, n
            MatOut(i,k) = x(i)
        end do
        b(k) = 0.0d0
    end do
end subroutine InvMat

! ---------------------------------------------------------------------------------------

! Find difference between two vectors
function Find_Diff(vec_a, vec_b)
    integer, intent(in) :: vec_a(:), vec_b(:)
    integer, dimension(size(vec_a) - size(vec_b)) :: Find_Diff

    integer :: i, j, count
    logical :: exist

    count = 0
    do i = 1, size(vec_a)
        exist = .false.
        do j = 1, size(vec_b)
            if(vec_a(i) == vec_b(j)) then
                exist = .true.
                exit
            end if
        end do

        if(exist == .false.) then
            count = count + 1
            Find_Diff(count) = vec_a(i)
        end if
    end do
end function Find_Diff

! ---------------------------------------------------------------------------------------

! Is two vector identical
function Equal_Vector(Va, Vb) result(equal)
    double precision, intent(in) :: Va(3)
    double precision, intent(in) :: Vb(3)

    logical :: equal

    if( abs(Va(1) - Vb(1)) < eps .and. &
        abs(Va(2) - Vb(2)) < eps .and. &
        abs(Va(3) - Vb(3)) < eps ) then
        equal = .true.
    else
        equal = .false.
    end if
end function Equal_Vector

! ---------------------------------------------------------------------------------------

! Return the location of the minimum in the section between start and end
function Find_Minimum(array, start, end) result(value)
    integer, dimension(1:), intent(in) :: array
    integer, intent(in) :: start
    integer, intent(in) :: end

    integer :: i, min, loc, value

    min = array(start)  ! Assume the first is the min
    loc = start         ! Record its position

    ! Start with next elements
    do i = start+1, end

        ! If array(i) less than the min?
        if(array(i) < min) then

            ! Record its position
            min = array(i)
            loc = i
        end if
    end do

    ! Return the position
    value = loc
end function Find_Minimum

! ---------------------------------------------------------------------------------------

! Swap the values of its two formal arguments
subroutine Swap(entity_a, entity_b)
    integer, intent(inout) :: entity_a
    integer, intent(inout) :: entity_b

    integer :: entity_t

    entity_t = entity_a
    entity_a = entity_b
    entity_b = entity_t
end subroutine Swap

! ---------------------------------------------------------------------------------------

! Receives an array() and sorts it into ascending order
subroutine Sort(array, size)
    integer, dimension(1:), intent(inout) :: array
    integer, intent(in) :: size

    integer :: i, loc

    ! Except for the last
    do i = 1, size - 1

        ! Find min from this to last
        loc = Find_Minimum(array, i, size)

        ! Swap this and the minimum
        call  Swap(array(i), array(loc))
    end do
end subroutine Sort

! ---------------------------------------------------------------------------------------

! Receives an array() and sorts it into ascending order
subroutine Sort2(array1, array2, size)
    integer, dimension(1:), intent(inout) :: array1
    integer, dimension(1:), intent(inout) :: array2
    integer, intent(in) :: size

    integer :: i, loc

    ! Except for the last
    do i = 1, size - 1

        ! Find min from this to last
        loc = Find_Minimum(array1, i, size)

        ! Swap this and the minimum
        call  Swap(array1(i), array1(loc))
        call  Swap(array2(i), array2(loc))
    end do
end subroutine Sort2

! ---------------------------------------------------------------------------------------

! The upper triangularization algorithm of Gauss method with full pivoting
function Triangularization(mat, eps, diag, kp, lp) result(flag)
    double precision, intent(in)    :: mat(:,:)
    double precision, intent(in)    :: eps
    double precision, intent(inout) :: diag(:,:)
    integer, intent(inout) :: kp(:)
    integer, intent(inout) :: lp(:)

    double precision :: po, t0
    integer :: i, j, k, lo, ko, flag, n

    n = ubound(mat, 1)

    diag = mat
    flag = 1
    k    = 1

    do while(flag == 1 .and. k < n)
        po = diag(k,k)
        lo = k
        ko = k

        do i = k, n
            do j = k, n
                if(abs(diag(i,j)) > abs(po)) then
                    po = diag(i,j)
                    lo = i
                    ko = j
                end if
            end do
        end do
        Lp(k) = lo
        Kp(k) = ko

        if(abs(po) < eps) then
            flag = 0
        else
            if(lo /= k) then
                do j = k, n
                    t0      = diag(k,j)
                    diag(k,j)  = diag(lo,j)
                    diag(lo,j) = t0
                end do
            end if
            if(ko /= k) then
                do i = 1, n
                    t0      = diag(i,k)
                    diag(i,k)  = diag(i,ko)
                    diag(i,ko) = t0
                end do
            end if
            do i = k + 1, n
                diag(i,k) = diag(i,k) / po
                do j = k + 1, n
                    diag(i,j) = diag(i,j) - diag(i,k) * diag(k,j)
                end do
            end do
            k = k + 1
        end if
    end do

    if(flag == 1 .and. abs(diag(n,n)) < eps) flag = 0
end function Triangularization

! ---------------------------------------------------------------------------------------

! Function calculating cross product of 3D vectors
function Cross(a, b)
    double precision, intent(in) :: a(3), b(3)
    double precision :: cross(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
end function Cross

! ---------------------------------------------------------------------------------------

! Subroutine Judging if four points are in the same place
subroutine IsPlanar(judge, nor_vec, n_point, position)
    double precision :: position(n_point*3), nor_vec(3), V(3,3), crit
    integer :: n_point
    logical :: judge

    judge = .false.
    if(n_point == 4) then
        V = reshape((/position(4:11:3) - position(1)*(/1,1,1/), &
                      position(5:12:3) - position(2)*(/1,1,1/), &
                      position(6:13:3) - position(3)*(/1,1,1/)/), shape(V))
        nor_vec = cross(V(2,:), V(3,:)-V(1,:))
        crit    = dot_product(nor_vec,V(1,:))

        if(abs(crit) < 1.0d-10) then
            judge = .True.
            return
        end if
    end if
end subroutine IsPlanar

! ---------------------------------------------------------------------------------------

! Function Returning Weight Masks for Gaussian Quadrature
subroutine GaussianMask(gp, weight, order)
    integer :: order
    double precision :: gp(order), weight(order)

    select case (order)
    case(1)
        gp     = 0.0d0
        weight = 2.0d0
    case(2)
        gp     = (/1.0d0 / sqrt(3.0d0), -1.0d0/sqrt(3.0d0)/)
        weight = (/1.0d0, 1.0d0/)
    case(3)
        gp     = (/sqrt(3.0d0/5.0d0), 0.0d0, -sqrt(3.0d0/5.0d0)/)
        weight = (/5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0/)
    end select
end subroutine GaussianMask

! ---------------------------------------------------------------------------------------

! Subroutine performing spacial rotations along a given axis
subroutine AxialRotate(PtIn, PtOut, axis, theta)
    double precision :: PtIn(3), PtOut(3), axis(3), theta
    double precision :: C, S, R(3,3)

    C    = cos(theta)
    S    = sin(theta)
    axis = axis/sqrt(dot_product(axis, axis))

    R = reshape((/axis(1)**2+(1.0d0-axis(1)**2)*C,     axis(1)*axis(2)*(1.0d0-C)-axis(3)*S, axis(1)*axis(3)*(1.0d0-C)+axis(2)*S, &
                  axis(1)*axis(2)*(1.0d0-C)+axis(3)*S, axis(2)**2+(1.0d0-axis(2)**2)*C,     axis(2)*axis(3)*(1.0d0-C)-axis(1)*S, &
                  axis(1)*axis(3)*(1.0d0-C)-axis(2)*S, axis(2)*axis(3)*(1.0d0-C)+axis(1)*S, axis(3)**2+(1.0d0-axis(3)**2)*C/),(/3,3/))
    PtOut = matmul(R, PtIn)
end subroutine AxialRotate

! ---------------------------------------------------------------------------------------

function Rotate_Vector(Vn, theta_a, theta_b) result(Q)
    double precision, intent(in) :: Vn(3), theta_a, theta_b
    
    double precision :: theta_x, theta_y, theta_z, gamma
    double precision :: Q(3,3), S(3,3), V1(3), V2(3)

    Q(:,:) = 0.0d0
    S(:,:) = 0.0d0

    call Director(Vn, V1, V2)

    theta_x = theta_a * V1(1) + theta_b * V2(1)
    theta_y = theta_a * V1(2) + theta_b * V2(2)
    theta_z = theta_a * V1(3) + theta_b * V2(3)

    gamma = sqrt(theta_x**2.0d0 + theta_y**2.0d0 + theta_z**2.0d0)

    Q(1,1) = 1.0d0
    Q(2,2) = 1.0d0
    Q(3,3) = 1.0d0

    if(abs(gamma) /= 0.0d0) then
        S(1,2) = -theta_z
        S(1,3) =  theta_y
        S(2,1) =  theta_z
        S(2,3) = -theta_x
        S(3,1) = -theta_y
        S(3,2) =  theta_x
        Q      =  Q + sin(gamma)/gamma*S + 0.5d0 * (sin(gamma/2.0d0)/(gamma/2.0d0))**2.0d0 * matmul(S,S)
    end if
end function Rotate_Vector

! ---------------------------------------------------------------------------------------

! Perform least-square approximation
subroutine LeastSquare(coeff, value, Ncoeff, Nval, sets)
    integer :: Ncoeff, Nval, sets
    double precision :: coeff(Ncoeff, sets), value(Nval, Ncoeff + sets)

    double precision :: ATA(Ncoeff, Ncoeff), ATY(Ncoeff)
    integer :: i

    ATA = matmul(transpose(value(:, 1:Ncoeff)), value(:, 1:Ncoeff))
    call InvMat(ATA, ATA, Ncoeff)

    do i = 1, sets
        ATY = matmul(transpose(value(:, 1:Ncoeff)), value(:, Ncoeff + i))
        coeff(:,i) = matmul(ATA, ATY)
    end do
end subroutine LeastSquare

! ---------------------------------------------------------------------------------------

end module Math_Kernel