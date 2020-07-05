!
! =============================================================================
!
! SATURN Educational Version
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
program SATURN_EDU

    use Math_Kernel
    use Data_Struct
    use DISP4
    use Solver

    implicit none

    ! define internal variables
    integer :: node_n, element_n            ! # of nodes, # of elements
    type(NodeType), allocatable :: node(:)  ! Node
    type(ElemType), allocatable :: elem(:)  ! Element
    type(PropType) :: prop                  ! Property
    type(ProbType) :: prob                  ! Problem

    integer :: dof(3)                       ! # of DOFs (total, free, fixed)
    integer, allocatable :: maxa(:)         ! Addresses of diagonal elements
    double precision, allocatable :: Kt(:)  ! Stiffness vector
    double precision, allocatable :: U(:)   ! Diplacement vector
    double precision, allocatable :: R(:)   ! Load vector
    integer :: n_eq                         ! Number of equations
    integer :: n_wk                         ! Number of elements below skyline of matrix
    integer :: m_hbw                        ! Maximum half bandwidth

    call Main

contains

! --------------------------------------------------------------------------------

! Main subroutine
subroutine Main()

    ! Define variables
    double precision :: ts_main, te_main, ts_assem, te_assem, ts_solve, te_solve

    call cpu_time(ts_main)

    ! Open files
    open(unit=1, file="SATURN.txt",     form="formatted")
    open(unit=2, file="SATURN_out.txt", form="formatted")
    open(unit=3, file="SATURN_res.txt", form="formatted")
    open(unit=4, file="SATURN_pos.txt", form="formatted")

    ! Read input file
    !call Read_Input_File
    call Set2DQuadEleRectDomain

    ! # of DOFs and assign equation numbers
    n_eq = Set_DOF_Number(node, dof)

    ! Column index
    call SKY_Index(node, elem, maxa, n_eq, n_wk, m_hbw)

    ! Allocate memory
    allocate(Kt(n_wk), U(dof(1)), R(dof(1)))

    ! Print information
    call Print_Information(prob, node, elem, prop, dof, n_eq, n_wk, m_hbw)

    ! Assemble total stiffness matrix
    write(0, "(a)"), " 1 - Assembling Stiffness and Load"
    call cpu_time(ts_assem)
    call Assemble_Kt(prob, node, elem, prop, Kt, maxa, n_eq)

    ! Assemble load vector
    call Assemble_Load(node, elem, prop, R)
    call cpu_time(te_assem)
    call cpu_time(ts_solve)
    
    ! slove linear system (find U in K*U=R)
    print *, "2/ Solving Linear System"
    call Solve_Equations
    call cpu_time(te_solve)

    ! calculate stress and print solutions
    print *, "3/ Printing Output Files"
    call Displacement_Stress

    print *, "4/ Completed !"; print * 
    print *, "Strain energy = ", 0.5d0 * dot_product(R(1:dof(2)), U(1:dof(2)))
    PRINT *, "Strain energy =   4.447735498387236E-008"
    ! deallocate memory
    deallocate(node, elem, Kt, R, U, maxa)
    call cpu_time(te_main)

    ! close files
    close(unit=1); close(unit=2); close(unit=3); close(unit=4)
    
    ! Compute time-consuming
    call Compute_Cost(te_main, ts_main, te_assem, ts_assem, te_solve, ts_solve)
end subroutine Main

! --------------------------------------------------------------------------------

! read input file
subroutine Read_Input_File()
    double precision :: thickness
    character :: bufs
    integer :: bufi, i

    ! read nodal information
    read(1,*) bufs; read(1,*) node_n; read(1,*) bufs
    allocate(Node(node_n)) ! allocate a node array
      
    do i=1, node_n
        read(1,*) bufi, Node(i).x(1:2), Node(i).bc(1:2), Node(i).pm(1:2)
    end do

    ! read elemental information
    read(1,*) bufs; read(1,*) element_n; read(1,*) bufs
    allocate(elem(element_n)) ! allocate a element array
    
    do i=1, element_n
        read(1,*) bufi, elem(i).cn(:), elem(i).q(:)
    end do

    ! read properties
    read(1,*) bufs; read(1,*) thickness
    prop.thick = thickness
    read(1,*) bufs; read(1,*) prop.young
    read(1,*) bufs; read(1,*) prop.poisson
end subroutine Read_Input_File

! --------------------------------------------------------------------------------

! # of DOFs and assign equation numbers to DOFs
! dof(1): # of total DOFs, dof(2): # of free DOFs, dof(3): # of fixed DOFs
function Set_DOF_Number(node, dof) result(n_eq)
    type(NodeType), intent(inout) :: node(:)
    integer, intent(out) :: dof(3)

    integer :: i, j, n_eq

    write(3, *) "EQUATION NUMBER"
    write(3, *) "---------------------"
    write(3, *) "    node   dof    eqn"

    dof(1) = size(node) * nDPN
    dof(2) = 0
    dof(3) = 0

    do i = 1, size(node)
        do j = 1, nDPN
            if (node(i)%bc(j) == 0) then
                dof(2) = dof(2) + 1
                node(i)%eq_n(j) = dof(2)
                write(3, "(3i7)") i, j, dof(2)
            else
                dof(3) = dof(3) + 1
                node(i)%eq_n(j) = dof(1) - dof(3) + 1
            end if
        end do
    end do
    write(3, *)

    n_eq = dof(2)
end function Set_DOF_Number

! --------------------------------------------------------------------------------

! Calculate column index
subroutine SKY_Index(node, elem, maxa, n_eq, n_wk, m_hbw)
    integer, allocatable, intent(out) :: maxa(:)
    type(NodeType), intent(in) :: node(:)
    type(ElemType), intent(in) :: elem(:)
    integer, intent(in)  :: n_eq
    integer, intent(out) :: n_wk
    integer, intent(out) :: m_hbw

    integer, allocatable :: column_h(:)
    integer, allocatable :: a_index(:)
    integer :: i, j, k, n_elem
 n_elem = size(elem)

    ! Allocate maxa array
    allocate(column_h(n_eq))
    allocate(a_index(nDPN*nNPE))
    allocate(maxa(n_eq+1))

    column_h(:) = 0

    do i = 1, n_elem

        ! Assemblage index
        do j = 1, nDPN
            do k = 1, nNPE
                a_index(nNPE*j+k-nNPE) = node(elem(i)%cn(k))%eq_n(j)
            end do
        end do

        ! Column height
        do k = 1, nDPN * nNPE
            do j = 1, nDPN * nNPE
                if(a_index(j) <= n_eq .and. a_index(k) <= n_eq) then
                    if(a_index(j) < a_index(k)) then
                        if(a_index(k)-a_index(j) > column_h(a_index(k))) then
                            column_h(a_index(k)) = a_index(k) - a_index(j)
                        end if
                    end if
                end if
            end do
        end do
    end do

    ! maxa array
    do i = 1, n_eq + 1
        maxa(i) = 0
    end do

    maxa(1) = 1
    maxa(2) = 2
    m_hbw   = 0
    if(n_eq > 1) then
        do i = 2, n_eq
            if(column_h(i) > m_hbw) m_hbw = column_h(i)
            maxa(i + 1) = maxa(i) + column_h(i) + 1
        end do
    end if

    m_hbw = m_hbw + 1
    n_wk  = maxa(n_eq + 1) - maxa(1)

    ! Deallocate
    deallocate(column_h)
    deallocate(a_index)
end subroutine SKY_Index

! --------------------------------------------------------------------------------

! Assemble total stiffness matrix
subroutine Assemble_Kt(prob, node, elem, prop, Kt, maxa, n_eq)
    type(NodeType), intent(in) :: node(:)
    type(ElemType), intent(in) :: elem(:)
    type(PropType), intent(in) :: prop
    type(ProbType), intent(in) :: prob
    double precision, intent(out) :: Kt(:)
    integer, intent(in) :: maxa(:)
    integer, intent(in) :: n_eq

    double precision, allocatable :: Ke(:,:)
    integer, allocatable :: a_index(:)
    double precision :: eNode(nNPE, nDIM)
    integer :: i, j, k, address

    allocate(Ke(nDPN*nNPE, nDPN*nNPE))
    allocate(a_index(nDPN*nNPE))

    Kt(:) = 0.0d0

    ! Assemble stiffness matrix
    do i = 1, size(elem)

        ! Nodal position of elements
        do j = 1, nNPE
            do k = 1, nDIM
                eNode(j, k) = node(elem(i)%cn(j))%x(k)
            end do
        end do

        ! Planestress element stiffness matrix
        call Plane_Stiffness(prop.young, prop.poisson, prop.thick, eNode, Ke)

        ! Print all element stiffness
        !call Print_Matrix(i, Ke)

        ! Assemblage index
        do j = 1, nDPN
            do k = 1, nNPE
                a_index(nNPE*j+k-nNPE) = node(elem(i)%cn(k))%eq_n(j)
            end do
        end do

        ! Assemble stiffness matrix
        do j = 1, nDPN * nNPE
            do k = 1, nDPN * nNPE
                if(a_index(j) <= n_eq .and. a_index(k) <= n_eq) then
                    if(a_index(j) <= a_index(k)) then
                        address = maxa(a_index(k)) + a_index(k) - a_index(j)
                        Kt(address) = Kt(address) + Ke(j, k)
                    end if
                end if
            end do
        end do
    end do

    ! Deallocate memory
    deallocate(Ke, a_index)
end subroutine Assemble_Kt

! --------------------------------------------------------------------------------

! assemble load vector
subroutine Assemble_Load(node, elem, prop, R)
    type(NodeType),   intent(in)  :: node(:)
    type(ElemType),   intent(in)  :: elem(:)
    type(PropType),   intent(in)  :: prop
    double precision, intent(out) :: R(:)

    double precision :: eNode(nNPE, nDIM)
    double precision :: NodalLoad(nDPN*nNPE)    ! Equivalent nodal load
    integer :: i, j, k

    R(:) = 0

    ! Assemble load vector for nodal load
    do i = 1, size(node)
        do j = 1, nDPN
            R(node(i)%eq_n(j)) = node(i)%pm(j)
        end do
    end do

    ! Assemble load vector for body force
    do i = 1, size(elem)

        ! Nodal position of elements
        do j = 1, nNPE
            do k = 1, nDIM
                eNode(j,k) = node(elem(i)%cn(j))%x(k)
            end do
        end do

        ! Calculate equivalent nodal load from body force
        call Plane_Load(eNode, elem(i).q, NodalLoad)

        ! Assemble load vector
        do j = 1, nDPN
            do k = 1, nNPE
                R(Node(elem(i).cn(k)).eq_n(j)) &
                    = R(Node(elem(i).cn(k)).eq_n(j)) + NodalLoad(nNPE*j+k-nNPE)
            end do
        end do
    end do
end subroutine Assemble_Load

! --------------------------------------------------------------------------------

! love linear equations
subroutine Solve_Equations()
    integer :: i

    U(:) = 0.0d0; i=3
    U(1:dof(2)) = R(1:dof(2))
    call Skyline_COLSOL(Kt(1:maxa(dof(2)+1)-1), U(1:dof(2)), maxa(1:dof(2)+1), &
        dof(2),  maxa(dof(2)+1)-1, dof(2)+1, 1, i)
    call Skyline_COLSOL(Kt(1:maxa(dof(2)+1)-1), U(1:dof(2)), maxa(1:dof(2)+1), &
        dof(2),  maxa(dof(2)+1)-1, dof(2)+1, 2, i)
end subroutine Solve_Equations

! --------------------------------------------------------------------------------

! calculate stress and print solutions
subroutine Displacement_Stress()
    double precision :: eNode(nNPE,nDIM)                ! nodal position of element
    double precision :: displace(nNPE*nDPN)             ! nodal displacement vector of element
    double precision :: Stress(nNPE,3)                  ! Sxx, Syy, Sxy in Gauss points or nodal positions (4)
    double precision :: scale_factor, max_pos, max_disp ! scaling factor
    integer :: i, j, k
   
    ! print strain energy
    write(3,"(a17, E14.6)") "STRAIN ENERGY = ", 0.5d0 * dot_product(R(1:dof(2)), U(1:dof(2)))
    write(4,*) element_n

    ! set scaling factor for the plotting
    max_pos  = 0.0d0
    max_disp = 0.0d0
    do i=1, node_n
        if( max_disp < dabs(U(Node(i).eq_n(1))) ) then
            max_disp = dabs(U(Node(i).eq_n(1)))
            max_pos  = sqrt(Node(i).x(1)*Node(i).x(1)+Node(i).x(2)*Node(i).x(2))
        end if

        if( max_disp < dabs(U(Node(i).eq_n(2))) ) then
            max_disp = dabs(U(Node(i).eq_n(2)))
            max_pos  = sqrt(Node(i).x(1)*Node(i).x(1)+Node(i).x(2)*Node(i).x(2))
        end if
    end do
      
    scale_factor = (1.2d0 * max_pos - max_pos) / max_disp        ! 1.2 * max_pos = (scale_factor * max_disp + max_pos)

    write(4,"(E14.6)"), scale_factor

    ! print nodal displacement
    write(3,*)
    write(3,*) "DISPLACEMENT "
    write(3,*) "------------------------------"
    write(3,*) "  Node      Dx         Dy     "
   
    do i=1, node_n
        write(3,"(1x,i4,2x,2(1P,E11.3))") i, U(Node(i).eq_n(1)), U(Node(i).eq_n(2))
    end do
    write(3,*)

    do i=1, element_n
        ! nodal position of element
        do j=1, nNPE
            do k=1, nDIM
                eNode(j, k) = Node( elem(i).cn(j) ).x(k)
            end do
        end do
      
        ! displacement vector of element 
        do j=1, nDPN
            do k=1, nNPE
                displace(nNPE * j + k - nNPE) = U( Node( elem(i).cn(k) ).eq_n(j) )
            end do
        end do

        ! calculate stress of element
        call Plane_Stress(prop.young, prop.poisson, eNode, displace, Stress)
      
        ! print stress
        write(3,"(a21,i4)") " STRESS of ELEMENT : ", i
        write(3,*) "----------------------------------------------" 
        write(3,*) " Position       Sxx        Syy        Sxy     "
      
        do j=1, nNPE
            write(3,"(1x,i4,5x, 3x,3(1P,E11.3))") j, Stress(j,:)
        end do
        write(3,*)
        
        ! print deformed shape and stress for post processing by using MATLAB 
        write(4,"(1x,28(1P,E13.5))") eNode(1,:), displace(1), displace(5), Stress(1,:),&
                                      eNode(2,:), displace(2), displace(6), Stress(2,:),&
                                      eNode(3,:), displace(3), displace(7), Stress(3,:),&
                                      eNode(4,:), displace(4), displace(8), Stress(4,:)
    end do
end subroutine Displacement_Stress

! --------------------------------------------------------------------------------

! print matrix
subroutine Print_Matrix(M)
    double precision, intent(in) :: M(:,:)
    integer :: i

    write(2,*) "---------------------------"
    do i=1, nNPE * nDPN
        write(2,"(8E12.4)" ) M(i,:)
    end do
    write(2,*)
end subroutine Print_Matrix

! --------------------------------------------------------------------------------

! Print information
subroutine Print_Information(prob, node, elem, prop, dof, n_eq, n_wk, m_hbw)
    type(NodeType), intent(in) :: node(:)
    type(ElemType), intent(in) :: elem(:)
    type(PropType), intent(in) :: prop
    type(ProbType), intent(in) :: prob
    integer, intent(in) :: dof(3)
    integer, intent(in) :: n_eq
    integer, intent(in) :: n_wk
    integer, intent(in) :: m_hbw

    write(0, "(a)"), " [SATURN Edu. Ver.]"
    write(0, "(a)")
    write(0, "(a)"), " ===================================================================="
    write(0, "(a)")
    write(0, "(a$)"), " Problem name: "//trim(prob%name)
    write(0, "(a$)"), ", Domain size: "//trim(adjustl(Int2Str(prob%n_domain)))
    write(0, "(a$)"), " Problem : "//trim(adjustl(prob%name))
    write(0, "(a)"), " Young's modulus: "//trim(adjustl(ES2Str(prop%Young, 9,3)))//&
                     ", Poisson ratio: "//trim(adjustl(Dble2Str(prop%Poisson, 8,4)))//&
                     ", thickness: "//trim(adjustl(Dble2Str(prop%thick, 8,4)))
    write(0, "(a)"), " # total elements: "//trim(adjustl(Int2Str(size(elem))))//&
                     ", # total nodes: "//trim(adjustl(Int2Str(size(node))))
    write(0, "(a)"), " # total DOFs: "//trim(adjustl(Int2Str(dof(1))))//&
                     ", # free DOFs: "//trim(adjustl(Int2Str(dof(2))))//&
                     ", # fixed DOFs: "//trim(adjustl(Int2Str(dof(3))))
    write(0, "(a$)")," # of equations = "//trim(adjustl(Int2Str(n_eq)))//","
    write(0, "(a )")," # of matrix elements = "//trim(adjustl(Int2Str(n_wk)))
    write(0, "(a$)")," Max. half bandwidth = "//trim(adjustl(Int2Str(m_hbw)))//","
    write(0, "(a)"), " Required memory = "//trim(adjustl(Dble2Str(real(n_wk)*8.0d0/1024.0d0/1024.0d0, 10, 3)))//" MB"
    write(0, "(a)")
    write(0, "(a)"), " ===================================================================="
    write(0, "(a)")
end subroutine Print_Information

! --------------------------------------------------------------------------------

! set rectangular domain with triangular element
subroutine Set2DQuadEleRectDomain()
    double precision :: young                       ! Young's modulus
    double precision :: possion                     ! Possion ratio
    double precision :: bxforce,    byforce         ! body force
    double precision :: pxforce,    pyforce         ! point load
    double precision :: dxforce,    dyforce         ! distributed load
    double precision :: thick                       ! thinkness
    double precision :: x_width,    y_width         ! length
    double precision :: n_i_node,   n_j_node        ! the number of nodes in the x,y-direction
    double precision :: n_i_element, n_j_element
    double precision :: del_x,      del_y
    integer :: i, j, k, domainsize, numbering, x_fix_surf

    domainsize  = 10

    young       =  1.7472e7             ! Young's modulus
    possion     =  0.3                  ! Possion ratio
    bxforce     =  0.0d0                ! x-direction body force
    byforce     = -1.0d0                ! y-direction body force
    pxforce     =  0.0d0                ! x-direction point load
    pyforce     =  0.0d0                ! y-direction point load
    dxforce     =  0.0d0                ! x-direction distributed load
    dyforce     =  0.0d0                ! y-direction distributed load
    thick       =  1.00d0               ! thickness
    x_width     =  1.0d0                ! x length
    y_width     =  1.0d0                ! y length
    n_i_node    =  domainsize + 1       ! the number of nodes in the x-direction
    n_j_node    =  domainsize + 1       ! the number of nodes in the y-direction
    n_i_element =  domainsize
    n_j_element =  domainsize
    del_x       =  x_width / (n_i_node - 1)
    del_y       =  y_width / (n_j_node - 1)

    node_n      = n_i_node * n_j_node
    element_n   = n_i_element * n_j_element

    allocate(Node(node_n))         ! allocate a node array 
    allocate(elem(element_n))   ! allocate a element array
    
    do i = 1, node_n
        Node(i).bc(:) = 0
        Node(i).eq_n(:) = 0
        Node(i).pm(:) = 0.0d0
    end do
     
    do i = 1, element_n
        elem(i).q(:) = 0.0d0
    end do
    
    ! set nodal position
    do i = 1, n_j_node
        do j = 1, n_i_node
            numbering = n_i_node * (j - 1) + i
            node(numbering).x(1) = del_x * (i - 1)
            node(numbering).x(2) = del_y * (j - 1)
        end do
    end do

    ! set connectivity
    do j = 1, n_j_element
        do i = 1, n_i_element
            numbering = n_i_element * (j - 1) + i

            elem(numbering).cn(1)	= n_i_node * (j - 1) + i
            elem(numbering).cn(2)	= n_i_node * (j - 1) + i + 1
            elem(numbering).cn(3)	= n_i_node * j + i + 1
            elem(numbering).cn(4)	= n_i_node * j + i
        end do
    end do

    ! imposing boundary condition
    x_fix_surf = 1								! left side B.C. : u, v = 1
	do j = 1, n_j_node
		numbering = n_i_node * (j - 1) + x_fix_surf
        
		Node(numbering).bc(1) = 1
		Node(numbering).bc(2) = 1
	end do

	! set body force
    elem(:).q(1) = bxforce
    elem(:).q(2) = byforce

    ! set thickness
    prop.thick = thick

	! set properties
	prop.young		= young
	prop.poisson	= possion

end subroutine Set2DQuadEleRectDomain

! --------------------------------------------------------------------------------

! Print computational cost
subroutine Compute_Cost(m_end, m_start, a_end, a_start, s_end, s_start)
    double precision, intent(in) :: m_end, m_start, a_end, a_start, s_end, s_start

    double precision :: cost_m, cost_a, cost_s

    cost_m = m_end - m_start
    cost_a = a_end - a_start
    cost_s = s_end - s_start

    write(0, "(a)")
    write(0, "(a)"), " [Computational cost]"
    write(0, "(a, a15, a)"), " 1. Total     : ", &
        trim(adjustl(Dble2Str((cost_m), 7, 2)))//"[sec], ", &
        trim(adjustl(Dble2Str((cost_m/60.0d0), 7, 2)))//"[min], "//&
        trim(adjustl(Dble2Str((cost_m/60.0d0/60.0d0), 6, 2)))//"[hr]"
    write(0, "(a, a15, a)"), " 2. Assembling: ", &
        trim(adjustl(Dble2Str((cost_a), 7, 2)))//"[sec], ", &
        trim(adjustl(Dble2Str((cost_a/60.0d0), 7, 2)))//"[min], "//&
        trim(adjustl(Dble2Str((cost_a/60.0d0/60.0d0), 6, 2)))//"[hr]"
    write(0, "(a, a15, a)"), " 3. Solving   : ", &
        trim(adjustl(Dble2Str((cost_s), 7, 2)))//"[sec], ", &
        trim(adjustl(Dble2Str((cost_s/60.0d0), 7, 2)))//"[min], "//&
        trim(adjustl(Dble2Str((cost_s/60.0d0/60.0d0), 6, 2)))//"[hr]"
    write(0, "(a)")
end subroutine Compute_Cost

! --------------------------------------------------------------------------------

end program SATURN_EDU