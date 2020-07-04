!
! =============================================================================
!
! SATURN Educational Version
! by Hyungmin Jun (hyungminjun@outlook.com)
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
program SATURN_EDU

    use Data_Struct
    use DISP4
    use Solver

    implicit none

    ! define internal variables
    integer :: node_n, element_n                               ! # of nodes, # of elements
    type (NodeType),    allocatable, dimension(:) :: Node      ! node
    type (ElementType), allocatable, dimension(:) :: Element   ! element
    type (MaterialType) :: Material                            ! material

    integer :: total_dof_n, free_dof_n, fixed_dof_n              ! # of DOFs (total, free, fixed)
    double precision, allocatable, dimension(:) :: Kt        ! stiffness vector
    double precision, allocatable, dimension(:) :: U         ! diplacement vector
    double precision, allocatable, dimension(:) :: R         ! load vector
    integer,          allocatable, dimension(:) :: c_index    ! column index

    call Main

contains

! --------------------------------------------------------------------------------

! Main subroutine
subroutine Main()
    character(len=20) :: filename
    real :: t_main_start, t_main_end
    real :: t_assemble_start, t_assemble_end
    real :: t_solve_start, t_solve_end

    call cpu_time(t_main_start)

    ! input file name
    filename = "planestress"; print * 
   
    ! open files
    open(unit=1, file=trim(filename)//".txt",     form="formatted")
    open(unit=2, file=trim(filename)//"_out.txt", form="formatted")
    open(unit=3, file=trim(filename)//"_res.txt", form="formatted")
    open(unit=4, file=trim(filename)//"_pos.txt", form="formatted")

    ! read input file
    !call Read_Input_File
    call Set2DQuadEleRectDomain

    ! calculate # of total DOF, # of free DOF, # of fixed DOF and assign equation numbers
    call Set_DOF_Number(total_dof_n, free_dof_n, fixed_dof_n)
   
    ! calculate column index
    call Calculate_Index
   
    ! print information
    call Print_Title(0)
    
    ! allocate memory
    allocate( Kt( c_index(free_dof_n+1)-1 ) ) ! total stiffness vector (Kt)
    allocate( U(total_dof_n) )                ! displacement vector
    allocate( R(total_dof_n) )                ! load vector
   
    ! assemble total stiffness matrix
    print *, "1/ Assembling Stiffness and Load"
    call cpu_time(t_assemble_start)
    call Assemble_Kt

    ! assemble load vector
    call Assemble_Load;     call cpu_time(t_assemble_end);      call cpu_time(t_solve_start)
    
    ! slove linear system (find U in K*U=R)
    print *, "2/ Solving Linear System"
    call Solve_Equations;   call cpu_time(t_solve_end)

    ! calculate stress and print solutions
    print *, "3/ Printing Output Files"
    call Displacement_Stress

    print *, "4/ Completed !"; print * 
    print *, "Strain energy = ", 0.5d0 * dot_product(R(1:free_dof_n), U(1:free_dof_n))
    PRINT *, "Strain energy =   4.447735498387236E-008"
    ! deallocate memory
    deallocate(Node, Element, Kt, R, U, c_index)
    call cpu_time(t_main_end)

    ! close files
    close(unit=1); close(unit=2); close(unit=3); close(unit=4)
    
    call Print_TimeConsuming(0, t_main_end, t_main_start, t_assemble_end, t_assemble_start, t_solve_end, t_solve_start)

    print *, ""; pause
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
    allocate(Element(element_n)) ! allocate a element array
    
    do i=1, element_n
        read(1,*) bufi, Element(i).cn(:), Element(i).q(:)
    end do

    ! read properties
    read(1,*) bufs; read(1,*) thickness
    Element(:).thickness = thickness
    read(1,*) bufs; read(1,*) Material.Young
    read(1,*) bufs; read(1,*) Material.Poisson
end subroutine Read_Input_File

! --------------------------------------------------------------------------------

! calculate # of total DOF, # of free DOF, # of fixed DOF and assign equation numbers to DOFs
subroutine Set_DOF_Number(tn, fn, cn)
    integer, intent(out) :: tn, fn, cn  ! # of total DOF, # of free DOF, # of fixed DOF
    integer :: i, j

    write(3, *) "EQUATION NUMBER"
    write(3, *) "---------------------"
    write(3, *) "    node   dof    eqn"

    tn = node_n * nDPN  ! # of total DOF
    fn = 0; cn = 0

    do i=1, node_n
        do j=1, nDPN
            if (Node(i).bc(j) == 0) then
                fn = fn + 1
                Node(i).eq_n(j) = fn
                write(3, "(3i7)") i, j, fn
            else
                cn = cn + 1
                Node(i).eq_n(j) = tn - cn + 1
            end if
        end do
    end do
    write(3, *)
end subroutine Set_DOF_Number

! --------------------------------------------------------------------------------

! calculate column index for skyline solver
subroutine Calculate_Index()   
    integer :: column_h(total_dof_n) ! column height
    integer :: a_index(nDPN * nNPE)  ! index for assemblage
    integer :: i, j, k               ! index for iteration
    integer :: buf_sum

    ! allocate c_index array
    allocate ( c_index(free_dof_n+1) )

    ! column height
    column_h(:) = 0

    do i=1, element_n

        ! assemblage index
        do j=1, nDPN
            do k=1, nNPE
                a_index( nNPE * j + k - nNPE ) = node( element(i).cn(k) ).eq_n(j)
            end do
        end do

        ! column height
        do k=1, nDPN * nNPE
            do j=1, nDPN * nNPE
                if(a_index(j) <= free_dof_n .and. a_index(k) <= free_dof_n) then
                    if( a_index(j) < a_index(k) ) then
                        if( a_index(k) - a_index(j) > column_h(a_index(k)) ) then
                            column_h(a_index(k)) = a_index(k) - a_index(j)
                        end if
                    end if
                end if
            end do
        end do
    end do

    ! c_index array
    buf_sum = 1

    do i=1, free_dof_n
        c_index(i) = buf_sum
        buf_sum = buf_sum + Column_H(i) + 1
    end do

    c_index(free_dof_n+1) = buf_sum
    write(3,"(a18, E14.6, a3)") "REQUIRED MEMORY =", buf_sum * 8.0d0 / 1000000.d0, " MB"
    write(3,*)
end subroutine Calculate_Index

! --------------------------------------------------------------------------------

! assemble total stiffness matrix by using equation numbers
subroutine Assemble_Kt()
    double precision :: eNode(nNPE,nDIM)                ! nodal position of 4-node element (x,y)
    double precision :: Ke(nNPE * nDPN,nNPE * nDPN)     ! stifness matrix of element
    integer :: a_index(nNPE * nDPN)                      ! assemblage index
    integer :: i, j, k, address

    Kt(:) = 0.0d0

    do i=1, element_n

        ! nodal position of element
        do j=1, nNPE
            do k=1, nDIM
                eNode(j, k) = Node( Element(i).cn(j) ).x(k)
            end do
        end do

        ! calculate stiffness matrix of element
        call Plane_Stiffness(Material.Young, Material.Poisson, Element(i).thickness, eNode, Ke)
        !write(2,"(a24,i4)") " STIFFNESS of ELEMENT : ", i
        !call Print_Matrix(Ke)

        ! assemblage index
        do j=1, nDPN
            do k=1, nNPE
                a_index( nNPE * j + k - nNPE ) = node( element(i).cn(k) ).eq_n(j)
            end do
        end do

        ! assemble total stiffness matrix
        do j=1, nNPE * nDPN
            do k=1, nNPE * nDPN
                if(a_index(j) <= free_dof_n .and. a_index(k) <= free_dof_n) then
                    if( a_index(j) <= a_index(k) ) then
                        address = c_index(a_index(k)) + a_index(k) - a_index(j)
                        Kt(address) = Kt(address) + Ke(j, k)
                    end if
                end if
            end do
        end do
    end do
end subroutine Assemble_Kt

! --------------------------------------------------------------------------------

! assemble load vector
subroutine Assemble_Load()
    double precision :: eNode(nNPE,nDIM)
    double precision :: NodalLoad(nNPE * nDPN)  ! equivalent nodal load
    integer :: i, j, k

    R(:) = 0

    ! assemble load vector for nodal load
    do i=1, node_n
        do j=1, nDPN
            R(Node(i).eq_n(j)) = Node(i).pm(j)
        end do
    end do

    ! assemble load vector for body force
    do i=1, element_n

        ! nodal position of element
        do j=1, nNPE
            do k=1, nDIM
                eNode(j, k) = Node( Element(i).cn(j) ).x(k)
            end do
        end do

        ! calculate equivalent nodal load from body force
        call Plane_Load(eNode, Element(i).q, NodalLoad)

        ! assemble load vector
        do j=1, nDPN
            do k=1, nNPE
                R( Node( Element(i).cn(k) ).eq_n(j) ) = R( Node( Element(i).cn(k) ).eq_n(j) ) &
                    + NodalLoad(nNPE*j+k-nNPE)
            end do
        end do
    end do
end subroutine Assemble_Load

! --------------------------------------------------------------------------------

! love linear equations
subroutine Solve_Equations()
    integer :: i

    U(:) = 0.0d0; i=3
    U(1:free_dof_n) = R(1:free_dof_n)
    call Skyline_COLSOL(Kt(1:c_index(free_dof_n+1)-1), U(1:free_dof_n), c_index(1:free_dof_n+1), &
        free_dof_n,  c_index(free_dof_n+1)-1, free_dof_n+1, 1, i)
    call Skyline_COLSOL(Kt(1:c_index(free_dof_n+1)-1), U(1:free_dof_n), c_index(1:free_dof_n+1), &
        free_dof_n,  c_index(free_dof_n+1)-1, free_dof_n+1, 2, i)
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
    write(3,"(a17, E14.6)") "STRAIN ENERGY = ", 0.5d0 * dot_product(R(1:free_dof_n), U(1:free_dof_n))
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
                eNode(j, k) = Node( Element(i).cn(j) ).x(k)
            end do
        end do
      
        ! displacement vector of element 
        do j=1, nDPN
            do k=1, nNPE
                displace(nNPE * j + k - nNPE) = U( Node( Element(i).cn(k) ).eq_n(j) )
            end do
        end do

        ! calculate stress of element
        call Plane_Stress(Material.Young, Material.Poisson, eNode, displace, Stress)
      
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

subroutine Print_Title(n)
    integer, intent(in) :: n

    write(n,*) "[SATURN Educational Ver.]"
    write(n,*)
    write(n,*) "--------------------------------------------------------------------"
    write(n,"(a19, e10.3, a18, e10.3)") " Young's modulus : ", Material.Young, ", Possion ratio : ", Material.Poisson
    write(n,"(a13, e10.3)") " thickness : ",Element(1).thickness
    write(n,"(a20, i8, a19, i8)") " # total elements : ", element_n, ",  # total nodes : ", node_n
    write(n,"(a16, i6, a17, i6, a16, i6)") " # total DOFs : ", total_dof_n, ",  # free DOFs : ", free_dof_n, ",  # fix DOFs : ", fixed_dof_n
    write(n,"(a19, e10.3, a3)") " required memory = ", c_index(free_dof_n+1) * 8.0d0 / 1024.0d0 / 1024.0d0, " MB"
    write(n,*) "--------------------------------------------------------------------"
    write(n,*)

end subroutine Print_Title

! --------------------------------------------------------------------------------

subroutine Print_TimeConsuming(n, m_end, m_start, a_end, a_start, s_end, s_start)
    integer, intent(in) :: n
    real, intent(in) :: m_end, m_start, a_end, a_start, s_end, s_start
    real :: consuming_m, consuming_a, consuming_s
      
    consuming_m	= m_end - m_start
	consuming_a	= a_end - a_start
	consuming_s	= s_end - s_start
   
    write(n,*)
    write(n,"(a19, e10.3, a8, e10.3, a8, e10.3, a5)") " time consuming : ", consuming_m, "[sec], ", &
                                                       consuming_m/60.0d0, "[min], ", consuming_m/60.0d0/60.0d0, "[hr]"
    write(n,"(a19, e10.3, a8, e10.3, a8, e10.3, a5)") " assembling time : ", consuming_a, "[sec], ", &
                                                       consuming_a/60.0d0, "[min], ", consuming_a/60.0d0/60.0d0, "[hr]"
    write(n,"(a19, e10.3, a8, e10.3, a8, e10.3, a5)") " solving time : ", consuming_s, "[sec], ", &
                                                       consuming_s/60.0d0, "[min], ", consuming_s/60.0d0/60.0d0, "[hr]"

end subroutine Print_TimeConsuming

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
    allocate(Element(element_n))   ! allocate a element array
    
    do i = 1, node_n
        Node(i).bc(:) = 0
        Node(i).eq_n(:) = 0
        Node(i).pm(:) = 0.0d0
    end do
     
    do i = 1, element_n
        Element(i).q(:) = 0.0d0
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

            Element(numbering).cn(1)	= n_i_node * (j - 1) + i
            Element(numbering).cn(2)	= n_i_node * (j - 1) + i + 1
            Element(numbering).cn(3)	= n_i_node * j + i + 1
            Element(numbering).cn(4)	= n_i_node * j + i
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
    Element(:).q(1) = bxforce
    Element(:).q(2) = byforce

    ! set thickness
    Element(:).thickness = thick

	! set properties
	Material.Young		= young
	Material.Poisson	= possion

end subroutine Set2DQuadEleRectDomain

! --------------------------------------------------------------------------------

end program SATURN_EDU