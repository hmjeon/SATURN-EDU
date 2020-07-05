!
! =============================================================================
!
! Solver - Skyline COLSOL
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
module Solver

   public Skyline_COLSOL

contains

! --------------------------------------------------------------------------------

! Skyline solver, a(n_nwk), x(nn), maxa(nnm)
subroutine Skyline_COLSOL(a, b, x, maxa, nn, n_wk, nnm, kkk, iout)
    double precision, intent(inout) :: a(:)
    double precision, intent(in)    :: b(:)
    double precision, intent(out)   :: x(:)
    integer, intent(in) :: maxa(:)
    integer, intent(in) :: nn
    integer, intent(in) :: n_wk
    integer, intent(in) :: nnm
    integer, intent(in) :: kkk
    integer, intent(in) :: iout

    double precision :: cc, c, d
    integer :: n, kn, kl, ku, kh, k, ic, klt, j, ki, nd, kk, l

    x = b

    ! Perform l*d*l(t) factorization of stiffness matrix
    if (kkk - 2) 40, 150, 150

40  do 140 n = 1, nn
        kn = maxa(n)
        kl = kn + 1
        ku = maxa(n + 1) - 1
        kh = ku - kl
        if (kh) 110, 90, 50
50      k = n - kh
        ic = 0
        klt = ku
        do 80 j = 1, kh
            ic = ic + 1
            klt = klt - 1
            ki = maxa(k)
            nd = maxa(k + 1) - ki - 1
            if (nd) 80, 80, 60
60          kk = min0(ic, nd)
            c = 0.0d0
            do 70 l = 1, kk
70          c = c + a(ki + l) * a(klt + l)
            a(klt) = a(klt) - c
80      k = k + 1
90      k = n
        d = 0.0d0
        do 100 kk = kl, ku
            k = k - 1
            ki = maxa(k)
            c = a(kk) / a(ki)
            d = d + c * a(kk)
100     a(kk) = c
        a(kn) = a(kn) - d
110     if (a(kn)) 120, 120, 140
120     write(iout, 2000) n, a(kn)
        write(*, 2000) n, a(kn)
        stop
        go to 800
140 continue
    go to 900

    ! Reduce right-hand-side load vector
150 do 180 n = 1, nn
        kl = maxa(n) + 1
        ku = maxa(n + 1) - 1
        if (ku - kl) 180, 160, 160
160     k = n
        c = 0.0d0
        do 170 kk = kl, ku
            k = k - 1
170     c = c + a(kk) * x(k)
        x(n) = x(n) - c
180 continue

    ! Back-substitute
    do 200 n = 1, nn
        k = maxa(n)
200 x(n) = x(n) / a(k)
    if (nn==1) go to 900
    n = nn
    do 230 l = 2, nn
        kl = maxa(n) + 1
        ku = maxa(n + 1) - 1
        if (ku - kl) 230, 210, 210
210     k = n
        do 220 kk = kl, ku
            k = k - 1
220     x(k) = x(k) - a(kk) * x(n)
230 n = n - 1
    go to 900
800 stop
900 return
2000 format (//' Stop - Stiffness matrix not positive definite',//, &
        ' Nonpositive pivot for equation ',i8,//, ' Pivot = ',e20.12 )
end subroutine Skyline_COLSOL

! --------------------------------------------------------------------------------

end module Solver