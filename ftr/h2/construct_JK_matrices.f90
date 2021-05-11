subroutine construct_JK_matirices(Nb, Nb4, Jmatrix, Kmatrix, Pmatrix, &
& twoint, twoidx)
implicit none
integer :: Nb, Nb4
real(8) :: Jmatrix(Nb, Nb)
real(8) :: Kmatrix(Nb, Nb)
real(8) :: Pmatrix(Nb, Nb)
real(8) :: twoint(Nb4)
integer :: twoidx(Nb4, 4)
integer :: i, j, ij, c, d, cd, ijcd, pair_index
real(8) :: pref1, tmp1, tmp2

Jmatrix = 0.0d+0
Kmatrix = 0.0d+0

do ijcd = 1, Nb4
    i = twoidx(ijcd, 1)
    j = twoidx(ijcd, 2)
    ij = pair_index(i, j)
    c = twoidx(ijcd, 3)
    d = twoidx(ijcd, 4)
    cd = pair_index(c, d)
    pref1 = 1.0d+0
    if (i .eq. j) pref1 = pref1 * 5.0d-1
    if (c .eq. d) pref1 = pref1 * 5.0d-1
    if (ij .eq. cd) pref1 = pref1 * 5.0d-1
    tmp1 = pref1 * twoint(ijcd)
    tmp2 = tmp1 * 5.0d-1


    ! J matrix
    Jmatrix(i, j) = Jmatrix(i, j) + tmp1 * Pmatrix(c, d) * 2.0d+0
    Jmatrix(j, i) = Jmatrix(j, i) + tmp1 * Pmatrix(c, d) * 2.0d+0
    Jmatrix(c, d) = Jmatrix(c, d) + tmp1 * Pmatrix(i, j) * 2.0d+0
    Jmatrix(d, c) = Jmatrix(d, c) + tmp1 * Pmatrix(i, j) * 2.0d+0

    ! K matrix
    Kmatrix(i, c) = Kmatrix(i, c) + tmp2 * Pmatrix(j, d)
    Kmatrix(j, c) = Kmatrix(j, c) + tmp2 * Pmatrix(i, d)
    Kmatrix(i, d) = Kmatrix(i, d) + tmp2 * Pmatrix(j, c)
    Kmatrix(j, d) = Kmatrix(j, d) + tmp2 * Pmatrix(i, c)
    Kmatrix(c, i) = Kmatrix(c, i) + tmp2 * Pmatrix(d, j)
    Kmatrix(c, j) = Kmatrix(c, j) + tmp2 * Pmatrix(d, i)
    Kmatrix(d, i) = Kmatrix(d, i) + tmp2 * Pmatrix(c, j)
    Kmatrix(d, j) = Kmatrix(d, j) + tmp2 * Pmatrix(c, i)

enddo

end subroutine
