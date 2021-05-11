subroutine calculate_two_integral(Nb, Nb4, max_cnt_lngth, &
& twoint, twoidx, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz)
implicit none
! Arguments
integer :: Nb, Nb4, max_cnt_lngth
integer :: orbital_center(Nb)
real(8) :: orbital_exponent(max_cnt_lngth, Nb)
real(8) :: contr_coefficient(max_cnt_lngth, Nb)
integer :: contraction_length(Nb)
real(8) :: xyz(max_cnt_lngth, 2)
real(8) :: norm_constant(max_cnt_lngth, Nb)
real(8) :: twoint(Nb4)
integer :: twoidx(Nb4, 4)

! Local variables
integer :: i, j, k, l, c, d, e, f, ij, cd, ijcd
integer :: pair_index
real(8) :: xi, yi, zi, xj, yj, zj, xc, yc, zc, xd, yd, zd, xp, yp, zp, xq, yq, zq
real(8) :: norm1, expo1, cont1, norm2, expo2, cont2
real(8) :: norm3, expo3, cont3, norm4, expo4, cont4
real(8) :: pref1, pref2, pref3, pref4, pref5, R_AB2, R_CD2, R_PQ2
real(8) :: pi, gamma, delta, eta, twopi25, Fgamma

! calculate overlap, kinetic, nucattr

pi = dacos(-1.0d+0)
twopi25 = 2.0d+0 * pi ** 2 * sqrt(pi)

twoint = 0.0d+0
! 1st index
do i = 1, Nb
    xi = xyz(1, orbital_center(i))
    yi = xyz(2, orbital_center(i))
    zi = xyz(3, orbital_center(i))
    do k = 1, contraction_length(i)
        norm1 = norm_constant(k, i)
        expo1 = orbital_exponent(k, i)
        cont1 = contr_coefficient(k, i)

        ! 2nd index
        do j = 1, i
            ij = pair_index(i, j)
            xj = xyz(1, orbital_center(j))
            yj = xyz(2, orbital_center(j))
            zj = xyz(3, orbital_center(j))
            R_AB2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            do l = 1, contraction_length(j)
                norm2 = norm_constant(l, j)
                expo2 = orbital_exponent(l, j)
                cont2 = contr_coefficient(l, j)

                gamma = expo1 + expo2
                xp = (expo1 * xi + expo2 * xj) / gamma
                yp = (expo1 * yi + expo2 * yj) / gamma
                zp = (expo1 * zi + expo2 * zj) / gamma
                pref1 = (expo1 * expo2) / gamma
                pref2 = norm1 * norm2 * cont1 * cont2 * dexp(-pref1 * R_AB2)

                ! 3rd index
                do c = 1, Nb
                    xc = xyz(1, orbital_center(c))
                    yc = xyz(2, orbital_center(c))
                    zc = xyz(3, orbital_center(c))
                    do e = 1, contraction_length(c)
                        norm3 = norm_constant(e, c)
                        expo3 = orbital_exponent(e, c)
                        cont3 = contr_coefficient(e, c)

                        ! 4th index
                        do d = 1, c
                            cd = pair_index(c, d)
                            if (cd .gt. ij) cycle
                            ijcd = pair_index(ij, cd)
                            xd = xyz(1, orbital_center(d))
                            yd = xyz(2, orbital_center(d))
                            zd = xyz(3, orbital_center(d))
                            R_CD2 = (xc - xd) ** 2 + (yc - yd) ** 2 + (zc - zd) ** 2
                            do f = 1, contraction_length(d)
                                norm4 = norm_constant(f, d)
                                expo4 = orbital_exponent(f, d)
                                cont4 = contr_coefficient(f, d)

                                delta = expo3 + expo4
                                eta = gamma + delta
                                xq = (expo3 * xc + expo4 * xd) / delta
                                yq = (expo3 * yc + expo4 * yd) / delta
                                zq = (expo3 * zc + expo4 * zd) / delta
                                R_PQ2 = (xp - xq) ** 2 + (yp - yq) ** 2 + (zp - zq) ** 2
                                pref3 = (expo3 * expo4) / delta
                                pref4 = norm3 * norm4 * cont3 * cont4 * dexp(-pref3 * R_CD2)
                                pref5 = 1.0d+0 / (gamma * delta * sqrt(eta))

                                twoint(ijcd) = twoint(ijcd) + &
                                & twopi25 * pref2 * pref4 * pref5 * Fgamma(gamma * delta / eta * R_PQ2)

                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo
enddo

! setup coulomb repulsion integrals indices
do i = 1, Nb
    do j = 1, i
        ij = pair_index(i, j)
        do c = 1, Nb
            do d = 1, c
                cd = pair_index(c, d)
                if (cd .gt. ij) cycle
                ijcd = pair_index(ij, cd)
                twoidx(ijcd, 1) = i
                twoidx(ijcd, 2) = j
                twoidx(ijcd, 3) = c
                twoidx(ijcd, 4) = d
            enddo
        enddo
    enddo
enddo

end subroutine
