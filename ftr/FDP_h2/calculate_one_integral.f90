subroutine calculate_one_integral(Nb, max_cnt_lngth, num_of_atom, &
& overlap, kinetic, nucattr, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz, Zcharge)
implicit none
! Arguments
integer :: Nb, max_cnt_lngth, num_of_atom
real(8) :: overlap(Nb, Nb)
real(8) :: kinetic(Nb, Nb)
real(8) :: nucattr(Nb, Nb)
integer :: orbital_center(Nb)
real(8) :: orbital_exponent(max_cnt_lngth, Nb)
real(8) :: contr_coefficient(max_cnt_lngth, Nb)
real(8) :: Zcharge(num_of_atom)
integer :: contraction_length(Nb)
real(8) :: xyz(max_cnt_lngth, 2)
real(8) :: norm_constant(max_cnt_lngth, Nb)

! Local variables
integer :: i, j, k, l, m
real(8) :: xi, yi, zi, xj, yj, zj, xp, yp, zp, xpc, ypc, zpc
real(8) :: norm1, expo1, cont1, norm2, expo2, cont2, tmp1, tmp2, tmp3, gamma
real(8) :: pref1, pref2, pref3, pref4, R_AB2, R_PC2
real(8) :: pi, Fgamma

! calculate overlap, kinetic, nucattr

pi = dacos(-1.0d+0)

do i = 1, Nb
    xi = xyz(1, orbital_center(i))
    yi = xyz(2, orbital_center(i))
    zi = xyz(3, orbital_center(i))

    do k = 1, contraction_length(i)
        norm1 = norm_constant(k, i)
        expo1 = orbital_exponent(k, i)
        cont1 = contr_coefficient(k, i)

        do j = 1, Nb
            xj = xyz(1, orbital_center(j))
            yj = xyz(2, orbital_center(j))
            zj = xyz(3, orbital_center(j))

            R_AB2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            tmp1 = 0.0d+0
            tmp2 = 0.0d+0
            tmp3 = 0.0d+0

            do l = 1, contraction_length(j)
                norm2 = norm_constant(l, j)
                expo2 = orbital_exponent(l, j)
                cont2 = contr_coefficient(l, j)

                gamma = expo1 + expo2
                xp = (expo1 * xi + expo2 * xj) / gamma
                yp = (expo1 * yi + expo2 * yj) / gamma
                zp = (expo1 * zi + expo2 * zj) / gamma
                pref1 = (expo1 * expo2) / gamma
                pref2 = (pi / gamma) ** 1.5d+0 * dexp(-pref1 * R_AB2)
                pref3 = norm1 * norm2 * cont1 * cont2
                pref4 = 2.0d+0 * pi / gamma * dexp(-pref1 * R_AB2)

                ! overlap
                tmp1 = tmp1 + pref3 * pref2

                ! kinetic
                tmp2 = tmp2 + pref3 * pref2 * pref1 * (3.0d+0 - 2.0d+0 * pref1 * R_AB2)

                ! nucleus attraction
                do m = 1, num_of_atom
                    xpc = xp - xyz(1, m)
                    ypc = yp - xyz(2, m)
                    zpc = zp - xyz(3, m)
                    R_PC2 = xpc ** 2 + ypc **2 + zpc ** 2
                    tmp3 = tmp3 - pref3 * pref4 * Zcharge(m) * Fgamma(gamma * R_PC2)
                enddo

            enddo

            overlap(j, i) = overlap(j, i) + tmp1
            kinetic(j, i) = kinetic(j, i) + tmp2
            nucattr(j, i) = nucattr(j, i) + tmp3

        enddo
    enddo
enddo

end subroutine
