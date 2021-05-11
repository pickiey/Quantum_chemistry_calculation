program FD_polarizability_h2
implicit none

integer :: num_of_atom, num_of_elec
integer :: nb, nb2, nb4
real(8) :: pi, ddot
real(8) :: nucleus_repulsion
integer :: contraction_length(4)
real(8) :: xyz(3, 2)
integer :: orbital_center(4)
real(8) :: orbital_exponent(3, 4)
real(8) :: contr_coefficient(3, 4)
real(8) :: norm_constant(3, 4)
real(8) :: overlap(4, 4)
real(8) :: kinetic(4, 4)
real(8) :: nucattr(4, 4)
real(8) :: dipole(4, 4, 3)
real(8) :: hcore(4, 4)
real(8) :: molecular_orbital(4, 4)
real(8) :: pmatrix(4, 4)
real(8) :: fockmatrix(4, 4)
real(8) :: jmatrix(4, 4)
real(8) :: kmatrix(4, 4)
real(8) :: wmatrix(4, 4)
real(8) :: orbital_energy(4)
real(8) :: twoint(55)
integer :: twoidx(55, 4)
real(8) :: zcharge(2)
integer :: i, j
real(8) :: tmp, work(15)
integer :: lwork, info
real(8) :: old_energy, new_energy
logical :: scfsuccess
real(8) :: d_epsilon

num_of_atom = 2
num_of_elec = num_of_atom
scfsuccess = .false.
nb = 4
nb2 = nb * (nb + 1) / 2
nb4 = nb2 * (nb2 + 1) / 2
pi = dacos(-1.0d+0)
lwork = 15
nucleus_repulsion = 1.0d+0 / 1.322808d+0
zcharge(1) = 1.0d+0
zcharge(2) = 1.0d+0
overlap = 0.0d+0
kinetic = 0.0d+0
nucattr = 0.0d+0
dipole  = 0.0d+0

write(*, *) 'input d_epsilon'
read(*, *) d_epsilon

! set up xyz = geometry
xyz = 0.0d+0
xyz(3, 1) = 0.661404d+0
xyz(3, 2) = -0.661404d+0

! set up orbital_center
orbital_center(1) = 1; orbital_center(2) = 1
orbital_center(3) = 2; orbital_center(4) = 2

! set up contraction_length
contraction_length(1) = 3; contraction_length(2) = 1
contraction_length(3) = 3; contraction_length(4) = 1

! set up orbital_exponent
orbital_exponent = 0.0d+0
orbital_exponent(1, 1) = 18.7311370d+0
orbital_exponent(2, 1) = 2.8253937d+0
orbital_exponent(3, 1) = 0.6401217d+0
orbital_exponent(1, 2) = 0.1612778d+0
orbital_exponent(:, 3) = orbital_exponent(:, 1)
orbital_exponent(:, 4) = orbital_exponent(:, 2)

! set up contr_coefficient
contr_coefficient = 0.0d+0
contr_coefficient(1, 1) = 0.03349460d+0
contr_coefficient(2, 1) = 0.23472695d+0
contr_coefficient(3, 1) = 0.81375733d+0
contr_coefficient(1, 2) = 1.0d+0
contr_coefficient(:, 3) = contr_coefficient(:, 1)
contr_coefficient(:, 4) = contr_coefficient(:, 2)

! calculate norm_constant
do i = 1, nb
    do j = 1, contraction_length(i)
        norm_constant(j, i) = (2.0d+0 * orbital_exponent(j, i) / pi) ** 0.75d+0
    enddo
enddo

! calculate overlap, kinetic, nucattr
call calculate_one_integral(nb, 3, num_of_atom, &
& overlap, kinetic, nucattr, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz, zcharge)

! calculate dipole
call calculate_dipole_integral(nb, 3, num_of_atom, &
& dipole, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz, zcharge)

hcore = kinetic + nucattr - d_epsilon * dipole(:, :, 3)

! calculate coulomb repulsion integrals
call calculate_two_integral(nb, nb4, 3, &
& twoint, twoidx, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz)

! initialize density
pmatrix = 0.0d+0
pmatrix(1, 1) = 5.0d-1
pmatrix(2, 2) = 5.0d-1
pmatrix(3, 3) = 5.0d-1
pmatrix(4, 4) = 5.0d-1

! construct fock, j, and k matrices
call construct_jk_matirices(nb, nb4, jmatrix, kmatrix, &
& pmatrix, twoint, twoidx)
fockmatrix = hcore + jmatrix - kmatrix
molecular_orbital = fockmatrix
wmatrix = overlap
tmp = ddot(nb ** 2, hcore + fockmatrix, 1, pmatrix, 1)
new_energy = tmp * 5.0d-1 + nucleus_repulsion

! start scf calculation
do i = 1, 50 !!! up to 50 iterations

    ! write scf information
    write(6, '(a, i5)') ' iteration = ', i
    write(6, '(a, f20.10, a)') ' hartree-fock energy is ', &
    & new_energy, ' (hartree)'

    old_energy = new_energy

    ! solve f*c = e*s*c
    call dsygv(1, 'v', 'u', nb, molecular_orbital, nb, wmatrix, &
    & nb, orbital_energy, work, lwork, info)

    ! construct density matrix
    pmatrix = matmul(molecular_orbital(:, 1:num_of_elec/2) &
    & , transpose(molecular_orbital(:, 1:num_of_elec/2)))
    pmatrix = pmatrix * 2.0d+0

    ! construct fock, j, and k matrices
    call construct_jk_matirices(nb, nb4, jmatrix, kmatrix, &
    & pmatrix, twoint, twoidx)
    fockmatrix = hcore + jmatrix - kmatrix

    ! check brillouin's condition
    wmatrix = matmul(matmul(fockmatrix, pmatrix), overlap)
    wmatrix = wmatrix - matmul(matmul(overlap, pmatrix), fockmatrix)
    tmp = ddot(nb ** 2, wmatrix, 1, wmatrix, 1)
    write(6, '(a, f20.10)') ' residual = ', sqrt(tmp)

    molecular_orbital = fockmatrix
    wmatrix = overlap

    ! calculate hartree-fock energy
    new_energy = ddot(nb ** 2, hcore + fockmatrix, 1, pmatrix, 1) * 5.0d-1 &
    & + nucleus_repulsion

    if (dabs(old_energy - new_energy) .le. 1.0d-10 .and. tmp .le. 1.0d-10) then
        scfsuccess = .true.
        exit
    endif

enddo
! end of scf calculation

if (scfsuccess) then
    call dsygv(1, 'v', 'u', nb, molecular_orbital, nb, wmatrix, nb, &
    & orbital_energy, work, lwork, info)
    write(6, *)
    write(6, *) 'hartree-fock scf calculation converged.'
    write(6, '(a, f20.10, a)') ' hartree-fock energy is ', &
    & new_energy, ' (hartree)'
    write(6, '(a, f20.10, a)') ' dz is ', &
    & ddot(nb ** 2, pmatrix, 1, dipole(:, :, 3), 1), ' (a.u.)'
    write(6, '(a, f20.10, a)') ' alpha_zz is ', &
    & ddot(nb ** 2, pmatrix, 1, dipole(:, :, 3), 1) / d_epsilon, ' (a.u.)'
else
    write(6, *) 'hartree-fock scf calculation not converged.'
endif

! punch information
open(10, file='checkpoint.out')
write(10, *) num_of_elec, nb, nb * (nb + 1) * (nb ** 2 + nb * 2) / 8
write(10, *) fockmatrix
write(10, *) molecular_orbital
write(10, *) dipole
write(10, *) twoint
write(10, *) twoint
close(10)

end program
