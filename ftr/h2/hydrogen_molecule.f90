program hydrogen_molecule
implicit none

integer :: num_of_atom, num_of_elec
integer :: Nb, Nb2, Nb4
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
real(8) :: Hcore(4, 4)
real(8) :: molecular_orbital(4, 4)
real(8) :: Pmatrix(4, 4)
real(8) :: Fockmatrix(4, 4)
real(8) :: Jmatrix(4, 4)
real(8) :: Kmatrix(4, 4)
real(8) :: Wmatrix(4, 4)
real(8) :: orbital_energy(4)
real(8) :: twoint(55)
integer :: twoidx(55, 4)
real(8) :: Zcharge(2)
integer :: i, j
real(8) :: tmp, work(15)
integer :: lwork, info
real(8) :: old_energy, new_energy
logical :: SCFsuccess

num_of_atom = 2
num_of_elec = num_of_atom
SCFsuccess = .false.
Nb = 4
Nb2 = Nb * (Nb + 1) / 2
Nb4 = Nb2 * (Nb2 + 1) / 2
pi = dacos(-1.0d+0)
lwork = 15
nucleus_repulsion = 1.0d+0 / 1.322808d+0
Zcharge(1) = 1.0d+0
Zcharge(2) = 1.0d+0
overlap = 0.0d+0
kinetic = 0.0d+0
nucattr = 0.0d+0

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
do i = 1, Nb
    do j = 1, contraction_length(i)
        norm_constant(j, i) = (2.0d+0 * orbital_exponent(j, i) / pi) ** 0.75d+0
    enddo
enddo

! calculate overlap, kinetic, nucattr
call calculate_one_integral(Nb, 3, num_of_atom, &
& overlap, kinetic, nucattr, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz, Zcharge)

Hcore = kinetic + nucattr

! calculate coulomb repulsion integrals
call calculate_two_integral(Nb, Nb4, 3, &
& twoint, twoidx, orbital_exponent, contr_coefficient, &
& contraction_length, orbital_center, norm_constant, xyz)

! initialize density
Pmatrix = 0.0d+0
Pmatrix(1, 1) = 5.0d-1
Pmatrix(2, 2) = 5.0d-1
Pmatrix(3, 3) = 5.0d-1
Pmatrix(4, 4) = 5.0d-1

! construct Fock, J, and K matrices
call construct_JK_matirices(Nb, Nb4, Jmatrix, Kmatrix, &
& Pmatrix, twoint, twoidx)
Fockmatrix = Hcore + Jmatrix - Kmatrix
molecular_orbital = Fockmatrix
Wmatrix = overlap
tmp = ddot(Nb ** 2, Hcore + Fockmatrix, 1, Pmatrix, 1)
new_energy = tmp * 5.0d-1 + nucleus_repulsion

! start SCF calculation
do i = 1, 50 !!! up to 50 iterations

    ! write SCF information
    write(6, '(a, i5)') ' Iteration = ', i
    write(6, '(a, f20.10, a)') ' Hartree-Fock energy is ', &
    & new_energy, ' (Hartree)'

    old_energy = new_energy

    ! solve F*C = e*S*C
    call dsygv(1, 'V', 'U', Nb, molecular_orbital, Nb, Wmatrix, &
    & Nb, orbital_energy, work, lwork, info)

    ! construct density matrix
    Pmatrix = matmul(molecular_orbital(:, 1:num_of_elec/2) &
    & , transpose(molecular_orbital(:, 1:num_of_elec/2)))
    Pmatrix = Pmatrix * 2.0d+0

    ! construct Fock, J, and K matrices
    call construct_JK_matirices(Nb, Nb4, Jmatrix, Kmatrix, &
    & Pmatrix, twoint, twoidx)
    Fockmatrix = Hcore + Jmatrix - Kmatrix

    ! check Brillouin's condition
    Wmatrix = matmul(matmul(Fockmatrix, Pmatrix), overlap)
    Wmatrix = Wmatrix - matmul(matmul(overlap, Pmatrix), Fockmatrix)
    tmp = ddot(Nb ** 2, Wmatrix, 1, Wmatrix, 1)
    write(6, '(a, f20.10)') ' Residual = ', sqrt(tmp)

    molecular_orbital = Fockmatrix
    Wmatrix = overlap

    ! calculate Hartree-Fock energy
    new_energy = ddot(Nb ** 2, Hcore + Fockmatrix, 1, Pmatrix, 1) * 5.0d-1 &
    & + nucleus_repulsion

    if (dabs(old_energy - new_energy) .le. 1.0d-10 .and. tmp .le. 1.0d-10) then
        SCFsuccess = .true.
        exit
    endif

enddo
! end of SCF calculation

if (SCFsuccess) then
    call dsygv(1, 'V', 'U', Nb, molecular_orbital, Nb, Wmatrix, Nb, &
    & orbital_energy, work, lwork, info)
    write(6, *)
    write(6, *) 'Hartree-Fock SCF calculation converged.'
    write(6, '(a, f20.10, a)') ' Hartree-Fock energy is ', &
    & new_energy, ' (Hartree)'
else
    write(6, *) 'Hartree-Fock SCF calculation not converged.'
endif

! punch information
open(10, file='checkpoint.out')
write(10, *) num_of_elec, Nb, Nb * (Nb + 1) * (Nb ** 2 + Nb * 2) / 8
write(10, *) Fockmatrix
write(10, *) molecular_orbital
write(10, *) twoint
write(10, *) twoint
close(10)

end program
