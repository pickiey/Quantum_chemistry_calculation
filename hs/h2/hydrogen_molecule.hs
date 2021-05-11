import Numeric.LinearAlgebra


num_of_atom, num_of_elec  :: Int
Nb, Nb2, Nb4              :: Int
pi, ddot                  :: Double
nucleus_repulsion         :: Double
contraction_length(4)     ::
xyz(3><2)                 :: Matrix Double
orbital_center            :: Vector Int
orbital_exponent(3><4)    :: Matrix Double
contr_coefficient(3><4)   :: Matrix Double
norm_constant(3><4)       :: Matrix Double
overlap(4><4)             :: Matrix Double
kinetic(4><4)             :: Matrix Double
nucattr(4><4)             :: Matrix Double
Hcore(4><4)               :: Matrix Double
molecular_orbital(4><4)   :: Matrix Double
Pmatrix(4><4)             :: Matrix Double
Fockmatrix(4><4)          :: Matrix Double
Jmatrix(4><4)             :: Matrix Double
Kmatrix(4><4)             :: Matrix Double
Wmatrix(4><4)             :: Matrix Double
orbital_energy(4)         ::
twoint(55)                ::
twoidx(55><4)             :: Matrix Int
Zcharge(2)                :: Double
i, j                      :: Int
tmp, work(15)             ::
lwork, info               :: Int
old_energy, new_energy    :: Double
SCFsuccess                :: Bool

num_of_atom = 2
num_of_elec = num_of_atom
SCFsuccess = false
Nb = 4
Nb2 = Nb * (Nb + 1) / 2
Nb4 = Nb2 * (Nb2 + 1) / 2
pi = dacos (-1.0e0)
lwork = 15
nucleus_repulsion = 1.0e0 / 1.322808e0
Zcharge(1) = 1.0e0
Zcharge(2) = 1.0e0
overlap = 0.0e0
kinetic = 0.0e0
nucattr = 0.0e0

-- set up xyz = geometry
xyz = 0.0d+0
xyz(3, 1) = 0.661404d+0
xyz(3, 2) = -0.661404d+0

-- set up orbital_center
orbital_center(1) = 1; orbital_center(2) = 1
orbital_center(3) = 2; orbital_center(4) = 2

-- set up contraction_length
contraction_length(1) = 3; contraction_length(2) = 1
contraction_length(3) = 3; contraction_length(4) = 1

-- set up orbital_exponent
orbital_exponent = 0.0d+0
orbital_exponent(1, 1) = 18.7311370e0
orbital_exponent(2, 1) = 2.8253937e0
orbital_exponent(3, 1) = 0.6401217e0
orbital_exponent(1, 2) = 0.1612778e0
orbital_exponent(:, 3) = orbital_exponent(:, 1)
orbital_exponent(:, 4) = orbital_exponent(:, 2)

-- set up contr_coefficient
contr_coefficient = 0.0e0
contr_coefficient(1, 1) = 0.03349460e0
contr_coefficient(2, 1) = 0.23472695e0
contr_coefficient(3, 1) = 0.81375733e0
contr_coefficient(1, 2) = 1.0e0
contr_coefficient(:, 3) = contr_coefficient(:, 1)
contr_coefficient(:, 4) = contr_coefficient(:, 2)


























main = do
  var <- 初期値
  hydmol var
  print $ 結果


hydmol :: a -> IO()
hydmol var = loop
  loop = do
    var更新処理
    loop
