real(8) function Fgamma(x)
implicit none
real(8) x, tmp, pi

pi = dacos(-1.0d+0)
if (x .lt. 0.0d+0) call abort
tmp = x
if (x .lt. 1.0d-10) tmp = 1.0d-10
tmp = erf(dsqrt(tmp)) * 0.5d+0 * dsqrt(pi / tmp)
Fgamma = tmp

end function
