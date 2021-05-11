integer function pair_index(i, j)
implicit none
integer i, j, ii, jj
ii = max(i, j)
jj = min(i, j)

pair_index = ii * (ii + 1) / 2 - (ii - jj)

end function
