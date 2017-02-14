! Computes the matrix for the finite-element linear system
! Dr. John T. Hwang
! Sicheng He
! June, 2016

subroutine getmtx(num_nodes, num_elems, num_cons, nnz, &
     E, nodes, elems, areas, cons, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, num_cons, nnz, E, nodes, elems, areas, cons
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems, areas
  !f2py depend(num_cons) cons
  !f2py depend(nnz) data, rows, cols

  ! Input
  integer, intent(in) :: num_nodes, num_elems, num_cons, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: areas(num_elems)
  integer, intent(in) :: cons(num_cons)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: Telem(2, 6), Kelem(2, 2)
  double precision :: xyz1(3), xyz2(3), A, L, cos_xyz(3)
  double precision :: data6x6(6, 6)
  integer :: rows6x6(6, 6), cols6x6(6, 6)
  integer :: ones11(6, 6), ones12(6, 6), ones21(6, 6), ones22(6, 6)
  integer :: k, k1, k2, ielem, icons, index, ind1, ind2

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  do k2 = 1, 6
     do k1 = 1, 6
        rows6x6(k1, k2) = mod(k1-1, 3)
        cols6x6(k1, k2) = mod(k2-1, 3)
     end do
  end do

  ones11(:, :) = 0
  ones12(:, :) = 0
  ones21(:, :) = 0
  ones22(:, :) = 0

  ones11(1:3, 1:3) = 1
  ones12(1:3, 4:6) = 1
  ones21(4:6, 1:3) = 1
  ones22(4:6, 4:6) = 1

  Telem(:, :) = 0.

  data(:) = 0.
  rows(:) = 0
  cols(:) = 0

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))
     A = areas(ielem)

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:3) = cos_xyz
     Telem(2, 4:6) = cos_xyz

     data6x6(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data6x6 = data6x6 * E * A / L

     ind1 = 3 * (elems(ielem, 1)-1)
     ind2 = 3 * (elems(ielem, 2)-1)

     do k1 = 1, 6
        do k2 = 1, 6
           index = index + 1
           data(index) = data(index) + data6x6(k1, k2)
           rows(index) = rows(index) + rows6x6(k1, k2) + 1
           cols(index) = cols(index) + cols6x6(k1, k2) + 1

           rows(index) = rows(index) + ones11(k1, k2) * ind1
           cols(index) = cols(index) + ones11(k1, k2) * ind1

           rows(index) = rows(index) + ones12(k1, k2) * ind1
           cols(index) = cols(index) + ones12(k1, k2) * ind2

           rows(index) = rows(index) + ones21(k1, k2) * ind2
           cols(index) = cols(index) + ones21(k1, k2) * ind1

           rows(index) = rows(index) + ones22(k1, k2) * ind2
           cols(index) = cols(index) + ones22(k1, k2) * ind2
        end do
     end do
  end do

  do icons = 1, num_cons
     do k = 1, 3
        index = index + 1
        data(index) = 1.
        rows(index) = 3 * (cons(icons)-1) + k
        cols(index) = 3 * num_nodes + 3 * (icons-1) + k

        index = index + 1
        data(index) = 1.
        cols(index) = 3 * (cons(icons)-1) + k
        rows(index) = 3 * num_nodes + 3 * (icons-1) + k
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

end subroutine getmtx



subroutine getmtx2D(num_nodes, num_elems, num_cons, nnz, &
     E, nodes, elems, areas, cons, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, num_cons, nnz, E, nodes, elems, areas, cons
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems, areas
  !f2py depend(num_cons) cons
  !f2py depend(nnz) data, rows, cols

  ! Input
  integer, intent(in) :: num_nodes, num_elems, num_cons, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 2)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: areas(num_elems)
  integer, intent(in) :: cons(num_cons)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: Telem(2, 4), Kelem(2, 2)
  double precision :: xy1(2), xy2(2), A, L, cos_xy(2)
  double precision :: data4x4(4, 4)
  integer :: rows4x4(4, 4), cols4x4(4, 4)
  integer :: ones11(4, 4), ones12(4, 4), ones21(4, 4), ones22(4, 4)
  integer :: k, k1, k2, ielem, icons, index, ind1, ind2

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  do k2 = 1, 4
     do k1 = 1, 4
        rows4x4(k1, k2) = mod(k1-1, 2)
        cols4x4(k1, k2) = mod(k2-1, 2)
     end do
  end do

  ones11(:, :) = 0
  ones12(:, :) = 0
  ones21(:, :) = 0
  ones22(:, :) = 0

  ones11(1:2, 1:2) = 1
  ones12(1:2, 3:4) = 1
  ones21(3:4, 1:2) = 1
  ones22(3:4, 3:4) = 1

  Telem(:, :) = 0.

  data(:) = 0.
  rows(:) = 0
  cols(:) = 0

  index = 0
  do ielem = 1, num_elems
     xy1 = nodes(elems(ielem, 1), :)
     xy2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xy2 - xy1, xy2 - xy1))
     A = areas(ielem)

     cos_xy = (xy2 - xy1) / L
     Telem(1, 1:2) = cos_xy
     Telem(2, 3:4) = cos_xy

     data4x4(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data4x4 = data4x4 * E * A / L

     ind1 = 2 * (elems(ielem, 1)-1)
     ind2 = 2 * (elems(ielem, 2)-1)

     do k1 = 1, 4
        do k2 = 1, 4
           index = index + 1
           data(index) = data(index) + data4x4(k1, k2)
           rows(index) = rows(index) + rows4x4(k1, k2) + 1
           cols(index) = cols(index) + cols4x4(k1, k2) + 1

           rows(index) = rows(index) + ones11(k1, k2) * ind1
           cols(index) = cols(index) + ones11(k1, k2) * ind1

           rows(index) = rows(index) + ones12(k1, k2) * ind1
           cols(index) = cols(index) + ones12(k1, k2) * ind2

           rows(index) = rows(index) + ones21(k1, k2) * ind2
           cols(index) = cols(index) + ones21(k1, k2) * ind1

           rows(index) = rows(index) + ones22(k1, k2) * ind2
           cols(index) = cols(index) + ones22(k1, k2) * ind2
        end do
     end do
  end do

  do icons = 1, num_cons
     do k = 1, 2
        index = index + 1
        data(index) = 1.
        rows(index) = 2 * (cons(icons)-1) + k
        cols(index) = 2 * num_nodes + 2 * (icons-1) + k

        index = index + 1
        data(index) = 1.
        cols(index) = 2 * (cons(icons)-1) + k
        rows(index) = 2 * num_nodes + 2 * (icons-1) + k
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

  !print*, "data", data

end subroutine getmtx2D





subroutine getmtx2(num_nodes, num_elems, nnz, &
     E, nodes, elems, areas, data)

  implicit none

  !f2py intent(in) num_nodes, num_elems, nnz, E, nodes, elems, areas
  !f2py intent(out) data
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems, areas
  !f2py depend(nnz) data

  ! Input
  integer, intent(in) :: num_nodes, num_elems, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: areas(num_elems)

  ! Output
  double precision, intent(out) :: data(nnz)

  ! Working
  double precision :: Telem(2, 6), Kelem(2, 2)
  double precision :: xyz1(3), xyz2(3), A, L, cos_xyz(3)
  double precision :: data6x6(6, 6)
  integer :: k1, k2, ielem, index

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  Telem(:, :) = 0.

  data(36 * num_elems+1 : nnz) = 1. ! NOTICE: this may cause a possible bug: it will give the right boundary condition but the K matrix is disturbed!

  !data(:) = 1.
  !print*, "data(1:2)", data(1:2), "data(:)", data(:)

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))
     A = areas(ielem)

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:3) = cos_xyz
     Telem(2, 4:6) = cos_xyz

     data6x6(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data6x6 = data6x6 * E * A / L

     do k1 = 1, 6
        do k2 = 1, 6
           index = index + 1
           data(index) = data(index) + data6x6(k1, k2)
        end do
     end do
  end do

  if (index .ne. 36 * num_elems) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  !print*, "data", data, "E, A, L", E, A, L

end subroutine getmtx2









subroutine getmtx2D2(num_nodes, num_elems, nnz, &
     E, nodes, elems, areas, data)

  implicit none

  !f2py intent(in) num_nodes, num_elems, nnz, E, nodes, elems, areas
  !f2py intent(out) data
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems, areas
  !f2py depend(nnz) data

  ! Input
  integer, intent(in) :: num_nodes, num_elems, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 2)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: areas(num_elems)

  ! Output
  double precision, intent(out) :: data(nnz)

  ! Working
  double precision :: Telem(2, 4), Kelem(2, 2)
  double precision :: xyz1(2), xyz2(2), A, L, cos_xyz(2)
  double precision :: data4x4(4, 4)
  integer :: k1, k2, ielem, index

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  Telem(:, :) = 0.

  data(16 * num_elems+1 : nnz) = 1. ! NOTICE: this may cause a possible bug: it will give the right boundary condition but the K matrix is disturbed!

  !data(:) = 1.
  !print*, "data(1:2)", data(1:2), "data(:)", data(:)

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))
     A = areas(ielem)

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:2) = cos_xyz
     Telem(2, 3:4) = cos_xyz

     data4x4(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data4x4 = data4x4 * E * A / L

     do k1 = 1, 4
        do k2 = 1, 4
           index = index + 1
           data(index) = data(index) + data4x4(k1, k2)
        end do
     end do
  end do

  if (index .ne. 16 * num_elems) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  !print*, "data", data, "E, A, L", E, A, L

end subroutine getmtx2D2


















subroutine getresder(num_nodes, num_elems, num_aug, nnz, &
     E, nodes, elems, disp_aug, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, num_aug, nnz, E, nodes, elems, disp_aug
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems
  !f2py depend(num_aug) disp_aug
  !f2py depend(nnz) data, rows, cols

  ! Input
  integer, intent(in) :: num_nodes, num_elems, num_aug, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 3)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: disp_aug(num_aug)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: Telem(2, 6), Kelem(2, 2)
  double precision :: xyz1(3), xyz2(3), L, cos_xyz(3)
  double precision :: data6x6(6, 6)
  integer :: rows6x6(6, 6), cols6x6(6, 6)
  integer :: ones11(6, 6), ones12(6, 6), ones21(6, 6), ones22(6, 6)
  integer :: k1, k2, ielem, index, ind1, ind2, ind_aug

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  do k2 = 1, 6
     do k1 = 1, 6
        rows6x6(k1, k2) = mod(k1-1, 3)
        cols6x6(k1, k2) = mod(k2-1, 3)
     end do
  end do

  ones11(:, :) = 0
  ones12(:, :) = 0
  ones21(:, :) = 0
  ones22(:, :) = 0

  ones11(1:3, 1:3) = 1
  ones12(1:3, 4:6) = 1
  ones21(4:6, 1:3) = 1
  ones22(4:6, 4:6) = 1

  Telem(:, :) = 0.

  data(:) = 0.
  rows(:) = 0

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:3) = cos_xyz
     Telem(2, 4:6) = cos_xyz

     data6x6(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data6x6 = data6x6 * E / L

     ind1 = 3 * (elems(ielem, 1)-1)
     ind2 = 3 * (elems(ielem, 2)-1)

     do k1 = 1, 6
        do k2 = 1, 6
           index = index + 1

           ind_aug = 0
           data(index) = data(index) + data6x6(k1, k2)
           rows(index) = rows(index) + rows6x6(k1, k2) + 1
           ind_aug = ind_aug + cols6x6(k1, k2) + 1

           rows(index) = rows(index) + ones11(k1, k2) * ind1
           ind_aug = ind_aug + ones11(k1, k2) * ind1

           rows(index) = rows(index) + ones12(k1, k2) * ind1
           ind_aug = ind_aug + ones12(k1, k2) * ind2

           rows(index) = rows(index) + ones21(k1, k2) * ind2
           ind_aug = ind_aug + ones21(k1, k2) * ind1

           rows(index) = rows(index) + ones22(k1, k2) * ind2
           ind_aug = ind_aug + ones22(k1, k2) * ind2

           data(index) = data(index) * disp_aug(ind_aug)

           cols(index) = ielem
        end do
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

end subroutine getresder




subroutine getresder2D(num_nodes, num_elems, num_aug, nnz, &
     E, nodes, elems, disp_aug, data, rows, cols)

  implicit none

  !f2py intent(in) num_nodes, num_elems, num_aug, nnz, E, nodes, elems, disp_aug
  !f2py intent(out) data, rows, cols
  !f2py depend(num_nodes) nodes
  !f2py depend(num_elems) elems
  !f2py depend(num_aug) disp_aug
  !f2py depend(nnz) data, rows, cols

  ! Input
  integer, intent(in) :: num_nodes, num_elems, num_aug, nnz
  double precision, intent(in) :: E, nodes(num_nodes, 2)
  integer, intent(in) :: elems(num_elems, 2)
  double precision, intent(in) :: disp_aug(num_aug)

  ! Output
  double precision, intent(out) :: data(nnz)
  integer, intent(out) :: rows(nnz), cols(nnz)

  ! Working
  double precision :: Telem(2, 4), Kelem(2, 2)
  double precision :: xyz1(2), xyz2(2), L, cos_xyz(2)
  double precision :: data4x4(4, 4)
  integer :: rows4x4(4, 4), cols4x4(4, 4)
  integer :: ones11(4, 4), ones12(4, 4), ones21(4, 4), ones22(4, 4)
  integer :: k1, k2, ielem, index, ind1, ind2, ind_aug

  Kelem(1, 1) =  1.
  Kelem(1, 2) = -1.
  Kelem(2, 1) = -1.
  Kelem(2, 2) =  1.

  do k2 = 1, 4
     do k1 = 1, 4
        rows4x4(k1, k2) = mod(k1-1, 2)
        cols4x4(k1, k2) = mod(k2-1, 2)
     end do
  end do

  ones11(:, :) = 0
  ones12(:, :) = 0
  ones21(:, :) = 0
  ones22(:, :) = 0

  ones11(1:2, 1:2) = 1
  ones12(1:2, 3:4) = 1
  ones21(3:4, 1:2) = 1
  ones22(3:4, 3:4) = 1

  Telem(:, :) = 0.

  data(:) = 0.
  rows(:) = 0

  index = 0
  do ielem = 1, num_elems
     xyz1 = nodes(elems(ielem, 1), :)
     xyz2 = nodes(elems(ielem, 2), :)
     L = sqrt(dot_product(xyz2 - xyz1, xyz2 - xyz1))

     cos_xyz = (xyz2 - xyz1) / L
     Telem(1, 1:2) = cos_xyz
     Telem(2, 3:4) = cos_xyz

     data4x4(:, :) = matmul(matmul(transpose(Telem), Kelem), Telem)
     data4x4 = data4x4 * E / L

     ind1 = 2 * (elems(ielem, 1)-1)
     ind2 = 2 * (elems(ielem, 2)-1)

     do k1 = 1, 4
        do k2 = 1, 4
           index = index + 1

           ind_aug = 0
           data(index) = data(index) + data4x4(k1, k2)
           rows(index) = rows(index) + rows4x4(k1, k2) + 1
           ind_aug = ind_aug + cols4x4(k1, k2) + 1

           rows(index) = rows(index) + ones11(k1, k2) * ind1
           ind_aug = ind_aug + ones11(k1, k2) * ind1

           rows(index) = rows(index) + ones12(k1, k2) * ind1
           ind_aug = ind_aug + ones12(k1, k2) * ind2

           rows(index) = rows(index) + ones21(k1, k2) * ind2
           ind_aug = ind_aug + ones21(k1, k2) * ind1

           rows(index) = rows(index) + ones22(k1, k2) * ind2
           ind_aug = ind_aug + ones22(k1, k2) * ind2

           data(index) = data(index) * disp_aug(ind_aug)

           cols(index) = ielem
        end do
     end do
  end do

  if (index .ne. nnz) then
     print *, 'Error in getmtx: did not reach end of nnz vectors'
  end if

  rows(:) = rows(:) - 1
  cols(:) = cols(:) - 1

end subroutine getresder2D

