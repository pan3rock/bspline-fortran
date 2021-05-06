subroutine db1ink_wrapper(x, nx, fcn, kx, iknot, tx, bcoef, iflag) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : db1ink
    integer(c_int), intent(in) :: nx
    integer(c_int), intent(in) :: kx
    real(c_double), intent(in) :: x(nx)
    real(c_double), intent(in) :: fcn(nx)
    integer(c_int), intent(in) :: iknot
    real(c_double), intent(inout) :: tx(nx + kx)
    real(c_double), intent(out) :: bcoef(nx)
    integer(c_int), intent(out) :: iflag

    call db1ink(x, nx, fcn, kx, iknot, tx, bcoef, iflag)
end subroutine db1ink_wrapper

subroutine db1val_wrapper(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : db1val
    integer(kind=c_int), intent(in) :: idx
    integer(kind=c_int), intent(in) :: nx
    integer(kind=c_int), intent(in) :: kx
    real(kind=c_double), intent(in) :: xval
    real(kind=c_double), intent(in) :: tx(nx + kx)
    real(kind=c_double), intent(in) :: bcoef(nx)
    real(kind=c_double), intent(out) :: f
    integer(kind=c_int), intent(out) :: iflag
    integer(kind=c_int), intent(inout) :: inbvx
    real(kind=c_double), intent(inout) :: w0(3 * kx)

    logical :: extrap = .true.
    call db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

end subroutine db1val_wrapper
