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
    real(kind=c_double), dimension(*), intent(in) :: tx
    real(kind=c_double), intent(in) :: bcoef(nx)
    real(kind=c_double), intent(out) :: f
    integer(kind=c_int), intent(out) :: iflag
    integer(kind=c_int), intent(inout) :: inbvx
    real(kind=c_double), intent(inout) :: w0(3 * kx)

    logical :: extrap = .true.
    call db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

end subroutine db1val_wrapper


subroutine dbspvn_wrapper(t,n,jhigh,k,x,vnikx,ileft,iflag)
    use iso_c_binding
    use bspline_sub_module, only : dbspvn, dintrv
    use bspline_kinds_module, only: wp, ip
    real(kind=c_double), intent(in) :: t(n + k)
    integer(kind=c_int), intent(in) :: n
    integer(kind=c_int), intent(in) :: jhigh
    integer(kind=c_int), intent(in) :: k
    real(kind=c_double), intent(in) :: x
    real(kind=c_double), intent(out) :: vnikx(k)
    integer(kind=c_int), intent(out) :: ileft
    integer(kind=c_int), intent(out) :: iflag

    logical :: extrap = .true.
    real(kind=c_double) :: xt
    integer(kind=c_int) :: mflag
    integer(kind=c_int) :: inbv = 1
    integer(kind=c_int) :: index = 1
    real(kind=c_double), dimension(2 * k) :: work
    integer(kind=c_int) :: iwork

    if (extrap) then
        if (x<t(1_ip)) then
            xt = t(1_ip)
        else if(x>t(n+k)) then
            xt = t(n+k)
        else
            xt = x
        end if
    else
        xt = x
    end if

    call dintrv(t, n+1, xt, inbv, ileft, mflag)
    if (xt<t(k)) then
        iflag = 404_ip  ! dbvalu - x is not greater than or equal to t(k)
        return
    end if

    if (mflag/=0_ip) then

        if (xt>t(ileft)) then
            iflag = 405_ip  ! dbvalu - x is not less than or equal to t(n+1)
            return
        end if

        do
            if (ileft==k) then
                iflag = 406_ip  ! dbvalu - a left limiting value cannot be obtained at t(k)
                return
            end if
            ileft = ileft - 1_ip
            if (xt/=t(ileft)) exit
        end do

    end if

    call dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork,iflag)
end subroutine dbspvn_wrapper