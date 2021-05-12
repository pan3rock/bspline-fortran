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

    logical :: extrap = .false.
    call db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)
end subroutine db1val_wrapper

subroutine dbvalu_wrapper(t,a,n,k,ideriv,x,inbv,work,iflag,val) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : dbvalu
    implicit none
    integer(c_int), intent(in) :: n
    integer(c_int), intent(in) :: k
    real(c_double), intent(in) :: t(n + k)
    real(c_double), intent(in) :: a(n)
    integer(c_int), intent(in) :: ideriv
    real(c_double), intent(in) :: x
    integer(c_int), intent(out) :: inbv
    real(c_double), intent(inout) :: work(3 * k)
    integer(c_int), intent(out) :: iflag
    real(c_double), intent(out) :: val

    logical :: extrap = .false.
    call dbvalu(t,a,n,k,ideriv,x,inbv,work,iflag,val,extrap)
end subroutine dbvalu_wrapper


subroutine db1spvn_wrapper(t,n,jhigh,k,x,vnikx,ileft,iflag) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : dbspvn, dintrv
    use bspline_kinds_module, only: wp, ip
    integer(kind=c_int), intent(in) :: n
    integer(kind=c_int), intent(in) :: k
    real(kind=c_double), intent(in) :: t(n + k)
    integer(kind=c_int), intent(in) :: jhigh
    real(kind=c_double), intent(in) :: x
    real(kind=c_double), intent(out) :: vnikx(k)
    integer(kind=c_int), intent(out) :: ileft
    integer(kind=c_int), intent(out) :: iflag

    real(kind=c_double) :: xt
    integer(kind=c_int) :: mflag
    integer(kind=c_int) :: inbv = 1
    integer(kind=c_int) :: index = 1
    real(kind=c_double), dimension(2 * k) :: work
    integer(kind=c_int) :: iwork

    logical :: extrap = .false.
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
end subroutine db1spvn_wrapper


subroutine db1fqad_wrapper(cproc,tx,bcoef,nx,kx,idx,x1,x2,tol,f,iflag,w0) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : db1fqad, b1fqad_func
    type(c_funptr), intent(in), value :: cproc
    integer(c_int), intent(in) :: nx
    integer(c_int), intent(in) :: kx
    real(c_double), intent(in) :: tx(nx + kx)
    real(c_double), intent(in) :: bcoef(nx)
    integer(c_int), intent(in) :: idx
    real(c_double), intent(in) :: x1
    real(c_double), intent(in) :: x2
    real(c_double), intent(in) :: tol
    real(c_double), intent(out) :: f
    integer(c_int), intent(out) :: iflag
    real(c_double), intent(inout) :: w0(3 * kx)

    procedure(b1fqad_func), pointer :: proc
    call c_f_procpointer(cproc, proc)

    call db1fqad(proc,tx,bcoef,nx,kx,idx,x1,x2,tol,f,iflag,w0)

end subroutine db1fqad_wrapper


subroutine dintrv_wrapper(xt,lxt,xx,ilo,ileft,mflag) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : dintrv
    implicit none

    integer(c_int),intent(in) :: lxt    
    real(c_double), intent(in) :: xt(lxt)     
    real(c_double),intent(in) :: xx     
    integer(c_int),intent(inout) :: ilo    
    integer(c_int),intent(out) :: ileft  
    integer(c_int),intent(out) :: mflag  

    logical :: extrap = .false.
    call dintrv(xt,lxt,xx,ilo,ileft,mflag,extrap)                                             
end subroutine dintrv_wrapper

subroutine dbknot_wrapper(x,n,k,t) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : dbknot
    implicit none

    integer(c_int),intent(in) :: n
    integer(c_int),intent(in) :: k
    real(c_double),intent(in) :: x(n)
    real(c_double),intent(out) :: t(n + k)

    call dbknot(x,n,k,t)
end subroutine dbknot_wrapper

subroutine dbtpcf_wrapper(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag) bind(c)
    use iso_c_binding
    use bspline_sub_module, only : dbtpcf
    implicit none

    integer(c_int),intent(in) :: n
    integer(c_int),intent(in) :: nf
    integer(c_int),intent(in) :: ldf
    integer(c_int),intent(in) :: k
    real(c_double),intent(in) :: x(n)
    real(c_double),intent(in) :: fcn(ldf, nf)
    real(c_double),intent(in) :: t(n + k)
    real(c_double),intent(out) :: bcoef(nf, n)
    real(c_double),intent(out) :: work(2 * k * (n + 1))
    integer(c_int),intent(out) :: iflag

    call dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag)
end subroutine dbtpcf_wrapper