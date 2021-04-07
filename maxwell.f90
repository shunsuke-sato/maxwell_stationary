module global_variables
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  integer,parameter :: nx = 512
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: Lx = 200d-9/0.529d-10, dx = Lx/nx
  real(8),parameter :: clight = 137.035999084d0
  real(8),parameter :: cm = 0.529177210903e-8
  real(8),parameter :: gamma = 2.0d0/ev
  integer,parameter :: nw = 128
  real(8),parameter :: wi = 20d0/ev, wf = 80d0/ev, dw = (wf-wi)/nw
  real(8),parameter :: volume = 112.0394d0

  integer,parameter :: nene = 4
  real(8) :: Eex_surface = 0d0 !(0.1d0/ev)/volume
  complex(8) :: zeps_w(0:nw,nene)
  real(8) :: Eene(nene), Temp(nene)


  real(8) :: Eex_x(0:nx),Temp_x(0:nx)
  complex(8) :: zeps_x(0:nx)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none
  integer :: iw
  real(8) :: ww, reflectivity, transmission
  real(8) :: ref_gs(0:nw), tra_gs(0:nw)
  real(8) :: ref_pp(0:nw), tra_pp(0:nw)



  call setup_dielectric_function

  Eex_surface = 0d0 !(0.1d0/ev)/volume
  call calc_energy_distribution

  open(30,file='ref_trans_gs.out')
  do iw = 0, nw
    ww = wi + dw*iw
    call setup_dielectric_function_in_real_space(iw)
    call calc_maxwell(iw, reflectivity, transmission)
    ref_gs(iw) = reflectivity
    tra_gs(iw) = transmission
    write(30,"(999e26.16e3)")ww,reflectivity,transmission
  end do
  close(30)


  Eex_surface = (0.5d0/ev)/volume
  call calc_energy_distribution

  open(30,file='ref_trans_pp.out')
  do iw = 0, nw
    ww = wi + dw*iw
    call setup_dielectric_function_in_real_space(iw)
    call calc_maxwell(iw, reflectivity, transmission)
    ref_pp(iw) = reflectivity
    tra_pp(iw) = transmission
    write(30,"(999e26.16e3)")ww,reflectivity,transmission
  end do
  close(30)
  

  open(30,file='transient_abs.out')
  do iw = 0, nw
    ww = wi + dw*iw
    write(30,"(999e26.16e3)")ww,-log(tra_pp(iw)/tra_gs(iw))
  end do
  close(30)

!  zeps = 1d0
!  zeps_dx2 = 1d0
!  omega = 1d0
!  kx = omega/clight
!
!  Ez(nx) = 1d0
!  dEz_dx(nx) = zi*kx
!
!  do ix = 0, nx-1
!     Ez_0  = Ez(nx-ix)
!     dEz_dx_0 = dEz_dx(nx-ix)
!
!     zk1 = dEz_dx_0
!     zl1 = -omega**2/clight**2*zeps(nx-ix)*Ez_0
!
!     zk2 = dEz_dx_0 -0.5d0*dx*zl1
!     zl2 = -omega**2/clight**2*zeps_dx2(nx-ix)*(Ez_0-0.5d0*dx*zk1)
!
!     zk3 = dEz_dx_0 -0.5d0*dx*zl2
!     zl3 = -omega**2/clight**2*zeps_dx2(nx-ix)*(Ez_0-0.5d0*dx*zk2)
!
!     zk4 = dEz_dx_0 -1d0*dx*zl3
!     zl4 = -omega**2/clight**2*zeps(nx-ix-1)*(Ez_0-1d0*dx*zk3)
!
!     Ez(nx-ix-1) = Ez_0 - dx/6d0*(zk1+2d0*zk2+2d0*zk3+zk4)
!     dEz_dx(nx-ix-1) = dEz_dx_0 - dx/6d0*(zl1+2d0*zl2+2d0*zl3+zl4)     
!     
!  end do
!
!
!  open(20,file="zEz.out")
!  do ix = 0, nx
!     write(20,"(99e26.16e3)")dx*ix,Ez(ix)
!  end do
!  close(20)

  
end program main
!-------------------------------------------------------------------------------
subroutine calc_maxwell(iw, reflectivity, transmission)
  use global_variables
  implicit none
  integer,intent(in) :: iw
  real(8),intent(out) :: reflectivity, transmission
  real(8) :: ww, k0
  complex(8) :: Ez, dEz_dx
  complex(8) :: Ez_0, dEz_dx_0
  complex(8) :: zeps_mid
  complex(8) :: zk1, zk2, zk3, zk4
  complex(8) :: zl1, zl2, zl3, zl4
  integer :: ix
  complex(8) :: za, zr, zt


  ww = wi + dw*iw
  k0 = ww/clight

  Ez = 1d0
  dEz_dx = zi*k0*Ez
  zt = Ez

  do ix = 0, nx-1
    zk1 = -k0**2*zeps_x(nx-ix)*Ez
    zl1 = dEz_dx

    zeps_mid = 0.5d0*(zeps_x(nx-ix)+zeps_x(nx-ix-1))
    zk2 = -k0**2*zeps_mid*(Ez-0.5d0*dx*zl1)
    zl2 = dEz_dx -0.5d0*dx*zk1

    zk3 = -k0**2*zeps_mid*(Ez-0.5d0*dx*zl2)
    zl3 = dEz_dx -0.5d0*dx*zk2

    zk4 = -k0**2*zeps_x(nx-ix-1)*(Ez-dx*zl3)
    zl4 = dEz_dx -dx*zk3

    dEz_dx = dEz_dx -dx/6d0*(zk1+2d0*zk2+2d0*zk3+zk4)
    Ez = Ez -dx/6d0*(zl1+2d0*zl2+2d0*zl3+zl4)

  end do

  za = 0.5d0*(Ez+dEz_dx/(zi*k0))
  zr = 0.5d0*(Ez-dEz_dx/(zi*k0))

  reflectivity = abs(zr/za)**2
  transmission = abs(zt/za)**2

end subroutine calc_maxwell
!-------------------------------------------------------------------------------
subroutine setup_dielectric_function
  use global_variables
  implicit none
  integer,parameter :: nt = 15503
  real(8) :: jt(0:nt,3), tt(0:nt), Act(0:nt,3), Act0
  real(8) :: ww, xx, window, dt
  integer :: it, iw, it_t
  complex(8) :: zsigma, zj, zE0, zeps
  character(256) :: cfilename(nene)
  real(8) :: ss
  integer :: iene

  Temp(1) = 300d0
  Eene(1) = -2.09855190d0/volume
  cfilename(1) ='../Te_300K/td.general/total_current'

  Temp(2) = 1000d0
  Eene(2) = -2.09838665d0/volume
  cfilename(2) ='../Te_1000K/td.general/total_current'

  Temp(3) = 2000d0
  Eene(3) = -2.09784157d0/volume
  cfilename(3) ='../Te_2000K/td.general/total_current'

  Temp(4) = 13000d0
  Eene(4) = -2.06769810d0/volume
  cfilename(4) ='../Te_13000K/td.general/total_current'
  
  ss = Eene(1)
  Eene = Eene -ss

  do iene = 1, nene

    open(20,file=trim(cfilename(iene)))
    read(20,*); read(20,*)
    read(20,*); read(20,*)
    do it = 0, nt
      read(20,*)it_t,tt(it), jt(it,1:3)
    end do
    dt = tt(1)-tt(0)

    Act0 = -6.851799983950d-01/137.035999679d0

    do iw = 0, nw
      ww = wi + dw*iw
      zj = 0d0
      do it = 0, nt
!        xx = (tt(it)-tt(0))/(tt(nt)-tt(0))
!        window = 1d0 -3d0*xx**2 + 2d0*xx**3
        window = exp(-gamma*(tt(it)-tt(0)))
        zj = zj + exp(zi*ww*(tt(it)-tt(0))-0.5d0*dt)*window*jt(it,3)
      end do
      zj = zj*dt/volume
      zE0 = Act0*exp(zi*ww*tt(0))
      zsigma = zj/zE0

      zeps_w(iw,iene) = 1d0 + 4d0*pi*zi*zsigma/ww

    end do

  end do


end subroutine setup_dielectric_function
!-------------------------------------------------------------------------------
subroutine calc_energy_distribution
  use global_variables
  implicit none
  complex(8) :: zeps_IR
  real(8) :: ww, k0, ss, xx
  complex(8) :: zk,zt,zalpha, zbeta, zE0
  integer,parameter :: Np = nene
  real(8) :: xn(0:Np-1),yn(0:Np-1),an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
  integer :: ix

  zeps_IR = -66.405d0 + zi*44.255d0
  ww = 1.6d0/ev
  k0 = ww/clight
  zk = k0 * sqrt(zeps_IR)
  zt = 2d0*k0*zk/(2d0*k0*zk*cos(zk*Lx) - zI*(k0**2+zk**2)*sin(zk*Lx))
  zalpha = -2d0*k0*(k0+zk)*exp(-zi*zk*Lx)/(&
    -(k0+zk)**2*exp(-zi*zk*Lx)+(k0-zk)**2*exp(zi*zk*Lx) )
  zbeta  = 2d0*k0*(k0-zk)*exp(zi*zk*Lx)/(&
    -(k0+zk)**2*exp(-zi*zk*Lx)+(k0-zk)**2*exp(zi*zk*Lx) )

  do ix = 0, nx
    xx = dx*ix
    zE0 = zalpha*exp(zi*zk*xx) + zbeta*exp(-zi*zk*xx)
    Eex_x(ix) = abs(zE0)**2
  end do

  ss = Eex_x(0)
  Eex_x = Eex_x*Eex_surface/ss

  xn(0:np-1) = Eene(1:np)
  yn(0:np-1) = Temp(1:np)
  call spline(Np,xn,yn,an,bn,cn,dn)
  

  do ix = 0, nx
    xx = dx*ix
    call calc_signal(Np,xn,yn,an,bn,cn,dn,Eex_x(ix),ss)
    Temp_x(ix) = ss
  end do

  open(30,file='Eex_Te_dist.out')
  do ix = 0, nx
    xx = dx*ix
    write(30,"(999e26.16e3)")xx, Eex_x(ix), Temp_x(ix)
  end do
  close(30)

end subroutine calc_energy_distribution
!-------------------------------------------------------------------------------
subroutine setup_dielectric_function_in_real_space(iw)
  use global_variables
  implicit none
  integer,intent(in) :: iw
  integer,parameter :: Np = nene
  real(8) :: xn(0:Np-1)
  real(8) :: yn_r(0:Np-1),an_r(0:Np-2),bn_r(0:Np-2),cn_r(0:Np-2),dn_r(0:Np-2)
  real(8) :: yn_i(0:Np-1),an_i(0:Np-2),bn_i(0:Np-2),cn_i(0:Np-2),dn_i(0:Np-2)
  integer :: iene, ix
  real(8) :: sr, si, xx


  do iene = 1, nene
    xn(iene-1)    = Eene(iene)
    yn_r(iene-1)  = real( zeps_w(iw,iene))
    yn_i(iene-1)  = aimag(zeps_w(iw,iene))
  end do

  call spline(Np,xn,yn_r,an_r,bn_r,cn_r,dn_r)
  call spline(Np,xn,yn_i,an_i,bn_i,cn_i,dn_i)

  zeps_x = 0d0

  do ix = 0, nx
    xx = dx*ix
    call calc_signal(Np,xn,yn_r,an_r,bn_r,cn_r,dn_r,Eex_x(ix),sr)
    call calc_signal(Np,xn,yn_i,an_i,bn_i,cn_i,dn_i,Eex_x(ix),si)
    zeps_x(ix) = sr + zi*si
  end do

end subroutine setup_dielectric_function_in_real_space
!-------------------------------------------------------------------------------
subroutine spline(Np,xn,yn,an,bn,cn,dn)
  implicit none
  integer,intent(in) :: Np
  real(8),intent(in) :: xn(0:Np-1),yn(0:Np-1)
  real(8),intent(out) :: an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
  integer :: i,j,k,Npm2
  real(8) :: dxn(0:Np-1),dyn(0:Np-1),u(1:Np-2),v(1:Np-2),Amat(1:Np-2,1:Np-2)
  real(8) :: Amat_t(1:Np-2,1:Np-2),Amat_chk(1:Np-2,1:Np-2)
  real(8) :: ss 
! for lapack
  integer :: LWORK
  integer, allocatable :: IPIV(:) ! dimension N
  real(8), allocatable :: WORK(:) ! dimension LWORK
  integer :: info

  Npm2 = Np-2
  LWORK = Npm2*Npm2*6
  allocate(IPIV(Npm2),WORK(Npm2))


  do i = 0,Np-2
    dxn(i) = xn(i+1) - xn(i)
    dyn(i) = yn(i+1) - yn(i)
  end do

  do i = 1,Npm2
    v(i) = 6d0*(dyn(i)/dxn(i) - dyn(i-1)/dxn(i-1))
  end do


  Amat = 0d0
  Amat(1,1) = 2d0*(dxn(1) + dxn(0))
  Amat(1,2) = dxn(1)
  do i = 2,Npm2-1
    Amat(i,i+1) = dxn(i)
    Amat(i,i) = 2d0*(dxn(i)+dxn(i-1))
    Amat(i,i-1) = dxn(i-1)
  end do
  Amat(Npm2,Npm2) = 2d0*(dxn(Npm2)+dxn(Npm2-1))
  Amat(Npm2,Npm2-1) = dxn(Npm2-1)

! inverse matrix problem
  Amat_t = Amat


  call DGETRF(Npm2, Npm2, Amat_t, Npm2, IPIV, info)  ! factorize
  call DGETRI(Npm2, Amat_t, Npm2, IPIV, WORK, LWORK, info)  ! inverse

!  check inverse matrix problem
!  do i = 1,Npm2
!    do j = 1,Npm2
!      ss = 0d0
!      do k = 1,Npm2
!        ss = ss + Amat(i,k)*Amat_t(k,j)
!      end do
!      Amat_chk(i,j) = ss
!    end do
!  end do
!
!  do i = 1,Npm2
!    write(*,'(999e16.6e3)')(Amat_chk(i,j),j=1,Npm2)
!  end do
!
!  stop

  do i = 1,Npm2
    u(i) = sum(Amat_t(i,:)*v(:))
!    write(*,*)u(i),v(i)
  end do
!  stop


! for b
  bn(0) = 0d0
  bn(1:Np-2) = 0.5d0*u(1:Np-2)
! for a
  do i = 0,Npm2-1
    an(i) = (u(i+1) -2d0*bn(i))/(6d0*dxn(i))
  end do
  an(Npm2) = (0d0 -2d0*bn(Npm2))/(6d0*dxn(Npm2))
! for d
  dn(0:Npm2) = yn(0:Npm2)
! for c
  cn(0) = dyn(0)/dxn(0) - dxn(0)*(u(0+1)+2d0*0d0)/6d0
  do i = 1,Npm2-1
    cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*u(i))/6d0
  end do
  cn(Npm2) = dyn(Npm2)/dxn(Npm2) - dxn(Npm2)*(0d0+2d0*u(Npm2))/6d0

  return
end subroutine spline
!-------------------------------------------------------------------------------
subroutine calc_signal(Np,xn,yn,an,bn,cn,dn,xx,signal)
  implicit none
  integer,intent(in) :: Np
  real(8),intent(in) :: xn(0:Np-1),yn(0:Np-1)
  real(8),intent(in) :: an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
  real(8),intent(in) :: xx
  real(8),intent(out) :: signal
  integer :: ip
  real(8) :: dx, xx_abs
  

  xx_abs = abs(xx)

  if(xx_abs < xn(0))then
    write(*,*)"Error: xx<xn(0)"
    stop
  end if

  if(xx_abs > xn(Np-1))then
    write(*,*)"Error: xx>xn(Np-1)"
    stop
  end if

  do ip = 1,Np-1
    if(xx_abs <= xn(ip))exit
  end do
  ip = ip - 1
  dx = xx_abs - xn(ip)

  signal = an(ip)*dx**3 + bn(ip)*dx**2 + cn(ip)*dx + dn(ip)
  signal = sign(signal,xx)

end subroutine calc_signal
