program main
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  integer,parameter :: nx = 512
  real(8),parameter :: Lx = 200d-9/0.529d-10, dx = Lx/nx
  real(8),parameter :: clight = 137.035999084d0
  complex(8) :: Ez(0:nx), dEz_dx(0:nx), zeps(0:nx), zeps_dx2(0:nx)
  complex(8) :: Ez_0, dEz_dx_0
  complex(8) :: zk1, zk2, zk3, zk4
  complex(8) :: zl1, zl2, zl3, zl4
  real(8) :: omega, kx
  integer :: ix


  zeps = 1d0
  zeps_dx2 = 1d0
  omega = 1d0
  kx = omega/clight

  Ez(nx) = 1d0
  dEz_dx(nx) = zi*kx

  do ix = 0, nx-1
     Ez_0  = Ez(nx-ix)
     dEz_dx_0 = dEz_dx(nx-ix)

     zk1 = dEz_dx_0
     zl1 = -omega**2/clight**2*zeps(nx-ix)*Ez_0

     zk2 = dEz_dx_0 -0.5d0*dx*zl1
     zl2 = -omega**2/clight**2*zeps_dx2(nx-ix)*(Ez_0-0.5d0*dx*zk1)

     zk3 = dEz_dx_0 -0.5d0*dx*zl2
     zl3 = -omega**2/clight**2*zeps_dx2(nx-ix)*(Ez_0-0.5d0*dx*zk2)

     zk4 = dEz_dx_0 -1d0*dx*zl3
     zl4 = -omega**2/clight**2*zeps(nx-ix-1)*(Ez_0-1d0*dx*zk3)

     Ez(nx-ix-1) = Ez_0 - dx/6d0*(zk1+2d0*zk2+2d0*zk3+zk4)
     dEz_dx(nx-ix-1) = dEz_dx_0 - dx/6d0*(zl1+2d0*zl2+2d0*zl3+zl4)     
     
  end do


  open(20,file="zEz.out")
  do ix = 0, nx
     write(20,"(99e26.16e3)")dx*ix,Ez(ix)
  end do
  close(20)

  
end program main
