subroutine read_parameters(Lx,Ly,beta, &
  n_bins,mc_sweeps,Kx,Ky,Kz) &
  bind(c,name='read_parameters')

  use,intrinsic :: iso_c_binding,only:c_char,c_int,c_double,c_long
  implicit none

  integer(kind=c_int),    intent(inout) :: Lx
  integer(kind=c_int),    intent(inout) :: Ly

  real(kind=c_double),    intent(inout) :: beta
  integer(kind=c_int),    intent(inout) :: n_bins
  integer(kind=c_long),   intent(inout) :: mc_sweeps

  real(kind=c_double),    intent(inout) :: Kx
  real(kind=c_double),    intent(inout) :: Ky
  real(kind=c_double),    intent(inout) :: Kz

  namelist /Lattice/ Lx,Ly
  namelist /Simulation/ beta,n_bins,mc_sweeps
  namelist /Hamiltonian/ Kx,Ky,Kz

  open(unit=100,file='parameters',status='old')
  read(unit=100,nml=Lattice)
  read(unit=100,nml=Simulation)
  read(unit=100,nml=Hamiltonian)
  close(unit=100)
  
endsubroutine read_parameters

