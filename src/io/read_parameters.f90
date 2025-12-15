subroutine read_parameters(N,beta, &
  n_bins,mc_sweeps,alpha) &
  bind(c,name='read_parameters')

  use,intrinsic :: iso_c_binding,only:c_char,c_int,c_double,c_long
  implicit none

  integer(kind=c_int),    intent(inout) :: N

  real(kind=c_double),    intent(inout) :: beta
  integer(kind=c_int),    intent(inout) :: n_bins
  integer(kind=c_long),   intent(inout) :: mc_sweeps

  real(kind=c_double),    intent(inout) :: alpha

  namelist /Lattice/ N
  namelist /Simulation/ beta,n_bins,mc_sweeps
  namelist /Hamiltonian/ alpha

  open(unit=100,file='parameters',status='old')
  read(unit=100,nml=Lattice)
  read(unit=100,nml=Simulation)
  read(unit=100,nml=Hamiltonian)
  close(unit=100)
  
endsubroutine read_parameters

