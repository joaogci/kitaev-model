subroutine read_parameters_replica(n_exchange_steps,n_exchange_frequency,n_replica_walkers) &
  bind(c,name='read_parameters_replica')

  use,intrinsic :: iso_c_binding,only:c_char,c_int,c_double,c_long
  implicit none

  integer(kind=c_int),    intent(inout) :: n_exchange_steps
  integer(kind=c_int),    intent(inout) :: n_exchange_frequency
  integer(kind=c_int),    intent(inout) :: n_replica_walkers

  namelist /Replica_Exchange/ n_exchange_steps,n_exchange_frequency,n_replica_walkers

  open(unit=100,file='parameters',status='old')
  read(unit=100,nml=Replica_Exchange)
  close(unit=100)
  
endsubroutine read_parameters_replica

