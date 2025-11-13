program driver
  use ESMF
  use NUOPC
  !use nuopc_comp   ! <-- contains NUOPC_COMPONENT_ATM, _OCN, _CPL
  implicit none

  type(ESMF_GridComp) :: atm, ocn, connector
  integer :: rc

  ! Initialize ESMF
  call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, rc=rc)

  !------------------------------------------------------
  ! Create components
  !------------------------------------------------------
  atm       = ESMF_GridCompCreate(name="ATM", rc=rc)
  ocn       = ESMF_GridCompCreate(name="OCN", rc=rc)
  connector = ESMF_GridCompCreate(name="ATM-OCN Connector", rc=rc)

  !------------------------------------------------------
  ! Derive NUOPC component type
  !------------------------------------------------------
  call NUOPC_CompDerive(atm,       NUOPC_COMPONENT_ATM, rc=rc)
  call NUOPC_CompDerive(ocn,       NUOPC_COMPONENT_OCN, rc=rc)
  call NUOPC_CompDerive(connector, NUOPC_COMPONENT_CPL, rc=rc)

  !------------------------------------------------------
  ! Set services
  !------------------------------------------------------
  call NUOPC_CompSetServices(atm,       userRoutine=ATM_SetServices,       rc=rc)
  call NUOPC_CompSetServices(ocn,       userRoutine=OCN_SetServices,       rc=rc)
  call NUOPC_CompSetServices(connector, userRoutine=Connector_SetServices, rc=rc)

  !------------------------------------------------------
  ! Add components to driver and run
  !------------------------------------------------------
  call NUOPC_DriverAddComp(connector, atm, rc=rc)
  call NUOPC_DriverAddComp(connector, ocn, rc=rc)
  call NUOPC_DriverRun(connector, rc=rc)

  ! Finalize
  call ESMF_Finalize(rc=rc)
end program driver

