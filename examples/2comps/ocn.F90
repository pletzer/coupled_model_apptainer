subroutine OCN_SetServices(comp, rc)
  use ESMF
  use NUOPC
  implicit none
  type(ESMF_GridComp), intent(inout) :: comp
  integer, intent(out) :: rc

  rc = ESMF_SUCCESS
  call NUOPC_CompSpecialize(comp, SetServices, &
       specLabel='Initialize', userRoutine=OCN_Initialize, rc=rc)
end subroutine OCN_SetServices

subroutine OCN_Initialize(comp, rc)
  use ESMF
  use NUOPC
  implicit none
  type(ESMF_GridComp), intent(inout) :: comp
  integer, intent(out) :: rc
  type(ESMF_Field) :: fld
  type(ESMF_Grid)  :: grid

  rc = ESMF_SUCCESS

  ! Create dummy 1x1 grid
  grid = ESMF_GridCreateNoPeriDim(ESMF_KIND_R8, dims=(/1,1/), rc=rc)

  ! Export temperature
  fld = ESMF_FieldCreate(grid, name='temperature', typekind=ESMF_TYPEKIND_R8, rc=rc)
  call NUOPC_Advertise(comp, export=fld, rc=rc)

  ! Import heat flux
  fld = ESMF_FieldCreate(grid, name='heat_flux', typekind=ESMF_TYPEKIND_R8, rc=rc)
  call NUOPC_Advertise(comp, import=fld, rc=rc)
end subroutine OCN_Initialize

