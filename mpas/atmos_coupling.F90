! ###########################################################################################
!> \file atmos_coupling.F90
!> Procedures for coupling the MPAS dynamical core to the CCPP Physics.
!>
! ###########################################################################################
module atmos_coupling_mod
  use mpas_kind_types, only : strKIND, RKIND
  use ufs_mpas_io,     only : domain_ptr,dyn_mpas_exchange_halo
  use ufs_mpas_io,     only : lucats,lumatch,luseas,lutype,li,albd,slmo,sfem,sfz0,therin,scfx,sfhc
  
  implicit none

  public :: ufs_physics_to_mpas
  public :: ufs_mpas_to_physics
  public :: ufs_microphysics_to_mpas
  public :: ufs_mpas_to_microphysics
  public :: ufs_mpas_grid_to_physics
  public :: ufs_mpas_sfc_to_physics
  public :: ufs_mpas_landuse_update
  public :: ufs_mpas_gwd_to_physics
  public :: ufs_mpas_reference_pressure
  
contains
  !> #########################################################################################
  !> Procedure to populate CCPP data containers with MPAS pool data.
  !> Called BEFORE CCPP Radiation and Physics Groups.
  !>
  !> Analogous to MPAS_to_physics in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> This procedure accesses MPAS data using MPAS native procedures and stores the data
  !> locally in the data-containers defined above. The MPAS "state" is then translated to the
  !> CCPP "state" needed by the physics.
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_physics(physics_state, surface_state, radiation)
    use GFS_typedefs,         only : GFS_statein_type, GFS_sfcprop_type, GFS_radtend_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use atm_core,             only : atm_compute_output_diagnostics
    use mpas_constants,       only : gravity, rvord

    ! Arguments
    type(GFS_statein_type),   intent(inout) :: physics_state
    type(GFS_sfcprop_type),   intent(inout) :: surface_state
    type(GFS_radtend_type),   intent(inout) :: radiation

    ! Locals
    type(mpas_pool_type), pointer :: state_pool, diag_pool, mesh_pool, diag_phys
    integer :: iCol, iLay, iTracer, ithread
    integer, pointer :: nCellsSolve, num_scalars, nVertLevels, index_qv
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND) :: rho1, rho2, tem1, tem2, theta, fzm_p, fzp_p, z0, z1, z2, w1, w2,rho_a
    real(kind=RKIND), pointer :: qv(:,:), qc(:,:), qr(:,:), qi(:,:), qs(:,:), qg(:,:)
    real(kind=RKIND), pointer :: ux(:,:), uy(:,:), theta_m(:,:), rho_zz(:,:), zgrid(:,:), zz(:,:)
    real(kind=RKIND), pointer :: exner(:,:), tracers(:,:,:), pressure_b(:,:), pressure_p(:,:)
    real(kind=RKIND), pointer :: w(:,:), surface_pressure(:),  prsi(:,:), prsl(:,:), rho(:,:)
    real(RKIND), pointer :: sfc_albedo(:),sfc_emiss(:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics',  diag_phys)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_pool,  'uReconstructZonal',      ux)
    call mpas_pool_get_array(diag_pool,  'uReconstructMeridional', uy)
    call mpas_pool_get_array(state_pool, 'scalars',                tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',                      w, timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'exner',                  exner)
    call mpas_pool_get_array(mesh_pool,  'zgrid',                  zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',                     zz)
    call mpas_pool_get_array(state_pool, 'theta_m',                theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',                 rho_zz,  timeLevel=1)
    call mpas_pool_get_array(diag_pool,  'pressure_base',          pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',             pressure_p)
    call mpas_pool_get_array(diag_pool,  'surface_pressure',       surface_pressure)
    call mpas_pool_get_array(diag_phys,  'sfc_albedo',             sfc_albedo)
    call mpas_pool_get_array(diag_phys,  'sfc_emiss' ,             sfc_emiss)
 
    ! Local variables
    allocate(prsl(nCellsSolve, nVertLevels))
    allocate(prsi(nCellsSolve, nVertLevels + 1))
    allocate(rho( nCellsSolve, nVertLevels))
    
    ! Copy fields from MPAS data containers to physics data containers.
    ! [k, i] -> [i, k]
    ! Retain bottom-up convention
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             ! Scalars (tracer,layer,col) -> (col,layer,tracer)
             do iTracer = 1,num_scalars
                physics_state % qgrs(iCol,iLay,iTracer) = max(0._RKIND, tracers(iTracer,iLay,iCol))
             end do

             ! Air denisty (rho) (TODO: Pass to CCPP Physics)
             rho(iCol,iLay) = zz(iLay,iCol) * rho_zz(iLay,iCol)

             ! Potential temperature (theta_m -> theta)
             theta = theta_m(iLay,iCol) / (1._RKIND + rvord * tracers(index_qv,iLay,iCol))

             ! Air temperature (theta -> temp)
             physics_state % tgrs(iCol,iLay)   = theta*exner(iLay,iCol)

             ! Winds at grid center
             physics_state % ugrs(iCol,iLay)   = ux(iLay,iCol)
             physics_state % vgrs(iCol,iLay)   = uy(iLay,iCol)

             ! Layer geopotential and height.
             physics_state % phil(iCol,iLay)   = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay,iCol))*gravity !(m -> m2/s2)
             physics_state % zgrid(iCol,iLay)  = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay,iCol))

             ! Level geopotential and height.
             physics_state % phii(iCol,iLay)   = zgrid(iLay,iCol)*gravity !(m -> m2/s2)
             physics_state % zigrid(iCol,iLay) = zgrid(iLay,iCol)

             ! Layer thickness.
             physics_state % dzgrid(iCol,iLay) = zgrid(iLay+1,iCol) - zgrid(iLay,iCol)

             ! Exner funciton
             physics_state % prslk(iCol,iLay)  = exner(iLay,iCol)

             ! MPAS provides vertical velocity at interfaces, compute layer mean, and convert from w -> omega
             physics_state % vvl(iCol,iLay) =  -0.5*(w(iLay,iCol) + w(iLay+1,iCol))*rho(iCol,iLay)*gravity

             ! Pressure (non-hydrostatic)
             prsl(iCol,iLay) = pressure_p(iLay,iCol) + pressure_b(iLay,iCol)
          end do
          do iLay = nVertLevels,nVertLevels+1
             physics_state % phii(iCol,iLay)     = zgrid(iLay,iCol)*gravity !(m -> m2/s2)
             physics_state % zigrid(iCol,iLay)   = zgrid(iLay,iCol)
             physics_state % dzgrid(iCol,iLay-1) = zgrid(iLay,iCol) - zgrid(iLay-1,iCol)
          end do

          ! Set surface temperature to lowest level temperature (revisit for coupling)
          theta = theta_m(1,iCol) / (1._RKIND + rvord * tracers(index_qv,1,iCol))
          surface_state % tsfc(iCol) = theta*exner(1,iCol)
          !
          ! Surface radiative properties
          ! DJS2026:
          ! For UFS-FV3, these fields come from set_emis() and set_alb() in
          ! GFS_radiation_surface.F90. The surface albedo also has spectral dependencies (nIR/uvvus)
          !
          ! For UFS-MPAS, these fields come from the MPAS sfc_input pool, which can be updated (daily)
          ! by calling ufs_mpas_landuse_update. The MPAS surface albedo is broadband,
          ! so we need to asign the same albedo for all channels in GFS_radiation_surface.
          radiation % semis(iCol) = sfc_emiss(iCol)
          radiation % salb(iCol)  = sfc_albedo(iCol)
       end do
    end do

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = rho_zz(1,iCol) * zz(1,iCol) * (1. + tracers(index_qv,1,iCol))
          rho2 = rho_zz(2,iCol) * zz(2,iCol) * (1. + tracers(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*(zgrid(2,iCol) - zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + pressure_b(1,iCol)
       end do
    end do

    ! Interpolation of pressure and temperature from layer-center to layer-interface
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 2,nVertLevels
             tem1 = 1./(zgrid(iLay+1,iCol)-zgrid(iLay-1,iCol))
             fzm_p = (zgrid(iLay,  iCol)-zgrid(iLay-1,iCol)) * tem1
             fzp_p = (zgrid(iLay+1,iCol)-zgrid(iLay,  iCol)) * tem1
             physics_state % tgri(iCol,iLay) = fzm_p*physics_state % tgrs(iCol,iLay) + fzp_p*physics_state % tgrs(iCol,iLay-1)
             prsi(iCol,iLay) = fzm_p*prsl(iCol,iLay) + fzp_p*prsl(iCol,iLay-1)
          enddo
       enddo
    enddo

    ! Interpolation of pressure and temperature to the top-of-the-model
    iLay = nVertLevels + 1
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          z0 = zgrid(iLay,iCol)
          z1 = 0.5*(zgrid(iLay  ,iCol)+zgrid(iLay-1,iCol))
          z2 = 0.5*(zgrid(iLay-1,iCol)+zgrid(iLay-2,iCol))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.-w1
          physics_state % tgri(iCol,iLay) = w1*physics_state % tgrs(iCol,iLay-1) + w2*physics_state % tgrs(iCol,iLay-2)
          prsi(iCol,iLay) = exp(w1*log(prsl(iCol,iLay-1))+w2*log(prsl(iCol,iLay-2)))
       end do
    end do

    ! Recalculate the pressure and temperature  at the surface as an extrapolation of
    ! the pressures in the 2 layers above the surface.
    iLay = 1
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          z0 = zgrid(iLay,iCol)
          z1 = 0.5*(zgrid(iLay  ,iCol)+zgrid(iLay+1,iCol))
          z2 = 0.5*(zgrid(iLay+1,iCol)+zgrid(iLay+2,iCol))
          w1 = (z0-z2)/(z1-z2)
          w2 = 1.-w1
          physics_state % tgri(iCol,iLay) = w1*physics_state % tgrs(iCol,iLay) + w2*physics_state % tgrs(iCol,iLay+1)
          prsi(iCol,iLay) = w1*prsl(iCol,iLay)+w2*prsl(iCol,iLay+1)
       end do
    end do

    ! Calculation of the hydrostatic pressure
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! Pressure at layer-interfaces
          iLay = nVertLevels + 1
          physics_state % prsi(iCol,iLay) = prsi(iCol,iLay)
          do iLay = nVertLevels,1,-1
             rho_a = rho(iCol,iLay) / (1.+tracers(index_qv,iLay,iCol))
             physics_state % prsi(iCol,iLay)  = physics_state % prsi(iCol,iLay+1) + &
                  gravity*rho(iCol,iLay)*physics_state % dzgrid(iCol,iLay)
          end do
          ! Pressure at layer-centers
          do iLay = nVertLevels,1,-1
             physics_state % prsl(iCol,iLay) = 0.5*(physics_state % prsi(iCol,iLay+1) + physics_state % prsi(iCol,iLay) )
          end do
          ! Pressure layer thickness
          do iLay = 1,nVertLevels
             physics_state % dp(iCol,iLay) = physics_state % prsi(iCol,iLay) - physics_state % prsi(iCol,iLay+1)
          end do
          ! Pressure difference across layer-centers
          physics_state % dp(iCol,1) = physics_state % prsi(iCol,1) - physics_state % prsi(iCol,2)
          do iLay = 2,nVertLevels
             physics_state % dp(iCol,iLay) = physics_state % prsi(iCol,iLay) - physics_state % prsi(iCol,iLay+1)
             physics_state % dpc(iCol,iLay)  = physics_state % prsl(iCol,iLay-1) - physics_state % prsl(iCol,iLay)
          end do
          physics_state % dpc(iCol,1)  = physics_state % prsi(iCol,1) - physics_state % prsl(iCol,1)
          ! Surface pressure
          physics_state % pgr(iCol) = physics_state % prsi(iCol,1)
       end do
    end do

    ! Housekeeping
    deallocate (prsl)
    deallocate (prsi)
    deallocate (rho)
    nullify (mesh_pool)
    nullify (state_pool)
    nullify (diag_pool)

  end subroutine ufs_mpas_to_physics

  !> #########################################################################################
  !> Procedure to update MPAS prognostic tendencies with physics (CCPP) tendencies.
  !> Called AFTER physics, BEFORE calling dynamics (current timestep).
  !>
  !> - Update scalar tendencies "scalars_tend" from MPAS "tend" pool
  !> - Convert from potential temerature (theta) to modified potential temperature (theta_m).
  !> - Update "theta_m" tendency from MPAS "tend" pool
  !> - Compute winds on grid edges
  !> - Update "ru" tendency from MPAS "tend" pool
  !>
  !> Analogous to phys_get_tend in physics/mpas_atmphys_todynamics.F
  !> Here, instead of updating the state with physics tendencies from the MPAS "tend_phys"
  !> pool, we use tendencies from the CCPP Physics data containers.
  !>
  !> #########################################################################################
  subroutine ufs_physics_to_mpas(physics_state)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_constants,     only : rv, rgas, gravity

    ! Arguments
    type(GFS_stateout_type), intent(in) :: physics_state

    ! Locals
    type(mpas_pool_type),  pointer :: state_pool, mesh_pool, tend_pool, diag_pool, tend_phys
    real(kind=RKIND), pointer :: mass(:,:), mass_edge(:,:), exner(:,:), theta_m(:,:), zgrid(:,:), zz(:,:)
    real(kind=RKIND), pointer :: pressure_b(:,:), pressure_p(:,:), tend_th_phys(:,:)
    real(kind=RKIND), pointer :: tend_theta_phys(:,:), tend_theta_dyn(:,:)
    real(kind=RKIND), pointer :: tend_u_phys(:,:), tend_ru_dyn(:,:)
    real(kind=RKIND), pointer :: tend_uzonal(:,:), tend_umerid(:,:)
    real(kind=RKIND), pointer :: scalars(:,:,:), tend_scalars_phys(:,:,:), tend_scalars_dyn(:,:,:)
    real(kind=RKIND), pointer :: surface_pressure(:)
    integer, pointer :: nCells, nCellsSolve, num_scalars, nVertLevels, nEdges, nEdgesSolve
    integer, pointer :: index_qv => null()
    integer, pointer :: index_qc => null()
    integer, pointer :: index_qi => null()
    integer, pointer :: index_qr => null()
    integer, pointer :: index_qs => null()
    integer, pointer :: index_qg => null()
    integer, pointer :: index_nc => null()
    integer, pointer :: index_ni => null()
    integer, pointer :: index_nifa => null()
    integer, pointer :: index_nwfa => null()
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer :: iCol,iLay,ithread,iScalar
    real(kind=RKIND):: coeff, tem1, tem2, rho1, rho2
    logical :: debug=.false.
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_physics_to_mpas'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nEdges',               nEdges)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nEdgesSolve',          nEdgesSolve)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state',        state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',         mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend',         tend_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',         diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend_physics', tend_phys)

    ! Get MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(mesh_pool,  'nCells',      nCells)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'index_qc',    index_qc)
    call mpas_pool_get_dimension(state_pool, 'index_qi',    index_qi)
    call mpas_pool_get_dimension(state_pool, 'index_qr',    index_qr)
    call mpas_pool_get_dimension(state_pool, 'index_qs',    index_qs)
    call mpas_pool_get_dimension(state_pool, 'index_qg',    index_qg)
    call mpas_pool_get_dimension(state_pool, 'index_nc',    index_nc)
    call mpas_pool_get_dimension(state_pool, 'index_ni',    index_ni)
    call mpas_pool_get_dimension(state_pool, 'index_nifa',  index_nifa)
    call mpas_pool_get_dimension(state_pool, 'index_nwfa',  index_nwfa)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool,'theta_m',          theta_m,1)
    call mpas_pool_get_array(state_pool,'scalars',          scalars,1)
    call mpas_pool_get_array(state_pool,'rho_zz',           mass,   1)
    call mpas_pool_get_array(mesh_pool, 'zz',               zz)
    call mpas_pool_get_array(mesh_pool, 'zgrid',            zgrid)
    call mpas_pool_get_array(diag_pool, 'rho_edge',         mass_edge)
    call mpas_pool_get_array(diag_pool, 'surface_pressure', surface_pressure)
    call mpas_pool_get_array(diag_pool, 'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool, 'pressure_p',       pressure_p)
    call mpas_pool_get_array(diag_pool, 'exner',            exner)
    call mpas_pool_get_array(tend_phys, 'tend_uzonal',      tend_uzonal)
    call mpas_pool_get_array(tend_phys, 'tend_umerid',      tend_umerid)

    ! Allocate/initialize local variables
    allocate(tend_th_phys(nVertLevels,nCellsSolve+1))
    allocate(tend_theta_phys(nVertLevels,nCellsSolve+1))
    allocate(tend_scalars_phys(num_scalars, nVertLevels,nCellsSolve+1))
    allocate(tend_u_phys(nVertLevels,nEdges+1))
    tend_th_phys(:,:)        = 0._RKIND
    tend_theta_phys(:,:)     = 0._RKIND
    tend_scalars_phys(:,:,:) = 0._RKIND
    tend_u_phys(:,:)         = 0._RKIND

    !> GJF: Add accumulated tendencies from the physics group

    !> #####################################################################################
    !> 1) Update MPAS "scalar" tendencies.
    !> #####################################################################################

    ! Specific humidity
    do ithread=1,nThreads
      do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        do iLay = 1,nVertLevels 
           tend_scalars_phys(index_qv,iLay,iCol) = tend_scalars_phys(index_qv,iLay,iCol) + &
                physics_state % dqdt(iCol,iLay,index_qv)*mass(iLay,iCol)
        end do
      end do
    end do

    ! Liquid cloud water
    if(associated(index_qc)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_qc,iLay,iCol) = tend_scalars_phys(index_qc,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_qc)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Ice cloud water
    if(associated(index_qi)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_qi,iLay,iCol) = tend_scalars_phys(index_qi,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_qi)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Rain water
    if(associated(index_qr)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_qr,iLay,iCol) = tend_scalars_phys(index_qr,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_qr)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Snow
    if(associated(index_qs)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_qs,iLay,iCol) = tend_scalars_phys(index_qs,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_qs)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Graupel
    if(associated(index_qg)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_qg,iLay,iCol) = tend_scalars_phys(index_qg,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_qg)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Liquid number concentration
    if(associated(index_nc)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_nc,iLay,iCol) = tend_scalars_phys(index_nc,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_nc)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Ice number concentration
    if(associated(index_ni)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_ni,iLay,iCol) = tend_scalars_phys(index_ni,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_ni)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! "Ice friendly" aerosol number concentration
    if(associated(index_nifa)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_nifa,iLay,iCol) = tend_scalars_phys(index_nifa,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_nifa)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! "Water friendly" aerosol number concentration
    if(associated(index_nwfa)) then
      do ithread=1,nThreads
        do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels 
             tend_scalars_phys(index_nwfa,iLay,iCol) = tend_scalars_phys(index_nwfa,iLay,iCol) + &
                  physics_state % dqdt(iCol,iLay,index_nwfa)*mass(iLay,iCol)
          end do
        end do
      end do
    end if

    ! Update MPAS tendencies (scalars)
    call mpas_pool_get_array(tend_pool, 'scalars_tend', tend_scalars_dyn)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
             do iScalar = 1, num_scalars
                tend_scalars_dyn(iScalar,iLay,iCol) = tend_scalars_dyn(iScalar,iLay,iCol) + tend_scalars_phys(iScalar,iLay,iCol)
             end do
          end do
       end do
    end do

    !> #####################################################################################
    !> 2) Update MPAS tendency "theta_m"
    !> #####################################################################################
    ! Add accumulated temperature tendencies from the physics group, convert from temperature
    ! to rho*potential temperature (T -> rtheta).
    do ithread=1,nThreads
       do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             tend_th_phys(iLay,iCol) = tend_th_phys(iLay,iCol) + (physics_state % dtdt(iCol,iLay)/exner(iLay,iCol))*mass(iLay,iCol)
          end do
       end do
    end do

    ! Convert from the tendency of potential temperature to the tendency of the  modified
    ! potential temperature (rtheta -> rtheta_m).
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
             coeff = (1. + rv/rgas * scalars(index_qv,iLay,iCol))
             tend_th_phys(iLay,iCol) = coeff * tend_th_phys(iLay,iCol) + rv/rgas * theta_m(iLay,iCol) * tend_scalars_phys(index_qv,iLay,iCol) / coeff
             tend_theta_phys(iLay,iCol) = tend_theta_phys(iLay,iCol) + tend_th_phys(iLay,iCol)
          end do
       end do
    end do

    ! Update MPAS tendency (theta_m)
    call mpas_pool_get_array(tend_pool, 'theta_m', tend_theta_dyn)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1, nVertLevels
               tend_theta_dyn(iLay,iCol) = tend_theta_dyn(iLay,iCol) + tend_theta_phys(iLay,iCol)
          end do
       end do
    end do

    !> #####################################################################################
    !> 3) Update MPAS tendency "ru"
    !> #####################################################################################

    ! First, need to regrid from grid-centers (physics) to grid-edges (dynamics)...
    ! Copy CCPP accumulated physics tendencies to MPAS scratch arrays.
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
           do iLay = 1, nVertLevels
              tend_uzonal(iLay,iCol) = physics_state % dudt(iCol,iLay)
              tend_umerid(iLay,iCol) = physics_state % dvdt(iCol,iLay)
           end do
       end do
    end do

    ! Next, update the halo points of the scratch arrays.
    call dyn_mpas_exchange_halo('tend_uzonal',debug)
    call dyn_mpas_exchange_halo('tend_umerid',debug)

    ! Finally, compute wind tendency at grid-edges.
    call tend_toEdges(mesh_pool, tend_uzonal, tend_umerid, tend_u_phys)

    ! Update MPAS tendency (ru)
    call mpas_pool_get_array(tend_pool, 'u', tend_ru_dyn)
    do iCol = 1,nEdgesSolve
       do iLay = 1, nVertLevels
          tend_ru_dyn(iLay,iCol) = tend_ru_dyn(iLay,iCol) + tend_u_phys(iLay,iCol)*mass_edge(iLay,iCol)
       end do
    end do

    !> #####################################################################################
    !> Diagnostics
    !> #####################################################################################
    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = mass(1,iCol) * zz(1,iCol) * (1. + scalars(index_qv,1,iCol))
          rho2 = mass(2,iCol) * zz(2,iCol) * (1. + scalars(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*( zgrid(2,iCol) -  zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + &
               pressure_b(1,iCol)
       enddo
    enddo

    ! Housekeeping
    deallocate(tend_th_phys)
    deallocate(tend_theta_phys)
    deallocate(tend_scalars_phys)
    deallocate(tend_u_phys)
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)
    nullify (tend_phys)

  end subroutine ufs_physics_to_mpas

  !> #########################################################################################
  !> Procedure to compute diabatic heating tendency from microphysics and store in MPAS pool.
  !> Called AFTER microphysics, BEFORE calling dynamics (next timestep)
  !>
  !> This routine is where we compute the diabatic heating rate, "rt_diabatic_tend", due to the
  !> microphysics.
  !> Follows microphysics_to_MPAS from src/core_atmosphere/physics/mpas_atmphys_interface.F.
  !> compute the diabatic heating rate, "rt_diabatic_tend" and save in MPAS memory for use by
  !> the dynamics at the subsequent time step in dynamics/mpas_atm_time_integration.F.
  !> Additionally, update any other fields needed by the dynamics (e.g., theta_m, rtheta_p)
  !>
  !> #########################################################################################
  subroutine ufs_microphysics_to_mpas(physics_state)
    use GFS_typedefs,       only : GFS_stateout_type
    use mpas_derived_types, only : mpas_pool_type
    use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_array
    use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_config
    use mpas_constants,     only : gravity, rvord, rv, rgas, p0, cp

    ! Arguments
    type(GFS_stateout_type),     intent(in   ) :: physics_state

    ! Locals
    type(mpas_pool_type), pointer :: diag_pool, mesh_pool, state_pool, tend_pool
    integer :: iCol, ithread, iLay, iTracer
    integer, pointer :: nCellsSolve
    integer, pointer :: index_qv => null()
    integer, pointer :: num_scalars, nVertLevels
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND) :: rho1, rho2, tem1, tem2, coeff, rcv, theta_dyn
    real(kind=RKIND), pointer :: config_dt
    real(kind=RKIND), pointer :: tracers(:,:,:), rt_diabatic_tend(:,:), rho_zz(:,:), theta_m(:,:)
    real(kind=RKIND), pointer :: zz(:,:), zgrid(:,:), exner(:,:), exner_b(:,:), rtheta_b(:,:), theta(:,:)
    real(kind=RKIND), pointer :: rtheta_p(:,:), pressure_b(:,:), pressure_p(:,:), surface_pressure(:), dtheta_dt_mp(:,:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_microphysics_to_mpas'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'tend',  tend_pool)

    ! Model timestep
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_dt', config_dt)

    ! MPAS dimensions
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(state_pool, 'scalars',          tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'rho_zz',           rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m',          theta_m, timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zgrid',            zgrid)
    call mpas_pool_get_array(mesh_pool,  'zz',               zz)
    call mpas_pool_get_array(diag_pool,  'pressure_base',    pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p',       pressure_p)
    call mpas_pool_get_array(diag_pool,  'surface_pressure', surface_pressure)
    call mpas_pool_get_array(diag_pool,  'exner',            exner)
    call mpas_pool_get_array(diag_pool,  'exner_base',       exner_b)
    call mpas_pool_get_array(diag_pool,  'rtheta_p',         rtheta_p)
    call mpas_pool_get_array(diag_pool,  'rtheta_base',      rtheta_b)
    call mpas_pool_get_array(tend_pool,  'rt_diabatic_tend', rt_diabatic_tend)
    call mpas_pool_get_array(diag_pool,  'dtheta_dt_mp',     dtheta_dt_mp)

    ! The MPAS version of microphysics schemes update the state within the dynamics;
    ! for CCPP/UFS, we will need to update the state variables here for use by the dynamics.
    ! Also, Update water vapor, cloud liquid water, rain mixing ratios, modified potential temperature,
    ! and potential temperature heating rate from microphysics
    allocate(theta(nVertLevels, nCellsSolve))
    rcv = rgas/(cp-rgas)
    do ithread=1,nThreads
      do iCol=cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
        do iLay = 1,nVertLevels
          
          ! Initialize diabatic heating tendency (initially fill with theta_m before updating)
          rt_diabatic_tend(iLay,iCol) = theta_m(iLay,iCol)
          
          ! Update potential temperature (theta) with microphysics tendency
          coeff = (1._RKIND + rvord * tracers(index_qv,iLay,iCol))
          theta_dyn = theta_m(ilay,iCol)/coeff
          theta(iLay,iCol) = theta_dyn + config_dt * (physics_state % dtdt(iCol,iLay) / exner(iLay,iCol))
          
          ! Scalars (col,layer,tracer) -> (tracer,layer,col)
          do iTracer = 1,num_scalars
            tracers(iTracer,iLay,iCol) = max(0._RKIND, tracers(iTracer,iLay,iCol) + config_dt * physics_state % dqdt(iCol,iLay,iTracer))
          end do
          
          ! update the virtual temperature coefficient with updated qv
          coeff = (1._RKIND + rvord * tracers(index_qv,iLay,iCol))
          ! Modified potential temperature (theta ->theta_m)
          theta_m(iLay,iCol) = theta(iLay,iCol)*coeff
          
          ! Now compute diabatic heating due to microphsyics, save for next time step
          rt_diabatic_tend(iLay,iCol) = (theta_m(iLay,iCol) - rt_diabatic_tend(iLay,iCol)) / config_dt
          
          ! Save the straight theta tendency due to microphysics
          dtheta_dt_mp(iLay,iCol) =  (theta(iLay,iCol) - theta_dyn) / config_dt
          
          ! Density weighted perturbation potential temperature
          rtheta_p(iLay,iCol) = rho_zz(iLay,iCol) * theta_m(iLay,iCol) - rtheta_b(iLay,iCol)
          
          ! Exner function
          exner(iLay,iCol) = (zz(iLay,iCol)*(rgas/P0)*(rtheta_p(iLay,iCol)+rtheta_b(iLay,iCol)))**rcv

          ! Perturbation pressure
          pressure_p(iLay,iCol) = zz(iLay,iCol)*rgas*(exner(iLay,iCol)*rtheta_p(iLay,iCol) + &
                                    (exner(iLay,iCol)-exner_b(iLay,iCol))*rtheta_b(iLay,iCol))  
        end do
      end do
    end do
    
    ! write(*,*) 'num_scalars',num_scalars
    ! if (associated(index_qv)) write(*,*) 'mean max/min ten qv',sum(tracers(index_qv,:,:)) / real(size(tracers(index_qv,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qv)), config_dt*minval(physics_state % ten_q(:,:,index_qv))
    ! if (associated(index_qc)) write(*,*) 'mean max/min ten qc',sum(tracers(index_qc,:,:)) / real(size(tracers(index_qc,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qc)), config_dt*minval(physics_state % ten_q(:,:,index_qc)) 
    ! if (associated(index_qi)) write(*,*) 'mean max/min ten qi',sum(tracers(index_qi,:,:)) / real(size(tracers(index_qi,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qi)), config_dt*minval(physics_state % ten_q(:,:,index_qi))
    ! if (associated(index_qr)) write(*,*) 'mean max/min ten qr',sum(tracers(index_qr,:,:)) / real(size(tracers(index_qr,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qr)), config_dt* minval(physics_state % ten_q(:,:,index_qr))
    ! if (associated(index_qs)) write(*,*) 'mean max/min ten qs',sum(tracers(index_qs,:,:)) / real(size(tracers(index_qs,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qs)), config_dt*minval(physics_state % ten_q(:,:,index_qs))
    ! if (associated(index_qg)) write(*,*) 'mean max/min ten qg',sum(tracers(index_qg,:,:)) / real(size(tracers(index_qg,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_qg)), config_dt*minval(physics_state % ten_q(:,:,index_qg))
    ! if (associated(index_nc)) write(*,*) 'mean max/min ten nc',sum(tracers(index_nc,:,:)) / real(size(tracers(index_nc,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_nc)), config_dt*minval(physics_state % ten_q(:,:,index_nc))
    ! if (associated(index_ni)) write(*,*) 'mean max/min ten ni',sum(tracers(index_ni,:,:)) / real(size(tracers(index_ni,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_ni)), config_dt*minval(physics_state % ten_q(:,:,index_ni))
    ! if (associated(index_nr)) write(*,*) 'mean max/min ten nr',sum(tracers(index_nr,:,:)) / real(size(tracers(index_nr,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_nr)), config_dt*minval(physics_state % ten_q(:,:,index_nr))
    ! if (associated(index_ns)) write(*,*) 'mean max/min ten ns',sum(tracers(index_ns,:,:)) / real(size(tracers(index_ns,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_ns)), config_dt*minval(physics_state % ten_q(:,:,index_ns))
    ! if (associated(index_ng)) write(*,*) 'mean max/min ten ng',sum(tracers(index_ng,:,:)) / real(size(tracers(index_ng,:,:))),config_dt*maxval(physics_state %ten_q(:,:,index_ng)), config_dt*minval(physics_state % ten_q(:,:,index_ng))
    ! if (associated(index_nifa)) write(*,*) 'mean max/min ten nifa',sum(tracers(index_nifa,:,:)) / real(size(tracers(index_nifa,:,:))), config_dt*maxval(physics_state % ten_q(:,:,index_nifa)), config_dt*minval(physics_state % ten_q(:,:,index_nifa))
    ! if (associated(index_nwfa)) write(*,*) 'mean max/min ten nwfa',sum(tracers(index_nwfa,:,:)) /real(size(tracers(index_nwfa,:,:))),config_dt*maxval(physics_state % ten_q(:,:,index_nwfa)), config_dt*minval(physics_state % ten_q(:,:,index_nwfa))
    

    ! Calculation of the surface pressure using hydrostatic assumption down to the surface.
    ! (from mpas_atmphys_interface.F:MPAS_to_physics())
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          tem1 = zgrid(2,iCol) - zgrid(1,iCol)
          tem2 = zgrid(3,iCol) - zgrid(2,iCol)
          rho1 = rho_zz(1,iCol) * zz(1,iCol) * (1. + tracers(index_qv,1,iCol))
          rho2 = rho_zz(2,iCol) * zz(2,iCol) * (1. + tracers(index_qv,2,iCol))
          surface_pressure(iCol) = 0.5*gravity*(zgrid(2,iCol) - zgrid(1,iCol)) &
               * (rho1 - 0.5*(rho2-rho1)*tem1/(tem1+tem2))
          surface_pressure(iCol) = surface_pressure(iCol) + pressure_p(1,iCol) + pressure_b(1,iCol)
       end do
    end do

    ! Housekeeping
    nullify (state_pool)
    nullify (mesh_pool)
    nullify (diag_pool)

  end subroutine ufs_microphysics_to_mpas

  !> #########################################################################################
  !> Procedure to update physics (CCPP) state using updated MPAS state.
  !> Called AFTER dynamics, BEFORE microphysics.
  !> 
  !> Analogous to microphysics_from_MPAS in src/core_atmosphere/physics/mpas_atmphys_interface.F
  !>
  !> #########################################################################################
  subroutine ufs_mpas_to_microphysics(physics_state, physics_statein)
    use GFS_typedefs,         only : GFS_stateout_type, GFS_statein_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_array, mpas_pool_get_dimension
    use mpas_constants,       only : rvord, gravity

    ! Arguments
    type(GFS_stateout_type), intent(inout) :: physics_state
    type(GFS_statein_type),  intent(inout) :: physics_statein

    ! Locals
    type(mpas_pool_type), pointer :: state_pool, diag_pool, mesh_pool
    integer :: ithread, iCol, iTracer, iLay
    integer, pointer :: num_scalars, nVertLevels, nCellsSolve, index_qv
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(kind=RKIND), pointer :: rho_zz(:,:), theta_m(:,:), zz(:,:), zgrid(:,:), exner(:,:)
    real(kind=RKIND), pointer :: tracers(:,:,:), rho(:,:), pressure_b(:,:), pressure_p(:,:), w(:,:)
    real(kind=RKIND) :: theta, pres
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_to_microphysics'

    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'state', state_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag',  diag_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)

    call mpas_pool_get_dimension(mesh_pool,  'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(mesh_pool,  'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(state_pool, 'num_scalars', num_scalars)
    call mpas_pool_get_dimension(state_pool, 'index_qv',    index_qv)

    call mpas_pool_get_array(state_pool, 'rho_zz',       rho_zz,  timeLevel=1)
    call mpas_pool_get_array(state_pool, 'theta_m',      theta_m, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'scalars',      tracers, timeLevel=1)
    call mpas_pool_get_array(state_pool, 'w',            w,       timeLevel=1)
    call mpas_pool_get_array(mesh_pool,  'zz',           zz)
    call mpas_pool_get_array(mesh_pool,  'zgrid',        zgrid)
    call mpas_pool_get_array(diag_pool,  'exner',        exner)
    call mpas_pool_get_array(diag_pool,  'pressure_base',pressure_b)
    call mpas_pool_get_array(diag_pool,  'pressure_p'   ,pressure_p)

    allocate(rho(nCellsSolve, nVertLevels))
    ! Update fields needed by microphysics...
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          do iLay = 1,nVertLevels
             ! Scalars (tracer,layer,col) -> (col,layer,tracer)
             do iTracer = 1,num_scalars
                physics_state % gq0(iCol,iLay,iTracer) = max(0._RKIND, tracers(iTracer,iLay,iCol))
             end do

             ! Air denisty (rho) (TODO: Pass to CCPP Physics)
             rho(iCol,iLay) = zz(iLay,iCol) * rho_zz(iLay,iCol)

             ! Potential temperature (theta_m -> theta)
             theta = theta_m(iLay,iCol) / (1._RKIND + rvord * max(0._RKIND,tracers(index_qv,iLay,iCol)))

             ! Air temperature (theta -> t)
             physics_state % gt0(iCol,iLay) = theta*exner(iLay,iCol)

             ! Vertical velocity (w -> omega)
             physics_statein % vvl(iCol,iLay) = -0.5*(w(iLay,iCol) + w(iLay+1,iCol))*rho(iCol,iLay)*gravity
          end do
       end do
    end do
    deallocate(rho)
    nullify(diag_pool)
    nullify(mesh_pool)
    nullify(state_pool)

    ! DJS:  Update hydrostatic pressure after dynamics, before MP?
    
    ! GJF: Remove microphysics heating from state before calling microphysics. This is done
    ! at line 3317 of mpas_atm_time_integration.F/atm_recover_large_step_variables_work.
 
  end subroutine ufs_mpas_to_microphysics

!> #########################################################################################
!> Procedure to transfer MPAS grid information to physics DDTs.
!>
!> #########################################################################################
  subroutine ufs_mpas_grid_to_physics(physics_grid)
    use GFS_typedefs,         only : GFS_grid_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config
    use mpas_constants,       only : pii
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    ! Arguments
    type(GFS_grid_type),      intent(inout) :: physics_grid

    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool
    integer :: i, ierr, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(RKIND), pointer :: lat(:), lon(:), area(:), meshDensity(:)
    real(RKIND), pointer :: nominalMinDc, config_len_disp
    real(RKIND)          :: rad2deg
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_grid_to_physics'

    ierr = 0
    rad2deg = 180.0_RKIND/pii
    
    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    
    call mpas_pool_get_array(mesh_pool,  'latCell',     lat)
    call mpas_pool_get_array(mesh_pool,  'lonCell',     lon)
    call mpas_pool_get_array(mesh_pool,  'areaCell',    area)
    call mpas_pool_get_array(mesh_pool,  'meshDensity', meshDensity)
    
    ! (from mpas_atm_core.F/atm_core_init Determine horizontal length scale used by horizontal diffusion and 3-d divergence damping
    nullify(nominalMinDc)
    call mpas_pool_get_array(mesh_pool, 'nominalMinDc', nominalMinDc)

    nullify(config_len_disp)
    call mpas_pool_get_config(domain_ptr % blocklist % configs, 'config_len_disp', config_len_disp)

    ! If config_len_disp was specified as a valid value, use that
    if (config_len_disp > 0.0_RKIND) then
      ! But if nominalMinDc was available in the input file and is different, print a warning
      if (nominalMinDc > 0.0_RKIND .and. abs(nominalMinDc - config_len_disp) > 1.0e-6_RKIND * config_len_disp) then
        call mpas_log_write(subname // ' WARNING: nominalMinDc was read from input file as a positive value ($r) that differs', &
                            realArgs=[nominalMinDc], messageType=MPAS_LOG_WARN)
        call mpas_log_write(subname // ' WARNING: from the specified config_len_disp value ($r)', &
                            realArgs=[config_len_disp], messageType=MPAS_LOG_WARN)
      end if
      nominalMinDc = config_len_disp
    ! Otherwise, try to use nominalMinDc
    else
      if (nominalMinDc > 0.0_RKIND) then
        call mpas_log_write('Setting config_len_disp to $r based on nominalMinDc value in input file', realArgs=[nominalMinDc])
          config_len_disp = nominalMinDc
      else
         call mpas_log_write(subname // ' ERROR: Both config_len_disp and nominalMinDc are <= 0.0.', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write(subname // ' ERROR: Please either specify config_len_disp in the &nhyd_model namelist group,', &
                             messageType=MPAS_LOG_ERR)
         call mpas_log_write(subname // ' ERROR: or use an input file that provides a valid value for the nominalMinDc variable.', &
                             messageType=MPAS_LOG_ERR)
        ierr = 1
      end if
    end if
    if (ierr/=0)  call mpas_log_write(subname // ' ERROR: Call to ufs_mpas_grid_to_physics() failed', messageType=MPAS_LOG_CRIT)

    do ithread = 1,nThreads
       do i = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          physics_grid % xlat(i)   = lat(i)
          physics_grid % xlon(i)   = lon(i)
          physics_grid % xlat_d(i) = physics_grid % xlat(i) * rad2deg
          physics_grid % xlon_d(i) = physics_grid % xlon(i) * rad2deg
         ! write(*,*) 'ithread, i, xlat, xlon, xlat_d, xlon_d',ithread, i, physics_grid % xlat(i),physics_grid % xlon(i),physics_grid % xlat_d(i),physics_grid % xlon_d(i)
          physics_grid % sinlat(i) = sin(physics_grid % xlat(i))
          physics_grid % coslat(i) = sqrt(1.0_RKIND - physics_grid % sinlat(i) * physics_grid % sinlat(i))
          physics_grid % area(i)   = area(i)
          !formula for dx comes from mpas_atmphys_driver_gwdo.F instead of sqrt(area) as in FV3
          physics_grid % dx(i)     = config_len_disp / meshDensity(i)**0.25
       end do
    end do
    
  end subroutine ufs_mpas_grid_to_physics

!> #########################################################################################
!> Procedure to transfer MPAS information to physics srfprop DDT
!>
!> #########################################################################################
  subroutine ufs_mpas_sfc_to_physics(physics_sfcprop)
    use GFS_typedefs,         only : GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array

    ! Arguments
    type(GFS_sfcprop_type),      intent(inout) :: physics_sfcprop
    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool, sfc_input
    integer :: i, ierr, iCol, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:), landmask(:)
    integer, pointer :: ivgtyp(:)
    real(RKIND), pointer :: sst(:), snow(:), tmn(:), albbck(:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_sfc_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',  mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)

    !using fv3atm_sfc_io.F90/Sfc_io_transfer() as a template; mpas_init_atm_static.F from MPAS-model for syntax
    call mpas_pool_get_array(sfc_input, 'landmask',  landmask)
    call mpas_pool_get_array(sfc_input, 'ivgtyp',    ivgtyp)
    call mpas_pool_get_array(sfc_input, 'sst',       sst)
    call mpas_pool_get_array(sfc_input, 'snow',      snow)
    call mpas_pool_get_array(sfc_input, 'tmn' ,      tmn)
    call mpas_pool_get_array(sfc_input, 'sfc_albbck',albbck)
    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          physics_sfcprop % slmsk(iCol) = landmask(iCol)
          ! Vegetation type. Valid because the MPAS mesh uses
          ! mminlu = MODIFIED_IGBP_MODIS_NOAH and the namelist sets ivegsrc = 1 (IGBP),
          physics_sfcprop % vtype(iCol) = ivgtyp(iCol)
          physics_sfcprop % tsfco(iCol) = sst(iCol)
          physics_sfcprop % weasd(iCol) = snow(iCol)
          physics_sfcprop % tg3(iCol)   = tmn(iCol)
          !zorl - z0/znt in MPAS, read in as sfz0; landuse_init_forMPAS not used yet; diag_physics pool, z0 - not initialized yet
          !alvsf, alvwf, alnsf, alnwf - MPAS doesn't split into visible/nir and strong/weak coszen dependency; set these to the value that we have (background snow-free albedo of surface)?
          physics_sfcprop % alvsf(iCol) = albbck(iCol)
          physics_sfcprop % alvwf(iCol) = albbck(iCol)
          physics_sfcprop % alnsf(iCol) = albbck(iCol)
          physics_sfcprop % alnwf(iCol) = albbck(iCol)
       end do
    end do

  end subroutine ufs_mpas_sfc_to_physics

  !> #########################################################################################
  !> Procedure to populate CCPP data container with MPAS pool data.
  !> Initial condition fields needed by the CCPP GWD parameterization.
  !> 
  !> #########################################################################################
  subroutine ufs_mpas_gwd_to_physics(control,surface)
    use GFS_typedefs,         only : GFS_control_type
    use GFS_typedefs,         only : GFS_sfcprop_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
    ! Arguments
    type(GFS_control_type), intent(in   ) :: control
    type(GFS_sfcprop_type), intent(inout) :: surface
    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool, sfc_input
    integer :: i, ierr, iCol, ithread
    integer, pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    real(RKIND), pointer :: var2d(:), con(:), oa1(:), oa2(:), oa3(:), oa4(:)
    real(RKIND), pointer :: ol1(:), ol2(:), ol3(:), ol4(:)
    real(RKIND), pointer :: var2dss(:), conss(:), oa1ss(:), oa2ss(:), oa3ss(:), oa4ss(:)
    real(RKIND), pointer :: ol1ss(:), ol2ss(:), ol3ss(:), ol4ss(:)
    real(RKIND), pointer :: var2dls(:), conls(:), oa1ls(:), oa2ls(:), oa3ls(:), oa4ls(:)
    real(RKIND), pointer :: ol1ls(:), ol2ls(:), ol3ls(:), ol4ls(:)
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_gwd_to_physics'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',      mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input', sfc_input)

    ! DJS: Do we need these? Only used when coupling to YSU scheme?
    call mpas_pool_get_array(sfc_input, 'var2d', var2d)
    call mpas_pool_get_array(sfc_input, 'con',   con)
    call mpas_pool_get_array(sfc_input, 'oa1',   oa1)
    call mpas_pool_get_array(sfc_input,	'oa2',   oa2)
    call mpas_pool_get_array(sfc_input,	'oa3',   oa3)
    call mpas_pool_get_array(sfc_input,	'oa4',   oa4)
    call mpas_pool_get_array(sfc_input, 'ol1',   ol1)
    call mpas_pool_get_array(sfc_input, 'ol2',   ol2)
    call mpas_pool_get_array(sfc_input, 'ol3',   ol3)
    call mpas_pool_get_array(sfc_input, 'ol4',   ol4)
    !
    if (control%gwd_opt==3 .or. control%gwd_opt==33 .or. &
        control%gwd_opt==2 .or. control%gwd_opt==22 ) then
       call mpas_pool_get_array(sfc_input, 'var2dss', var2dss)
       call mpas_pool_get_array(sfc_input, 'conss',   conss)
       call mpas_pool_get_array(sfc_input, 'oa1ss',   oa1ss)
       call mpas_pool_get_array(sfc_input, 'oa2ss',   oa2ss)
       call mpas_pool_get_array(sfc_input, 'oa3ss',   oa3ss)
       call mpas_pool_get_array(sfc_input, 'oa4ss',   oa4ss)
       call mpas_pool_get_array(sfc_input, 'ol1ss',   ol1ss)
       call mpas_pool_get_array(sfc_input, 'ol2ss',   ol2ss)
       call mpas_pool_get_array(sfc_input, 'ol3ss',   ol3ss)
       call mpas_pool_get_array(sfc_input, 'ol4ss',   ol4ss)
       call mpas_pool_get_array(sfc_input, 'var2dls', var2dls)
       call mpas_pool_get_array(sfc_input, 'conls',   conls)
       call mpas_pool_get_array(sfc_input, 'oa1ls',   oa1ls)
       call mpas_pool_get_array(sfc_input, 'oa2ls',   oa2ls)
       call mpas_pool_get_array(sfc_input, 'oa3ls',   oa3ls)
       call mpas_pool_get_array(sfc_input, 'oa4ls',   oa4ls)
       call mpas_pool_get_array(sfc_input, 'ol1ls',   ol1ls)
       call mpas_pool_get_array(sfc_input, 'ol2ls',   ol2ls)
       call mpas_pool_get_array(sfc_input, 'ol3ls',   ol3ls)
       call mpas_pool_get_array(sfc_input, 'ol4ls',   ol4ls)
    end if

    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! DJS: Do we need these? Only used when coupling to YSU scheme?
          surface % hprime(iCol,1) = var2d(iCol)
          surface % oc(iCol)       = con(iCol)
          surface % oa4(iCol,1)    = oa1(iCol)
          surface % oa4(iCol,2)    = oa2(iCol)
          surface % oa4(iCol,3)    = oa3(iCol)
          surface % oa4(iCol,4)    = oa4(iCol)
          surface % clx(iCol,1)    = ol1(iCol)
          surface % clx(iCol,2)    = ol2(iCol)
          surface % clx(iCol,3)    = ol3(iCol)
          surface % clx(iCol,4)    = ol4(iCol)
          ! DJS: Where are these set? Do we need them?
          !      drag_suite_run doesn't use these...
          !      GFS_GWD_generic_pre sets them to zero, since mntvar(11:14) is initialized to zero.
          !      Not using GFS_GWD_generic_pre in UFS-MPAS, instead these are initializaed to zero elsewhere.

          ! UGWPv1 requires these fields, TODO: do the other options need them?
          if (control%gwd_opt==2) then
             surface % gamma(iCol) = 0.0
             surface % sigma(iCol) = 0.0
             surface % theta(iCol) = 0.0
             surface % elvmax(iCol) = 0.0
          end if
          if (control%gwd_opt==3 .or. control%gwd_opt==33 .or. &
              control%gwd_opt==2 .or. control%gwd_opt==22 ) then
             surface % hprime(iCol,1) = var2dls(iCol)
             surface % oc(iCol)       = conls(iCol)
             surface % oa4(iCol,1)    = oa1ls(iCol)
             surface % oa4(iCol,2)    = oa2ls(iCol)
             surface % oa4(iCol,3)    = oa3ls(iCol)
             surface % oa4(iCol,4)    = oa4ls(iCol)
             surface % clx(iCol,1)    = ol1ls(iCol)
             surface % clx(iCol,2)    = ol2ls(iCol)
             surface % clx(iCol,3)    = ol3ls(iCol)
             surface % clx(iCol,4)    = ol4ls(iCol)
             surface % varss(iCol)    = var2dss(iCol)
             surface % ocss(iCol)     = conss(iCol)
             surface % oa4ss(iCol,1)  = oa1ss(iCol)
             surface % oa4ss(iCol,2)  = oa2ss(iCol)
             surface % oa4ss(iCol,3)  = oa3ss(iCol)
             surface % oa4ss(iCol,4)  = oa4ss(iCol)
             surface % clxss(iCol,1)  = ol1ss(iCol)
             surface % clxss(iCol,2)  = ol2ss(iCol)
             surface % clxss(iCol,3)  = ol3ss(iCol)
             surface % clxss(iCol,4)  = ol4ss(iCol)
          end if
       end do
    end do

  end subroutine ufs_mpas_gwd_to_physics

  !> #########################################################################################
  !> Build the dycore-neutral reference pressure profile that UGWPv1 expects in ak/bk.
  !>
  !> FV3 passes its hybrid sigma-pressure coefficients to control_initialize(); MPAS is on a
  !> terrain-following height coordinate and has no such coefficients, so cires_ugwpv1_init
  !> would otherwise dereference null pointers at cires_ugwpv1_module.F90:252:
  !>
  !>     pmb(k) = ak(k) + pref*bk(k)        ! pref = 1.e5
  !>     zkm(k) = -hpskm*alog(pmb(k)/pref)
  !>
  !> Setting bk = 0 and ak = a reference layer pressure makes that evaluate pmb(k) = ak(k).
  !> The reference used here is the global mean of the MPAS base-state pressure, which is
  !> MPAS's own notion of a reference profile rather than an invented one.
  !>
  !> The mean is taken with a global reduction on purpose. pressure_base is (nVertLevels,
  !> nCells) because the terrain-following coordinate makes the base state terrain-dependent,
  !> so a per-rank mean or a single arbitrary column would change with the MPI decomposition
  !> and break decomposition-independence. Terrain spread in layer height is under 2% above
  !> 15 km, so the choice only perturbs UGWP's launch level by a level or two.
  !>
  !> Ordered surface -> TOA, matching both the MPAS bottom-up convention and the
  !> "ak -pa bk-dimensionless from surf" comment in cires_ugwpv1_module.F90.
  !> Called once, before MPAS_initialize().
  !> #########################################################################################
  subroutine ufs_mpas_reference_pressure(levs, ak, bk)
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_derived_types,   only : MPAS_LOG_CRIT
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array
    use mpas_dmpar,           only : mpas_dmpar_sum_real_array, mpas_dmpar_sum_int
    use mpas_log,             only : mpas_log_write

    ! Arguments
    integer,          intent(in)  :: levs
    real(kind=RKIND), intent(out) :: ak(levs+1)   !< reference pressure at layer centres (Pa)
    real(kind=RKIND), intent(out) :: bk(levs+1)   !< zero; MPAS is not on a hybrid coordinate

    ! Locals
    type(mpas_pool_type), pointer :: mesh_pool, diag_pool
    integer, pointer :: nCellsSolve, nVertLevels
    real(kind=RKIND), pointer :: pressure_base(:,:)
    real(kind=RKIND), allocatable :: localSum(:), globalSum(:)
    integer :: iCol, iLay, nCellsGlobal
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_reference_pressure'

    ! Access MPAS data pools
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh', mesh_pool)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag', diag_pool)
    call mpas_pool_get_dimension(mesh_pool, 'nCellsSolve', nCellsSolve)
    call mpas_pool_get_dimension(mesh_pool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(diag_pool, 'pressure_base', pressure_base)

    if (nVertLevels /= levs) then
       call mpas_log_write(subname // ' ERROR: nVertLevels does not match levs', &
                           messageType=MPAS_LOG_CRIT)
    end if
    if (.not. associated(pressure_base)) then
       call mpas_log_write(subname // ' ERROR: diag pressure_base is not associated; ' // &
                           'MPAS base state must be initialized before this call', &
                           messageType=MPAS_LOG_CRIT)
    end if

    allocate(localSum(levs), globalSum(levs))
    localSum(:) = 0.0_RKIND

    do iCol = 1, nCellsSolve
       do iLay = 1, levs
          localSum(iLay) = localSum(iLay) + pressure_base(iLay,iCol)
       end do
    end do

    call mpas_dmpar_sum_real_array(domain_ptr % dminfo, levs, localSum, globalSum)
    call mpas_dmpar_sum_int(domain_ptr % dminfo, nCellsSolve, nCellsGlobal)

    if (nCellsGlobal <= 0) then
       call mpas_log_write(subname // ' ERROR: global cell count is not positive', &
                           messageType=MPAS_LOG_CRIT)
    end if

    ak(1:levs) = globalSum(1:levs) / real(nCellsGlobal, RKIND)
    ! Model lid. cires_ugwpv1_module.F90:248 documents the top of the profile as zero
    ! pressure; UGWPv1 reads only 1:levs, so this element exists to satisfy the ak(levs+1)
    ! declaration in its interface.
    ak(levs+1) = 0.0_RKIND
    bk(:)      = 0.0_RKIND

    call mpas_log_write(subname // ': MPAS reference pressure profile from base state, ' // &
                        '$i global cells', intArgs=[nCellsGlobal])
    call mpas_log_write(subname // ': surface ak = $r Pa, model-top layer ak = $r Pa', &
                        realArgs=[ak(1), ak(levs)])

    deallocate(localSum, globalSum)
    nullify (mesh_pool)
    nullify (diag_pool)

  end subroutine ufs_mpas_reference_pressure

  !> #########################################################################################
  !> Procedure to populate MPAS diag_phys pool with CCPP data.
  !>
  !> The fields from diag_phys are allocated by MPAS, populated with CCPP Physics data, and
  !> written to output. So essentially we copy physics arrays back into MPAS memory to use
  !> MPAS's native output functionality.
  !> #########################################################################################
  subroutine ufs_mpas_phys_diag(control,radiation,diagnostics,tbd)
    use GFS_typedefs,         only : GFS_control_type
    use GFS_typedefs,         only : GFS_radtend_type
    use GFS_typedefs,         only : GFS_diag_type
    use GFS_typedefs,         only : GFS_tbd_type
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config

    ! Arguments
    type(GFS_control_type), intent(in) :: control
    type(GFS_radtend_type), intent(in) :: radiation
    type(GFS_diag_type),    intent(in) :: diagnostics
    type(GFS_tbd_type),     intent(in) :: tbd

    ! Locals
    type(mpas_pool_type), pointer :: diag_phys
    real(RKIND), pointer :: swdnb(:),swdnbc(:),swupb(:),swupbc(:)
    real(RKIND), pointer :: lwdnb(:),lwdnbc(:),lwupb(:),lwupbc(:)
    real(RKIND), pointer :: re_cloud(:,:),re_ice(:,:),re_snow(:,:)
    real(RKIND), pointer :: sfc_albedo(:),sfc_emiss(:)
    real(RKIND), pointer :: refl10cm(:,:)
    real(RKIND), pointer :: rainc(:),rainnc(:),frainnc(:),snownc(:),graupelnc(:)
    real(RKIND), pointer :: raincv(:),rainncv(:),snowncv(:),graupelncv(:)
    real(RKIND), pointer :: dusfcg(:),dvsfcg(:),dusfc_ls(:),dvsfc_ls(:)
    real(RKIND), pointer :: dusfc_bl(:),dvsfc_bl(:),dusfc_ss(:),dvsfc_ss(:)
    real(RKIND), pointer :: dusfc_fd(:),dvsfc_fd(:)
    real(RKIND), pointer :: dtaux3d(:,:), dtauy3d(:,:)
    real(RKIND), pointer :: dtaux3d_ls(:,:), dtauy3d_ls(:,:), dtaux3d_ss(:,:), dtauy3d_ss(:,:)
    real(RKIND), pointer :: dtaux3d_fd(:,:), dtauy3d_fd(:,:), dtaux3d_bl(:,:), dtauy3d_bl(:,:)
    integer,     pointer :: nThreads, cellSolveThreadStart(:), cellSolveThreadEnd(:)
    integer :: iCol, ithread
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_phys_diag'

    ! Get openMP information
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'nThreads',             nThreads)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadStart', cellSolveThreadStart)
    call mpas_pool_get_dimension(domain_ptr % blocklist % dimensions,  'cellSolveThreadEnd',   cellSolveThreadEnd)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics',  diag_phys)

    ! Grab fields from MPAS pools
    call mpas_pool_get_array(diag_phys,'swdnb'     , swdnb     )
    call mpas_pool_get_array(diag_phys,'swdnbc'    , swdnbc    )
    call mpas_pool_get_array(diag_phys,'swupb'     , swupb     )
    call mpas_pool_get_array(diag_phys,'swupbc'    , swupbc    )
    call mpas_pool_get_array(diag_phys,'lwdnb'     , lwdnb     )
    call mpas_pool_get_array(diag_phys,'lwdnbc'    , lwdnbc    )
    call mpas_pool_get_array(diag_phys,'lwupb'     , lwupb     )
    call mpas_pool_get_array(diag_phys,'lwupbc'    , lwupbc    )
    call mpas_pool_get_array(diag_phys,'refl10cm'  , refl10cm  )
    call mpas_pool_get_array(diag_phys,'rainc'     , rainc     )
    call mpas_pool_get_array(diag_phys,'rainnc'    , rainnc    )
    call mpas_pool_get_array(diag_phys,'frainnc'   , frainnc   )
    call mpas_pool_get_array(diag_phys,'snownc'    , snownc    )
    call mpas_pool_get_array(diag_phys,'graupelnc' , graupelnc )
    call mpas_pool_get_array(diag_phys,'raincv'    , raincv    )
    call mpas_pool_get_array(diag_phys,'rainncv'   , rainncv   )
    call mpas_pool_get_array(diag_phys,'snowncv'   , snowncv   )
    call mpas_pool_get_array(diag_phys,'graupelncv', graupelncv)
    call mpas_pool_get_array(diag_phys,'re_cloud'  , re_cloud  )
    call mpas_pool_get_array(diag_phys,'re_ice'    , re_ice    )
    call mpas_pool_get_array(diag_phys,'re_snow'   , re_snow   )
    call mpas_pool_get_array(diag_phys,'sfc_albedo', sfc_albedo)
    call mpas_pool_get_array(diag_phys,'sfc_emiss' , sfc_emiss )
    ! UFS GWD diagnostics are conditionally allocated.
    if (control % ldiag_ugwp .or. control % do_ugwp_v1) then
       call mpas_pool_get_array(diag_phys,'dusfcg'    , dusfcg    )
       call mpas_pool_get_array(diag_phys,'dvsfcg'    , dvsfcg    )
       call mpas_pool_get_array(diag_phys,'dusfc_ls'  , dusfc_ls  )
       call mpas_pool_get_array(diag_phys,'dvsfc_ls'  , dvsfc_ls  )
       call mpas_pool_get_array(diag_phys,'dusfc_bl'  , dusfc_bl  )
       call mpas_pool_get_array(diag_phys,'dvsfc_bl'  , dvsfc_bl  )
       call mpas_pool_get_array(diag_phys,'dusfc_ss'  , dusfc_ss  )
       call mpas_pool_get_array(diag_phys,'dvsfc_ss'  , dvsfc_ss  )
       call mpas_pool_get_array(diag_phys,'dusfc_fd'  , dusfc_fd  )
       call mpas_pool_get_array(diag_phys,'dvsfc_fd'  , dvsfc_fd  )
       call mpas_pool_get_array(diag_phys,'dtaux3d'   , dtaux3d   )
       call mpas_pool_get_array(diag_phys,'dtauy3d'   , dtauy3d   )
       call mpas_pool_get_array(diag_phys,'dtaux3d_ls', dtaux3d_ls)
       call mpas_pool_get_array(diag_phys,'dtauy3d_ls', dtauy3d_ls)
       call mpas_pool_get_array(diag_phys,'dtaux3d_ss', dtaux3d_ss)
       call mpas_pool_get_array(diag_phys,'dtauy3d_ss', dtauy3d_ss)
       call mpas_pool_get_array(diag_phys,'dtaux3d_bl', dtaux3d_bl)
       call mpas_pool_get_array(diag_phys,'dtauy3d_bl', dtauy3d_bl)
       call mpas_pool_get_array(diag_phys,'dtaux3d_fd', dtaux3d_fd)
       call mpas_pool_get_array(diag_phys,'dtauy3d_fd', dtauy3d_fd)
    end if

    do ithread = 1,nThreads
       do iCol = cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread)
          ! Radiation fluxes at surface
          swdnb(iCol)  = radiation%sfcfsw(iCol)%dnfxc
          swdnbc(iCol) = radiation%sfcfsw(iCol)%dnfx0
          swupb(iCol)  = radiation%sfcfsw(iCol)%upfxc
          swupbc(iCol) = radiation%sfcfsw(iCol)%upfx0
          lwdnb(iCol)  = radiation%sfcflw(iCol)%dnfxc
          lwdnbc(iCol) = radiation%sfcflw(iCol)%dnfx0
          lwupb(iCol)  = radiation%sfcflw(iCol)%upfxc
          lwupbc(iCol) = radiation%sfcflw(iCol)%upfx0
          ! Reflectivity
          refl10cm(:,iCol) = diagnostics%refl_10cm(iCol,:)
          ! Instantaneous precipitation
          raincv(iCol)     = diagnostics%rain(iCol)
          rainncv(iCol)    = diagnostics%rainc(iCol)
          snowncv(iCol)    = diagnostics%snow(iCol)
          graupelncv(iCol) = diagnostics%graupel(iCol)
          ! Accumulated precipitation
          rainc(iCol)      = diagnostics%cnvprcp(iCol)
          rainnc(iCol)     = diagnostics%totprcp(iCol)
          frainnc(iCol)    = diagnostics%totice(iCol)
          snownc(iCol)     = diagnostics%totsnw(iCol)
          graupelnc(iCol)  = diagnostics%totgrp(iCol)
          ! Hydrometeor effective radii
          re_cloud(:,iCol) = tbd%phy_f3d(iCol,:,control%nleffr)
          re_ice(:,iCol)   = tbd%phy_f3d(iCol,:,control%nieffr)
          re_snow(:,iCol)  = tbd%phy_f3d(iCol,:,control%nseffr)
          ! Surface radiative properties
          sfc_albedo(iCol) = radiation%sfalb(iCol)
          sfc_emiss(iCol)  = radiation%semis(iCol)
          ! Gravity-wave physics diagnostics (conditionally allocated)
          if (control % ldiag_ugwp .or. control % do_ugwp_v1) then
             dusfcg(iCol)       = diagnostics%dusfcg(iCol)
             dvsfcg(iCol)       = diagnostics%dvsfcg(iCol)
             dusfc_ls(iCol)     = diagnostics%du_ogwcol(iCol)
             dvsfc_ls(iCol)     = diagnostics%dv_ogwcol(iCol)
             dusfc_bl(iCol)     = diagnostics%du_oblcol(iCol)
             dvsfc_bl(iCol)     = diagnostics%dv_oblcol(iCol)
             dusfc_ss(iCol)     = diagnostics%du_osscol(iCol)
             dvsfc_ss(iCol)     = diagnostics%dv_osscol(iCol)
             dusfc_fd(iCol)     = diagnostics%du_ofdcol(iCol)
             dvsfc_fd(iCol)     = diagnostics%dv_ofdcol(iCol)
             dtaux3d(:,iCol)    = diagnostics%dudt_gw(iCol,:)
             dtauy3d(:,iCol)    = diagnostics%dvdt_gw(iCol,:)
             dtaux3d_ls(:,iCol) = diagnostics%dudt_ogw(iCol,:)
             dtauy3d_ls(:,iCol) = diagnostics%dvdt_ogw(iCol,:)
             dtaux3d_ss(:,iCol) = diagnostics%dudt_oss(iCol,:)
             dtauy3d_ss(:,iCol) = diagnostics%dvdt_oss(iCol,:)
             dtaux3d_bl(:,iCol) = diagnostics%dudt_obl(iCol,:)
             dtauy3d_bl(:,iCol) = diagnostics%dvdt_obl(iCol,:)
             dtaux3d_fd(:,iCol) = diagnostics%dudt_ofd(iCol,:)
             dtauy3d_fd(:,iCol) = diagnostics%dvdt_ofd(iCol,:)
          end if
       end do
    end do
  end subroutine ufs_mpas_phys_diag

  !> ########################################################################################
  !> This routine computes physics surface properties using data from the LANDUSE.TBL (read in
  !> during model initialization), and MPAS grid information.
  !> Called once per calendar day.
  !> ########################################################################################
  subroutine ufs_mpas_landuse_update(julday)
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_subpool, mpas_pool_get_dimension
    use mpas_log,             only : mpas_log_write
    use mpas_derived_types,   only : MPAS_LOG_ERR, MPAS_LOG_WARN, MPAS_LOG_CRIT
    use mpas_pool_routines,   only : mpas_pool_get_array, mpas_pool_get_config
    implicit none
    ! Arguments
    integer, intent(in) :: julday
    ! Locals
    type(mpas_pool_type), pointer :: sfc_input, diag_phys, mesh
    real(kind=RKIND), pointer :: latCell(:),snoalb(:),snowc(:),xice(:),albbck(:),embck(:),&
         xicem(:),xland(:),z0(:),mavail(:),sfc_albedo(:),sfc_emiss(:),thc(:),ust(:),znt(:)
    integer, pointer :: nCells,isice,iswater,ivgtyp(:),landmask(:)
    logical,pointer :: config_frac_seaice, config_sfc_albedo
    integer :: iCell,isn,is
    real(kind=RKIND) :: xice_threshold
    character(len=*), parameter :: subname = 'atmos_coupling::ufs_mpas_landuse_update'

    call mpas_log_write(subname //'   Updating MPAS surface properting', messageType=MPAS_LOG_WARN)

    ! Are these really physics flags? (i.e. Do the scheme need theses? Or just here for the surface init?)
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_frac_seaice', config_frac_seaice)
    call mpas_pool_get_config( domain_ptr % blocklist % configs, 'config_sfc_albedo' , config_sfc_albedo)

    ! Access MPAS data pools.
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'diag_physics', diag_phys)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'mesh',         mesh)
    call mpas_pool_get_subpool(domain_ptr % blocklist % structs, 'sfc_input',    sfc_input)

    ! Dimensions
    call mpas_pool_get_dimension(mesh,'nCells',nCells)

    ! Arrays
    call mpas_pool_get_array(mesh,     'latCell'   ,latCell    )
    call mpas_pool_get_array(sfc_input,'isice_lu'  ,isice      )
    call mpas_pool_get_array(sfc_input,'iswater_lu',iswater    )
    call mpas_pool_get_array(sfc_input,'landmask'  ,landmask   )
    call mpas_pool_get_array(sfc_input,'ivgtyp'    ,ivgtyp     )
    call mpas_pool_get_array(sfc_input,'snoalb'    ,snoalb     )
    call mpas_pool_get_array(sfc_input,'snowc'     ,snowc      )
    call mpas_pool_get_array(sfc_input,'xice'      ,xice       )
    call mpas_pool_get_array(sfc_input,'xland'     ,xland      )
    call mpas_pool_get_array(sfc_input,'sfc_albbck',albbck     )
    call mpas_pool_get_array(diag_phys,'sfc_emibck',embck      )
    call mpas_pool_get_array(diag_phys,'mavail'    ,mavail     )
    call mpas_pool_get_array(diag_phys,'sfc_albedo',sfc_albedo )
    call mpas_pool_get_array(diag_phys,'sfc_emiss' ,sfc_emiss  )
    call mpas_pool_get_array(diag_phys,'thc'       ,thc        )
    call mpas_pool_get_array(diag_phys,'ust'       ,ust        )
    call mpas_pool_get_array(diag_phys,'xicem'     ,xicem      )
    call mpas_pool_get_array(diag_phys,'z0'        ,z0         )
    call mpas_pool_get_array(diag_phys,'znt'       ,znt        )

    if(.not. config_frac_seaice) then
       xice_threshold = 0.5_RKIND
    elseif(config_frac_seaice) then
       xice_threshold = 0.02
    endif

    do iCell = 1, nCells
       !finds the season as function of julian day (summer=1, winter=2): summer in the Northern
       !Hemisphere is defined between March 15th and September 15th (winter otherwise).
       isn = 1
       if(julday.lt.105 .or. julday.ge.288) isn=2
       if(latCell(iCell) .lt. 0.) isn=3-isn

       is = ivgtyp(iCell)

       !set no data points to water:
       if(is.eq.0) is = iswater
       if(.not. config_sfc_albedo) albbck(iCell) = albd(is,isn)/100.
       sfc_albedo(iCell) = albbck(iCell)

       if(snowc(iCell) .gt. 0.5) then
          albbck(iCell) = albd(isice,isn) / 100.
          if(config_sfc_albedo) then
             sfc_albedo(iCell) = snoalb(iCell)
          else
             sfc_albedo(iCell) = albbck(iCell) / (1+scfx(is,isn))
          endif
       endif
       thc(iCell)    = therin(is,isn) / 100.
       z0(iCell)     = sfz0(is,isn) / 100.
       znt(iCell)    = z0(icell)
       embck(iCell)  = sfem(is,isn)
       sfc_emiss(iCell) = embck(iCell)

       if(associated(mavail)) mavail(iCell) = slmo(is,isn)
       if(associated(ust)) ust(iCell) = 0.0001

       !set sea-ice points to land with ice/snow surface properties:
       if(xice(iCell) .ge. xice_threshold) then
          albbck(iCell) = albd(isice,isn) / 100.
          embck(iCell)  = sfem(isice,isn)
          if(config_frac_seaice) then
             !0.08 is the albedo over open water.
             !0.98 is the emissivity over open water.
             sfc_albedo(iCell) = xice(iCell)*albbck(iCell) + (1-xice(iCell))*0.08
             sfc_emiss(iCell)  = xice(iCell)*embck(iCell)  + (1-xice(iCell))*0.98
          else
             sfc_albedo(iCell) = albbck(iCell)
             sfc_emiss(iCell)  = embck(iCell)
          endif
          thc(iCell) = therin(isice,isn) / 100.
          z0(icell)  = sfz0(isice,isn) / 100.
          znt(iCell) = z0(iCell)

          if(associated(mavail)) mavail(iCell) = slmo(isice,isn)
       endif
    end do ! nCells

    call mpas_log_write(subname //'   Finished updating MPAS surface properties ', messageType=MPAS_LOG_WARN)
  end subroutine ufs_mpas_landuse_update

  !> ########################################################################################
  !> Interpolate wind-tendencies from centers to edges of grid-cells.
  !>
  !> See MPAS-Model/src/core_atmosphere/physics/mpas_atmphys_todynamics.F: tend_toEdges
  !> ########################################################################################
  subroutine tend_toEdges(mesh,Ux_tend,Uy_tend,U_tend)
    use mpas_derived_types,   only : mpas_pool_type
    use mpas_pool_routines,   only : mpas_pool_get_dimension, mpas_pool_get_array

    type(mpas_pool_type),intent(in):: mesh
    real(kind=RKIND),intent(in),dimension(:,:),target:: Ux_tend,Uy_tend
    real(kind=RKIND),intent(out),dimension(:,:):: U_tend

    ! locals
    integer:: iCell,iEdge,k,j
    integer:: cell1, cell2
    integer,pointer:: nCells,nCellsSolve,nEdges
    integer,dimension(:,:),pointer:: cellsOnEdge
    real(kind=RKIND),dimension(:,:),pointer:: east,north,edgeNormalVectors

    ! Grab dimensions
    call mpas_pool_get_dimension(mesh,'nCells',nCells)
    call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
    call mpas_pool_get_dimension(mesh,'nEdges',nEdges)

    ! Grad arrays
    call mpas_pool_get_array(mesh,'east',east)
    call mpas_pool_get_array(mesh,'north',north)
    call mpas_pool_get_array(mesh,'edgeNormalVectors',edgeNormalVectors)
    call mpas_pool_get_array(mesh,'cellsOnEdge',cellsOnEdge)

    do iEdge = 1, nEdges
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       U_tend(:,iEdge) =  Ux_tend(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell1)   &
                                                 +  edgeNormalVectors(2,iEdge) * east(2,cell1)   &
                                                 +  edgeNormalVectors(3,iEdge) * east(3,cell1))  &
                        + Uy_tend(:,cell1) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell1)  &
                                                 +  edgeNormalVectors(2,iEdge) * north(2,cell1)  &
                                                 +  edgeNormalVectors(3,iEdge) * north(3,cell1)) &
                        + Ux_tend(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * east(1,cell2)   &
                                                 +  edgeNormalVectors(2,iEdge) * east(2,cell2)   &
                                                 +  edgeNormalVectors(3,iEdge) * east(3,cell2))  &
                        + Uy_tend(:,cell2) * 0.5 * (edgeNormalVectors(1,iEdge) * north(1,cell2)  &
                                                 +  edgeNormalVectors(2,iEdge) * north(2,cell2)  &
                                                 +  edgeNormalVectors(3,iEdge) * north(3,cell2))
    end do
  end subroutine tend_toEdges
  
end module atmos_coupling_mod
