module scheme

    use global, only: SCHEME_NAME_LENGTH

    use global_vars, only : imx
    use global_vars, only : jmx
    use global_vars, only : kmx

    use global_vars, only : mu_ref
    use global_vars, only : n_var
    use global_vars, only : F_p
    use global_vars, only : G_p
    use global_vars, only : H_p
    use global_vars, only : mass_residue
    use global_vars, only : x_mom_residue
    use global_vars, only : y_mom_residue
    use global_vars, only : z_mom_residue
    use global_vars, only : energy_residue
    use global_vars, only : TKE_residue
    use global_vars, only : omega_residue
    use global_vars, only : turbulence

    use global_vars, only : scheme_name

    use utils, only: alloc, dealloc, dmsg
!    use grid, only: imx, jmx, kmx
!    use state, only: n_var, mu_ref
    use face_interpolant, only: setup_interpolant_scheme, &
            destroy_interpolant_scheme
    use van_leer, only: &
            setup_scheme_van_leer => setup_scheme, &
            destroy_scheme_van_leer => destroy_scheme, &
            compute_fluxes_van_leer => compute_fluxes, &
            get_residue_van_leer => get_residue, &
            F_van_leer => F, &
            G_van_leer => G, &
            H_van_leer => H, &
            residue_van_leer => residue
    use ausm, only: &
            setup_scheme_ausm => setup_scheme, &
            destroy_scheme_ausm => destroy_scheme, &
            compute_fluxes_ausm => compute_fluxes, &
            get_residue_ausm => get_residue, &
            F_ausm => F, &
            G_ausm => G, &
            H_ausm => H, &
            residue_ausm => residue
    use ldfss0, only: &
            setup_scheme_ldfss0=> setup_scheme, &
            destroy_scheme_ldfss0=> destroy_scheme, &
            compute_fluxes_ldfss0 => compute_fluxes, &
            get_residue_ldfss0 => get_residue, &
            F_ldfss0 => F, &
            G_ldfss0 => G, &
            H_ldfss0 => H, &
            residue_ldfss0 => residue
!    use source, only: setup_source, destroy_source
!   use hlle, only: &
!           setup_scheme_hlle => setup_scheme, &
!           destroy_scheme_hlle => destroy_scheme, &
!           get_residue_hlle => get_residue
    include "turbulence_models/include/scheme/import_module.inc"

    implicit none
    private

!    character(len=SCHEME_NAME_LENGTH) :: scheme_name
!    real, public, dimension(:, :, :, :), pointer :: F_p, G_p, H_p
!    real, public, dimension(:, :, :), pointer :: mass_residue, x_mom_residue, &
!                             y_mom_residue, z_mom_residue, energy_residue

    include "turbulence_models/include/scheme/variable_deceleration.inc" 

    ! Public members
    public :: scheme_name
    public :: setup_scheme
    public :: destroy_scheme
    public :: compute_fluxes
    public :: compute_residue

    contains

        subroutine setup_scheme()
            implicit none

            call setup_interpolant_scheme()

            select case (scheme_name)
                case ("van_leer")
                    call setup_scheme_van_leer()
                    F_p => F_van_leer
                    G_p => G_van_leer
                    H_p => H_van_leer
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_van_leer(:, :, :, 5)
                    include "turbulence_models/include/scheme/van_leer_setup.inc" 
                case ("ausm")
                    call setup_scheme_ausm()
                    F_p => F_ausm
                    G_p => G_ausm
                    H_p => H_ausm
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ausm(:, :, :, 5)
                    include "turbulence_models/include/scheme/ausm_setup.inc" 
                case ("ldfss0")
                    call setup_scheme_ldfss0()
                    F_p => F_ldfss0
                    G_p => G_ldfss0
                    H_p => H_ldfss0
                    mass_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 1)
                    x_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 2)
                    y_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 3)
                    z_mom_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 4)
                    energy_residue(1:imx-1, 1:jmx-1, 1:kmx-1) => residue_ldfss0(:, :, :, 5)
                    include "turbulence_models/include/scheme/ldfss0_setup.inc" 
!               case ("hlle")
!                   call setup_scheme_hlle()
                case default
                    call dmsg(5, 'scheme', 'setup_scheme', &
                            'Scheme not recognized.')
                    stop
            end select

 !           if (mu_ref /= 0.0) then
  !              call setup_source()
  !          end if

        end subroutine setup_scheme

        subroutine deallocate_memory()

            implicit none

            nullify(F_p)
            nullify(G_p)
            nullify(H_p)
            nullify(mass_residue)
            nullify(x_mom_residue)
            nullify(y_mom_residue)
            nullify(z_mom_residue)
            nullify(energy_residue)
            include "turbulence_models/include/scheme/deallocate_memory.inc" 

        end subroutine deallocate_memory

        subroutine destroy_scheme()

            implicit none

            select case (scheme_name)
                case ("van_leer")
                    call destroy_scheme_van_leer
                case ("ausm")
                    call destroy_scheme_ausm()
                case ("ldfss0")
                    call destroy_scheme_ldfss0()
!               case ("hlle")
!                   call destroy_scheme_hlle()
                case default
                    call dmsg(5, 'scheme', 'destroy_scheme', &
                            'Scheme not recognized.')
                    stop
            end select
            
!            if (mu_ref /= 0.0) then
!                call destroy_source()
!            end if

            call destroy_interpolant_scheme()
            call deallocate_memory()

        end subroutine destroy_scheme

        subroutine compute_fluxes
        
            implicit none

            select case (scheme_name)
                case ("van_leer")
                    call compute_fluxes_van_leer()
                case ("ausm")
                    call compute_fluxes_ausm()
                case ("ldfss0")
                    call compute_fluxes_ldfss0()
!               case ("hlle")
!                   call compute_fluxes_hlle()
                case default
                    call dmsg(5, 'scheme', 'compute_fluxes', &
                            'Scheme not recognized.')
                    stop
            end select
            
        end subroutine compute_fluxes

        subroutine compute_residue()
            
            implicit none
            
            select case (scheme_name)
                case ("van_leer")
                    call get_residue_van_leer()
                case ("ausm")
                    call get_residue_ausm()
                case ("ldfss0")
                    call get_residue_ldfss0()
!               case ("hlle")
!                   call get_residue_hlle()
                case default
                    call dmsg(5, 'scheme', 'compute_residue', &
                            'Scheme not recognized.')
                    stop
            end select
        
        end subroutine compute_residue

end module scheme
