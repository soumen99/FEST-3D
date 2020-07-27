module diffusivity
! Setup, destroy, calculate molecular diffusivity


use vartypes
use viscosity     , only : mu
use copy_bc       , only : copy1
use utils         , only : alloc


implicit none
  private
  real(wp), dimension(:, :, :), allocatable, target     :: diff
  public :: setup_viscosity
  public :: calculate_viscosity
  public :: diff

  contains

    subroutine calculate_diffusivity(qp, scheme, flow, bc, dims)
      !< Calculate molecular and turbulent viscosity
      implicit none
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(boundarytype), intent(in) :: bc
      !< boundary conditions and fixed values
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center


    integer i,j,k
    real(wp) :: pressure
    real(wp) :: density
    real(wp) :: T
    integer :: imx, jmx, kmx

    imx = dims%imx
    jmx = dims%jmx
    kmx= dims%kmx

    ! calculate molecular diffusivity
    if (flow%diff_ref /=0.) then
        select case (trim(flow%diff_variation))
        case ('constant')

            continue

            case ('Einstein-Stokes law')
                do k = 0,kmx
                    do j = 0,jmx
                        do i = 0,imx
                            density = qp(i,j,k,1)
                            pressure = qp(i,j,k,5)
                            T = pressure/(density*flow%R_gas)
                            diff(i,j,k) = ((flow%diff_ref)*(flow%mu_ref)*T)/(flow%T_ref*mu(i,j,k))
                        end do
                    end do
                end do

                case DEFAULT
            print*,"diffusivity variation not recognized:"
            print*, "   found '",trim(flow%diff_variation),"'"
            print*, "accepted values: 1) Einstein-Stokes law"
            print*, "                 2) constant"
            Fatal_error


        end select
    end if


    if(any(isnan(diff))then
        Fatal_error
      end if

    end subroutine calculate_diffusivity




    subroutine setup_diffusivity(scheme,flow, dims)
      !< Allocate and pointer for molecular and turbulent viscosity
      implicit none
      type(extent), intent(in) :: dims
      !< Extent of the domain:imx,jmx,kmx
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(flowtype), intent(in) :: flow
      !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
      integer :: imx, jmx, kmx

      imx = dims%imx
      jmx = dims%jmx
      kmx = dims%kmx
      !setup_molecular_diffusivity()
      if (flow%diff_ref/=0.) then
        call alloc(diff, -2, imx+2, -2, jmx+2, -2, kmx+2)
        diff = flow%diff_ref !intialize
      end if

    end subroutine setup_diffusivity


end module diffusivity



