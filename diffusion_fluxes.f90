  module scalar_diffusion

  #include"error.h"
  use vartypes
  use gradients  , only : gradphi_x
  use gradients  , only : gradphi_y
  use gradients  , only : gradphi_z
  use diffusivity, only : diff
  use gradient_diffusion, only : sc
  use viscosity,   only : mu
  use viscosity,   only : mu_t
  use utils      , only : alloc
  implicit none
  private

  integer :: imx, jmx, kmx

  public :: compute_viscous_fluxes

  contains

    subroutine compute_diffusion_fluxes(F, G, H, qp, cells, Ifaces, Jfaces, Kfaces, scheme, flow, dims)
      !< Call to all diffusion flux subroutine based on


        implicit none
        type(schemetype), intent(in) :: scheme
        !< finite-volume Schemes
        type(flowtype), intent(in) :: flow
        !< Information about fluid flow: freestream-speed, ref-viscosity,etc.
        type(extent), intent(in) :: dims
        !< Extent of the domain:imx,jmx,kmx
        real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
        real(wp), dimension(:, :, :, :), intent(inout) :: F, G, H
        type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
        !< Input cell quantities: volume
        type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
        !< Store face quantites for I faces
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
        !< Store face quantites for J faces
        type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
        !< Store face quantites for K faces
        integer, dimension(3) :: flagsi=(/1,0,0/)
        integer, dimension(3) :: flagsj=(/0,1,0/)
        integer, dimension(3) :: flagsk=(/0,0,1/)

        imx = dims%imx
        jmx = dims%jmx
        kmx = dims%kmx

        call compute_diffusion_fluxes_scalar(F, qp, cells, Ifaces, flagsi, scheme, flow, dims)
        call compute_diffusion_fluxes_scalar(G, qp, cells, Jfaces, flagsj, scheme, flow, dims)
        call compute_diffusion_fluxes_scalar(H, qp, cells, Kfaces, flagsk,scheme, flow, dims)


        if (any(isnan(G))) then
              Fatal_error
            end if
            if (any(isnan(F))) then
              Fatal_error
            end if
            if (any(isnan(H))) then
              Fatal_error
            end if

    end subroutine compute_viscous_fluxes


    subroutine compute_diffusion_fluxes_scalar(F, qp, cells, faces, flags, dims)
      !< Compute viscous fluxes for scalar transport equation
      implicit none
      real(wp), dimension(:, :, :, :), intent(inout) :: F
      !< Flux array
      type(schemetype), intent(in) :: scheme
      !< finite-volume Schemes
      type(extent), intent(in) :: dims
      !< Extent of the domain: imx,jmx,kmx
      integer, dimension(3), intent(in) :: flags
      !< flags for direction switch
      real(wp), dimension(-2:dims%imx+2, -2:dims%jmx+2, -2:dims%kmx+2, 1:dims%n_var), intent(in) :: qp
      !< Store primitive variable at cell center
      type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
      !< Input cell quantities: volume
      type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
      !< Face quantities: area and unit normal
      ! local variables
      real(wp) :: rhoface
      real(wp) :: normal_comp
      real(wp) :: d_LR
      real(wp) :: delx
      real(wp) :: dely
      real(wp) :: delz
      real(wp) :: mut_f
      real(wp) :: diff_f
      real(wp) :: delphi
      real(wp) :: nx
      real(wp) :: ny
      real(wp) :: nz
      real(wp) :: area
      integer :: i, j, k
      integer :: ii, jj, kk
      !--- scalar variable requirement ---!
      real(wp) :: dphidx, dphidy, dphidz



      ii = flags(1)
      jj = flags(2)
      kk = flags(3)

      !---------------------------------------------------------------------
      ! Calculating the scalar diffusive fluxes at the faces
      !--------------------------------------------------------------------
      do k = 1, dims%kmx - 1 + kk
       do j = 1, dims%jmx - 1 + jj
        do i = 1, dims%imx - 1 + ii

          !--- FACE Gradients ---!
          ! Gradients at face as average of gradients at cell centres
          dphidx = 0.5 * (gradphi_x(i-ii, j-jj, k-kk) + gradphi_x(i, j, k))
          dphidy = 0.5 * (gradphi_y(i-ii, j-jj, k-kk) + gradphi_y(i, j, k))
          dphidz = 0.5 * (gradphi_z(i-ii, j-jj, k-kk) + gradphi_z(i, j, k))

          !--- For ODD-EVEN coupling error ---!
          ! distance between cell center of adjacent cell for the i,j,k face
          delx = cells(i, j, k)%centerx - cells(i-ii, j-jj, k-kk)%centerx
          dely = cells(i, j, k)%centery - cells(i-ii, j-jj, k-kk)%centery
          delz = cells(i, j, k)%centerz - cells(i-ii, j-jj, k-kk)%centerz

          d_LR = sqrt(delx*delx + dely*dely + delz*delz)

          ! difference in state across face
          delphi = qp(i, j, k, n_var) - qp(i-ii, j-jj, k-kk, n_var)

          !normal_comp   = ( delta(phi) - (grad(phi).dot.delR) )/magnitudeR
          !new grad(phi) =  grad(phi) + correction(normal_comp.dot.delR/magnitudeR)
          normal_comp = (delphi - (dphidx*delx + dphidy*dely + dphidz*delz))/d_LR
          dphidx       =  dphidx + (normal_comp * delx / d_LR)
          dphidy       =  dphidy + (normal_comp * dely / d_LR)
          dphidz       =  dphidz + (normal_comp * delz / d_LR)
          !--- end of ODD-EVEN coupling correction ---!

          diff_f = 0.5*(diff(i-ii, j-jj, k-kk)+ diff(i, j, k))
          if(trim(scheme%turbulence)/='none') then
            mut_f = 0.5*(mu_t(i-ii, j-jj, k-kk) + mu_t(i, j, k))
          else
            mut_f = 0.0
          end if


          ! calling some element from memory and keep them handy for calculation
          nx    = faces(i,j,k)%nx
          ny    = faces(i,j,k)%ny
          nz    = faces(i,j,k)%nz
          area  = faces(i,j,k)%A

          ! adding viscous fluxes to stored convective flux
          F(i, j, k, dims%n_var) = F(i, j, k, dims%n_var) - area*((diff_f+mut_f/sc)*(dphidx*nx + dphidy*ny + dphidz*nz))

        end do
       end do
      end do
    end subroutine compute_diffusion_fluxes_scalar

