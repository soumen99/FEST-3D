! include in paralle module, send_recv_i subroutine

      select case (turbulence)

          case ("none")
              !include nothing
              continue

          case ("sst")
              include "turbulence_models/sst/parallel/recv_imin.inc"

          case DEFAULT
              call dmsg(5, "turbulence_models/include/parallel", "recv_imin", &
                        "ERROR: Turbulence model not recognised")
              STOP

      end select