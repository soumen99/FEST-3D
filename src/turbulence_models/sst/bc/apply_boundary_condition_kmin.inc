!module boundary_conditions
!TODO Index 6 and 7 are wrong , check logrith definintion in main BC file

                    if ((bc_kmn(i, j) .and. current_condition) .eq. sst_BC_FIX_tk) then
                        call fix_tk("kmin", i * cell_ind, bc_imn_fix_values(i, j, 6))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. sst_BC_FIX_tw) then
                        call fix_tw("kmin", i * cell_ind, bc_imn_fix_values(i, j, 7))
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. sst_BC_COPY_tk) then
                        call copy_tk("kmin", i * cell_ind)
                    else if ((bc_kmn(i, j) .and. current_condition) .eq. sst_BC_COPY_tw) then
                        call copy_tw("kmin", i * cell_ind)
                    end if
