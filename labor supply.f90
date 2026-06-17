!##############################################################################
! PROGRAM LaborSupplyCollocation
!
! ## Extension of the life cycle model by endogenous labor supply
!    — Collocation method (Mongey Case I)
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
! Translated from MATLAB labor_supply.m collocation solver.
!##############################################################################
module collocation_globals

    use toolbox

    implicit none

    ! =========================================================================
    ! Model parameters
    ! =========================================================================
    integer, parameter :: JJ = 80       ! number of ages
    integer, parameter :: JR = 45       ! retirement age
    integer, parameter :: NP = 2        ! number of permanent productivity states
    integer, parameter :: NS = 7        ! number of persistent shock states
    integer, parameter :: NA = 200      ! number of asset grid points

    ! household preferences
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam  = 1d0 - 1d0/gamma
    real*8, parameter :: beta  = 0.98d0
    real*8, parameter :: nu    = 0.335d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

    ! asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 200d0
    real*8, parameter :: a_grow = 0.05d0

    ! net prices
    real*8, parameter :: r_val = 0.04d0
    real*8, parameter :: w_val = 1.0d0

    ! =========================================================================
    ! Collocation dimensions
    ! =========================================================================
    integer, parameter :: nx   = NA + 1         ! no. of collocation nodes (= no. basis fns)
    integer, parameter :: nxd  = nx             ! distribution grid = collocation grid
    integer, parameter :: ndim = nx * NS        ! total combined state dimension

    integer, parameter :: is_initial = 4        ! middle eta state

    ! =========================================================================
    ! Grids and stochastic processes
    ! =========================================================================
    real*8 :: a(0:NA)                    ! asset grid
    real*8 :: eta(NS)                    ! persistent shock values
    real*8 :: pi_mat(NS, NS)             ! shock transition matrix
    real*8 :: theta(NP)                  ! permanent productivity
    real*8 :: dist_theta(NP)             ! initial distribution over theta
    real*8 :: eff(JJ)                    ! age-efficiency profile
    real*8 :: psi(JJ+1)                  ! conditional survival probabilities
    real*8 :: pen(JJ)                    ! pension payments

    ! =========================================================================
    ! Precomputed arrays
    ! =========================================================================
    real*8 :: cash_arr(JJ, 0:NA, NP, NS)    ! (1+r)*a(ia) + pen(j)
    real*8 :: wage_arr(JJ, NP, NS)          ! w*eff(j)*theta(ip)*eta(is)

    ! =========================================================================
    ! Solution arrays  (age, asset_node, productivity, shock)
    ! =========================================================================
    real*8 :: V(JJ, 0:NA, NP, NS), EV(JJ, 0:NA, NP, NS)
    real*8 :: c(JJ, 0:NA, NP, NS), l(JJ, 0:NA, NP, NS)
    real*8 :: aplus(JJ, 0:NA, NP, NS)
    real*8 :: residual(JJ), ev_error(JJ)

    ! =========================================================================
    ! Distribution
    ! =========================================================================
    real*8 :: phi_dist(0:NA, NP, NS, JJ)    ! mass at (ia, ip, is, ij)
    real*8 :: mass(JJ, NP), mass_error(JJ)

    ! =========================================================================
    ! Cohort aggregates
    ! =========================================================================
    real*8 :: c_coh(JJ), y_coh(JJ), l_coh(JJ), h_coh(JJ)
    real*8 :: a_coh(JJ), v_coh(JJ)
    real*8 :: pen_coh(JJ), income_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_l(JJ), cv_h(JJ), corr_hl(JJ)
    real*8 :: frac_bor(JJ)

    ! =========================================================================
    ! Communication variables for optimisation objective
    ! =========================================================================
    integer :: ij_com, ip_com, is_com   ! current age / productivity / shock
    real*8  :: gold_coefEV(0:NA)        ! next-period EV coefficients for current (ip,is)

    ! =========================================================================
    ! Function approximation space for the asset dimension
    ! =========================================================================
    type(func_space) :: fspace_a   ! linear spline space on a(0:NA)

    ! =========================================================================
    ! PEV operator  (pi ⊗ I_nx), constructed once in initialisation
    ! =========================================================================
    type(spmat) :: PEV

contains


    !##########################################################################
    ! FUNCTION utility
    !
    ! CRRA utility of the consumption-leisure bundle.
    !##########################################################################
    function utility(cons, lab) result(u)

        implicit none
        real*8, intent(in) :: cons, lab
        real*8 :: u, flow, c_help, l_help

        c_help = max(cons, 1d-10)
        l_help = min(max(lab, 0d0), 1d0-1d-10)
        flow   = c_help**nu * (1d0-l_help)**(1d0-nu)

        if(abs(egam) < 1d-12)then
            u = log(flow)
        else
            u = flow**egam/egam
        endif

    end function utility


    !##########################################################################
    ! FUNCTION bellman_max
    !
    ! Bellman right-hand side for goldenx maximisation.
    !   u(c(aplus), l(aplus)) + beta*psi(j+1) * EV_{j+1}(aplus)
    !
    ! ielem is the element index (1-based) passed by goldenx.
    !   ia = ielem - 1   gives the asset-grid index.
    !
    ! Module communication variables:
    !   ij_com, ip_com, is_com  — current state
    !   gold_coefEV(0:NA)       — next-period EV coefficients
    !##########################################################################
    function bellman_max(ap, ielem) result(val)

        implicit none
        real*8,  intent(in) :: ap
        integer, intent(in) :: ielem
        real*8  :: val, available, wage_val, cons, lab, cont
        integer :: ia

        ia = ielem - 1
        available = cash_arr(ij_com, ia, ip_com, is_com)
        wage_val  = wage_arr(ij_com, ip_com, is_com)

        ! labour supply
        if(ij_com < JR)then
            lab = nu + (1d0-nu)*(ap - available)/wage_val
            lab = min(max(lab, 0d0), 1d0-1d-10)
        else
            lab = 0d0
        endif

        ! consumption
        cons = max(available + wage_val*lab - ap, 1d-10)

        ! continuation value via linear interpolation (hot path — direct call)
        cont = linint_Gen(ap, a, gold_coefEV)

        ! Bellman RHS (maximised by goldenx)
        val = utility(cons, lab) + beta*psi(ij_com+1)*cont

    end function bellman_max

end module collocation_globals


!##############################################################################
!##############################################################################
! MAIN PROGRAM
!##############################################################################
!##############################################################################
program LaborSupplyCollocation

    use collocation_globals

    implicit none

    ! initialise everything
    call initialize()

    ! start the clock
    call tic()

    ! solve household problem via collocation
    call solve_collocation()

    ! calculate distribution over state space (Mongey)
    call get_distribution_mongey()

    ! aggregate individual decisions
    call aggregation_mongey()

    ! stop the clock
    call toc()

    ! write output
    call output()

    ! make plots (policy functions + cohort aggregates)
    call make_plots()

    ! close files
    close(21)

contains


    !##########################################################################
    ! SUBROUTINE initialize
    !
    ! Sets up grids, stochastic processes, precomputed arrays, and the PEV
    ! sparse operator.
    !##########################################################################
    subroutine initialize()

        implicit none
        integer :: ij, ia, ip, is

        ! ---- survival probabilities ----
        psi = (/1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, &
                0.99906d0, 0.99908d0, 0.99906d0, 0.99907d0, 0.99901d0, &
                0.99899d0, 0.99896d0, 0.99893d0, 0.99890d0, 0.99887d0, &
                0.99886d0, 0.99878d0, 0.99871d0, 0.99862d0, 0.99853d0, &
                0.99841d0, 0.99835d0, 0.99819d0, 0.99801d0, 0.99785d0, &
                0.99757d0, 0.99735d0, 0.99701d0, 0.99676d0, 0.99650d0, &
                0.99614d0, 0.99581d0, 0.99555d0, 0.99503d0, 0.99471d0, &
                0.99435d0, 0.99393d0, 0.99343d0, 0.99294d0, 0.99237d0, &
                0.99190d0, 0.99137d0, 0.99085d0, 0.99000d0, 0.98871d0, &
                0.98871d0, 0.98721d0, 0.98612d0, 0.98462d0, 0.98376d0, &
                0.98226d0, 0.98062d0, 0.97908d0, 0.97682d0, 0.97514d0, &
                0.97250d0, 0.96925d0, 0.96710d0, 0.96330d0, 0.95965d0, &
                0.95619d0, 0.95115d0, 0.94677d0, 0.93987d0, 0.93445d0, &
                0.92717d0, 0.91872d0, 0.91006d0, 0.90036d0, 0.88744d0, &
                0.87539d0, 0.85936d0, 0.84996d0, 0.82889d0, 0.81469d0, &
                0.79705d0, 0.78081d0, 0.76174d0, 0.74195d0, 0.72155d0, &
                0.00000d0/)

        ! ---- age-efficiency profile ----
        eff(1:JR-1) = &
            (/1.0000d0, 1.0719d0, 1.1438d0, 1.2158d0, 1.2842d0, 1.3527d0, &
              1.4212d0, 1.4897d0, 1.5582d0, 1.6267d0, 1.6952d0, 1.7217d0, &
              1.7438d0, 1.7748d0, 1.8014d0, 1.8279d0, 1.8545d0, 1.8810d0, &
              1.9075d0, 1.9341d0, 1.9606d0, 1.9623d0, 1.9640d0, 1.9658d0, &
              1.9675d0, 1.9692d0, 1.9709d0, 1.9726d0, 1.9743d0, 1.9760d0, &
              1.9777d0, 1.9700d0, 1.9623d0, 1.9546d0, 1.9469d0, 1.9392d0, &
              1.9315d0, 1.9238d0, 1.9161d0, 1.9084d0, 1.9007d0, 1.8354d0, &
              1.7701d0, 1.7048d0/)
        eff(JR:JJ) = 0d0

        ! ---- pension ----
        pen = 0d0
        pen(JR:JJ) = 0.5d0*sum(eff)/dble(JR-1)*0.33d0

        ! ---- permanent productivity ----
        dist_theta = 1d0/dble(NP)
        theta(1)   = exp(-sqrt(sigma_theta))
        theta(2)   = exp( sqrt(sigma_theta))

        ! ---- persistent shock process (Rouwenhorst) ----
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi_mat)
        eta = exp(eta)

        ! ---- asset grid (growing) ----
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! ---- precompute cash-on-hand and wage ----
        do ij = 1, JJ
            do ip = 1, NP
                do is = 1, NS
                    wage_arr(ij, ip, is) = w_val*eff(ij)*theta(ip)*eta(is)
                    do ia = 0, NA
                        cash_arr(ij, ia, ip, is) = (1d0+r_val)*a(ia) + pen(ij)
                    enddo
                enddo
            enddo
        enddo

        ! ---- define spline function space on asset grid (CompEcon style) ----
        fspace_a = fundef(a, 1)    ! linear B-spline, breaks = a(0:NA)

        ! ---- construct PEV = kron(pi, I_nx) (sparse) ----
        PEV = sp_kron(sparse(pi_mat), sp_eye(nx))

        ! ---- open output file ----
        open(21, file='output_collocation.out')

    end subroutine initialize


    !##########################################################################
    ! SUBROUTINE solve_collocation
    !
    ! Backward-induction collocation solver.
    !
    ! For each age (JJ down to 1) and each productivity type ip:
    !   1. Compute Bellman RHS at every (ia, is) via golden-section search.
    !   2. Set V = rhsV, EV = PEV * V  (since Phi = I at collocation nodes).
    !##########################################################################
    subroutine solve_collocation()

        implicit none
        integer  :: ij, ia, ip, is
        real*8   :: ap_low_vec(nx), ap_high_vec(nx)
        real*8   :: ap_opt_vec(nx), fopt_vec(nx)
        real*8   :: Vvec(ndim), EVvec(ndim)
        real*8   :: rhsV(0:NA, NS)
        real*8   :: c_pol(0:NA, NS), l_pol(0:NA, NS), ap_pol(0:NA, NS)
        real*8   :: g1(ndim), g2(ndim)
        integer  :: max_ip, max_is

        ! =====================================================================
        !  Final period  (j = JJ)
        ! =====================================================================
        ij = JJ
        max_ip = 1          ! retirement: only one "effective" ip
        max_is = 1
        do ip = 1, max_ip
            do is = 1, max_is
                do ia = 0, NA
                    aplus(ij, ia, ip, is) = 0d0
                    l(ij, ia, ip, is)     = 0d0
                    c(ij, ia, ip, is)     = max(cash_arr(ij, ia, ip, is), 1d-10)
                    V(ij, ia, ip, is)     = utility(c(ij, ia, ip, is), 0d0)
                enddo
            enddo
            ! broadcast to all shock states
            do is = 1, NS
                aplus(ij, :, ip, is) = 0d0
                c(ij, :, ip, is)     = c(ij, :, 1, 1)
                l(ij, :, ip, is)     = 0d0
                V(ij, :, ip, is)     = V(ij, :, 1, 1)
            enddo
        enddo

        ! broadcast to all productivity types (retirement: same for all ip)
        do ip = 2, NP
            do is = 1, NS
                aplus(ij, :, ip, is) = aplus(ij, :, 1, is)
                c(ij, :, ip, is)     = c(ij, :, 1, is)
                l(ij, :, ip, is)     = l(ij, :, 1, is)
                V(ij, :, ip, is)     = V(ij, :, 1, is)
            enddo
        enddo

        ! EV  via PEV operator
        do ip = 1, NP
            Vvec = reshape(V(ij, :, ip, :), (/ndim/))
            EVvec = sp_matmul(PEV, Vvec)
            EV(ij, :, ip, :) = reshape(EVvec, (/nx, NS/))
        enddo

        residual(ij)  = 0d0
        ev_error(ij)  = 0d0

        ! =====================================================================
        !  Backward induction  (j = JJ-1 down to 1)
        ! =====================================================================
        do ij = JJ-1, 1, -1

            if(ij >= JR)then
                max_ip = 1
                max_is = 1
            else
                max_ip = NP
                max_is = NS
            endif

            do ip = 1, max_ip

                ! ---- compute Bellman RHS via golden-section search ----
                rhsV  = 0d0
                c_pol = 0d0
                l_pol = 0d0
                ap_pol = 0d0

                ij_com = ij
                ip_com = ip

                do is = 1, max_is
                    is_com = is

                    ! next-period EV coefficients for this (ip, is)
                    gold_coefEV(:) = EV(ij+1, :, ip, is)

                    ! build bounds arrays for all asset nodes at once
                    do ia = 0, NA
                        ap_low_vec(ia+1) = 0d0
                        ap_high_vec(ia+1) = cash_arr(ij, ia, ip, is) &
                            + wage_arr(ij, ip, is)*(1d0-1d-10) - 1d-10
                        ap_high_vec(ia+1) = min(a_u, max(0d0, ap_high_vec(ia+1)))
                    enddo

                    ! vectorised golden-section maximisation
                    call goldenx(bellman_max, ap_low_vec, ap_high_vec, &
                                 ap_opt_vec, fopt_vec)

                    ! unpack results
                    ap_pol(:, is) = ap_opt_vec
                    rhsV(:, is)   = fopt_vec

                    ! recover policies at the optima
                    do ia = 0, NA
                        call recover_policies(ap_opt_vec(ia+1), ij, ia, ip, is, &
                                              c_pol(ia, is), l_pol(ia, is))
                    enddo
                enddo

                ! ---- Newton step (trivial because Phi = I) ----
                ! c = rhsV,  ce = PEV * rhsV
                Vvec = reshape(rhsV, (/ndim/))
                EVvec = sp_matmul(PEV, Vvec)

                ! ---- store ----
                V(ij, :, ip, 1:max_is)  = rhsV(:, 1:max_is)
                EV(ij, :, ip, 1:max_is) = reshape(EVvec, (/nx, max_is/))
                c(ij, :, ip, 1:max_is)  = c_pol(:, 1:max_is)
                l(ij, :, ip, 1:max_is)  = l_pol(:, 1:max_is)
                aplus(ij, :, ip, 1:max_is) = ap_pol(:, 1:max_is)

                ! ---- residual check ----
                g1 = Vvec - reshape(rhsV, (/ndim/))
                g2 = EVvec - sp_matmul(PEV, Vvec)
                residual(ij) = max(residual(ij), maxval(abs(g1)), maxval(abs(g2)))

                ! ---- EV error (eqn 2 check) ----
                ev_error(ij) = max(ev_error(ij), maxval(abs(g2)))

                ! ---- fill retirement states ----
                if(ij >= JR)then
                    do is = 1, NS
                        aplus(ij, :, 1, is) = ap_pol(:, 1)
                        c(ij, :, 1, is)     = c_pol(:, 1)
                        l(ij, :, 1, is)     = l_pol(:, 1)
                        V(ij, :, 1, is)     = rhsV(:, 1)
                        EV(ij, :, 1, is)    = EV(ij, :, 1, 1)
                    enddo
                endif
            enddo

            ! ---- broadcast retirement solution to all productivity types ----
            if(ij >= JR)then
                do ip = 2, NP
                    do is = 1, NS
                        aplus(ij, :, ip, is) = aplus(ij, :, 1, is)
                        c(ij, :, ip, is)     = c(ij, :, 1, is)
                        l(ij, :, ip, is)     = l(ij, :, 1, is)
                        V(ij, :, ip, is)     = V(ij, :, 1, is)
                        EV(ij, :, ip, is)    = EV(ij, :, 1, is)
                    enddo
                enddo
            endif

            write(*,'(a,i3,a)')'Age: ', ij, ' DONE!'
        enddo

    end subroutine solve_collocation


    !##########################################################################
    ! SUBROUTINE recover_policies
    !
    ! Computes consumption and labour from a given aplus.
    !##########################################################################
    subroutine recover_policies(ap, ij, ia, ip, is, cons, lab)

        implicit none
        real*8,  intent(in)  :: ap
        integer, intent(in)  :: ij, ia, ip, is
        real*8,  intent(out) :: cons, lab
        real*8 :: available, wage_val

        available = cash_arr(ij, ia, ip, is)
        wage_val  = wage_arr(ij, ip, is)

        if(ij < JR)then
            lab = nu + (1d0-nu)*(ap - available)/wage_val
            lab = min(max(lab, 0d0), 1d0-1d-10)
        else
            lab = 0d0
        endif

        cons = max(available + wage_val*lab - ap, 1d-10)

    end subroutine recover_policies


    !##########################################################################
    ! SUBROUTINE get_distribution_mongey
    !
    ! Forward distribution using Mongey's Q-matrix approach.
    !   L_{j+1} = Q_j' * L_j
    ! where Q_j encodes both asset-interpolation weights (QX) and the
    ! exogenous shock transition (pi).
    !##########################################################################
    subroutine get_distribution_mongey()

        implicit none
        integer      :: ij, ip, is, ia
        real*8       :: L_vec(ndim), L_next(ndim), ap_vec(ndim)
        type(spmat)  :: QX_sp, Q_sp
        integer      :: idx

        ! ---- zero out distribution ----
        phi_dist = 0d0
        mass     = 0d0

        ! ---- initial distribution (age 1) ----
        ! all mass at a(0) in the middle eta state
        do ip = 1, NP
            L_vec = 0d0
            L_vec((is_initial-1)*nxd + 1) = dist_theta(ip)
            phi_dist(:, ip, :, 1) = reshape(L_vec, (/nxd, NS/))
        enddo

        mass(1, :) = sum(sum(phi_dist(:, :, :, 1), dim=1), dim=2)
        mass_error(1) = abs(sum(mass(1, :)) - 1d0)

        ! ---- iterate forward ----
        do ij = 1, JJ-1
            do ip = 1, NP

                ! extract aplus for age ij, productivity ip
                do is = 1, NS
                    do ia = 0, NA
                        idx = (is-1)*nxd + ia + 1
                        ap_vec(idx) = aplus(ij, ia, ip, is)
                    enddo
                enddo

                ! clamp to grid bounds
                ap_vec = min(max(ap_vec, a(0)), a(NA))

                ! sparse linear-spline basis at aplus points (CSR, CompEcon style)
                QX_sp = funbas(fspace_a, ap_vec)

                ! build Mongey transition matrix  Q = dprod(QX, pi)
                call sp_build_mongey_Q(QX_sp, pi_mat, Q_sp)
                call sp_free(QX_sp)

                ! current distribution vector (Mongey stacking)
                L_vec = reshape(phi_dist(:, ip, :, ij), (/ndim/))

                ! advance:  L_{j+1} = Q' * L_j
                L_next = sp_transpose_matmul(Q_sp, L_vec)

                ! store
                phi_dist(:, ip, :, ij+1) = reshape(L_next, (/nxd, NS/))

                ! free sparse matrix
                call sp_free(Q_sp)
            enddo

            ! mass check
            mass(ij+1, :) = sum(sum(phi_dist(:, :, :, ij+1), dim=1), dim=2)
            mass_error(ij+1) = abs(sum(mass(ij+1, :)) - 1d0)
        enddo

    end subroutine get_distribution_mongey


    !##########################################################################
    ! SUBROUTINE aggregation_mongey
    !
    ! Aggregates cohort statistics from the distribution and policy functions.
    !##########################################################################
    subroutine aggregation_mongey()

        implicit none
        integer :: ij, ia, ip, is
        real*8  :: hmat(0:NA, NS)
        real*8  :: csq, ysq, lsq, hsq, hls, denom, Lmass
        real*8  :: cval, lval, hval, yval

        c_coh  = 0d0; y_coh  = 0d0; l_coh  = 0d0; h_coh  = 0d0
        a_coh  = 0d0; v_coh  = 0d0
        cv_c   = 0d0; cv_y   = 0d0; cv_l   = 0d0; cv_h   = 0d0
        corr_hl = 0d0; frac_bor = 0d0

        do ij = 1, JJ

            csq = 0d0; ysq = 0d0; lsq = 0d0; hsq = 0d0; hls = 0d0

            ! precompute productivity for this age
            do is = 1, NS
                do ia = 0, NA
                    hmat(ia, is) = eff(ij)*eta(is)
                enddo
            enddo

            do ip = 1, NP
                do is = 1, NS
                    do ia = 0, NA
                        Lmass = phi_dist(ia, ip, is, ij)
                        if(Lmass < 1d-16) cycle

                        cval = c(ij, ia, ip, is)
                        lval = l(ij, ia, ip, is)
                        hval = hmat(ia, is)*theta(ip)
                        yval = hval*lval

                        ! cohort means
                        c_coh(ij) = c_coh(ij) + cval*Lmass
                        y_coh(ij) = y_coh(ij) + yval*Lmass
                        l_coh(ij) = l_coh(ij) + lval*Lmass
                        h_coh(ij) = h_coh(ij) + hval*Lmass
                        a_coh(ij) = a_coh(ij) + a(ia)*Lmass
                        v_coh(ij) = v_coh(ij) + V(ij, ia, ip, is)*Lmass

                        ! second moments
                        csq = csq + cval**2*Lmass
                        ysq = ysq + yval**2*Lmass
                        lsq = lsq + lval**2*Lmass
                        hsq = hsq + hval**2*Lmass
                        hls = hls + hval*lval*Lmass

                        ! borrowing constrained
                        if(aplus(ij, ia, ip, is) < 1d-6) &
                            frac_bor(ij) = frac_bor(ij) + Lmass
                    enddo
                enddo
            enddo

            ! pension and income
            pen_coh(ij)    = pen(ij)*sum(mass(ij, :))
            income_coh(ij) = y_coh(ij) + r_val*a_coh(ij)

            ! coefficients of variation
            cv_c(ij) = sqrt(max(csq - c_coh(ij)**2, 0d0))/max(c_coh(ij), 1d-10)
            cv_y(ij) = sqrt(max(ysq - y_coh(ij)**2, 0d0))/max(y_coh(ij), 1d-10)
            cv_l(ij) = sqrt(max(lsq - l_coh(ij)**2, 0d0))/max(l_coh(ij), 1d-10)
            cv_h(ij) = sqrt(max(hsq - h_coh(ij)**2, 0d0))/max(h_coh(ij), 1d-10)

            ! correlation hours-productivity
            denom = sqrt(max(hsq - h_coh(ij)**2, 0d0)) &
                  * sqrt(max(lsq - l_coh(ij)**2, 0d0))
            corr_hl(ij) = (hls - h_coh(ij)*l_coh(ij))/max(denom, 1d-10)
        enddo

        frac_bor(JJ) = 1d0

    end subroutine aggregation_mongey


    !##########################################################################
    ! SUBROUTINE output
    !
    ! Writes cohort aggregates and diagnostics to output_collocation.out.
    !##########################################################################
    subroutine output()

        implicit none
        integer :: ij, ia, ip, is, iamax(JJ), ages(JJ)

        ages = 20 + (/(ij, ij=1,JJ)/)

        ! ---- check maximum grid point used ----
        iamax = 0
        do ij = 1, JJ
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        if(phi_dist(ia, ip, is, ij) > 1d-8) iamax(ij) = ia
                    enddo
                enddo
            enddo
        enddo

        ! ---- policy / Newton diagnostics ----
        write(21, '(a)') ' AGE    RESIDUAL     EV_ERROR      MIN_C       MIN_L       MAX_L      MAX_APLUS'
        do ij = 1, JJ
            write(21, '(i4, 2es13.4, 2f12.4, f12.4, f13.4)') &
                ages(ij), residual(ij), ev_error(ij), &
                minval(c(ij, :, :, :)), minval(l(ij, :, :, :)), &
                maxval(l(ij, :, :, :)), maxval(aplus(ij, :, :, :))
        enddo
        write(21, '(a/)') '--------------------------------------------------------------------'

        ! ---- cohort aggregates ----
        write(21, '(a,a)') ' IJ      CONS     HOURS  EARNINGS    INCOME      PENS    ASSETS', &
            '     CV(C)     CV(L)     CV(Y)     VALUE     IAMAX'
        do ij = 1, JJ
            write(21, '(i3,10f10.3,i10)') ij, c_coh(ij), l_coh(ij), y_coh(ij), &
                income_coh(ij), pen_coh(ij), a_coh(ij), &
                cv_c(ij), cv_l(ij), cv_y(ij), v_coh(ij), iamax(ij)
        enddo
        write(21, '(a/)') '--------------------------------------------------------------------'

        ! ---- distribution diagnostic ----
        write(21, '(a)') ' AGE     MASS_ERROR    FRAC_BOR'
        do ij = 1, JJ
            write(21, '(i4, es14.4, f12.4)') ages(ij), mass_error(ij), frac_bor(ij)
        enddo

    end subroutine output


    !##########################################################################
    ! SUBROUTINE make_plots
    !
    ! Plots policy functions and cohort aggregates, matched to the MATLAB
    ! make_policy_plots / make_distribution_plots in labor_supply.m.
    !##########################################################################
    subroutine make_plots()

        implicit none
        integer  :: ages(JJ), pick(5), ip, is, k, j
        real*8   :: xplot(0:NA)
        character(len=20) :: leg

        ages = 20 + (/(j, j=1,JJ)/)
        ip = 1
        is = is_initial

        ! representative ages: 1, 20, JR, 60, JJ
        pick = (/1, 20, JR, 60, JJ/)
        xplot = a

        ! ---- Savings policy (aplus) ----
        do k = 1, size(pick)
            j = pick(k)
            write(leg, '("Age ", i3)') ages(j)
            call plot(a(0:NA), aplus(j, 0:NA, ip, is), legend=trim(leg))
        enddo
        call execplot(xlabel='Assets a', ylabel='Next-period assets a''')

        ! ---- Consumption policy ----
        do k = 1, size(pick)
            j = pick(k)
            write(leg, '("Age ", i3)') ages(j)
            call plot(a(0:NA), c(j, 0:NA, ip, is), legend=trim(leg))
        enddo
        call execplot(xlabel='Assets a', ylabel='Consumption c')

        ! ---- Labour policy ----
        do k = 1, size(pick)
            j = pick(k)
            write(leg, '("Age ", i3)') ages(j)
            call plot(a(0:NA), l(j, 0:NA, ip, is), legend=trim(leg))
        enddo
        call execplot(xlabel='Assets a', ylabel='Hours l')

        ! ---- Value function ----
        do k = 1, size(pick)
            j = pick(k)
            write(leg, '("Age ", i3)') ages(j)
            call plot(a(0:NA), V(j, 0:NA, ip, is), legend=trim(leg))
        enddo
        call execplot(xlabel='Assets a', ylabel='Value function V')

        ! ---- Cohort means: consumption & earnings+pension ----
        call plot(dble(ages), c_coh, legend='Consumption')
        call plot(dble(ages), y_coh + pen_coh, legend='Earnings + pension')
        call execplot(xlabel='Age j', ylabel='Mean')

        ! ---- Cohort labour ----
        call plot(dble(ages), l_coh)
        call execplot(xlabel='Age j', ylabel='Hours')

        ! ---- Cohort assets ----
        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

        ! ---- Coefficients of variation ----
        call plot(dble(ages), cv_c, legend='Consumption')
        call plot(dble(ages), cv_y, legend='Earnings')
        call plot(dble(ages), cv_l, legend='Hours')
        call execplot(xlabel='Age j', ylabel='Coefficient of Variation')

        ! ---- Variance decomposition of earnings ----
        call plot(dble(ages), cv_y, legend='Earnings')
        call plot(dble(ages), cv_l, legend='Hours')
        call plot(dble(ages), cv_h, legend='Productivity')
        call plot(dble(ages), corr_hl, legend='Corr(H, l)')
        call execplot(xlabel='Age j', ylabel='Variance Decomposition')

        ! ---- Fraction borrowing constrained ----
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1))
        call execplot(xlabel='Age j', ylabel='Frac. Borrowing Constrained')

    end subroutine make_plots

end program LaborSupplyCollocation
