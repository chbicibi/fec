! module mod_prob1
!     use interface
!     use util
!     use individual
!     use problem
!     use mod_problem
!     implicit none

!     private
!     public :: TProb1

!     type, extends(TProblem) :: TProb1
!         integer :: nvar, nobj, ncon
!         real(8), allocatable :: ubounds(:), lbounds(:)

!         contains

!         generic :: initialize => initialize_prob

!         procedure :: initialize_prob => initialize_prob1
!         procedure :: call_indiv => call_prob1
!         procedure :: scaling
!     end type TProb1

!     contains

!     ! ============================================================================
!     ! Problem Class
!     ! ============================================================================

!     subroutine initialize_prob1(this, nvar, nobj, ncon, ubounds, lbounds)
!         implicit none
!         class(TProb1), intent(inout) :: this
!         integer, intent(in) :: nvar, nobj, ncon
!         real(8), intent(in) :: ubounds(:), lbounds(:)

!         this%nvar = nvar
!         this%nobj = nobj
!         this%ncon = ncon
!         this%ubounds = ubounds
!         this%lbounds = lbounds
!     end subroutine initialize_prob1


!     ! ############################################################################
!     ! 目的関数呼び出しサブルーチン
!     ! ############################################################################
!     subroutine call_prob1(this, indiv)
!         implicit none
!         class(TProb1), intent(inout) :: this
!         class(TIndiv), intent(inout) :: indiv

!         call TNK(indiv%dvariables, indiv%objectives, indiv%constraints)

!         ! call obj_fun(this%scaling(indiv%dvariables), this%nobj, indiv%objectives, this%ncon, indiv%constraints)
!         indiv%feasible = all(indiv%constraints > 0)
!     end subroutine call_prob1

!     function scaling(this, variables) result(res)
!         implicit none
!         class(TProb1), intent(in) :: this
!         real(8), intent(in) :: variables(:)
!         real(8), allocatable :: res(:)

!         res = this%lbounds + (this%ubounds - this%lbounds) * variables
!     end function scaling
! end module mod_prob1


include "util.f08"
include "interface.f08"


module mod_nsga2
    implicit none

    type :: TIndividual
        integer :: idx
        double precision, allocatable :: x(:)
        double precision, allocatable :: f(:)
        double precision, allocatable :: g(:)
        double precision, allocatable :: fitness(:)
        double precision, allocatable :: crowding(:)
        integer :: rank
        logical :: is_eval = .false.
    end type TIndividual

    type :: TPopulation
        integer :: i_gen
        type(TIndividual), allocatable :: ind(:)
    end type TPopulation

contains
end module mod_nsga2


module mod_problem
    implicit none
    contains

    subroutine problem_OSY(x, f, g)
        double precision, intent(in) :: x(:)
        double precision, intent(out), allocatable :: f(:)
        double precision, intent(out), allocatable :: g(:)
    ! print *, 'problem_OSY', size(x)
    ! print *, 'problem_OSY', allocated(x)

    f = [-sum([25d0, 1d0, 1d0, 1d0, 1d0] * (x(:5) - [2d0, 2d0, 1d0, 4d0, 1d0]) ** 2), &
         sum(x ** 2)]

    g = [-2d0 + x(1) + x(2),              &
          6d0 - x(1) - x(2),              &
          2d0 + x(1) - x(2),              &
          2d0 - x(1) + 3d0 * x(2),        &
          4d0 - (x(3) - 3d0) ** 2 - x(4), &
         -4d0 + (x(5) - 3d0) ** 2 + x(6)]
         print *, f
    end subroutine problem_OSY
end module mod_problem


program main
    use util
    use mod_nsga2
    use mod_problem
    implicit none

    ! type(TNSGA2C) :: optimizer
    ! type(TProb1) :: problem

    integer :: n_x, n_f, n_g, n_ind, i_gen, n_gen
    character(32) :: tf, ts, tc, tm
    real(8), allocatable :: lbounds(:), ubounds(:)
    integer :: unit, start_time, i, accumulator

    type(TPopulation), allocatable :: pop(:)

    ! 前処理
    call set_randomseed


    ! パラメータ読み込み
    ! open(newunit=unit, file="inputopt.txt", status="old")
    !     read(unit, *) n_x        ! ... number of design variables
    !     read(unit, *) n_f        ! ... number of objective functions
    !     read(unit, *) n_g        ! ... number of constraint functions
    !     read(unit, *) n_ind      ! ... population size
    !     read(unit, *)             ! ... neighborhood size
    !     read(unit, *) n_gen      ! ... ngen: total geration number
    !     read(unit, *) tf          ! ... test function type
    !     read(unit, *) ts          ! ... selection method type
    !     read(unit, *) tc          ! ... crossover method type
    !     read(unit, *) tm          ! ... mutation method type
    ! close(unit)
    n_gen = 1
    n_ind = 10
    n_x = 6
    n_f = 2
    n_g = 6


    ! 設計範囲設定
    ! allocate(lbounds(nvar))
    ! allocate(ubounds(nvar))
    ! open(newunit=unit, file="bounds.txt", status="old")
    !     do i = 1, nvar
    !         read(unit, *) lbounds(i), ubounds(i)
    !     end do
    ! close(unit)


    ! 初期化
    ! call optimizer%initialize(nx=nvar, m=nobj, N=npop, f=tf, selection=ts, crossover=tc, mutation=tm)
    ! optimizer%elite_preservation = .true.
    ! optimizer%sharing = .true.
    ! optimizer%dup_rejection = .false.
    ! optimizer%history_preservation = .true.
    ! optimizer%show_log = .true.

    ! 初期集団生成
    allocate(pop(n_gen))
    allocate(pop(1)%ind(n_ind))
    do i = 1, n_ind
        pop(1)%ind(i)%idx = i
        allocate(pop(1)%ind(i)%x(6), source=random_array(n_x))
        ! pop(1)%ind(i)%x = random_array(n_x)
        ! print *, pop(1)%ind(i)%x
    end do
    accumulator = n_ind

    ! 評価
    call evaluate(pop(1)%ind, problem_OSY)

    do i = 1, size(pop(1)%ind)
        print *, pop(1)%ind(i)%f
    end do


        ! 最適化ループ

        ! do i_gen = 1, n_gen
        !         eval(pop(i_gen))

    ! 用意されているテスト関数を使用する場合は次の2行をコメントアウトしてください
    ! call problem%initialize(nvar, nobj, ncon, ubounds, lbounds) ! 問題オブジェクト初期化
    ! call optimizer%set_problem(problem)                         ! 問題オブジェクト設定

    ! 計算
    ! call system_clock(start_time)      ! 開始時間記録
    ! call optimizer%prepare_calculation ! 初期個体生成
    ! call optimizer%run(ngen)           ! 最適化
    ! call elapsed_time(start_time)      ! 実行時間表示


    ! 計算結果出力
    ! call optimizer%save_result("result_all.csv", elite="all", feasible="all")          ! 最終世代の全個体
    ! call optimizer%save_result("result_elite.csv", elite="only", feasible="all")       ! 最終世代の非劣解個体
    ! call optimizer%save_result("result_feasible.csv", elite="all", feasible="only")    ! 最終世代の全個体 (実行可能解のみ)
    ! call optimizer%save_history("history_all.csv", elite="all", feasible="all")        ! 全世代の全個体
    ! call optimizer%save_history("history_elite.csv", elite="only", feasible="all")     ! 全世代の非劣解個体
    ! call optimizer%save_history("history_feasible.csv", elite="only", feasible="only") ! 全世代の非劣解個体 (実行可能解のみ)
    ! call optimizer%save_history("history_infeasible.csv", elite="all", feasible="not") ! 全世代の全個体 (制約違反解のみ)

    contains

    subroutine evaluate(ind, proc)
        type(TIndividual), intent(inout) :: ind(:)
        integer :: i
        interface
            subroutine proc(x, f, g)
                double precision, intent(in) :: x(:)
                double precision, intent(out), allocatable :: f(:)
                double precision, intent(out), allocatable :: g(:)
            end subroutine
        end interface

        ! procedure :: proc
        print *, size(ind)

        do i = 1, size(ind)
            if (ind(i)%is_eval) continue

            print *, size(ind(i)%x)

            call proc(ind(i)%x, ind(i)%f, ind(i)%g)
            ind(i)%is_eval = .true.
        end do
    end subroutine evaluate

    subroutine calc_rank(ind)
        obj = reshape([(pop%ind(i)%f(n_f), i = 1, size(pop))], &
                      [n_f, n_ind])

        pop%rank = rank_pareto(obj)

        if (.true.) then

          ! pop%crowding = crowding_distance(objectives, pop%rank)
          ! pop%crowding = min(crowding_distance(objectives, pop%rank), 1d0)
          pop%crowding = tanh(crowding_distance(objectives, pop%rank))
          pop%fitness = a * (1.0d0 - a) ** (pop%rank - pop%crowding)
          ! pop%fitness = 1.0d0 / (pop%rank + 1 - pop%crowding)
        else
          pop%fitness = a * (1.0d0 - a) ** (pop%rank - 1)
          ! pop%fitness = 1.0d0 / pop%rank
        end if
    end subroutine calc_rank

    function rank_pareto_con2_proc(val, con, mask, proc) result(rank)
        real(8), intent(in) :: val(:, :), con(:, :)
        logical, intent(in) :: mask(:)
        procedure(func_2d_2d_1l_l) :: proc
        integer, allocatable :: rank(:)
        integer :: n, i, j

        n = size(val, dim=2)
        rank = rank_pareto_core([((proc(val(:, [i, j]), con(:, [i, j]), mask([i, j])), &
                                   j = 1, n), i = 1, n)], n)
    end function rank_pareto_con2_proc

    function rank_pareto_core(dominated, n) result(rank)
        logical, intent(in) :: dominated(n, n)
        integer, intent(in) :: n
        integer, allocatable :: rank(:)
        integer, allocatable :: num_dominated(:)
        logical, allocatable :: front_mask(:)
        integer :: rank_no, i

        rank = filled(n, 0)
        ! num_dominated = [(count(dominated(n*i+1:n*(i+1))), i = 0, n - 1)]
        allocate(num_dominated(n), source=[(count(dominated(:, i)), i = 1, n)])

        do rank_no = 1, n
            front_mask = rank == 0 .and. num_dominated == 0
            where (front_mask) rank = rank_no
            if (all(rank > 0)) exit
            do i = 1, n
                ! if (rank(i) == 0) num_dominated(i) = num_dominated(i) - &
                !                     count(front_mask .and. dominated(n*(i-1)+1:n*i))
                if (rank(i) == 0) num_dominated(i) = num_dominated(i) - &
                    count(front_mask .and. dominated(:, i))
            end do
        end do
    end function rank_pareto_core
end program main
