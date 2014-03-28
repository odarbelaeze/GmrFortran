module modsample

    use modtypes
    use modhelpers

    implicit none

contains

    subroutine sample_create(r_, s_, species_, st_, stat_)

        double precision, allocatable, dimension(:,:), intent(inout) :: r_
        double precision, allocatable, dimension(:,:), intent(inout) :: s_
        integer, allocatable, dimension(:), intent(inout) :: species_
        type(sampleTraits), intent(in) :: st_

        integer, optional, intent(out) :: stat_
        integer :: err = 0

        integer :: i, j, k, n, id

        if (present(stat_)) stat_ = 0

        if (allocated(r_)) deallocate(r_, stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 1
            print *, "*** r_: Deallocation request denied"
        end if

        if (allocated(s_)) deallocate(s_, stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 2
            print *, "*** s_: Deallocation request denied"
        end if

        if (allocated(species_)) deallocate(species_, stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 3
            print *, "*** species_: Deallocation request denied"
        end if

        n = st_ % w * st_ % h * st_ % l + st_ % nElectrons

        allocate(r_(3, n), stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 4
            print *, "*** r_: Allocation request denied"
            return
        endif

        allocate(s_(3, n), stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 5
            print *, "*** s_: Allocation request denied"
            return
        endif

        allocate(species_(n), stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 6
            print *, "*** species_: Allocation request denied"
            return
        endif

        id = 1

        do i = 1, st_ % w, 1
            do j = 1, st_ % h, 1
                do k = 1, st_ % l, 1
                    r_(:, id) = (/ i, j, k /) + (/ 0.0d0, 0.0d0, 0.0d0 /)
                    call helpers_randomspin(s_(:, id))
                    species_(id) = SP_ION
                    id = id + 1
                end do
            end do
        end do

        call random_number(r_(:, id:))

        do i = 1, st_ % nElectrons, 1
            r_(:, id) = r_(:, id) * (/ st_ % w, st_ % h, st_ % l /)
            call helpers_randomspin(s_(:, id))
            species_(id) = SP_ELECTRON
            id = id + 1;
        end do

    end subroutine sample_create


    subroutine sample_destroy(r_, s_, species_, stat_)
        double precision, allocatable, dimension(:,:), intent(inout) :: r_
        double precision, allocatable, dimension(:,:), intent(inout) :: s_
        integer, allocatable, dimension(:), intent(inout) :: species_

        integer, optional :: stat_
        integer :: err = 0

        if (present(stat_)) stat_ = 0

        if (allocated(r_)) deallocate(r_, stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 1
            print *, "*** r_: Deallocation request denied"
        end if

        if (allocated(s_)) deallocate(s_, stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 2
            print *, "*** s_: Deallocation request denied"
        end if

        if (allocated(species_)) deallocate(species_, stat=err)
        if (err /= 0) then
            if (present(stat_)) stat_ = 2
            print *, "*** species_: Deallocation request denied"
        end if

    end subroutine sample_destroy

    subroutine sample_updatenbh(innb_, inbh_, ennb_, enbh_, &
                                r_, species_, st_, it_, stat_)

        integer, allocatable, dimension(:), intent(inout) :: innb_
        integer, allocatable, dimension(:,:), intent(inout) :: inbh_
        integer, allocatable, dimension(:), intent(inout) :: ennb_
        integer, allocatable, dimension(:,:), intent(inout) :: enbh_
        double precision, dimension(:,:), intent(in) :: r_
        integer, dimension(:), intent(in) :: species_
        type(sampleTraits), intent(in) :: st_
        type(interactionTraits), intent(in) :: it_
        integer, optional, intent(inout) :: stat_

        double precision, dimension(3) :: dims
        double precision :: distance
        integer :: inbhsize, enbhsize
        integer :: i, j, eid, iid
        integer :: err = 0

        if (present(stat_)) stat_ = 0

        if (.not. allocated(r_) .or. .not. allocated(species_)) then
            if (present(stat_)) stat_ = 1
            return
        end if

        if (.not. allocated(innb_)) then
            allocate(innb_(size(r_, dim=2)), stat=err)
            if (err /= 0) then
                print *, "***innb_: Allocation request denied"
                if (present(stat_)) stat_ = 2
                return
            end if
        end if

        if (.not. allocated(ennb_)) then
            allocate(ennb_(size(r_, dim=2)), stat=err)
            if (err /= 0) then
                print *, "***ennb_: Allocation request denied"
                if (present(stat_)) stat_ = 2
                return
            end if
        end if

        ennb_ = 0
        innb_ = 0
        dims = (/ st_ % w, st_ % l, st_ % h /)

        do i = 1, size(r_, dim=2), 1
            do j = 1, size(r_, dim=2), 1
                distance = helpers_distance(r_(:,i), r_(:,j), dims)
                if (species_(i) == SP_ELECTRON .and. & 
                    species_(j) == SP_ELECTRON) then
                    if (distance < it_ % eeCutOff) ennb_(i) = ennb_(i) + 1
                else if (species_(i) == SP_ION .and. & 
                         species_(j) == SP_ION) then
                    if (distance < it_ % iiCutOff) innb_(i) = innb_(i) + 1
                else
                    if (distance < it_ % eiCutOff) then
                        if (species_(j) == SP_ELECTRON) ennb_(i) = ennb_(i) + 1
                        if (species_(j) == SP_ION) innb_(i) = innb_(i) + 1
                    end if
                end if
            end do
        end do

        if (.not. allocated(inbh_) .or. &
            size(inbh_, dim=1) < maxval(innb_) .or. &
            size(inbh_, dim=2) /= size(innb_, dim=1)) then
            
            if (allocated(inbh_)) deallocate(inbh_, stat=err)
            if (err /= 0) then
                print *, "***inbh_: Deallocation request denied"
                if (present(stat_)) stat_ = 3
                return
            end if

            inbhsize = 2

            do while (inbhsize < maxval(innb_))
                inbhsize = inbhsize * 2
            end do

            allocate(inbh_(inbhsize, size(innb_, dim=1)), stat=err)
            if (err /= 0) then
                print *, "***inbh_: Allocation request denied"
                if (present(stat_)) stat_ = 3
                return
            end if
        end if

        if (.not. allocated(enbh_) .or. &
            size(enbh_, dim=1) < maxval(ennb_) .or. &
            size(enbh_, dim=2) /= size(ennb_, dim=1)) then
            
            if (allocated(enbh_)) deallocate(enbh_, stat=err)
            if (err /= 0) then
                print *, "***enbh_: Deallocation request denied"
                if (present(stat_)) stat_ = 4
                return
            end if

            enbhsize = 2

            do while (enbhsize < maxval(ennb_))
                enbhsize = enbhsize * 2
            end do

            allocate(enbh_(enbhsize, size(ennb_, dim=1)), stat=err)
            if (err /= 0) then
                print *, "***enbh_: Allocation request denied"
                if (present(stat_)) stat_ = 4
                return
            end if
        end if

        do i = 1, size(r_, dim=2), 1
            eid = 1
            iid = 1
            do j = 1, size(r_, dim=2), 1
                distance = helpers_distance(r_(:,i), r_(:,j), dims)
                if (species_(i) == SP_ELECTRON .and. & 
                    species_(j) == SP_ELECTRON) then
                    if (distance < it_ % eeCutOff) then
                        enbh_(eid, i) = j
                        eid = eid + 1
                    end if
                else if (species_(i) == SP_ION .and. & 
                         species_(j) == SP_ION) then
                    if (distance < it_ % iiCutOff) then
                        inbh_(iid, i) = j
                        iid = iid + 1
                    end if
                else
                    if (distance < it_ % eiCutOff) then
                        if (species_(j) == SP_ELECTRON) then
                            enbh_(eid, i) = j
                            eid = eid + 1
                        end if
                        if (species_(j) == SP_ION) then
                            inbh_(iid, i) = j
                            iid = iid + 1
                        end if
                    end if
                end if
            end do
        end do

    end subroutine sample_updatenbh


    subroutine sample_destoynbh(innb_, inbh_, ennb_, enbh_, stat_)

        integer, allocatable, dimension(:), intent(inout) :: innb_
        integer, allocatable, dimension(:,:), intent(inout) :: inbh_
        integer, allocatable, dimension(:), intent(inout) :: ennb_
        integer, allocatable, dimension(:,:), intent(inout) :: enbh_
        integer, optional, intent(out) :: stat_ = 0

        integer :: err = 0

        if (allocated(innb_)) deallocate(innb_, stat=err)
        if (err /= 0) then
            print *, "***innb_: Deallocation request denied"
            if (present(stat_)) stat_ = 1
        end if

        if (allocated(inbh_)) deallocate(inbh_, stat=err)
        if (err /= 0) then
            print *, "***inbh_: Deallocation request denied"
            if (present(stat_)) stat_ = 1
        end if

        if (allocated(ennb_)) deallocate(ennb_, stat=err)
        if (err /= 0) then
            print *, "***ennb_: Deallocation request denied"
            if (present(stat_)) stat_ = 1
        end if

        if (allocated(enbh_)) deallocate(enbh_, stat=err)
        if (err /= 0) then
            print *, "***enbh_: Deallocation request denied"
            if (present(stat_)) stat_ = 1
        end if
        
    end subroutine sample_destoynbh

end module modsample
