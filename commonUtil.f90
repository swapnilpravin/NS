!-------------------------------------------------------------------------
! Common Fortran functions and routines
!-------------------------------------------------------------------------
!
! MODULE: commonUtil
!
!> @author
!> Swapnil Pravin <swapnilpravin@gmail.com>
!
! DESCRIPTION:
!> A set of Fortran functions and routines for common utility tasks
!  related to file IO, array operations etc.
!
! REVISION HISTORY:
! 03 Feb 2014 - Added doxygen setup
!
!-------------------------------------------------------------------------

module commonUtil
    
    implicit none

contains

    !---------------------------------------------------------------------
    ! DESCRIPTION
    !> Returns the number of lines in a file
    !
    !> @param[in] input_file: ( character(len=30) ) file name 
    !> @param[out] N: ( integer ) number of lines
    !---------------------------------------------------------------------

    subroutine number_of_lines(input_file, N)
        implicit none
        ! dummy vars
        character (len=*), intent(in) :: input_file
        integer, intent(out) :: N
        
        ! local vars
        integer :: temp=0, IOstatus=0, i=0

        open(unit=10, file=input_file, status="old", action="read")
        do
            read(10,*,iostat=IOstatus) temp
            if (iostatus < 0) then
                exit
            else
                i=i+1
            endif
        enddo
            
        close(10)
        N = i
    end subroutine
    
    !---------------------------------------------------------------------
    ! DESCRIPTION
    !> Load data from file into a 1D array. Each line in the file should 
    !  contain one data element
    !
    !> @param[in] filename: ( character(len=*) ) file name 
    !> @param[out] number_of_lines: ( integer ) number of lines in file
    !> @param[out] output_array: ( double(:) ) 1D double array in which
    !                           data is to be filled.
    !---------------------------------------------------------------------

    subroutine load1DArrayFromFile(filename, output_array)

        ! Dummy vars
        character (len=*), intent(in) :: filename
        !integer, intent(in) :: number_of_lines
        real, dimension(:), allocatable, intent(inout) :: output_array
        
        ! local vars
        integer :: i
        integer :: N ! Number of lines in filename

        call number_of_lines(filename, N)
        allocate(output_array(N))

        open(unit=10, file=filename, status='old', action='read')
        do i=1,N
            read(10,*) output_array(i)
        enddo
        
        close(10)
    end subroutine
    
    !---------------------------------------------------------------------
    ! DESCRIPTION
    !> Write 1D double array to file, one element per line
    !
    !> @param[in] array: ( double(:) ) array to write
    !> @param[out] filename: ( character(len=*) ) file name 
    !---------------------------------------------------------------------

    subroutine write1DArrayToFile(array, filename)
         
        !Dummy vars
        character (len=*) :: filename
        double precision, dimension(:) :: array
        ! local vars
        integer :: N
        integer :: i
        
        open(unit=10, file=filename, status='replace', action='write')
        
        N = size(array)
        do i=1,N
            write(10,*) array(i)
        enddo
        close(10)
        print *, filename, " written to disk."
    end subroutine write1DArrayToFile

    subroutine load2DArrayFromFile(filename, output_array)

        ! Dummy vars
        character (len=*) :: filename
        real, dimension(:,:), allocatable, intent(out) :: output_array

        ! local vars
        integer :: i,j
        integer :: Nrows, Ncols


        call number_of_lines(filename, Nrows)

    end subroutine

    !---------------------------------------------------------------------
    ! DESCRIPTION
    !> Write 2D double array to file, one row per line
    !
    !> @param[in] array: ( double(:,:) ) array to write
    !> @param[out] filename: ( character(len=*) ) file name 
    !---------------------------------------------------------------------

    subroutine write2DArrayToFile(array, filename)
        
        !Dummy vars
        character (len=*) :: filename
        double precision, dimension(:,:) :: array
        ! local vars
        integer, dimension(1:2) :: N
        integer :: i, j
        
        open(unit=10, file=filename, status='replace', action='write')
        
        N = shape(array)
        do i=1,N(1)
            write(10,*) ( array(i,j), j = 1,N(2) )
        enddo
        close(10)
        print *, filename, " written to disk."
    end subroutine write2DArrayToFile
        
    !---------------------------------------------------------------------
    ! DESCRIPTION
    !> Create an array with n elements between xmin and xmax (including both)
    !
    !> @param[in] xmin: ( double ) first element of array
    !> @param[in] xmax: ( double ) last element of array    
    !> @param[in] n: ( integer ) number of elements in the array
    !> @param[out] x: ( double(n) ) array to fill the elements in
    !---------------------------------------------------------------------

    subroutine linspace(xmin, xmax, n, x)

        ! dummy vars
        real, intent(in) :: xmin
        real, intent(in) :: xmax
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: x
        
        ! local vars
        integer :: i
        real :: dx
        
        dx = (xmax - xmin) / (n-1)

        do i=1,n
            x(i) = xmin + dx * (i-1)
        enddo
    end subroutine linspace
	
	!---------------------------------------------------------------------
	! DESCRIPTION
	!> call this before generating random numbers for a matrix
	!
	!       REAL :: r(5,5)
    !       CALL init_random_seed()
    !       CALL RANDOM_NUMBER(r)
	!---------------------------------------------------------------------
	subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
    end subroutine init_random_seed

end module commonUtil
