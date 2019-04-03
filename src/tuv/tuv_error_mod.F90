module tuv_error_mod
  implicit none

contains

  subroutine tuv_error_fatal(msg)
#ifdef CPRNAG
   ! NAG does not provide these as intrinsics, but it does provide modules
   ! that implement commonly used POSIX routines.
    use f90_unix_proc, only: abort
#endif

    character(len=*), intent(in) :: msg

    write(*,*) 'ERROR: '//trim(msg)
    call abort
    
  end subroutine tuv_error_fatal
end module tuv_error_mod
