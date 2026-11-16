module fvtkhtg_kinds_m
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
   implicit none
   private

   ! Fundamental integer and floating-point kind parameters
   integer, parameter, public :: id_k = int64
   integer, parameter, public :: i64_k = int64
   integer, parameter, public :: i32_k = int32
   integer, parameter, public :: byte_k = int8
   integer, parameter, public :: r64_k = real64
   integer, parameter, public :: r32_k = real32

   ! Status and error codes
   integer(i32_k), parameter, public :: HTG_OK = 0_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_INVALID_ARG = 1_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_OUT_OF_BOUNDS = 2_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_INVALID_STATE = 3_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_TOPOLOGY = 4_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_INDEXING = 5_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_MASK = 6_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_INTERFACE = 7_i32_k
   integer(i32_k), parameter, public :: HTG_ERR_NOT_FOUND = 8_i32_k

   ! Global node indexing mode constants
   integer(i32_k), parameter, public :: HTG_INDEX_UNSET = 0_i32_k
   integer(i32_k), parameter, public :: HTG_INDEX_IMPLICIT = 1_i32_k
   integer(i32_k), parameter, public :: HTG_INDEX_EXPLICIT = 2_i32_k

   ! Field data type constants
   integer(i32_k), parameter, public :: HTG_TYPE_REAL64 = 1_i32_k
   integer(i32_k), parameter, public :: HTG_TYPE_REAL32 = 2_i32_k
   integer(i32_k), parameter, public :: HTG_TYPE_INT64 = 3_i32_k
   integer(i32_k), parameter, public :: HTG_TYPE_INT32 = 4_i32_k
   integer(i32_k), parameter, public :: HTG_TYPE_BIT = 5_i32_k

   ! Sentinel constants
   integer(id_k), parameter, public :: HTG_INVALID_ID = -1_id_k
   integer(id_k), parameter, public :: HTG_LEAF_ELDER_SENTINEL = -1_id_k
   integer(id_k), parameter, public :: HTG_UINT32_MAX = 4294967295_id_k

   ! Diagnostic status container
   type, public :: htg_status_t
      integer(i32_k) :: code = HTG_OK
      character(len=256) :: message = ''
   contains
      procedure :: ok => status_ok
      procedure :: clear => status_clear
      procedure :: set => status_set
   end type htg_status_t

   public :: htg_set_ok
   public :: htg_fail

contains

   pure logical function status_ok(self) result(is_ok)
      class(htg_status_t), intent(in) :: self
      is_ok = (self%code == HTG_OK)
   end function status_ok

   pure subroutine status_clear(self)
      class(htg_status_t), intent(inout) :: self
      self%code = HTG_OK
      self%message = ''
   end subroutine status_clear

   pure subroutine status_set(self, code, msg)
      class(htg_status_t), intent(inout) :: self
      integer(i32_k), intent(in) :: code
      character(len=*), intent(in), optional :: msg
      self%code = code
      if (present(msg)) then
         self%message = msg
      else
         self%message = ''
      end if
   end subroutine status_set

   pure subroutine htg_set_ok(status)
      type(htg_status_t), intent(out), optional :: status
      if (present(status)) then
         status%code = HTG_OK
         status%message = ''
      end if
   end subroutine htg_set_ok

   pure subroutine htg_fail(code, msg, status)
      integer(i32_k), intent(in) :: code
      character(len=*), intent(in) :: msg
      type(htg_status_t), intent(out), optional :: status
      if (present(status)) then
         status%code = code
         status%message = msg
      else
         error stop msg
      end if
   end subroutine htg_fail

end module fvtkhtg_kinds_m
