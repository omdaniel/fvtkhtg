module fvtkhtg_bit_array_m
   use fvtkhtg_kinds_m, only: byte_k, i32_k, id_k, htg_status_t, &
      HTG_ERR_INVALID_ARG, HTG_ERR_OUT_OF_BOUNDS, htg_fail, htg_set_ok
   implicit none
   private

   integer(id_k), parameter :: MIN_BYTE_GROWTH = 8_id_k

   ! MSB-first packed bit array: bit `id` lives in `bytes(id / 8)` at bit position
   ! `7 - mod(id, 8)`. Trailing unused bits in the final active byte are always 0.
   type, public :: htg_bit_array_t
      integer(id_k) :: num_bits = 0_id_k
      integer(id_k) :: byte_capacity = 0_id_k
      integer(byte_k), allocatable :: bytes(:)
   contains
      procedure :: initialize => bit_array_initialize
      procedure :: clear => bit_array_clear
      procedure :: resize => bit_array_resize
      procedure :: set_number_of_values => bit_array_resize
      procedure :: get_number_of_values => bit_array_get_number_of_values
      procedure :: get_number_of_bytes => bit_array_get_number_of_bytes
      procedure :: get_byte_capacity => bit_array_get_byte_capacity
      procedure :: set_value => bit_array_set_value
      procedure :: insert_value => bit_array_insert_value
      procedure :: insert_next_value => bit_array_insert_next_value
      procedure :: get_value => bit_array_get_value
      procedure :: fill => bit_array_fill
      procedure :: get_raw_bytes => bit_array_get_raw_bytes
   end type htg_bit_array_t

contains

   pure function active_byte_count(nbits) result(nbytes)
      integer(id_k), intent(in) :: nbits
      integer(id_k) :: nbytes
      if (nbits <= 0_id_k) then
         nbytes = 0_id_k
      else
         nbytes = (nbits + 7_id_k) / 8_id_k
      end if
   end function active_byte_count

   pure subroutine zero_trailing_bits(self)
      type(htg_bit_array_t), intent(inout) :: self
      integer(id_k) :: active_bytes, last_idx
      integer(i32_k) :: rem, pos

      if (self%num_bits <= 0_id_k) return
      active_bytes = active_byte_count(self%num_bits)
      rem = int(mod(self%num_bits, 8_id_k), i32_k)
      if (rem /= 0_i32_k) then
         last_idx = active_bytes - 1_id_k
         do pos = 0_i32_k, 7_i32_k - rem
            self%bytes(last_idx) = ibclr(self%bytes(last_idx), pos)
         end do
      end if
   end subroutine zero_trailing_bits

   pure subroutine bit_array_clear(self)
      class(htg_bit_array_t), intent(inout) :: self
      if (allocated(self%bytes)) deallocate(self%bytes)
      self%num_bits = 0_id_k
      self%byte_capacity = 0_id_k
   end subroutine bit_array_clear

   pure subroutine bit_array_initialize(self, num_bits, init_value, status)
      class(htg_bit_array_t), intent(inout) :: self
      integer(id_k), intent(in), optional :: num_bits
      logical, intent(in), optional :: init_value
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: target_bits, nbytes
      logical :: fill_val

      target_bits = 0_id_k
      if (present(num_bits)) target_bits = num_bits

      fill_val = .false.
      if (present(init_value)) fill_val = init_value

      if (target_bits < 0_id_k) then
         call self%clear()
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_bit_array_t%initialize: num_bits must be non-negative', status)
         return
      end if

      call self%clear()
      if (target_bits > 0_id_k) then
         nbytes = active_byte_count(target_bits)
         allocate(self%bytes(0_id_k : nbytes - 1_id_k))
         self%num_bits = target_bits
         self%byte_capacity = nbytes
         if (fill_val) then
            self%bytes = -1_byte_k
            call zero_trailing_bits(self)
         else
            self%bytes = 0_byte_k
         end if
      end if

      call htg_set_ok(status)
   end subroutine bit_array_initialize

   pure subroutine bit_array_resize(self, num_bits, status)
      class(htg_bit_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: num_bits
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: old_active, needed_bytes, new_cap
      integer(byte_k), allocatable :: tmp(:)

      if (num_bits < 0_id_k) then
         call htg_fail(HTG_ERR_INVALID_ARG, &
            'htg_bit_array_t%resize: num_bits must be non-negative', status)
         return
      end if

      old_active = active_byte_count(self%num_bits)
      needed_bytes = active_byte_count(num_bits)

      if (needed_bytes > self%byte_capacity) then
         new_cap = max(needed_bytes, max(MIN_BYTE_GROWTH, 2_id_k * self%byte_capacity))
         allocate(tmp(0_id_k : new_cap - 1_id_k))
         tmp = 0_byte_k
         if (old_active > 0_id_k .and. allocated(self%bytes)) then
            tmp(0_id_k : old_active - 1_id_k) = self%bytes(0_id_k : old_active - 1_id_k)
         end if
         call move_alloc(tmp, self%bytes)
         self%byte_capacity = new_cap
      end if

      if (num_bits < self%num_bits) then
         if (needed_bytes < old_active .and. allocated(self%bytes)) then
            self%bytes(needed_bytes : old_active - 1_id_k) = 0_byte_k
         end if
         self%num_bits = num_bits
         call zero_trailing_bits(self)
      else
         self%num_bits = num_bits
      end if

      call htg_set_ok(status)
   end subroutine bit_array_resize

   pure function bit_array_get_number_of_values(self) result(nbits)
      class(htg_bit_array_t), intent(in) :: self
      integer(id_k) :: nbits
      nbits = self%num_bits
   end function bit_array_get_number_of_values

   pure function bit_array_get_number_of_bytes(self) result(nbytes)
      class(htg_bit_array_t), intent(in) :: self
      integer(id_k) :: nbytes
      nbytes = active_byte_count(self%num_bits)
   end function bit_array_get_number_of_bytes

   pure function bit_array_get_byte_capacity(self) result(cap)
      class(htg_bit_array_t), intent(in) :: self
      integer(id_k) :: cap
      cap = self%byte_capacity
   end function bit_array_get_byte_capacity

   pure subroutine bit_array_set_value(self, id, val, status)
      class(htg_bit_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: id
      logical, intent(in) :: val
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: byte_idx
      integer(i32_k) :: bit_pos

      if (id < 0_id_k .or. id >= self%num_bits) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_bit_array_t%set_value: bit id out of bounds', status)
         return
      end if

      byte_idx = id / 8_id_k
      bit_pos = 7_i32_k - int(mod(id, 8_id_k), i32_k)
      if (val) then
         self%bytes(byte_idx) = ibset(self%bytes(byte_idx), bit_pos)
      else
         self%bytes(byte_idx) = ibclr(self%bytes(byte_idx), bit_pos)
      end if

      call htg_set_ok(status)
   end subroutine bit_array_set_value

   pure subroutine bit_array_insert_value(self, id, val, status)
      class(htg_bit_array_t), intent(inout) :: self
      integer(id_k), intent(in) :: id
      logical, intent(in) :: val
      type(htg_status_t), intent(out), optional :: status

      type(htg_status_t) :: resize_status
      integer(id_k) :: byte_idx
      integer(i32_k) :: bit_pos

      if (id < 0_id_k) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_bit_array_t%insert_value: bit id must be non-negative', status)
         return
      end if

      if (id >= self%num_bits) then
         call self%resize(id + 1_id_k, resize_status)
         if (.not. resize_status%ok()) then
            call htg_fail(resize_status%code, trim(resize_status%message), status)
            return
         end if
      end if

      byte_idx = id / 8_id_k
      bit_pos = 7_i32_k - int(mod(id, 8_id_k), i32_k)
      if (val) then
         self%bytes(byte_idx) = ibset(self%bytes(byte_idx), bit_pos)
      else
         self%bytes(byte_idx) = ibclr(self%bytes(byte_idx), bit_pos)
      end if

      call htg_set_ok(status)
   end subroutine bit_array_insert_value

   pure subroutine bit_array_insert_next_value(self, val, inserted_id, status)
      class(htg_bit_array_t), intent(inout) :: self
      logical, intent(in) :: val
      integer(id_k), intent(out), optional :: inserted_id
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: next_id

      next_id = self%num_bits
      call self%insert_value(next_id, val, status)
      if (present(status)) then
         if (.not. status%ok()) then
            if (present(inserted_id)) inserted_id = -1_id_k
            return
         end if
      end if
      if (present(inserted_id)) inserted_id = next_id
   end subroutine bit_array_insert_next_value

   logical function bit_array_get_value(self, id, status) result(val)
      class(htg_bit_array_t), intent(in) :: self
      integer(id_k), intent(in) :: id
      type(htg_status_t), intent(out), optional :: status

      integer(id_k) :: byte_idx
      integer(i32_k) :: bit_pos

      val = .false.
      if (id < 0_id_k .or. id >= self%num_bits) then
         call htg_fail(HTG_ERR_OUT_OF_BOUNDS, &
            'htg_bit_array_t%get_value: bit id out of bounds', status)
         return
      end if

      byte_idx = id / 8_id_k
      bit_pos = 7_i32_k - int(mod(id, 8_id_k), i32_k)
      val = btest(self%bytes(byte_idx), bit_pos)
      call htg_set_ok(status)
   end function bit_array_get_value

   pure subroutine bit_array_fill(self, val)
      class(htg_bit_array_t), intent(inout) :: self
      logical, intent(in) :: val

      integer(id_k) :: nbytes

      if (self%num_bits <= 0_id_k) return
      nbytes = active_byte_count(self%num_bits)
      if (val) then
         self%bytes(0_id_k : nbytes - 1_id_k) = -1_byte_k
         call zero_trailing_bits(self)
      else
         self%bytes(0_id_k : nbytes - 1_id_k) = 0_byte_k
      end if
   end subroutine bit_array_fill

   subroutine bit_array_get_raw_bytes(self, ptr)
      class(htg_bit_array_t), target, intent(inout) :: self
      integer(byte_k), pointer, intent(out) :: ptr(:)

      integer(id_k) :: nbytes

      nbytes = active_byte_count(self%num_bits)
      if (nbytes <= 0_id_k .or. .not. allocated(self%bytes)) then
         ptr => null()
      else
         ptr(0_id_k : nbytes - 1_id_k) => self%bytes(0_id_k : nbytes - 1_id_k)
      end if
   end subroutine bit_array_get_raw_bytes

end module fvtkhtg_bit_array_m
