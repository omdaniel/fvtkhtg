module test_foundations_probe_m
   use, intrinsic :: iso_fortran_env, only: error_unit
   use fvtkhtg_kinds_m, only: byte_k, i32_k, id_k
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   use fvtkhtg_cell_data_m, only: htg_cell_data_t
   implicit none
   private

   integer(byte_k), parameter, public :: PROBE_BYTE_0X80 = ibset(0_byte_k, 7)
   public :: lldb_checkpoint

contains

   subroutine lldb_checkpoint(bits, cd)
      type(htg_bit_array_t), intent(in) :: bits
      type(htg_cell_data_t), intent(in) :: cd

      if (bits%num_bits /= 10_id_k .or. &
          bits%bytes(0) /= -127_byte_k .or. &
          bits%bytes(1) /= PROBE_BYTE_0X80 .or. &
          cd%num_arrays /= 1_i32_k) then
         write(error_unit, '(a)') 'ASSERTION FAILED: lldb_checkpoint state mismatch'
         error stop 1
      end if
   end subroutine lldb_checkpoint

end module test_foundations_probe_m

program test_foundations
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use fvtkhtg_kinds_m, only: byte_k, i32_k, i64_k, id_k, r32_k, r64_k, &
      htg_status_t, HTG_OK, HTG_ERR_INVALID_ARG, HTG_ERR_OUT_OF_BOUNDS, &
      HTG_ERR_INVALID_STATE, HTG_ERR_TOPOLOGY, HTG_ERR_INDEXING, &
      HTG_ERR_MASK, HTG_ERR_INTERFACE, HTG_ERR_NOT_FOUND, &
      HTG_INDEX_UNSET, HTG_INDEX_IMPLICIT, HTG_INDEX_EXPLICIT, &
      HTG_TYPE_REAL64, HTG_TYPE_REAL32, HTG_TYPE_INT64, HTG_TYPE_INT32, &
      HTG_TYPE_BIT, HTG_INVALID_ID, HTG_LEAF_ELDER_SENTINEL, HTG_UINT32_MAX, &
      htg_set_ok, htg_fail
   use fvtkhtg_bit_array_m, only: htg_bit_array_t
   use fvtkhtg_cell_data_m, only: htg_field_array_t, htg_cell_data_t
   use test_foundations_probe_m, only: lldb_checkpoint
   implicit none

   integer(byte_k), parameter :: BYTE_0X80 = ibset(0_byte_k, 7)

   call test_kinds_and_status()
   call test_bit_array_msb_patterns()
   call test_bit_array_resize_and_insert()
   call test_field_array_all_types()
   call test_cell_data_container_and_stability()
   call test_negative_error_paths()
   call run_lldb_inspection_fixture()

   write(output_unit, '(a)') 'PASS: all foundations unit tests succeeded.'

contains

   subroutine assert_true(cond, msg)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) then
         write(error_unit, '(a,a)') 'ASSERTION FAILED: ', msg
         error stop 1
      end if
   end subroutine assert_true

   subroutine test_kinds_and_status()
      type(htg_status_t) :: st

      call assert_true(storage_size(0_byte_k) == 8, 'byte_k must be 8 bits')
      call assert_true(storage_size(0_i32_k) == 32, 'i32_k must be 32 bits')
      call assert_true(storage_size(0_i64_k) == 64, 'i64_k must be 64 bits')
      call assert_true(storage_size(0_id_k) == 64, 'id_k must be 64 bits')
      call assert_true(storage_size(0.0_r32_k) == 32, 'r32_k must be 32 bits')
      call assert_true(storage_size(0.0_r64_k) == 64, 'r64_k must be 64 bits')

      call assert_true(HTG_OK == 0_i32_k, 'HTG_OK == 0')
      call assert_true(HTG_ERR_INVALID_ARG == 1_i32_k, 'HTG_ERR_INVALID_ARG == 1')
      call assert_true(HTG_ERR_OUT_OF_BOUNDS == 2_i32_k, 'HTG_ERR_OUT_OF_BOUNDS == 2')
      call assert_true(HTG_ERR_INVALID_STATE == 3_i32_k, 'HTG_ERR_INVALID_STATE == 3')
      call assert_true(HTG_ERR_TOPOLOGY == 4_i32_k, 'HTG_ERR_TOPOLOGY == 4')
      call assert_true(HTG_ERR_INDEXING == 5_i32_k, 'HTG_ERR_INDEXING == 5')
      call assert_true(HTG_ERR_MASK == 6_i32_k, 'HTG_ERR_MASK == 6')
      call assert_true(HTG_ERR_INTERFACE == 7_i32_k, 'HTG_ERR_INTERFACE == 7')
      call assert_true(HTG_ERR_NOT_FOUND == 8_i32_k, 'HTG_ERR_NOT_FOUND == 8')

      call assert_true(HTG_INDEX_UNSET == 0_i32_k, 'HTG_INDEX_UNSET == 0')
      call assert_true(HTG_INDEX_IMPLICIT == 1_i32_k, 'HTG_INDEX_IMPLICIT == 1')
      call assert_true(HTG_INDEX_EXPLICIT == 2_i32_k, 'HTG_INDEX_EXPLICIT == 2')

      call assert_true(HTG_INVALID_ID == -1_id_k, 'HTG_INVALID_ID == -1')
      call assert_true(HTG_LEAF_ELDER_SENTINEL == -1_id_k, 'HTG_LEAF_ELDER_SENTINEL == -1')
      call assert_true(HTG_UINT32_MAX == 4294967295_id_k, 'HTG_UINT32_MAX == 2^32 - 1')

      call assert_true(st%ok(), 'default status is ok')
      call st%set(HTG_ERR_TOPOLOGY, 'cycle detected')
      call assert_true(.not. st%ok(), 'status after set error is not ok')
      call assert_true(st%code == HTG_ERR_TOPOLOGY, 'status code matches')
      call assert_true(trim(st%message) == 'cycle detected', 'status message matches')

      call st%clear()
      call assert_true(st%ok() .and. len_trim(st%message) == 0, 'status clear resets state')

      call htg_fail(HTG_ERR_MASK, 'mask fault', st)
      call assert_true(st%code == HTG_ERR_MASK, 'htg_fail populates status code')
      call assert_true(trim(st%message) == 'mask fault', 'htg_fail populates message')

      call htg_set_ok(st)
      call assert_true(st%ok() .and. len_trim(st%message) == 0, 'htg_set_ok clears status')
   end subroutine test_kinds_and_status

   subroutine test_bit_array_msb_patterns()
      type(htg_bit_array_t), target :: ba
      type(htg_bit_array_t) :: ba_copy
      integer(byte_k), pointer :: raw(:)
      type(htg_status_t) :: st
      integer(id_k) :: i

      call ba%initialize(status=st)
      call assert_true(st%ok(), 'empty bit_array initialize ok')
      call assert_true(ba%get_number_of_values() == 0_id_k, 'empty num_bits == 0')
      call assert_true(ba%get_number_of_bytes() == 0_id_k, 'empty num_bytes == 0')
      call ba%get_raw_bytes(raw)
      call assert_true(.not. associated(raw), 'empty raw pointer is null')

      ! Single byte MSB-first verification: bit 0 -> 0x80 (-128_int8), bit 7 -> 0x01 (1_int8)
      call ba%initialize(8_id_k, .false., st)
      call assert_true(st%ok(), '8-bit initialize ok')
      call ba%get_raw_bytes(raw)
      call assert_true(associated(raw), '8-bit raw pointer associated')
      call assert_true(lbound(raw, 1) == 0 .and. ubound(raw, 1) == 0, '8-bit raw bounds 0:0')
      call assert_true(raw(0) == 0_byte_k, 'initial byte is 0')

      call ba%set_value(0_id_k, .true., st)
      call assert_true(st%ok() .and. raw(0) == BYTE_0X80, 'bit 0 sets MSB 0x80 (-128_int8)')
      call assert_true(ba%get_value(0_id_k, st) .and. st%ok(), 'get_value(0) is true')

      call ba%set_value(0_id_k, .false., st)
      call ba%set_value(7_id_k, .true., st)
      call assert_true(st%ok() .and. raw(0) == 1_byte_k, 'bit 7 sets LSB 0x01 (1_int8)')

      call ba%set_value(1_id_k, .true., st)
      call assert_true(raw(0) == 65_byte_k, 'bits 1 and 7 yield 0x41 (65_int8)')

      ! Cross-byte boundaries: 7 -> 8 -> 63 -> 64 -> 65
      call ba%initialize(66_id_k, .false., st)
      call assert_true(st%ok(), '66-bit initialize ok')
      call assert_true(ba%get_number_of_bytes() == 9_id_k, '66 bits = 9 active bytes')
      call ba%set_value(7_id_k, .true., st)
      call ba%set_value(8_id_k, .true., st)
      call ba%set_value(63_id_k, .true., st)
      call ba%set_value(64_id_k, .true., st)
      call ba%set_value(65_id_k, .true., st)

      call ba%get_raw_bytes(raw)
      call assert_true(lbound(raw, 1) == 0 .and. ubound(raw, 1) == 8, '66-bit raw bounds 0:8')
      call assert_true(raw(0) == 1_byte_k, 'byte 0 has bit 7 -> 0x01')
      call assert_true(raw(1) == BYTE_0X80, 'byte 1 has bit 8 -> 0x80')
      do i = 2_id_k, 6_id_k
         call assert_true(raw(i) == 0_byte_k, 'intermediate bytes 2..6 are 0')
      end do
      call assert_true(raw(7) == 1_byte_k, 'byte 7 has bit 63 -> 0x01')
      call assert_true(raw(8) == -64_byte_k, 'byte 8 has bits 64,65 -> 0xC0 (-64_int8)')

      ! Non-byte-multiple init_value=.true. must keep trailing slack bits zeroed
      call ba%initialize(10_id_k, .true., st)
      call ba%get_raw_bytes(raw)
      call assert_true(raw(0) == -1_byte_k, '10-bit true init byte 0 == 0xFF')
      call assert_true(raw(1) == -64_byte_k, '10-bit true init byte 1 == 0xC0 (trailing 6 bits 0)')

      ! Intrinsic copy preserves 0: lower bound
      ba_copy = ba
      call assert_true(lbound(ba_copy%bytes, 1) == 0, 'intrinsic copy preserves 0: lbound')
      call assert_true(ba_copy%bytes(0) == -1_byte_k .and. ba_copy%bytes(1) == -64_byte_k, &
         'intrinsic copy preserves exact packed bytes')
   end subroutine test_bit_array_msb_patterns

   subroutine test_bit_array_resize_and_insert()
      type(htg_bit_array_t), target :: ba
      integer(byte_k), pointer :: raw(:)
      type(htg_status_t) :: st
      integer(id_k) :: id_out, i

      call ba%initialize(13_id_k, .true., st)
      call ba%get_raw_bytes(raw)
      ! 13 bits: byte 0 = 0xFF (-1), byte 1 = 1111_1000 = 0xF8 (-8_int8)
      call assert_true(raw(0) == -1_byte_k .and. raw(1) == -8_byte_k, &
         '13-bit fill true has 3 trailing zero bits (0xF8)')

      ! Shrink to 5 bits: byte 0 must become 1111_1000 (-8_int8), byte 1 zeroed
      call ba%resize(5_id_k, st)
      call assert_true(st%ok() .and. ba%get_number_of_values() == 5_id_k, 'shrunk to 5 bits')
      call ba%get_raw_bytes(raw)
      call assert_true(size(raw) == 1 .and. raw(0) == -8_byte_k, &
         'shrunk 5-bit byte 0 cleared bits 5..7 to 0 (0xF8)')

      ! Grow to 20 bits via set_number_of_values: bits 0..4 stay true, bits 5..19 must be false
      call ba%set_number_of_values(20_id_k, st)
      call assert_true(st%ok() .and. ba%get_number_of_values() == 20_id_k, 'grown to 20 bits')
      call assert_true(ba%get_byte_capacity() >= 8_id_k, 'geometric byte capacity >= 8')
      do i = 0_id_k, 4_id_k
         call assert_true(ba%get_value(i, st), 'preserved bit 0..4 is true')
      end do
      do i = 5_id_k, 19_id_k
         call assert_true(.not. ba%get_value(i, st), 'newly grown bit 5..19 is false')
      end do

      ! fill(.false.) then sparse insert_value and insert_next_value
      call ba%fill(.false.)
      call ba%resize(0_id_k, st)
      call ba%insert_next_value(.true., id_out, st)
      call assert_true(st%ok() .and. id_out == 0_id_k, 'first insert_next_value -> id 0')
      call ba%insert_next_value(.false., status=st)
      call assert_true(st%ok() .and. ba%get_number_of_values() == 2_id_k, &
         'insert_next_value without optional inserted_id')
      call ba%insert_value(64_id_k, .true., st)
      call assert_true(st%ok() .and. ba%get_number_of_values() == 65_id_k, &
         'insert_value(64) grows num_bits to 65')
      call assert_true(ba%get_value(0_id_k), 'sparse insert bit 0 is true')
      call assert_true(.not. ba%get_value(1_id_k), 'sparse insert bit 1 is false')
      call assert_true(.not. ba%get_value(63_id_k), 'sparse insert bit 63 is false')
      call assert_true(ba%get_value(64_id_k), 'sparse insert bit 64 is true')
   end subroutine test_bit_array_resize_and_insert

   subroutine test_field_array_all_types()
      type(htg_field_array_t), target :: f_r64, f_r32, f_i64, f_i32, f_bit
      real(r64_k) :: tup_r64_in(0:2), tup_r64_out(0:2)
      real(r32_k) :: tup_r32_in(0:1), tup_r32_out(0:1)
      integer(i64_k) :: tup_i64_in(0:1), tup_i64_out(0:1)
      integer(i32_k) :: tup_i32_in(0:2), tup_i32_out(0:2)
      logical :: tup_bit_in(0:2), tup_bit_out(0:2)
      real(r64_k), pointer :: view_r64(:, :)
      real(r32_k), pointer :: view_r32(:, :)
      integer(i64_k), pointer :: view_i64(:, :)
      integer(i32_k), pointer :: view_i32(:, :)
      type(htg_bit_array_t), pointer :: view_bit
      integer(byte_k), pointer :: view_bit_bytes(:)
      type(htg_status_t) :: st

      ! 1. REAL64 3-component normals with amortized growth
      call f_r64%initialize('InterfaceNormals', HTG_TYPE_REAL64, 3_i32_k, 2_id_k, st)
      call assert_true(st%ok(), 'r64 3-comp initialize')
      tup_r64_in = [0.0_r64_k, 1.0_r64_k, 0.0_r64_k]
      call f_r64%set_tuple_r64(0_id_k, tup_r64_in, st)
      call f_r64%set_component_r64(1_id_k, 0_i32_k, 0.6_r64_k, st)
      call f_r64%set_component_r64(1_id_k, 1_i32_k, 0.8_r64_k, st)
      call f_r64%set_component_r64(1_id_k, 2_i32_k, 0.0_r64_k, st)
      call f_r64%ensure_tuples(10_id_k, st)
      call assert_true(st%ok() .and. f_r64%get_number_of_tuples() == 10_id_k, &
         'r64 ensure_tuples(10)')
      call assert_true(f_r64%get_tuple_capacity() >= 10_id_k, 'r64 tuple_capacity >= 10')
      call f_r64%get_tuple_r64(0_id_k, tup_r64_out, st)
      call assert_true(all(abs(tup_r64_out - tup_r64_in) <= 0.0_r64_k), &
         'r64 tuple 0 preserved across growth')
      call f_r64%get_view_r64(view_r64, st)
      call assert_true(st%ok() .and. associated(view_r64), 'r64 view associated')
      call assert_true(lbound(view_r64, 1) == 0 .and. ubound(view_r64, 1) == 2, &
         'r64 view dim 1 bounds 0:2')
      call assert_true(lbound(view_r64, 2) == 0 .and. ubound(view_r64, 2) == 9, &
         'r64 view dim 2 bounds 0:9')
      call assert_true(abs(view_r64(1, 1) - 0.8_r64_k) <= 0.0_r64_k .and. &
         abs(view_r64(2, 9)) <= 0.0_r64_k, 'r64 view values and zero-init')

      ! Shrink and re-grow clears truncated entries
      call f_r64%set_component_r64(9_id_k, 2_i32_k, 42.0_r64_k, st)
      call f_r64%set_number_of_tuples(3_id_k, st)
      call f_r64%ensure_tuples(10_id_k, st)
      call assert_true(abs(f_r64%get_component_r64(9_id_k, 2_i32_k, st)) <= 0.0_r64_k, &
         'shrunk then regrown tuple is zeroed')
      call f_r64%set_scalar_r64(4_id_k, 3.25_r64_k, st)
      call assert_true(abs(f_r64%get_scalar_r64(4_id_k, st) - 3.25_r64_k) <= 0.0_r64_k, &
         'r64 scalar convenience roundtrip')

      ! 2. REAL32 2-component field
      call f_r32%initialize('Velocity2D', HTG_TYPE_REAL32, 2_i32_k, 3_id_k, st)
      tup_r32_in = [1.5_r32_k, -2.5_r32_k]
      call f_r32%set_tuple_r32(2_id_k, tup_r32_in, st)
      call f_r32%get_tuple_r32(2_id_k, tup_r32_out, st)
      call assert_true(all(abs(tup_r32_out - tup_r32_in) <= 0.0_r32_k), 'r32 tuple roundtrip')
      call f_r32%set_scalar_r32(0_id_k, 7.5_r32_k, st)
      call assert_true(abs(f_r32%get_scalar_r32(0_id_k, st) - 7.5_r32_k) <= 0.0_r32_k, &
         'r32 scalar convenience roundtrip')
      call f_r32%get_view_r32(view_r32, st)
      call assert_true(lbound(view_r32, 1) == 0 .and. lbound(view_r32, 2) == 0, &
         'r32 view 0-based bounds')

      ! 3. INT64 2-component field
      call f_i64%initialize('GlobalPair', HTG_TYPE_INT64, 2_i32_k, 4_id_k, st)
      tup_i64_in = [5000000000_i64_k, -7_i64_k]
      call f_i64%set_tuple_i64(3_id_k, tup_i64_in, st)
      call f_i64%get_tuple_i64(3_id_k, tup_i64_out, st)
      call assert_true(all(tup_i64_out == tup_i64_in), 'i64 tuple roundtrip')
      call f_i64%set_scalar_i64(0_id_k, 99_i64_k, st)
      call assert_true(f_i64%get_scalar_i64(0_id_k, st) == 99_i64_k, 'i64 scalar roundtrip')
      call f_i64%get_view_i64(view_i64, st)
      call assert_true(lbound(view_i64, 1) == 0 .and. ubound(view_i64, 2) == 3, &
         'i64 view 0-based bounds')

      ! 4. INT32 3-component field
      call f_i32%initialize('MaterialIds', HTG_TYPE_INT32, 3_i32_k, 2_id_k, st)
      tup_i32_in = [10_i32_k, 20_i32_k, 30_i32_k]
      call f_i32%set_tuple_i32(1_id_k, tup_i32_in, st)
      call f_i32%get_tuple_i32(1_id_k, tup_i32_out, st)
      call assert_true(all(tup_i32_out == tup_i32_in), 'i32 tuple roundtrip')
      call f_i32%set_scalar_i32(0_id_k, -42_i32_k, st)
      call assert_true(f_i32%get_scalar_i32(0_id_k, st) == -42_i32_k, 'i32 scalar roundtrip')
      call f_i32%get_view_i32(view_i32, st)
      call assert_true(view_i32(2, 1) == 30_i32_k, 'i32 view matches component 2')

      ! 5. BIT 3-component field (3 * 3 = 9 bits -> 2 active bytes)
      call f_bit%initialize('Flags', HTG_TYPE_BIT, 3_i32_k, 3_id_k, st)
      tup_bit_in = [.true., .false., .true.]
      call f_bit%set_tuple_bit(0_id_k, tup_bit_in, st)
      call f_bit%set_scalar_bit(1_id_k, .true., st)
      call assert_true(f_bit%get_scalar_bit(1_id_k, st), 'bit scalar convenience roundtrip')
      call f_bit%set_scalar_bit(1_id_k, .false., st)
      call f_bit%set_component_bit(2_id_k, 2_i32_k, .true., st)
      call f_bit%get_tuple_bit(0_id_k, tup_bit_out, st)
      call assert_true(all(tup_bit_out .eqv. tup_bit_in), 'bit tuple 0 roundtrip')
      call assert_true(f_bit%get_component_bit(2_id_k, 2_i32_k, st), 'bit (2,2) -> bit 8 is true')
      call f_bit%get_view_bit(view_bit, st)
      call assert_true(associated(view_bit) .and. view_bit%get_number_of_values() == 9_id_k, &
         'bit view has 9 bits')
      call f_bit%get_view_bit_bytes(view_bit_bytes, st)
      call assert_true(lbound(view_bit_bytes, 1) == 0 .and. size(view_bit_bytes) == 2, &
         'bit byte view bounds 0:1')
      ! Bits 0 and 2 set in byte 0 -> 1010_0000 = -96_int8 (0xA0); bit 8 set in byte 1 -> -128_int8
      call assert_true(view_bit_bytes(0) == -96_byte_k .and. view_bit_bytes(1) == BYTE_0X80, &
         'bit field MSB byte packing matches 0xA0, 0x80')
   end subroutine test_field_array_all_types

   subroutine test_cell_data_container_and_stability()
      type(htg_cell_data_t), target :: cd
      type(htg_field_array_t), pointer :: normals_ptr, intercepts_ptr, p_idx
      type(htg_status_t) :: st
      character(len=32) :: extra_name
      integer(i32_k) :: k

      call cd%add_array('InterfaceNormals', HTG_TYPE_REAL64, 3_i32_k, 4_id_k, st)
      call assert_true(st%ok(), 'add InterfaceNormals')
      call cd%add_array('InterfaceIntercepts', HTG_TYPE_REAL64, 3_i32_k, 4_id_k, st)
      call assert_true(st%ok(), 'add InterfaceIntercepts')
      call cd%add_array('Density', HTG_TYPE_REAL64, 1_i32_k, 4_id_k, st)
      call assert_true(st%ok(), 'add Density')

      call cd%set_active_scalars('Density', st)
      call assert_true(st%ok() .and. cd%get_active_scalars() == 'Density', 'active scalars set')

      call cd%get_array('InterfaceNormals', normals_ptr, st)
      call assert_true(st%ok() .and. associated(normals_ptr), 'get_array InterfaceNormals')
      call normals_ptr%set_component_r64(0_id_k, 1_i32_k, -1.0_r64_k, st)

      ! Add 16 more arrays to force multiple array_capacity reallocations and verify normals_ptr
      ! remains completely valid and non-dangling!
      do k = 1_i32_k, 16_i32_k
         write(extra_name, '(a,i0)') 'ExtraField_', k
         call cd%add_array(trim(extra_name), HTG_TYPE_INT32, 1_i32_k, 4_id_k, st)
         call assert_true(st%ok(), 'add extra array')
      end do
      call assert_true(cd%get_number_of_arrays() == 19_i32_k, '19 arrays registered')
      call assert_true(normals_ptr%get_name() == 'InterfaceNormals', &
         'normals_ptr survived multiple container reallocations without dangling')
      call assert_true(abs(normals_ptr%get_component_r64(0_id_k, 1_i32_k, st) + 1.0_r64_k) <= &
         0.0_r64_k, 'normals_ptr data intact after container growth')

      ! Grow all 19 arrays simultaneously via ensure_tuples
      call cd%ensure_tuples(35_id_k, st)
      call assert_true(st%ok() .and. normals_ptr%get_number_of_tuples() == 35_id_k, &
         'cd%ensure_tuples grew all arrays')

      ! 0-based get_array_by_index and remove_array
      call cd%get_array_by_index(1_i32_k, intercepts_ptr, st)
      call assert_true(st%ok() .and. intercepts_ptr%get_name() == 'InterfaceIntercepts', &
         '0-based index 1 is InterfaceIntercepts')

      call cd%remove_array('Density', st)
      call assert_true(st%ok() .and. .not. cd%has_array('Density'), 'Density removed')
      call assert_true(.not. cd%has_active_scalars() .and. cd%get_active_scalars() == '', &
         'removing active scalar array clears active_scalars')

      call cd%get_array_by_index(2_i32_k, p_idx, st)
      call assert_true(st%ok() .and. p_idx%get_name() == 'ExtraField_1', &
         'subsequent arrays shifted cleanly after remove_array')

      ! Verify remove_array recomputes max num_tuples across remaining arrays
      call cd%add_array('LargeTemp', HTG_TYPE_REAL64, 1_i32_k, 120_id_k, st)
      call assert_true(st%ok() .and. cd%get_number_of_tuples() == 120_id_k, &
         'adding larger array updates cd num_tuples to 120')
      call cd%remove_array('LargeTemp', st)
      call assert_true(st%ok() .and. cd%get_number_of_tuples() == 35_id_k, &
         'removing largest array recomputes cd num_tuples back to 35')
      call assert_true(normals_ptr%get_name() == 'InterfaceNormals', &
         'normals_ptr still valid after add/remove')

      call cd%clear()
      call assert_true(cd%get_number_of_arrays() == 0_i32_k, 'cd%clear empties container')
   end subroutine test_cell_data_container_and_stability

   subroutine test_negative_error_paths()
      type(htg_bit_array_t) :: ba
      type(htg_field_array_t), target :: fa
      type(htg_cell_data_t), target :: cd
      type(htg_field_array_t), pointer :: fa_ptr
      real(r32_k), pointer :: bad_view(:, :)
      real(r64_k) :: bad_tuple(0:1)
      type(htg_status_t) :: st
      logical :: bval

      call ba%initialize(-5_id_k, .false., st)
      call assert_true(st%code == HTG_ERR_INVALID_ARG, 'negative num_bits rejected')

      call ba%initialize(4_id_k, .false., st)
      call ba%set_value(4_id_k, .true., st)
      call assert_true(st%code == HTG_ERR_OUT_OF_BOUNDS, 'set_value(num_bits) rejected')
      bval = ba%get_value(-1_id_k, st)
      call assert_true(.not. bval .and. st%code == HTG_ERR_OUT_OF_BOUNDS, &
         'get_value(-1) rejected')
      call ba%insert_value(-1_id_k, .true., st)
      call assert_true(st%code == HTG_ERR_OUT_OF_BOUNDS, 'insert_value(-1) rejected')

      ! Uninitialized field_array operations
      call fa%ensure_tuples(5_id_k, st)
      call assert_true(st%code == HTG_ERR_INVALID_STATE, 'uninitialized ensure_tuples rejected')

      call fa%initialize('   ', HTG_TYPE_REAL64, 1_i32_k, 1_id_k, st)
      call assert_true(st%code == HTG_ERR_INVALID_ARG, 'blank field name rejected')

      call fa%initialize('BadType', 99_i32_k, 1_i32_k, 1_id_k, st)
      call assert_true(st%code == HTG_ERR_INVALID_ARG, 'invalid data_type rejected')

      call fa%initialize('BadComp', HTG_TYPE_REAL64, 0_i32_k, 1_id_k, st)
      call assert_true(st%code == HTG_ERR_INVALID_ARG, '0 components rejected')

      call fa%initialize('ValidR64', HTG_TYPE_REAL64, 3_i32_k, 2_id_k, st)
      call assert_true(st%ok(), 'ValidR64 initialized')
      call fa%get_view_r32(bad_view, st)
      call assert_true(.not. associated(bad_view) .and. st%code == HTG_ERR_INVALID_STATE, &
         'wrong type view rejected')
      bad_tuple = [1.0_r64_k, 2.0_r64_k]
      call fa%set_tuple_r64(0_id_k, bad_tuple, st)
      call assert_true(st%code == HTG_ERR_INVALID_ARG, 'mismatched tuple slice length rejected')
      call fa%set_component_r64(2_id_k, 0_i32_k, 1.0_r64_k, st)
      call assert_true(st%code == HTG_ERR_OUT_OF_BOUNDS, 'tuple_id out of bounds rejected')
      call fa%set_component_r64(0_id_k, 3_i32_k, 1.0_r64_k, st)
      call assert_true(st%code == HTG_ERR_OUT_OF_BOUNDS, 'comp_id out of bounds rejected')

      ! CellData negative checks
      call cd%add_array('ValidR64', HTG_TYPE_REAL64, 1_i32_k, 2_id_k, st)
      call cd%add_array('ValidR64', HTG_TYPE_REAL64, 1_i32_k, 2_id_k, st)
      call assert_true(st%code == HTG_ERR_INVALID_ARG, 'duplicate array name rejected')

      call cd%get_array('Missing', fa_ptr, st)
      call assert_true(.not. associated(fa_ptr) .and. st%code == HTG_ERR_NOT_FOUND, &
         'missing get_array returns HTG_ERR_NOT_FOUND')

      call cd%remove_array('Missing', st)
      call assert_true(st%code == HTG_ERR_NOT_FOUND, 'missing remove_array returns NOT_FOUND')

      call cd%set_active_scalars('Missing', st)
      call assert_true(st%code == HTG_ERR_NOT_FOUND, &
         'missing set_active_scalars returns NOT_FOUND')

      call cd%get_array_by_index(1_i32_k, fa_ptr, st)
      call assert_true(.not. associated(fa_ptr) .and. st%code == HTG_ERR_OUT_OF_BOUNDS, &
         'out-of-bounds get_array_by_index rejected')
   end subroutine test_negative_error_paths

   subroutine run_lldb_inspection_fixture()
      type(htg_bit_array_t), target :: inspect_bits
      type(htg_cell_data_t), target :: inspect_cd
      type(htg_field_array_t), pointer :: normals
      type(htg_status_t) :: st

      ! Build deterministic bit pattern: 10 bits with bits 0, 7, 8 set
      ! Expected bytes(0) = 0x81 (-127_int8), bytes(1) = 0x80 (-128_int8)
      call inspect_bits%initialize(10_id_k, .false., st)
      call inspect_bits%set_value(0_id_k, .true., st)
      call inspect_bits%set_value(7_id_k, .true., st)
      call inspect_bits%set_value(8_id_k, .true., st)

      call inspect_cd%add_array('InterfaceNormals', HTG_TYPE_REAL64, 3_i32_k, 2_id_k, st)
      call inspect_cd%set_active_scalars('InterfaceNormals', st)
      call inspect_cd%get_array('InterfaceNormals', normals, st)
      call normals%set_tuple_r64(0_id_k, [1.0_r64_k, 0.0_r64_k, 0.0_r64_k], st)
      call normals%set_tuple_r64(1_id_k, [0.0_r64_k, 1.0_r64_k, 0.0_r64_k], st)

      call lldb_checkpoint(inspect_bits, inspect_cd)
   end subroutine run_lldb_inspection_fixture

end program test_foundations
