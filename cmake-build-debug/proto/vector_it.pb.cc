// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: vector_it.proto

#include "vector_it.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
extern PROTOBUF_INTERNAL_EXPORT_vector_5fit_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_vector_vector_5fit_2eproto;
namespace vector {
class vectorDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<vector> _instance;
} _vector_default_instance_;
class full_iter_vecDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<full_iter_vec> _instance;
} _full_iter_vec_default_instance_;
}  // namespace vector
static void InitDefaultsscc_info_full_iter_vec_vector_5fit_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::vector::_full_iter_vec_default_instance_;
    new (ptr) ::vector::full_iter_vec();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_full_iter_vec_vector_5fit_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_full_iter_vec_vector_5fit_2eproto}, {
      &scc_info_vector_vector_5fit_2eproto.base,}};

static void InitDefaultsscc_info_vector_vector_5fit_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::vector::_vector_default_instance_;
    new (ptr) ::vector::vector();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_vector_vector_5fit_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_vector_vector_5fit_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_vector_5fit_2eproto[2];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_vector_5fit_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_vector_5fit_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_vector_5fit_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::vector::vector, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::vector::vector, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::vector::vector, iter_),
  PROTOBUF_FIELD_OFFSET(::vector::vector, vec_value_),
  0,
  ~0u,
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::vector::full_iter_vec, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::vector::full_iter_vec, vec_t_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, sizeof(::vector::vector)},
  { 9, -1, sizeof(::vector::full_iter_vec)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::vector::_vector_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::vector::_full_iter_vec_default_instance_),
};

const char descriptor_table_protodef_vector_5fit_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\017vector_it.proto\022\006vector\")\n\006vector\022\014\n\004i"
  "ter\030\001 \002(\005\022\021\n\tvec_value\030\002 \003(\002\".\n\rfull_ite"
  "r_vec\022\035\n\005vec_t\030\002 \003(\0132\016.vector.vector"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_vector_5fit_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_vector_5fit_2eproto_sccs[2] = {
  &scc_info_full_iter_vec_vector_5fit_2eproto.base,
  &scc_info_vector_vector_5fit_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_vector_5fit_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_vector_5fit_2eproto = {
  false, false, descriptor_table_protodef_vector_5fit_2eproto, "vector_it.proto", 116,
  &descriptor_table_vector_5fit_2eproto_once, descriptor_table_vector_5fit_2eproto_sccs, descriptor_table_vector_5fit_2eproto_deps, 2, 0,
  schemas, file_default_instances, TableStruct_vector_5fit_2eproto::offsets,
  file_level_metadata_vector_5fit_2eproto, 2, file_level_enum_descriptors_vector_5fit_2eproto, file_level_service_descriptors_vector_5fit_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_vector_5fit_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_vector_5fit_2eproto)), true);
namespace vector {

// ===================================================================

class vector::_Internal {
 public:
  using HasBits = decltype(std::declval<vector>()._has_bits_);
  static void set_has_iter(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000001) ^ 0x00000001) != 0;
  }
};

vector::vector(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  vec_value_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:vector.vector)
}
vector::vector(const vector& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_),
      vec_value_(from.vec_value_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  iter_ = from.iter_;
  // @@protoc_insertion_point(copy_constructor:vector.vector)
}

void vector::SharedCtor() {
  iter_ = 0;
}

vector::~vector() {
  // @@protoc_insertion_point(destructor:vector.vector)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void vector::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void vector::ArenaDtor(void* object) {
  vector* _this = reinterpret_cast< vector* >(object);
  (void)_this;
}
void vector::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void vector::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const vector& vector::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_vector_vector_5fit_2eproto.base);
  return *internal_default_instance();
}


void vector::Clear() {
// @@protoc_insertion_point(message_clear_start:vector.vector)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  vec_value_.Clear();
  iter_ = 0;
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* vector::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  _Internal::HasBits has_bits{};
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // required int32 iter = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 8)) {
          _Internal::set_has_iter(&has_bits);
          iter_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated float vec_value = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 21)) {
          ptr -= 1;
          do {
            ptr += 1;
            _internal_add_vec_value(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<float>(ptr));
            ptr += sizeof(float);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<21>(ptr));
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedFloatParser(_internal_mutable_vec_value(), ptr, ctx);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  _has_bits_.Or(has_bits);
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* vector::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:vector.vector)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required int32 iter = 1;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(1, this->_internal_iter(), target);
  }

  // repeated float vec_value = 2;
  for (int i = 0, n = this->_internal_vec_value_size(); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteFloatToArray(2, this->_internal_vec_value(i), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:vector.vector)
  return target;
}

size_t vector::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:vector.vector)
  size_t total_size = 0;

  // required int32 iter = 1;
  if (_internal_has_iter()) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
        this->_internal_iter());
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated float vec_value = 2;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_vec_value_size());
    size_t data_size = 4UL * count;
    total_size += 1 *
                  ::PROTOBUF_NAMESPACE_ID::internal::FromIntSize(this->_internal_vec_value_size());
    total_size += data_size;
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void vector::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:vector.vector)
  GOOGLE_DCHECK_NE(&from, this);
  const vector* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<vector>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:vector.vector)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:vector.vector)
    MergeFrom(*source);
  }
}

void vector::MergeFrom(const vector& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:vector.vector)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  vec_value_.MergeFrom(from.vec_value_);
  if (from._internal_has_iter()) {
    _internal_set_iter(from._internal_iter());
  }
}

void vector::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:vector.vector)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void vector::CopyFrom(const vector& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:vector.vector)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool vector::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void vector::InternalSwap(vector* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  vec_value_.InternalSwap(&other->vec_value_);
  swap(iter_, other->iter_);
}

::PROTOBUF_NAMESPACE_ID::Metadata vector::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

class full_iter_vec::_Internal {
 public:
};

full_iter_vec::full_iter_vec(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  vec_t_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:vector.full_iter_vec)
}
full_iter_vec::full_iter_vec(const full_iter_vec& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      vec_t_(from.vec_t_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:vector.full_iter_vec)
}

void full_iter_vec::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_full_iter_vec_vector_5fit_2eproto.base);
}

full_iter_vec::~full_iter_vec() {
  // @@protoc_insertion_point(destructor:vector.full_iter_vec)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void full_iter_vec::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void full_iter_vec::ArenaDtor(void* object) {
  full_iter_vec* _this = reinterpret_cast< full_iter_vec* >(object);
  (void)_this;
}
void full_iter_vec::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void full_iter_vec::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const full_iter_vec& full_iter_vec::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_full_iter_vec_vector_5fit_2eproto.base);
  return *internal_default_instance();
}


void full_iter_vec::Clear() {
// @@protoc_insertion_point(message_clear_start:vector.full_iter_vec)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  vec_t_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* full_iter_vec::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .vector.vector vec_t = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_vec_t(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<18>(ptr));
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* full_iter_vec::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:vector.full_iter_vec)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .vector.vector vec_t = 2;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_vec_t_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(2, this->_internal_vec_t(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:vector.full_iter_vec)
  return target;
}

size_t full_iter_vec::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:vector.full_iter_vec)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .vector.vector vec_t = 2;
  total_size += 1UL * this->_internal_vec_t_size();
  for (const auto& msg : this->vec_t_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void full_iter_vec::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:vector.full_iter_vec)
  GOOGLE_DCHECK_NE(&from, this);
  const full_iter_vec* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<full_iter_vec>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:vector.full_iter_vec)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:vector.full_iter_vec)
    MergeFrom(*source);
  }
}

void full_iter_vec::MergeFrom(const full_iter_vec& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:vector.full_iter_vec)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  vec_t_.MergeFrom(from.vec_t_);
}

void full_iter_vec::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:vector.full_iter_vec)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void full_iter_vec::CopyFrom(const full_iter_vec& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:vector.full_iter_vec)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool full_iter_vec::IsInitialized() const {
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(vec_t_)) return false;
  return true;
}

void full_iter_vec::InternalSwap(full_iter_vec* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  vec_t_.InternalSwap(&other->vec_t_);
}

::PROTOBUF_NAMESPACE_ID::Metadata full_iter_vec::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace vector
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::vector::vector* Arena::CreateMaybeMessage< ::vector::vector >(Arena* arena) {
  return Arena::CreateMessageInternal< ::vector::vector >(arena);
}
template<> PROTOBUF_NOINLINE ::vector::full_iter_vec* Arena::CreateMaybeMessage< ::vector::full_iter_vec >(Arena* arena) {
  return Arena::CreateMessageInternal< ::vector::full_iter_vec >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
