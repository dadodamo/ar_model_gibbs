// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: scalar_it.proto

#include "scalar_it.pb.h"

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
extern PROTOBUF_INTERNAL_EXPORT_scalar_5fit_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_scalar_scalar_5fit_2eproto;
namespace scalar {
class scalarDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<scalar> _instance;
} _scalar_default_instance_;
class full_scalar_itDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<full_scalar_it> _instance;
} _full_scalar_it_default_instance_;
}  // namespace scalar
static void InitDefaultsscc_info_full_scalar_it_scalar_5fit_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::scalar::_full_scalar_it_default_instance_;
    new (ptr) ::scalar::full_scalar_it();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_full_scalar_it_scalar_5fit_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_full_scalar_it_scalar_5fit_2eproto}, {
      &scc_info_scalar_scalar_5fit_2eproto.base,}};

static void InitDefaultsscc_info_scalar_scalar_5fit_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::scalar::_scalar_default_instance_;
    new (ptr) ::scalar::scalar();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_scalar_scalar_5fit_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_scalar_scalar_5fit_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_scalar_5fit_2eproto[2];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_scalar_5fit_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_scalar_5fit_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_scalar_5fit_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  PROTOBUF_FIELD_OFFSET(::scalar::scalar, _has_bits_),
  PROTOBUF_FIELD_OFFSET(::scalar::scalar, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::scalar::scalar, iter_),
  PROTOBUF_FIELD_OFFSET(::scalar::scalar, value_),
  1,
  0,
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::scalar::full_scalar_it, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::scalar::full_scalar_it, scalar_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, 7, sizeof(::scalar::scalar)},
  { 9, -1, sizeof(::scalar::full_scalar_it)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::scalar::_scalar_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::scalar::_full_scalar_it_default_instance_),
};

const char descriptor_table_protodef_scalar_5fit_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\017scalar_it.proto\022\006scalar\"%\n\006scalar\022\014\n\004i"
  "ter\030\001 \002(\005\022\r\n\005value\030\002 \002(\001\"0\n\016full_scalar_"
  "it\022\036\n\006scalar\030\001 \003(\0132\016.scalar.scalar"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_scalar_5fit_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_scalar_5fit_2eproto_sccs[2] = {
  &scc_info_full_scalar_it_scalar_5fit_2eproto.base,
  &scc_info_scalar_scalar_5fit_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_scalar_5fit_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_scalar_5fit_2eproto = {
  false, false, descriptor_table_protodef_scalar_5fit_2eproto, "scalar_it.proto", 114,
  &descriptor_table_scalar_5fit_2eproto_once, descriptor_table_scalar_5fit_2eproto_sccs, descriptor_table_scalar_5fit_2eproto_deps, 2, 0,
  schemas, file_default_instances, TableStruct_scalar_5fit_2eproto::offsets,
  file_level_metadata_scalar_5fit_2eproto, 2, file_level_enum_descriptors_scalar_5fit_2eproto, file_level_service_descriptors_scalar_5fit_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_scalar_5fit_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_scalar_5fit_2eproto)), true);
namespace scalar {

// ===================================================================

class scalar::_Internal {
 public:
  using HasBits = decltype(std::declval<scalar>()._has_bits_);
  static void set_has_iter(HasBits* has_bits) {
    (*has_bits)[0] |= 2u;
  }
  static void set_has_value(HasBits* has_bits) {
    (*has_bits)[0] |= 1u;
  }
  static bool MissingRequiredFields(const HasBits& has_bits) {
    return ((has_bits[0] & 0x00000003) ^ 0x00000003) != 0;
  }
};

scalar::scalar(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:scalar.scalar)
}
scalar::scalar(const scalar& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      _has_bits_(from._has_bits_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&value_, &from.value_,
    static_cast<size_t>(reinterpret_cast<char*>(&iter_) -
    reinterpret_cast<char*>(&value_)) + sizeof(iter_));
  // @@protoc_insertion_point(copy_constructor:scalar.scalar)
}

void scalar::SharedCtor() {
  ::memset(reinterpret_cast<char*>(this) + static_cast<size_t>(
      reinterpret_cast<char*>(&value_) - reinterpret_cast<char*>(this)),
      0, static_cast<size_t>(reinterpret_cast<char*>(&iter_) -
      reinterpret_cast<char*>(&value_)) + sizeof(iter_));
}

scalar::~scalar() {
  // @@protoc_insertion_point(destructor:scalar.scalar)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void scalar::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void scalar::ArenaDtor(void* object) {
  scalar* _this = reinterpret_cast< scalar* >(object);
  (void)_this;
}
void scalar::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void scalar::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const scalar& scalar::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_scalar_scalar_5fit_2eproto.base);
  return *internal_default_instance();
}


void scalar::Clear() {
// @@protoc_insertion_point(message_clear_start:scalar.scalar)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    ::memset(&value_, 0, static_cast<size_t>(
        reinterpret_cast<char*>(&iter_) -
        reinterpret_cast<char*>(&value_)) + sizeof(iter_));
  }
  _has_bits_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* scalar::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
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
      // required double value = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 17)) {
          _Internal::set_has_value(&has_bits);
          value_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
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

::PROTOBUF_NAMESPACE_ID::uint8* scalar::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:scalar.scalar)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = _has_bits_[0];
  // required int32 iter = 1;
  if (cached_has_bits & 0x00000002u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt32ToArray(1, this->_internal_iter(), target);
  }

  // required double value = 2;
  if (cached_has_bits & 0x00000001u) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(2, this->_internal_value(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:scalar.scalar)
  return target;
}

size_t scalar::RequiredFieldsByteSizeFallback() const {
// @@protoc_insertion_point(required_fields_byte_size_fallback_start:scalar.scalar)
  size_t total_size = 0;

  if (_internal_has_value()) {
    // required double value = 2;
    total_size += 1 + 8;
  }

  if (_internal_has_iter()) {
    // required int32 iter = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
        this->_internal_iter());
  }

  return total_size;
}
size_t scalar::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:scalar.scalar)
  size_t total_size = 0;

  if (((_has_bits_[0] & 0x00000003) ^ 0x00000003) == 0) {  // All required fields are present.
    // required double value = 2;
    total_size += 1 + 8;

    // required int32 iter = 1;
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
        this->_internal_iter());

  } else {
    total_size += RequiredFieldsByteSizeFallback();
  }
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void scalar::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:scalar.scalar)
  GOOGLE_DCHECK_NE(&from, this);
  const scalar* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<scalar>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:scalar.scalar)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:scalar.scalar)
    MergeFrom(*source);
  }
}

void scalar::MergeFrom(const scalar& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:scalar.scalar)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  cached_has_bits = from._has_bits_[0];
  if (cached_has_bits & 0x00000003u) {
    if (cached_has_bits & 0x00000001u) {
      value_ = from.value_;
    }
    if (cached_has_bits & 0x00000002u) {
      iter_ = from.iter_;
    }
    _has_bits_[0] |= cached_has_bits;
  }
}

void scalar::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:scalar.scalar)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void scalar::CopyFrom(const scalar& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:scalar.scalar)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool scalar::IsInitialized() const {
  if (_Internal::MissingRequiredFields(_has_bits_)) return false;
  return true;
}

void scalar::InternalSwap(scalar* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  swap(_has_bits_[0], other->_has_bits_[0]);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(scalar, iter_)
      + sizeof(scalar::iter_)
      - PROTOBUF_FIELD_OFFSET(scalar, value_)>(
          reinterpret_cast<char*>(&value_),
          reinterpret_cast<char*>(&other->value_));
}

::PROTOBUF_NAMESPACE_ID::Metadata scalar::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

class full_scalar_it::_Internal {
 public:
};

full_scalar_it::full_scalar_it(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  scalar_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:scalar.full_scalar_it)
}
full_scalar_it::full_scalar_it(const full_scalar_it& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      scalar_(from.scalar_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:scalar.full_scalar_it)
}

void full_scalar_it::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_full_scalar_it_scalar_5fit_2eproto.base);
}

full_scalar_it::~full_scalar_it() {
  // @@protoc_insertion_point(destructor:scalar.full_scalar_it)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void full_scalar_it::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void full_scalar_it::ArenaDtor(void* object) {
  full_scalar_it* _this = reinterpret_cast< full_scalar_it* >(object);
  (void)_this;
}
void full_scalar_it::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void full_scalar_it::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const full_scalar_it& full_scalar_it::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_full_scalar_it_scalar_5fit_2eproto.base);
  return *internal_default_instance();
}


void full_scalar_it::Clear() {
// @@protoc_insertion_point(message_clear_start:scalar.full_scalar_it)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  scalar_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* full_scalar_it::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated .scalar.scalar scalar = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_scalar(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<10>(ptr));
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

::PROTOBUF_NAMESPACE_ID::uint8* full_scalar_it::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:scalar.full_scalar_it)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated .scalar.scalar scalar = 1;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_scalar_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(1, this->_internal_scalar(i), target, stream);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:scalar.full_scalar_it)
  return target;
}

size_t full_scalar_it::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:scalar.full_scalar_it)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated .scalar.scalar scalar = 1;
  total_size += 1UL * this->_internal_scalar_size();
  for (const auto& msg : this->scalar_) {
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

void full_scalar_it::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:scalar.full_scalar_it)
  GOOGLE_DCHECK_NE(&from, this);
  const full_scalar_it* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<full_scalar_it>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:scalar.full_scalar_it)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:scalar.full_scalar_it)
    MergeFrom(*source);
  }
}

void full_scalar_it::MergeFrom(const full_scalar_it& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:scalar.full_scalar_it)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  scalar_.MergeFrom(from.scalar_);
}

void full_scalar_it::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:scalar.full_scalar_it)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void full_scalar_it::CopyFrom(const full_scalar_it& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:scalar.full_scalar_it)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool full_scalar_it::IsInitialized() const {
  if (!::PROTOBUF_NAMESPACE_ID::internal::AllAreInitialized(scalar_)) return false;
  return true;
}

void full_scalar_it::InternalSwap(full_scalar_it* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  scalar_.InternalSwap(&other->scalar_);
}

::PROTOBUF_NAMESPACE_ID::Metadata full_scalar_it::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace scalar
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::scalar::scalar* Arena::CreateMaybeMessage< ::scalar::scalar >(Arena* arena) {
  return Arena::CreateMessageInternal< ::scalar::scalar >(arena);
}
template<> PROTOBUF_NOINLINE ::scalar::full_scalar_it* Arena::CreateMaybeMessage< ::scalar::full_scalar_it >(Arena* arena) {
  return Arena::CreateMessageInternal< ::scalar::full_scalar_it >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
