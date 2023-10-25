// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: scalar_it.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_scalar_5fit_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_scalar_5fit_2eproto

#include <limits>
#include <string>

#include <google/protobuf/port_def.inc>
#if PROTOBUF_VERSION < 3014000
#error This file was generated by a newer version of protoc which is
#error incompatible with your Protocol Buffer headers. Please update
#error your headers.
#endif
#if 3014000 < PROTOBUF_MIN_PROTOC_VERSION
#error This file was generated by an older version of protoc which is
#error incompatible with your Protocol Buffer headers. Please
#error regenerate this file with a newer version of protoc.
#endif

#include <google/protobuf/port_undef.inc>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/arena.h>
#include <google/protobuf/arenastring.h>
#include <google/protobuf/generated_message_table_driven.h>
#include <google/protobuf/generated_message_util.h>
#include <google/protobuf/metadata_lite.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>  // IWYU pragma: export
#include <google/protobuf/extension_set.h>  // IWYU pragma: export
#include <google/protobuf/unknown_field_set.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
#define PROTOBUF_INTERNAL_EXPORT_scalar_5fit_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_scalar_5fit_2eproto {
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTableField entries[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::AuxiliaryParseTableField aux[]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::ParseTable schema[2]
    PROTOBUF_SECTION_VARIABLE(protodesc_cold);
  static const ::PROTOBUF_NAMESPACE_ID::internal::FieldMetadata field_metadata[];
  static const ::PROTOBUF_NAMESPACE_ID::internal::SerializationTable serialization_table[];
  static const ::PROTOBUF_NAMESPACE_ID::uint32 offsets[];
};
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_scalar_5fit_2eproto;
namespace scalar {
class full_scalar_it;
class full_scalar_itDefaultTypeInternal;
extern full_scalar_itDefaultTypeInternal _full_scalar_it_default_instance_;
class scalar;
class scalarDefaultTypeInternal;
extern scalarDefaultTypeInternal _scalar_default_instance_;
}  // namespace scalar
PROTOBUF_NAMESPACE_OPEN
template<> ::scalar::full_scalar_it* Arena::CreateMaybeMessage<::scalar::full_scalar_it>(Arena*);
template<> ::scalar::scalar* Arena::CreateMaybeMessage<::scalar::scalar>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace scalar {

// ===================================================================

class scalar PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:scalar.scalar) */ {
 public:
  inline scalar() : scalar(nullptr) {}
  virtual ~scalar();

  scalar(const scalar& from);
  scalar(scalar&& from) noexcept
    : scalar() {
    *this = ::std::move(from);
  }

  inline scalar& operator=(const scalar& from) {
    CopyFrom(from);
    return *this;
  }
  inline scalar& operator=(scalar&& from) noexcept {
    if (GetArena() == from.GetArena()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }

  inline const ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance);
  }
  inline ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
  }

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return GetMetadataStatic().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return GetMetadataStatic().reflection;
  }
  static const scalar& default_instance();

  static inline const scalar* internal_default_instance() {
    return reinterpret_cast<const scalar*>(
               &_scalar_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(scalar& a, scalar& b) {
    a.Swap(&b);
  }
  inline void Swap(scalar* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(scalar* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline scalar* New() const final {
    return CreateMaybeMessage<scalar>(nullptr);
  }

  scalar* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<scalar>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const scalar& from);
  void MergeFrom(const scalar& from);
  PROTOBUF_ATTRIBUTE_REINITIALIZES void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  const char* _InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) final;
  ::PROTOBUF_NAMESPACE_ID::uint8* _InternalSerialize(
      ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  inline void SharedCtor();
  inline void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(scalar* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "scalar.scalar";
  }
  protected:
  explicit scalar(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_scalar_5fit_2eproto);
    return ::descriptor_table_scalar_5fit_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kValueFieldNumber = 2,
    kIterFieldNumber = 1,
  };
  // required double value = 2;
  bool has_value() const;
  private:
  bool _internal_has_value() const;
  public:
  void clear_value();
  double value() const;
  void set_value(double value);
  private:
  double _internal_value() const;
  void _internal_set_value(double value);
  public:

  // required int32 iter = 1;
  bool has_iter() const;
  private:
  bool _internal_has_iter() const;
  public:
  void clear_iter();
  ::PROTOBUF_NAMESPACE_ID::int32 iter() const;
  void set_iter(::PROTOBUF_NAMESPACE_ID::int32 value);
  private:
  ::PROTOBUF_NAMESPACE_ID::int32 _internal_iter() const;
  void _internal_set_iter(::PROTOBUF_NAMESPACE_ID::int32 value);
  public:

  // @@protoc_insertion_point(class_scope:scalar.scalar)
 private:
  class _Internal;

  // helper for ByteSizeLong()
  size_t RequiredFieldsByteSizeFallback() const;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::internal::HasBits<1> _has_bits_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  double value_;
  ::PROTOBUF_NAMESPACE_ID::int32 iter_;
  friend struct ::TableStruct_scalar_5fit_2eproto;
};
// -------------------------------------------------------------------

class full_scalar_it PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:scalar.full_scalar_it) */ {
 public:
  inline full_scalar_it() : full_scalar_it(nullptr) {}
  virtual ~full_scalar_it();

  full_scalar_it(const full_scalar_it& from);
  full_scalar_it(full_scalar_it&& from) noexcept
    : full_scalar_it() {
    *this = ::std::move(from);
  }

  inline full_scalar_it& operator=(const full_scalar_it& from) {
    CopyFrom(from);
    return *this;
  }
  inline full_scalar_it& operator=(full_scalar_it&& from) noexcept {
    if (GetArena() == from.GetArena()) {
      if (this != &from) InternalSwap(&from);
    } else {
      CopyFrom(from);
    }
    return *this;
  }

  inline const ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet& unknown_fields() const {
    return _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance);
  }
  inline ::PROTOBUF_NAMESPACE_ID::UnknownFieldSet* mutable_unknown_fields() {
    return _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
  }

  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* descriptor() {
    return GetDescriptor();
  }
  static const ::PROTOBUF_NAMESPACE_ID::Descriptor* GetDescriptor() {
    return GetMetadataStatic().descriptor;
  }
  static const ::PROTOBUF_NAMESPACE_ID::Reflection* GetReflection() {
    return GetMetadataStatic().reflection;
  }
  static const full_scalar_it& default_instance();

  static inline const full_scalar_it* internal_default_instance() {
    return reinterpret_cast<const full_scalar_it*>(
               &_full_scalar_it_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    1;

  friend void swap(full_scalar_it& a, full_scalar_it& b) {
    a.Swap(&b);
  }
  inline void Swap(full_scalar_it* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(full_scalar_it* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline full_scalar_it* New() const final {
    return CreateMaybeMessage<full_scalar_it>(nullptr);
  }

  full_scalar_it* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<full_scalar_it>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const full_scalar_it& from);
  void MergeFrom(const full_scalar_it& from);
  PROTOBUF_ATTRIBUTE_REINITIALIZES void Clear() final;
  bool IsInitialized() const final;

  size_t ByteSizeLong() const final;
  const char* _InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) final;
  ::PROTOBUF_NAMESPACE_ID::uint8* _InternalSerialize(
      ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const final;
  int GetCachedSize() const final { return _cached_size_.Get(); }

  private:
  inline void SharedCtor();
  inline void SharedDtor();
  void SetCachedSize(int size) const final;
  void InternalSwap(full_scalar_it* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "scalar.full_scalar_it";
  }
  protected:
  explicit full_scalar_it(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_scalar_5fit_2eproto);
    return ::descriptor_table_scalar_5fit_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kScalarFieldNumber = 1,
  };
  // repeated .scalar.scalar scalar = 1;
  int scalar_size() const;
  private:
  int _internal_scalar_size() const;
  public:
  void clear_scalar();
  ::scalar::scalar* mutable_scalar(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::scalar::scalar >*
      mutable_scalar();
  private:
  const ::scalar::scalar& _internal_scalar(int index) const;
  ::scalar::scalar* _internal_add_scalar();
  public:
  const ::scalar::scalar& scalar(int index) const;
  ::scalar::scalar* add_scalar();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::scalar::scalar >&
      scalar() const;

  // @@protoc_insertion_point(class_scope:scalar.full_scalar_it)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::scalar::scalar > scalar_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  friend struct ::TableStruct_scalar_5fit_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// scalar

// required int32 iter = 1;
inline bool scalar::_internal_has_iter() const {
  bool value = (_has_bits_[0] & 0x00000002u) != 0;
  return value;
}
inline bool scalar::has_iter() const {
  return _internal_has_iter();
}
inline void scalar::clear_iter() {
  iter_ = 0;
  _has_bits_[0] &= ~0x00000002u;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 scalar::_internal_iter() const {
  return iter_;
}
inline ::PROTOBUF_NAMESPACE_ID::int32 scalar::iter() const {
  // @@protoc_insertion_point(field_get:scalar.scalar.iter)
  return _internal_iter();
}
inline void scalar::_internal_set_iter(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _has_bits_[0] |= 0x00000002u;
  iter_ = value;
}
inline void scalar::set_iter(::PROTOBUF_NAMESPACE_ID::int32 value) {
  _internal_set_iter(value);
  // @@protoc_insertion_point(field_set:scalar.scalar.iter)
}

// required double value = 2;
inline bool scalar::_internal_has_value() const {
  bool value = (_has_bits_[0] & 0x00000001u) != 0;
  return value;
}
inline bool scalar::has_value() const {
  return _internal_has_value();
}
inline void scalar::clear_value() {
  value_ = 0;
  _has_bits_[0] &= ~0x00000001u;
}
inline double scalar::_internal_value() const {
  return value_;
}
inline double scalar::value() const {
  // @@protoc_insertion_point(field_get:scalar.scalar.value)
  return _internal_value();
}
inline void scalar::_internal_set_value(double value) {
  _has_bits_[0] |= 0x00000001u;
  value_ = value;
}
inline void scalar::set_value(double value) {
  _internal_set_value(value);
  // @@protoc_insertion_point(field_set:scalar.scalar.value)
}

// -------------------------------------------------------------------

// full_scalar_it

// repeated .scalar.scalar scalar = 1;
inline int full_scalar_it::_internal_scalar_size() const {
  return scalar_.size();
}
inline int full_scalar_it::scalar_size() const {
  return _internal_scalar_size();
}
inline void full_scalar_it::clear_scalar() {
  scalar_.Clear();
}
inline ::scalar::scalar* full_scalar_it::mutable_scalar(int index) {
  // @@protoc_insertion_point(field_mutable:scalar.full_scalar_it.scalar)
  return scalar_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::scalar::scalar >*
full_scalar_it::mutable_scalar() {
  // @@protoc_insertion_point(field_mutable_list:scalar.full_scalar_it.scalar)
  return &scalar_;
}
inline const ::scalar::scalar& full_scalar_it::_internal_scalar(int index) const {
  return scalar_.Get(index);
}
inline const ::scalar::scalar& full_scalar_it::scalar(int index) const {
  // @@protoc_insertion_point(field_get:scalar.full_scalar_it.scalar)
  return _internal_scalar(index);
}
inline ::scalar::scalar* full_scalar_it::_internal_add_scalar() {
  return scalar_.Add();
}
inline ::scalar::scalar* full_scalar_it::add_scalar() {
  // @@protoc_insertion_point(field_add:scalar.full_scalar_it.scalar)
  return _internal_add_scalar();
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::scalar::scalar >&
full_scalar_it::scalar() const {
  // @@protoc_insertion_point(field_list:scalar.full_scalar_it.scalar)
  return scalar_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace scalar

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_scalar_5fit_2eproto
