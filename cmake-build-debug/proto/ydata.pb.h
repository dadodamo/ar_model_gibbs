// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: ydata.proto

#ifndef GOOGLE_PROTOBUF_INCLUDED_ydata_2eproto
#define GOOGLE_PROTOBUF_INCLUDED_ydata_2eproto

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
#define PROTOBUF_INTERNAL_EXPORT_ydata_2eproto
PROTOBUF_NAMESPACE_OPEN
namespace internal {
class AnyMetadata;
}  // namespace internal
PROTOBUF_NAMESPACE_CLOSE

// Internal implementation detail -- do not use these members.
struct TableStruct_ydata_2eproto {
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
extern const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_ydata_2eproto;
namespace y_data {
class full_y;
class full_yDefaultTypeInternal;
extern full_yDefaultTypeInternal _full_y_default_instance_;
class vector;
class vectorDefaultTypeInternal;
extern vectorDefaultTypeInternal _vector_default_instance_;
}  // namespace y_data
PROTOBUF_NAMESPACE_OPEN
template<> ::y_data::full_y* Arena::CreateMaybeMessage<::y_data::full_y>(Arena*);
template<> ::y_data::vector* Arena::CreateMaybeMessage<::y_data::vector>(Arena*);
PROTOBUF_NAMESPACE_CLOSE
namespace y_data {

// ===================================================================

class vector PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:y_data.vector) */ {
 public:
  inline vector() : vector(nullptr) {}
  virtual ~vector();

  vector(const vector& from);
  vector(vector&& from) noexcept
    : vector() {
    *this = ::std::move(from);
  }

  inline vector& operator=(const vector& from) {
    CopyFrom(from);
    return *this;
  }
  inline vector& operator=(vector&& from) noexcept {
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
  static const vector& default_instance();

  static inline const vector* internal_default_instance() {
    return reinterpret_cast<const vector*>(
               &_vector_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    0;

  friend void swap(vector& a, vector& b) {
    a.Swap(&b);
  }
  inline void Swap(vector* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(vector* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline vector* New() const final {
    return CreateMaybeMessage<vector>(nullptr);
  }

  vector* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<vector>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const vector& from);
  void MergeFrom(const vector& from);
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
  void InternalSwap(vector* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "y_data.vector";
  }
  protected:
  explicit vector(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_ydata_2eproto);
    return ::descriptor_table_ydata_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kVecValueFieldNumber = 2,
  };
  // repeated double vec_value = 2;
  int vec_value_size() const;
  private:
  int _internal_vec_value_size() const;
  public:
  void clear_vec_value();
  private:
  double _internal_vec_value(int index) const;
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      _internal_vec_value() const;
  void _internal_add_vec_value(double value);
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      _internal_mutable_vec_value();
  public:
  double vec_value(int index) const;
  void set_vec_value(int index, double value);
  void add_vec_value(double value);
  const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
      vec_value() const;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
      mutable_vec_value();

  // @@protoc_insertion_point(class_scope:y_data.vector)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedField< double > vec_value_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  friend struct ::TableStruct_ydata_2eproto;
};
// -------------------------------------------------------------------

class full_y PROTOBUF_FINAL :
    public ::PROTOBUF_NAMESPACE_ID::Message /* @@protoc_insertion_point(class_definition:y_data.full_y) */ {
 public:
  inline full_y() : full_y(nullptr) {}
  virtual ~full_y();

  full_y(const full_y& from);
  full_y(full_y&& from) noexcept
    : full_y() {
    *this = ::std::move(from);
  }

  inline full_y& operator=(const full_y& from) {
    CopyFrom(from);
    return *this;
  }
  inline full_y& operator=(full_y&& from) noexcept {
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
  static const full_y& default_instance();

  static inline const full_y* internal_default_instance() {
    return reinterpret_cast<const full_y*>(
               &_full_y_default_instance_);
  }
  static constexpr int kIndexInFileMessages =
    1;

  friend void swap(full_y& a, full_y& b) {
    a.Swap(&b);
  }
  inline void Swap(full_y* other) {
    if (other == this) return;
    if (GetArena() == other->GetArena()) {
      InternalSwap(other);
    } else {
      ::PROTOBUF_NAMESPACE_ID::internal::GenericSwap(this, other);
    }
  }
  void UnsafeArenaSwap(full_y* other) {
    if (other == this) return;
    GOOGLE_DCHECK(GetArena() == other->GetArena());
    InternalSwap(other);
  }

  // implements Message ----------------------------------------------

  inline full_y* New() const final {
    return CreateMaybeMessage<full_y>(nullptr);
  }

  full_y* New(::PROTOBUF_NAMESPACE_ID::Arena* arena) const final {
    return CreateMaybeMessage<full_y>(arena);
  }
  void CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) final;
  void CopyFrom(const full_y& from);
  void MergeFrom(const full_y& from);
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
  void InternalSwap(full_y* other);
  friend class ::PROTOBUF_NAMESPACE_ID::internal::AnyMetadata;
  static ::PROTOBUF_NAMESPACE_ID::StringPiece FullMessageName() {
    return "y_data.full_y";
  }
  protected:
  explicit full_y(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  private:
  static void ArenaDtor(void* object);
  inline void RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena* arena);
  public:

  ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadata() const final;
  private:
  static ::PROTOBUF_NAMESPACE_ID::Metadata GetMetadataStatic() {
    ::PROTOBUF_NAMESPACE_ID::internal::AssignDescriptors(&::descriptor_table_ydata_2eproto);
    return ::descriptor_table_ydata_2eproto.file_level_metadata[kIndexInFileMessages];
  }

  public:

  // nested types ----------------------------------------------------

  // accessors -------------------------------------------------------

  enum : int {
    kVecTFieldNumber = 1,
  };
  // repeated .y_data.vector vec_t = 1;
  int vec_t_size() const;
  private:
  int _internal_vec_t_size() const;
  public:
  void clear_vec_t();
  ::y_data::vector* mutable_vec_t(int index);
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::y_data::vector >*
      mutable_vec_t();
  private:
  const ::y_data::vector& _internal_vec_t(int index) const;
  ::y_data::vector* _internal_add_vec_t();
  public:
  const ::y_data::vector& vec_t(int index) const;
  ::y_data::vector* add_vec_t();
  const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::y_data::vector >&
      vec_t() const;

  // @@protoc_insertion_point(class_scope:y_data.full_y)
 private:
  class _Internal;

  template <typename T> friend class ::PROTOBUF_NAMESPACE_ID::Arena::InternalHelper;
  typedef void InternalArenaConstructable_;
  typedef void DestructorSkippable_;
  ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::y_data::vector > vec_t_;
  mutable ::PROTOBUF_NAMESPACE_ID::internal::CachedSize _cached_size_;
  friend struct ::TableStruct_ydata_2eproto;
};
// ===================================================================


// ===================================================================

#ifdef __GNUC__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif  // __GNUC__
// vector

// repeated double vec_value = 2;
inline int vector::_internal_vec_value_size() const {
  return vec_value_.size();
}
inline int vector::vec_value_size() const {
  return _internal_vec_value_size();
}
inline void vector::clear_vec_value() {
  vec_value_.Clear();
}
inline double vector::_internal_vec_value(int index) const {
  return vec_value_.Get(index);
}
inline double vector::vec_value(int index) const {
  // @@protoc_insertion_point(field_get:y_data.vector.vec_value)
  return _internal_vec_value(index);
}
inline void vector::set_vec_value(int index, double value) {
  vec_value_.Set(index, value);
  // @@protoc_insertion_point(field_set:y_data.vector.vec_value)
}
inline void vector::_internal_add_vec_value(double value) {
  vec_value_.Add(value);
}
inline void vector::add_vec_value(double value) {
  _internal_add_vec_value(value);
  // @@protoc_insertion_point(field_add:y_data.vector.vec_value)
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
vector::_internal_vec_value() const {
  return vec_value_;
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >&
vector::vec_value() const {
  // @@protoc_insertion_point(field_list:y_data.vector.vec_value)
  return _internal_vec_value();
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
vector::_internal_mutable_vec_value() {
  return &vec_value_;
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedField< double >*
vector::mutable_vec_value() {
  // @@protoc_insertion_point(field_mutable_list:y_data.vector.vec_value)
  return _internal_mutable_vec_value();
}

// -------------------------------------------------------------------

// full_y

// repeated .y_data.vector vec_t = 1;
inline int full_y::_internal_vec_t_size() const {
  return vec_t_.size();
}
inline int full_y::vec_t_size() const {
  return _internal_vec_t_size();
}
inline void full_y::clear_vec_t() {
  vec_t_.Clear();
}
inline ::y_data::vector* full_y::mutable_vec_t(int index) {
  // @@protoc_insertion_point(field_mutable:y_data.full_y.vec_t)
  return vec_t_.Mutable(index);
}
inline ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::y_data::vector >*
full_y::mutable_vec_t() {
  // @@protoc_insertion_point(field_mutable_list:y_data.full_y.vec_t)
  return &vec_t_;
}
inline const ::y_data::vector& full_y::_internal_vec_t(int index) const {
  return vec_t_.Get(index);
}
inline const ::y_data::vector& full_y::vec_t(int index) const {
  // @@protoc_insertion_point(field_get:y_data.full_y.vec_t)
  return _internal_vec_t(index);
}
inline ::y_data::vector* full_y::_internal_add_vec_t() {
  return vec_t_.Add();
}
inline ::y_data::vector* full_y::add_vec_t() {
  // @@protoc_insertion_point(field_add:y_data.full_y.vec_t)
  return _internal_add_vec_t();
}
inline const ::PROTOBUF_NAMESPACE_ID::RepeatedPtrField< ::y_data::vector >&
full_y::vec_t() const {
  // @@protoc_insertion_point(field_list:y_data.full_y.vec_t)
  return vec_t_;
}

#ifdef __GNUC__
  #pragma GCC diagnostic pop
#endif  // __GNUC__
// -------------------------------------------------------------------


// @@protoc_insertion_point(namespace_scope)

}  // namespace y_data

// @@protoc_insertion_point(global_scope)

#include <google/protobuf/port_undef.inc>
#endif  // GOOGLE_PROTOBUF_INCLUDED_GOOGLE_PROTOBUF_INCLUDED_ydata_2eproto
