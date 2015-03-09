// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_SMALLOBJECT_HH
#define DUNE_FUNCTIONS_COMMON_SMALLOBJECT_HH

#include <utility>

namespace Dune {
namespace Functions {


/**
 * \brief A wrapper providing small object optimization
 *
 * \tparam Base Base class type of wrapped objects
 * \tparam bufferSize Size of small object buffer
 *
 * This class encapsulates small object optimization for polymorphic types.
 * The type of objects passed to the constructor must be derived from the
 * base class type Base.
 *
 * If the size of the derived type fits into the static buffer, then the
 * wrapped object is stored there, otherwise it is allocated dynamically.
 *
 * In order to make the copy constructor work for polymorphic types,
 * Base must provide clone() and clone(void*). The former should return
 * a pointer to a dynamically allocated clone, while the latter
 * should call the appropriate placement-new with the passed pointer.
 *
 * Similarly the polymorphic type has to implement move(void*).
 * This should call placement-new and can std::move all the
 * data but leave the object in a valid and probably unusable state.
 */
template<class Base, size_t bufferSize>
class SmallObject
{
public:

  template<class Derived>
  SmallObject(Derived&& derived)
  {
    if (sizeof(Derived)<bufferSize)
      p_ = new (buffer_) Derived(std::forward<Derived>(derived));
    else
      p_ = new Derived(std::forward<Derived>(derived));
  }

  SmallObject(SmallObject&& other)
  {
    moveToWrappedObject(std::move(other));
  }

  SmallObject(const SmallObject& other)
  {
    copyToWrappedObject(other);
  }

  ~SmallObject()
  {
    destroyWrappedObject();
  }

  SmallObject& operator=(const SmallObject& other)
  {
    destroyWrappedObject();
    copyToWrappedObject(other);
    return *this;
  }

  SmallObject& operator=(SmallObject&& other)
  {
    destroyWrappedObject();
    moveToWrappedObject(std::move(other));
    return *this;
  }

  bool bufferUsed() const
  {
    return ((void*) (p_) == (void*)(&buffer_));
  }

  const Base& get() const
  {
    return *p_;
  }

  Base& get()
  {
    return *p_;
  }

private:

  void destroyWrappedObject()
  {
    if (bufferUsed())
      p_->~Base();
    else
      delete p_;
  }

  void moveToWrappedObject(SmallObject&& other)
  {
    if (other.bufferUsed())
      p_ = other.p_->move(buffer_);
    else
    {
      // We don't need to check for &other_!=this, because you can't
      // have an rvalue to *this and call it's assignment/constructor
      // at the same time. (Despite trying to shot yourself in the foot
      // with std::move explicitly.)

      // Take ownership of allocated object
      p_ = other.p_;

      // Leave pointer in a clear state to avoid double freeing it.
      other.p_ = 0;
    }
  }

  void copyToWrappedObject(const SmallObject& other)
  {
    if (&other!=this)
    {
      if (other.bufferUsed())
        p_ = other.p_->clone(buffer_);
      else
        p_ = other.p_->clone();
    }
  }

  Base* p_;
  char buffer_[bufferSize];
};


} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_SMALLOBJECT_HH
