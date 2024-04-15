// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_BACKENDS_ISTL_SINGLEROWCOLMATRIX_HH
#define DUNE_FUNCTIONS_BACKENDS_ISTL_SINGLEROWCOLMATRIX_HH

#include <cassert>
#include <type_traits>

#include <dune/common/ftraits.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/std/type_traits.hh>

namespace Dune { namespace Functions {
namespace Impl
{
  template <class M, class = void>
  struct SizeType
  {
    using type = std::size_t;
  };

  template <class M>
  struct SizeType<M,std::void_t<typename M::size_type>>
  {
    using type = typename M::size_type;
  };

  template <class M, class = void>
  struct BlockType
  {
    using type = void;
  };

  template <class M>
  struct BlockType<M,std::void_t<typename M::block_type>>
  {
    using type = typename M::block_type;
  };

} // end namespace Impl


//! Matrix-wrapper representing a matrix with a single column. Used in hierarchic block structure.
template <class Matrix>
class SingleColumnMatrix
    : public Matrix
{
  using Base = Matrix;

  static inline constexpr Dune::index_constant<0> zero_ = {};

  template <class M>
  using DynamicIndexAccessible = decltype(std::declval<M>()[0u][zero_]);

public:
  using size_type = typename Impl::SizeType<Base>::type;
  using block_type = typename Impl::BlockType<Base>::type;
  using field_type = typename Base::field_type;

private:
  template <size_type I0 = 0, class Mat>
  static constexpr auto row_range (const Mat& mat)
  {
    if constexpr (Dune::Std::is_detected_v<DynamicIndexAccessible, Mat>)
      return Dune::IntegralRange<size_type>{I0,Mat::N()};
    else
      return Dune::StaticIntegralRange<size_type,Mat::N(),I0>{};
  }

public:
  using Matrix::Matrix;
  using Matrix::operator=;

  const Matrix& matrix () const
  {
    return static_cast<const Matrix&>(*this);
  }

  Matrix& matrix ()
  {
    return static_cast<Matrix&>(*this);
  }

  //! Provide fixed single column resize by `resize(.,.)` method
  template <class M = Matrix,
    decltype(std::declval<M>().resize(0u,0u), bool{}) = true>
  void resize (std::size_t r, std::size_t)
  {
    matrix().resize(r, 1);
  }

  //! Provide fixed single column resize by `setSize(.,.)` method
  template <class M = Matrix,
    decltype(std::declval<M>().setSize(0u,0u), bool{}) = true>
  void setSize (std::size_t r, std::size_t)
  {
    matrix().setSize(r, 1);
  }

  //! y = A x
  template <class X, class Y>
  void mv (const X& x, Y& y) const
  {
    assert(this->N() == y.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].mv(x, y[i]);
    });
  }

  //! y = A^T x
  template <class X, class Y>
  void mtv (const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    assert(this->N() > 0);
    (*this)[zero_][zero_].mtv(x[zero_], y);
    Dune::Hybrid::forEach(row_range<1>(*this), [&](auto i) -> void {
      (*this)[i][zero_].umtv(x[i], y);
    });
  }

  //! y += A x
  template <class X, class Y>
  void umv (const X& x, Y& y) const
  {
    assert(this->N() == y.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].umv(x, y[i]);
    });
  }

  //! y -= A x
  template <class X, class Y>
  void mmv (const X& x, Y& y) const
  {
    assert(this->N() == y.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].mmv(x, y[i]);
    });
  }

  //! y += alpha A x
  template <class X, class Y>
  void usmv (const field_type& alpha, const X& x, Y& y) const
  {
    assert(this->N() == y.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].usmv(alpha, x, y[i]);
    });
  }

  //! y += A^T x
  template <class X, class Y>
  void umtv (const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].umtv(x[i], y);
    });
  }

  //! y -= A^T x
  template <class X, class Y>
  void mmtv (const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].mmtv(x[i], y);
    });
  }

  //! y += alpha A^T x
  template <class X, class Y>
  void usmtv (const field_type& alpha, const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].usmtv(alpha, x[i], y);
    });
  }

  //! y += A^H x
  template <class X, class Y>
  void umhv (const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].umhv(x[i], y);
    });
  }

  //! y -= A^H x
  template <class X, class Y>
  void mmhv (const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].mmhv(x[i], y);
    });
  }

  //! y += alpha A^H x
  template <class X, class Y>
  void usmhv (const field_type& alpha, const X& x, Y& y) const
  {
    assert(this->N() == x.size());
    Dune::Hybrid::forEach(row_range(*this), [&](auto i) -> void {
      (*this)[i][zero_].usmhv(alpha, x[i], y);
    });
  }
};

template <class M>
struct IsSingleColumnMatrix : std::false_type {};

template <class M>
struct IsSingleColumnMatrix<SingleColumnMatrix<M>> : std::true_type {};



//! Matrix-wrapper representing a matrix with a single row. Used in hierarchic block structure.
template <class Matrix>
class SingleRowMatrix
    : public Matrix
{
  using Base = Matrix;

  static inline constexpr Dune::index_constant<0> zero_ = {};

  template <class M>
  using DynamicIndexAccessible = decltype(std::declval<M>()[zero_][0u]);

public:
  using block_type = typename Impl::BlockType<Base>::type;
  using size_type = typename Impl::SizeType<Base>::type;
  using field_type = typename Base::field_type;

private:
  template <size_type I0 = 0, class Mat>
  static constexpr auto col_range (const Mat& mat)
  {
    if constexpr (Dune::Std::is_detected_v<DynamicIndexAccessible,Mat>)
      return Dune::IntegralRange<size_type>{I0,Mat::M()};
    else
      return Dune::StaticIntegralRange<size_type,Mat::M(),I0>{};
  }

public:
  using Matrix::Matrix;
  using Matrix::operator=;

  const Matrix& matrix () const
  {
    return static_cast<const Matrix&>(*this);
  }

  Matrix& matrix ()
  {
    return static_cast<Matrix&>(*this);
  }

  template <class M = Matrix,
    decltype(std::declval<M>().resize(0u,0u), bool{}) = true>
  void resize (std::size_t, std::size_t c)
  {
    matrix().resize(1, c);
  }

  template <class M = Matrix,
    decltype(std::declval<M>().setSize(0u,0u), bool{}) = true>
  void setSize (std::size_t, std::size_t c)
  {
    matrix().setSize(1, c);
  }

  //! y = A x
  template <class X, class Y>
  void mv (const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    (*this)[zero_][zero_].mv(x[zero_], y);
    Dune::Hybrid::forEach(col_range<1>(*this), [&](auto j) -> void {
      (*this)[zero_][j].umv(x[j], y);
    });
  }

  //! y = A^T x
  template <class X, class Y>
  void mtv (const X& x, Y& y) const
  {
    assert(this->M() == y.size());
    assert(this->M() > 0);
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].mtv(x, y[j]);
    });
  }

  //! y += A x
  template <class X, class Y>
  void umv(const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].umv(x[j], y);
    });
  }

  //! y -= A x
  template <class X, class Y>
  void mmv (const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].mmv(x[j], y);
    });
  }

  //! y += alpha A x
  template <class X, class Y>
  void usmv (const field_type& alpha, const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].usmv(alpha, x[j], y);
    });
  }

  //! y += A^T x
  template <class X, class Y>
  void umtv (const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].umv(x[j], y);
    });
  }

  //! y -= A^T x
  template <class X, class Y>
  void mmtv (const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].mmtv(x[j], y);
    });
  }

  //! y += alpha A^T x
  template <class X, class Y>
  void usmtv (const field_type& alpha, const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].usmtv(alpha, x[j], y);
    });
  }

  //! y += A^H x
  template <class X, class Y>
  void umhv (const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].umhv(x[j], y);
    });
  }

  //! y -= A^H x
  template <class X, class Y>
  void mmhv (const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].mmhv(x[j], y);
    });
  }

  //! y += alpha A^H x
  template <class X, class Y>
  void usmhv (const field_type& alpha, const X& x, Y& y) const
  {
    assert(this->M() == x.size());
    Dune::Hybrid::forEach(col_range(*this), [&](auto j) -> void {
      (*this)[zero_][j].usmtv(alpha, x[j], y);
    });
  }
};

template <class M>
struct IsSingleRowMatrix : std::false_type {};

template <class M>
struct IsSingleRowMatrix<SingleRowMatrix<M>> : std::true_type {};

} // end namespace Functions

template <class Matrix>
class FieldTraits<Functions::SingleColumnMatrix<Matrix>>
    : public FieldTraits<Matrix>
{};

template <class Matrix>
class FieldTraits<Functions::SingleRowMatrix<Matrix>>
    : public FieldTraits<Matrix>
{};

} // end namespace Dune

#endif // DUNE_FUNCTIONS_BACKENDS_ISTL_SINGLEROWCOLMATRIX_HH
