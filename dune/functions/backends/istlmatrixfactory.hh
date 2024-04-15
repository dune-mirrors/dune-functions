// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_BACKENDS_ISTL_MATRIXFACTORY_HH
#define DUNE_FUNCTIONS_BACKENDS_ISTL_MATRIXFACTORY_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/indices.hh>

#include <dune/functions/backends/singlerowcolmatrix.hh>
#include <dune/functions/functionspacebases/containerdescriptors.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>

namespace Dune::Functions {
namespace ContainerDescriptors {

template<class SparseFactory>
struct ISTLMatrixFactory
    : public SparseFactory
{
  void operator() (const Unknown&, const Unknown&) const
  {
    DUNE_THROW(Dune::NotImplemented, "Cannot create a matrix. The container descriptor is unknown.");
  }

  template <class CD>
  void operator() (const Unknown&, const CD&) const
  {
    DUNE_THROW(Dune::NotImplemented, "Cannot create a matrix. The container descriptor is unknown.");
  }

  template <class CD>
  void operator() (const CD&, const Unknown&) const
  {
    DUNE_THROW(Dune::NotImplemented, "Cannot create a matrix. The container descriptor is unknown.");
  }

private: // specializations for row or column descriptor it not type-uniform

  template<std::size_t m, std::size_t n, class Row, class Col>
  auto make_multitype_matrix (const Row& row, const Col& col) const
  {
    auto make_row = [&](auto i) {
      return unpackIntegerSequence([&](auto... jj) {
        using RowType = Dune::MultiTypeBlockVector<decltype((*this)(row[i],col[jj]))...>;
        return RowType{(*this)(row[i],col[jj])...};
      }, std::make_index_sequence<n>());
    };

    return unpackIntegerSequence([&](auto... ii) {
      using MatrixType = Dune::MultiTypeBlockMatrix<decltype(make_row(ii))...>;
      return MatrixType{make_row(ii)...};
    }, std::make_index_sequence<m>());
  }

public:

  template<class... V, class... W>
  auto operator() (const Tuple<V...>& row, const Tuple<W...>& col) const
  {
    return make_multitype_matrix<sizeof...(V),sizeof...(W)>(row,col);
  }

  template<class... V, class W, std::size_t n>
  auto operator() (const Tuple<V...>& row, const Array<W,n>& col) const
  {
    return make_multitype_matrix<sizeof...(V),n>(row,col);
  }

  template<class V, std::size_t m, class... W>
  auto operator() (const Array<V,m>& row, const Tuple<W...>& col) const
  {
    return make_multitype_matrix<m,sizeof...(W)>(row,col);
  }

  template<class... V, class W, std::size_t n>
  auto operator() (const Tuple<V...>& row, const UniformArray<W,n>& col) const
  {
    return make_multitype_matrix<sizeof...(V),n>(row,col);
  }

  template<class V, std::size_t m, class... W>
  auto operator() (const UniformArray<V,m>& row, const Tuple<W...>& col) const
  {
    return make_multitype_matrix<m,sizeof...(W)>(row,col);
  }

private:

  template<std::size_t m, class Row, class Flat>
  auto make_multitype_singlecol_matrix (const Row& row, const Flat& col) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      using ColMatrixType = Dune::MultiTypeBlockMatrix<
        Dune::MultiTypeBlockVector<decltype((*this)(row[ii],col))>...>;
      using MatrixType = Dune::Functions::SingleColumnMatrix<ColMatrixType>;
      return MatrixType{
        Dune::MultiTypeBlockVector<decltype((*this)(row[ii],col))>{(*this)(row[ii],col)}...
      };
    }, std::make_index_sequence<m>());
  }

public:

  template<class... V, std::size_t n>
  auto operator() (const Tuple<V...>& row, const UniformArray<Value,n>& col) const
  {
    make_multitype_singlecol_matrix<sizeof...(V)>(row,col);
  }

  template<class... V>
  auto operator() (const Tuple<V...>& row, const UniformVector<Value>& col) const
  {
    make_multitype_singlecol_matrix<sizeof...(V)>(row,col);
  }

private:

  template<std::size_t n, class Flat, class Col>
  auto make_multitype_singlerow_matrix (const Flat& row, const Col& col) const
  {
    return unpackIntegerSequence([&](auto... ii) {
      using RowType = Dune::MultiTypeBlockVector<decltype((*this)(row,col[ii]))...>;
      using RowMatrixType = Dune::MultiTypeBlockMatrix<RowType>;
      using MatrixType = Dune::Functions::SingleRowMatrix<RowMatrixType>;
      return MatrixType{RowType{(*this)(row,col[ii])...}};
    }, std::make_index_sequence<n>());
  }

public:

  template<std::size_t m, class... W>
  auto operator() (const UniformArray<Value,m>& row, const Tuple<W...>& col) const
  {
    return make_multitype_singlerow_matrix<sizeof...(W)>(row,col);
  }

  template<class... W>
  auto operator() (const UniformVector<Value>& row, const Tuple<W...>& col) const
  {
    return make_multitype_singlerow_matrix<sizeof...(W)>(row,col);
  }

private: // specializations for row or column type-uniform descriptors

  template<class Row, class Col>
  auto make_matrix (const Row& row, const Col& col) const
  {
    using VW = decltype((*this)(row[0u],col[0u]));
    Dune::Matrix<VW> mat(row.size(),col.size());
    for (std::size_t i = 0; i < row.size(); ++i)
      for (std::size_t j = 0; j < col.size(); ++j)
        mat[i][j] = (*this)(row[i],col[j]);
    return mat;
  }

public:

  template<class V, std::size_t m, class W, std::size_t n>
  auto operator() (const Array<V,m>& row, const Array<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W>
  auto operator() (const Array<V,m>& row, const Vector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W, std::size_t n>
  auto operator() (const Array<V,m>& row, const UniformArray<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W>
  auto operator() (const Array<V,m>& row, const UniformVector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W, std::size_t n>
  auto operator() (const Vector<V>& row, const Array<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W>
  auto operator() (const Vector<V>& row, const Vector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W, std::size_t n>
  auto operator() (const Vector<V>& row, const UniformArray<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W>
  auto operator() (const Vector<V>& row, const UniformVector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W, std::size_t n>
  auto operator() (const UniformArray<V,m>& row, const Array<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W>
  auto operator() (const UniformArray<V,m>& row, const Vector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W, std::size_t n>
  auto operator() (const UniformArray<V,m>& row, const UniformArray<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, std::size_t m, class W>
  auto operator() (const UniformArray<V,m>& row, const UniformVector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W, std::size_t n>
  auto operator() (const UniformVector<V>& row, const Array<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W>
  auto operator() (const UniformVector<V>& row, const Vector<W>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W, std::size_t n>
  auto operator() (const UniformVector<V>& row, const UniformArray<W,n>& col) const
  {
    return make_matrix(row,col);
  }

  template<class V, class W>
  auto operator() (const UniformVector<V>& row, const UniformVector<W>& col) const
  {
    return make_matrix(row,col);
  }

private:

  template<class Row, class Flat>
  auto make_singlecol_matrix (const Row& row, const Flat& col) const
  {
    using ColMatrixType = Dune::Matrix<decltype((*this)(row[0u],col))>;
    using MatrixType = Dune::Functions::SingleColumnMatrix<ColMatrixType>;
    MatrixType mat(row.size(),1);
    for (std::size_t i = 0; i < row.size(); ++i)
      mat[i][0] = (*this)(row[0u],col);
    return mat;
  }

public:

  template<class V, std::size_t m, std::size_t n>
  auto operator() (const Array<V,m>& row, const UniformArray<Value,n>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V, std::size_t m>
  auto operator() (const Array<V,m>& row, const UniformVector<Value>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V, std::size_t n>
  auto operator() (const Vector<V>& row, const UniformArray<Value,n>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V>
  auto operator() (const Vector<V>& row, const UniformVector<Value>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V, std::size_t m, std::size_t n>
  auto operator() (const UniformArray<V,m>& row, const UniformArray<Value,n>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V, std::size_t m>
  auto operator() (const UniformArray<V,m>& row, const UniformVector<Value>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V, std::size_t n>
  auto operator() (const UniformVector<V>& row, const UniformArray<Value,n>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

  template<class V>
  auto operator() (const UniformVector<V>& row, const UniformVector<Value>& col) const
  {
    return make_singlecol_matrix(row,col);
  }

private:

  template<class Flat, class Col>
  auto make_singlerow_matrix (const Flat& row, const Col& col) const
  {
    using RowMatrixType = Dune::Matrix<decltype((*this)(row,col[0u]))>;
    using MatrixType = Dune::Functions::SingleRowMatrix<RowMatrixType>;
    MatrixType mat(1,col.size());
    for (std::size_t j = 0; j < col.size(); ++j)
      mat[0][j] = (*this)(row,col[j]);
    return mat;
  }

public:

  template<std::size_t m, class W, std::size_t n>
  auto operator() (const UniformArray<Value,m>& row, const Array<W,n>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<std::size_t m, class W>
  auto operator() (const UniformArray<Value,m>& row, const Vector<W>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<std::size_t m, class W, std::size_t n>
  auto operator() (const UniformArray<Value,m>& row, const UniformArray<W,n>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<std::size_t m, class W>
  auto operator() (const UniformArray<Value,m>& row, const UniformVector<W>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<class W, std::size_t n>
  auto operator() (const UniformVector<Value>& row, const Array<W,n>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<class W>
  auto operator() (const UniformVector<Value>& row, const Vector<W>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<class W, std::size_t n>
  auto operator() (const UniformVector<Value>& row, const UniformArray<W,n>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

  template<class W>
  auto operator() (const UniformVector<Value>& row, const UniformVector<W>& col) const
  {
    return make_singlerow_matrix(row,col);
  }

public: // flat matrices

  using SparseFactory::operator();
};

template <class T>
struct ISTLSparseMatrixFactory
{
  auto operator() (const UniformVector<Value>& row, const UniformVector<Value>& col) const
  {
    return Dune::BCRSMatrix<T>{}; //(row.size(),col.size());
  }

  template <std::size_t m, std::size_t n>
  auto operator() (const UniformVector<UniformArray<Value,m>>& row, const UniformVector<UniformArray<Value,n>>& col) const
  {
    return Dune::BCRSMatrix<Dune::FieldMatrix<T,m,n>>{}; //(row.size(), col.size());
  }
};

struct ISTLSparsityPatternFactory
{
  auto operator() (const UniformVector<Value>& row, const UniformVector<Value>& col) const
  {
    return Dune::MatrixIndexSet(row.size(),col.size());
  }

  template <std::size_t m, std::size_t n>
  auto operator() (const UniformVector<UniformArray<Value,m>>& row, const UniformVector<UniformArray<Value,n>>& col) const
  {
    return Dune::MatrixIndexSet(row.size(), col.size());
  }
};

} // end namespace ContainerDescriptors


// Construct an istl matrix type compatible with the container descriptors
template<class SparseMatrixFactory = ContainerDescriptors::ISTLSparseMatrixFactory<double>,
         class RowDescriptor, class ColDescriptor>
auto istlMatrixFactory (const RowDescriptor& row, const ColDescriptor& col)
{
  auto factory = ContainerDescriptors::ISTLMatrixFactory<SparseMatrixFactory>{};
  return factory(row, col);
}

} // end namespace Dune::Functions

#endif // DUNE_FUNCTIONS_BACKENDS_ISTL_MATRIXFACTORY_HH
