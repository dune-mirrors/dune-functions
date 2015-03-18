// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/functions/common/functionfromcallable.hh>

namespace Dune {
namespace Functions {

namespace Imp {

struct AllTrueBitSetVector
{
  struct AllTrueBitSet
  {
      bool test(int i) const { return true; }
  } allTrue_;

  const AllTrueBitSet& operator[](int i) const
  {
      return allTrue_;
  }
};



template<class C>
struct FlatIndexContainerAccess
{
  template<class T>
  static void setEntry(C& x, int i, const T& xi)
  {
    x = xi;
  }

  static auto getEntry(const C& x, int i) -> decltype(x)
  {
    return x;
  }

  static int size(const C& x)
  {
    return 1;
  }
};

template<class K, int n>
struct FlatIndexContainerAccess<typename Dune::FieldVector<K, n> >
{
  typedef typename Dune::FieldVector<K, n> RT;

  template<class RFT>
  static void setEntry(RT& x, int i, const RFT& xi)
  {
    x[i] = xi;
  }

  static auto getEntry(const RT& x, int i) -> decltype(x[i])
  {
    return x[i];
  }

  static int size(const RT& x)
  {
    return n;
  }
};

template<class K, int n, int m>
struct FlatIndexContainerAccess<typename Dune::FieldMatrix<K, n, m> >
{
  typedef typename Dune::FieldMatrix<K, n, m> RT;

  template<class RFT>
  static void setEntry(RT& x, int i, const RFT& xi)
  {
    x[i/m][i%m] = xi;
  }

  static auto getEntry(const RT& x, int i) -> decltype(x[i/m][i%m])
  {
    return x[i/m][i%m];
  }

  static int size(const RT& x)
  {
    return n*m;
  }
};

template<class K>
struct FlatIndexContainerAccess<typename std::vector<K> >
{
  typedef typename std::vector<K> RT;

  template<class RFT>
  static void setEntry(RT& x, int i, const RFT& xi)
  {
    x[i] = xi;
  }

  static auto getEntry(const RT& x, int i) -> decltype(x[i])
  {
    return x[i];
  }

  static int size(const RT& x)
  {
      return x.size();
  }
};

template<int k, class Alloc>
struct FlatIndexContainerAccess<Dune::BitSetVectorConstReference<k, Alloc> >
{
  typedef typename Dune::BitSetVectorConstReference<k, Alloc> RT;

  template<class RFT>
  static void setEntry(RT& x, int i, const RFT& xi)
  {
    DUNE_THROW(Dune::Exception, "Can't modify Dune::BitSetVector<k>::const_reference");
  }

  static bool getEntry(const RT& x, int i)
  {
    return x.test(i);
  }

  static int size(const RT& x)
  {
      return x.size();
  }
};

template<>
struct FlatIndexContainerAccess<AllTrueBitSetVector::AllTrueBitSet>
{
  template<class RFT>
  static void setEntry(AllTrueBitSetVector::AllTrueBitSet& x, int i, const RFT& xi)
  {
    DUNE_THROW(Dune::Exception, "Can't modify AllTrueBitSet");
  }

  static bool getEntry(const AllTrueBitSetVector::AllTrueBitSet& x, int i)
  {
    return true;
  }

  static int size(const AllTrueBitSetVector::AllTrueBitSet& x)
  {
    DUNE_THROW(Dune::NotImplemented, "AllTrueBitSet has no size");
  }
};




} // namespace Imp




/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatIndexContainerAccess.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 * \param bitVector A vector with flags marking ald DOFs that should be interpolated
 */
template <class B, class C, class F, class BV>
void interpolate(const B& basis, C& coeff, F&& f, BV&& bitVector)
{
  using GridView = typename B::GridView;
  using Element = typename GridView::template Codim<0>::Entity;

  using FiniteElement = typename B::LocalView::Tree::FiniteElement;
  using FunctionBaseClass = typename Dune::LocalFiniteElementFunctionBase<FiniteElement>::type;

  using LocalBasisRange = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
  using LocalDomain = typename Element::Geometry::LocalCoordinate;

  using CoefficientBlock = typename std::decay<decltype(coeff[0])>::type;
  using BitVectorBlock = typename std::decay<decltype(bitVector[0])>::type;

  auto&& gridView = basis.gridView();

  // Make a grid function supporting local evaluation out of f
  auto gf = makeGridViewFunction(std::forward<F>(f), gridView);

  // Obtain a local view of f
  auto localF = localFunction(gf);

  // Note that we capture j by reference. Hence we can switch
  // the selected component later on by modifying j. Maybe we
  // should avoid this naughty statefull lambda hack in favor
  // of a separate helper class.
  int j=0;
  auto localFj = [&](const LocalDomain& x){
    using FunctionRange = typename std::decay<decltype(localF(LocalDomain(0)))>::type;
    return Imp::FlatIndexContainerAccess<FunctionRange>::getEntry(localF(x), j);
  };

  using FunctionFromCallable = typename Dune::Functions::FunctionFromCallable<LocalBasisRange(LocalDomain), decltype(localFj), FunctionBaseClass>;

  auto basisIndexSet = basis.indexSet();
  coeff.resize(basisIndexSet.size());

  auto processed = Dune::BitSetVector<1>(basisIndexSet.size(), false);
  auto interpolationValues = std::vector<LocalBasisRange>();

  auto localView = basis.localView();
  auto localIndexSet = basisIndexSet.localIndexSet();

  auto blockSize = Imp::FlatIndexContainerAccess<CoefficientBlock>::size(coeff[0]);

  for (const auto& e : elements(gridView))
  {
    localView.bind(e);
    localIndexSet.bind(localView);
    localF.bind(e);

    const auto& fe = localView.tree().finiteElement();

    // check if all components have already been processed
    bool allProcessed = true;
    for (size_t i=0; i<fe.localBasis().size(); ++i)
    {
      // if index was already processed we don't need to do further checks
      auto index = localIndexSet.index(i)[0];
      if (processed[index][0])
        continue;

      // if index was not processed, check if any entry is marked for interpolation
      auto&& bitVectorBlock = bitVector[index];
      for(int k=0; k<blockSize; ++k)
      {
        if (Imp::FlatIndexContainerAccess<BitVectorBlock>::getEntry(bitVectorBlock,k))
        {
          allProcessed = false;
          break;
        }
      }
    }

    if (not(allProcessed))
    {
      // We loop over j defined above and thus over the components of the
      // range type of localF.
      for(j=0; j<blockSize; ++j)
      {

        // We cannot use localFj directly because interpolate requires a Dune::VirtualFunction like interface
        fe.localInterpolation().interpolate(FunctionFromCallable(localFj), interpolationValues);
        for (size_t i=0; i<fe.localBasis().size(); ++i)
        {
          size_t index = localIndexSet.index(i)[0];
          auto interpolateHere = Imp::FlatIndexContainerAccess<BitVectorBlock>::getEntry(bitVector[index],j);

          if (not(processed[index][0]) and interpolateHere)
            Imp::FlatIndexContainerAccess<CoefficientBlock>::setEntry(coeff[index], j, interpolationValues[i]);
        }
      }
      for (size_t i=0; i<fe.localBasis().size(); ++i)
        processed[localIndexSet.index(i)[0]][0] = true;
    }
  }
}



/**
 * \brief Interpolate given function in discrete function space
 *
 * Notice that this will only work if the range type of f and
 * the block type of coeff are compatible and supported by
 * FlatIndexContainerAccess.
 *
 * \param basis Global function space basis of discrete function space
 * \param coeff Coefficient vector to represent the interpolation
 * \param f Function to interpolate
 */
template <class B, class C, class F>
void interpolate(const B& basis, C& coeff, F&& f)
{
  interpolate (basis, coeff, f, Imp::AllTrueBitSetVector());
}


} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_INTERPOLATE_HH
