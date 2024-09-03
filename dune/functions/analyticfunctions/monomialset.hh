// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_ANALYTICFUNCTIONS_MONOMIALSET_HH
#define DUNE_FUNCTIONS_ANALYTICFUNCTIONS_MONOMIALSET_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>



namespace Dune::Functions {

  namespace Impl {

    // Computes y = [1, x, x^2, x^3, ..., x^maxOrder] with x a single coordinate
    template<int maxOrder, class DomainFieldType, class RangeType>
    void computePowers(const DomainFieldType& x, RangeType& y)
    {
      if constexpr(maxOrder >= 0)
      {
        y[0] = 1;
        for(auto k : Dune::range(maxOrder))
          y[k+1] = y[k]*x;
      }
    }

  } // namespace Impl


  /**
   * \brief Set of all monomials up to degree \p maxOrder as vector valued functions.
   *
   * \ingroup FunctionImplementations
   *
   * \tparam RangeFieldType scalar type.
   * \tparam dim Domain dimension.
   * \tparam maxOrder Maximal monomial order.
   *
   * This vector valued function contains the dim-variate monomials up to
   * order maxOrder as components and models the
   * \ref Concept::DifferentiableFunction<Range(Domain)> concept.
   *
   * This class is a primary template with functions for documentation purposes only.
   * This is currently specialized for dim=1, dim=2, and dim=3 only.
   */
  template<class RangeFieldType, int dimension, int maxOrder>
  struct MonomialSet{
    static constexpr int dim = dimension;
    static constexpr int size = Dune::binomial(maxOrder + dim, dim);
    static_assert(dimension <= 3 and dimension >= 1, "This Primary Template should never be compiled");
    /**
     * \brief Return array of monomial exponents with shape `size x dim`
     *
     * The k-the entry of the returned array is the exponent
     * multiindex of the monomial corresponding the the k-th
     * component of the function.
     */
    static constexpr std::array<std::array<std::size_t,dim>,size> exponents();

    /**
    * \brief Return array of monomial evaluations
    *
    * The k-the entry of the returned array is the value
    * of the monomial corresponding the the k-th
    * entry of the return value of exponents().
    * Note that the return type depends on \p dim and might change in the future.
    * \tparam DomainFieldType The scalar type of the domain.
    */
    template<class DomainFieldType>
    constexpr Dune::FieldVector<RangeFieldType,size> operator()(const Dune::FieldVector<DomainFieldType,dim>& x) const;

    /**
    * \brief Set of all first order derivatives of monomials up to degree \p maxOrder as vector of vector valued functions.
    */
    struct Derivative
    {
      /**
      * \brief Return array of arrays of monomial derivatives
      *
      * The ith component of the k-the entry of the returned structure
      * is the derivative in direction i of the monomial corresponding
      * the the k-th entry of the return value of exponents().
      * Note that the return type depends on \p dim and might change in the future.
      * \tparam DomainFieldType The scalar type of the domain.
      */
      template<class DomainFieldType>
      constexpr Dune::FieldVector<RangeFieldType,size> operator()(const Dune::FieldVector<DomainFieldType,dim>& x) const;
      /**
      * \brief Set of all second order derivatives of monomials up to degree \p maxOrder as vector of matrix valued functions.
      */
      struct Hessian
      {
        /**
        * \brief Return array of Matrices of monomial derivatives
        *
        * The (i,j)th component of the k-the entry of the returned structure
        * is the derivative in direction (i,j) of the monomial corresponding
        * the the k-th entry of the return value of exponents().
        * Note that the return type depends on \p dim and might change in the future.
        * \tparam DomainFieldType The scalar type of the domain.
        */
        template<class DomainFieldType>
        constexpr Dune::FieldVector<RangeFieldType,size> operator()(const Dune::FieldVector<DomainFieldType,dim>& x) const;
      };
      /**
      * \brief Construct the Hessian object from a Derivative
      */
      constexpr friend auto derivative(const Derivative & d)
      {
        return Hessian{};
      }

    };
    /**
    * \brief Construct the Derivative object from a MonomialSet
    */
    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }


  };

  #ifndef DOXYGEN
  // Specialization for dim = 1
  template<class RangeFieldType, int maxOrder>
  struct MonomialSet<RangeFieldType, 1, maxOrder>
  {
    static constexpr int dim = 1;
    static constexpr int size = maxOrder+1;

    /**
     * \brief Return array of monomial exponent with shape `size x 1`
     *
     * The k-the entry of the returned array is the exponent
     * multiindex of the monomial corresponding the the k-th
     * component of the function.
     */
    static constexpr auto exponents()
    {
      auto p = std::array<std::array<std::size_t,1>,size>{};
      for(auto k : Dune::range(size))
        p[k][0] = k;
      return p;
    }

    template<class DomainFieldType>
    constexpr auto operator()(const Dune::FieldVector<DomainFieldType,1>& x) const
    {
      auto y = Dune::FieldVector<RangeFieldType,size>{};
      Impl::computePowers<maxOrder>(x[0], y);
      return y;
    }

    struct Derivative
    {
      template<class DomainFieldType>
      constexpr auto operator()(const Dune::FieldVector<DomainFieldType,1>& x) const
      {
        auto xPowers = Dune::FieldVector<RangeFieldType,size>{};
        Impl::computePowers<maxOrder-1>(x[0], xPowers);

        auto y = Dune::FieldMatrix<RangeFieldType,size,1>{};
        for(auto order : Dune::range(1, maxOrder+1))
          y[order][0] = order*xPowers[order-1];
        return y;
      }

      struct Hessian
      {
        template<class DomainFieldType>
        constexpr auto operator()(const Dune::FieldVector<DomainFieldType,1>& x) const
        {
          auto xPowers = std::array<RangeFieldType,maxOrder+1>{};
          Impl::computePowers<maxOrder-2>(x[0],xPowers);

          auto y = Dune::FieldVector<Dune::FieldMatrix<RangeFieldType,1,1>,size>{};
          for(auto order : Dune::range(maxOrder+1))
            if (order-1 > 0)
              y[order][0][0] = order*(order-1)*xPowers[order-2];
          return y;
        }
      };

      constexpr friend auto derivative(const Derivative& d)
      {
        return Hessian{};
      }

    };

    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };



  // Specialization for dim = 2
  template<class RangeFieldType, int maxOrder>
  struct MonomialSet<RangeFieldType, 2, maxOrder>
  {
    static constexpr int dim = 2;
    static constexpr int size = (maxOrder+1)*(maxOrder+2)/2;

    /**
     * \brief Return array of monomial exponents with shape `size x 2`
     *
     * The k-the entry of the returned array is the exponent
     * multiindex of the monomial corresponding the the k-th
     * component of the function.
     */
    static constexpr auto exponents()
    {
      auto p = std::array<std::array<std::size_t,2>,size>{};
      std::size_t index = 0;
      for(auto order : Dune::range(maxOrder+1))
      {
        for(auto k : Dune::range(order+1))
        {
          p[index][0] = order-k;
          p[index][1] = k;
          ++index;
        }
      }
      return p;
    }

    template<class DomainFieldType>
    constexpr auto operator()(const Dune::FieldVector<DomainFieldType,2>& x) const
    {
      auto xPowers = std::array<Dune::FieldVector<RangeFieldType,size>,dim>{};
      for(auto j : Dune::range(dim))
        Impl::computePowers<maxOrder>(x[j], xPowers[j]);

      auto y = Dune::FieldVector<RangeFieldType,size>{};
      std::size_t index = 0;
      for(auto order : Dune::range(maxOrder+1))
      {
        for(auto k : Dune::range(order+1))
        {
          y[index] = xPowers[0][order-k]*xPowers[1][k];
          ++index;
        }
      }
      return y;
    }

    struct Derivative
    {
      template<class DomainFieldType>
      constexpr auto operator()(const Dune::FieldVector<DomainFieldType,2>& x) const
      {
        auto xPowers = std::array<Dune::FieldVector<RangeFieldType,size>,dim>{};
        for(auto j : Dune::range(dim))
          Impl::computePowers<maxOrder-1>(x[j], xPowers[j]);

        auto y = Dune::FieldMatrix<RangeFieldType,size,2>{};
        std::size_t index = 0;
        for(auto order : Dune::range(maxOrder+1))
        {
          for(auto k : Dune::range(order+1))
          {
            if (order-k > 0)
              y[index][0] = (order-k)*xPowers[0][order-k-1]*xPowers[1][k];
            if (k > 0)
              y[index][1] = k*xPowers[0][order-k]*xPowers[1][k-1];
            ++index;
          }
        }
        return y;
      }

      struct Hessian
      {
        template<class DomainFieldType>
        constexpr auto operator()(const Dune::FieldVector<DomainFieldType,2>& x) const
        {
          auto xPowers = std::array<std::array<RangeFieldType,maxOrder+1>,dim>{};
          for(auto j : Dune::range(dim))
            Impl::computePowers<maxOrder-2>(x[j], xPowers[j]);

          auto y = Dune::FieldVector<Dune::FieldMatrix<RangeFieldType,2,2>,size>{};
          std::size_t index = 0;
          for(auto order : Dune::range(maxOrder+1))
          {
            for(auto k : Dune::range(order+1))
            {
              if (order-k > 1)
                y[index][0][0] = (order-k-1)*(order-k)*xPowers[0][order-k-2]*xPowers[1][k];
              if (k > 0 and order-k > 0){
                auto mixed = k*(order-k)*xPowers[0][order-k-1]*xPowers[1][k-1];
                y[index][0][1]= mixed;
                y[index][1][0]= mixed;
              }
              if (k > 1)
                y[index][1][1] = k*(k-1)*xPowers[0][order-k]*xPowers[1][k-2];

              ++index;
            }
          }
          return y;
        }
      };

      constexpr friend auto derivative(const Derivative & d)
      {
        return Hessian{};
      }

    };

    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };


    // Specialization for dim = 3
  template<class RangeFieldType, int maxOrder>
  struct MonomialSet<RangeFieldType, 3, maxOrder>
  {
    static constexpr int dim = 3;
    static constexpr int size = Dune::binomial(maxOrder + dim, dim);

    /**
     * \brief Return array of monomial exponents of shape `size x 3`
     *
     * The k-the entry of the returned array is the exponent
     * multiindex of the monomial corresponding the the k-th
     * component of the function. Note, that this ordering is tensorbased, i.e. x,y,z,xx,xy,yy,xz,yz,zz, ...
     */
    static constexpr auto exponents()
    {
      auto p = std::array<std::array<std::size_t,3>,size>{};
      std::size_t index = 0;
      for(auto order : Dune::range(maxOrder+1))
      {
        for(auto k : Dune::range(order+1))
        {
          for (auto l : Dune::range(order-k+1))
          {
            p[index][0] = order-k-l;
            p[index][1] = l;
            p[index][2] = k;
            ++index;
          }

        }
      }
      return p;
    }

    template<class DomainFieldType>
    constexpr auto operator()(const Dune::FieldVector<DomainFieldType,3>& x) const
    {
      auto xPowers = std::array<std::array<RangeFieldType,maxOrder+1>,dim>{};
      for(auto j : Dune::range(dim))
        Impl::computePowers<maxOrder>(x[j], xPowers[j]);

      auto y = Dune::FieldVector<RangeFieldType,size>{};
      std::size_t index = 0;
      for(auto order : Dune::range(maxOrder+1))
      {
        for(auto k : Dune::range(order+1))
        {
          for (auto l : Dune::range(order-k+1))
          {
            y[index] = xPowers[0][order-k-l]*xPowers[1][l]*xPowers[2][k];
            ++index;
          }
        }
      }
      return y;
    }

    struct Derivative
    {
      template<class DomainFieldType>
      constexpr auto operator()(const Dune::FieldVector<DomainFieldType,3>& x) const
      {
        auto xPowers = std::array<std::array<RangeFieldType,maxOrder+1>,dim>{};
        for(auto j : Dune::range(dim))
        {
          xPowers[j][0] = 1.0;
          for(auto k : Dune::range(maxOrder))
            xPowers[j][k+1] = xPowers[j][k]*x[j];
        }

        auto y = Dune::FieldMatrix<RangeFieldType,size,3>{};
        std::size_t index = 0;
        for(auto order : Dune::range(maxOrder+1))
        {
          for(auto k : Dune::range(order+1))
          {
            for (auto l : Dune::range(order-k+1))
            {
              if (order-k-l > 0)
                y[index][0] = (order-k-l)*xPowers[0][order-k-l-1]*xPowers[1][l]*xPowers[2][k];
              if (l > 0)
                y[index][1] = l*xPowers[0][order-k-l]*xPowers[1][l-1]*xPowers[2][k];
              if (k > 0)
                y[index][2] = k*xPowers[0][order-k-l]*xPowers[1][l]*xPowers[2][k-1];
              ++index;
            }
          }
        }
        return y;
      }

      struct Hessian
      {
        template<class DomainFieldType>
        constexpr auto operator()(const Dune::FieldVector<DomainFieldType,3>& x) const
        {
          auto xPowers = std::array<std::array<RangeFieldType,maxOrder+1>,dim>{};
          for(auto j : Dune::range(dim))
            Impl::computePowers<maxOrder>(x[j], xPowers[j]);

          auto y = Dune::FieldVector<Dune::FieldMatrix<RangeFieldType,3,3>,size>{};
          std::size_t index = 0;
          for(auto order : Dune::range(maxOrder+1))
          {
            for(auto k : Dune::range(order+1))
            {
              for (auto l : Dune::range(order-k+1))
              {
                // xx
                if (order-k-l-1 > 0)
                  y[index][0][0] = (order-k-l)*(order-k-l-1)*xPowers[0][order-k-l-2]*xPowers[1][l]*xPowers[2][k];
                // xy and yx
                if (order-k-l > 0 and l > 0){
                  y[index][0][1] = (order-k-l)*l*xPowers[0][order-k-l-1]*xPowers[1][l-1]*xPowers[2][k];
                  y[index][1][0] = y[index][0][1];
                }
                // yy
                if (l-1 > 0)
                  y[index][1][1] = l*(l-1)*xPowers[0][order-k-l]*xPowers[1][l-2]*xPowers[2][k];
                // xz and zx
                if (k > 0 and order-k-l > 0)
                {
                  y[index][0][2] = (order-k-l)*k*xPowers[0][order-k-l-1]*xPowers[1][l]*xPowers[2][k-1];
                  y[index][2][0] = y[index][0][2];
                }
                // yz
                if (l > 0 and k > 0)
                {
                  y[index][1][2] = l*k*xPowers[0][order-k-l]*xPowers[1][l-1]*xPowers[2][k-1];
                  y[index][2][1] = y[index][1][2];
                }
                // zz
                if (k-1 > 0)
                  y[index][2][2] = (k-1)*k*xPowers[0][order-k-l]*xPowers[1][l]*xPowers[2][k-2];
                ++index;
              }
            }
          }
          return y;
        }

      };

      constexpr friend auto derivative(const Derivative & d)
      {
        return Hessian{};
      }

    };

    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };
  #endif
} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_MONOMIALSET_HH
