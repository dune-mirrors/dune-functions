#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DOUBLEPIOLA_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_DOUBLEPIOLA_HH
/** \brief Transforms shape function values and derivatives from reference
 * element coordinates to world coordinates using the double contravariant Piola
 * transform
 *
 *  This transformation preserves normal traces of vector fields.
 *  It is therefore the canonical transformation for tensorvalued
 * H(div)-conforming finite elements.
 *
 * See for example:
 *  Aznaran, Francis & Kirby, Robert & Farrell, Patrick. (2021). Transformations
 * for Piola-mapped elements.
 */
#include "dune/common/fmatrix.hh"
#include "dune/common/fvector.hh"
// #include "dune/functions/functionspacebases/arnoldwintherbasis.hh"
#include <cstddef>
namespace Dune {

template <class K, class S, int l, int m, int n,
          std::enable_if_t<IsNumber<S>::value, int> = 0>
auto operator*(FieldMatrix<FieldVector<K, l>, n, m>const& mat, S const &scalar)
{
  auto result = mat;
  for(auto i = 0u; i < m; ++i)
    for(auto j = 0u; j < n; ++i)
      result[i][j] *= scalar;
  return result;
}

namespace Functions::Impl {
  struct DoubleContravariantPiolaTransformator {
      /** \brief Double Piola-transform a set of shape-function values
       *
       * \param[in,out] values The values to be Piola-transformed
       */
      template <typename Values, typename LocalCoordinate, typename Geometry>
      static auto apply(Values &values, const LocalCoordinate &xi,
                        const Geometry &geometry) {
        auto jacobian = geometry.jacobian(xi);
        auto integrationElement = geometry.integrationElement(xi);
        assert(values[0].N() == values[0].M());
        assert(values[0].N() == jacobian.N());
        assert(jacobian.M() == jacobian.N());

        for (auto &value : values) {
          // auto tmp = value;
          // value = 0;
          // for (std::size_t k = 0; k < jacobian.N(); k++)
          //   for (std::size_t l = 0; l < jacobian.N(); l++)
          //     for (auto&& [jacobian_k_i, i] : sparseRange(jacobian[k]))
          //       for (auto&& [jacobian_l_j, j] : sparseRange(jacobian[l]))
          //         value[k][l] += jacobian_k_i * tmp[i][j] * jacobian_l_j;

          value = jacobian * value * transpose(jacobian);
          value = value/(integrationElement * integrationElement);
        }
        return;
      }

      /** \brief Piola-transform a set of shape-function derivatives
       *
       * \param[in,out] gradients The shape function derivatives to be
       * Piola-transformed
       *
       * \bug The current implementation works only for affine geometries.
       *   The Piola transformation for non-affine geometries requires
       *   second derivatives of the geometry, which we don't get
       *   from the dune-grid Geometry interface.
       */
      template <typename Gradients, typename LocalCoordinate, typename Geometry>
      static auto applyJacobian(Gradients &gradients, const LocalCoordinate &xi,
                                const Geometry &geometry) {
        apply(gradients, xi, geometry);
        return;

        // auto jacobian = geometry.jacobian(xi);
        // auto integrationElement = geometry.integrationElement(xi);
        // assert(gradients[0].N() == gradients[0].M());
        // assert(gradients[0].N() == jacobian.N());
        // assert(jacobian.M() == jacobian.N());

        // // We exploit that the result of the piola transform does only depend on
        // // the symmetric part of the input and is itself symmetric
        // for (auto &value : gradients) {
        //   auto tmp = value;
        //   value = 0;
        //   for (std::size_t k = 0; k < jacobian.N(); k++)
        //     for (std::size_t l = 0; l < jacobian.N(); l++)
        //       for (std::size_t m = 0; m < value[k][l].size(); ++m) {
        //         for (auto &&[jacobian_k_i, i] : sparseRange(jacobian[k]))
        //           for (auto &&[jacobian_l_j, j] : sparseRange(jacobian[l]))
        //             value[k][l][m] +=
        //                 jacobian_k_i * tmp[i][j][m] * jacobian_l_j;
        //         value[k][l][m] = value[k][l][m] *
        //                          (1. / integrationElement * integrationElement);
        //       }
        // }
      }

      /** \brief Wrapper around a callable that applies the inverse Piola
       * transform
       *
       * The LocalInterpolation implementations in dune-localfunctions expect
       * local-valued functions, but the ones dune-functions expect
       * global-valued ones.  Therefore, we need to stuff the inverse Piola
       * transform between dune-functions and dune-localfunctions, and this is
       * what this class does.
       */
      template <class Function, class LocalCoordinate, class Element>
      class LocalValuedFunction {
          const Function &f_;
          const Element &element_;

        public:
          LocalValuedFunction(const Function &f, const Element &element)
              : f_(f), element_(element) {}

          auto operator()(const LocalCoordinate &xi) const {
            auto globalValue = f_(xi);

            // Apply the inverse Piola transform
            auto jacobianInverse = element_.geometry().jacobianInverse(xi);
            auto integrationElement =
                element_.geometry().integrationElement(xi);

            globalValue = jacobianInverse * globalValue * transpose(jacobianInverse);

            globalValue *= integrationElement * integrationElement;

            return globalValue;
          }
      };
  };
} // namespace Functions::Impl
} // namespace Dune
#endif