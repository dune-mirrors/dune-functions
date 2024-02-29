#ifndef DUNE_FUNCTIONS_PYTHON_RANGETYPE_HH
#define DUNE_FUNCTIONS_PYTHON_RANGETYPE_HH

#include <tuple>

#include <dune/python/common/dynvector.hh>
#include <dune/python/common/fvector.hh>
#include <dune/python/common/tuplevector.hh>

#include <dune/typetree/nodetags.hh>

namespace Dune
{

  namespace Python
  {

    /**
     * \brief Register container types for the range type of discrete functions
     * to the pythong bindings by calling the corresponding register functions
     * recursively.
    */
    template <class R>
    struct RegisterRange
    {
      static_assert(std::is_arithmetic_v<R>);

      // nothing to register, as K is a basic type
      static void apply(pybind11::module scope) {}
    };

    template <class K, std::size_t n>
    struct RegisterRange<std::array<K,n>>
    {
      static void apply(pybind11::module scope)
      {
        RegisterRange<K>::apply(scope);
      }
    };

    template <class K, int n>
    struct RegisterRange<Dune::FieldVector<K,n>>
    {
      static void apply(pybind11::module scope)
      {
        RegisterRange<K>::apply(scope);
        registerFieldVector<K,n>(scope);
      }
    };

    template <class K>
    struct RegisterRange<std::vector<K>>
    {
      static void apply(pybind11::module scope)
      {
        RegisterRange<K>::apply(scope);
      }
    };

    template <class K>
    struct RegisterRange<Dune::DynamicVector<K>>
    {
      static void apply(pybind11::module scope)
      {
        RegisterRange<K>::apply(scope);
        registerDynamicVector<K>(scope);
      }
    };

    template <class... K>
    struct RegisterRange<Dune::TupleVector<K...>>
    {
      static void apply(pybind11::module scope)
      {
        (RegisterRange<K>::apply(scope),...);
        registerTupleVector<K...>(scope);
      }
    };


    /**
     * \brief Data type to use as range type of discrete functions.
     *
     * The data type is constructed by mapping the hierarchic basis tree
     * into a container type. Power nodes are thereby mapped into FieldVector
     * or DynamicVector and composite nodes into TupleVector types. From
     * the leaf node of the basis the range type is extracted from the local
     * basis range type.
     *
     * \tparam Node    The type of the basis node.
     * \tparam NodeTag Tag representing the type of the node.
     **/
    template <class Node, class NodeTag = typename Node::NodeTag>
    struct RangeType;

    template <class Node>
    struct RangeType<Node, Dune::TypeTree::LeafNodeTag>
    {
      using FiniteElement = typename Node::FiniteElement;
      using LocalBasis = typename FiniteElement::Traits::LocalBasisType;
      using type = typename LocalBasis::Traits::RangeType;
    };

    template <class Node>
    struct RangeType<Node, Dune::TypeTree::PowerNodeTag>
    {
      using child_type = typename RangeType<typename Node::ChildType>::type;
      using type = Dune::FieldVector<child_type, Node::degree()>;
    };

    template <class Node>
    struct RangeType<Node, Dune::TypeTree::DynamicPowerNodeTag>
    {
      using child_type = typename RangeType<typename Node::ChildType>::type;
      using type = Dune::DynamicVector<child_type>;
    };

    template <class Node>
    struct RangeType<Node, Dune::TypeTree::CompositeNodeTag>
    {
      template <class> struct make_type;
      template <class... Children>
      struct make_type<std::tuple<Children...>>
      {
        using type = Dune::TupleVector<typename RangeType<Children>::type...>;
      };
      using type = typename make_type<typename Node::ChildType>::type;
    };

  } // end namespace Python
} // end namespace Dune

#endif // DUNE_FUNCTIONS_PYTHON_RANGETYPE_HH
