
#include <dune/typetree2/typetree.hh>
#include <dune/typetree/pairtraversal.hh>

#include <iostream>

struct Counter
{
  Counter()
    : _id(_ids++)
  {
    std::cout << "Constructed id = " << id() << std::endl;
  }

  Counter(const Counter& rhs)
    : _id(_ids++)
  {
    rhs.assert_valid();
    std::cout << "Copy-Constructed id = " << id() << " from id = " << rhs.id() << std::endl;
  }

  Counter(Counter&& rhs)
    : _id(rhs._id)
  {
    rhs.assert_valid();
    rhs._id = -1;
    std::cout << "Move-Constructed id = " << id() << std::endl;
  }

  ~Counter()
  {
    std::cout << "Destructed id = " << _id << std::endl;
  }

  Counter& operator=(const Counter& rhs)
  {
    rhs.assert_valid();
    assert_valid();
    std::cout << "Assigned id = " << id() << " from id = " << rhs.id() << std::endl;
    return *this;
  }

  Counter& operator=(Counter&& rhs)
  {
    assert_valid();
    rhs.assert_valid();
    std::cout << "Move-Assigned id = " << id() << " from id = " << rhs.id() << std::endl;
    rhs._id = -1;
    return *this;
  }

  int id() const
  {
    assert_valid();
    return _id;
  }

  void assert_valid() const
  {
    assert(_id != -1);
  }

  int _id;
  static int _ids;
};

int Counter::_ids = 0;

struct SimpleLeaf
  : public Dune::TypeTree2::LeafNode
  , public Counter
{
  static const char* name()
  {
    return "SimpleLeaf";
  }

  using Super = Dune::TypeTree2::LeafNode;
  using Super::Super;
};

struct SimpleLeafDerived
  : public SimpleLeaf
{
  static const char* name()
  {
    return "SimpleLeafDerived";
  }

  using Super = SimpleLeaf;
  using Super::Super;
};

template<class T, std::size_t k>
struct SimplePower
  : public Dune::TypeTree2::StaticNonUniformTypeTree<T,k>
  , public Counter
{
  static const char* name()
  {
    return "SimplePower";
  }

  using Super = Dune::TypeTree2::StaticNonUniformTypeTree<T,k>;
  using Super::Super;
};

template<class... Children>
struct SimpleComposite
  : public Dune::TypeTree2::VariadicNonUniformTypeTree<Children...>
  , public Counter
{
  static const char* name()
  {
    return "SimpleComposite";
  }

  using Super = Dune::TypeTree2::VariadicNonUniformTypeTree<Children...>;
  using Super::Super;
};

template<class T>
struct SimpleDynamicPower
  : public Dune::TypeTree2::DynamicNonUniformTypeTree<T>
  , public Counter
{
  static const char* name()
  {
    return "SimpleDynamicPower";
  }

  using Super = Dune::TypeTree2::DynamicNonUniformTypeTree<T>;
  using Super::Super;
};


struct TreePrinter
  : public Dune::TypeTree::TreeVisitor
  , public Dune::TypeTree::DynamicTraversal
{

  template<class T, class TreePath>
  void leaf(const T& t, TreePath treePath) const
  {
    pre(t,treePath);
  }

  template<class T, class TreePath>
  void pre(const T& t, TreePath treePath) const
  {
    for (std::size_t i = 0; i < treePath.size(); ++i)
      std::cout << "  ";
    std::cout << t.name() << " " << t.id() << std::endl;
  }
};




struct PairPrinter
  : public Dune::TypeTree::TreePairVisitor
  , public Dune::TypeTree::DynamicTraversal
{

  template<class T1, class T2, class TreePath>
  void leaf(const T1& t1, const T2& t2, TreePath treePath) const
  {
    pre(t1,t2,treePath);
  }

  template<class T1, class T2, class TreePath>
  void pre(const T1& t1, const T2& t2, TreePath treePath) const
  {
    for (std::size_t i = 0; i < treePath.size(); ++i)
      std::cout << "  ";
    std::cout << t1.name() << " " << t1.id() << "      " << t2.name() << " " << t2.id() << std::endl;
  }
};
