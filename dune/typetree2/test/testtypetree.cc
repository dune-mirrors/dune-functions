#include "config.h"

#include <dune/common/classname.hh>

#include "typetreetestswitch.hh"

#if TEST_TYPETREE_INVALID

int main()
{
  return 0;
}

#else

#include "typetreetestutility.hh"


int main(int argc, char** argv)
{

  // basic tests

  // leaf node
  SimpleLeaf sl1;

  typedef SimplePower<SimpleLeaf,3> SP1;
  SP1 sp1_1{sl1, sl1, sl1};

  SimpleLeaf sl2;
  SP1 sp1_2(Dune::index_constant<3>{}, sl2);
  SP1 sp1_2a(Dune::index_constant<3>{}, sl2);

  typedef SimpleComposite<SimpleLeaf,SP1,SimpleLeaf> SC1;
  SC1 sc1_1(sl1,sp1_2,sl2);

  TreePrinter treePrinter;
  Dune::TypeTree::applyToTree(const_cast<const SC1&>(sc1_1),treePrinter);

  typedef SimpleComposite<SimpleLeaf,SimpleLeaf,SimpleLeaf> SC2;
  SC2 sc2(sl1,sl1,sl1);

  typedef SimpleComposite<SimpleLeaf,SP1,SimpleLeaf,SC1> SVC1;
  SVC1 svc1_1(sl1,sp1_1,sl2,sc1_1);

  typedef SimpleDynamicPower<SVC1> SDP1;
  SDP1 sdp_1(2, svc1_1);
  Dune::TypeTree::applyToTree(sdp_1,treePrinter);

  SP1 sp1_3(SimpleLeaf(),SimpleLeaf(),sl1);
  Dune::TypeTree::applyToTree(sp1_3,TreePrinter());

  SVC1 svc1_2(SimpleLeaf(),SP1(sp1_2),sl2,const_cast<const SC1&>(sc1_1));

  typedef SimpleComposite<SimpleLeaf,SC2,SimpleLeaf,SC1> SVC2;
  SVC2 svc2_1(sl1,sc2,sl2,sc1_1);

  Dune::TypeTree::applyToTreePair(svc1_2,svc2_1,PairPrinter());

  typedef SimpleDynamicPower<SimpleLeaf> SDP;
  SDP sdp(2, sl1);

  // Test valid and invalid child access. Invalid access should be caught at compile time
  auto const _0 = Dune::TypeTree::index_constant<0>();
  auto const _1 = Dune::TypeTree::index_constant<1>();
  auto const _2 = Dune::TypeTree::index_constant<2>();

  // 1: valid access
  auto x1 = Dune::TypeTree::child(sp1_1, _0);
#ifdef FAILURE2
  // 2: invalid access (too few children)
  {
    auto const _3 = Dune::TypeTree::index_constant<3>();
    auto x2 = Dune::TypeTree::child(sp1_1, _3);
  }
#endif
#ifdef FAILURE3
  // 3: invalid access (child has no children)
  auto x3 = Dune::TypeTree::child(sp1_1, _0, _0);
#endif

  // 4: valid access
  auto x4 = Dune::TypeTree::child(sc1_1, _1, _2);
#ifdef FAILURE5
  // 5: invalid access (too few children)
  {
    auto const _3 = Dune::TypeTree::index_constant<3>();
    auto x5 = Dune::TypeTree::child(sc1_1, _3);
  }
#endif
#ifdef FAILURE6
  // 6: invalid access (child has no children)
  auto x6 = Dune::TypeTree::child(sc1_1, _0, _0);
#endif

  return 0;
}

#endif
