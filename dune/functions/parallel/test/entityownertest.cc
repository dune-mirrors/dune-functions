#include <config.h>

#include <iostream>
#include <set>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/parallel/entityowner.hh>
#include <dune/grid/yaspgrid.hh>

template <class Arg, class... Args>
std::string concat (Arg&& arg, Args&&... args)
{
  std::stringstream ss;
  ss << std::forward<Arg>(arg);
  if constexpr(sizeof...(args) > 0)
    ss << concat(std::forward<Args>(args)...);
  return ss.str();
}

template <class... Args>
void print (Args&&... args)
{
  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance();
  std::cout << concat("[", mpiHelper.rank(), "]  ", std::forward<Args>(args)...) << std::flush;
}


template <class GridView>
void test (std::string name, Dune::TestSuite& testSuite, GridView const& gridView)
{
  Dune::TestSuite subTestSuite{name};

  print(name, ":\n");
  const int rank = gridView.comm().rank();

  // auto const& idSet = gridView.grid().globalIdSet();
  // using IdType = typename GridView::Grid::GlobalIdSet::IdType;
  // std::set<IdType> localEntities;

  // // exchange all entities
  // for (auto const& e : elements(gridView))
  // {
  //   for (int codim = 0; codim <= GridView::dimension; ++codim)
  //     for (unsigned int i = 0; i < e.subEntities(codim); ++i)
  //       localEntities.insert(idSet.subId(e,i,codim));
  // }

  // std::vector<std::vector<IdType>> distributedEntities(gridView.comm().size());

  // {
  //   // 1. fill with local data
  //   distributedEntities[rank].insert(distributedEntities[rank].begin(),
  //                                   localEntities.begin(), localEntities.end());

  //   // 2. send local data to all other ranks
  //   for (int r = 0; r < gridView.comm().size(); ++r) {
  //     if (r != rank)
  //       gridView.comm().isend(distributedEntities[rank], r, 12345);
  //   }

  //   // 3. receiv data from all other ranks
  //   for (int r = 0; r < gridView.comm().size(); ++r) {
  //     if (r != rank)
  //       gridView.comm().recv(distributedEntities[r], r, 12345);
  //   }
  // }


  Dune::Functions::EntityOwner entityOwner{gridView};
  for (auto const& e : elements(gridView))
  {
    Dune::Hybrid::forEach(Dune::StaticIntegralRange<int,GridView::dimension+1>{}, [&](auto codim)
    {
      for (unsigned int i = 0; i < e.subEntities(codim); ++i)
      {
        Dune::PartitionType pt = e.template subEntity<codim>(i).partitionType();
        int ownerRank = entityOwner.owner(e,i,codim);
        if (pt == Dune::InteriorEntity)
          subTestSuite.check(ownerRank == rank);
        else if (pt == Dune::OverlapEntity || pt == Dune::FrontEntity || pt == Dune::GhostEntity)
          subTestSuite.check(ownerRank != rank);
        // else if (pt == Dune::BorderEntity) {
        //   auto const& entitySet = distributedEntities[ownerRank];
        //   auto it = std::find(entitySet.begin(), entitySet.end(), idSet.subId(e,i,codim));
        //   subTestSuite.check(it != entitySet.end());
        // }
      }
    });
  }

  testSuite.subTest(subTestSuite);
}


int main (int argc, char** argv)
{
  Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite testSuite{"EntityOwnerTest"};

  Dune::YaspGrid<2> grid2d0o({1.0,1.0}, {4,4}, 0u, 0); // overlap = 0
  test("YaspGrid<2> overlap=0", testSuite, grid2d0o.leafGridView());
  Dune::YaspGrid<2> grid2d1o({1.0,1.0}, {8,8}, 0u, 1); // overlap = 1
  test("YaspGrid<2> overlap=1", testSuite, grid2d1o.leafGridView());
  Dune::YaspGrid<2> grid2d2o({1.0,1.0}, {16,16}, 0u, 2); // overlap = 2
  test("YaspGrid<2> overlap=2", testSuite, grid2d2o.leafGridView());

  Dune::YaspGrid<3> grid3d0o({1.0,1.0,1.0}, {8,8,8}, 0u, 0); // overlap = 0
  test("YaspGrid<3> overlap=0", testSuite, grid3d0o.leafGridView());
  Dune::YaspGrid<3> grid3d1o({1.0,1.0,1.0}, {16,16,16}, 0u, 1); // overlap = 1
  test("YaspGrid<3> overlap=1", testSuite, grid3d1o.leafGridView());
  Dune::YaspGrid<3> grid3d2o({1.0,1.0,1.0}, {32,32,32}, 0u, 2); // overlap = 2
  test("YaspGrid<3> overlap=2", testSuite, grid3d2o.leafGridView());

  return testSuite.exit();
}