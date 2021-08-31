// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <functional>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/functions/common/referencehelper.hh>

#include <dune/common/test/testsuite.hh>



class CopyCounter
{
public:
  CopyCounter() : count_(0) {}
  CopyCounter(std::size_t count) : count_(count) {}
  CopyCounter(const CopyCounter& other) :
    count_(other.count_ + 1)
  {}

  auto& getCount() {
    return count_;
  }

  const auto& getCount() const {
    return count_;
  }

  void setCount(std::size_t count) {
    count_ = count;
  }

private:
  mutable std::size_t count_;
};


int main (int argc, char *argv[]) try
{
  // Set up MPI, if available
  Dune::MPIHelper::instance(argc, argv);

  Dune::TestSuite suite;

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with mutable l-value");
    CopyCounter c;
    Dune::Functions::resolveRef(c).setCount(42);
    Dune::Functions::resolveRef(c).getCount();
    subSuite.check(Dune::Functions::resolveRef(c).getCount() == 42, "Checking resolveRef");
    subSuite.check(not Dune::Functions::IsReferenceWrapper_v<decltype(c)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with const l-value");
    const CopyCounter c(42);
    Dune::Functions::resolveRef(c).getCount();
    subSuite.check(Dune::Functions::resolveRef(c).getCount() == 42, "Checking resolveRef");
    subSuite.check(not Dune::Functions::IsReferenceWrapper_v<decltype(c)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with mutable reference_wrapper of mutable l-value");
    CopyCounter c;
    auto c_ref = std::ref(c);
    Dune::Functions::resolveRef(c_ref).setCount(42);
    Dune::Functions::resolveRef(c_ref).getCount();
    subSuite.check(Dune::Functions::resolveRef(c_ref).getCount() == 42, "Checking resolveRef");
    subSuite.check(Dune::Functions::IsReferenceWrapper_v<decltype(c_ref)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with const reference_wrapper of mutable l-value");
    CopyCounter c;
    const auto c_ref = std::ref(c);
    Dune::Functions::resolveRef(c_ref).setCount(42);
    Dune::Functions::resolveRef(c_ref).getCount();
    subSuite.check(Dune::Functions::resolveRef(c_ref).getCount() == 42, "Checking resolveRef");
    subSuite.check(Dune::Functions::IsReferenceWrapper_v<decltype(c_ref)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with mutable reference_wrapper of const l-value");
    const CopyCounter c(42);
    auto c_ref = std::ref(c);
    Dune::Functions::resolveRef(c_ref).getCount();
    subSuite.check(Dune::Functions::resolveRef(c_ref).getCount() == 42, "Checking resolveRef");
    subSuite.check(Dune::Functions::IsReferenceWrapper_v<decltype(c_ref)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with const reference_wrapper of const l-value");
    const CopyCounter c(42);
    const auto c_ref = std::ref(c);
    Dune::Functions::resolveRef(c_ref).getCount();
    subSuite.check(Dune::Functions::resolveRef(c_ref).getCount() == 42, "Checking resolveRef");
    subSuite.check(Dune::Functions::IsReferenceWrapper_v<decltype(c_ref)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with const reference_wrapper of const l-value (via std::cref)");
    CopyCounter c(42);
    auto c_ref = std::cref(c);
    Dune::Functions::resolveRef(c_ref).getCount();
    subSuite.check(Dune::Functions::resolveRef(c_ref).getCount() == 42, "Checking resolveRef");
    subSuite.check(Dune::Functions::IsReferenceWrapper_v<decltype(c_ref)>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  suite.subTest([]() {
    Dune::TestSuite subSuite("Checking with const reference_wrapper r-value of mutable l-value");
    CopyCounter c;
    Dune::Functions::resolveRef(std::ref(c)).setCount(42);
    Dune::Functions::resolveRef(std::ref(c)).getCount();
    subSuite.check(Dune::Functions::resolveRef(std::ref(c)).getCount() == 42, "Checking resolveRef");
    subSuite.check(Dune::Functions::IsReferenceWrapper_v<decltype(std::ref(c))>, "Checking IsReferenceWrapper_v");
    return subSuite;
  }());

  return suite.exit();
}
catch (std::exception& e) {
  std::cout << e.what() << std::endl;
  return 1;
}
