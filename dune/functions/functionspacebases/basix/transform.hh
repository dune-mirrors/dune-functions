
void evaluateFunction (const Domain& x, std::vector<Range>& out) const
{
  out.resize(size());

  using GlobalRange = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange()>>;
#if USE_OUT_VECTOR_FOR_MDSPAN
  GlobalRange _global{&out[0][0], out.size()};
#else
  globalEvaluationBuffer_.resize(out.size() * dimRange());
  GlobalRange _global{globalEvaluationBuffer_.data(), out.size()};
#endif

  localEvaluationBuffer_.resize(out.size() * LocalBasis::dimRange());
  using LocalRange = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,LocalBasis::dimRange()>>;
  LocalRange _local{localEvaluationBuffer_.data(), out.size()};
  lb_->evaluateFunction(x, _local);

  auto J = geometry_->jacobian(x);
  auto K = geometry_->jacobianInverse(x);
  auto detJ = geometry_->integrationElement(x);

  using Jacobian = Std::mdspan<typename Geometry::ctype, Std::extents<std::size_t,Geometry::coorddimension,Geometry::mydimension>>;
  Jacobian _J{&J[0][0]};

  using JacobianInverse = Std::mdspan<typename Geometry::ctype, Std::extents<std::size_t,Geometry::mydimension,Geometry::coorddimension>>;
  JacobianInverse _K{&K[0][0]};

  auto map = lb_->basix().template map_fn<GlobalRange, LocalRange, Jacobian, JacobianInverse>();
  map(_global, _local, _J, detJ, _K);


#if !USE_OUT_VECTOR_FOR_MDSPAN
  // copy the output values back into the output variable
  for (std::size_t i = 0; i < out.size(); ++i)
    for (std::size_t j = 0; j < dimRange(); ++j)
      out[i][j] = _global(i,j);
#endif
}


/// \brief Evaluate all shape function jacobians in a point x
void evaluateJacobian (const Domain& x, std::vector<Jacobian>& out) const
{
  out.resize(size());

  using GlobalJacobian = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange(),dimDomain>>;
#if USE_OUT_VECTOR_FOR_MDSPAN
  GlobalPartials _global{&out[0][0], out.size()};
#else
  globalEvaluationBuffer_.resize(out.size() * dimRange());
  GlobalPartials _global{globalEvaluationBuffer_.data(), out.size()};
#endif

  localEvaluationBuffer_.resize(out.size() * LocalBasis::dimRange() * dimDomain);
  using LocalJacobian = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,LocalBasis::dimRange(),dimDomain>>;
  LocalJacobian _local{localEvaluationBuffer_.data(), out.size()};
  lb_->evaluateJacobian(x, _local);

#if !USE_OUT_VECTOR_FOR_MDSPAN
  // copy the output values back into the output variable
  for (std::size_t i = 0; i < out.size(); ++i)
    for (std::size_t j = 0; j < dimRange(); ++j)
      for (std::size_t k = 0; k < dimDomain; ++k)
        out[i][j][k] = _global(i,j,k);
#endif
}


/// \brief Evaluate all shape function partial derivatives with given orders in a point x
void partial (const std::array<unsigned int,dimDomain>& order,
              const Domain& x, std::vector<Range>& out) const
{
  out.resize(size());

  using GlobalPartials = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,dimRange()>>;
#if USE_OUT_VECTOR_FOR_MDSPAN
  GlobalPartials _global{&out[0][0], out.size()};
#else
  globalEvaluationBuffer_.resize(out.size() * dimRange());
  GlobalPartials _global{globalEvaluationBuffer_.data(), out.size()};
#endif

  localEvaluationBuffer_.resize(out.size() * LocalBasis::dimRange());
  using LocalPartials = Std::mdspan<F, Std::extents<std::size_t,Std::dynamic_extent,LocalBasis::dimRange()>>;
  LocalPartials _local{localEvaluationBuffer_.data(), out.size()};
  lb_->partial(order, x, _local);

  // TODO: Transformation needs to be implemented.
  DUNE_THROW(Dune::NotImplemented, "Transform not yet implemented.");

#if !USE_OUT_VECTOR_FOR_MDSPAN
  // copy the output values back into the output variable
  for (std::size_t i = 0; i < out.size(); ++i)
    for (std::size_t j = 0; j < dimRange(); ++j)
      out[i][j] = _global(i,j);
#endif
}