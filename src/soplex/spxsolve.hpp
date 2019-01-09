
  template <class R>
  typename SPxSolverBase<R>::Status SPxSolverBase<R>::getPrimalSol(VectorBase<R>& p_vector) const
  {

    if (!isInitialized())
      {
        /* exit if presolving/simplifier cleared the problem */
        if (status() == NO_PROBLEM)
          return status();
        throw SPxStatusException("XSOLVE06 Not Initialized");
      }
    if (rep() == ROW)
      p_vector = coPvec();
    else
      {
        const typename SPxBasisBase<R>::Desc& ds = this->desc();

        for (int i = 0; i < this->nCols(); ++i)
          {
            switch (ds.colStatus(i))
              {
              case SPxBasisBase<R>::Desc::P_ON_LOWER :
                p_vector[i] = SPxLPBase<R>::lower(i);
                break;
              case SPxBasisBase<R>::Desc::P_ON_UPPER :
              case SPxBasisBase<R>::Desc::P_FIXED :
                p_vector[i] = SPxLPBase<R>::upper(i);
                break;
              case SPxBasisBase<R>::Desc::P_FREE :
                p_vector[i] = 0;
                break;
              case SPxBasisBase<R>::Desc::D_FREE :
              case SPxBasisBase<R>::Desc::D_ON_UPPER :
              case SPxBasisBase<R>::Desc::D_ON_LOWER :
              case SPxBasisBase<R>::Desc::D_ON_BOTH :
              case SPxBasisBase<R>::Desc::D_UNDEFINED :
                break;
              default:
                throw SPxInternalCodeException("XSOLVE07 This should never happen.");
              }
          }
        for (int j = 0; j < dim(); ++j)
          {
            if (this->baseId(j).isSPxColId())
              p_vector[ this->number(SPxColId(this->baseId(j))) ] = fVec()[j];
          }
      }
    return status();
  }

