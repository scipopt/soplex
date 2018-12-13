
  /// checks a solution for feasibility
  template <class R>
  bool SPxMainSM<R>::checkSolution(SPxLPBase<R>& lp, DVectorBase<R> sol)
  {
    for(int i = lp.nRows()-1; i >= 0; --i)
      {
        const SVectorBase<R>& row = lp.rowVector(i);
        R activity = 0;

        for(int k = 0; k < row.size(); k++)
          activity += row.value(k)*sol[row.index(k)];

        if(!GE(activity, lp.lhs(i), feastol()) || !LE(activity, lp.rhs(i), feastol()))
          return false;
      }

    return true;
  }

