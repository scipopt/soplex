// General specializations
/// solves real LP during iterative refinement
template <class R>
typename SPxSolverBase<R>::Status SoPlexBase<R>::_solveRealForRational(bool fromscratch, VectorBase<R>& primal, VectorBase<R>& dual, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols, bool& returnedBasis)
{
  assert(_isConsistent());

  assert(_solver.nRows() == numRowsRational());
  assert(_solver.nCols() == numColsRational());
  assert(primal.dim() == numColsRational());
  assert(dual.dim() == numRowsRational());

  typename SPxSolverBase<R>::Status result = SPxSolverBase<R>::UNKNOWN;

#ifndef SOPLEX_MANUAL_ALT
  if( fromscratch || !_hasBasis )
    _enableSimplifierAndScaler();
  else
    _disableSimplifierAndScaler();
#else
  _disableSimplifierAndScaler();
#endif

  // start timing
  _statistics->syncTime->start();

  // if preprocessing is applied, we need to restore the original LP at the end
  SPxLPRational* rationalLP = 0;
  if( _simplifier != 0 || _scaler != 0 )
    {
      spx_alloc(rationalLP);
      rationalLP = new (rationalLP) SPxLPRational(_solver);
    }

  // if preprocessing is applied, the basis may change, hence invalidate the rational basis factorization; if no
  if( _simplifier != 0 || _scaler != 0 )
    _rationalLUSolver.clear();

  // stop timing
  _statistics->syncTime->stop();

  returnedBasis = false;
  try
    {
      // apply problem simplification
      typename SPxSimplifier<R>::Result simplificationStatus = SPxSimplifier<R>::OKAY;
      if( _simplifier != 0 )
        {
          // do not remove bounds of boxed variables or sides of ranged rows if bound flipping is used
          bool keepbounds = intParam(SoPlexBase<R>::RATIOTESTER) == SoPlexBase<R>::RATIOTESTER_BOUNDFLIPPING;
          simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlexBase<R>::EPSILON_ZERO), realParam(SoPlexBase<R>::FPFEASTOL), realParam(SoPlexBase<R>::FPOPTTOL), keepbounds);
        }

      // apply scaling after the simplification
      if( _scaler != 0 && simplificationStatus == SPxSimplifier<R>::OKAY )
        _scaler->scale(_solver, false);

      // run the simplex method if problem has not been solved by the simplifier
      if( simplificationStatus == SPxSimplifier<R>::OKAY )
        {
          MSG_INFO1( spxout, spxout << std::endl );

          _solveRealLPAndRecordStatistics();

          MSG_INFO1( spxout, spxout << std::endl );
        }

      ///@todo move to private helper methods
      // evaluate status flag
      if( simplificationStatus == SPxSimplifier<R>::INFEASIBLE )
        result = SPxSolverBase<R>::INFEASIBLE;
      else if( simplificationStatus == SPxSimplifier<R>::DUAL_INFEASIBLE )
        result = SPxSolverBase<R>::INForUNBD;
      else if( simplificationStatus == SPxSimplifier<R>::UNBOUNDED )
        result = SPxSolverBase<R>::UNBOUNDED;
      else if( simplificationStatus == SPxSimplifier<R>::VANISHED || simplificationStatus == SPxSimplifier<R>::OKAY )
        {
          result = simplificationStatus == SPxSimplifier<R>::VANISHED ? SPxSolverBase<R>::OPTIMAL : _solver.status();

          ///@todo move to private helper methods
          // process result
          switch( result )
            {
            case SPxSolverBase<R>::OPTIMAL:
              // unsimplify if simplifier is active and LP is solved to optimality; this must be done here and not at solution
              // query, because we want to have the basis for the original problem
              if( _simplifier != 0 )
                {
                  assert(!_simplifier->isUnsimplified());
                  assert(simplificationStatus == SPxSimplifier<R>::VANISHED || simplificationStatus == SPxSimplifier<R>::OKAY);

                  bool vanished = simplificationStatus == SPxSimplifier<R>::VANISHED;

                  // get solution vectors for transformed problem
                  DVectorReal tmpPrimal(vanished ? 0 : _solver.nCols());
                  DVectorReal tmpSlacks(vanished ? 0 : _solver.nRows());
                  DVectorReal tmpDual(vanished ? 0 : _solver.nRows());
                  DVectorReal tmpRedCost(vanished ? 0 : _solver.nCols());

                  if( !vanished )
                    {
                      assert(_solver.status() == SPxSolverBase<R>::OPTIMAL);

                      _solver.getPrimalSol(tmpPrimal);
                      _solver.getSlacks(tmpSlacks);
                      _solver.getDualSol(tmpDual);
                      _solver.getRedCostSol(tmpRedCost);

                      // unscale vectors
                      if( _scaler != 0 )
                        {
                          _scaler->unscalePrimal(_solver, tmpPrimal);
                          _scaler->unscaleSlacks(_solver, tmpSlacks);
                          _scaler->unscaleDual(_solver, tmpDual);
                          _scaler->unscaleRedCost(_solver, tmpRedCost);
                        }

                      // get basis of transformed problem
                      basisStatusRows.reSize(_solver.nRows());
                      basisStatusCols.reSize(_solver.nCols());
                      _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr(), basisStatusRows.size(), basisStatusCols.size());
                    }

                  ///@todo catch exception
                  _simplifier->unsimplify(tmpPrimal, tmpDual, tmpSlacks, tmpRedCost, basisStatusRows.get_ptr(), basisStatusCols.get_ptr());

                  // store basis for original problem
                  basisStatusRows.reSize(numRowsRational());
                  basisStatusCols.reSize(numColsRational());
                  _simplifier->getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr(), basisStatusRows.size(), basisStatusCols.size());
                  returnedBasis = true;

                  primal = _simplifier->unsimplifiedPrimal();
                  dual = _simplifier->unsimplifiedDual();
                }
              // if the original problem is not in the solver because of scaling, we also need to store the basis
              else
                {
                  _solver.getPrimalSol(primal);
                  _solver.getDualSol(dual);

                  // unscale vectors
                  if( _scaler != 0 )
                    {
                      _scaler->unscalePrimal(_solver, primal);
                      _scaler->unscaleDual(_solver, dual);
                    }

                  // get basis of transformed problem
                  basisStatusRows.reSize(_solver.nRows());
                  basisStatusCols.reSize(_solver.nCols());
                  _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr(), basisStatusRows.size(), basisStatusCols.size());
                  returnedBasis = true;
                }
              break;

            case SPxSolverBase<R>::ABORT_CYCLING:
              if( _simplifier == 0 && boolParam(SoPlexBase<R>::ACCEPTCYCLING) )
                {
                  _solver.getPrimalSol(primal);
                  _solver.getDualSol(dual);

                  // unscale vectors
                  if( _scaler != 0 )
                    {
                      _scaler->unscalePrimal(_solver, primal);
                      _scaler->unscaleDual(_solver, dual);
                    }

                  // get basis of transformed problem
                  basisStatusRows.reSize(_solver.nRows());
                  basisStatusCols.reSize(_solver.nCols());
                  _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr(), basisStatusRows.size(), basisStatusCols.size());
                  returnedBasis = true;
                  result = SPxSolverBase<R>::OPTIMAL;
                }
              break;
            case SPxSolverBase<R>::ABORT_TIME:
            case SPxSolverBase<R>::ABORT_ITER:
            case SPxSolverBase<R>::ABORT_VALUE:
            case SPxSolverBase<R>::REGULAR:
            case SPxSolverBase<R>::RUNNING:
            case SPxSolverBase<R>::UNBOUNDED:
              break;
            case SPxSolverBase<R>::INFEASIBLE:
              // if simplifier is active we cannot return a Farkas ray currently
              if( _simplifier != 0 )
                break;

              // return Farkas ray as dual solution
              _solver.getDualfarkas(dual);

              // unscale vectors
              if( _scaler != 0 )
                _scaler->unscaleDual(_solver, dual);

              // if the original problem is not in the solver because of scaling, we also need to store the basis
              basisStatusRows.reSize(_solver.nRows());
              basisStatusCols.reSize(_solver.nCols());
              _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr(), basisStatusRows.size(), basisStatusCols.size());
              returnedBasis = true;
              break;

            case SPxSolverBase<R>::INForUNBD:
            case SPxSolverBase<R>::SINGULAR:
            default:
              break;
            }
        }
    }
  catch( ... )
    {
      MSG_INFO1( spxout, spxout << "Exception thrown during floating-point solve.\n" );
      result = SPxSolverBase<R>::ERROR;
    }

  // restore original LP if necessary
  if( _simplifier != 0 || _scaler != 0 )
    {
      assert(rationalLP != 0);
      _solver.loadLP((SPxLPBase<R>)(*rationalLP));
      rationalLP->~SPxLPRational();
      spx_free(rationalLP);
      if( _hasBasis )
        _solver.setBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());
    }

  return result;
}

/// solves real LP with recovery mechanism
template <class R>
typename SPxSolverBase<R>::Status SoPlexBase<R>::_solveRealStable(bool acceptUnbounded, bool acceptInfeasible, VectorBase<R>& primal, VectorBase<R>& dual, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols, bool& returnedBasis, const bool forceNoSimplifier)
{
  typename SPxSolverBase<R>::Status result = SPxSolverBase<R>::UNKNOWN;

  bool fromScratch = false;
  bool solved = false;
  bool solvedFromScratch = false;
  bool initialSolve = true;
  bool increasedMarkowitz = false;
  bool relaxedTolerances = false;
  bool tightenedTolerances = false;
  bool switchedScaler = false;
  bool switchedSimplifier = false;
  bool switchedRatiotester = false;
  bool switchedPricer = false;
  bool turnedoffPre = false;

  Real markowitz = _slufactor.markowitz();
  int ratiotester = intParam(SoPlexBase<R>::RATIOTESTER);
  int pricer = intParam(SoPlexBase<R>::PRICER);
  int simplifier = intParam(SoPlexBase<R>::SIMPLIFIER);
  int scaler = intParam(SoPlexBase<R>::SCALER);
  int type = intParam(SoPlexBase<R>::ALGORITHM);

  if( forceNoSimplifier )
    setIntParam(SoPlexBase<R>::SIMPLIFIER, SoPlexBase<R>::SIMPLIFIER_OFF);

  while( true )
    {
      assert(!increasedMarkowitz || GE(_slufactor.markowitz(), 0.9));

      result = _solveRealForRational(fromScratch, primal, dual, basisStatusRows, basisStatusCols, returnedBasis);

      solved = (result == SPxSolverBase<R>::OPTIMAL)
        || (result == SPxSolverBase<R>::INFEASIBLE && acceptInfeasible)
        || (result == SPxSolverBase<R>::UNBOUNDED && acceptUnbounded);

      if( solved )
        break;

      //         if( _isSolveStopped() )
      //            break;

      if( initialSolve )
        {
          MSG_INFO1( spxout, spxout << "Numerical troubles during floating-point solve." << std::endl );
          initialSolve = false;
        }

      if( !turnedoffPre
          && (intParam(SoPlexBase<R>::SIMPLIFIER) != SoPlexBase<R>::SIMPLIFIER_OFF || intParam(SoPlexBase<R>::SCALER) != SoPlexBase<R>::SCALER_OFF) )
        {
          MSG_INFO1( spxout, spxout << "Turning off preprocessing." << std::endl );

          turnedoffPre = true;

          setIntParam(SoPlexBase<R>::SCALER, SoPlexBase<R>::SCALER_OFF);
          setIntParam(SoPlexBase<R>::SIMPLIFIER, SoPlexBase<R>::SIMPLIFIER_OFF);

          fromScratch = true;
          _solver.reLoad();
          solvedFromScratch = true;
          continue;
        }

      setIntParam(SoPlexBase<R>::SCALER, scaler);
      setIntParam(SoPlexBase<R>::SIMPLIFIER, simplifier);

      if( !increasedMarkowitz )
        {
          MSG_INFO1( spxout, spxout << "Increasing Markowitz threshold." << std::endl );

          _slufactor.setMarkowitz(0.9);
          increasedMarkowitz = true;
          try
            {
              _solver.factorize();
              continue;
            }
          catch( ... )
            {
              MSG_DEBUG( std::cout << std::endl << "Factorization failed." << std::endl );
            }
        }

      if( !solvedFromScratch )
        {
          MSG_INFO1( spxout, spxout << "Solving from scratch." << std::endl );

          fromScratch = true;
          _solver.reLoad();

          solvedFromScratch = true;
          continue;
        }

      setIntParam(SoPlexBase<R>::RATIOTESTER, ratiotester);
      setIntParam(SoPlexBase<R>::PRICER, pricer);

      if( !switchedScaler )
        {
          MSG_INFO1( spxout, spxout << "Switching scaling." << std::endl );

          if( scaler == int(SoPlexBase<R>::SCALER_OFF) )
            setIntParam(SoPlexBase<R>::SCALER, SoPlexBase<R>::SCALER_BIEQUI);
          else
            setIntParam(SoPlexBase<R>::SCALER, SoPlexBase<R>::SCALER_OFF);

          fromScratch = true;
          _solver.reLoad();

          solvedFromScratch = true;
          switchedScaler = true;
          continue;
        }

      if( !switchedSimplifier && !forceNoSimplifier )
        {
          MSG_INFO1( spxout, spxout << "Switching simplification." << std::endl );

          if( simplifier == int(SoPlexBase<R>::SIMPLIFIER_OFF) )
            setIntParam(SoPlexBase<R>::SIMPLIFIER, SoPlexBase<R>::SIMPLIFIER_AUTO);
          else
            setIntParam(SoPlexBase<R>::SIMPLIFIER, SoPlexBase<R>::SIMPLIFIER_OFF);

          fromScratch = true;
          _solver.reLoad();

          solvedFromScratch = true;
          switchedSimplifier = true;
          continue;
        }

      setIntParam(SoPlexBase<R>::SIMPLIFIER, SoPlexBase<R>::SIMPLIFIER_OFF);

      if( !relaxedTolerances )
        {
          MSG_INFO1( spxout, spxout << "Relaxing tolerances." << std::endl );

          setIntParam(SoPlexBase<R>::ALGORITHM, SoPlexBase<R>::ALGORITHM_PRIMAL);
          _solver.setDelta((_solver.feastol() * 1e3 > 1e-3) ? 1e-3 : (_solver.feastol() * 1e3));
          relaxedTolerances = _solver.feastol() >= 1e-3;
          solvedFromScratch = false;
          continue;
        }

      if( !tightenedTolerances && result != SPxSolverBase<R>::INFEASIBLE )
        {
          MSG_INFO1( spxout, spxout << "Tightening tolerances." << std::endl );

          setIntParam(SoPlexBase<R>::ALGORITHM, SoPlexBase<R>::ALGORITHM_DUAL);
          _solver.setDelta(_solver.feastol() * 1e-3 < 1e-9 ? 1e-9 : _solver.feastol() * 1e-3);
          tightenedTolerances = _solver.feastol() <= 1e-9;
          solvedFromScratch = false;
          continue;
        }

      setIntParam(SoPlexBase<R>::ALGORITHM, type);

      if( !switchedRatiotester )
        {
          MSG_INFO1( spxout, spxout << "Switching ratio test." << std::endl );

          _solver.setType(_solver.type() == SPxSolverBase<R>::LEAVE ? SPxSolverBase<R>::ENTER : SPxSolverBase<R>::LEAVE);
          if( _solver.ratiotester() != (SPxRatioTester<R>*)&_ratiotesterTextbook )
            setIntParam(SoPlexBase<R>::RATIOTESTER, RATIOTESTER_TEXTBOOK);
          else
            setIntParam(SoPlexBase<R>::RATIOTESTER, RATIOTESTER_FAST);
          switchedRatiotester = true;
          solvedFromScratch = false;
          continue;
        }

      if( !switchedPricer )
        {
          MSG_INFO1( spxout, spxout << "Switching pricer." << std::endl );

          _solver.setType(_solver.type() == SPxSolverBase<R>::LEAVE ? SPxSolverBase<R>::ENTER : SPxSolverBase<R>::LEAVE);
          if( _solver.pricer() != (SPxPricer<R>*)&_pricerDevex )
            setIntParam(SoPlexBase<R>::PRICER, PRICER_DEVEX);
          else
            setIntParam(SoPlexBase<R>::PRICER, PRICER_STEEP);
          switchedPricer = true;
          solvedFromScratch = false;
          continue;
        }

      MSG_INFO1( spxout, spxout << "Giving up." << std::endl );

      break;
    }

  _solver.setFeastol(realParam(SoPlexBase<R>::FPFEASTOL));
  _solver.setOpttol(realParam(SoPlexBase<R>::FPOPTTOL));
  _slufactor.setMarkowitz(markowitz);

  setIntParam(SoPlexBase<R>::RATIOTESTER, ratiotester);
  setIntParam(SoPlexBase<R>::PRICER, pricer);
  setIntParam(SoPlexBase<R>::SIMPLIFIER, simplifier);
  setIntParam(SoPlexBase<R>::SCALER, scaler);
  setIntParam(SoPlexBase<R>::ALGORITHM, type);

  return result;
}
