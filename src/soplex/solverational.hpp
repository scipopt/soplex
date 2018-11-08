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

/// solves current problem with iterative refinement and recovery mechanism
template <class R>
void SoPlexBase<R>::_performOptIRStable(
                                        SolRational& sol,
                                        bool acceptUnbounded,
                                        bool acceptInfeasible,
                                        int minRounds,
                                        bool& primalFeasible,
                                        bool& dualFeasible,
                                        bool& infeasible,
                                        bool& unbounded,
                                        bool& stoppedTime,
                                        bool& stoppedIter,
                                        bool& error)
{
  // start rational solving timing
  _statistics->rationalTime->start();

  primalFeasible = false;
  dualFeasible = false;
  infeasible = false;
  unbounded = false;
  stoppedTime = false;
  stoppedIter = false;
  error = false;

  // set working tolerances in floating-point solver
  _solver.setFeastol(realParam(SoPlexBase<R>::FPFEASTOL));
  _solver.setOpttol(realParam(SoPlexBase<R>::FPOPTTOL));

  // declare vectors and variables
  typename SPxSolverBase<R>::Status result = SPxSolverBase<R>::UNKNOWN;

  _modLower.reDim(numColsRational(), false);
  _modUpper.reDim(numColsRational(), false);
  _modLhs.reDim(numRowsRational(), false);
  _modRhs.reDim(numRowsRational(), false);
  _modObj.reDim(numColsRational(), false);

  DVectorBase<R> primalReal(numColsRational());
  DVectorBase<R> dualReal(numRowsRational());

  Rational boundsViolation;
  Rational sideViolation;
  Rational redCostViolation;
  Rational dualViolation;
  Rational primalScale;
  Rational dualScale;
  Rational maxScale;

  // solve original LP
  MSG_INFO1( spxout, spxout << "Initial floating-point solve . . .\n" );

  if( _hasBasis )
    {
      assert(_basisStatusRows.size() == numRowsRational());
      assert(_basisStatusCols.size() == numColsRational());
      _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
      _hasBasis = (_solver.basis().status() > SPxBasisBase<R>::NO_PROBLEM);
    }

  for( int r = numRowsRational() - 1; r >= 0; r-- )
    {
      assert(_solver.maxRowObj(r) == 0.0);
    }

  _statistics->rationalTime->stop();
  result = _solveRealStable(acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols, _hasBasis);

  // evaluate result
  switch( result )
    {
    case SPxSolverBase<R>::OPTIMAL:
      MSG_INFO1( spxout, spxout << "Floating-point optimal.\n" );
      break;
    case SPxSolverBase<R>::INFEASIBLE:
      MSG_INFO1( spxout, spxout << "Floating-point infeasible.\n" );
      // the floating-point solve returns a Farkas ray if and only if the simplifier was not used, which is exactly
      // the case when a basis could be returned
      if( _hasBasis )
        {
          sol._dualFarkas = dualReal;
          sol._hasDualFarkas = true;
        }
      else
        sol._hasDualFarkas = false;
      infeasible = true;
      return;
    case SPxSolverBase<R>::UNBOUNDED:
      MSG_INFO1( spxout, spxout << "Floating-point unbounded.\n" );
      unbounded = true;
      return;
    case SPxSolverBase<R>::ABORT_TIME:
      stoppedTime = true;
      return;
    case SPxSolverBase<R>::ABORT_ITER:
      stoppedIter = true;
      return;
    default:
      error = true;
      return;
    }

  _statistics->rationalTime->start();

  // store floating-point solution of original LP as current rational solution and ensure that solution vectors have
  // right dimension; ensure that solution is aligned with basis
  sol._primal.reDim(numColsRational(), false);
  sol._slacks.reDim(numRowsRational(), false);
  sol._dual.reDim(numRowsRational(), false);
  sol._redCost.reDim(numColsRational(), false);
  sol._isPrimalFeasible= true;
  sol._isDualFeasible = true;

  for( int c = numColsRational() - 1; c >= 0; c-- )
    {
      typename SPxSolverBase<R>::VarStatus& basisStatusCol = _basisStatusCols[c];

      if( basisStatusCol == SPxSolverBase<R>::ON_LOWER )
        sol._primal[c] = lowerRational(c);
      else if( basisStatusCol == SPxSolverBase<R>::ON_UPPER )
        sol._primal[c] = upperRational(c);
      else if( basisStatusCol == SPxSolverBase<R>::FIXED )
        {
          // it may happen that lower and upper are only equal in the R LP but different in the rational LP; we do
          // not check this to avoid rational comparisons, but simply switch the basis status to the lower bound; this
          // is necessary, because for fixed variables any reduced cost is feasible
          sol._primal[c] = lowerRational(c);
          basisStatusCol = SPxSolverBase<R>::ON_LOWER;
        }
      else if( basisStatusCol == SPxSolverBase<R>::ZERO )
        sol._primal[c] = 0;
      else
        sol._primal[c] = primalReal[c];
    }
  _rationalLP->computePrimalActivity(sol._primal, sol._slacks);

  int dualSize = 0;
  for( int r = numRowsRational() - 1; r >= 0; r-- )
    {
      typename SPxSolverBase<R>::VarStatus& basisStatusRow = _basisStatusRows[r];

      // it may happen that left-hand and right-hand side are different in the rational, but equal in the R LP,
      // leading to a fixed basis status; this is critical because rows with fixed basis status are ignored in the
      // computation of the dual violation; to avoid rational comparisons we do not check this but simply switch to
      // the left-hand side status
      if( basisStatusRow == SPxSolverBase<R>::FIXED )
        basisStatusRow = SPxSolverBase<R>::ON_LOWER;

      {
        sol._dual[r] = dualReal[r];
        if( dualReal[r] != 0.0 )
          dualSize++;
      }
    }
  // we assume that the objective function vector has less nonzeros than the reduced cost vector, and so multiplying
  // with -1 first and subtracting the dual activity should be faster than adding the dual activity and negating
  // afterwards
  _rationalLP->getObj(sol._redCost);
  _rationalLP->subDualActivity(sol._dual, sol._redCost);

  // initial scaling factors are one
  primalScale = _rationalPosone;
  dualScale = _rationalPosone;

  // control progress
  Rational maxViolation;
  Rational bestViolation = _rationalPosInfty;
  const Rational violationImprovementFactor = 0.9;
  const Rational errorCorrectionFactor = 1.1;
  Rational errorCorrection = 2;
  int numFailedRefinements = 0;

  // store basis status in case solving modified problem failed
  DataArray< typename SPxSolverBase<R>::VarStatus > basisStatusRowsFirst;
  DataArray< typename SPxSolverBase<R>::VarStatus > basisStatusColsFirst;

  // refinement loop
  const bool maximizing = (intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE);
  const int maxDimRational = numColsRational() > numRowsRational() ? numColsRational() : numRowsRational();
  SolRational factorSol;
  bool factorSolNewBasis = true;
  int lastStallRefinements = 0;
  int nextRatrecRefinement = 0;
  do
    {
      // decrement minRounds counter
      minRounds--;

      MSG_DEBUG( std::cout << "Computing primal violations.\n" );

      // compute violation of bounds
      boundsViolation = 0;
      for( int c = numColsRational() - 1; c >= 0; c-- )
        {
          // lower bound
          assert((lowerRational(c) > _rationalNegInfty) == _lowerFinite(_colTypes[c]));
          if( _lowerFinite(_colTypes[c]) )
            {
              if( lowerRational(c) == 0 )
                {
                  _modLower[c] = sol._primal[c];
                  _modLower[c] *= -1;
                  if( _modLower[c] > boundsViolation )
                    boundsViolation = _modLower[c];
                }
              else
                {
                  _modLower[c] = lowerRational(c);
                  _modLower[c] -= sol._primal[c];
                  if( _modLower[c] > boundsViolation )
                    boundsViolation = _modLower[c];
                }
            }

          // upper bound
          assert((upperRational(c) < _rationalPosInfty) == _upperFinite(_colTypes[c]));
          if( _upperFinite(_colTypes[c]) )
            {
              if( upperRational(c) == 0 )
                {
                  _modUpper[c] = sol._primal[c];
                  _modUpper[c] *= -1;
                  if( _modUpper[c] < -boundsViolation )
                    boundsViolation = -_modUpper[c];
                }
              else
                {
                  _modUpper[c] = upperRational(c);
                  _modUpper[c] -= sol._primal[c];
                  if( _modUpper[c] < -boundsViolation )
                    boundsViolation = -_modUpper[c];
                }
            }
        }

      // compute violation of sides
      sideViolation = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
        {
          const typename SPxSolverBase<R>::VarStatus& basisStatusRow = _basisStatusRows[r];

          // left-hand side
          assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
          if( _lowerFinite(_rowTypes[r]) )
            {
              if( lhsRational(r) == 0 )
                {
                  _modLhs[r] = sol._slacks[r];
                  _modLhs[r] *= -1;
                }
              else
                {
                  _modLhs[r] = lhsRational(r);
                  _modLhs[r] -= sol._slacks[r];
                }

              if( _modLhs[r] > sideViolation )
                sideViolation = _modLhs[r];
              // if the activity is feasible, but too far from the bound, this violates complementary slackness; we
              // count it as side violation here
              else if( basisStatusRow == SPxSolverBase<R>::ON_LOWER && _modLhs[r] < -sideViolation )
                sideViolation = -_modLhs[r];
            }

          // right-hand side
          assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));
          if( _upperFinite(_rowTypes[r]) )
            {
              if( rhsRational(r) == 0 )
                {
                  _modRhs[r] = sol._slacks[r];
                  _modRhs[r] *= -1;
                }
              else
                {
                  _modRhs[r] = rhsRational(r);
                  _modRhs[r] -= sol._slacks[r];
                }

              if( _modRhs[r] < -sideViolation )
                sideViolation = -_modRhs[r];
              // if the activity is feasible, but too far from the bound, this violates complementary slackness; we
              // count it as side violation here
              else if( basisStatusRow == SPxSolverBase<R>::ON_UPPER && _modRhs[r] > sideViolation )
                sideViolation = _modRhs[r];
            }
        }

      MSG_DEBUG( std::cout << "Computing dual violations.\n" );

      // compute reduced cost violation
      redCostViolation = 0;
      for( int c = numColsRational() - 1; c >= 0; c-- )
        {
          if( _colTypes[c] == RANGETYPE_FIXED )
            continue;

          const typename SPxSolverBase<R>::VarStatus& basisStatusCol = _basisStatusCols[c];
          assert(basisStatusCol != SPxSolverBase<R>::FIXED);

          if( ((maximizing && basisStatusCol != SPxSolverBase<R>::ON_LOWER) || (!maximizing && basisStatusCol != SPxSolverBase<R>::ON_UPPER))
              && sol._redCost[c] < -redCostViolation )
            {
              MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
                         << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
                         << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
                         << ", sol._redCost[c] = " << rationalToString(sol._redCost[c])
                         << "\n" );
              redCostViolation = -sol._redCost[c];
            }

          if( ((maximizing && basisStatusCol != SPxSolverBase<R>::ON_UPPER) || (!maximizing && basisStatusCol != SPxSolverBase<R>::ON_LOWER))
              && sol._redCost[c] > redCostViolation )
            {
              MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
                         << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
                         << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
                         << ", sol._redCost[c] = " << rationalToString(sol._redCost[c])
                         << "\n" );
              redCostViolation = sol._redCost[c];
            }
        }

      // compute dual violation
      dualViolation = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
        {
          if( _rowTypes[r] == RANGETYPE_FIXED )
            continue;

          const typename SPxSolverBase<R>::VarStatus& basisStatusRow = _basisStatusRows[r];
          assert(basisStatusRow != SPxSolverBase<R>::FIXED);

          if( ((maximizing && basisStatusRow != SPxSolverBase<R>::ON_LOWER) || (!maximizing && basisStatusRow != SPxSolverBase<R>::ON_UPPER))
              && sol._dual[r] < -dualViolation )
            {
              MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                         << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                         << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                         << ", sol._dual[r] = " << rationalToString(sol._dual[r])
                         << "\n" );
              dualViolation = -sol._dual[r];
            }

          if( ((maximizing && basisStatusRow != SPxSolverBase<R>::ON_UPPER) || (!maximizing && basisStatusRow != SPxSolverBase<R>::ON_LOWER))
              && sol._dual[r] > dualViolation )
            {
              MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                         << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                         << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                         << ", sol._dual[r] = " << rationalToString(sol._dual[r])
                         << "\n" );
              dualViolation = sol._dual[r];
            }
        }

      _modObj = sol._redCost;

      // output violations; the reduced cost violations for artificially introduced slack columns are actually violations of the dual multipliers
      MSG_INFO1( spxout, spxout
                 << "Max. bound violation = " << rationalToString(boundsViolation) << "\n"
                 << "Max. row violation = " << rationalToString(sideViolation) << "\n"
                 << "Max. reduced cost violation = " << rationalToString(redCostViolation) << "\n"
                 << "Max. dual violation = " << rationalToString(dualViolation) << "\n" );

      MSG_DEBUG( spxout
                 << std::fixed << std::setprecision(2) << std::setw(10)
                 << "Progress table: "
                 << std::setw(10) << _statistics->refinements << " & "
                 << std::setw(10) << _statistics->iterations << " & "
                 << std::setw(10) << _statistics->solvingTime->time() << " & "
                 << std::setw(10) << _statistics->rationalTime->time() << " & "
                 << std::setw(10) << rationalToString(boundsViolation > sideViolation ? boundsViolation : sideViolation, 2) << " & "
                 << std::setw(10) << rationalToString(redCostViolation > dualViolation ? redCostViolation : dualViolation, 2) << "\n");

      // terminate if tolerances are satisfied
      primalFeasible = (boundsViolation <= _rationalFeastol && sideViolation <= _rationalFeastol);
      dualFeasible = (redCostViolation <= _rationalOpttol && dualViolation <= _rationalOpttol);
      if( primalFeasible && dualFeasible )
        {
          if( minRounds < 0 )
            {
              MSG_INFO1( spxout, spxout << "Tolerances reached.\n" );
              break;
            }
          else
            {
              MSG_INFO1( spxout, spxout << "Tolerances reached but minRounds forcing additional refinement rounds.\n" );
            }
        }

      // terminate if some limit is reached
      if( _isSolveStopped(stoppedTime, stoppedIter) )
        break;

      // check progress
      maxViolation = boundsViolation;
      if( sideViolation > maxViolation )
        maxViolation = sideViolation;
      if( redCostViolation > maxViolation )
        maxViolation = redCostViolation;
      if( dualViolation > maxViolation )
        maxViolation = dualViolation;
      bestViolation *= violationImprovementFactor;
      if( maxViolation > bestViolation )
        {
          MSG_INFO2( spxout, spxout << "Failed to reduce violation significantly.\n" );
          numFailedRefinements++;
        }
      else
        bestViolation = maxViolation;

      // decide whether to perform rational reconstruction and/or factorization
      bool performRatfac = boolParam(SoPlexBase<R>::RATFAC)
        && lastStallRefinements >= intParam(SoPlexBase<R>::RATFAC_MINSTALLS) && _hasBasis && factorSolNewBasis;
      bool performRatrec = boolParam(SoPlexBase<R>::RATREC)
        && (_statistics->refinements >= nextRatrecRefinement || performRatfac);

      // attempt rational reconstruction
      errorCorrection *= errorCorrectionFactor;
      if( performRatrec && maxViolation > 0 )
        {
          MSG_INFO1( spxout, spxout << "Performing rational reconstruction . . .\n" );

          maxViolation *= errorCorrection; // only used for sign check later
          maxViolation.invert();
          if( _reconstructSolutionRational(sol, _basisStatusRows, _basisStatusCols, maxViolation) )
            {
              MSG_INFO1( spxout, spxout << "Tolerances reached.\n" );
              primalFeasible = true;
              dualFeasible = true;
              break;
            }
          nextRatrecRefinement = int(_statistics->refinements * RealParam(SoPlexBase<R>::RATREC_FREQ)) + 1;
          MSG_DEBUG( spxout << "Next rational reconstruction after refinement " << nextRatrecRefinement << ".\n" );
        }

      // solve basis systems exactly
      if( performRatfac && maxViolation > 0 )
        {
          MSG_INFO1( spxout, spxout << "Performing rational factorization . . .\n" );

          bool optimal;
          _factorizeColumnRational(sol, _basisStatusRows, _basisStatusCols, stoppedTime, stoppedIter, error, optimal);
          factorSolNewBasis = false;

          if( stoppedTime )
            {
              MSG_INFO1( spxout, spxout << "Stopped rational factorization.\n" );
            }
          else if( error )
            {
              // message was already printed; reset error flag and continue without factorization
              error = false;
            }
          else if( optimal )
            {
              MSG_INFO1( spxout, spxout << "Tolerances reached.\n" );
              primalFeasible = true;
              dualFeasible = true;
              break;
            }
          else if( boolParam(SoPlexBase<R>::RATFACJUMP) )
            {
              MSG_INFO1( spxout, spxout << "Jumping to exact basic solution.\n" );
              minRounds++;
              continue;
            }
        }

      // start refinement

      // compute primal scaling factor; limit increase in scaling by tolerance used in floating point solve
      maxScale = primalScale;
      maxScale *= _rationalMaxscaleincr;

      primalScale = boundsViolation > sideViolation ? boundsViolation : sideViolation;
      if( primalScale < redCostViolation )
        primalScale = redCostViolation;
      assert(primalScale >= 0);

      if( primalScale > 0 )
        {
          primalScale.invert();
          if( primalScale > maxScale )
            primalScale = maxScale;
        }
      else
        primalScale = maxScale;

      if( boolParam(SoPlexBase<R>::POWERSCALING) )
        primalScale.powRound();

      // apply scaled bounds
      if( primalScale <= 1 )
        {
          if( primalScale < 1 )
            primalScale = 1;
          for( int c = numColsRational() - 1; c >= 0; c-- )
            {
              if( _lowerFinite(_colTypes[c]) )
                {
                  if( _modLower[c] <= _rationalNegInfty )
                    _solver.changeLower(c, -realParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeLower(c, R(_modLower[c]));
                }
              if( _upperFinite(_colTypes[c]) )
                {
                  if( _modUpper[c] >= _rationalPosInfty )
                    _solver.changeUpper(c, RealParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeUpper(c, R(_modUpper[c]));
                }
            }
        }
      else
        {
          MSG_INFO2( spxout, spxout << "Scaling primal by " << rationalToString(primalScale) << ".\n" );

          for( int c = numColsRational() - 1; c >= 0; c-- )
            {
              if( _lowerFinite(_colTypes[c]) )
                {
                  _modLower[c] *= primalScale;
                  if( _modLower[c] <= _rationalNegInfty )
                    _solver.changeLower(c, -realParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeLower(c, R(_modLower[c]));
                }
              if( _upperFinite(_colTypes[c]) )
                {
                  _modUpper[c] *= primalScale;
                  if( _modUpper[c] >= _rationalPosInfty )
                    _solver.changeUpper(c, RealParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeUpper(c, R(_modUpper[c]));
                }
            }
        }

      // apply scaled sides
      assert(primalScale >= 1);
      if( primalScale == 1 )
        {
          for( int r = numRowsRational() - 1; r >= 0; r-- )
            {
              if( _lowerFinite(_rowTypes[r]) )
                {
                  if( _modLhs[r] <= _rationalNegInfty )
                    _solver.changeLhs(r, -realParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeLhs(r, R(_modLhs[r]));
                }
              if( _upperFinite(_rowTypes[r]) )
                {
                  if( _modRhs[r] >= _rationalPosInfty )
                    _solver.changeRhs(r, RealParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeRhs(r, R(_modRhs[r]));
                }
            }
        }
      else
        {
          for( int r = numRowsRational() - 1; r >= 0; r-- )
            {
              if( _lowerFinite(_rowTypes[r]) )
                {
                  _modLhs[r] *= primalScale;
                  if( _modLhs[r] <= _rationalNegInfty )
                    _solver.changeLhs(r, -realParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeLhs(r, R(_modLhs[r]));
                }
              if( _upperFinite(_rowTypes[r]) )
                {
                  _modRhs[r] *= primalScale;
                  if( _modRhs[r] >= _rationalPosInfty )
                    _solver.changeRhs(r, RealParam(SoPlexBase<R>::INFTY));
                  else
                    _solver.changeRhs(r, R(_modRhs[r]));
                }
            }
        }

      // compute dual scaling factor; limit increase in scaling by tolerance used in floating point solve
      maxScale = dualScale;
      maxScale *= _rationalMaxscaleincr;

      dualScale = redCostViolation > dualViolation ? redCostViolation : dualViolation;
      assert(dualScale >= 0);

      if( dualScale > 0 )
        {
          dualScale.invert();
          if( dualScale > maxScale )
            dualScale = maxScale;
        }
      else
        dualScale = maxScale;

      if( boolParam(SoPlexBase<R>::POWERSCALING) )
        dualScale.powRound();

      if( dualScale > primalScale )
        dualScale = primalScale;

      if( dualScale < 1 )
        dualScale = 1;
      else
        {
          MSG_INFO2( spxout, spxout << "Scaling dual by " << rationalToString(dualScale) << ".\n" );

          // perform dual scaling
          ///@todo remove _modObj and use dualScale * sol._redCost directly
          _modObj *= dualScale;
        }

      // apply scaled objective function
      for( int c = numColsRational() - 1; c >= 0; c-- )
        {
          if( _modObj[c] >= _rationalPosInfty )
            _solver.changeObj(c, RealParam(SoPlexBase<R>::INFTY));
          else if( _modObj[c] <= _rationalNegInfty )
            _solver.changeObj(c, -realParam(SoPlexBase<R>::INFTY));
          else
            _solver.changeObj(c, R(_modObj[c]));
        }
      for( int r = numRowsRational() - 1; r >= 0; r-- )
        {
          Rational newRowObj;
          if( _rowTypes[r] == RANGETYPE_FIXED )
            _solver.changeRowObj(r, 0.0);
          else
            {
              newRowObj = sol._dual[r];
              newRowObj *= dualScale;
              if( newRowObj >= _rationalPosInfty )
                _solver.changeRowObj(r, -realParam(SoPlexBase<R>::INFTY));
              else if( newRowObj <= _rationalNegInfty )
                _solver.changeRowObj(r, RealParam(SoPlexBase<R>::INFTY));
              else
                _solver.changeRowObj(r, -Real(newRowObj));
            }
        }

      MSG_INFO1( spxout, spxout << "Refined floating-point solve . . .\n" );

      // ensure that artificial slack columns are basic and inequality constraints are nonbasic; otherwise we may end
      // up with dual violation on inequality constraints after removing the slack columns; do not change this in the
      // floating-point solver, though, because the solver may require its original basis to detect optimality
      if( _slackCols.num() > 0 && _hasBasis )
        {
          int numOrigCols = numColsRational() - _slackCols.num();
          assert(_slackCols.num() <= 0 || boolParam(SoPlexBase<R>::EQTRANS));
          for( int i = 0; i < _slackCols.num(); i++ )
            {
              int row = _slackCols.colVector(i).index(0);
              int col = numOrigCols + i;

              assert(row >= 0);
              assert(row < numRowsRational());

              if( _basisStatusRows[row] == SPxSolverBase<R>::BASIC && _basisStatusCols[col] != SPxSolverBase<R>::BASIC )
                {
                  _basisStatusRows[row] = _basisStatusCols[col];
                  _basisStatusCols[col] = SPxSolverBase<R>::BASIC;
                  _rationalLUSolver.clear();
                }
            }
        }

      // load basis
      if( _hasBasis && _solver.basis().status() < SPxBasisBase<R>::REGULAR )
        {
          MSG_DEBUG( spxout << "basis (status = " << _solver.basis().status() << ") desc before set:\n"; _solver.basis().desc().dump() );
          _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
          MSG_DEBUG( spxout << "basis (status = " << _solver.basis().status() << ") desc after set:\n"; _solver.basis().desc().dump() );

          _hasBasis = _solver.basis().status() > SPxBasisBase<R>::NO_PROBLEM;
          MSG_DEBUG( spxout << "setting basis in solver " << (_hasBasis ? "successful" : "failed") << " (3)\n" );
        }

      // solve modified problem
      int prevIterations = _statistics->iterations;
      _statistics->rationalTime->stop();
      result = _solveRealStable(acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols, _hasBasis, primalScale > 1e20 || dualScale > 1e20);

      // count refinements and remember whether we moved to a new basis
      _statistics->refinements++;
      if( _statistics->iterations <= prevIterations )
        {
          lastStallRefinements++;
          _statistics->stallRefinements++;
        }
      else
        {
          factorSolNewBasis = true;
          lastStallRefinements = 0;
          _statistics->pivotRefinements = _statistics->refinements;
        }

      // evaluate result; if modified problem was not solved to optimality, stop refinement
      switch( result )
        {
        case SPxSolverBase<R>::OPTIMAL:
          MSG_INFO1( spxout, spxout << "Floating-point optimal.\n" );
          break;
        case SPxSolverBase<R>::INFEASIBLE:
          MSG_INFO1( spxout, spxout << "Floating-point infeasible.\n" );
          sol._dualFarkas = dualReal;
          sol._hasDualFarkas = true;
          infeasible = true;
          _solver.clearRowObjs();
          return;
        case SPxSolverBase<R>::UNBOUNDED:
          MSG_INFO1( spxout, spxout << "Floating-point unbounded.\n" );
          unbounded = true;
          _solver.clearRowObjs();
          return;
        case SPxSolverBase<R>::ABORT_TIME:
          stoppedTime = true;
          return;
        case SPxSolverBase<R>::ABORT_ITER:
          stoppedIter = true;
          _solver.clearRowObjs();
          return;
        default:
          error = true;
          _solver.clearRowObjs();
          return;
        }

      _statistics->rationalTime->start();

      // correct primal solution and align with basis
      MSG_DEBUG( std::cout << "Correcting primal solution.\n" );

      int primalSize = 0;
      Rational primalScaleInverse = primalScale;
      primalScaleInverse.invert();
      _primalDualDiff.clear();
      for( int c = numColsRational() - 1; c >= 0; c-- )
        {
          // force values of nonbasic variables to bounds
          typename SPxSolverBase<R>::VarStatus& basisStatusCol = _basisStatusCols[c];

          if( basisStatusCol == SPxSolverBase<R>::ON_LOWER )
            {
              if( sol._primal[c] != lowerRational(c) )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = lowerRational(c);
                  _primalDualDiff.value(i) -= sol._primal[c];
                  sol._primal[c] = lowerRational(c);
                }
            }
          else if( basisStatusCol == SPxSolverBase<R>::ON_UPPER )
            {
              if( sol._primal[c] != upperRational(c) )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = upperRational(c);
                  _primalDualDiff.value(i) -= sol._primal[c];
                  sol._primal[c] = upperRational(c);
                }
            }
          else if( basisStatusCol == SPxSolverBase<R>::FIXED )
            {
              // it may happen that lower and upper are only equal in the R LP but different in the rational LP; we
              // do not check this to avoid rational comparisons, but simply switch the basis status to the lower
              // bound; this is necessary, because for fixed variables any reduced cost is feasible
              basisStatusCol = SPxSolverBase<R>::ON_LOWER;
              if( sol._primal[c] != lowerRational(c) )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = lowerRational(c);
                  _primalDualDiff.value(i) -= sol._primal[c];
                  sol._primal[c] = lowerRational(c);
                }
            }
          else if( basisStatusCol == SPxSolverBase<R>::ZERO )
            {
              if( sol._primal[c] != 0 )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = sol._primal[c];
                  _primalDualDiff.value(i) *= -1;
                  sol._primal[c] = 0;
                }
            }
          else
            {
              if( primalReal[c] == 1.0 )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = primalScaleInverse;
                  sol._primal[c] += _primalDualDiff.value(i);
                }
              else if( primalReal[c] == -1.0 )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = primalScaleInverse;
                  _primalDualDiff.value(i) *= -1;
                  sol._primal[c] += _primalDualDiff.value(i);
                }
              else if( primalReal[c] != 0.0 )
                {
                  int i = _primalDualDiff.size();
                  _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                  _primalDualDiff.add(c);
                  _primalDualDiff.value(i) = primalReal[c];
                  _primalDualDiff.value(i) *= primalScaleInverse;
                  sol._primal[c] += _primalDualDiff.value(i);
                }
            }

          if( sol._primal[c] != 0 )
            primalSize++;
        }

      // update or recompute slacks depending on which looks faster
      if( _primalDualDiff.size() < primalSize )
        {
          _rationalLP->addPrimalActivity(_primalDualDiff, sol._slacks);
#ifndef NDEBUG
#ifdef SOPLEX_WITH_GMP
          {
            DVectorRational activity(numRowsRational());
            _rationalLP->computePrimalActivity(sol._primal, activity);
            assert(sol._slacks == activity);
          }
#endif
#endif
        }
      else
        _rationalLP->computePrimalActivity(sol._primal, sol._slacks);
      const int numCorrectedPrimals = _primalDualDiff.size();

      // correct dual solution and align with basis
      MSG_DEBUG( std::cout << "Correcting dual solution.\n" );

#ifndef NDEBUG
      {
        // compute reduced cost violation
        DVectorRational debugRedCost(numColsRational());
        debugRedCost = DVectorRational(_realLP->maxObj());
        debugRedCost *= -1;
        _rationalLP->subDualActivity(DVectorRational(dualReal), debugRedCost);

        Rational debugRedCostViolation = 0;
        for( int c = numColsRational() - 1; c >= 0; c-- )
          {
            if( _colTypes[c] == RANGETYPE_FIXED )
              continue;

            const SPxSolverBase<R>::VarStatus& basisStatusCol = _basisStatusCols[c];
            assert(basisStatusCol != SPxSolverBase<R>::FIXED);

            if( ((maximizing && basisStatusCol != SPxSolverBase<R>::ON_LOWER) || (!maximizing && basisStatusCol != SPxSolverBase<R>::ON_UPPER))
                && debugRedCost[c] < -debugRedCostViolation )
              {
                MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
                           << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
                           << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
                           << ", obj[c] = " << _realLP->obj(c)
                           << ", debugRedCost[c] = " << rationalToString(debugRedCost[c])
                           << "\n" );
                debugRedCostViolation = -debugRedCost[c];
              }

            if( ((maximizing && basisStatusCol != SPxSolverBase<R>::ON_UPPER) || (!maximizing && basisStatusCol != SPxSolverBase<R>::ON_LOWER))
                && debugRedCost[c] > debugRedCostViolation )
              {
                MSG_DEBUG( std::cout << "basisStatusCol = " << basisStatusCol
                           << ", lower tight = " << bool(sol._primal[c] <= lowerRational(c))
                           << ", upper tight = " << bool(sol._primal[c] >= upperRational(c))
                           << ", obj[c] = " << _realLP->obj(c)
                           << ", debugRedCost[c] = " << rationalToString(debugRedCost[c])
                           << "\n" );
                debugRedCostViolation = debugRedCost[c];
              }
          }

        // compute dual violation
        Rational debugDualViolation = 0;
        Rational debugBasicDualViolation = 0;
        for( int r = numRowsRational() - 1; r >= 0; r-- )
          {
            if( _rowTypes[r] == RANGETYPE_FIXED )
              continue;

            const SPxSolverBase<R>::VarStatus& basisStatusRow = _basisStatusRows[r];
            assert(basisStatusRow != SPxSolverBase<R>::FIXED);

            Rational val = (-dualScale * sol._dual[r]) - Rational(dualReal[r]);

            if( ((maximizing && basisStatusRow != SPxSolverBase<R>::ON_LOWER) || (!maximizing && basisStatusRow != SPxSolverBase<R>::ON_UPPER))
                && val > debugDualViolation )
              {
                MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                           << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                           << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                           << ", dualReal[r] = " << rationalToString(val)
                           << ", dualReal[r] = " << dualReal[r]
                           << "\n" );
                debugDualViolation = val;
              }

            if( ((maximizing && basisStatusRow != SPxSolverBase<R>::ON_UPPER) || (!maximizing && basisStatusRow != SPxSolverBase<R>::ON_LOWER))
                && val < -debugDualViolation )
              {
                MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                           << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                           << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                           << ", dualReal[r] = " << rationalToString(val)
                           << ", dualReal[r] = " << dualReal[r]
                           << "\n" );
                debugDualViolation = -val;
              }

            if( basisStatusRow == SPxSolverBase<R>::BASIC && spxAbs(val) > debugBasicDualViolation )
              {
                MSG_DEBUG( std::cout << "basisStatusRow = " << basisStatusRow
                           << ", lower tight = " << bool(sol._slacks[r] <= lhsRational(r))
                           << ", upper tight = " << bool(sol._slacks[r] >= rhsRational(r))
                           << ", dualReal[r] = " << rationalToString(val)
                           << ", dualReal[r] = " << dualReal[r]
                           << "\n" );
                debugBasicDualViolation = spxAbs(val);
              }
          }

        if( debugRedCostViolation > _solver.opttol() || debugDualViolation > _solver.opttol() || debugBasicDualViolation > 1e-9 )
          {
            MSG_WARNING( spxout, spxout << "Warning: floating-point dual solution with violation "
                         << rationalToString(debugRedCostViolation) << " / "
                         << rationalToString(debugDualViolation) << " / "
                         << rationalToString(debugBasicDualViolation)
                         << " (red. cost, dual, basic).\n" );
          }
      }
#endif

      Rational dualScaleInverseNeg = dualScale;
      dualScaleInverseNeg.invert();
      dualScaleInverseNeg *= -1;
      _primalDualDiff.clear();
      dualSize = 0;
      for( int r = numRowsRational() - 1; r >= 0; r-- )
        {
          typename SPxSolverBase<R>::VarStatus& basisStatusRow = _basisStatusRows[r];

          // it may happen that left-hand and right-hand side are different in the rational, but equal in the R LP,
          // leading to a fixed basis status; this is critical because rows with fixed basis status are ignored in the
          // computation of the dual violation; to avoid rational comparisons we do not check this but simply switch
          // to the left-hand side status
          if( basisStatusRow == SPxSolverBase<R>::FIXED )
            basisStatusRow = SPxSolverBase<R>::ON_LOWER;

          {
            if( dualReal[r] != 0 )
              {
                int i = _primalDualDiff.size();
                _ensureDSVectorRationalMemory(_primalDualDiff, maxDimRational);
                _primalDualDiff.add(r);
                _primalDualDiff.value(i) = dualReal[r];
                _primalDualDiff.value(i) *= dualScaleInverseNeg;
                sol._dual[r] -= _primalDualDiff.value(i);

                dualSize++;
              }
            else
              {
                // we do not check whether the dual value is nonzero, because it probably is; this gives us an
                // overestimation of the number of nonzeros in the dual solution
                dualSize++;
              }
          }
        }

      // update or recompute reduced cost values depending on which looks faster; adding one to the length of the
      // dual vector accounts for the objective function vector
      if( _primalDualDiff.size() < dualSize + 1 )
        {
          _rationalLP->addDualActivity(_primalDualDiff, sol._redCost);
#ifndef NDEBUG
          {
            DVectorRational activity(_rationalLP->maxObj());
            activity *= -1;
            _rationalLP->subDualActivity(sol._dual, activity);
          }
#endif
        }
      else
        {
          // we assume that the objective function vector has less nonzeros than the reduced cost vector, and so multiplying
          // with -1 first and subtracting the dual activity should be faster than adding the dual activity and negating
          // afterwards
          _rationalLP->getObj(sol._redCost);
          _rationalLP->subDualActivity(sol._dual, sol._redCost);
        }
      const int numCorrectedDuals = _primalDualDiff.size();

      if( numCorrectedPrimals + numCorrectedDuals > 0 )
        {
          MSG_INFO2( spxout, spxout << "Corrected " << numCorrectedPrimals << " primal variables and " << numCorrectedDuals << " dual values.\n" );
        }
    }
  while( true );

  // correct basis status for restricted inequalities
  if( _hasBasis )
    {
      for( int r = numRowsRational() - 1; r >= 0; r-- )
        {
          assert((lhsRational(r) == rhsRational(r)) == (_rowTypes[r] == RANGETYPE_FIXED));
          if( _rowTypes[r] != RANGETYPE_FIXED && _basisStatusRows[r] == SPxSolverBase<R>::FIXED )
            _basisStatusRows[r] = (maximizing == (sol._dual[r] < 0))
              ? SPxSolverBase<R>::ON_LOWER
              : SPxSolverBase<R>::ON_UPPER;
        }
    }

  // compute objective function values
  assert(sol._isPrimalFeasible == sol._isDualFeasible);
  if( sol._isPrimalFeasible)
    {
      sol._objVal = sol._primal * _rationalLP->maxObj();
      if( intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MINIMIZE )
        sol._objVal *= -1;
    }

  // set objective coefficients for all rows to zero
  _solver.clearRowObjs();

  // stop rational solving time
  _statistics->rationalTime->stop();
}
