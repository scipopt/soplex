

/// Returns unscaled lower bound of column \p i.
template<>
R SPxLPBase<R>::lowerUnscaled(int i) const
{
  assert(i >= 0 && i < nCols());
  if( _isScaled )
    return lp_scaler->lowerUnscaled(*this, i);
  else
    return LPColSetBase<R>::lower(i);
}

