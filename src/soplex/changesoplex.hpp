  template <class R>
void SPxSolverBase<R>::changeLower(int i, const R& newLower, bool scale)
{
  if( newLower != (scale ? lowerUnscaled<R>(i) : lower(i)) )
   {
      R oldLower = lower(i);
      // This has to be done before calling changeLowerStatus() because that is calling
      // basis.dualColStatus() which calls lower() and needs the changed value.
      SPxLPBase<R>::changeLower(i, newLower, scale);

      if (SPxBasisBase<R>::status() > SPxBasisBase<R>::NO_PROBLEM)
      {
         changeLowerStatus(i, lower(i), oldLower);
         unInit();
      }
   }
}
