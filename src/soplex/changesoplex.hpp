template <class R>
void SPxSolverBase<R>::changeUpper(int i, const R& newUpper, bool scale)
{
  if( newUpper != (scale ? this->upperUnscaled(i) : this->upper(i)) )
    {
      R oldUpper = this->upper(i);
      SPxLPBase<R>::changeUpper(i, newUpper, scale);

      if (SPxBasisBase<R>::status() > SPxBasisBase<R>::NO_PROBLEM)
        {
          changeUpperStatus(i, this->upper(i), oldUpper);
          unInit();
        }
    }
}
