//  LAST EDIT: Thu Oct 12 14:54:23 2000 by Tobias Pfender (opt32!bzfpfend) 
//-----------------------------------------------------------------------------
// name:     OSI Interface for SOPLEX
// author:   Tobias Pfender
//           Konrad-Zuse-Zentrum Berlin (Germany)
//           email: pfender@zib.de
// date:     01/17/2002
//-----------------------------------------------------------------------------
// Copyright (C) 2002, Tobias Pfender, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifdef COIN_USE_SPX
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iostream>

#include <string>

#include "OsiSpxSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiPackedMatrix.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif

void OsiSpxSolverInterface::printBounds()
{
  int nc = getNumCols();
  int nr = getNumRows();
  const char * s = getRowSense();
  const double * b = getRightHandSide();
  const double * rng = getRowRange();
  const double * cl = getColLower();
  const double * cu = getColUpper();
  const double * rl = getRowLower();
  const double * ru = getRowUpper();
  
  std::cout << "ncols=" << nc << ", nrows=" << nr;
  std::cout << std::endl << "sns=";
  int i;
  for( i = 0; i < nr; ++i )
    std::cout << " " << s[i];
  std::cout << std::endl << "rhs=";
  for( i = 0; i < nr; ++i )
    std::cout << " " << b[i];
  std::cout << std::endl << "rng=";
  for( i = 0; i < nr; ++i )
    std::cout << " " << rng[i];
  std::cout << std::endl << "cl =";
  for( i = 0; i < nc; ++i )
    std::cout << " " << cl[i];
  std::cout << std::endl << "cu =";
  for( i = 0; i < nc; ++i )
    std::cout << " " << cu[i];
  std::cout << std::endl << "rl =";
  for( i = 0; i < nr; ++i )
    std::cout << " " << rl[i];
  std::cout << std::endl << "ru =";
  for( i = 0; i < nr; ++i )
    std::cout << " " << ru[i];
  std::cout << std::endl;
}

//--------------------------------------------------------------------------
void OsiSpxSolverInterfaceUnitTest( const std::string & mpsDir )
{
  // Test default constructor
  {
    OsiSpxSolverInterface m;
    assert( m.rowsense_==NULL );
    assert( m.rhs_==NULL );
    assert( m.rowrange_==NULL );
    assert( m.colsol_==NULL );
    assert( m.rowsol_==NULL );
    assert( m.matrixByRow_==NULL );
    assert( m.matrixByCol_==NULL );
    assert( m.getApplicationData() == NULL );
    int i=2346;
    m.setApplicationData(&i);
    assert( *((int *)(m.getApplicationData())) == i );
  }

  {    
    OsiRelFltEq eq;
    OsiSpxSolverInterface m;
    std::string fn = mpsDir+"exmip1";
    m.readMps(fn.c_str(),"mps");
    int ad = 13579;
    m.setApplicationData(&ad);
    assert( *((int *)(m.getApplicationData())) == ad );

    {
      assert( m.getNumCols()==8 );
      const OsiPackedMatrix * colCopy = m.getMatrixByCol();
      assert( colCopy->getNumCols() == 8 );
      assert( colCopy->getMajorDim() == 8 );
      assert( colCopy->getNumRows() == 5 );
      assert( colCopy->getMinorDim() == 5 );
      assert (colCopy->getVectorLengths()[7] == 2 );
      OsiPackedMatrix revColCopy;
      revColCopy.reverseOrderedCopyOf(*colCopy);
      OsiPackedMatrix rev2ColCopy;      
      rev2ColCopy.reverseOrderedCopyOf(revColCopy);
      assert( rev2ColCopy.getNumCols() == 8 );
      assert( rev2ColCopy.getMajorDim() == 8 );
      assert( rev2ColCopy.getNumRows() == 5 );
      assert( rev2ColCopy.getMinorDim() == 5 );
      assert( rev2ColCopy.getVectorLengths()[7] == 2 );
    }
    
    {
      OsiSpxSolverInterface im;    
      assert( im.getNumCols() == 0 ); 
    }
    
    // Test copy constructor and assignment operator
    {
      OsiSpxSolverInterface lhs;
      {      
        assert( *((int *)(m.getApplicationData())) == ad );
        OsiSpxSolverInterface im(m);   
        assert( *((int *)(im.getApplicationData())) == ad );

        OsiSpxSolverInterface imC1(im);
        assert( imC1.getNumCols() == im.getNumCols() );
        assert( imC1.getNumRows() == im.getNumRows() );   
        assert( *((int *)(imC1.getApplicationData())) == ad ); 
        
        //im.setModelPtr(m);
        
        OsiSpxSolverInterface imC2(im);
        assert( imC2.getNumCols() == im.getNumCols() );
        assert( imC2.getNumRows() == im.getNumRows() );  
        assert( *((int *)(imC2.getApplicationData())) == ad ); 
        
        lhs=imC2;
      }
      // Test that lhs has correct values even though rhs has gone out of scope

      assert( lhs.getNumCols() == m.getNumCols() );
      assert( lhs.getNumRows() == m.getNumRows() );      
      assert( *((int *)(lhs.getApplicationData())) == ad );
    }
    
    // Test clone
    {
      OsiSpxSolverInterface soplexSi(m);
      OsiSolverInterface * siPtr = &soplexSi;
      OsiSolverInterface * siClone = siPtr->clone();
      OsiSpxSolverInterface * soplexClone = dynamic_cast<OsiSpxSolverInterface*>(siClone);
      assert( soplexClone != NULL );
      assert( soplexClone->getNumRows() == soplexSi.getNumRows() );
      assert( soplexClone->getNumCols() == m.getNumCols() );
      
      assert( *((int *)(soplexClone->getApplicationData())) == ad );
      delete siClone;
    }
   
    // test infinity
    {
      OsiSpxSolverInterface si;
      assert( eq( si.getInfinity(), soplex::infinity ) );
    }     
    
    // Test setting solution
    {
      OsiSpxSolverInterface m1(m);
      int i;

      double * cs = new double[m1.getNumCols()];
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        cs[i] = i + .5;
      m1.setColSolution(cs);
      for ( i = 0;  i < m1.getNumCols();  i++ ) 
        assert(m1.getColSolution()[i] == i + .5);
      
      double * rs = new double[m1.getNumRows()];
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        rs[i] = i - .5;
      m1.setRowPrice(rs);
      for ( i = 0;  i < m1.getNumRows();  i++ ) 
        assert(m1.getRowPrice()[i] == i - .5);

      delete [] cs;
      delete [] rs;
    }
    
    // Test fraction Indices
    {
      OsiSpxSolverInterface fim;
      std::string fn = mpsDir+"exmip1";
      fim.readMps(fn.c_str(),"mps");
      //fim.setModelPtr(m);
      // exmip1.mps has 2 integer variables with index 2 & 3
      assert(  fim.isContinuous(0) );
      assert(  fim.isContinuous(1) );
      assert( !fim.isContinuous(2) );
      assert( !fim.isContinuous(3) );
      assert(  fim.isContinuous(4) );
      
      assert( !fim.isInteger(0) );
      assert( !fim.isInteger(1) );
      assert(  fim.isInteger(2) );
      assert(  fim.isInteger(3) );
      assert( !fim.isInteger(4) );
      
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert(  fim.isBinary(2) );
      assert(  fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert( !fim.isIntegerNonBinary(2) );
      assert( !fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );

      
      // Test fractionalIndices
      {
	// Set a solution vector
	double * cs = new double[fim.getNumCols()];
	for ( int i = 0;  i < fim.getNumCols();  cs[i++] = 0.0 );
	cs[2] = 2.9;
	cs[3] = 3.0;
	fim.setColSolution(cs);

        OsiVectorInt fi = fim.getFractionalIndices();
        assert( fi.size() == 1 );
        assert( fi[0]==2 );
        
        // Set integer variables very close to integer values
        cs[2] = 5 + .00001/2.;
        cs[3] = 8 - .00001/2.;
	fim.setColSolution(cs);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 0 );
        
        // Set integer variables close, but beyond tolerances
        cs[2] = 5 + .00001*2.;
        cs[3] = 8 - .00001*2.;
	fim.setColSolution(cs);
        fi = fim.getFractionalIndices(1e-5);
        assert( fi.size() == 2 );
        assert( fi[0]==2 );
        assert( fi[1]==3 );

	delete [] cs;
      }
     
      // Change data so column 2 & 3 are integerNonBinary
      fim.setColUpper(2, 5);
      fim.setColUpper(3, 6.0);
      assert( !fim.isBinary(0) );
      assert( !fim.isBinary(1) );
      assert( !fim.isBinary(2) );
      assert( !fim.isBinary(3) );
      assert( !fim.isBinary(4) );
      
      assert( !fim.isIntegerNonBinary(0) );
      assert( !fim.isIntegerNonBinary(1) );
      assert(  fim.isIntegerNonBinary(2) );
      assert(  fim.isIntegerNonBinary(3) );
      assert( !fim.isIntegerNonBinary(4) );
    }
    
    // Test apply cuts method
    {      
      OsiSpxSolverInterface im(m);
      OsiCuts cuts;
      
      // Generate some cuts 
      {
        // Get number of rows and columns in model
        int nr=im.getNumRows();
        int nc=im.getNumCols();
        assert( nr == 5 );
        assert( nc == 8 );
        
        // Generate a valid row cut from thin air
        int c;
        {
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *el = new double[nc];
          for (c=0;c<nc;c++) el[c]=((double)c)*((double)c);
          
          OsiRowCut rc;
          rc.setRow(nc,inx,el);
          rc.setLb(-100.);
          rc.setUb(100.);
          rc.setEffectiveness(22);
          
          cuts.insert(rc);
          delete[]el;
          delete[]inx;
        }
        
        // Generate valid col cut from thin air
        {
          const double * soplexColLB = im.getColLower();
          const double * soplexColUB = im.getColUpper();
          int *inx = new int[nc];
          for (c=0;c<nc;c++) inx[c]=c;
          double *lb = new double[nc];
          double *ub = new double[nc];
          for (c=0;c<nc;c++) lb[c]=soplexColLB[c]+0.001;
          for (c=0;c<nc;c++) ub[c]=soplexColUB[c]-0.001;
          
          OsiColCut cc;
          cc.setLbs(nc,inx,lb);
          cc.setUbs(nc,inx,ub);
          
          cuts.insert(cc);
          delete [] ub;
          delete [] lb;
          delete [] inx;
        }
        
        {
          // Generate a row and column cut which have are ineffective
          OsiRowCut * rcP= new OsiRowCut;
          rcP->setEffectiveness(-1.);
          cuts.insert(rcP);
          assert(rcP==NULL);
          
          OsiColCut * ccP= new OsiColCut;
          ccP->setEffectiveness(-12.);
          cuts.insert(ccP);
          assert(ccP==NULL);
        }
        {
          //Generate inconsistent Row cut
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          rc.setLb(3.);
          rc.setUb(4.);
          assert(!rc.consistent());
          cuts.insert(rc);
        }
        {
          //Generate inconsistent col cut
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={-10};
          double el[ne]={2.5};
          cc.setUbs(ne,inx,el);
          assert(!cc.consistent());
          cuts.insert(cc);
        }
        {
          // Generate row cut which is inconsistent for model m
          OsiRowCut rc;
          const int ne=1;
          int inx[ne]={10};
          double el[ne]={2.5};
          rc.setRow(ne,inx,el);
          assert(rc.consistent());
          assert(!rc.consistent(im));
          cuts.insert(rc);
        }
        {
          // Generate col cut which is inconsistent for model m
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={30};
          double el[ne]={2.0};
          cc.setLbs(ne,inx,el);
          assert(cc.consistent());
          assert(!cc.consistent(im));
          cuts.insert(cc);
        }
        {
          // Generate col cut which is infeasible
          OsiColCut cc;
          const int ne=1;
          int inx[ne]={0};
          double el[ne]={2.0};
          cc.setUbs(ne,inx,el);
          cc.setEffectiveness(1000.);
          assert(cc.consistent());
          assert(cc.consistent(im));
          assert(cc.infeasible(im));
          cuts.insert(cc);
        }
      }
      assert(cuts.sizeRowCuts()==4);
      assert(cuts.sizeColCuts()==5);
      
      OsiSolverInterface::ApplyCutsReturnCode rc = im.applyCuts(cuts);
      assert( rc.getNumIneffective() == 2 );
      assert( rc.getNumApplied() == 2 );
      assert( rc.getNumInfeasible() == 1 );
      assert( rc.getNumInconsistentWrtIntegerModel() == 2 );
      assert( rc.getNumInconsistent() == 2 );
      assert( cuts.sizeCuts() == rc.getNumIneffective() +
        rc.getNumApplied() +
        rc.getNumInfeasible() +
        rc.getNumInconsistentWrtIntegerModel() +
        rc.getNumInconsistent() );
    }
    {    
      OsiSpxSolverInterface soplexSi(m);
      int nc = soplexSi.getNumCols();
      int nr = soplexSi.getNumRows();
      const double * cl = soplexSi.getColLower();
      const double * cu = soplexSi.getColUpper();
      const double * rl = soplexSi.getRowLower();
      const double * ru = soplexSi.getRowUpper();

      assert( nc == 8 );
      assert( nr == 5 );
      assert( eq(cl[0],2.5) );
      assert( eq(cl[1],0.0) );
      assert( eq(cu[1],4.1) );
      assert( eq(cu[2],1.0) );

      assert( eq(rl[0],2.5) );
      assert( eq(rl[4],3.0) );
      assert( eq(ru[1],2.1) );
      assert( eq(ru[4],15.0) );
      
      double newCs[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
      soplexSi.setColSolution(newCs);
      const double * cs = soplexSi.getColSolution();
      assert( eq(cs[0],1.0) );
      assert( eq(cs[7],8.0) );
      {
        OsiSpxSolverInterface solnSi(soplexSi);
        const double * cs = solnSi.getColSolution();
        assert( eq(cs[0],1.0) );
        assert( eq(cs[7],8.0) );
      }

      assert( !eq(cl[3],1.2345) );
      soplexSi.setColLower( 3, 1.2345 );
      assert( eq(soplexSi.getColLower()[3],1.2345) );
      
      assert( !eq(cu[4],10.2345) );
      soplexSi.setColUpper( 4, 10.2345 );
      assert( eq(soplexSi.getColUpper()[4],10.2345) );

      assert( eq(soplexSi.getObjValue(),0.0) );

      assert( eq( soplexSi.getObjCoefficients()[0],  1.0) );
      assert( eq( soplexSi.getObjCoefficients()[1],  0.0) );
      assert( eq( soplexSi.getObjCoefficients()[2],  0.0) );
      assert( eq( soplexSi.getObjCoefficients()[3],  0.0) );
      assert( eq( soplexSi.getObjCoefficients()[4],  2.0) );
      assert( eq( soplexSi.getObjCoefficients()[5],  0.0) );
      assert( eq( soplexSi.getObjCoefficients()[6],  0.0) );
      assert( eq( soplexSi.getObjCoefficients()[7], -1.0) );
    }
    
    // Test getMatrixByRow method
    { 
      const OsiSpxSolverInterface si(m);
      const OsiPackedMatrix * smP = si.getMatrixByRow();
      //const OsiPackedMatrix * osmP = dynamic_cast(const OsiSpxPackedMatrix*)(smP);
      //assert( osmP!=NULL );
      
      OsiRelFltEq eq;
      const double * ev = smP->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const int * mi = smP->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = smP->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    
      
      assert( smP->getMajorDim() == 5 ); 
      assert( smP->getNumElements() == 14 );
      
    }
    //--------------
    // Test rowsense, rhs, rowrange, getMatrixByRow
    {
      OsiSpxSolverInterface lhs;
      {     
        OsiSpxSolverInterface siC1(m);     
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );
	assert( siC1.colsol_!=NULL );
	assert( siC1.rowsol_!=NULL );
        assert( siC1.matrixByRow_==NULL );

        const char   * siC1rs  = siC1.getRowSense();
        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        
        const double * siC1rhs = siC1.getRightHandSide();
        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        
        const double * siC1rr  = siC1.getRowRange();
        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        
        const OsiPackedMatrix * siC1mbr = siC1.getMatrixByRow();
        assert( siC1mbr != NULL );
        
        const double * ev = siC1mbr->getElements();
        assert( eq(ev[0],   3.0) );
        assert( eq(ev[1],   1.0) );
        assert( eq(ev[2],  -2.0) );
        assert( eq(ev[3],  -1.0) );
        assert( eq(ev[4],  -1.0) );
        assert( eq(ev[5],   2.0) );
        assert( eq(ev[6],   1.1) );
        assert( eq(ev[7],   1.0) );
        assert( eq(ev[8],   1.0) );
        assert( eq(ev[9],   2.8) );
        assert( eq(ev[10], -1.2) );
        assert( eq(ev[11],  5.6) );
        assert( eq(ev[12],  1.0) );
        assert( eq(ev[13],  1.9) );
        
        const int * mi = siC1mbr->getVectorStarts();
        assert( mi[0]==0 );
        assert( mi[1]==5 );
        assert( mi[2]==7 );
        assert( mi[3]==9 );
        assert( mi[4]==11 );
        assert( mi[5]==14 );
        
        const int * ei = siC1mbr->getIndices();
        assert( ei[0]  ==  0 );
        assert( ei[1]  ==  1 );
        assert( ei[2]  ==  3 );
        assert( ei[3]  ==  4 );
        assert( ei[4]  ==  7 );
        assert( ei[5]  ==  1 );
        assert( ei[6]  ==  2 );
        assert( ei[7]  ==  2 );
        assert( ei[8]  ==  5 );
        assert( ei[9]  ==  3 );
        assert( ei[10] ==  6 );
        assert( ei[11] ==  0 );
        assert( ei[12] ==  4 );
        assert( ei[13] ==  7 );    
        
        assert( siC1mbr->getMajorDim() == 5 ); 
        assert( siC1mbr->getNumElements() == 14 );
        

        assert( siC1rs  == siC1.getRowSense() );
        assert( siC1rhs == siC1.getRightHandSide() );
        assert( siC1rr  == siC1.getRowRange() );

        // Change SOPLEX Model by adding free row
        OsiRowCut rc;
        rc.setLb(-DBL_MAX);
        rc.setUb( DBL_MAX);
        OsiCuts cuts;
        cuts.insert(rc);
        siC1.applyCuts(cuts);
             
        // Since model was changed, test that cached
        // data is now freed.
        assert( siC1.rowrange_==NULL );
        assert( siC1.rowsense_==NULL );
        assert( siC1.rhs_==NULL );
	assert( siC1.colsol_==NULL );
	assert( siC1.rowsol_==NULL );
	assert( siC1.matrixByRow_==NULL );
	assert( siC1.matrixByCol_==NULL );
        
        siC1rs  = siC1.getRowSense();
        siC1rhs = siC1.getRightHandSide();
        siC1rr  = siC1.getRowRange();

        assert( siC1rs[0]=='G' );
        assert( siC1rs[1]=='L' );
        assert( siC1rs[2]=='E' );
        assert( siC1rs[3]=='R' );
        assert( siC1rs[4]=='R' );
        assert( siC1rs[5]=='N' );

        assert( eq(siC1rhs[0],2.5) );
        assert( eq(siC1rhs[1],2.1) );
        assert( eq(siC1rhs[2],4.0) );
        assert( eq(siC1rhs[3],5.0) );
        assert( eq(siC1rhs[4],15.) ); 
        assert( eq(siC1rhs[5],0.0) ); 

        assert( eq(siC1rr[0],0.0) );
        assert( eq(siC1rr[1],0.0) );
        assert( eq(siC1rr[2],0.0) );
        assert( eq(siC1rr[3],5.0-1.8) );
        assert( eq(siC1rr[4],15.0-3.0) );
        assert( eq(siC1rr[5],0.0) );
    
        lhs=siC1;
      }
      // Test that lhs has correct values even though siC1 has gone out of scope    
      assert( lhs.rowrange_==NULL );
      assert( lhs.rowsense_==NULL );
      assert( lhs.rhs_==NULL ); 
      assert( lhs.colsol_!=NULL );
      assert( lhs.rowsol_!=NULL );
      assert( lhs.matrixByRow_==NULL ); 
      assert( lhs.matrixByCol_==NULL ); 
      
      const char * lhsrs  = lhs.getRowSense();
      assert( lhsrs[0]=='G' );
      assert( lhsrs[1]=='L' );
      assert( lhsrs[2]=='E' );
      assert( lhsrs[3]=='R' );
      assert( lhsrs[4]=='R' );
      assert( lhsrs[5]=='N' );
      
      const double * lhsrhs = lhs.getRightHandSide();
      assert( eq(lhsrhs[0],2.5) );
      assert( eq(lhsrhs[1],2.1) );
      assert( eq(lhsrhs[2],4.0) );
      assert( eq(lhsrhs[3],5.0) );
      assert( eq(lhsrhs[4],15.) ); 
      assert( eq(lhsrhs[5],0.0) ); 
      
      const double *lhsrr  = lhs.getRowRange();
      assert( eq(lhsrr[0],0.0) );
      assert( eq(lhsrr[1],0.0) );
      assert( eq(lhsrr[2],0.0) );
      assert( eq(lhsrr[3],5.0-1.8) );
      assert( eq(lhsrr[4],15.0-3.0) );
      assert( eq(lhsrr[5],0.0) );      
      
      const OsiPackedMatrix * lhsmbr = lhs.getMatrixByRow();
      assert( lhsmbr != NULL );       
      const double * ev = lhsmbr->getElements();
      assert( eq(ev[0],   3.0) );
      assert( eq(ev[1],   1.0) );
      assert( eq(ev[2],  -2.0) );
      assert( eq(ev[3],  -1.0) );
      assert( eq(ev[4],  -1.0) );
      assert( eq(ev[5],   2.0) );
      assert( eq(ev[6],   1.1) );
      assert( eq(ev[7],   1.0) );
      assert( eq(ev[8],   1.0) );
      assert( eq(ev[9],   2.8) );
      assert( eq(ev[10], -1.2) );
      assert( eq(ev[11],  5.6) );
      assert( eq(ev[12],  1.0) );
      assert( eq(ev[13],  1.9) );
      
      const int * mi = lhsmbr->getVectorStarts();
      assert( mi[0]==0 );
      assert( mi[1]==5 );
      assert( mi[2]==7 );
      assert( mi[3]==9 );
      assert( mi[4]==11 );
      assert( mi[5]==14 );
      
      const int * ei = lhsmbr->getIndices();
      assert( ei[0]  ==  0 );
      assert( ei[1]  ==  1 );
      assert( ei[2]  ==  3 );
      assert( ei[3]  ==  4 );
      assert( ei[4]  ==  7 );
      assert( ei[5]  ==  1 );
      assert( ei[6]  ==  2 );
      assert( ei[7]  ==  2 );
      assert( ei[8]  ==  5 );
      assert( ei[9]  ==  3 );
      assert( ei[10] ==  6 );
      assert( ei[11] ==  0 );
      assert( ei[12] ==  4 );
      assert( ei[13] ==  7 );    
      
      int md = lhsmbr->getMajorDim();
      assert(  md == 6 ); 
      assert( lhsmbr->getNumElements() == 14 );
    }
    
    //--------------    
  }

    
  // Do common solverInterface testing by calling the
  // base class testing method.
  {
    OsiSpxSolverInterface m;
    OsiSolverInterfaceCommonUnitTest(&m, mpsDir);
  }
}
#endif
