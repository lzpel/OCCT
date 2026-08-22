// Created on: 1998-07-02
// Created by: Stephanie HUMEAU
// Copyright (c) 1998-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <GeomFill_GuideTrihedronPlan.hxx>

#include <Adaptor3d_Curve.hxx>
#include <ElCLib.hxx>
#include <Geom_Plane.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomFill_Frenet.hxx>
#include <GeomFill_PlanFunc.hxx>
#include <GeomFill_TrihedronLaw.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec.hxx>
#include <IntCurveSurface_HInter.hxx>
#include <IntCurveSurface_IntersectionPoint.hxx>
#include <NCollection_Sequence.hxx>
#include <math_FunctionRoot.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>
#include <Precision.hxx>
#include <Standard_Type.hxx>

IMPLEMENT_STANDARD_RTTIEXT(GeomFill_GuideTrihedronPlan, GeomFill_TrihedronWithGuide)

// cadrum fix A: the guide reaching this law is often a BRepAdaptor_CompCurve that has
// lost the periodicity of its underlying closed curve (a closed wire is split into
// unit-parameter edge spans). Detect closure geometrically so seam handling works for
// both truly periodic adaptors and closed composite ones.
static bool cadrumGuideClosed(const occ::handle<Adaptor3d_Curve>& theGuide, double& thePeriod)
{
  thePeriod = theGuide->LastParameter() - theGuide->FirstParameter();
  if (theGuide->IsPeriodic())
    return true;
  const gp_Pnt aP1 = theGuide->Value(theGuide->FirstParameter());
  const gp_Pnt aP2 = theGuide->Value(theGuide->LastParameter());
  return aP1.Distance(aP2) < 1.e-4 * thePeriod;
}

// #include <gp_Trsf2d.hxx>
// #include <Bnd_Box2d.hxx>

#ifdef OCCT_DEBUG
static void TracePlan(const occ::handle<Geom_Surface>& /*Plan*/)
{
  std::cout << "Pas d'intersection Guide/Plan" << std::endl;
}
#endif

//=================================================================================================

static void InGoodPeriod(const double Prec, const double Period, double& Current)
{
  double Diff = Current - Prec;
  int    nb   = (int)std::trunc(Diff / Period);
  Current -= nb * Period;
  Diff = Current - Prec;
  if (Diff > Period / 2)
  {
    Current -= Period;
  }
  else if (Diff < -Period / 2)
  {
    Current += Period;
  }
}

//=================================================================================================

GeomFill_GuideTrihedronPlan::GeomFill_GuideTrihedronPlan(
  const occ::handle<Adaptor3d_Curve>& theGuide)
    : X(1, 1),
      XTol(1, 1),
      Inf(1, 1),
      Sup(1, 1),
      myStatus(GeomFill_PipeOk)
{
  myCurve.Nullify();
  myGuide = theGuide; // guide
  myTrimG = theGuide;
  myNbPts = 20; // nb points pour calculs
  Pole =
    new (NCollection_HArray2<gp_Pnt2d>)(1, 1, 1, myNbPts); // tab pr stocker Pprime (pt sur guide)
  frenet = new (GeomFill_Frenet)();
  XTol.Init(1.e-6);
  XTol(1) = myGuide->Resolution(1.e-6);
}

//=======================================================================
// function : Init
// purpose  : calcule myNbPts points sur la courbe guide (<=> normale)
//=======================================================================
void GeomFill_GuideTrihedronPlan::Init()
{
  myStatus = GeomFill_PipeOk;
  double f = myCurve->FirstParameter();
  double l = myCurve->LastParameter();
  double DeltaG;

  frenet->SetCurve(myCurve);
  DeltaG = (myGuide->LastParameter() - myGuide->FirstParameter()) / 2;

  Inf(1) = myGuide->FirstParameter() - DeltaG;
  Sup(1) = myGuide->LastParameter() + DeltaG;

  if (!myGuide->IsPeriodic())
  {
    myTrimG = myGuide->Trim(myGuide->FirstParameter() - DeltaG / 100,
                            myGuide->LastParameter() + DeltaG / 100,
                            DeltaG * 1.e-7);
  }
  else
  {
    myTrimG = myGuide;
  }
  DeltaG /= 3;

  double     aGPeriod  = 0.;
  const bool isGClosed = cadrumGuideClosed(myGuide, aGPeriod);

  // cadrum fix C: one plane/guide solve at parameter theT. The branch is chosen by
  // parameter continuity with thePrevW (spatially for the very first sample); a
  // closed guide is cut twice by every normal plane and spatial selection can hop
  // to the antipodal branch near the seam.
  auto aSolveAt = [&](const double theT,
                      const bool   theHasPrev,
                      const double thePrevW,
                      double&      theW) -> bool {
    gp_Pnt aP;
    gp_Vec aTg, aNm, aBn;
    myCurve->D0(theT, aP);
    frenet->D0(theT, aTg, aNm, aBn);
    occ::handle<Geom_Plane>          aPlan = new (Geom_Plane)(aP, aTg);
    occ::handle<GeomAdaptor_Surface> aPl   = new (GeomAdaptor_Surface)(aPlan);
    IntCurveSurface_HInter           anInt;
    anInt.Perform(myTrimG, aPl);
    if (anInt.NbPoints() == 0)
      return false;
    bool   aFirst = true;
    double aBest  = 0.;
    if (!theHasPrev)
    {
      double aDmin = 0.;
      for (int jj = 1; jj <= anInt.NbPoints(); jj++)
      {
        const double aD = aP.Distance(anInt.Point(jj).Pnt());
        if (aFirst || aD < aDmin)
        {
          aDmin  = aD;
          aBest  = anInt.Point(jj).W();
          aFirst = false;
        }
      }
    }
    else
    {
      for (int jj = 1; jj <= anInt.NbPoints(); jj++)
      {
        double aW = anInt.Point(jj).W();
        if (isGClosed)
          InGoodPeriod(thePrevW, aGPeriod, aW);
        if (aFirst || std::abs(aW - thePrevW) < std::abs(aBest - thePrevW))
        {
          aBest  = aW;
          aFirst = false;
        }
      }
    }
    theW = aBest;
    return true;
  };

  // Seed pass at myNbPts uniform samples.
  NCollection_Sequence<gp_Pnt2d> aSeq;
  bool   aHasPrev = false;
  double aPrevW   = 0.;
  for (int ii = 1; ii <= myNbPts; ii++)
  {
    double t = double(myNbPts - ii) * f + double(ii - 1) * l;
    t /= (myNbPts - 1);
    double aW = 0.;
    if (aSolveAt(t, aHasPrev, aPrevW, aW))
    {
      aSeq.Append(gp_Pnt2d(t, aW));
      aHasPrev = true;
      aPrevW   = aW;
    }
    else if (!aHasPrev)
    {
      myStatus = GeomFill_PlaneNotIntersectGuide;
      aSeq.Append(gp_Pnt2d(t, myGuide->FirstParameter()));
    }
    else
    {
      aSeq.Append(gp_Pnt2d(t, aPrevW)); // keep the branch value; refined below
    }
  }

  // cadrum fix C: refine adaptively until neighbouring samples differ by less than
  // DeltaG in guide parameter, so evaluation-time seeds interpolate inside one
  // branch. Bounded to keep degenerate inputs from exploding the table.
  const int aMaxPts   = 20 * myNbPts;
  bool      aInserted = true;
  while (aInserted && aSeq.Length() < aMaxPts)
  {
    aInserted = false;
    for (int ii = 1; ii < aSeq.Length() && aSeq.Length() < aMaxPts; ii++)
    {
      const gp_Pnt2d a1 = aSeq(ii);
      const gp_Pnt2d a2 = aSeq(ii + 1);
      if (std::abs(a2.Y() - a1.Y()) <= DeltaG || (a2.X() - a1.X()) < 1.e-9 * (l - f))
        continue;
      const double aTm = 0.5 * (a1.X() + a2.X());
      double       aWm = 0.;
      if (!aSolveAt(aTm, true, a1.Y(), aWm))
        aWm = 0.5 * (a1.Y() + a2.Y());
      aSeq.InsertAfter(ii, gp_Pnt2d(aTm, aWm));
      aInserted = true;
      ii++;
    }
  }

  Pole = new (NCollection_HArray2<gp_Pnt2d>)(1, 1, 1, aSeq.Length());
  for (int ii = 1; ii <= aSeq.Length(); ii++)
    Pole->SetValue(1, ii, aSeq(ii));
}

//=================================================================================================

bool GeomFill_GuideTrihedronPlan::SetCurve(const occ::handle<Adaptor3d_Curve>& C)
{
  myCurve = C;
  if (!myCurve.IsNull())
  {
    Init();
  }
  return true;
}

//=================================================================================================

occ::handle<Adaptor3d_Curve> GeomFill_GuideTrihedronPlan::Guide() const
{
  return myGuide;
}

//=================================================================================================

bool GeomFill_GuideTrihedronPlan::D0(const double Param,
                                     gp_Vec&      Tangent,
                                     gp_Vec&      Normal,
                                     gp_Vec&      BiNormal)
{
  gp_Pnt P, Pprime;
  //  gp_Vec To;

  myCurve->D0(Param, P);

  frenet->D0(Param, Tangent, Normal, BiNormal);

  // initialisation de la recherche
  InitX(Param);

  int Iter = 50;

  // fonction dont il faut trouver la racine : G(W)-Pl(U,V)=0
  GeomFill_PlanFunc E(P, Tangent, myGuide);

  // resolution
  // cadrum fix C: polish the interpolated seed inside its own table bracket. The
  // member bounds spanned twice the guide range and let the solver capture the
  // antipodal branch of a closed guide. Should the local search fail, retry from
  // the other periodic image of the seed; as a last resort keep the seed itself -
  // the refined table already makes it accurate to within DeltaG.
  double aGPeriod = 0.;
  const bool   isGClosed = cadrumGuideClosed(myGuide, aGPeriod);
  const double aPad      = 0.05 * aGPeriod;
  math_FunctionRoot Result(E, X(1), XTol(1), std::min(myWLo, X(1)) - aPad, std::max(myWHi, X(1)) + aPad, Iter);
  bool   aRootDone = Result.IsDone();
  double aRootVal  = aRootDone ? Result.Root() : 0.;
  if (!aRootDone && isGClosed)
  {
    const double aXAlt = (X(1) - myGuide->FirstParameter() > aGPeriod / 2.) ? X(1) - aGPeriod : X(1) + aGPeriod;
    math_FunctionRoot ResultAlt(E, aXAlt, XTol(1), aXAlt - aPad, aXAlt + aPad, Iter);
    if (ResultAlt.IsDone())
    {
      aRootDone = true;
      aRootVal  = ResultAlt.Root();
    }
  }
  if (!aRootDone)
  {
    aRootDone = true;
    aRootVal  = X(1);
  }

  if (aRootDone)
  {
    double Res = aRootVal;
    //      R = Result.Root();    // solution

    Pprime = myTrimG->Value(Res); // pt sur courbe guide
    gp_Vec n(P, Pprime);          // vecteur definissant la normale du triedre

    Normal   = n.Normalized();
    BiNormal = Tangent.Crossed(Normal);
    BiNormal.Normalize();
  }
  else
  { // Erreur...
#ifdef OCCT_DEBUG
    std::cout << "D0 :";
    // plan ortho a la trajectoire pour determiner Pprime
    occ::handle<Geom_Plane> Plan = new (Geom_Plane)(P, Tangent);
    TracePlan(Plan);
#endif
    myStatus = GeomFill_PlaneNotIntersectGuide;
    return false;
  }

  return true;
}

//=======================================================================
// function : D1
// purpose  : calculation of trihedron and first derivative
//=======================================================================
bool GeomFill_GuideTrihedronPlan::D1(const double Param,
                                     gp_Vec&      Tangent,
                                     gp_Vec&      DTangent,
                                     gp_Vec&      Normal,
                                     gp_Vec&      DNormal,
                                     gp_Vec&      BiNormal,
                                     gp_Vec&      DBiNormal)
{
  //  return false;
  gp_Pnt P, PG;
  gp_Vec To, TG;

  // triedre de frenet sur la trajectoire
  myCurve->D1(Param, P, To);
  frenet->D1(Param, Tangent, DTangent, Normal, DNormal, BiNormal, DBiNormal);

  // tolerance sur E
  int Iter = 50;

  // fonction dont il faut trouver la racine : G(W)-Pl(U,V)=0
  InitX(Param);
  GeomFill_PlanFunc E(P, Tangent, myGuide);

  // resolution
  // cadrum fix C: polish the interpolated seed inside its own table bracket. The
  // member bounds spanned twice the guide range and let the solver capture the
  // antipodal branch of a closed guide. Should the local search fail, retry from
  // the other periodic image of the seed; as a last resort keep the seed itself -
  // the refined table already makes it accurate to within DeltaG.
  double aGPeriod = 0.;
  const bool   isGClosed = cadrumGuideClosed(myGuide, aGPeriod);
  const double aPad      = 0.05 * aGPeriod;
  math_FunctionRoot Result(E, X(1), XTol(1), std::min(myWLo, X(1)) - aPad, std::max(myWHi, X(1)) + aPad, Iter);
  bool   aRootDone = Result.IsDone();
  double aRootVal  = aRootDone ? Result.Root() : 0.;
  if (!aRootDone && isGClosed)
  {
    const double aXAlt = (X(1) - myGuide->FirstParameter() > aGPeriod / 2.) ? X(1) - aGPeriod : X(1) + aGPeriod;
    math_FunctionRoot ResultAlt(E, aXAlt, XTol(1), aXAlt - aPad, aXAlt + aPad, Iter);
    if (ResultAlt.IsDone())
    {
      aRootDone = true;
      aRootVal  = ResultAlt.Root();
    }
  }
  if (!aRootDone)
  {
    aRootDone = true;
    aRootVal  = X(1);
  }

  if (aRootDone)
  {
    double Res = aRootVal;
    //      R = Result.Root();    // solution
    myTrimG->D1(Res, PG, TG);
    gp_Vec n(P, PG), dn; // vecteur definissant la normale du triedre
    double Norm = n.Magnitude();
    if (Norm < 1.e-12)
    {
      Norm = 1.0;
    }
    n /= Norm;

    Normal   = n;
    BiNormal = Tangent.Crossed(Normal);

    // derivee premiere du triedre
    double dedx, dedt, dtg_dt;
    E.Derivative(Res, dedx);
    E.DEDT(Res, To, DTangent, dedt);
    dtg_dt = -dedt / dedx;

    /*      double h=1.e-7, e, etg, etc;
          E.Value(Res, e);
          E.Value(Res+h, etg);
          if ( std::abs( (etg-e)/h - dedx) > 1.e-4) {
        std::cout << "err :" <<  (etg-e)/h - dedx << std::endl;
          }
          gp_Pnt pdbg;
          gp_Vec td, nb, bnb;
          myCurve->D0(Param+h, pdbg);
          frenet->D0(Param+h,td, nb, bnb);

          GeomFill_PlanFunc Edeb(pdbg, td, myGuide);
          Edeb.Value(Res, etc);
          if ( std::abs( (etc-e)/h - dedt) > 1.e-4) {
        std::cout << "err :" <<  (etc-e)/h - dedt << std::endl;
          } */

    dn.SetLinearForm(dtg_dt, TG, -1, To);

    DNormal.SetLinearForm(-(n * dn), n, dn);
    DNormal /= Norm;
    DBiNormal.SetLinearForm(Tangent.Crossed(DNormal), DTangent.Crossed(Normal));
  }
  else
  { // Erreur...
#ifdef OCCT_DEBUG
    std::cout << "D1 :";
    // plan ortho a la trajectoire
    occ::handle<Geom_Plane> Plan = new (Geom_Plane)(P, Tangent);
    TracePlan(Plan);
#endif
    myStatus = GeomFill_PlaneNotIntersectGuide;
    return false;
  }

  return true;
}

//=======================================================================
// function : D2
// purpose  : calculation of trihedron and derivatives
//=======================================================================
bool GeomFill_GuideTrihedronPlan::D2(const double Param,
                                     gp_Vec&      Tangent,
                                     gp_Vec&      DTangent,
                                     gp_Vec&      D2Tangent,
                                     gp_Vec&      Normal,
                                     gp_Vec&      DNormal,
                                     gp_Vec&      D2Normal,
                                     gp_Vec&      BiNormal,
                                     gp_Vec&      DBiNormal,
                                     gp_Vec&      D2BiNormal)
{
  //  gp_Pnt P, PG;
  gp_Pnt P;
  //  gp_Vec To,DTo,TG,DTG;
  gp_Vec To, DTo;

  myCurve->D2(Param, P, To, DTo);

  // triedre de Frenet sur la trajectoire
  frenet->D2(Param,
             Tangent,
             DTangent,
             D2Tangent,
             Normal,
             DNormal,
             D2Normal,
             BiNormal,
             DBiNormal,
             D2BiNormal);

  return false;
}

//=================================================================================================

occ::handle<GeomFill_TrihedronLaw> GeomFill_GuideTrihedronPlan::Copy() const
{
  occ::handle<GeomFill_GuideTrihedronPlan> copy = new (GeomFill_GuideTrihedronPlan)(myGuide);
  copy->SetCurve(myCurve);
  return copy;
}

//=================================================================================================

GeomFill_PipeError GeomFill_GuideTrihedronPlan::ErrorStatus() const
{
  return myStatus;
}

//=======================================================================
// function : NbIntervals
// purpose  : Version provisoire : Il faut tenir compte du guide
//=======================================================================
int GeomFill_GuideTrihedronPlan::NbIntervals(const GeomAbs_Shape S) const
{
  int           Nb;
  GeomAbs_Shape tmpS;
  switch (S)
  {
    case GeomAbs_C0:
      tmpS = GeomAbs_C1;
      break;
    case GeomAbs_C1:
      tmpS = GeomAbs_C2;
      break;
    case GeomAbs_C2:
      tmpS = GeomAbs_C3;
      break;
    default:
      tmpS = GeomAbs_CN;
  }

  Nb = myCurve->NbIntervals(tmpS);
  return Nb;
}

//=================================================================================================

void GeomFill_GuideTrihedronPlan::Intervals(NCollection_Array1<double>& TT,
                                            const GeomAbs_Shape         S) const
{
  GeomAbs_Shape tmpS;
  switch (S)
  {
    case GeomAbs_C0:
      tmpS = GeomAbs_C1;
      break;
    case GeomAbs_C1:
      tmpS = GeomAbs_C2;
      break;
    case GeomAbs_C2:
      tmpS = GeomAbs_C3;
      break;
    default:
      tmpS = GeomAbs_CN;
  }
  myCurve->Intervals(TT, tmpS);
}

//=================================================================================================

void GeomFill_GuideTrihedronPlan::SetInterval(const double First, const double Last)
{
  myTrimmed = myCurve->Trim(First, Last, Precision::Confusion());
}

//=================================================================================================

void GeomFill_GuideTrihedronPlan::GetAverageLaw(gp_Vec& ATangent,
                                                gp_Vec& ANormal,
                                                gp_Vec& ABiNormal)
{
  int    ii;
  double t, Delta = (myCurve->LastParameter() - myCurve->FirstParameter()) / 20.001;

  ATangent.SetCoord(0., 0., 0.);
  ANormal.SetCoord(0., 0., 0.);
  ABiNormal.SetCoord(0., 0., 0.);
  gp_Vec T, N, B;

  for (ii = 1; ii <= 20; ii++)
  {
    t = myCurve->FirstParameter() + (ii - 1) * Delta;
    D0(t, T, N, B);
    ATangent += T;
    ANormal += N;
    ABiNormal += B;
  }
  ATangent /= 20;
  ANormal /= 20;
  ABiNormal /= 20;
}

//=================================================================================================

bool GeomFill_GuideTrihedronPlan::IsConstant() const
{
  if ((myCurve->GetType() == GeomAbs_Line) && (myGuide->GetType() == GeomAbs_Line))
  {
    double Angle;
    Angle = myCurve->Line().Angle(myGuide->Line());
    if ((Angle < 1.e-12) || ((2 * M_PI - Angle) < 1.e-12))
    {
      return true;
    }
  }

  return false;
}

//=================================================================================================

bool GeomFill_GuideTrihedronPlan::IsOnlyBy3dCurve() const
{
  return false;
}

//=======================================================================
// function : Origine
// purpose  : Nothing!!
//=======================================================================
void GeomFill_GuideTrihedronPlan::Origine(const double, const double) {}

//==================================================================
// Function : InitX
// Purpose : recherche par interpolation d'une valeur initiale
//==================================================================
void GeomFill_GuideTrihedronPlan::InitX(const double Param)
{

  int    Ideb = 1, Ifin = Pole->RowLength(), Idemi;
  double Valeur, t1, t2;

  Valeur = Pole->Value(1, Ideb).X();
  if (Param == Valeur)
  {
    Ifin = Ideb + 1;
  }

  Valeur = Pole->Value(1, Ifin).X();
  if (Param == Valeur)
  {
    Ideb = Ifin - 1;
  }

  while (Ideb + 1 != Ifin)
  {
    Idemi  = (Ideb + Ifin) / 2;
    Valeur = Pole->Value(1, Idemi).X();
    if (Valeur < Param)
    {
      Ideb = Idemi;
    }
    else
    {
      if (Valeur > Param)
      {
        Ifin = Idemi;
      }
      else
      {
        Ideb = Idemi;
        Ifin = Ideb + 1;
      }
    }
  }

  t1          = Pole->Value(1, Ideb).X();
  t2          = Pole->Value(1, Ifin).X();
  double diff = t2 - t1;
  if (diff > 1.e-7)
  {
    double b = (Param - t1) / diff, a = (t2 - Param) / diff;

    X(1) = Pole->Value(1, Ideb).Coord(2) * a + Pole->Value(1, Ifin).Coord(2) * b; // param guide
  }
  else
  {
    X(1) = (Pole->Value(1, Ideb).Coord(2) + Pole->Value(1, Ifin).Coord(2)) / 2;
  }
  // cadrum fix C: remember the local bracket the seed was interpolated from, so the
  // evaluation-time root search can be confined to it.
  myWLo = std::min(Pole->Value(1, Ideb).Coord(2), Pole->Value(1, Ifin).Coord(2));
  myWHi = std::max(Pole->Value(1, Ideb).Coord(2), Pole->Value(1, Ifin).Coord(2));

  // Fold a seed interpolated past the seam of a closed guide back into the evaluable
  // range (composite adaptors of a closed wire report IsPeriodic()=false while still
  // wrapping at LastParameter); shift the bracket together with the seed.
  double aGPeriod = 0.;
  if (cadrumGuideClosed(myGuide, aGPeriod))
  {
    const double aFolded = ElCLib::InPeriod(X(1), myGuide->FirstParameter(), myGuide->LastParameter());
    myWLo += aFolded - X(1);
    myWHi += aFolded - X(1);
    X(1) = aFolded;
  }
}
