// Copyright (c) 2026 OPEN CASCADE SAS
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

// Fillets on two edges that share no vertex, where the stripes hold a different
// number of surface data. ChFi3d_Builder::PerformMoreSurfdata indexed the first
// stripe with an index computed for the stripe at the vertex.

#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GProp_GProps.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

#include <gtest/gtest.h>

namespace
{
//! Returns the edge of the shape whose middle is at the given point.
TopoDS_Edge findEdgeAt(const TopoDS_Shape& theShape, const gp_Pnt& theMiddle)
{
  for (TopExp_Explorer anExp(theShape, TopAbs_EDGE); anExp.More(); anExp.Next())
  {
    const TopoDS_Edge& anEdge = TopoDS::Edge(anExp.Current());
    BRepAdaptor_Curve  aCurve(anEdge);
    const gp_Pnt aMiddle = aCurve.Value(0.5 * (aCurve.FirstParameter() + aCurve.LastParameter()));
    if (aMiddle.Distance(theMiddle) < 1.e-6)
    {
      return anEdge;
    }
  }
  return TopoDS_Edge();
}

//! Volume of the shape.
double volumeOf(const TopoDS_Shape& theShape)
{
  GProp_GProps aProps;
  BRepGProp::VolumeProperties(theShape, aProps);
  return aProps.Mass();
}

//! A 22.5 square whose left side bows out into an arc, extruded 5 along Z.
TopoDS_Shape makeArcPrism()
{
  const gp_Pnt aP0(0.0, 0.0, 0.0), aP1(22.5, 0.0, 0.0), aP2(22.5, 22.5, 0.0), aP3(0.0, 22.5, 0.0);
  BRepBuilderAPI_MakeWire aWire;
  aWire.Add(BRepBuilderAPI_MakeEdge(aP0, aP1).Edge());
  aWire.Add(BRepBuilderAPI_MakeEdge(aP1, aP2).Edge());
  aWire.Add(BRepBuilderAPI_MakeEdge(aP2, aP3).Edge());
  aWire.Add(
    BRepBuilderAPI_MakeEdge(GC_MakeArcOfCircle(aP3, gp_Pnt(-7.5, 11.25, 0.0), aP0).Value()).Edge());
  const TopoDS_Face aFace = BRepBuilderAPI_MakeFace(aWire.Wire()).Face();
  return BRepPrimAPI_MakePrism(aFace, gp_Vec(0.0, 0.0, 5.0)).Shape();
}

//! Volume removed by a fillet of the given radius on the given edges, which must succeed.
double removedBy(const TopoDS_Shape& theShape, double theRadius, const TopoDS_Edge& theEdge)
{
  BRepFilletAPI_MakeFillet aFillet(theShape);
  aFillet.Add(theRadius, theEdge);
  aFillet.Build();
  EXPECT_TRUE(aFillet.IsDone());
  return volumeOf(theShape) - volumeOf(aFillet.Shape());
}
} // namespace

// The straight edge is one stripe with one surface data; the arc edge is another stripe with
// more. Filleting both at once must succeed and remove exactly what each fillet alone removes.
TEST(BRepFilletAPI_MakeFilletTest, DisjointStripesWithDifferentSurfDataCounts)
{
  const TopoDS_Shape aPrism    = makeArcPrism();
  const TopoDS_Edge  aStraight = findEdgeAt(aPrism, gp_Pnt(22.5, 11.25, 0.0));
  const TopoDS_Edge  anArc     = findEdgeAt(aPrism, gp_Pnt(-7.5, 11.25, 5.0));
  ASSERT_FALSE(aStraight.IsNull());
  ASSERT_FALSE(anArc.IsNull());
  const double aRemovedStraight = removedBy(aPrism, 2.5, aStraight);
  const double aRemovedArc      = removedBy(aPrism, 2.5, anArc);

  BRepFilletAPI_MakeFillet aFillet(aPrism);
  aFillet.Add(2.5, aStraight);
  aFillet.Add(2.5, anArc);
  ASSERT_NO_THROW(aFillet.Build());
  ASSERT_TRUE(aFillet.IsDone());
  EXPECT_TRUE(BRepCheck_Analyzer(aFillet.Shape()).IsValid());
  EXPECT_NEAR(volumeOf(aFillet.Shape()), volumeOf(aPrism) - aRemovedStraight - aRemovedArc, 1.e-3);
}
