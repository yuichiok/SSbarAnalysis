
/*------------------------------------------------------------------------------
VectorTools.cpp
 Created : 2022-09-14  okugawa
------------------------------------------------------------------------------*/

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#include "VectorTools.hh"

VectorTools::VectorTools() {}

VectorTools::VectorTools( ROOT::Math::XYZVector vec3, Float_t E )
{
    vector3.SetCoordinates(vec3.x(),vec3.y(),vec3.z());
    vector4.SetCoordinates(vec3.x(),vec3.y(),vec3.z(),E);
}

VectorTools::VectorTools( Float_t x, Float_t y, Float_t z, Float_t E )
{
    vector3.SetCoordinates(x,y,z);
    vector4.SetCoordinates(x,y,z,E);
}

void VectorTools::SetCoordinates(Float_t x, Float_t y, Float_t z, Float_t E)
{
    vector3.SetCoordinates(x,y,z);
    vector4.SetCoordinates(x,y,z,E);
}

// Common tools
Float_t VectorTools::GetThetaBetween( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 )
{
  Float_t acol =  -2.0;
  Float_t dot = v1.Dot(v2);
  acol = std::acos( dot/(v1.R()*v2.R()) );
/*
  ROOT::Math::XYZVector v = v1.Cross(v2);
  Float_t sinth  = v.R()/(v1.R()*v2.R());
  Float_t direct = v1.Px()*v2.Px() * v1.Py()*v2.Py() * v1.Pz()*v2.Pz();
  sinth = (direct>0) ? sinth : -sinth;
  acol  = std::asin( sinth );
*/
//   acol = std::asin( (v.R()/(v1.R()*v2.R())) );// * v2.Mag()/(v1+v2).Mag(); 
  return acol;
}

ROOT::Math::XYZVector     VectorTools::v3() { return vector3; }
ROOT::Math::PxPyPzEVector VectorTools::v4() { return vector4; }