
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
Float_t VectorTools::GetCosBetween( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 )
{
  Float_t dot     = v1.Dot(v2);
  Float_t cos_btw = dot/(v1.R()*v2.R());
  return cos_btw;
}

Float_t VectorTools::GetSinACol( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 )
{
  ROOT::Math::XYZVector cross     = v1.Cross(v2);
  Float_t sin_acol = cross.R()/(v1.R()*v2.R());
  return sin_acol;
}

Float_t VectorTools::GetThetaBetween( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 )
{
  Float_t dot = v1.Dot(v2);
  return std::acos( dot/(v1.R()*v2.R()) );
}


ROOT::Math::XYZVector     VectorTools::v3() { return vector3; }
ROOT::Math::PxPyPzEVector VectorTools::v4() { return vector4; }