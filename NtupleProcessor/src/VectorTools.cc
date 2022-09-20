
/*------------------------------------------------------------------------------
VectorTools.cpp
 Created : 2022-09-14  okugawa
------------------------------------------------------------------------------*/

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include "../include/VectorTools.hh"

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
Float_t VectorTools::GetSinacol( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 )
{
  Float_t sinacol =  -2.0;
  ROOT::Math::XYZVector v = v1.Cross(v2);
  sinacol = (v.R()/(v1.R()*v2.R()));// * v2.Mag()/(v1+v2).Mag(); 
  return sinacol;
}

ROOT::Math::XYZVector     VectorTools::v3() { return vector3; }
ROOT::Math::PxPyPzEVector VectorTools::v4() { return vector4; }