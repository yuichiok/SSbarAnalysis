
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

ROOT::Math::XYZVector     VectorTools::v3() { return vector3; }
ROOT::Math::PxPyPzEVector VectorTools::v4() { return vector4; }