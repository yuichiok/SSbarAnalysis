#ifndef GUARD_VectorTools_h
#define GUARD_VectorTools_h

/*------------------------------------------------------------------------------
   VectorTools
 Created : 2022-09-14  okugawa
 Main class of VectorTool program.
------------------------------------------------------------------------------*/

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

class VectorTools
{
  
  public:
    VectorTools();
    VectorTools( ROOT::Math::XYZVector vec3, Float_t E );
    VectorTools( Float_t x, Float_t y, Float_t z, Float_t E );
    virtual ~VectorTools() {};

    virtual void SetCoordinates(Float_t x, Float_t y, Float_t z, Float_t E);

    // return vectors
    virtual ROOT::Math::XYZVector     v3 ();
    virtual ROOT::Math::PxPyPzEVector v4 ();

  private:
    ROOT::Math::XYZVector     vector3;
    ROOT::Math::PxPyPzEVector vector4;

};

#endif