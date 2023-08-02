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
    // Constructor and Destructor
      VectorTools();
      VectorTools( ROOT::Math::XYZVector vec3, Float_t E );
      VectorTools( Float_t x, Float_t y, Float_t z, Float_t E );
      virtual ~VectorTools() {};

    // Set values
      virtual void SetCoordinates ( Float_t x, Float_t y, Float_t z, Float_t E );

    // Common Tools
      static Float_t GetCosBetween ( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 );
      static Float_t GetThetaBetween ( ROOT::Math::XYZVector v1, ROOT::Math::XYZVector v2 );

    // return private members
      virtual ROOT::Math::XYZVector     v3 ();
      virtual ROOT::Math::PxPyPzEVector v4 ();

  private:
    // Private members
      ROOT::Math::XYZVector     vector3;
      ROOT::Math::PxPyPzEVector vector4;

};

#endif