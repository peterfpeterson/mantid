// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <complex>
#include <fstream>
#include <iterator>
#include <list>
#include <map>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include "MantidKernel/Exception.h"

#include "MantidGeometry/Math/Triple.h"
#include "MantidGeometry/Objects/CSGObject.h"
#include "MantidGeometry/Objects/Rules.h"
#include "MantidGeometry/Surfaces/BaseVisit.h"
#include "MantidGeometry/Surfaces/Line.h"
#include "MantidGeometry/Surfaces/Surface.h"
#include "MantidKernel/Matrix.h"
#include "MantidKernel/V3D.h"

#ifdef ENABLE_OPENCASCADE
// Opencascade defines _USE_MATH_DEFINES without checking whether it is already
// used.
// Undefine it here before we include the headers to avoid a warning
#ifdef _MSC_VER
#undef _USE_MATH_DEFINES
#ifdef M_SQRT1_2
#undef M_SQRT1_2
#endif
#endif

#include "MantidKernel/WarningSuppressions.h"
GNU_DIAG_OFF("conversion")
GNU_DIAG_OFF("cast-qual")
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <TopoDS_Shape.hxx>
GNU_DIAG_ON("conversion")
GNU_DIAG_ON("cast-qual")
#endif

namespace Mantid::Geometry {
using Kernel::V3D;

Intersection::Intersection(std::unique_ptr<Rule> Ix, std::unique_ptr<Rule> Iy)
    : Rule(), A(std::move(Ix)), B(std::move(Iy))
/**
  Intersection constructor from two Rule ptrs.
  - Sets A,B's parents to *this
  - Allowed to control Ix and Iy
  @param Ix :: Rule A
  @param Iy :: Rule B
*/
{
  if (A)
    A->setParent(this);
  if (B)
    B->setParent(this);
}

Intersection::Intersection(Rule *Parent, std::unique_ptr<Rule> Ix, std::unique_ptr<Rule> Iy)
    : Rule(Parent), A(std::move(Ix)), B(std::move(Iy))
/**
  Intersection constructor from two Rule ptrs
  - Sets A,B's parents to *this.
  @param Parent :: Set the Rule::Parent pointer
  @param Ix :: Rule A
  @param Iy :: Rule B
*/
{
  if (A)
    A->setParent(this);
  if (B)
    B->setParent(this);
}

Intersection::Intersection(const Intersection &Iother)
    : Rule(), A(), B()
/**
  Copy constructor:
  Does a clone on the sub-tree below
  @param Iother :: Intersection to copy
*/
{
  if (Iother.A) {
    A = Iother.A->clone();
    A->setParent(this);
  }
  if (Iother.B) {
    B = Iother.B->clone();
    B->setParent(this);
  }
}

Intersection &Intersection::operator=(const Intersection &Iother)
/**
  Assignment operator :: Does a deep copy
  of the leaves of Iother.
  @param Iother :: object to copy
  @return *this
*/
{
  if (this != &Iother) {
    Rule::operator=(Iother);
    // first create new copy all fresh
    if (Iother.A) {
      A = Iother.A->clone();
      A->setParent(this);
    }
    if (Iother.B) {
      B = Iother.B->clone();
      B->setParent(this);
    }
  }
  return *this;
}

Intersection *Intersection::doClone() const
/**
  Virtual copy constructor
  @return new Intersection(this)
*/
{
  return new Intersection(*this);
}

std::unique_ptr<Intersection> Intersection::clone() const
/**
  Virtual copy constructor
  @return new Intersection(this)
*/
{
  return std::unique_ptr<Intersection>(doClone());
}

int Intersection::isComplementary() const
/**
  Determine is the rule has complementary
  sub components
  @retval 1 :: A side
  @retval -1 :: B side
  @retval 0 :: no complement
*/
{
  if (A && A->isComplementary())
    return 1;
  if (B && B->isComplementary())
    return -1;

  return 0;
}

void Intersection::setLeaves(std::unique_ptr<Rule> aR, std::unique_ptr<Rule> bR)
/**
  Replaces a both with a rule.
  No deletion is carried out but sets the parents.
  @param aR :: Rule on the left
  @param bR :: Rule on the right
*/
{
  A = std::move(aR);
  B = std::move(bR);
  if (A)
    A->setParent(this);
  if (B)
    B->setParent(this);
}

void Intersection::setLeaf(std::unique_ptr<Rule> nR, const int side)
/**
  Replaces a leaf with a rule.
  Calls delete on previous leaf.
  @param nR :: new rule
  @param side :: side to use
  - 0 == LHS
  - 1 == RHS
*/
{
  if (side) {
    B = std::move(nR);
    if (B)
      B->setParent(this);
  } else {
    A = std::move(nR);
    if (A)
      A->setParent(this);
  }
}

int Intersection::findLeaf(const Rule *R) const
/**
  Finds out if the Rule is the
  same as the leaves
  @param R :: Rule pointer to compare
  @retval 0 / 1 for LHS / RHS leaf
  @retval -1 :: neither leaf
*/
{
  if (A.get() == R)
    return 0;
  if (B.get() == R)
    return 1;
  return -1;
}

Rule *Intersection::findKey(const int KeyN)
/**
  Finds the leaf with the surface number KeyN
  @param KeyN :: Number to search for
  @retval 0 :: no leaf with that key number availiable
  @retval Rule* if an appropiate leaf is found
*/
{
  Rule *PtrOut = (A) ? A->findKey(KeyN) : nullptr;
  if (PtrOut)
    return PtrOut;
  return (B) ? B->findKey(KeyN) : nullptr;
}

std::string Intersection::display() const
/**
  Displaces a bracket wrapped object
  @return Bracketed string
*/
{
  std::string out;
  if (!A || !B)
    throw std::runtime_error("Intersection::display incomplete type");
  if (A->type() == -1)
    out = "(" + A->display() + ")";
  else
    out = A->display();

  out += " ";

  if (B->type() == -1)
    out += "(" + B->display() + ")";
  else
    out += B->display();
  return out;
}

std::string Intersection::displayAddress() const
/**
  Debug function that converts the
  the intersection ion space delimited unit
  to denote intersection.
  @return ( leaf leaf )
*/
{
  std::stringstream cx;
  cx << " [ " << this;
  if (A && B)
    cx << " ] (" + A->displayAddress() + " " + B->displayAddress() + ") ";
  else if (A)
    cx << " ] (" + A->displayAddress() + " 0x0 ) ";
  else if (B)
    cx << " ] ( 0x0 " + B->displayAddress() + ") ";
  else
    cx << " ] ( 0x0 0x0 ) ";
  return cx.str();
}

bool Intersection::isValid(const Kernel::V3D &Vec) const
/**
  Calculates if Vec is within the object
  @param Vec :: Point to test
  @retval 1 ::  Vec is within object
  @retval 0 :: Vec is outside object.
*/
{
  if (A && B) {
    return (A->isValid(Vec) && B->isValid(Vec));
  }
  return false;
}

bool Intersection::isValid(const std::map<int, int> &MX) const
/**
  Use MX to determine if the surface truth etc is
  valie
  @param MX :: map of key + logical value XOR sign
  @retval 1 ::  Both sides are valid
  @retval 0 :: Either side is invalid.
*/
{
  if (!A || !B)
    return false;
  return A->isValid(MX) && B->isValid(MX);
}

int Intersection::simplify()
/**
  Union simplification::
    - -S S simplify to True.
    - S S simplify to S
    - -S -S simplifies to -S
    @retval 1 if clauses removed (not top)
    @retval -1 replacement of this intersection is required
                by leaf 0
    @retval 0 if no work to do.
*/
{
  return 0;
}

/**
 * find the common bounding box with the two childs of intersection
 * @param xmax :: Maximum value for the bounding box in x direction
 * @param ymax :: Maximum value for the bounding box in y direction
 * @param zmax :: Maximum value for the bounding box in z direction
 * @param xmin :: Minimum value for the bounding box in x direction
 * @param ymin :: Minimum value for the bounding box in y direction
 * @param zmin :: Minimum value for the bounding box in z direction
 */
void Intersection::getBoundingBox(double &xmax, double &ymax, double &zmax, double &xmin, double &ymin, double &zmin) {
  double Axmax, Aymax, Azmax, Axmin, Aymin, Azmin;
  double Bxmax, Bymax, Bzmax, Bxmin, Bymin, Bzmin;
  Axmax = Bxmax = xmax;
  Aymax = Bymax = ymax;
  Azmax = Bzmax = zmax;
  Axmin = Bxmin = xmin;
  Aymin = Bymin = ymin;
  Azmin = Bzmin = zmin;
  A->getBoundingBox(Axmax, Aymax, Azmax, Axmin, Aymin, Azmin);
  B->getBoundingBox(Bxmax, Bymax, Bzmax, Bxmin, Bymin, Bzmin);
  xmax = (Axmax < Bxmax) ? Axmax : Bxmax;
  xmin = (Axmin > Bxmin) ? Axmin : Bxmin;
  ymax = (Aymax < Bymax) ? Aymax : Bymax;
  ymin = (Aymin > Bymin) ? Aymin : Bymin;
  zmax = (Azmax < Bzmax) ? Azmax : Bzmax;
  zmin = (Azmin > Bzmin) ? Azmin : Bzmin;
}

#ifdef ENABLE_OPENCASCADE
/**
 * Analyze intersection
 * @return the resulting TopoDS_Shape
 */
TopoDS_Shape Intersection::analyze() {
  TopoDS_Shape left = A->analyze();
  TopoDS_Shape right = B->analyze();
  BRepAlgoAPI_Common comm(left, right);
  return comm.Shape();
}
#endif

// -------------------------------------------------------------
//         UNION
//---------------------------------------------------------------

Union::Union(Rule *Parent, std::unique_ptr<Rule> Ix, std::unique_ptr<Rule> Iy)
    : Rule(Parent), A(std::move(Ix)), B(std::move(Iy))
/**
  Union constructor from two Rule ptrs.
  - Sets A,B's parents to *this
  - Allowed to control Ix and Iy
  @param Parent :: Rule that is the parent to this
  @param Ix :: Rule A
  @param Iy :: Rule B
*/
{
  if (A)
    A->setParent(this);
  if (B)
    B->setParent(this);
}

Union::Union(std::unique_ptr<Rule> Ix, std::unique_ptr<Rule> Iy)
    : Rule(), A(std::move(Ix)), B(std::move(Iy))
/**
  Union constructor from two Rule ptrs
  - Sets A,B's parents to *this.
  - Allowed to control Ix and Iy
  @param Ix :: Rule A
  @param Iy :: Rule B
*/
{
  if (A)
    A->setParent(this);
  if (B)
    B->setParent(this);
}

Union::Union(const Union &Iother)
    : Rule(Iother), A(), B()
/**
  Copy constructor:
  Does a clone on the sub-tree below
  @param Iother :: Union to copy
*/
{
  if (Iother.A) {
    A = Iother.A->clone();
    A->setParent(this);
  }
  if (Iother.B) {
    B = Iother.B->clone();
    B->setParent(this);
  }
}

Union &Union::operator=(const Union &Iother)
/**
  Assignment operator :: Does a deep copy
  of the leaves of Iother.
  @param Iother :: Union to assign to it.
  @return this union (copy).
*/
{
  if (this != &Iother) {
    Rule::operator=(Iother);
    // first create new copy all fresh
    if (Iother.A) {
      A = Iother.A->clone();
      A->setParent(this);
    }
    if (Iother.B) {
      B = Iother.B->clone();
      B->setParent(this);
    }
  }
  return *this;
}

Union *Union::doClone() const
/**
  Clone allows deep virtual coping
  @return new Union copy.
*/
{
  return new Union(*this);
}

std::unique_ptr<Union> Union::clone() const
/**
  Clone allows deep virtual coping
  @return new Union copy.
*/
{
  return std::unique_ptr<Union>(doClone());
}

void Union::setLeaf(std::unique_ptr<Rule> nR, const int side)
/**
  Replaces a leaf with a rule.
  Calls delete on previous leaf.
  @param nR :: new rule
  @param side :: side to use
  - 0 == LHS
  - 1 == RHS
*/
{
  if (side) {
    B = std::move(nR);
    if (B)
      B->setParent(this);
  } else {
    A = std::move(nR);
    if (A)
      A->setParent(this);
  }
  return;
}

void Union::setLeaves(std::unique_ptr<Rule> aR, std::unique_ptr<Rule> bR)
/**
   Replaces a both with a rule.
  No deletion is carried out but sets the parents.
  @param aR :: Rule on the left
  @param bR :: Rule on the right
*/
{
  A = std::move(aR);
  B = std::move(bR);
  if (A)
    A->setParent(this);
  if (B)
    B->setParent(this);
}

int Union::findLeaf(const Rule *R) const
/**
  Finds out if the Rule is the
  same as the leaves
  @param R :: Rule pointer to compare
  @retval 0 / 1 for LHS / RHS leaf
  @retval -1 :: neither leaf
*/
{
  if (A.get() == R)
    return 0;
  if (B.get() == R)
    return 1;
  return -1;
}

Rule *Union::findKey(const int KeyN)
/**
  Finds the leaf with the surface number KeyN
  @param KeyN :: Number to search for
  @retval 0 :: no leaf with that key number availiable
  @retval Rule* if an appropiate leaf is found
*/
{

  Rule *PtrOut = (A) ? A->findKey(KeyN) : nullptr;
  if (PtrOut)
    return PtrOut;
  return (B) ? B->findKey(KeyN) : nullptr;
}

int Union::isComplementary() const
/**
  Determine is the rule has complementary
  sub components
  @retval 1 :: A side
  @retval -1 :: B side
  @retval 0 :: no complement
*/
{
  if (A && A->isComplementary())
    return 1;
  if (B && B->isComplementary())
    return -1;

  return 0;
}

int Union::simplify()
/**
  Union simplification::
    - -S S simplify to True.
    - S S simplify to S
    - -S -S simplifies to -S
    @retval 1 if clauses removed (not top)
    @retval -1 replacement of this intersection is required
                by leaf 0
    @retval 0 if no work to do.
*/
{
  return 0;
}

bool Union::isValid(const Kernel::V3D &Vec) const
/**
  Calculates if Vec is within the object
  @param Vec :: Point to test
  @retval  1 ::  Vec is within object
  @retval 0 :: Vec is outside object.
*/
{
  return (A && A->isValid(Vec)) || (B && B->isValid(Vec));
}

bool Union::isValid(const std::map<int, int> &MX) const
/**
  Use MX to determine if the surface truth etc is
  valie
  @param MX :: map of key + logical value XOR sign
  @retval 1 :: if either side is valid
  @retval 0 :: Neither side is valid
*/
{
  return (A && A->isValid(MX)) || (B && B->isValid(MX));
}

std::string Union::display() const
/**
  Display the union in the form
  (N:M) where N,M are the downward
  rules
  @return bracket string
*/
{
  std::string out;
  if (!A || !B)
    throw std::runtime_error("Union::display incomplete type");
  if (A->type() == 1)
    out = "(" + A->display() + ")";
  else
    out = A->display();

  out += " : ";

  if (B->type() == 1)
    out += "(" + B->display() + ")";
  else
    out += B->display();

  return out;
}

std::string Union::displayAddress() const
/**
  Returns the memory address as
  a string. Displays addresses of leaves
  @return String of address
*/
{
  std::stringstream cx;
  cx << " [ " << this;

  if (A && B)
    cx << " ] (" + A->displayAddress() + " : " + B->displayAddress() + ") ";
  else if (A)
    cx << " ] (" + A->displayAddress() + " : 0x0 ) ";
  else if (B)
    cx << " ] ( 0x0 : " + B->displayAddress() + ") ";
  else
    cx << " ] ( 0x0 : 0x0 ) ";
  return cx.str();
}

/**
 * gets the bounding box for the Union Rule
 * @param xmax :: Maximum value for the bounding box in x direction
 * @param ymax :: Maximum value for the bounding box in y direction
 * @param zmax :: Maximum value for the bounding box in z direction
 * @param xmin :: Minimum value for the bounding box in x direction
 * @param ymin :: Minimum value for the bounding box in y direction
 * @param zmin :: Minimum value for the bounding box in z direction
 */
void Union::getBoundingBox(double &xmax, double &ymax, double &zmax, double &xmin, double &ymin, double &zmin) {
  double Axmax, Aymax, Azmax, Axmin, Aymin, Azmin;
  double Bxmax, Bymax, Bzmax, Bxmin, Bymin, Bzmin;
  Axmax = Bxmax = xmax;
  Aymax = Bymax = ymax;
  Azmax = Bzmax = zmax;
  Axmin = Bxmin = xmin;
  Aymin = Bymin = ymin;
  Azmin = Bzmin = zmin;
  A->getBoundingBox(Axmax, Aymax, Azmax, Axmin, Aymin, Azmin);
  B->getBoundingBox(Bxmax, Bymax, Bzmax, Bxmin, Bymin, Bzmin);
  xmax = (Axmax > Bxmax) ? Axmax : Bxmax;
  xmin = (Axmin < Bxmin) ? Axmin : Bxmin;
  ymax = (Aymax > Bymax) ? Aymax : Bymax;
  ymin = (Aymin < Bymin) ? Aymin : Bymin;
  zmax = (Azmax > Bzmax) ? Azmax : Bzmax;
  zmin = (Azmin < Bzmin) ? Azmin : Bzmin;
}

#ifdef ENABLE_OPENCASCADE
TopoDS_Shape Union::analyze() {
  TopoDS_Shape left = A->analyze();
  TopoDS_Shape right = B->analyze();
  BRepAlgoAPI_Fuse fuse(left, right);
  return fuse.Shape();
}
#endif

// -------------------------------------------------------------
//         SURF KEYS
//---------------------------------------------------------------

SurfPoint::SurfPoint()
    : Rule(), m_key(), keyN(0), sign(1)
/**
  Constructor with null key/number
*/
{}

SurfPoint *SurfPoint::doClone() const
/**
  Clone constructor
  @return new(*this)
*/
{
  return new SurfPoint(*this);
}

std::unique_ptr<SurfPoint> SurfPoint::clone() const
/**
 Clone constructor
 @return new(*this)
 */
{
  return std::unique_ptr<SurfPoint>(doClone());
}

void SurfPoint::setLeaf(std::unique_ptr<Rule> nR, const int /*unused*/)
/**
  Replaces a leaf with a rule.
  This REQUIRES that nR is of type SurfPoint
  @param nR :: new rule
  @param :: ignored
*/
{
  // std::cerr<<"Calling SurfPoint setLeaf"<<'\n';

  const auto *newX = dynamic_cast<const SurfPoint *>(nR.get());
  if (newX)
    *this = *newX;
}

void SurfPoint::setLeaves(std::unique_ptr<Rule> aR, std::unique_ptr<Rule> /*unused*/)
/**
  Replaces a leaf with a rule.
  This REQUIRES that nR is of type SurfPoint
  @param aR :: new rule
*/
{
  // std::cerr<<"Calling SurfPoint setLeaf"<<'\n';
  const auto *newX = dynamic_cast<const SurfPoint *>(aR.get());
  if (newX)
    *this = *newX;
}

int SurfPoint::findLeaf(const Rule *A) const
/**
  Determines if this rule is a particular leaf value
  uses memory address to compare.
  @param A :: Rule to check
  @return 0 if true and -1 if false.
*/
{
  return (this == A) ? 0 : -1;
}

Rule *SurfPoint::findKey(const int KeyNum)
/**
  Finds the leaf with the surface number KeyN
  @param KeyNum :: Number to search for
  @retval 0 :: no leaf with that key number availiable
  @retval Rule* if an appropiate leaf is found
*/
{
  return (KeyNum == keyN) ? this : nullptr;
}

void SurfPoint::setKeyN(const int Ky)
/**
  Sets the key and the sign
  @param Ky :: key value (+/- is used for sign)
*/
{
  sign = (Ky < 0) ? -1 : 1;
  keyN = sign * Ky;
}

void SurfPoint::setKey(const std::shared_ptr<Surface> &Spoint)
/**
  Sets the key pointer. The class takes ownership.
  @param Spoint :: new key values
*/
{
  m_key = Spoint;
}

int SurfPoint::simplify()
/**
  Impossible to simplify a simple
  rule leaf. Therefore returns 0
  @return 0
*/
{
  return 0;
}

bool SurfPoint::isValid(const Kernel::V3D &Pt) const
/**
  Determines if a point  is valid.
  @param Pt :: Point to test
  @retval 1 :: Pt is the +ve side of the
  surface or on the surface
  @retval 0 :: Pt is on the -ve side of the surface
*/
{
  if (m_key) {
    return (m_key->side(Pt) * sign) >= 0;
  }
  return false;
}

bool SurfPoint::isValid(const std::map<int, int> &MX) const
/**
  Use MX to determine if the surface truth etc is
  valid
  @param MX :: map of key + logical value XOR sign
  @return MX.second if key found or 0
*/
{
  auto lx = MX.find(keyN);
  if (lx == MX.end())
    return false;
  const int rtype = (lx->second) ? 1 : -1;
  return (rtype * sign) >= 0;
}

std::string SurfPoint::display() const
/**
  Returns the signed surface number as
  a string.
  @return string of the value
*/
{
  std::stringstream cx;
  cx << sign * keyN;
  return cx.str();
}

std::string SurfPoint::displayAddress() const
/**
  Returns the memory address as
  a string.
  @return memory address as string
*/
{
  std::stringstream cx;
  cx << this;
  return cx.str();
}

/**
 * gets the bounding box for the surface object held by SurfPoint
 * @param xmax :: Maximum value for the bounding box in x direction
 * @param ymax :: Maximum value for the bounding box in y direction
 * @param zmax :: Maximum value for the bounding box in z direction
 * @param xmin :: Minimum value for the bounding box in x direction
 * @param ymin :: Minimum value for the bounding box in y direction
 * @param zmin :: Minimum value for the bounding box in z direction
 */
void SurfPoint::getBoundingBox(double &xmax, double &ymax, double &zmax, double &xmin, double &ymin, double &zmin) {
  if (this->sign < 1) // If the object sign is positive then include
    m_key->getBoundingBox(xmax, ymax, zmax, xmin, ymin, zmin);
  else { // if the object sign is negative then get the complement
    std::vector<V3D> listOfPoints;
    double gXmax, gYmax, gZmax, gXmin, gYmin, gZmin;
    gXmax = xmax;
    gYmax = ymax;
    gZmax = zmax;
    gXmin = xmin;
    gYmin = ymin;
    gZmin = zmin;
    m_key->getBoundingBox(gXmax, gYmax, gZmax, gXmin, gYmin, gZmin);
    if (!((xmax <= gXmax && xmax >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
      listOfPoints.emplace_back(xmax, ymax, zmax);
    if (!((xmin <= gXmax && xmin >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
      listOfPoints.emplace_back(xmin, ymax, zmax);
    if (!((xmin <= gXmax && xmin >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
      listOfPoints.emplace_back(xmin, ymax, zmin);
    if (!((xmax <= gXmax && xmax >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
      listOfPoints.emplace_back(xmax, ymax, zmin);
    if (!((xmin <= gXmax && xmin >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
      listOfPoints.emplace_back(xmin, ymin, zmin);
    if (!((xmax <= gXmax && xmax >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
      listOfPoints.emplace_back(xmax, ymin, zmin);
    if (!((xmax <= gXmax && xmax >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
      listOfPoints.emplace_back(xmax, ymin, zmax);
    if (!((xmin <= gXmax && xmin >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
      listOfPoints.emplace_back(xmin, ymin, zmax);

    // group box inside input box
    if (((gXmax <= xmax && gXmax >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
        (gXmax != xmax || gYmax != ymax || gZmax != zmax))
      listOfPoints.emplace_back(gXmax, gYmax, gZmax);
    if (((gXmin <= xmax && gXmin >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
        (gXmin != xmin || gYmax != ymax || gZmax != zmax))
      listOfPoints.emplace_back(gXmin, gYmax, gZmax);
    if (((gXmin <= xmax && gXmin >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
        (gXmin != xmin || gYmax != ymax || gZmin != zmin))
      listOfPoints.emplace_back(gXmin, gYmax, gZmin);
    if (((gXmax <= xmax && gXmax >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
        (gXmax != xmax || gYmax != ymax || gZmin != zmin))
      listOfPoints.emplace_back(gXmax, gYmax, gZmin);
    if (((gXmin <= xmax && gXmin >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
        (gXmin != xmin || gYmin != ymin || gZmin != zmin))
      listOfPoints.emplace_back(gXmin, gYmin, gZmin);
    if (((gXmax <= xmax && gXmax >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
        (gXmax != xmax || gYmin != ymin || gZmin != zmin))
      listOfPoints.emplace_back(gXmax, gYmin, gZmin);
    if (((gXmax <= xmax && gXmax >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
        (gXmax != xmax || gYmin != ymin || gZmax != zmax))
      listOfPoints.emplace_back(gXmax, gYmin, gZmax);
    if (((gXmin <= xmax && gXmin >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
        (gXmin != xmin || gYmin != ymin || gZmax != zmax))
      listOfPoints.emplace_back(gXmin, gYmin, gZmax);

    if (!listOfPoints.empty()) {
      xmin = ymin = zmin = DBL_MAX;
      xmax = ymax = zmax = -DBL_MAX;
      for (std::vector<V3D>::const_iterator it = listOfPoints.begin(); it != listOfPoints.end(); ++it) {
        //			std::cout<<(*it)<<'\n';
        if ((*it)[0] < xmin)
          xmin = (*it)[0];
        if ((*it)[1] < ymin)
          ymin = (*it)[1];
        if ((*it)[2] < zmin)
          zmin = (*it)[2];
        if ((*it)[0] > xmax)
          xmax = (*it)[0];
        if ((*it)[1] > ymax)
          ymax = (*it)[1];
        if ((*it)[2] > zmax)
          zmax = (*it)[2];
      }
    }
  }
}

#ifdef ENABLE_OPENCASCADE
TopoDS_Shape SurfPoint::analyze() {
  // Check for individual type of surfaces
  TopoDS_Shape Result = m_key->createShape();
  if (sign > 0)
    Result.Complement();
  if (m_key->className() == "Plane") {
    // build a box
    gp_Pnt p(-1000.0, -1000.0, -1000.0);
    TopoDS_Shape world = BRepPrimAPI_MakeBox(p, 2000.0, 2000.0, 2000.0).Shape();
    return BRepAlgoAPI_Common(world, Result);
  }
  return Result;
}
#endif
//----------------------------------------
//       COMPOBJ
//----------------------------------------

CompObj::CompObj()
    : Rule(), objN(0), key(nullptr)
/**
  Constructor
*/
{}

CompObj *CompObj::doClone() const
/**
  Clone of this
  @return new copy of this
*/
{
  return new CompObj(*this);
}

std::unique_ptr<CompObj> CompObj::clone() const
/**
 Clone of this
 @return new copy of this
*/
{
  return std::unique_ptr<CompObj>(doClone());
}

void CompObj::setObjN(const int Ky)
/**
  Sets the object Number
  @param Ky :: key value
*/
{
  objN = Ky;
}

void CompObj::setObj(CSGObject *val)
/**
  Sets the object
  @param val :: Object value
*/
{
  key = val;
}

void CompObj::setLeaf(std::unique_ptr<Rule> aR, const int /*unused*/)
/**
  Replaces a leaf with a rule.
  This REQUIRES that aR is of type SurfPoint
  @param aR :: new rule
  @param :: Null side point
*/
{
  const auto *newX = dynamic_cast<const CompObj *>(aR.get());
  // Make a copy
  if (newX)
    *this = *newX;
}

void CompObj::setLeaves(std::unique_ptr<Rule> aR, std::unique_ptr<Rule> oR)
/**
  Replaces a leaf with a rule.
  This REQUIRES that aR is of type CompObj
  @param aR :: new rule
  @param oR :: Null other rule
*/
{
  (void)oR; // Avoid compiler warning

  const auto *newX = dynamic_cast<const CompObj *>(aR.get());
  if (newX)
    *this = *newX;
}

Rule *CompObj::findKey(const int i)
/**
  This is a complementary object and we dont
  search into CompObjs. If that is needed
  then the CompObj should be removed first
  @param i :: Null index key
  @return 0
*/
{
  (void)i; // Avoid compiler warning
  return nullptr;
}

int CompObj::findLeaf(const Rule *A) const
/**
  Check to see if this is a copy of a given Rule
  @param A :: Rule Ptr to find
  @retval 0 on success -ve on failuire
*/
{
  return (this == A) ? 0 : -1;
}

bool CompObj::isValid(const Kernel::V3D &Pt) const
/**
  Determines if a point  is valid.
  Checks to see if the point is valid in the object
  and returns ture if it is not valid.
  @param  Pt :: Point to test
  @retval not valid in the object
  @retval true is no object is set
*/
{
  if (key)
    return !(key->isValid(Pt));
  return true;
}

bool CompObj::isValid(const std::map<int, int> &SMap) const
/**
  Determines if a point  is valid.
  @param  SMap :: map of surface number and true values
  @return :: status
*/
{
  return (key) ? !(key->isValid(SMap)) : true;
}

int CompObj::simplify()
/**
  Impossible to simplify a simple
  rule leaf. Therefore returns 0
  @return 0
*/
{
  return 0;
}

std::string CompObj::display() const
/**
  Displays the object as \#number
  @return string component
*/
{
  std::stringstream cx;
  cx << "#" << objN;
  return cx.str();
}

std::string CompObj::displayAddress() const
/**
  Returns the memory address as
  a string.
  @return memory address as string
*/
{
  std::stringstream cx;
  cx << this;
  return cx.str();
}

/**
 * gets the bounding box for the Complementary of the object
 * @param xmax :: Maximum value for the bounding box in x direction
 * @param ymax :: Maximum value for the bounding box in y direction
 * @param zmax :: Maximum value for the bounding box in z direction
 * @param xmin :: Minimum value for the bounding box in x direction
 * @param ymin :: Minimum value for the bounding box in y direction
 * @param zmin :: Minimum value for the bounding box in z direction
 */
void CompObj::getBoundingBox(double &xmax, double &ymax, double &zmax, double &xmin, double &ymin, double &zmin) {
  /// To calculate the bounding box is formed by the all the points in input box
  /// that are not in object bounding box
  /// and all the points in object bounding box that are in the input box
  std::vector<V3D> listOfPoints;
  double gXmax, gYmax, gZmax, gXmin, gYmin, gZmin;
  gXmax = xmax;
  gYmax = ymax;
  gZmax = zmax;
  gXmin = xmin;
  gYmin = ymin;
  gZmin = zmin;
  key->getBoundingBox(gXmax, gYmax, gZmax, gXmin, gYmin, gZmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmax, ymax, zmax);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmin, ymax, zmax);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmin, ymax, zmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmax, ymax, zmin);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmin, ymin, zmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmax, ymin, zmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmax, ymin, zmax);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmin, ymin, zmax);

  // object box inside input box
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmax != xmax || gYmax != ymax || gZmax != zmax))
    listOfPoints.emplace_back(gXmax, gYmax, gZmax);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmin != xmin || gYmax != ymax || gZmax != zmax))
    listOfPoints.emplace_back(gXmin, gYmax, gZmax);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmin != xmin || gYmax != ymax || gZmin != zmin))
    listOfPoints.emplace_back(gXmin, gYmax, gZmin);
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmax != xmax || gYmax != ymax || gZmin != zmin))
    listOfPoints.emplace_back(gXmax, gYmax, gZmin);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmin != xmin || gYmin != ymin || gZmin != zmin))
    listOfPoints.emplace_back(gXmin, gYmin, gZmin);
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmax != xmax || gYmin != ymin || gZmin != zmin))
    listOfPoints.emplace_back(gXmax, gYmin, gZmin);
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmax != xmax || gYmin != ymin || gZmax != zmax))
    listOfPoints.emplace_back(gXmax, gYmin, gZmax);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmin != xmin || gYmin != ymin || gZmax != zmax))
    listOfPoints.emplace_back(gXmin, gYmin, gZmax);

  if (!listOfPoints.empty()) {
    xmin = ymin = zmin = DBL_MAX;
    xmax = ymax = zmax = -DBL_MAX;
    for (std::vector<V3D>::const_iterator it = listOfPoints.begin(); it != listOfPoints.end(); ++it) {
      //			std::cout<<(*it)<<'\n';
      if ((*it)[0] < xmin)
        xmin = (*it)[0];
      if ((*it)[1] < ymin)
        ymin = (*it)[1];
      if ((*it)[2] < zmin)
        zmin = (*it)[2];
      if ((*it)[0] > xmax)
        xmax = (*it)[0];
      if ((*it)[1] > ymax)
        ymax = (*it)[1];
      if ((*it)[2] > zmax)
        zmax = (*it)[2];
    }
  }
}

#ifdef ENABLE_OPENCASCADE
TopoDS_Shape CompObj::analyze() {
  TopoDS_Shape Result = const_cast<Rule *>(key->topRule())->analyze();
  Result.Complement();
  return Result;
}
#endif
// -----------------------------------------------
// BOOLVALUE
// -----------------------------------------------

BoolValue::BoolValue()
    : Rule(), status(-1)
/**
  Constructor
*/
{}

BoolValue &BoolValue::operator=(const BoolValue &A)
/**
  Assignment operator
  @param A :: BoolValue to copy
  @return *this
 */
{
  if (this != &A) {
    Rule::operator=(A);
    status = A.status;
  }
  return *this;
}

BoolValue *BoolValue::doClone() const
/**
  Clone constructor
  @return new(*this)
 */
{
  return new BoolValue(*this);
}

std::unique_ptr<BoolValue> BoolValue::clone() const
/**
  Clone constructor
  @return new(*this)
*/
{
  return std::unique_ptr<BoolValue>(doClone());
}

void BoolValue::setLeaf(std::unique_ptr<Rule> aR, const int /*unused*/)
/**
  Replaces a leaf with a rule.
  This REQUIRES that aR is of type SurfPoint
  @param aR :: new rule

  secont argument ignored ( Null side point)
*/
{
  // std::cerr<<"Calling BoolValue setLeaf"<<'\n';
  const auto *newX = dynamic_cast<const BoolValue *>(aR.get());
  if (newX)
    *this = *newX;
}

void BoolValue::setLeaves(std::unique_ptr<Rule> aR, std::unique_ptr<Rule> oR)
/**
  Replaces a leaf with a rule.
  This REQUIRES that aR is of type SurfPoint
  @param aR :: new rule
  @param oR :: Null other rule
*/
{
  (void)oR; // Avoid compiler warning
  // std::cerr<<"Calling BoolValue setLeaves"<<'\n';
  const auto *newX = dynamic_cast<const BoolValue *>(aR.get());
  if (newX)
    *this = *newX;
}

int BoolValue::findLeaf(const Rule *A) const { return (this == A) ? 0 : -1; }

bool BoolValue::isValid(const Kernel::V3D &pt) const
/**
  Determines if a point  is valid.
  @param pt :: Point to test
  @return status
*/
{
  (void)pt; // Avoid compiler warning
  return status > 0;
}

bool BoolValue::isValid(const std::map<int, int> &map) const
/**
  Determines if a point  is valid.
  @param map :: map of surface number and true values
  @return :: status
*/
{
  (void)map; // Avoid compiler warning
  return status > 0;
}

int BoolValue::simplify()
/**
  Bool value is always in simplest form
  @return 0
*/
{
  return 0;
}

std::string BoolValue::display() const
/**
  @return string of
   - "true" "false" or "unknown"
*/

{
  switch (status) {
  case 1:
    return " True ";
  case 0:
    return " False ";
  }
  return " Unknown ";
}

std::string BoolValue::displayAddress() const
/**
  Returns the memory address as
  a string.
  @return string of Address
*/
{
  std::stringstream cx;
  cx << this;
  return cx.str();
}

/**
 * gets the bounding box for the BoolValue Rule
 * @param xmax :: Maximum value for the bounding box in x direction
 * @param ymax :: Maximum value for the bounding box in y direction
 * @param zmax :: Maximum value for the bounding box in z direction
 * @param xmin :: Minimum value for the bounding box in x direction
 * @param ymin :: Minimum value for the bounding box in y direction
 * @param zmin :: Minimum value for the bounding box in z direction
 */
void BoolValue::getBoundingBox(double &xmax, double &ymax, double &zmax, double &xmin, double &ymin, double &zmin) {
  // Returns what ever bounding box it gets
  (void)xmax; // Avoid compiler warning
  (void)ymax; // Avoid compiler warning
  (void)zmax; // Avoid compiler warning
  (void)xmin; // Avoid compiler warning
  (void)ymin; // Avoid compiler warning
  (void)zmin; // Avoid compiler warning
}

//----------------------------------------
//       COMPGRP
//----------------------------------------

CompGrp::CompGrp(Rule *Parent, std::unique_ptr<Rule> Cx)
    : Rule(Parent), A(std::move(Cx))
/**
  Constructor to build parent and complent tree
  @param Parent :: Rule that is the parent to this
  @param Cx :: Complementary object below
*/
{
  if (A)
    A->setParent(this);
}

CompGrp::CompGrp(const CompGrp &Cother)
    : Rule(Cother), A()
/**
  Standard copy constructor
  @param Cother :: CompGrp to copy
 */
{
  if (Cother.A) {
    A = std::unique_ptr<Rule>(Cother.A->clone());
    A->setParent(this);
  }
}

CompGrp &CompGrp::operator=(const CompGrp &Cother)
/**
  Standard assignment operator
  @param Cother :: CompGrp to copy
  @return *this
 */
{
  if (this != &Cother) {
    Rule::operator=(Cother);
    if (Cother.A) {
      A = std::unique_ptr<Rule>(Cother.A->clone());
      A->setParent(this);
    }
  }
  return *this;
}

CompGrp *CompGrp::doClone() const
/**
  Clone of this
  @return new copy of this
*/
{
  return new CompGrp(*this);
}

std::unique_ptr<CompGrp> CompGrp::clone() const
/**
  Clone of this
  @return new copy of this
*/
{
  return std::unique_ptr<CompGrp>(doClone());
}

void CompGrp::setLeaf(std::unique_ptr<Rule> nR, const int side)
/**
  Replaces a leaf with a rule.
  No deletion is carried out
  @param nR :: new rule
  @param side :: side to use
*/
{
  (void)side; // Avoid compiler warning
  A = std::move(nR);
  if (A)
    A->setParent(this);
}

void CompGrp::setLeaves(std::unique_ptr<Rule> aR, std::unique_ptr<Rule> oR)
/**
  Replaces a leaf with a rule.
  No deletion is carried out but sets the parents.
  @param aR :: new rule
  @param oR :: Null other rule
*/
{
  (void)oR; // Avoid compiler warning
  A = std::move(aR);
  if (A)
    A->setParent(this);
}

Rule *CompGrp::findKey(const int i)
/**
  This is a complementary object and we dont
  search into CompGrps. If that is needed
  then the CompGrp should be removed first
  @param i :: Null index key
  @return 0
*/
{
  (void)i; // Avoid compiler warning
  return nullptr;
}

int CompGrp::findLeaf(const Rule *R) const
/**
  Check to see if this is a copy of a given Rule
  @param R :: Rule Ptr to find
  @retval 0 on success -ve on failuire
*/
{
  return (A.get() == R) ? 0 : -1;
}

bool CompGrp::isValid(const Kernel::V3D &Pt) const
/**
  Determines if a point  is valid.
  Checks to see if the point is valid in the object
  and returns ture if it is not valid.
  Note the complementary reverse in the test.
  @param  Pt :: Point to test
  @retval not valid in the object
  @retval true is no object is set
*/
{
  // Note:: if isValid is true then return 0:
  if (A)
    return !(A->isValid(Pt));
  return true;
}

bool CompGrp::isValid(const std::map<int, int> &SMap) const
/**
  Determines if a point  is valid.
  @param  SMap :: map of surface number and true values
  @return :: status
*/
{
  // Note:: if isValid is true then return 0:
  if (A)
    return !A->isValid(SMap);
  return true;
}

int CompGrp::simplify()
/**
  Impossible to simplify a simple
  rule leaf. Therefore returns 0
  @return 0
*/
{
  return 0;
}

std::string CompGrp::display() const
/**
  Displays the object as \#number
  @return string component
*/
{
  std::stringstream cx;
  if (A)
    cx << "#( " << A->display() << " )";
  return cx.str();
}

std::string CompGrp::displayAddress() const
/**
  Returns the memory address as
  a string.
  @return memory address as string
*/
{
  std::stringstream cx;
  cx << "#( [" << this << "] ";
  if (A)
    cx << A->displayAddress();
  else
    cx << "0x0";
  cx << " ) ";
  return cx.str();
}

/**
 * gets the bounding box for the complement of the group
 * @param xmax :: Maximum value for the bounding box in x direction
 * @param ymax :: Maximum value for the bounding box in y direction
 * @param zmax :: Maximum value for the bounding box in z direction
 * @param xmin :: Minimum value for the bounding box in x direction
 * @param ymin :: Minimum value for the bounding box in y direction
 * @param zmin :: Minimum value for the bounding box in z direction
 */
void CompGrp::getBoundingBox(double &xmax, double &ymax, double &zmax, double &xmin, double &ymin, double &zmin) {
  /// To calculate the bounding box is formed by the all the points in input box
  /// that are not in group bounding box
  /// and all the points in group bounding box that are in the input box
  std::vector<V3D> listOfPoints;
  double gXmax, gYmax, gZmax, gXmin, gYmin, gZmin;
  gXmax = xmax;
  gYmax = ymax;
  gZmax = zmax;
  gXmin = xmin;
  gYmin = ymin;
  gZmin = zmin;
  A->getBoundingBox(gXmax, gYmax, gZmax, gXmin, gYmin, gZmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmax, ymax, zmax);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmin, ymax, zmax);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmin, ymax, zmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymax <= gYmax && ymax >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmax, ymax, zmin);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmin, ymin, zmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmin <= gZmax && zmin >= gZmin)))
    listOfPoints.emplace_back(xmax, ymin, zmin);
  if (!((xmax <= gXmax && xmax >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmax, ymin, zmax);
  if (!((xmin <= gXmax && xmin >= gXmin) && (ymin <= gYmax && ymin >= gYmin) && (zmax <= gZmax && zmax >= gZmin)))
    listOfPoints.emplace_back(xmin, ymin, zmax);

  // group box inside input box
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmax != xmax || gYmax != ymax || gZmax != zmax))
    listOfPoints.emplace_back(gXmax, gYmax, gZmax);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmin != xmin || gYmax != ymax || gZmax != zmax))
    listOfPoints.emplace_back(gXmin, gYmax, gZmax);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmin != xmin || gYmax != ymax || gZmin != zmin))
    listOfPoints.emplace_back(gXmin, gYmax, gZmin);
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmax <= ymax && gYmax >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmax != xmax || gYmax != ymax || gZmin != zmin))
    listOfPoints.emplace_back(gXmax, gYmax, gZmin);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmin != xmin || gYmin != ymin || gZmin != zmin))
    listOfPoints.emplace_back(gXmin, gYmin, gZmin);
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmin <= zmax && gZmin >= zmin)) &&
      (gXmax != xmax || gYmin != ymin || gZmin != zmin))
    listOfPoints.emplace_back(gXmax, gYmin, gZmin);
  if (((gXmax <= xmax && gXmax >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmax != xmax || gYmin != ymin || gZmax != zmax))
    listOfPoints.emplace_back(gXmax, gYmin, gZmax);
  if (((gXmin <= xmax && gXmin >= xmin) && (gYmin <= ymax && gYmin >= ymin) && (gZmax <= zmax && gZmax >= zmin)) &&
      (gXmin != xmin || gYmin != ymin || gZmax != zmax))
    listOfPoints.emplace_back(gXmin, gYmin, gZmax);

  if (!listOfPoints.empty()) {
    xmin = ymin = zmin = DBL_MAX;
    xmax = ymax = zmax = -DBL_MAX;
    for (std::vector<V3D>::const_iterator it = listOfPoints.begin(); it != listOfPoints.end(); ++it) {
      //			std::cout<<(*it)<<'\n';
      if ((*it)[0] < xmin)
        xmin = (*it)[0];
      if ((*it)[1] < ymin)
        ymin = (*it)[1];
      if ((*it)[2] < zmin)
        zmin = (*it)[2];
      if ((*it)[0] > xmax)
        xmax = (*it)[0];
      if ((*it)[1] > ymax)
        ymax = (*it)[1];
      if ((*it)[2] > zmax)
        zmax = (*it)[2];
    }
  }
}

#ifdef ENABLE_OPENCASCADE
TopoDS_Shape CompGrp::analyze() {
  TopoDS_Shape Result = A->analyze();
  Result.Complement();
  return Result;
}

TopoDS_Shape BoolValue::analyze() { return TopoDS_Shape(); }
#endif
} // namespace Mantid::Geometry
