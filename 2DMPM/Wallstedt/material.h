// Philip Wallstedt 2004-2009
#ifndef PATCH_H
#define PATCH_H

#include <map>
#include <vector>
#include <list>
#include "tensor.h"
#include "patch.h"
using namespace std;

// Rather than re-construct the shape and gradient weights within 
// each interpolation function we find them once per time step and store
// them in partContribs.
struct partContribs{
   struct portion{
      int idx; // node index
      double weight,gradx,grady; // weights
      void setPortion(int a,double b,double c,double d){idx=a;weight=b;gradx=c;grady=d;}
   };
   int Npor; // how many nodes this particle contributes to
   vector<portion>portionArray; // contribution list size set in post-setup configuration
   const portion&operator[](int l)const{return portionArray[l];}
};

/////// managed Array /////////////////////////
class material;

// All arrays used in the patch are managed by this class.
// An array is declared once (in the patch) and after that its
// size is adjusted together with all the other arrays of
// its group.  This makes it easy to add new variables to
// a patch, which we do all the time.

// The array manager creates and manages groups of arrays.  All arrays in a group
// have the same size, but do not contain the same type of data.

// These are vectors of pointers to managed arrays.  Every time
// a new array is declared, its address is added to one of these
// vectors as part of the array's constructor.
// These are only used during array construction.  Later on, the
// patch constructor will make copies of these, then clear the
// originals.
static vector<arraySC*>commonPartsPtr;
static vector<arraySC*>commonNodesPtr;


// These specialize each group of managed arrays and define the 
// constructor which adds each array's address to the managed group.
template<typename T>class partArray:public managedArray<T>{friend class patch;partArray(){commonPartsPtr.push_back(this);}};
template<typename T>class nodeArray:public managedArray<T>{friend class patch;nodeArray(){commonNodesPtr.push_back(this);}};

/////// material class ///////////////////////////

// A patch holds node and particle data for a particular spatial region of
// the problem.  Parallel and/or multi-resolution versions of GIMP
// need multiple patches that communicate with each other.
// Ultimately, patches may need to be created and destroyed as the
// solution progresses, though this capability is not currently used.
class material{
// const patch-wide values
public:
private:
   vector<arraySC*>parts;           // the managed group of particles
   int numberOfParts;               // Nparts
public:
   double Ymod;                     // Young's Modulus
   double dens;                     // Density
   double pois;                     // Poisson's ratio
   double damp;                     // Damping coefficient
   double load;                     // Some common load imposed on the problem
   double vwav;                     // The wave speed in the material = sqrt(Ymod/dens)
   double ppe;                      // Particles per edge.  2 ppe = 4 particles in a cell
// particle variables
   partArray<partContribs>    pCon;          // nodes to which this particle contributes, and weights for each
   partArray<Matrix2>         pF;            // deformation gradient
   partArray<Matrix2>         pGv;           // velocity gradient
   partArray<Matrix2>         pVS;           // volume*stress
   partArray<Vector2>         pX;            // initial position
   partArray<Vector2>         px;            // position
   partArray<Vector2>         pv;            // velocity
   partArray<Vector2>         pfe;           // external force
   partArray<Vector2>         pw;            // momentum
   partArray<Vector2>         pvI;           // velocity increment
   partArray<Vector2>         pxI;           // position increment
   partArray<double>          pm;            // mass
   partArray<double>          pV;            // Initial volume
   partArray<double>          pJ;            // Jacobian
// declared methods
   int copyPart(const patch&pch,int i);
   void resizeParts(int);
   void resizeNodes(int);
   patch(const int Nx,
         const int Ny,
         const double bx,
         const double by,
         const double ex,
         const double ey,
         const int Ng,
         const double th);
   patch&operator=(const patch&); // never defined, never used
public:
// inline methods
   int Npart()const{return numberOfParts;}
   int Nnode()const{return numberOfNodes;}
   int insertNode(const Vector2&v){
      resizeNodes(Nnode()+1);
      gx[Nnode()-1]=v;
      return Nnode()-1;
   }
   int appendPart(){
      int nidx;
      resizeParts(Npart()+1);
      nidx=Npart()-1;
      return nidx;
   }
   // These are all inlined and optimized by the compiler (I hope)
   // These functions find the cell within which a particle is contained
   // The lower left node of the cell is always defined as the "parent"
   // node of that cell.
   double rsx(const double x)const{return(x-regionBegin.x)/dx+Nghost;}
   double rsy(const double y)const{return(y-regionBegin.y)/dy+Nghost;}
   int rsi(const double x)const{return int(floor(rsx(x)));}
   int rsj(const double y)const{return int(floor(rsy(y)));}
   int inCell(const double x,const double y)const{ // find the cell into which x,y falls
      const int i=rsi(x);
      const int j=rsj(y);
      if(i>=0&&i<I&&j>=0&&j<J)return j*I+i;
      else throw"inCell:outside findable region";
   }
   int inCell(const int     p)const{return inCell(px[p].x,px[p].y);}
   int inCell(const Vector2&p)const{return inCell(    p.x,    p.y);}
   bool inRegion(const double x,const double y){return x>=regionBegin.x&&x<regionEnd.x&&y>=regionBegin.y&&y<regionEnd.y;}
   bool inRegion(const Vector2&p){return inRegion(p.x,p.y);}
   // This is a variation on finding the cell.  However, this is for
   // GIMP particles where portions of a particle may fall in up to
   // four cells.  Hence, we return the lowest and left-est node to
   // which the particle may contribute.
   int inCell9(const Vector2&p)const{
      const double&x=p.x;
      const double&y=p.y;
      const int i=rsi(x);
      const int j=rsj(y);
      if(!(i>=0&&i<I&&j>=0&&j<J))throw"inCell9:outside findable region";
      const double enx=rsx(x)-double(rsi(x));
      const double eny=rsy(y)-double(rsj(y));
      const int io=(enx<.5?i-1:i);
      const int jo=(eny<.5?j-1:j);
      return jo*I+io;
   }
};

template<typename T>
std::ostream& operator<<(std::ostream&os,managedArray<T>&a){
   for(int i=0;i<a.size();i++)os<<i<<'['<<a[i]<<"]\n";
   return os;
}

#endif







