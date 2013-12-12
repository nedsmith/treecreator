#include <windows.h>
#include <commctrl.h>
#include <fstream>
#include <tchar.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>

#include "resource.h"
#include "StdAfx.h"

#pragma comment(linker,"/manifestdependency:\"type='win32' name='Microsoft.Windows.Common-Controls' version='6.0.0.0' processorArchitecture='*' publicKeyToken='6595b64144ccf1df' language='*'\"")

#define TWOPI 6.2831
#define MINSTRETCH 0.5
#define MAXBRANCHES 25
#define MAXDIVS 50
#define HALFROOT3 0.86603
#define UPDIRECTION (vector(0,1,0))


using namespace std;

class vector;
class vertex;
class face;
class growPoint;
class Camera;
class FixedCamera;


static TCHAR szWindowClass[] = _T("win32app");
static TCHAR szTitle[] = _T("Tree Maker");

HINSTANCE hInst;
HWND hWnd1;


HWND levelSpin;
HWND angle1Spin;
HWND angle2Spin;
HWND heightSpin;
HWND divSpin;
HWND sizeSpin;
HWND thickSpin;
HWND bScaleSpin;
HWND randSpin;



// Define 3-vectors
class vector{
public:
double x,y,z;
vector () {x=0;y=0;z=0;};
vector(double xcoord, double ycoord, double zcoord) {x=xcoord; y=ycoord; z=zcoord;};
vector operator+(vector v2);
vector operator-(vector v2);
vector operator=(vector v2);
vector operator*(vector v2);
vector operator*(double scalar);
vector unit();
void makeUnit();
double mod();
};

//vector addition
vector vector::operator+(vector v2) {
	vector sum;
	sum.x = x+v2.x;
	sum.y = y+v2.y;
	sum.z = z+v2.z;
	return sum;
};

//vector subtraction
vector vector::operator-(vector v2) {
	vector sum;
	sum.x = x-v2.x;
	sum.y = y-v2.y;
	sum.z = z-v2.z;
	return sum;
};


//vector assignment
vector vector::operator=(vector v2) {
	x=v2.x;
	y=v2.y;
	z=v2.z;
	return *this;
};


// scalar multiplication
vector vector::operator*(double scalar) {
	vector scaled;
	scaled.x = x*scalar;
	scaled.y = y*scalar;
	scaled.z = z*scalar;
	return scaled;
}


// Cross product
vector vector::operator*(vector v2) {
vector cross;
cross.x = y*v2.z - z*v2.y;
cross.y = z*v2.x - x*v2.z;
cross.z = x*v2.y - y*v2.x;
return cross;
};


//scalar product (maybe dodgy - not a member of vector class)
double dot(vector v1,vector v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}



// returns the absolute magnitude of a vector
double vector::mod() {
	return sqrt(x*x+y*y+z*z);
};

// returns a unit version of the vector
vector vector::unit() {
	double modulus;
	modulus=this->mod();
	if (modulus<1e-5f) {return vector(1,0,0); }
	else {return (*this)*(1.0/modulus); }
};

// makes a vector a unit;
void vector::makeUnit() {
	*this=this->unit();
};


vector world1, world2, world3;	//global coordinates


// Vertex and Face Classes--------------------------------------------------------------------------------

//store vertex coordinates as:
typedef short int VXCOMPONENT; 
#define VXSCALE 100.0

// Polygon vertex
class vertex{
private:
	~vertex() {};
public:
	VXCOMPONENT comp1,comp2,comp3;
	short int id;
	vertex *next;
	vertex(double xcoord,double ycoord,double zcoord);
	vertex(vector vect);	
	vector pos(void);
	void setPos(vector vect);
	static int numVerts;
	static vertex *firstvert;
	static vertex *lastvert;
	static void updateIDs();
	static void delAll();
};

int vertex::numVerts=0;
vertex *vertex::firstvert=NULL;
vertex *vertex::lastvert=NULL;


// Polygon face
class face {
private:
	~face() {};
public:
	vertex	*a,*b,*c,*d;	//corners of the face in cyclic order
	bool quad;				//is it a quadrilateral? otherwise it's a triangle
	short int id;
	face *next;
	face(vertex *corner1, vertex *corner2, vertex *corner3);
	face(vertex *corner1, vertex *corner2, vertex *corner3, vertex *corner4);
	static int numFaces;
	static face *firstface;
	static face *lastface;
	static void delAll();
};

int face::numFaces=0;
face *face::firstface=NULL;
face *face::lastface=NULL;


// Vertex constructor - automatically makes a linked list
vertex::vertex(double xcoord,double ycoord,double zcoord) {
comp1 = (VXCOMPONENT) (xcoord * VXSCALE);
comp2 = (VXCOMPONENT) (ycoord * VXSCALE);
comp3 = (VXCOMPONENT) (zcoord * VXSCALE);
if (numVerts==0) {
	firstvert = this;
}
else {lastvert->next = this;}
lastvert = this;
next = NULL;
id=numVerts+1;
numVerts++;
};


vertex::vertex(vector vect) {
comp1 = (VXCOMPONENT) (vect.x * VXSCALE);
comp2 = (VXCOMPONENT) (vect.y * VXSCALE);
comp3 = (VXCOMPONENT) (vect.z * VXSCALE);
if (numVerts==0) {
	firstvert = this;
}
else {lastvert->next = this;}
lastvert = this;
next = NULL;
id=numVerts+1;
numVerts++;
};

// deletes all verticies
void vertex::delAll() {
	vertex *v1, *v2;
	v1=firstvert;
	while (v1!=NULL) {
		v2=v1;
		v1=v1->next;
		delete v2;
	}
	numVerts=0;
	firstvert=NULL;
	lastvert=NULL;
}


vector vertex::pos(void) {
	return vector(((double)comp1)/VXSCALE, ((double)comp2)/VXSCALE, ((double)comp3)/VXSCALE);
}


void vertex::setPos(vector vect) {
	comp1 = (VXCOMPONENT) (vect.x * VXSCALE);
	comp2 = (VXCOMPONENT) (vect.y * VXSCALE);
	comp3 = (VXCOMPONENT) (vect.z * VXSCALE);
}


// face constructor
face::face(vertex *corner1, vertex *corner2, vertex *corner3) {
a=corner1;
b=corner2;
c=corner3;
d=NULL;
quad=0;
if (numFaces==0) {
	firstface = this;}
else {lastface->next = this;}
lastface = this;
next = NULL;
id=numFaces+1;
numFaces++;
};

// and the quad version
face::face(vertex *corner1, vertex *corner2, vertex *corner3, vertex *corner4) {
a=corner1;
b=corner2;
c=corner3;
d=corner4;
quad=1;
if (numFaces==0) {
	firstface = this;}
else {lastface->next = this;}
lastface = this;
next = NULL;
id=numFaces+1;
numFaces++;
};


// deletes all faces
void face::delAll() {
	face *f1, *f2;
	f1=firstface;
	while (f1!=NULL) {
		f2=f1;
		f1=f1->next;
		delete f2;
	}
	numFaces=0;
	firstface=NULL;
	lastface=NULL;
}


// returns a random float between 5 and -5
double random() {
	return  (double)((rand()%1000)-500)*0.01;
};

// a random float between 0 and 1
double rand1() {
	return (double) (rand()%10000)*0.001;
};


// returns a random vector;
vector randVector() {
	return vector(random(),random(),random());
};



// GrowPoint Class ------------------------------------------------------------------------------------
// a point that a branch may grow from
class growPoint {
private:
	vector centrePoint;		// central point to grow from
	double radius;			// approximate radius

	vector growAxis;	// direction to grow along
	vector axis2;		// points from centrePoint to startVert
	vector axis3;		// these should be an orthonomal triple
	
	vector steady;		//   a unit vector perpendicular to growAxis
						// recording the direction of the last fork

	bool clockwise;

public:
	vertex *startVert;	//pointer to the first of a sequence of consecutive verts in the linked-list
	int divisions;		//where the number of such verts is divisions

	// constructors and destructor
	growPoint(vector middle, int divs, double thickness);
	growPoint(vector middle,vector ax1, vector ax2, vector ax3, int divs, vertex *v1, bool clocks, double thickness, vector axSteady);
	growPoint(vector middle,vector ax1, vector ax2, vector ax3, int divs, vertex *v1, bool clocks, double thickness);
	~growPoint();

	void grow(double length, vector newGrowAxis, double thickness);
	void grow(double length);
	void growRel(double relLength, double turnAngle, double bendAngle, double scale);
	void spike(double length);
	growPoint *fork(double length, vector branch1, vector branch2, int divs1, int divs2, double thickness);
	growPoint *forkRel(double relLength, double turnAngle, double branchAngle1, double branchAngle2, double scale);
	growPoint *forkRel(double relLength, double turnAngle, double branchAngle, double scale);
	void checkRadius();
	static void finish(growPoint * gp, double length);
	growPoint * invert();
	void randomise();
	void smooth();


};


//makes a growpoint with new verts in a circle facing upwards
growPoint::growPoint(vector middle, int divs, double thickness) {
	int i;
	double theta;
	vertex *v1;

	centrePoint=middle;
	growAxis=UPDIRECTION;
	axis2=UPDIRECTION*vector(1,0,0);
	if (axis2.mod()<1e-5) axis2=UPDIRECTION*vector(0,1,0);
	axis2.makeUnit();
	axis3=growAxis*axis2;
	steady=axis2;
	divisions=divs;
	clockwise=true;
	radius=thickness;

	v1 = new vertex(middle+axis2*thickness);
	for (i=1;i<divs;i++) {
		theta=TWOPI*(double)i/(double)divs;
		new vertex(middle+axis2*(thickness*cos(theta))+axis3*(thickness*sin(theta)));
	}

	startVert=v1;
};

//makes a growpoint with existing verts
growPoint::growPoint(vector middle,vector ax1, vector ax2, vector ax3, int divs, vertex *v1, bool clocks, double thickness, vector axSteady) {
	centrePoint=middle;
	growAxis=ax1;
	axis2=ax2;
	axis3=ax3;
	steady=axSteady;
	divisions=divs;
	startVert=v1;
	clockwise=clocks;
	radius=thickness;
};


//old version without steady vector
growPoint::growPoint(vector middle,vector ax1, vector ax2, vector ax3, int divs, vertex *v1, bool clocks, double thickness) {
	centrePoint=middle;
	growAxis=ax1;
	axis2=ax2;
	axis3=ax3;
	steady=ax2;
	divisions=divs;
	startVert=v1;
	clockwise=clocks;
	radius=thickness;
};


//destructor
growPoint::~growPoint() {
// (empty)
};

// grows along the predefined axis and turns to face a new direction
void growPoint::grow(double length, vector newGrowAxis, double thickness) {
	vertex *loVert1, *loVert2, *hiVert1, *hiVert2, *newStartVert;

	vector newMiddle,newAx2,newAx3,newSteady;
	vector norm, proj2, proj3;
	vector curveDir, noCurveDir; //directions along and perpendicular to the curve of the branch
	double stretch;
	double co1, co2, theta, theta2, direction;
	int i;


// calculate plane of join and new axes

	growAxis.makeUnit();
	newGrowAxis.makeUnit();
	norm=(growAxis+newGrowAxis).unit();
	noCurveDir=growAxis*newGrowAxis;
	if (noCurveDir.mod()<1e-5) {
		noCurveDir=axis2;
		curveDir=axis3;
	} else {
		noCurveDir.makeUnit();
		curveDir=(norm*noCurveDir).unit();
	}

	proj2=(norm*axis2)*norm;
	if (proj2.mod()<1e-5) {
		proj2=noCurveDir;
		proj3=curveDir;
	} else {
		proj2.makeUnit();
		proj3=(norm*proj2).unit();
	}




	newMiddle = centrePoint+growAxis*length;
	stretch=dot(growAxis,norm);
	if (stretch < MINSTRETCH) stretch=MINSTRETCH;

co1 = dot(proj2,noCurveDir);
co2 = dot(proj2,curveDir);


if (fabs(co1)<1e-5) {
	if (co2>0) {theta=TWOPI*0.25;}
	else {theta=-TWOPI*0.25f;}
} else {
	theta=atan(co2/co1);
	if (co1<0.0f) theta+=TWOPI*0.5;
}

co2=co2*stretch;


// draw new verticies and faces of updated growPoint

if (clockwise) {direction=1;}
else {direction=-1;}

newStartVert = new vertex(newMiddle+noCurveDir*(cos(theta)*thickness)+curveDir*(sin(theta)*thickness/stretch));

hiVert1 = newStartVert;
loVert1 = startVert;


for (i=1;i<divisions;i++) {
theta2=theta+direction*((double)i/(double)divisions)*TWOPI;

hiVert2=hiVert1;
loVert2=loVert1;
hiVert1 = new vertex(newMiddle+noCurveDir*(cos(theta2)*thickness)+curveDir*(sin(theta2)*thickness/stretch));
loVert1 = loVert1->next;

new face(loVert1,hiVert1,hiVert2,loVert2);
}

new face(startVert,newStartVert,hiVert1,loVert1);



	newAx2=(newGrowAxis*proj2)*newGrowAxis;
	if (newAx2.mod()<1e-5f) {
		newAx2=axis2;
		newAx3=axis3;
	} else {
		newAx2.makeUnit();
		newAx3=(newGrowAxis*newAx2).unit();
	}

	newSteady=(newGrowAxis*steady)*newGrowAxis;
	if (newSteady.mod()<1e-5) {
		newSteady=newAx2;
	} else {
		newSteady.makeUnit();
	}

centrePoint=newMiddle;
growAxis=newGrowAxis;
axis2=newAx2;
axis3=newAx3;
steady=newSteady;
startVert=newStartVert;
radius=thickness;

}


// overloaded version of grow that keeps existing direction and radius
void growPoint::grow(double length) {
	grow(length, growAxis, radius);
}



// alternative version of grow specifying two angles (works relative to steady vector)
void growPoint::growRel(double relLength, double turnAngle, double bendAngle, double scale) {
	vector steadyPerp = growAxis*steady;
	grow(radius*relLength,growAxis*cos(bendAngle)+steady*(sin(bendAngle)*cos(turnAngle))+steadyPerp*(sin(bendAngle)*sin(turnAngle)),radius*scale);
}



void growPoint::spike(double length) {
	vertex *apex, *v1, *v2;
	int i;

	apex = new vertex(centrePoint+growAxis*length*radius);
	v2=startVert;
	v1=startVert;
	for (i=1;i<divisions;i++) {
		v1=v1->next;
		new face(v1,v2,apex);
		v2=v1;		
	}
	new face(startVert,v2,apex);

}


void growPoint::finish(growPoint *gp, double length) {
	gp->spike(length);
	delete gp;
}


// growPoint::fork
// my forking routine - this is a big function
growPoint *growPoint::fork(double length, vector branch1, vector branch2, int divs1, int divs2, double thickness) {

vector norm1, norm2, norm3, middle, centre1, arcCentre, out, ax1, ax2, ax3, proj, branch1ax2;
vector turnPoint1, turnPoint2, forkPoint;
double share1, share2, shareFork, totDivs, minDivs;
bool evens;
double separation, radius1, radius2, radius3, dist1, dist2, dist3, cosec1, cosec2, cosec3, oldRadius;
double stretch1, stretch2, stretch3;
double theta, thetaInc, co1, co2, direction;
int shift;
bool forward;
int i, j;
vertex *vTop, *vBottom, *vEnd, *startBranch1, *startBranch2, *vLoEnd;
vertex *hiVert1, *loVert1, *hiVert2, *loVert2;
vertex *backSegment[MAXDIVS];
int backDivs;


// tidy up input
if (divs1<3) divs1=3;
if (divs2<3) divs2=3;
if (divs1+divs2 < divisions) divs1=divisions-divs2;
if (divs1+divisions < divs2) divs2=divisions+divs1;
if (divisions+divs2 < divs1) divs1=divisions+divs2;


growAxis.makeUnit();
branch1.makeUnit();
branch2.makeUnit();

if (clockwise) {direction=1;}
else {direction=-1;}

// Find normals to the 3 planes of intersection
norm1 = (growAxis+branch1);
norm2 = (growAxis+branch2);
norm3 = (branch2-branch1);
out=((norm2*norm1)*direction).unit();

ax1=(out*norm1).unit();
ax2=(norm2*out).unit();
ax3=(out*norm3).unit();

stretch1=dot(norm1.unit(),branch1);
stretch2=dot(norm2.unit(),branch2);
stretch3=dot(norm3.unit(),branch2);  //NB yes really

if (stretch1<MINSTRETCH) stretch1=MINSTRETCH;
if (stretch2<MINSTRETCH) stretch2=MINSTRETCH;
if (stretch3<MINSTRETCH) stretch3=MINSTRETCH;

// work out no. of divisions in each bounding arc
totDivs = (double)(divisions+divs1+divs2);
if ((int)totDivs%2==0) {evens=true;}
else {evens=false;}
share1=totDivs/2-(double)divs2;
share2=totDivs/2-(double)divs1;
shareFork=totDivs/2-(double)divisions;

minDivs=shareFork;
if (share1<minDivs) minDivs=share1;
if (share2<minDivs) minDivs=share2;


//determine overall size
separation = direction*thickness*minDivs/(double)divisions;

// calculate sizes of arcs


if (minDivs>0.25) {

	cosec1=share1/minDivs;
	cosec2=share2/minDivs;
	cosec3=shareFork/minDivs;

	radius1=separation*cosec1;
	radius2=separation*cosec2;
	radius3=separation*cosec3;

	dist1=radius1*sqrt(1-pow(cosec1,-2))/stretch1;
	dist2=radius2*sqrt(1-pow(cosec2,-2))/stretch2;
	dist3=radius3*sqrt(1-pow(cosec3,-2))/stretch3;

} else { // bug fix: for case when divs1+divs2=divisions

	// dummy values - don't use
	cosec1=999;
	cosec2=999;
	cosec3=999;

	radius1=direction*thickness*share1/(double)divisions;
	radius2=direction*thickness*share2/(double)divisions;
	radius3=direction*thickness*shareFork/(double)divisions;

	dist1=radius1/stretch1;
	dist2=radius2/stretch2;
	dist3=radius3/stretch3;

}


turnPoint1=ax1*((dist1+radius1)/stretch1);
turnPoint2=ax2*((dist2+radius2)/stretch2);
forkPoint =ax3*((dist3+radius3)/stretch3);

middle = centrePoint + growAxis*length - (turnPoint1+turnPoint2)*0.5;

//Now draw the verticies of the arcs

arcCentre=middle+ax1*dist1;
(minDivs<0.25) ? theta = 0 : theta=direction*asin(1.0f/cosec1);
thetaInc=direction*(TWOPI-2.0*direction*theta)/share1;

theta+=thetaInc;
startBranch1 = new vertex(arcCentre-ax1*(cos(theta)*radius1/stretch1)+out*(sin(theta)*radius1));
backSegment[0] = startBranch1;

for (i=1;i<=share1-1;i++) {	
	theta+=thetaInc;
	backSegment[i] = new vertex(arcCentre-ax1*(cos(theta)*radius1/stretch1)+out*(sin(theta)*radius1));
}
backDivs = i;

arcCentre=middle+ax3*dist3;
(minDivs<0.25) ? theta = 0 : theta=direction*asin(1.0f/cosec3);
thetaInc=direction*(TWOPI-2.0*direction*theta)/shareFork;

vBottom=vertex::lastvert;
if (evens) {
	startBranch2 = vertex::lastvert;	
} else {
	theta+=thetaInc/2;
	startBranch2 = new vertex(arcCentre-ax3*(cos(theta)*radius3/stretch3)-out*(sin(theta)*radius3));  //NB sign change: -out
}

for (i=1;i<=shareFork;i++) {	
	theta+=thetaInc;
	new vertex(arcCentre-ax3*(cos(theta)*radius3/stretch3)-out*(sin(theta)*radius3));
}

vTop = vertex::lastvert;

arcCentre=middle+ax2*dist2;
(minDivs<0.25) ? theta = 0 : theta=direction*asin(1.0f/cosec2);
thetaInc=direction*(TWOPI-2.0*direction*theta)/share2;

for (i=1;i<share2;i++) {	//NB not <=   (share1 etc take half-integer values)
	theta+=thetaInc;
	new vertex(arcCentre-ax2*(cos(theta)*radius2/stretch2)+out*(sin(theta)*radius2));
}

vEnd=vertex::lastvert;

//determine angle between fork and original branch axis

proj=(growAxis*out)*growAxis;
co1 = dot(proj,axis2);
co2 = dot(proj,axis3);

if (fabs(co1)<1e-5) {
	if (co2>0) {theta=TWOPI*0.25;}
	else {theta=-TWOPI*0.25;}
} else {
	theta=atan(co2/co1);
	if (co1<0.0f) theta+=TWOPI*0.5;
}


//join up arcs to original growPoint

loVert1=startVert;
loVert2=startVert;
shift=(int)(theta*(double)divisions/TWOPI);
for (i=1; i<divisions; i++) {
	loVert2=loVert2->next;
	if (i==shift) {loVert1=loVert2;}
}
vLoEnd=loVert2;
hiVert1=vTop;
forward=true;


for (i=0; i<divisions; i++) {

hiVert2=hiVert1;
loVert2=loVert1;

if (hiVert1==vEnd) {
	hiVert1=vBottom;
	forward=false;
	j=backDivs-2;
}
else if (hiVert1==startBranch1) {
	hiVert1=vTop;
	forward=true;
}
else if (forward) {hiVert1=hiVert1->next;}
else {hiVert1=backSegment[j]; j--;}

if (loVert1==vLoEnd) {loVert1=startVert;}
else {loVert1 = loVert1->next;}

new face(loVert1,loVert2, hiVert2, hiVert1);

}


//add extra triangle if divs add up to an odd number
if (evens==false) new face(vBottom,vEnd,startBranch2);



oldRadius=radius;

//guide existing growpoint along path of branch 2
startVert=startBranch2;
divisions=divs2;
growAxis=branch2;
radius=thickness*(double)divs2/(double)divisions;
centrePoint=middle+(turnPoint2+forkPoint)*0.5;
axis2=((branch2*((startBranch2->pos())-centrePoint))*branch2).unit();
axis3=(branch2*axis2).unit();
steady=((branch2*out)*branch2).unit();


// return a new growpoint along branch 1
centre1=middle+(turnPoint1+forkPoint)*0.5;
branch1ax2=((branch1*((startBranch1->pos())-centre1))*branch1).unit();
return new growPoint(centre1, branch1, branch1ax2, (branch1*branch1ax2).unit(), divs1, startBranch1, !clockwise, thickness*(double)divs1/(double)divisions, ((branch1*out)*branch1).unit());

} //-------------------------------------End growPoint::fork------------------------------------------------------------------------




//simpler version of fork taking three angles and returning two equal branches
growPoint *growPoint::forkRel(double relLength, double turnAngle, double branchAngle1, double branchAngle2, double scale) {

	vector branch1,branch2;
	vector steadyPerp;
	double co1, co2, co3;
	int newDivs;

	steadyPerp=growAxis*steady;

	co1=cos(branchAngle1);
	co2=sin(branchAngle1)*cos(turnAngle);
	co3=sin(branchAngle1)*sin(turnAngle);

	branch1=growAxis*co1+steady*co2+steadyPerp*co3;

	co1=cos(branchAngle2);
	co2=sin(branchAngle2)*cos(turnAngle);
	co3=sin(branchAngle2)*sin(turnAngle);

	branch2=growAxis*co1-steady*co2-steadyPerp*co3;

	if (scale<0.5) {newDivs=(divisions/2)+1;}
	else {newDivs = (int)((double)divisions * scale ) + 0;}
	if (newDivs<3) newDivs=3;

	return fork(relLength*radius, branch1, branch2, newDivs, newDivs, radius*scale*(double)divisions/(double)newDivs);

}

// older (symmetrical) version of forkRel
growPoint *growPoint::forkRel(double relLength, double turnAngle, double branchAngle, double scale) {
return forkRel(relLength, turnAngle, branchAngle, branchAngle, scale);
}



// updates the stored radius to the actual distance from the centre to the first vertex
// NB never do this
void growPoint::checkRadius() {
	radius=(centrePoint-(startVert->pos())).mod();
}


// creates a new grow point using the same verticies but facing the opposite direction
// returns pointer to the new grow point
growPoint * growPoint::invert() {
    growPoint * gp = new growPoint(centrePoint, growAxis*(-1.0), axis2, axis3*(-1.0), divisions, startVert, !clockwise, radius, steady);
    return gp;
}


void growPoint::randomise() {

double amount = 0.1;
vertex * v = startVert;
for (int i=0; i<divisions; i++) {
	v->setPos(v->pos() +  (v->pos()-centrePoint)*(rand1()*amount));
	v=v->next;
}
}


void growPoint::smooth() {
	
	vector oldPoints[MAXDIVS];
	vector newPoints[MAXDIVS];
	vertex * v;



	// store the old positions of verts along the ring
	v=startVert;
	for (int i=0; i<divisions; i++) {
		oldPoints[i]=v->pos();
		v=v->next;
	}


	for (int i=0; i<divisions; i++) {
		



	}

	// update the verts to new positions
	v=startVert;
	for (int i=0; i<divisions; i++) {
		v->setPos(newPoints[i]);
		v=v->next;
	}





}


// Global function
void SetWorldAxes(void) {
	world1=UPDIRECTION;
	world1.makeUnit();
	world2=world1*vector(1,0,0);
	if (world2.mod()<1e-5) world2=world1*vector(0,1,0);
	world2.makeUnit();
	world3=world1*world2;
}



// makes a cube with centre (midx,midy,midz) and side length edge
void MakeCube(double midx, double midy, double midz, double edge) {
double h=edge/(double)2.0;
vertex *v1,*v2,*v3,*v4,*v5,*v6,*v7,*v8;
v1 = new vertex(midx+h,midy+h,midz+h);
v2 = new vertex(midx+h,midy+h,midz-h);
v3 = new vertex(midx+h,midy-h,midz+h);
v4 = new vertex(midx+h,midy-h,midz-h);
v5 = new vertex(midx-h,midy+h,midz+h);
v6 = new vertex(midx-h,midy+h,midz-h);
v7 = new vertex(midx-h,midy-h,midz+h);
v8 = new vertex(midx-h,midy-h,midz-h);

new face(v1, v2, v4, v3);
new face(v5, v6, v8, v7);
new face(v1, v2, v6, v5);
new face(v3, v4, v8, v7);
new face(v1, v3, v7, v5);
new face(v2, v4, v8, v6);

};





// makes a triangular prism without endcaps
void MakePrism(vector startpos, vector length, double thickness) {
vector i1,endpos,perp,perp2;
double perpmod;
vertex *v1, *v2, *v3, *v4, *v5, *v6;

if (length.mod()<1e-5) length=vector(0,0,0.01f);  // if length is 0, make it up

endpos=startpos+length;
i1=vector(1,0,0);		//i1 is a coordinate unit vector
perp=i1*length;			//try and make perp perpendicular to length using cross products
perpmod=perp.mod();
if (perpmod<1e-5) {		//if perp is 0
	i1=vector(0,1,0);	//try again
	perp=i1*length;		//this always works
	perpmod=perp.mod();
}

perp2=length*perp;		//perp and perp2 are coordinate vectors in the endcap-plane

perp=perp*(thickness/perpmod);		//stretch perp and perp2 to the right length
perp2=perp2*(thickness/perp2.mod());

//Make the prism
v1 = new vertex(startpos+perp);
v2 = new vertex(startpos-(perp*0.5)+(perp2*HALFROOT3));
v3 = new vertex(startpos-(perp*0.5)-(perp2*HALFROOT3));
v4 = new vertex(endpos+perp);
v5 = new vertex(endpos-(perp*0.5)+(perp2*HALFROOT3));
v6 = new vertex(endpos-(perp*0.5)-(perp2*HALFROOT3));

new face(v1,v2,v5,v4);
new face(v2,v3,v6,v5);
new face(v3,v1,v4,v6);

};





// tree attribute list - all properties for drawing a tree
struct treeAttributes {
	double scale;
	double height;
	double segmentLength;
	int iterations;
	int radDivisions;
	double branchAng1;
	double branchAng2;
	double branchScale;
};







// creates a tree
void makeTree(vector position, treeAttributes * pTreeAtts) {

	growPoint *gp[MAXBRANCHES];
	int level[MAXBRANCHES];
	int i;
	int lastBranch=1;
	int minLevels=4;


double trunkLength = pTreeAtts->height;
double size = pTreeAtts->segmentLength;
double thickness = pTreeAtts->scale;
int maxLevels = pTreeAtts->iterations;
int trunkDivs = pTreeAtts->radDivisions;
double randomness = 0;



	for (i=0; i<MAXBRANCHES; i++) {
		gp[i]=NULL;
		level[i]=0;
	}

	gp[0] = new growPoint(position,trunkDivs,thickness);
	level[0]=0;
	gp[1] = gp[0]->forkRel(trunkLength,0,pTreeAtts->branchAng1, pTreeAtts->branchAng2, 0.6);
	level[1]=0;

	do {
			if ( level[lastBranch]<minLevels   || 
					(random()>-3   &&   level[lastBranch]<maxLevels)    ) {
						gp[lastBranch + 1] = gp[lastBranch]->forkRel(size*(10+randomness*random())*0.1,  0  ,pTreeAtts->branchAng1+ TWOPI*((random()*randomness*0.01)), pTreeAtts->branchAng2,pTreeAtts->branchScale);
				level[lastBranch+1] = ++level[lastBranch];
				lastBranch++;
			} else {
				growPoint::finish(gp[lastBranch],size);
				lastBranch--;
			}
	} while (lastBranch>=0);



}




// Opens a Save File dialogue box, then
// creates a .obj file containing all verticies and faces
// requires handle to the parent window
// NB makes an error if some verts don't exist
void makeObjFile(HWND hWnd) {
vertex *v;
face *f;
ofstream objfile;
OPENFILENAME ofn;

TCHAR filename[MAX_PATH] = _T("Mesh");

    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = hWnd;
    ofn.lpstrFilter = _T("Wavefront mesh (*.obj)\0*.obj\0");
    ofn.lpstrFile = filename;
    ofn.nMaxFile = MAX_PATH;
    ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT;
    ofn.lpstrDefExt = _T("obj");

    if(GetSaveFileName(&ofn))
    {

		objfile.open(filename);
		objfile << "o 1\n\n";

		// write vertices
		for (v=vertex::firstvert; v!=NULL; v=v->next) {
			objfile << "v " << v->comp1 << " " << v->comp2 << " " << v->comp3 << "\n";
		}

		objfile << "\nusemtl Default\n\n";

		// write faces
		for (f=face::firstface; f!=NULL; f=f->next) {
			objfile << "f " << f->a->id << " " << f->b->id << " " << f->c->id;
			if (f->quad==1) {objfile << " " << f->d->id; }
			objfile << "\n";
		}

		objfile.close();
	}
} // ------------------- end makeObjFile



// updates the ID numbers of all verticies (unused)
void vertex::updateIDs() {
	int i=1;
	vertex *v;

	for (v=firstvert; v!=NULL; v=v->next) {
		v->id=i;
		i++;
	}
}


// deletes all faces and verticies
void clearAll() {
	vertex::delAll();
	face::delAll();
}



class Camera {
protected:
	vector camPosition;

	vector viewAxis;
	vector axis2;
	vector axis3;


public:
	double xViewField;
	double yViewField;

	void moveTo (vector pos) {camPosition=pos; };
	void pointAt(vector target);
	vector proj(vector target);
	Camera() {moveTo(vector(0,0,-10)); pointAt(vector(0,0,0)); xViewField=1; yViewField=0.75; }
	~Camera() {};

};


class FixedCamera : public Camera {
private:
	vector centre;

	double theta;
	double phi;
	double distance;
	double angleStep;
	double zoomStep;

public:

	void update(void);

	void rotateLeft(void)  {theta-=angleStep; update(); }
	void rotateRight(void) {theta+=angleStep; update(); }
	void rotateUp(void)    {phi+=angleStep;			if (phi>3.14159) phi=3.14159;		update(); }
	void rotateDown(void)  {phi-=angleStep;			if (phi<0) phi=0;					update(); }
	void zoomIn(void)      {distance/=zoomStep;		if (distance<1e-5) distance=1e-5;	update(); }
	void zoomOut(void)     {distance*=zoomStep;		if (distance>1e5) distance=1e5;		update(); }

	FixedCamera();
	~FixedCamera() {};

};


// updates a fixed camera's position and rotation to match the data in FixedCamera
void FixedCamera::update(void) {
	viewAxis=world1*cos(phi) + world2*(sin(phi)*cos(theta)) + world3*(sin(phi)*sin(theta));
	axis2= world2*sin(theta) - world3*cos(theta);
	axis3=viewAxis*axis2;
	camPosition=centre - viewAxis * distance;
}


// Fixed Camera constructor
FixedCamera::FixedCamera() {
	centre = vector(0,0,0);
	theta = 0;
	phi = TWOPI/4;
	distance = 20;
	angleStep = 0.01;
	zoomStep = 1.1;
	update();
}



// point at a target
void Camera::pointAt(vector target) {
	viewAxis=target-camPosition;
	viewAxis.makeUnit();
	axis2=UPDIRECTION*viewAxis;
	axis2.makeUnit();
	axis3=viewAxis*axis2;
}


// projects the target vector into the camera's view
vector Camera::proj(vector target) {
vector relVect, result;
double scale;

relVect=target-camPosition;
scale=dot(relVect, viewAxis);
result=vector(1e10,0,0); //default result for off-screen vectors

if (scale>1e-5) {		// if in front of camera
	result = vector(dot(relVect,axis2)/scale,dot(relVect,axis3)/scale, 0);
}

return result;

}



class Graphic {

private:
	HWND hDrawWindow;
	RECT drawArea;
	HDC hdcBuffer;
	HBITMAP hbmBuffer;
	HGDIOBJ hOldBuffer;
	int width;
	int height;

public:

	FixedCamera *cam;

	Graphic(HWND hWnd);
	~Graphic();
	void chooseDrawArea();
	void render();
	void redraw();
	int getX(vector v);
	int getY(vector v);

};


// rescales the drawing area to match the window
void Graphic::chooseDrawArea(void) {
	RECT rc;

	GetClientRect(hDrawWindow, &rc);
	rc.bottom-=0;
	rc.top+=0;
	rc.left+=200;
	rc.right-=0;
	drawArea=rc;

	width=rc.right-rc.left;
	height=rc.bottom-rc.top;

}


// constructor
Graphic::Graphic(HWND hWnd) {
	HDC hdc;

	hDrawWindow = hWnd;
	chooseDrawArea();

	cam = new FixedCamera();
	
	hdc = GetDC(hWnd);
	hdcBuffer = CreateCompatibleDC(hdc);
	hbmBuffer = CreateCompatibleBitmap(hdc, 1600, 1200);
	hOldBuffer = SelectObject(hdcBuffer, hbmBuffer);
	ReleaseDC(hWnd, hdc);

}

// destructor
Graphic::~Graphic() {
	delete cam;
	SelectObject(hdcBuffer, hOldBuffer);
	DeleteDC(hdcBuffer);
	DeleteObject(hbmBuffer);
}



int Graphic::getX(vector v) {
return (int)( (double)(width)/2  + ((cam->proj(v)).x) * (double)(width) / ((cam->xViewField)*2));
}


int Graphic::getY(vector v) {
return (int)( (double)(height)/2  + ((cam->proj(v)).y) * (double)(width) / ((cam->xViewField)*2));
}

// for use out of WM_PAINT
void Graphic::render() {

HDC hdc;
HPEN hPen=NULL;
HBRUSH hBlackBrush;
HGDIOBJ hObjOld;
face *f;
RECT rcBuf;

rcBuf.left=0;
rcBuf.top=0;
rcBuf.right=width;
rcBuf.bottom=height;


hBlackBrush = CreateSolidBrush(RGB(0,0,0));
FillRect(hdcBuffer, &rcBuf, hBlackBrush);
DeleteObject(hBlackBrush);

hPen=CreatePen(PS_SOLID, 0, RGB(0,255,0));
hObjOld=SelectObject(hdcBuffer, hPen);

for (f=face::firstface; f->next != NULL && f!=face::lastface; f=f->next) {
	
MoveToEx(hdcBuffer, getX(f->a->pos()), getY(f->a->pos()), NULL);
LineTo(hdcBuffer, getX(f->b->pos()), getY(f->b->pos()));
LineTo(hdcBuffer, getX(f->c->pos()), getY(f->c->pos()));

#ifdef FULLFACERENDER
if (f->quad) {
LineTo(hdcBuffer, getX(f->d->pos()), getY(f->d->pos()));
}
LineTo(hdcBuffer, getX(f->a->pos()), getY(f->a->pos()));
#endif

}

SelectObject(hdcBuffer, hObjOld);

DeleteObject(hPen);

hdc = GetDC(hDrawWindow);
BitBlt(hdc, drawArea.left, drawArea.top, width, height, hdcBuffer, 0,0, SRCCOPY);
ReleaseDC(hDrawWindow, hdc);

}


// for use in WM_PAINT - doesn't render to buffer
void Graphic::redraw() {
	HDC hdc;
	PAINTSTRUCT ps;

	hdc = BeginPaint(hDrawWindow, &ps);
	BitBlt(hdc, drawArea.left, drawArea.top, width, height, hdcBuffer, 0,0, SRCCOPY);
	EndPaint(hDrawWindow, &ps);
}










// add an up-down control with buddy window and caption
//	compile with ComCtl32.lib
HWND WINAPI CreateSpinner(
	  HWND hWndParent,
	  HINSTANCE hInst,
	  int xpos, int ypos,
	  int min, int max,
	  int initial,
	  TCHAR caption[]
	  )
{

	      HWND hWndCaption = CreateWindowEx(
			  0,
			  WC_STATIC,
			  caption,
			  WS_CHILD | WS_VISIBLE | WS_OVERLAPPED | SS_LEFT,
			  xpos, ypos,
			  120, 21,
			  hWndParent,
			  0,          //identifier - change this
			  hInst,
			  NULL
		  );


		  HWND hWndBuddy = CreateWindowEx(
			  0,
			  WC_EDIT,
			  _T("7"),
			  WS_CHILD | WS_VISIBLE | WS_OVERLAPPED | ES_NUMBER | ES_RIGHT | WS_BORDER,
			  xpos+120, ypos,
			  40, 21,
			  hWndParent,
			  0,          //identifier - change this
			  hInst,
			  NULL
		  );

		  HWND hWndSpin = CreateWindowEx(
			  0,
			  UPDOWN_CLASS,
			  _T(""),
			  WS_CHILD | WS_VISIBLE | WS_OVERLAPPED | UDS_AUTOBUDDY | UDS_SETBUDDYINT,
			  xpos + 160, ypos,
			  20, 21,
			  hWndParent,
			  0,          //identifier - change this
			  hInst,
			  NULL
		  );

SendMessage(hWndSpin, UDM_SETRANGE, 0, MAKELONG((short)max, (short) min));
SendMessage(hWndSpin, UDM_SETPOS32, 0, (LPARAM)(INT) initial);


		  SetFocus(hWndSpin);

		  return hWndSpin;
}



// some global variables

Graphic *view1;
bool sceneChange;
treeAttributes *tree1;



// uploads data from the spinners into a treeAttributes struct
void GetAttributes() {
		tree1->iterations	= SendMessage(levelSpin, UDM_GETPOS, 0, 0);
		tree1->branchAng1	= (TWOPI/360.0)*SendMessage(angle1Spin, UDM_GETPOS, 0, 0);
		tree1->branchAng2	= (TWOPI/360.0)*SendMessage(angle2Spin, UDM_GETPOS, 0, 0);
		tree1->radDivisions = SendMessage(divSpin, UDM_GETPOS, 0, 0);
		tree1->branchScale	= 0.01 * SendMessage(bScaleSpin, UDM_GETPOS, 0, 0);

		double thick	= 0.01 * SendMessage(thickSpin, UDM_GETPOS, 0, 0);
		double length	= 0.01 * SendMessage(sizeSpin, UDM_GETPOS, 0, 0);
		tree1->scale			= thick;
		tree1->height			= 0.01 * SendMessage(heightSpin, UDM_GETPOS, 0, 0) / thick;
		tree1->segmentLength	= length / thick;
}






// forward definition of WndProc
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);






// WinMain starts here --------------------------------------------------------------------------------------------------------------------------
int WINAPI WinMain(HINSTANCE hInstance,
                   HINSTANCE hPrevInstance,
                   LPSTR lpCmdLine,
                   int nCmdShow)
{
    WNDCLASSEX wc;

    wc.cbSize = sizeof(WNDCLASSEX);
    wc.style          = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc    = WndProc;
    wc.cbClsExtra     = 0;
    wc.cbWndExtra     = 0;
    wc.hInstance      = hInstance;
    wc.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_APPLICATION));
    wc.hCursor        = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground  = (HBRUSH)(COLOR_WINDOW);
    wc.lpszMenuName   = MAKEINTRESOURCE(TOP_MENU);
    wc.lpszClassName  = szWindowClass;
    wc.hIconSm        = LoadIcon(wc.hInstance, MAKEINTRESOURCE(IDI_APPLICATION));

    if (!RegisterClassEx(&wc))
    { MessageBox(NULL, _T("Call to RegisterClassEx failed"),_T("Error"),NULL);
        return 1; }

    hInst = hInstance; 

    HWND hWnd = CreateWindow(
        szWindowClass,
        szTitle,
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT,
        800, 600,
        NULL,
        NULL,
        hInstance,
        NULL
    );

    if (!hWnd)
    {  MessageBox(NULL,_T("Call to CreateWindow failed!"),_T("Error"),NULL);
        return 1;  }


// initialise common controls - compile with Comctl32.lib
	INITCOMMONCONTROLSEX icce;
	icce.dwSize = sizeof(INITCOMMONCONTROLSEX);
	icce.dwICC = ICC_UPDOWN_CLASS | ICC_WIN95_CLASSES;
	if (!InitCommonControlsEx(&icce)) {	 
		MessageBox(hWnd, _T("Call to InitCommonControlsEx Failed"), _T("Error"), NULL); 
	}


// create spinners	
int indent = 5;
levelSpin	= CreateSpinner(hWnd,hInst, indent, 10,  1, 12,		5,		_T("Levels"));
angle1Spin	= CreateSpinner(hWnd,hInst, indent, 40,  0, 180,	30,		_T("Angle 1"));
angle2Spin	= CreateSpinner(hWnd,hInst, indent, 70,  0, 180,	20,		_T("Angle 2"));
heightSpin	= CreateSpinner(hWnd,hInst, indent, 100, 0, 10000,	600,	_T("Trunk length"));
divSpin		= CreateSpinner(hWnd,hInst, indent, 130, 3, 50,		11,		_T("Trunk divisions"));
sizeSpin	= CreateSpinner(hWnd,hInst, indent, 160, 0, 10000,	200,	_T("Branch length"));
thickSpin	= CreateSpinner(hWnd,hInst, indent, 190, 0, 1000,	40,		_T("Branch thickness"));
bScaleSpin	= CreateSpinner(hWnd,hInst, indent, 220, 0, 150,	80,		_T("Taper (%)"));
randSpin	= CreateSpinner(hWnd,hInst, indent, 250, 0, 32767,	0,		_T("Random seed"));

view1 = new Graphic(hWnd);
tree1 = new treeAttributes;
GetAttributes();
SetWorldAxes();
makeTree(vector(0,-12,0),tree1);


ShowWindow(hWnd, nCmdShow);
UpdateWindow(hWnd);

    // Message loop:
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0))
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    return (int) msg.wParam;
}




LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{


    switch (message)
    {
	case WM_CREATE:
		SetTimer(hWnd, 1, 0, NULL);
    break;
    case WM_DESTROY:

		delete tree1;
		delete view1;
		
		clearAll();// delete all verts & faces

        PostQuitMessage(0);
        break;
	case WM_TIMER:
		view1->cam->rotateLeft();
		view1->render();
		break;
	case WM_PAINT:	
		view1->redraw();
		break;
	case WM_SIZE:
		view1->chooseDrawArea(); // resize the main graphic
		break;
	case WM_COMMAND:
		switch (LOWORD(wParam)) {
			case FILE_EXIT:
				PostMessage(hWnd,WM_CLOSE,0,0);
			break;
			case EXPORT_OBJ:
				makeObjFile(hWnd);// export mesh as .obj
			break;
		}
		break;
	case WM_VSCROLL:

		// redraw the tree mesh when a parameter is changed
		clearAll();
		GetAttributes();
		srand(SendMessage(randSpin, UDM_GETPOS, 0, 0));
		try {
		makeTree(vector(0,-12,0),tree1);
		} catch (bad_alloc) {
			MessageBox(hWnd,_T("Out of memory"),_T("Error"),0);
		}
		break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
        break;
    }

    return 0;
}