#ifndef PHOTON_H
#define PHOTON_H

#include <math.h>
#include "/opt/root/include/TRandom3.h"

#include <eigen3/Eigen/Dense>

using Eigen::Vector3d;

class Photon {
	public:
		Photon() { rndgen = new TRandom3(0); my_rg = true; }
		Photon( TRandom3 *rg ) { rndgen = rg; my_rg = false; }
		~Photon() { if (my_rg) delete rndgen; }

        void SetDir (const double dx, const double dy, const double dz) {dir(0) = dx; dir(1) = dy; dir(2) = dz;}
        void SetDir (const Vector3d v) { dir = v; }

        const double GetDx() { return dir(0); }
        const double GetDy() { return dir(1); }
        const double GetDz() { return dir(2); }
        const Vector3d GetDir() { return dir; }
// for interfacing with ROOT TGeoNavigator
        const double *GetDirection() { direction[0]=dir(0); direction[1]=dir(1); direction[2]=dir(2); return direction;}
        const double *Move(const double *origin, double step);

		void Emission();
		void Rayleigh();
        void Mirror(Vector3d normal);
        void Diffuse(Vector3d normal);
        void Refraction(Vector3d normal, double n0, double n1);
        bool Border(Vector3d normal, double n0, double n1);
		void Print();

	private:
        Vector3d dir;			// Direction
        Vector3d old_dir;		// Previous direction
        double c, s, n;		// Temporary variables
        double s2, sc, cs, sp, sm, r;	// Temporary variables for Fresnel formulae
        Vector3d dir_n; 		// Temporary: normal component of dir

		void random_dir();
		TRandom3 *rndgen;
		bool my_rg;
// for interfacing with ROOT TGeoNavigator
        double direction[3];
        double point[3];
};

inline void Photon::random_dir()
{
	/* This is a variant of the algorithm for computing a random point
   * on the unit sphere; the algorithm is suggested in Knuth, v2,
   * 3rd ed, p136; and attributed to Robert E Knop, CACM, 13 (1970),
   * 326.
   * Taken from GNU scientific library v1.6
   */

	//Begin with the polar method for getting x,y inside a unit circle

	do {
		dir(0) = -1. + 2. * rndgen->Rndm();
		dir(1) = -1. + 2. * rndgen->Rndm();
		dir(2) = 0.;
		c = dir.squaredNorm();
    } while (c > 1.0);

	dir(2) = -1. + 2. * c;              // z uniformly distributed from -1 to 1
	s = 2. * sqrt (1. - c);         // factor to adjust x,y so that x^2+y^2 is equal to 1-z^2
	dir(0) *= s;
	dir(1) *= s;
}

// Isotropic emission
inline void Photon::Emission()
{
	random_dir();
}

// Perfect specular reflection
inline void Photon::Mirror(Vector3d normal)
{
	c = dir.dot(normal);
	dir -= normal*(c*2.);
}

// Perfect Lambertian reflection
inline void Photon::Diffuse(Vector3d normal)
{
    //register double norm2;

	do {
	    random_dir();
	    dir -= normal;
	} while (dir.squaredNorm() < 0.000001);
	dir.normalize();
}

// Rayleigh scattering
inline void Photon::Rayleigh()
{
	old_dir = dir;
	do {
		random_dir();
		c = dir.dot(old_dir);
	} while ( rndgen->Rndm()*2 > (c*c + 1.));
}

// Refraction
inline void Photon::Refraction(Vector3d normal, double n0, double n1)
{
	c = dir.dot(normal);
	dir_n = normal * c; // dir_n = normal component of dir
	dir -= dir_n; 	// dir = tangential component of dir
	dir *= n0/n1; 	// dir = tangential component of new direction
	dir += normal * sqrt(1.-dir.squaredNorm()); //  scaled normal gives normal component of new direction
}

// Interface between two dielectrics
inline bool Photon::Border(Vector3d normal, double n0, double n1)
// returns: true if crossed, false if reflected
// probabilities are calculated according to 
// Fresnel equations for non-polarized light
{
	n = n0/n1;
	c = dir.dot(normal);
	s2 = 1. - c*c;
	s = sqrt(s2);
	if (s*n > 1.) {
		Mirror(normal);	// Full internal reflection
		return false;
	}
	sc = s*sqrt(1. - n*n*s2);
	cs = c*s*n;

	sp = sc + cs; sp *= sp;
	sm = sc - cs; sm *= sm;

	r = sm/sp*(1.+(1.-sp)/(1.-sm))/2.;

	if (rndgen->Rndm() < r) {
		Mirror(normal);	// Reflected
		return false;
	}

	Refraction(normal, n0, n1);
	return true;
}

inline const double *Photon::Move(const double *origin, double step)
{
    Vector3d trans = dir*step;
    point[0] = origin[0] + trans(0);
    point[1] = origin[1] + trans(1);
    point[2] = origin[2] + trans(2);
    return point;
}

#endif
