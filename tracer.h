#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include <string>
#include "/opt/root/include/TRandom3.h"
#include "/opt/root/include/TGeoManager.h"
#include "/opt/root/include/TGeoNavigator.h"
#include "/opt/root/include/TGeoMedium.h"
#include "/opt/root/include/TGeoVolume.h"
#include "/opt/root/include/TGeoNode.h"

using std::vector;

struct basic_ovr { // basic overload
	float loss;
	float specular;
	float diffuse;
	int diff_model;
};

class Photon;

class Tracer {
    public:
        // status returned by OnePhoton()
        enum Status {
            outside = -1,  // initial point is outside the geometry
            pmrange =  -2, // photosensor number out of range
            maxsteps = 1,  // max steps reached
            escaped = 2,   // photon escaped
            absmedium = 3, // absorption by the medium
            absbndry = 4,  // absorption on a boundary
            sensor = 5     // entered a photosensor
        };
        // process on a boundary
        enum BndProcess {
            absorption = 0,
            specular = 1,
            diffuse = 2,
            fresnel= 3
        };
        // status returned by ProcessBoundary()
        enum BndStatus {
            reflected = -1,     // photon stays in the same volume
            absorbed = 0,       // photon is absorbed
            transmitted = 1     // photon passes to the next volume
        };

	public:
        Tracer();
        ~Tracer() {delete rndgen;}

		bool Init(char *geometry_path, char* override_path);
		bool ReadScanData(char *scan_path);
        bool Scan(int nevents, vector<vector<float> > &result);
        void SetLogging(bool log_on) {logging_on = log_on;}
        void Log1(std::string comment);
        void Log2(std::string comment);
        bool CheckConsistency12();
        void SwapNavigators();
        Status OnePhoton(float x, float y, float z);
        BndProcess RollTheDice(basic_ovr *ovr);
        BndStatus ProcessBoundary(basic_ovr *ovr);
		bool PhotonBomb(float x, float y, float z, int n);
		int CountPMs();
		vector <int> GetHitMap();
		void PrintOmap();
		int GetNmat() {return nmat;}
		void SetMaxSteps(int steps) {max_steps = steps;}

private:
		int npms; // total number of photosensors
		int nmat; // total number of materials
		int max_steps; // max number of steps to avoid infinite tracking
		TGeoManager *gm;
        TGeoNavigator *nav; // current navigator
        TGeoNavigator *nav2; // auxilary navigator

		vector < vector < basic_ovr* >* > omap; // override map
		vector<vector<float> > olist; // override list read from file
		vector<vector<float> > scanlist; // list of points to simulate
		TRandom3 *rndgen;
		Photon *ph;
		vector <int> hit_map;

        bool logging_on;
        int log_status;
        vector <std::string> st_mnem = {"absorption", "specular", "diffuse", "fresnel", "emission"};
        vector <std::string> bnd_mnem = {"reflected", "absorbed", "transmitted"};

		TGeoNode *nextnode;				// next node after potential boundary crossing
// properties of the current medium
		int mat;	  // material index
		float n;      // refractive index
		float abscoef; // absorption coefficient (1/MFP)
        float emission; // re-emisson probability
		float raylen; // Rayleigh MFP	
// properties of (potentially) next medium
		int mat_next; // material index		
		float n_next; // refractive index
};

#endif
