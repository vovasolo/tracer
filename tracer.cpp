#include <math.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "/opt/root/include/TRandom3.h"
#include "/opt/root/include/TGeoManager.h"
#include "/opt/root/include/TGeoNavigator.h"
#include "/opt/root/include/TGeoMedium.h"
#include "/opt/root/include/TGeoVolume.h"
#include "/opt/root/include/TGeoNode.h"

#include <eigen3/Eigen/Dense>

#include "photon.h"
#include "tracer.h"

#define FARAWAY 1.0e12

using Eigen::Vector3d;
using std::vector;

vector<vector<float> > readColumnData(std::ifstream *file)
{
    // The result of the read is placed in here
    vector<vector<float> > data;
    std::string   line;
    // Read one line at a time into the variable line:
    while(std::getline(*file, line)) {
        vector<float>   lineData;
        std::stringstream  lineStream(line);

        float value;
        while(lineStream >> value) 
            lineData.push_back(value);

        data.push_back(lineData);
    }
    return data;
}

void writeColumnData(vector<vector<float> > *data, std::ofstream &file)
{
    for (auto line = data->begin(); line != data->end(); ++line) {
        for (auto col = line->begin(); col != line->end(); ++col) {
            file << *col << " ";
        }
        file << std::endl;
    }
/* alternatively c++11
    for (auto const line : data) {
        for (auto const &col : line) {
            file << col << " ";
        }
        file << std::endl;
    }
*/
}

Tracer::Tracer()
{
    rndgen = new TRandom3(0);
    ph = new Photon(rndgen);
    logging_on = false;
}

void Tracer::Log1(std::string comment)
{
    if (!logging_on)
        return;
//    std::cout << "Volume:" << nav->GetCurrentVolume()->GetName() << std::endl;
    const double *log_pos = nav->GetCurrentPoint();
    std::cout << comment << " " << log_pos[0] << " " << log_pos[1] << " " << log_pos[2] << " -> ";// std::endl;
    const double *log_dir = nav->GetCurrentDirection();
    std::cout << log_dir[0] << " " << log_dir[1] << " " << log_dir[2] << std::endl;
}

void Tracer::Log2(std::string comment)
{
    if (!logging_on)
        return;
    const double *log_pos = nav2->GetCurrentPoint();
    std::cout << comment << " " << log_pos[0] << " " << log_pos[1] << " " << log_pos[2] << std::endl;
}

bool Tracer::CheckConsistency12()
{
    const double *pos1 = nav->GetCurrentPoint();
    const double *pos2 = nav2->GetCurrentPoint();
    double dx = pos2[0]-pos1[0]; double dy = pos2[1]-pos1[1]; double dz = pos2[2]-pos1[2];
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    if (dist > 1e-10) {
        std::cout << "deltapos: " << dx << " " << dy << " " << dz << " " << dist << std::endl;
        return false;
    }
    return true;
}
void Tracer::SwapNavigators()
{
    TGeoNavigator *tmp = nav;
    nav = nav2;
    nav2 = tmp;
}

// OnePhoton provides an elementary photon tracing loop
// photon is emitted in 4pi from (x, y, z)
// returns status from enum Status
Tracer::Status Tracer::OnePhoton(float x, float y, float z)
{
	int cnt = max_steps; 			// counter of steps

    nav->SetCurrentPoint(x, y, z);
// check if inside trackable volume
	nav->FindNode();
	if (nav->IsOutside()) 
        return outside; // report bad starting point
//let's begin
	ph->Emission();
    if (logging_on)
        std::cout << "*" << std::endl;
    log_status = 4;
    nav->SetCurrentDirection(ph->GetDirection());

jmp_next: // jump here if the photon moves to the next volume 
//    std::cout << "Current: " << nav->GetCurrentVolume()->GetName() << std::endl;
// get optical properties of the current medium
	TGeoMedium *medium = nav->GetCurrentVolume()->GetMedium();
	n = medium->GetParam(0);
	abscoef = medium->GetParam(2);
    emission = medium->GetParam(3);
    raylen = medium->GetParam(4);
// get material index of the current volume	
	mat = nav->GetCurrentVolume()->GetMaterial()->GetIndex();

jmp_same: // jump here if the photon stays in the same volume 
    if (--cnt < 0) return maxsteps; // check for max steps

    nav->FindNextBoundary();

//    std::cout << nextnode->GetVolume()->GetName() << std::endl;
//  logging
    Log1(st_mnem[log_status]);
//    nav->InspectState();

    double step = nav->GetStep(); 			// step is the distance to the next boundary

//check for absorption and Rayleigh on the way to the boundary
    double abspath = abscoef > 0 ? -log(1.-rndgen->Rndm())/abscoef : FARAWAY;	// sampled path to absorption
    double raypath = raylen > 0 ? -log(1.-rndgen->Rndm())*raylen : FARAWAY;		// sampled path to Rayleigh
    if (abspath < step || raypath < step) {		// if photon won't reach the boundary
    // check what gonna happen first
        if (abspath < raypath) {				// absorption - re-emission
            if (rndgen->Rndm() > emission) {
                if (logging_on) {
                    nav->SetStep(abspath);
                    nav->Step(kFALSE);
                    Log1("absmedium");
                }
                return absmedium;    // just absorption
            }
            nav->SetStep(abspath);
            ph->Emission();
        } else {                                  // Rayleigh
            nav->SetStep(raypath);
            ph->Rayleigh();
        }
        nav->Step(kFALSE); // move to the re-emission or Rayleigh position
        nav->SetCurrentDirection(ph->GetDirection());
        goto jmp_same; // continue in the same volume
	}

// photon will reach the boundary, however we don't know yet if it's going to cross it or not
// get material index of the next volume with the help of the auxilary navigator

    nav2->SetCurrentPoint(nav->GetCurrentPoint());
    nav2->SetCurrentDirection(nav->GetCurrentDirection()); // for calculation of normal
    nav2->FindNode();
    nextnode = nav2->FindNextBoundaryAndStep(); // go to the boundary, stop right AFTER crossing
    // check if we have reached World's end
    if (nav2->IsOutside()) {
        Log2("escaped2");
        return escaped;			// photon escaped
    }
//    std::cout << "Next: " << nextnode->GetVolume()->GetName() << std::endl;
    mat_next = nextnode->GetVolume()->GetMaterial()->GetIndex();
// check for overrides
    basic_ovr *ovr = omap[mat] ? omap[mat]->at(mat_next) : 0;
// process what happens on the boundary
// the returned status: -1 if reflected, 0 if absorbed, 1 if transmitted
    BndStatus status = ProcessBoundary(ovr);
//    if (logging_on)
//        std::cout << bnd_mnem[status+1] << std::endl;
    if (status == absorbed) {
        Log2("absorbed");
        return absbndry;
    }
    // restore fStep that was messed (shortened by 1e-6) by nav->Step() in ProcessBoundary()
    nav->SetStep(step);
    if (status == reflected) {  // staying in the same volume
        nav->Step(kTRUE, kFALSE);   // go to the boundary, stop right BEFORE crossing
        nav->SetCurrentDirection(ph->GetDirection());
        goto jmp_same;
    } else {		  // transmitted, continue in the next volume
        // check if we entered a photosensor
        const char* tit = nextnode->GetVolume()->GetTitle();
        if (tit[0] == 'P') {
            int pm_id = nextnode->GetNumber(); // sensor id is stored as node number in the geometry
            if (pm_id >= npms)                 // sanity check
                return pmrange;
            hit_map[pm_id]++;
            Log2("detected");
            return sensor;
        }
        // if not -- continue
        SwapNavigators(); // nav2 is already at the correct starting point -- let's use it
        nav->SetCurrentDirection(ph->GetDirection());

        // check if we have reached World's end
        // this should never nappen as nav2 already have checked
        if (nav->IsOutside()) {
            Log1("escaped1");
            return outside;			// photon escaped
        }
//        std::cout << nextnode->GetVolume()->GetName() << " ***** " << nxt1->GetVolume()->GetName() << std::endl;

        goto jmp_next;
    }
}

// random choice of a boundary process
Tracer::BndProcess Tracer::RollTheDice(basic_ovr *ovr)
{
    float dice = rndgen->Rndm();  // roll the dice
    float t = ovr->loss;
    if (dice < t) return absorption;
    if (dice < (t += ovr->specular)) return specular;
    if (dice < (t += ovr->diffuse)) return diffuse;
    return fresnel;
}

Tracer::BndStatus Tracer::ProcessBoundary(basic_ovr *ovr)
{
    const double *normal;
    Vector3d V3d_normal;
// process is Fresnel if override isn't defined
// otherwise it's randomly chosen according to the probabilities defined in the override
    BndProcess process = ovr ? RollTheDice(ovr) : fresnel;
    log_status = process;
    if (process != absorption) {      // we'll need the normal in every case except absorption
//        const double *pos = nav->GetCurrentPoint();
        normal = nav2->FindNormal();
//        nav->SetCurrentPoint(pos);
//        std::cout << "normal: " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
        V3d_normal << normal[0], normal[1], normal[2];
    }
    switch (process) {
        case absorption:                 // absorbed
            return absorbed;
        case specular:                 // specular reflection
            ph->Mirror(V3d_normal);
            return reflected;
        case diffuse:                 // diffuse scattering
            switch (ovr->diff_model) {
                case 0: //4Pi scattering
                    ph->Emission();  // transmitted if emission happened to the next volume and reflected otherwise
                    return ph->GetDir().dot(V3d_normal) > 0 ? transmitted : reflected;
                case 1: //2Pi lambertian, remaining in the same volume (back scattering)
                    ph->Diffuse(V3d_normal);
                    return reflected;
                case 2: //2Pi lambertian, scattering to the next volume
                    ph->Diffuse(-V3d_normal);
                    return transmitted;
            }
        case fresnel:                 // Fresnel
        // get refractive index of the next volume
            n_next = nextnode->GetVolume()->GetMedium()->GetParam(0);
            bool crossed = ph->Border(V3d_normal, n, n_next); // ???
            return crossed ? transmitted : reflected;
    }
}

bool Tracer::PhotonBomb(float x, float y, float z, int n)
{
	std::fill(hit_map.begin(), hit_map.end(), 0);
	for (int i=0; i<n; i++)
		OnePhoton(x, y, z);
	return true;	
}

bool Tracer::Scan(int nevents, vector<vector<float> > &result)
{
    for (int i=0; i<nevents && i<scanlist.size(); i++) {
        vector <float> evt = scanlist[i];
        PhotonBomb(evt[0], evt[1], evt[2], evt[3]);
        vector <float> hit(hit_map.begin(), hit_map.end());
        hit.push_back(evt[0]);
        hit.push_back(evt[1]);
        hit.push_back(evt[2]);
        hit.push_back(evt[3]);
        result.push_back(hit);
    }
    return true;
}


// olist format:
// From To Loss Specular Diffuse DiffuseModel(int)
// omap can be addressed as omap[mat_from][mat_to]
// omap is a vector of pointers to vectors of pointers to basic_ovr structure
// if there is no overrides from mat. X then omap[X] == 0
// if there are some overrides from mat. X but not to mat. Y then omap[X][Y] == 0

bool Tracer::Init(char *geometry_path, char* override_path)
{
	gm = TGeoManager::Import(geometry_path);

	nav = gm->GetCurrentNavigator();
	if (!nav)
		nav = gm->AddNavigator();
//    nav->GetCache()->BuildIdArray(); // enable navigation by node id
// secondary navigator for probing next volumes
    nav2 = gm->AddNavigator();

	npms = CountPMs();
	hit_map.resize(npms, 0);

// load overrides from file
	std::ifstream filein(override_path);
	olist = readColumnData(&filein);
	nmat = (int) olist[0][0]+0.1;
	for (int i=0; i<nmat; i++)
		omap.push_back(0);
	for (int i=1; i<olist.size(); i++) {
		int m_from = olist[i][0];
		int m_to = olist[i][1];
		if (omap[m_from] == 0) {
			omap[m_from] = new vector <basic_ovr*>;
			for (int j=0; j<nmat; j++)
				omap[m_from]->push_back(0);
		}

		basic_ovr *bo = new basic_ovr(); // initialize to 0
		bo->loss = olist[i][2]; // Loss
		bo->specular = olist[i][3]; // Specular
		bo->diffuse = olist[i][4]; // Diffuse
		bo->diff_model = olist[i][5]; // Diffuse Model

		omap[m_from]->at(m_to) = bo;
	}

	return true;	
}

int Tracer::CountPMs()
{
	TGeoVolume *top = gm->GetTopVolume();
	int nnodes = top->CountNodes();
	int n_pms = 0;
	for (int i=0; i<nnodes; i++) {
		TGeoNode *node = top->GetNode(i);
		if (!node) break;
		const char *tit = top->GetNode(i)->GetVolume()->GetTitle();
		if (tit[0] == 'P')
			n_pms++;
	} 
	return n_pms;
}

bool Tracer::ReadScanData(char *scan_path)
{
	std::ifstream filein(scan_path);
	scanlist = readColumnData(&filein);	
	return true;
}

void Tracer::PrintOmap()
{
	for (int i=0; i<nmat; i++) {
		if (omap[i]) {
			for (int j=0; j<nmat; j++)
				if (omap[i]->at(j)) {
					basic_ovr *bo = omap[i]->at(j);
					std::cout << i << " " << j << ": loss=" << bo->loss << " spec=" << bo->specular;
					std::cout << " diff=" << bo->diffuse << " model=" << bo->diff_model << std::endl;
			}
		}
	}
}
