//
// This Program is provided by Duke University and the authors as a service to the
// research community. It is provided without cost or restrictions, except for the
// User's acknowledgement that the Program is provided on an "As Is" basis and User
// understands that Duke University and the authors make no express or implied
// warranty of any kind.  Duke University and the authors specifically disclaim any
// implied warranty or merchantability or fitness for a particular purpose, and make
// no representations or warranties that the Program will not infringe the
// intellectual property rights of others. The User agrees to indemnify and hold
// harmless Duke University and the authors from and against any and all liability
// arising out of User's use of the Program.
//
// highMap.h
// Copyright 2005, Austin Eliazar, Ronald Parr, Duke University
//
// Header file for mapping code
// See lowMap.h for details (it is nearly identical)
//

#include "low.h"

extern PMapStarter highMap[H_MAP_WIDTH][H_MAP_HEIGHT];

// The nodes of the ancestry tree are stored here. Since each particle has a unique ID, we can quickly access the particles via their ID
// in this array. See the structure TAncestor for more details.
extern TAncestor h_particleID[H_ID_NUMBER];

// Our current set of particles being processed by the particle filter
extern TParticle h_particle[H_PARTICLE_NUMBER];
// We like to keep track of exactly how many particles we are currently using.
extern int h_cur_particles_used;


void HighInitializeFlags();
void HighInitializeWorldMap();

void HighResizeArray(TMapStarter *node, int deadID);
void HighDeleteObservation(short int x, short int y, short int node);
double HighComputeProb(int x, int y, double distance, int ID);

void HighAddTrace(double startx, double starty, double MeasuredDist, double theta, TAncestor *parent,  int addEnd);
double HighLineTrace(double startx, double starty, double theta, double MeasuredDist, int parentID);
