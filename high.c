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
// high.c
// Copyright 2005 Austin Eliazar, Ronald Parr, Duke University
//
// This file controls all of the major SLAM functions for the high level, in 
// hierarchcal SLAM. It is not used in non-hierarchical SLAM. 
// This code is very nearly an exact duplicate of low.c, and as such has not been
// heavily commented. See low.c for explinations.
//

#include "high.h"
#include "mt-rand.h"

// Threshold for culling particles.  x means that particles with prob. e^x worse
// then the best in the current round are culled
#define H_THRESH 12.0

// Number of passes to use to cull good particles
#define PASSES 9

// Maximum error we will allow for a trace in grid squares
#define MAX_TRACE_ERROR exp(-24.0/HIGH_VARIANCE)
#define WORST_POSSIBLE -10000000

struct TSample_struct {
  float x, y, theta, xG, yG, tG;
  double probability;
  int    parent;
};
typedef struct TSample_struct TSample;


// The number of iterations between writing out the map for video. 0 is off.
int H_VIDEO = 1;

int h_cleanID;
int h_availableID[H_ID_NUMBER];
// No. of children each particle gets
int h_children[H_PARTICLE_NUMBER];

TParticle h_savedParticle[H_PARTICLE_NUMBER];
int h_cur_saved_particles_used;

int h_curGeneration;
unsigned char h_map[H_MAP_WIDTH][H_MAP_HEIGHT];



double LogScorePosition(double x, double y, double theta, int parent, TSense sense)
{
  int i;
  double a, total;

  total = 0.0;
  for (i=0; i < SENSE_NUMBER; i++) {
    a = HighLineTrace(x, y, (sense[i].theta + theta), sense[i].distance, parent);
    total = total + log(MAX(MAX_TRACE_ERROR, a));
  }
  return total;
}



void HighLocalize(TPath *path, TSenseLog *obs)
{
  int i, j, k;
  int best, keepers, worst;
  int newchildren[H_SAMPLE_NUMBER];
  double moveAngle, threshold;
  double ftemp, total;
  TSample sample[H_SAMPLE_NUMBER];
  TPath *holdPath;
  TSenseLog *holdObs;

  holdPath = path;
  holdObs = obs;
  obs = obs->next;
 
 // Make particles
  j = 0;
  for (i=0; i < h_cur_particles_used; i++) {
    while (h_children[i] > 0) {
      h_children[i]--;
      sample[j].parent = i;
      sample[j].probability = 0.0;
      
      // Scatter them
      sample[j].xG = GAUSSIAN(0.8);
      sample[j].yG = GAUSSIAN(0.8);
      sample[j].tG = GAUSSIAN(0.025);
      sample[j].x = h_particle[i].x + sample[j].xG;
      sample[j].y = h_particle[i].y + sample[j].yG;
      sample[j].theta = h_particle[i].theta + sample[j].tG;
      j++;
    }
  }

  threshold = WORST_POSSIBLE;
  j = 0;
  while (path != NULL) {
    HighInitializeFlags();
    keepers = 0;
    best = 0;
    worst = 0;
    for (i=0; i < H_SAMPLE_NUMBER; i++) {
      if (sample[i].probability > threshold) {
	keepers++;
	// Move the particles one step
	moveAngle = sample[i].theta + path->T/2.0;
	sample[i].x = sample[i].x + (TURN_RADIUS * (cos(sample[i].theta + path->T) - cos(sample[i].theta))) +
	  (path->D * cos(moveAngle)) + (path->C * cos(moveAngle + M_PI/2));
	sample[i].y = sample[i].y + (TURN_RADIUS * (sin(sample[i].theta + path->T) - sin(sample[i].theta))) +
	  (path->D * sin(moveAngle)) + (path->C * sin(moveAngle + M_PI/2));
	sample[i].theta = sample[i].theta + path->T;
	
	// Score this step of the obs
	sample[i].probability = sample[i].probability + 
	                        LogScorePosition(sample[i].x, sample[i].y, sample[i].theta, 
						 h_particle[ sample[i].parent ].ancestryNode->ID, obs->sense);
	if (sample[i].probability > sample[best].probability)
	  best = i;
      }
      // Cull bad ones
      else 
	sample[i].probability = WORST_POSSIBLE;
    }

    fprintf(stderr, " ** %d  %.4f     %d\n", best, sample[best].probability, keepers);
    threshold = sample[best].probability - H_THRESH;
    j++;

    // Advance the path & observation
    path = path->next;
    obs = obs->next;
  }

  fprintf(stderr, "High level- Best of %d ", keepers);

  // Normalize
  total = 0.0;
  threshold = sample[best].probability;
  for (i=0; i < H_SAMPLE_NUMBER; i++) {
    if (sample[i].probability == WORST_POSSIBLE)
      sample[i].probability = 0.0;
    else {
      sample[i].probability = exp(sample[i].probability-threshold);
      total = total + sample[i].probability;
    }
  }

  for (i=0; i < H_SAMPLE_NUMBER; i++)
    sample[i].probability = sample[i].probability/total;

  // Count how many children each particle will get in next generation
  for (i = 0; i < H_SAMPLE_NUMBER; i++) 
    newchildren[i] = 0;

  i = j = 0;  // i = no. of survivors, j = no. of new samples
  while ((j < H_SAMPLE_NUMBER) && (i < H_PARTICLE_NUMBER)) {
    k = 0;
    ftemp = MTrandDec();
    while (ftemp > sample[k].probability) {
      ftemp = ftemp - sample[k].probability;
      k++;
    }
    if (newchildren[k] == 0)
      i++;
    newchildren[k]++;
    j++;
  }

  fprintf(stderr, "(%d kept ", i);

  // Now copy over new particles to savedParticles
  best = 0;
  k = 0; // pointer into saved particles
  for (i = 0; i < H_SAMPLE_NUMBER; i++)
    if (newchildren[i] > 0) {
      // We use the parent's x/y/t here because when we update the map, we want to go through 
      // each movement step again
      h_savedParticle[k].x =      h_particle[ sample[i].parent ].x + sample[i].xG;
      h_savedParticle[k].y =      h_particle[ sample[i].parent ].y + sample[i].yG;
      h_savedParticle[k].theta =  h_particle[ sample[i].parent ].theta + sample[i].tG;
      h_savedParticle[k].ancestryNode = h_particle[ sample[i].parent ].ancestryNode;
      h_savedParticle[k].probability = sample[i].probability;
      h_savedParticle[k].ancestryNode->numChildren++;
      h_children[k] = newchildren[i];

      if (h_savedParticle[k].probability > h_savedParticle[best].probability) 
	best = k;
      k++;
    }

  // This number records how many saved particles we are currently using, so that we can ignore anything beyond this
  // in later computations.
  h_cur_saved_particles_used = k;

  // We might need to continue generating children for particles, if we reach PARTICLE_NUMBER worth of distinct parents early
  // We renormalize over the chosen particles, and continue to sample from there.
  if (j < H_SAMPLE_NUMBER) {
    // Normalize particle probabilities. Note that they have already been exponentiated
    total = 0.0;
    for (i = 0; i < h_cur_saved_particles_used; i++) 
      total = total + h_savedParticle[i].probability;

    for (i=0; i < h_cur_saved_particles_used; i++)
      h_savedParticle[i].probability = h_savedParticle[i].probability/total;

    total = 0.0;
    for (i = 0; i < h_cur_saved_particles_used; i++) 
      total = total + h_savedParticle[i].probability;

    while (j < H_SAMPLE_NUMBER) {
      k = 0;
      ftemp = MTrandDec()*total;
      while (ftemp > (h_savedParticle[k].probability)) {
	ftemp = ftemp - h_savedParticle[k].probability;
	k++;
      }
      h_children[k]++;

      j++;
    }
  }
}



void HighAddToWorldModel(TPath *sourcePath, TSenseLog *sourceObs, int maxID)
{
  int i, ID;
  double moveAngle;
  TPath *path;
  TSenseLog *obs;

  path = sourcePath;
  obs = sourceObs->next;

  while (path != NULL) {
    HighInitializeFlags();
    for (ID=0; ID < maxID; ID++) {
      // Move the particle one step
      moveAngle = h_particle[ID].theta + (path->T/2.0);
      h_particle[ID].x = h_particle[ID].x + (TURN_RADIUS * (cos(h_particle[ID].theta + path->T) - cos(h_particle[ID].theta))) +
			 (path->D * cos(moveAngle)) + (path->C * cos(moveAngle + M_PI/2));
      h_particle[ID].y = h_particle[ID].y + (TURN_RADIUS * (sin(h_particle[ID].theta + path->T) - sin(h_particle[ID].theta))) +
			 (path->D * sin(moveAngle)) + (path->C * sin(moveAngle + M_PI/2));
      h_particle[ID].theta = h_particle[ID].theta + path->T;

      for (i=0; i < SENSE_NUMBER; i++) {
	// normalize readings relative to the pose of current assumed position
	HighAddTrace(h_particle[ID].x, h_particle[ID].y, obs->sense[i].distance, (obs->sense[i].theta + h_particle[ID].theta), 
		     h_particle[ID].ancestryNode, (obs->sense[i].distance < MAX_SENSE_RANGE));
      }
    }

    path = path->next;
    obs = obs->next;
  }

  HighInitializeFlags();
  if (obs != NULL) 
    for (ID=0; ID < maxID; ID++) {
      for (i=0; i < SENSE_NUMBER; i++) {
	// normalize readings relative to the pose of current assumed position
	HighAddTrace(h_particle[ID].x, h_particle[ID].y, obs->sense[i].distance, (obs->sense[i].theta + h_particle[ID].theta), 
		     h_particle[ID].ancestryNode, (obs->sense[i].distance < MAX_SENSE_RANGE));
      }
    }
}



//
// UpdateAncestry
//
// When a particle has been deemed fit to spawn :
//   a) determine how many children each particle gets. 
//   b) generate new particle IDs for thier children. only child inherits the parent ID
//   c) make a new ancestor entry for all particles spawning more than one child. includes :
//          - ID, generation, number of children, parent pointer, list of pointers to modified grid squares
//   d) add the particle's data to the map
//   e) prune 'dead' ancestry
//
void HighUpdateAncestry(TPath *path, TSenseLog *obs)
{
  int i, j;
  TAncestor *temp, *hold;
  TEntryList *entry, *workArray;
  TMapStarter *node;

  // Go through the current particle array, and prune out all particles that did not spawn any particles. We know that the particle
  // had to spawn samples, but those samples may not have become particles themselves, through not generating any samples for the 
  // next generation. Recurse up through there.
  for (i=0; i < h_cur_particles_used; i++) {
    temp = h_particle[i].ancestryNode;

    while (temp->numChildren == 0) {
      // Free up memory
      for (j=0; j < temp->total; j++)
	HighDeleteObservation(temp->mapEntries[j].x, temp->mapEntries[j].y, temp->mapEntries[j].node);

      free(temp->mapEntries);
      temp->mapEntries = NULL;

      // Recover the ID. 
      h_cleanID++;
      h_availableID[h_cleanID] = temp->ID;
      temp->generation = h_curGeneration;
      temp->ID = -42;

      hold = temp;
      temp = temp->parent;
      hold->parent = NULL;

      // Note the disappearance of this particle (may cause telescoping of particles, or outright deletion)
      temp->numChildren--;
    }
  }

  // Run through the particle IDs, checking for those that are in use. Of those, look at their parent.
  // Should the parent have only one child, give all of the altered map squares of the parent to the child, 
  // and dispose of the parent. 
  for (i = 0; i < H_ID_NUMBER-1; i++) {
    // These booleans mean (in order) the ID is in use, it has a parent (ie is not the root of the ancestry tree),
    // and that its parent has only one child (this ID) 
    if ((h_particleID[i].ID == i) && (h_particleID[i].parent != NULL) && (h_particleID[i].parent->numChildren == 1)) {
      while (h_particleID[i].parent->generation == -111)
	h_particleID[i].parent = h_particleID[i].parent->parent;

      // We need to collapse one of the branches on our ancestry tree. 
      // This entails moving all of the observed squares from the child up to the parent
      temp = h_particleID[i].parent;

      // We need to collapse one of the branches on our ancestry tree.
      // This entails moving all of the observed squares from the child up to the parent

      // Check to make sure that the parent's array is large enough to accomadate all of the entries of the child as well.
      if (temp->size < (temp->total + h_particleID[i].total)) {
	temp->size = (int)(ceil((temp->size + h_particleID[i].size)*1.75));
	workArray = (TEntryList *)malloc(sizeof(TEntryList)*temp->size);
	if (workArray == NULL) fprintf(stderr, "Malloc failed for workArray\n");


	for (j=0; j < temp->total; j++) {
	  workArray[j].x = temp->mapEntries[j].x;
	  workArray[j].y = temp->mapEntries[j].y;
	  workArray[j].node = temp->mapEntries[j].node;
	}
	// Note that temp->total hasn't changed- that will grow as the child's entries are added in
	free(temp->mapEntries);
	temp->mapEntries = workArray;
      }

      // Change all map entries of the child to have the ID of the parent
      // Also check to see if this entry supercedes an entry currently attributed to the parent.
      // Since collapses can merge all of the entries between the parent and the current child into the parent, this check is performed 
      // by comparing to see if the generation of the last observation (before the child's update) is at least as recent as parent's 
      // generation. If so, note that there is another "dead" entry in the observation array. It will be cleaned up later. If this puts 
      // the total number of used slot, minus the number of "dead", below the threshold, shrink the array (which cleans up the dead)
      entry = h_particleID[i].mapEntries;
      for (j=0; j < h_particleID[i].total; j++) {
	// Change the ID
	node = highMap[entry[j].x][entry[j].y];
	node->array[entry[j].node].ID = temp->ID;
	node->array[entry[j].node].source = temp->total;

	temp->mapEntries[temp->total].x = entry[j].x;
	temp->mapEntries[temp->total].y = entry[j].y;
	temp->mapEntries[temp->total].node = entry[j].node;
	temp->total++;

	// Check for pre-existing observation in the parent's list
	if (node->array[entry[j].node].parentGen >= temp->generation) {
	  node->array[entry[j].node].parentGen = -1;
	  node->dead++;
	}
      }

      // We do this in a second pass for a good reason. If there are more than one update for a given grid square which uses the child's
      // ID (as a consequence of an earlier collapse), then we want to make certain that the resizing doesn't take place until after all
      // entries have changed their ID appropriately.
      for (j=0; j < h_particleID[i].total; j++) {
	node = highMap[entry[j].x][entry[j].y];
	if ((int)((node->total - node->dead)*2.5) < node->size) 
	  HighResizeArray(node, -7);
      }

      // We're done with it- remove the array of updates from the child.
      free(entry);
      h_particleID[i].mapEntries = NULL;

      // Inherit the number of children
      temp->numChildren = h_particleID[i].numChildren;

      // Subtlety of the ancestry tree: since we only keep pointers up to the parent, we can't exactly change all of the
      // descendents of the child to now point to the parent. What we can do, however, is mark the change for later, and
      // update all of the ancestor particles in a single go, later. That will take a single O(P) pass
      h_particleID[i].generation = -111;
    }
  }

  // This is the step where we correct for redirections that arise from the collapse of a branch of the ancestry tree
  for (i=0; i < H_ID_NUMBER-1; i++) 
    if (h_particleID[i].ID == i) 
      while (h_particleID[i].parent->generation == -111) 
	h_particleID[i].parent = h_particleID[i].parent->parent;

  // Wipe the slate clean, so that we don't get confused by the mechinations of the previous changes
  // from deletes and merges. Updates can make thier own tables, as needed.
  //
  // When deleting, we could have a chain of two ancestors removed, such that the child relied on the 
  // entry of the parent in the ObsArray. However, as it stands, only the parent would be able to change
  // it's ObsArray entry to NULL. Then, when a new particle gets the ID of that dead child, we wouldn't 
  // know when the ObsArray entry came from a valid observation, or from the phantom of the dead child
  // -- We could get around this by realizing that the only way a child would ever look to a grid square, 
  // -- and thus would have the entry filled in, at this point would be for that particle having updated 
  // -- the square on a previous pass. That is, if the ObsArray entry points to a node with ID equal to 
  // -- this particle's ID, we know its safe (because the phantom only could point to a parent). If the 
  // -- node has any other ID, we could scrub it, and start as if it were NULL.
  //
  // -- So basically, the argument of above is able to extended into a general proof. Implement this later
  //

  j = 0;
  // Add the current savedParticles into the ancestry tree, and copy them over into the 'real' particle array
  for (i = 0; i < h_cur_saved_particles_used; i++) {
    while (h_savedParticle[i].ancestryNode->generation == -111) 
      h_savedParticle[i].ancestryNode = h_savedParticle[i].ancestryNode->parent;

    if (h_savedParticle[i].ancestryNode->numChildren == 1) {
      h_savedParticle[i].ancestryNode->generation = h_curGeneration;
      h_savedParticle[i].ancestryNode->numChildren = 0;
      h_particle[j].ancestryNode = h_savedParticle[i].ancestryNode;

      h_particle[j].x = h_savedParticle[i].x;
      h_particle[j].y = h_savedParticle[i].y;
      h_particle[j].theta = h_savedParticle[i].theta;
      h_particle[j].probability = h_savedParticle[i].probability;
      j++;
    }
    else if (h_savedParticle[i].ancestryNode->numChildren > 0) {
      temp = &(h_particleID[ h_availableID[h_cleanID] ]);
      temp->ID = h_availableID[h_cleanID];
      h_cleanID--;

      temp->parent = h_savedParticle[i].ancestryNode;
      temp->mapEntries = NULL;
      temp->total = 0;
      temp->size = 0;
      temp->generation = h_curGeneration;
      temp->numChildren = 0;
      temp->path = NULL;
      temp->seen = 0;

      if (h_cleanID < 0) {
	fprintf(stderr, " !!! Insufficient Number of Particle IDs : Abandon Ship !!!\n");
	h_cleanID = 0;
      }

      h_particle[i].ancestryNode = temp;
      h_particle[j].x = h_savedParticle[i].x;
      h_particle[j].y = h_savedParticle[i].y;
      h_particle[j].theta = h_savedParticle[i].theta;
      h_particle[j].probability = h_savedParticle[i].probability;
      j++;
    }
  }

  h_cur_particles_used = h_cur_saved_particles_used;

  HighAddToWorldModel(path, obs, h_cur_particles_used);

  // Clean up the ancestry particles which disappeared in branch collapses. Also, recover their IDs.
  for (i=0; i < H_ID_NUMBER-1; i++) 
    if (h_particleID[i].generation == -111) {
      h_particleID[i].generation = -1;
      h_particleID[i].numChildren = 0;
      h_particleID[i].parent = NULL;
      h_particleID[i].mapEntries = NULL;
      h_particleID[i].path = NULL;
      h_particleID[i].seen = 0;
      h_particleID[i].total = 0;
      h_particleID[i].size = 0;

      // Recover the ID. 
      h_cleanID++;
      h_availableID[h_cleanID] = i;
      h_particleID[i].ID = -3;
    }
}



//
// PrintMap
//
// name - name of map file
// particles - flag to indicate if particles should be shown
//
void HighPrintMap(char *name, TAncestor *parent)
{
  FILE *printFile;
  int x, y;
  int width, height;
  int startx, starty, lastx, lasty;
  char sysCall[128];
  double hit;

  width = H_MAP_WIDTH;
  height = H_MAP_HEIGHT;

  for(x=0; x < width; x++)
    for(y=0; y < height; y++)
      h_map[x][y] = 0;

  lastx = 0;
  lasty = 0;
  startx = width-1;
  starty = height-1;

  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      hit = HighComputeProb(x, y, 1.4, parent->ID);
      if (hit == UNKNOWN) 
	h_map[x][y] = 255;
      else {
	h_map[x][y] = (int) (230 - (hit * 230));
	if (x > lastx)
	  lastx = x;
	if (y > lasty)
	  lasty = y;
	if (x < startx)
	  startx = x;
	if (y < starty)
	  starty = y;
      }
    }
  }  

  sprintf(sysCall, "%s.ppm", name);
  printFile = fopen(sysCall, "w");
  fprintf(printFile, "P6\n # particles.ppm \n %d %d\n",
	  lastx-startx+1, lasty-starty+1);
  fprintf(printFile, "255\n");

  for (y = lasty; y >= starty; y--) 
    for (x = startx; x <= lastx; x++) {
      if (h_map[x][y] == 254) 
	fprintf(printFile, "%c%c%c", 255, 0, 0);
      else if (h_map[x][y] == 253) 
	fprintf(printFile, "%c%c%c", 0, 255, 200);
      else if (h_map[x][y] == 252) 
	fprintf(printFile, "%c%c%c", 255, 55, 55);
      else if (h_map[x][y] == 251) 
	fprintf(printFile, "%c%c%c", 50, 150, 255);
      else
	fprintf(printFile, "%c%c%c", h_map[x][y], h_map[x][y], h_map[x][y]);
    }
      
  fclose(printFile);
  sprintf(sysCall, "convert %s.ppm %s.png", name, name);
  system(sysCall);
  sprintf(sysCall, "chmod 666 %s.ppm", name);
  system(sysCall);
  sprintf(sysCall, "chmod 666 %s.png", name);
  system(sysCall);
  fprintf(stderr, "High map dumped to file\n");
}



void InitHighSlam()
{
  int i;

  // Initialize the worldMap
  HighInitializeWorldMap();

  // Initialize the ancestry and particles
  h_cleanID = H_ID_NUMBER - 2;    // ID_NUMBER-1 is being used as the root of the ancestry tree.

  // Initialize all of our unused ancestor particles to look unused.
  for (i = 0; i < H_ID_NUMBER; i++) {
    h_availableID[i] = i;

    h_particleID[i].total = 0;
    h_particleID[i].size = 0;
    h_particleID[i].generation = -1;;
    h_particleID[i].numChildren = 0;
    h_particleID[i].ID = -1;
    h_particleID[i].parent = NULL;
    h_particleID[i].mapEntries = NULL;
    h_particleID[i].path = NULL;
    h_particleID[i].seen = 0;
  }

  // Initialize the root of our ancestry tree.
  h_particleID[H_ID_NUMBER-1].generation = 0;
  h_particleID[H_ID_NUMBER-1].numChildren = 1;
  h_particleID[H_ID_NUMBER-1].ID = H_ID_NUMBER-1;
  h_particleID[H_ID_NUMBER-1].parent = NULL;
  h_particleID[H_ID_NUMBER-1].mapEntries = NULL;
  h_particleID[H_ID_NUMBER-1].size = 0;
  h_particleID[H_ID_NUMBER-1].total = 0;

  // Create all of our starting particles at the center of the map.
  for (i = 0; i < H_PARTICLE_NUMBER; i++) {
    h_particle[i].ancestryNode = &(h_particleID[H_ID_NUMBER-1]);
    h_particle[i].x = H_MAP_WIDTH / 2;
    h_particle[i].y = (H_MAP_HEIGHT / 2) + 100;
    h_particle[i].theta = 0.001;
    h_particle[i].probability = 0;
    h_children[i] = 0;
  }
  // We really only use the first particle, since they are all essentially the same.
  h_particle[0].probability = 1;
  h_cur_particles_used = 1;
  h_children[0] = H_SAMPLE_NUMBER;

  // We don't need to initialize the savedParticles, since Localization will create them for us, and they first are used in 
  // UpdateAncestry, which is called after Localization. This statement isn't necessary, then, but serves as a sort of placeholder 
  // when reading the code.
  h_cur_saved_particles_used = 0;

  h_curGeneration = 0;
}



void CloseHighSlam()
{
}



void HighSlam(TPath *path, TSenseLog *obs)
{
  int i, j;
  char name[16];

  HighInitializeFlags();

  if (h_curGeneration == 0) {
    HighAddToWorldModel(path, obs, 1);

    sprintf(name, "hmap00");
    HighPrintMap(name, h_particle[0].ancestryNode);
    sprintf(name, "rm hmap00.ppm");
    system(name);
  }
  else {
    // Localize off of the path
    HighLocalize(path, obs);

    HighUpdateAncestry(path, obs);

    if ((H_VIDEO) && (h_curGeneration % H_VIDEO == 0)) {
      sprintf(name, "hmap%.2d", (int) (h_curGeneration/H_VIDEO));
      j = 0;
      for (i = 0; i < h_cur_particles_used; i++)
	if (h_particle[i].probability > h_particle[j].probability)
	  j = i;

      HighPrintMap(name, h_particle[j].ancestryNode);
      sprintf(name, "rm hmap%.2d.ppm", (int) (h_curGeneration/H_VIDEO));
      system(name);
    }
  }

  h_curGeneration++;
  HighInitializeFlags();
}



