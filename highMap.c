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
// highMap.c
//
// Copyright 2005, Austin Eliazar, Ronald Parr, Duke University
//
// Code for generating and maintaining maps at the high level for hierarchical SLAM.
// Code is nearly identical to lowMap.c, and therefore isn't heavily commented. Go 
// there to see an explination of the functions.
//

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <math.h>
#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "highMap.h"

#define H_PRIOR_DIST 4.0
#define H_PRIOR (-1.0/(MAP_SCALE*8.0))


PMapStarter highMap[H_MAP_WIDTH][H_MAP_HEIGHT];
// The nodes of the ancestry tree are stored here. Since each particle has a unique ID, we can quickly access the particles via their ID
// in this array. See the structure TAncestor for more details.
TAncestor h_particleID[H_ID_NUMBER];
// Our current set of particles being processed by the particle filter
TParticle h_particle[H_PARTICLE_NUMBER];
// We like to keep track of exactly how many particles we are currently using.
int h_cur_particles_used;



void HighInitializeFlags()
{
  while (observationID > 0) {
    observationID--;
    flagMap[obsX[observationID]][obsY[observationID]] = 0;
    obsX[observationID] = 0;
    obsY[observationID] = 0;
  }
  observationID = 1;
}


//
// Initializes the highMap
// Always returns 0 to indicate that it was successful.
//
void HighInitializeWorldMap()
{
  int x, y;

  for (y=0; y < H_MAP_HEIGHT; y++)
    for (x=0; x < H_MAP_WIDTH; x++) {
      highMap[x][y] = NULL;
      flagMap[x][y] = 0;
    }

  for (x=0; x < AREA; x++) {
    obsX[x] = 0;
    obsY[x] = 0;
  }

  observationID = 1;
}


void HighDestroyMap()
{
  int x, y;

  // Get rid of the old map.
  for (y=0; y < H_MAP_HEIGHT; y++)
    for (x=0; x < H_MAP_WIDTH; x++) 
      while (highMap[x][y] != NULL) {
	free(highMap[x][y]->array);
	free(highMap[x][y]);
	highMap[x][y] = NULL;
      }
}



void HighResizeArray(TMapStarter *node, int deadID)
{
  short int i, j, ID, x, y;
  short int hash[H_ID_NUMBER];
  int source, last;
  TMapNode *temp;

  if (deadID >= 0)
    node->dead++;

  // Don't count the dead entries in computing the new size
  node->size = (int)(ceil((node->total - node->dead)*1.75));
  temp = (TMapNode *) malloc(sizeof(TMapNode)*node->size);
  if (temp == NULL) fprintf(stderr, "Malloc failed in expansion of arrays\n");

  for (i=0; i < H_ID_NUMBER; i++)
    hash[i] = -1;

  j = 0;
  for (i=0; i < node->total; i++) {
    if (node->array[i].ID == deadID) {
      // Denote that this has been removed already. Therefore, we won't try to remove it later.
      // We don't bother actually removing the source, since the only way that we can have a deadID is if we are in
      // the process of removing all updates that deadID has made.
      h_particleID[deadID].mapEntries[node->array[i].source].node = -1;
    }

    // This observation is the first one of this ID entered into the new array. Just copy it over, and note its position.
    else if (hash[node->array[i].ID] == -1) {
      // Copy the information into the new array.
      temp[j].ID = node->array[i].ID;
      temp[j].source = node->array[i].source;
      temp[j].parentGen = node->array[i].parentGen;
      temp[j].hits = node->array[i].hits;
      temp[j].distance = node->array[i].distance;

      // This entry is moving- alter its source to track it
      if (h_particleID[ temp[j].ID ].mapEntries[ temp[j].source ].node != i) {
	fprintf(stderr, "HUWAK !!%d (%d %d)\n", i, temp[j].ID, temp[j].source);
	exit(-1);
      }
      h_particleID[ temp[j].ID ].mapEntries[ temp[j].source ].node = j;

      // Note that an observation with this ID has already been entered into the new array, and where that was entered.
      hash[node->array[i].ID] = j;
      j++;
    }

    else if (node->array[i].distance > temp[hash[node->array[i].ID]].distance) {
      ID = node->array[i].ID;
      source = temp[hash[ID]].source;

      // Remove the source of the dead entry
      h_particleID[ID].total--;
      last = h_particleID[ID].total;

      h_particleID[ID].mapEntries[source].x = h_particleID[ID].mapEntries[last].x;
      h_particleID[ID].mapEntries[source].y = h_particleID[ID].mapEntries[last].y;
      h_particleID[ID].mapEntries[source].node = h_particleID[ID].mapEntries[last].node;

      // A source entry was moved. Make sure that the observation it links to notes the new source position.
      x = h_particleID[ID].mapEntries[source].x;
      y = h_particleID[ID].mapEntries[source].y;

      if ((highMap[x][y] == node) && (h_particleID[ID].mapEntries[source].node < i))
	temp[hash[ID]].source = source;
      else
	highMap[x][y]->array[ h_particleID[ID].mapEntries[source].node ].source = source;

      // Copy the more recent information into the slot previously held by the dead entry
      temp[hash[ID]].source = node->array[i].source;
      temp[hash[ID]].hits = node->array[i].hits;
      temp[hash[ID]].distance = node->array[i].distance;
      // We do not copy over the parentGen- we are inheriting it from the dead entry, since it was the predecessor
      // The ID does not need to be copied, since it was necessarily the same for both observations.

      // This entry is moving- alter its source to track it
      h_particleID[ID].mapEntries[ node->array[i].source ].node = hash[ID];
    }

    // There was already an entry for this ID. This new entry is an older form of the observation already recorded. Therefore, 
    // the new entry is dead, and should not be copied over, and it's source in the ancestry tree should be removed.
    else {
      // The new entry is an older form of the one already entered. We should inherit the new parentGen
      if (node->array[i].parentGen != -1)
	temp[hash[node->array[i].ID]].parentGen = node->array[i].parentGen;

      ID = node->array[i].ID;
      source = node->array[i].source;

      // Remove the source of the dead entry
      h_particleID[ID].total--;
      last = h_particleID[ID].total;

      if (last != source) {
	h_particleID[ID].mapEntries[source].x = h_particleID[ID].mapEntries[last].x;
	h_particleID[ID].mapEntries[source].y = h_particleID[ID].mapEntries[last].y;
	h_particleID[ID].mapEntries[source].node = h_particleID[ID].mapEntries[last].node;

	// A source entry was moved. Make sure that the observation it links to notes the new source position.
	x = h_particleID[ID].mapEntries[source].x;
	y = h_particleID[ID].mapEntries[source].y;

	if ((highMap[x][y] == node) && (h_particleID[ID].mapEntries[source].node <= i))
	  temp[hash[ID]].source = source;
	else
	  highMap[x][y]->array[ h_particleID[ID].mapEntries[source].node ].source = source;
      }
    }

  }

  // Note the new total, which should be the previous size minus the dead.
  node->total = j;
  // After completing this process, we have removed all dead entries.
  node->dead = 0;
  free(node->array);
  node->array = temp;
}



// When we add a new entry to workingArray, there is a chance
// that we will run into a dead entry. If so, we will need to delete the dead entry, by copying the last entry onto its
// location. We then need to recursively add the entry (that we just copied onto that spot) to the workingArray
static void AddToWorkingArray(int i, TMapStarter *node, short int workingArray[]) 
{
  int j, source, last;
  TEntryList *entries;

  // Keep an eye out for dead entries. They will be made apparent when two entries both have the same ID.
  if (workingArray[node->array[i].ID] == -1) 
    workingArray[node->array[i].ID] = i;

  else {
    // The node we are currently looking at is the dead one.
    if (node->array[i].distance < node->array[ workingArray[node->array[i].ID] ].distance) {
      // Otherwise, remove the source, then remove the entry. Follow with a recursive call.
      j = i;
      if (node->array[i].parentGen >= 0)
	node->array[ workingArray[node->array[i].ID] ].parentGen = node->array[i].parentGen;
    }

    // The previously entered entry is outdated. Replace it with this newer one.
    else {
      j = workingArray[node->array[i].ID];
      workingArray[node->array[i].ID] = i;
      if (node->array[j].parentGen >= 0)
	node->array[i].parentGen = node->array[j].parentGen;
    }
    
    // The node identified as "j" is dead. Remove its entry from the list of altered squares in the ancestor tree.
    h_particleID[node->array[j].ID].total--;

    entries = h_particleID[node->array[j].ID].mapEntries;
    source = node->array[j].source;
    last = h_particleID[node->array[j].ID].total;

    if (last != source) {
      entries[source].x = entries[last].x;
      entries[source].y = entries[last].y;
      entries[source].node = entries[last].node;

      // Somewhat confusing- we just removed an entry from the list of altered squares maintained by an ancestor particle (entries)
      // This means moving an entry from the end of that list to the spot which was vacated (entries[h_particleID[node->array[j].ID].total])
      // Therefore, the entry in the map corresponding to that last entry needs to point to the new entry.
      highMap[ entries[source].x ][ entries[source].y ]->array[ entries[source].node ].source = source;
    }    

    // Now remove the node itself
    node->total--;
    node->dead--;

    if (node->dead < 0) {
      fprintf(stderr, "DEATH! DESTRUCTION! %d : %d %d\n", node->dead, 
	      h_particleID[ node->array[0].ID ].mapEntries[ node->array[0].source ].x, 
	      h_particleID[ node->array[0].ID ].mapEntries[ node->array[0].source ].y);
      exit(-1);
    }

    if (j != node->total) {
      node->array[j].parentGen = node->array[node->total].parentGen;
      node->array[j].distance = node->array[node->total].distance;
      node->array[j].source = node->array[node->total].source;
      node->array[j].hits = node->array[node->total].hits;
      node->array[j].ID = node->array[node->total].ID;
      // We just moved the last entry in the list to position j. Update it's source entry in the ancestry tree to reflect its new position
      h_particleID[ node->array[j].ID ].mapEntries[ node->array[j].source ].node = j;

      // If the entry we just moved was in workingArray, we need to correct workingArray.
      // Also, we know that since it has been entered already, we don't need to enter it again
      if (workingArray[node->array[j].ID] == node->total)
	workingArray[node->array[j].ID] = j;
      else if (i != node->total) 
	// Final step- add this newly copied node to the working array (we don't want it skipped over)
	AddToWorkingArray(j, node, workingArray);
    }

  }
}



inline void HighBuildObservation(int x, int y, char usage)
{
  TAncestor *lineage;
  PAncestor stack[H_PARTICLE_NUMBER];
  short int workingArray[H_ID_NUMBER+1];
  int i, here, topStack;
  char flag;

  if (observationID >= AREA)
    fprintf(stderr, "cRoll over!\n");

  // Grab a slot
  flagMap[x][y] = observationID;
  obsX[observationID] = x;
  obsY[observationID] = y;
  observationID++;
  // Display ownership of this slot
  here = flagMap[x][y];

  // Initialize the slot and the ancestor particles
  for (i=0; i < H_ID_NUMBER; i++) 
    observationArray[here][i] = -1;

  for (i=0; i < H_ID_NUMBER; i++) {
    workingArray[i] = -1;
    h_particleID[i].seen = 0;
  }
  
  if (usage) {
    flag = 1;
    for (i=0; i < highMap[x][y]->total; i++) 
      if (highMap[x][y]->array[i].hits > 0)
	flag = 0;
  }

  // Fill in the areas of the array that made direct observations
  for (i=0; i < highMap[x][y]->total; i++) 
    AddToWorkingArray(i, highMap[x][y], workingArray);

  // Fill in the holes in the observation array, by using the value of their parents
  for (i=0; i < h_cur_particles_used; i++) {
    lineage = h_particle[i].ancestryNode;
    topStack = 0;

    // Eventually we will either get to an ancestor that we have already seen,
    // or we will hit the top of the tree (and thus its parent is NULL)
    // We never have to play with the root of the observation tree, because it has no parent
    while ((lineage != NULL) && (lineage->seen == 0)) {
      // put this ancestor on the stack to look at later
      stack[topStack] = lineage;
      topStack++;
      // Note that we already have seen this ancestor, for later lineage searches
      lineage->seen = 1;
      lineage = lineage->parent;  // Advance to this ancestor's parent
    }

    // Now trapse back down the stack, filling in each ancestor's info if need be
    while (topStack > 0) {
      topStack--;
      lineage = stack[topStack];
      // Try to fill in the holes of UNKNOWN. If the parent is also UNKNOWN, we know by construction 
      // that all of the other ancestors are also UNKNOWN, and thus the designation is correct
      if ((workingArray[lineage->ID] == -1) && (lineage->parent != NULL)) 
	workingArray[lineage->ID] = workingArray[lineage->parent->ID];
	if (workingArray[lineage->ID] == -1) 
	  flag = 0;
    }
  }

  if ((usage) && (flag)) 
    flagMap[x][y] = -2;
  else
    for (i=0; i < H_ID_NUMBER; i++) 
      observationArray[here][i] = workingArray[i];
}


//
// Finds the appropriate entry in the designated grid square, and then makes a duplicate of that entry
// modified according to the input.
//
void HighUpdateGridSquare(int x, int y, double distance, int hit, int parentID)
{
  TEntryList *tempEntry;
  int here, i;

  if (highMap[x][y] == NULL) {
    if (observationID >= AREA)
      fprintf(stderr, "dRoll over!\n");

    // Display ownership of this slot
    flagMap[x][y] = observationID;
    obsX[observationID] = x;
    obsY[observationID] = y;
    observationID++;

    highMap[x][y] = (TMapStarter *) malloc(sizeof(TMapStarter));
    if (highMap[x][y] == NULL) fprintf(stderr, "Malloc failed in creation of Map Starter at %d %d\n", x, y);
    highMap[x][y]->dead = 0;
    highMap[x][y]->total = 0;
    highMap[x][y]->size = 1;
    highMap[x][y]->array = (TMapNode *) malloc(sizeof(TMapNode));
    if (highMap[x][y]->array == NULL) fprintf(stderr, "Malloc failed in making initial map array for %d %d\n", x, y);

    // Initialize the slot
    for (i=0; i < H_ID_NUMBER; i++) 
      observationArray[flagMap[x][y]][i] = -1;
  }
  else if (flagMap[x][y] == 0) 
    HighBuildObservation(x, y, 0);

  here = observationArray[flagMap[x][y]][parentID];

  // If source->ID == parentID, just alter the source (no need for the obsArray)
  if ((here != -1) && (highMap[x][y]->array[here].ID == parentID)) {
    highMap[x][y]->array[here].hits = highMap[x][y]->array[here].hits + hit;
    highMap[x][y]->array[here].distance = highMap[x][y]->array[here].distance + distance;
  }
  // Use the lookup table to find the appropriate hits and distance 
  // already recorded by an ancestor, use source as parent. Dont forget 
  // to update the obsEntry.
  else {
    // We will be adding a new entry to the list- is there enough room?
    if (highMap[x][y]->size <= highMap[x][y]->total) {
      HighResizeArray(highMap[x][y], -71);
    }

    // Make all changes before incrementing highMap[x][y]->total, since it's used as an index

    // Update the obsArray
    observationArray[flagMap[x][y]][parentID] = highMap[x][y]->total;

    // Add an entry in to the list of altered map squares for this particle
    // First check to see if the size of that array is big enough to hold another entry
    if (h_particleID[parentID].size == 0) {
      h_particleID[parentID].size = 1;
      h_particleID[parentID].mapEntries = (TEntryList *) malloc(sizeof(TEntryList));
      if (h_particleID[parentID].mapEntries == NULL) fprintf(stderr, "Malloc failed in creation of entry list array\n");
    }
    else if (h_particleID[parentID].size <= h_particleID[parentID].total) {
      h_particleID[parentID].size = (int)(ceil(h_particleID[parentID].total*1.75));
      tempEntry = (TEntryList *) malloc(sizeof(TEntryList)*h_particleID[parentID].size);
      if (tempEntry == NULL) fprintf(stderr, "Malloc failed in expansion of entry list array\n");

      for (i=0; i < h_particleID[parentID].total; i++) {
	tempEntry[i].x = h_particleID[parentID].mapEntries[i].x;
	tempEntry[i].y = h_particleID[parentID].mapEntries[i].y;
	tempEntry[i].node = h_particleID[parentID].mapEntries[i].node;
      }

      free(h_particleID[parentID].mapEntries);
      h_particleID[parentID].mapEntries = tempEntry;
    }

    h_particleID[parentID].mapEntries[h_particleID[parentID].total].x = x;
    h_particleID[parentID].mapEntries[h_particleID[parentID].total].y = y;
    h_particleID[parentID].mapEntries[h_particleID[parentID].total].node = highMap[x][y]->total;

    i = highMap[x][y]->total;
    highMap[x][y]->array[i].source = h_particleID[parentID].total;
    highMap[x][y]->array[i].ID = parentID;
    h_particleID[parentID].total++;

    // Check to see if this square has been observed by an ancestor
    if (here == -1) {
      // Previously unknown; clean slate
      highMap[x][y]->array[i].hits = hit;
      highMap[x][y]->array[i].distance = distance + H_PRIOR_DIST;
      highMap[x][y]->array[i].parentGen = -2; 
    }
    else {
      // Include the pertinent info
      highMap[x][y]->array[i].hits = highMap[x][y]->array[here].hits + hit;
      highMap[x][y]->array[i].distance = distance + highMap[x][y]->array[here].distance;
      highMap[x][y]->array[i].parentGen = h_particleID[ highMap[x][y]->array[here].ID ].generation;
    }

    highMap[x][y]->total++;
  }
}



void HighDeleteObservation(short int x, short int y, short int node) {
  int total;

  // We already removed it in resizing the array
  if ((node == -1) || (highMap[x][y] == NULL))
    return;

  if (highMap[x][y]->total - highMap[x][y]->dead == 1) {
    free(highMap[x][y]->array);
    free(highMap[x][y]);
    highMap[x][y] = NULL;
    return;
  }

  // Look to shrink the array
  if ((int)((highMap[x][y]->total - 1 - highMap[x][y]->dead)*2.5) <= highMap[x][y]->size) {
    // Let resizing the array remove this entry
    HighResizeArray(highMap[x][y], highMap[x][y]->array[node].ID);
    if (highMap[x][y]->total == 0) {
      free(highMap[x][y]->array);
      free(highMap[x][y]);
      highMap[x][y] = NULL;
    }
    return;
  }

  highMap[x][y]->total--;
  total = highMap[x][y]->total;
  if (node != highMap[x][y]->total) {
    highMap[x][y]->array[node].hits      = highMap[x][y]->array[total].hits;
    highMap[x][y]->array[node].distance  = highMap[x][y]->array[total].distance;
    highMap[x][y]->array[node].ID        = highMap[x][y]->array[total].ID;
    highMap[x][y]->array[node].source    = highMap[x][y]->array[total].source;
    highMap[x][y]->array[node].parentGen = highMap[x][y]->array[total].parentGen;
    h_particleID[ highMap[x][y]->array[node].ID ].mapEntries[ highMap[x][y]->array[node].source ].node = node;
  }
}



//
// Input: x, y- location of a grid square
//        distance- the length of a line passing through the square
//        variance- an output, which will be filled the variance for errors in the lasercast that 
//                  should stop here.
// Output: returns the probability of trace of the given length through this square will be stopped by an obstacle
//
inline double HighComputeProbability(int x, int y, double distance, int parentID) 
{
  int here;

  if (highMap[x][y] == NULL) 
    return (1.0 - exp(H_PRIOR * distance));

  // If this grid square has been observed already this iteration, the last item in
  // the corresponding observation array will show ownership.
  // If that check fails, we have build the observation array entry for this square
  if (flagMap[x][y] == 0) 
    HighBuildObservation(x, y, 1);

  if (flagMap[x][y] == -2)
    return 0;

  here = observationArray[flagMap[x][y]][parentID];

  if (here == -1)
    return (1.0 - exp(H_PRIOR * distance));
  if (highMap[x][y]->array[here].hits == 0)
    return 0;
  return (1.0 - exp(-(highMap[x][y]->array[here].hits/highMap[x][y]->array[here].distance) * distance));
}



double HighComputeProb(int x, int y, double distance, int ID) 
{
  int i;

  if (highMap[x][y] == NULL) 
    return UNKNOWN;

  while (1) {
    for (i=0; i < highMap[x][y]->total; i++) {
      if (highMap[x][y]->array[i].ID == ID) {
	if (highMap[x][y]->array[i].hits == 0)
	  return 0;
	return (1.0 - exp(-(highMap[x][y]->array[i].hits/highMap[x][y]->array[i].distance) * distance));
      }
    }

    if (h_particleID[ID].parent == NULL)
      return UNKNOWN;
    else 
      ID = h_particleID[ID].parent->ID;
  }

  return UNKNOWN;
}




void HighAddTrace(double startx, double starty, double MeasuredDist, double theta, TAncestor *parent, int addEnd) {
  double overflow, slope; // Used for actually tracing the line
  int x, y, incX, incY, endx, endy;
  int xedge, yedge;       // Used in computing the midpoint. Recompensates for which edge of the square the line entered from
  double dx, dy;
  double distance, error, stdDist;
  double secant, cosecant;   // precomputed for speed

  secant = 1.0/fabs(cos(theta));
  cosecant = 1.0/fabs(sin(theta));

  distance = MIN(MeasuredDist, MAX_SENSE_RANGE);
  dx = (startx + (cos(theta) * distance));
  dy = (starty + (sin(theta) * distance));

  endx = (int) (dx);
  endy = (int) (dy);

  // Decide which x and y directions the line is travelling.
  if (startx > dx) {
    incX = -1;
    xedge = 1;
  }
  else {
    incX = 1;
    xedge = 0;
  }
  
  if (starty > dy) {
    incY = -1;
    yedge = 1;
  }
  else {
    incY = 1;
    yedge = 0;
  }

  if (fabs(startx - dx) > fabs(starty - dy)) {
    // The given starting point is non-integer. The line therefore starts at some point partially set in to the starting
    // square. Overflow starts at this offcenter amount, in order to make steps in the y direction at the right places.
    y = (int) (starty);
    overflow =  starty - y;
    if (incY == 1)
      overflow = 1.0 - overflow;
    slope = fabs(tan(theta));
    if (slope > 1.0) 
      slope = fabs((starty - dy) / (startx - dx));

    // The first square is a delicate thing, as we aren't doing a full square traversal in
    // either direction. So we figure out this strange portion of a step so that we can then
    // work off of the axes later. 
    // NOTE: we aren't computing the probability of this first step. Its a technical issue for
    // simplicity, and the odds of the sensor sitting on top of a solid object are sufficiently 
    // close to zero to ignore this tiny portion of a step. 
    error = fabs(((int)(startx)+incX+xedge)-startx);
    overflow = overflow - (slope*error);
    // The first step is actually in the y direction, due to the proximity of starty to the y axis. 
    if (overflow < 0.0) {
      y = y + incY;
      overflow = overflow + 1.0;
    }

    stdDist = slope*cosecant;

    for (x = (int) (startx) + incX; x != endx; x = x + incX) {
      overflow = overflow - slope;

      // Compute the distance travelled in this square
      if (overflow < 0.0)
	distance = (overflow+slope)*cosecant;
      else
	distance = stdDist;
      HighUpdateGridSquare(x, y, distance, 0, parent->ID);

      if (overflow < 0) {
	y = y + incY;
	distance = -overflow*cosecant;
	overflow = overflow + 1.0;
	HighUpdateGridSquare(x, y, distance, 0, parent->ID);
      }
    }

    if (addEnd) {
      if (incX < 0)
	distance = fabs((x+1) - dx)*secant;
      else
	distance = fabs(dx - x)*secant;
      HighUpdateGridSquare(endx, endy, distance, 1, parent->ID);
    }

  }

  else {
    x = (int) (startx);
    overflow = startx - x;
    if (incX == 1)
      overflow = 1.0 - overflow;
    slope = 1.0/fabs(tan(theta));

    // (See corresponding comments in the previous half of this function)
    error = fabs(((int)(starty)+incY+yedge)-starty);
    overflow = overflow - (error*slope);
    // The first step is actually in the y direction, due to the proximity of starty to the y axis. 
    if (overflow < 0.0) {
      x = x + incX;
      overflow = overflow + 1.0;
    }

    stdDist = slope*secant;

    for (y = (int) (starty) + incY; y != endy; y = y + incY) {
      overflow = overflow - slope;
      // Compute the distance travelled in this square. 
      if (overflow < 0)
	distance = (overflow+slope)*secant;
      else
	distance = stdDist;

      HighUpdateGridSquare(x, y, distance, 0, parent->ID);

      if (overflow < 0.0) {
	x = x + incX;
	distance = -overflow*secant;
	overflow = overflow + 1.0;
	HighUpdateGridSquare(x, y, distance, 0, parent->ID);
      }
    }

    if (addEnd) {
      if (incY < 0)
	distance = fabs(((y+1) - dy)/sin(theta));
      else
	distance = fabs((dy - y)/sin(theta));
      HighUpdateGridSquare(endx, endy, distance, 1, parent->ID);
    }
  }

}




//
// Inputs: x, y- starting point for the trace
//         theta- angle for the trace
//         measuredDist- the observed distance for this trace
//         parent- the most recent member of the ancestry for the particle being considered
//         hit- really an output, this will be filled with the total probability that this laser cast
//              hit an obstruction before reaching the maximum range of the sensor
// Output: The total evaluated probability for this laser cast (unnormalized). 
//
// Note that this trace automatically goes out to MAX_SENSE_RANGE, unless it is determined at some point 
// that any further trace has less than 0.01 probability of being reached, given the map.
//
double HighLineTrace(double startx, double starty, double theta, double MeasuredDist, int parentID) {
  double overflow, slope; // Used for actually tracing the line
  int x, y, incX, incY, endx, endy;
  double dx, dy;
  double totalProb; // Total probability that the line trace should have stopped before this step in the trace
  double eval;      // Total raw probability for the observation given this line trace through the map
  double prob, distance, error;
  double secant, cosecant;   // precomputed for speed
  double xblock, yblock;
  double xMotion, yMotion;
  double standardDist;

  eval = 0.0;
  totalProb = 1.0;
  secant = 1.0/fabs(cos(theta));
  cosecant = 1.0/fabs(sin(theta));

  distance = MIN(MeasuredDist+20.0, MAX_SENSE_RANGE);
  dx = (startx + (cos(theta) * distance));
  dy = (starty + (sin(theta) * distance));

  endx = (int) (dx);
  endy = (int) (dy);

  // Decide which x and y directions the line is travelling.
  if (startx > dx) {
    incX = -1;
    xblock = -startx;
  }
  else {
    incX = 1;
    xblock = 1.0-startx;
  }
  
  if (starty > dy) {
    incY = -1;
    yblock = -starty;
  }
  else {
    incY = 1;
    yblock = 1.0-starty;
  }
  
  if (fabs(startx - dx) > fabs(starty - dy)) {
    // The given starting point is non-integer. The line therefore starts at some point partially set in to the starting
    // square. Overflow starts at this offcenter amount, in order to make steps in the y direction at the right places.
    y = (int) (starty);
    overflow =  starty - y;
    if (incY == 1)
      overflow = 1.0 - overflow;
    slope = fabs(tan(theta));
    if (slope > 1.0) 
      slope = fabs((starty - dy) / (startx - dx));

    // The first square is a delicate thing, as we aren't doing a full square traversal in
    // either direction. So we figure out this strange portion of a step so that we can then
    // work off of the axes later. 
    // NOTE: we aren't computing the probability of this first step. Its a technical issue for
    // simplicity, and the odds of the sensor sitting on top of a solid object are sufficiently 
    // close to zero to ignore this tiny portion of a step. 
    dx = fabs((int)(startx)+xblock);
    dy = fabs(tan(theta)*dx);
    // The first step is actually in the y direction, due to the proximity of starty to the y axis. 
    if (overflow - dy < 0.0) {
      y = y + incY;
      overflow = (overflow - dy) + 1.0;
    }
    // Our first step is in fact in the x direction in this case. Set up for the overflow to 
    // be our starting offset plus this little extra we travel in the y direction.
    else 
      overflow = overflow - dy;

    standardDist = slope*cosecant;
    xMotion = -fabs(fabs(( ((int) (startx)) +xblock) * secant) - MeasuredDist);
    yMotion = -fabs(fabs((y+yblock) * cosecant) - MeasuredDist);

    for (x = (int) (startx) + incX; x != endx; x = x + incX) {
      overflow = overflow - slope;
      xMotion = xMotion + secant;
      if (overflow < 0.0) {
	error = fabs(yMotion);
	distance = (overflow+slope)*cosecant;
      }
      else {
	error = fabs(xMotion);
	distance = standardDist;
      }

      prob = totalProb * HighComputeProbability(x, y, distance, parentID);
      if (error < 20.0) 
	eval = eval + (prob * exp(-(error*error)/(2*HIGH_VARIANCE)));
      totalProb = totalProb - prob;
    
      if (overflow < 0.0) {
	y += incY;
	yMotion = yMotion + cosecant;
	error = fabs(xMotion);

	distance = -overflow*cosecant;
	overflow = overflow + 1.0;

	prob = totalProb * HighComputeProbability(x, y, distance, parentID);
	if (error < 20.0) 
	  eval = eval + (prob * exp(-(error*error)/(2*HIGH_VARIANCE)));
	totalProb = totalProb - prob;
      }
    
    }
  }

  else {
    x = (int) (startx);
    overflow = startx - x;
    if (incX == 1)
      overflow = 1.0 - overflow;
    slope = 1.0/fabs(tan(theta));

    // (See corresponding comments in the previous half of this function)
    dy = fabs((int)(starty)+yblock);
    dx = fabs(dy/tan(theta));
    // The first step is actually in the y direction, due to the proximity of starty to the y axis. 
    if (overflow - dx < 0) {
      x = x + incX;
      overflow = (overflow - dx) + 1.0;
    }
    // Our first step is in fact in the x direction in this case. Set up for the overflow to 
    // be our starting offset plus this little extra we travel in the y direction.
    else
      overflow = overflow - dx;

    standardDist = slope*secant;
    xMotion = -fabs(fabs((x+xblock) * secant) - MeasuredDist);
    yMotion = -fabs(fabs(( ((int) (starty)) +yblock) * cosecant) - MeasuredDist);

    for (y = (int) (starty) + incY; y != endy; y = y + incY) {
      yMotion = yMotion + cosecant;
      overflow = overflow - slope;

      if (overflow < 0.0) { 
	error = fabs(xMotion);
	distance = (overflow+slope)*secant;
      }
      else {
	error = fabs(yMotion);
	distance = standardDist;
      }

      prob = totalProb * HighComputeProbability(x, y, distance, parentID);
      if (error < 20.0) 
	eval = eval + (prob * exp(-(error*error)/(2*HIGH_VARIANCE)));
      totalProb = totalProb - prob;
    
      if (overflow < 0.0) {
	x += incX;
	xMotion = xMotion + secant;
	error = fabs(yMotion);

	distance = -overflow*secant;
	overflow = overflow + 1.0;

	prob = totalProb * HighComputeProbability(x, y, distance, parentID);
	if (error < 20.0) 
	  eval = eval + (prob * exp(-(error*error)/(2*HIGH_VARIANCE)));
	totalProb = totalProb - prob;
      }

    }
  }

  if (MeasuredDist >= MAX_SENSE_RANGE) 
    return (eval + totalProb);
  if (totalProb == 1)
    return 0;
  return (eval / (1.0 - totalProb));
}

