// These includes pull in interface definitions and utilities.
#include "ThisRobot.cpp"


// Set the maximum usuable distance for the laser range finder. This number is often less than the actual
// reliable distance for the specific LRF, because the laser casts 'scatter' at long distances.
#define MAX_SENSE_RANGE 7.95 * MAP_SCALE

// Some useful macros
#define SIGN(A) ((A) >= 0.0 ? (1.0) : (-1.0))
#define SQUARE(A) (((A) * (A)))
#define MAX(A,B) ((A) >= (B) ? (A) : (B))
#define MIN(A,B) ((A) >= (B) ? (B) : (A))



