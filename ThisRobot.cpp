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
// Copyright 2005 Austin Eliazar, Ronald Parr, Duke University
//

/*****************************************************
 *           PLAYBACK ROBOT
 * This is basically an empty file for
 * allowing the code to run offline,
 * with all of the direct calls to the robot 
 * doing nothing. This is to be used exclusively
 * in playback mode, and any attempts to do 
 * otherwise will probably give very boring results.
 *****************************************************/


// All base robot files are required to have this exact same definition. DO NOT CHANGE!
// This is the structure in which the SICK laser data is stored, so that the SLAM program
// can use it.

#include "ThisRobot.h"


TOdo odometry;


int InitializeThisRobot(int argc, char *argv[]) {
  return 0;
}

int ConnectOdometry(int argc, char *argv[]) {
  return 0;
}

int ConnectLaser(int argc, char *argv[]) {
  return 0;
}

int ConnectDrive(int argc, char *argv[]) {
  return 0;
}

// This is never going to be called if we are truly running this offline. If for some reason this is run
// not in playback mode, this function wont do anything.
void GetSensation(TSense &sense) {
  return;
}


void GetOdometry(TOdo &odometry) {
  return;
}


void Drive(double speed, double turn) {
  return;
}
