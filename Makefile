#
# Makefile for SLAM under RedHat Linux
#
# Copyright 2002, Austin Eliazar, Ronald Parr & Duke University
#
#
# Tested with:
#
# - Redhat 8.0
#
# Need to modify to include:
#
# - Generalizing and testing on other platforms
# - Make for pure localization
#


CC = g++
#CFLAGS = -g -I. -I/usr/local/gcc280/lib/g++-include
CFLAGS = -g -I. -O3 -Wall -I/usr/local/lib/g++-include
# For profiling the code : 
#CFLAGS += -pg

#LDFLAGS =  -lnsl -lnls -lsocket
LDFLAGS = -lpthread

SRC = mt-rand.o ThisRobot.o basic.o map.o lowMap.o low.o highMap.o high.o slam.o

slam : $(SRC)
	$(CC) $(CFLAGS) -o slam $(SRC) $(LDFLAGS)

slam.o : slam.cpp high.h
	$(CC) $(CFLAGS) -c slam.cpp

high.o : high.c high.h highMap.h
	$(CC) $(CFLAGS) -c high.c

highMap.o : highMap.c highMap.h low.h 
	$(CC) $(CFLAGS) -c highMap.c

low.o : low.c low.h lowMap.h
	$(CC) $(CFLAGS) -c low.c

lowMap.o : lowMap.c lowMap.h map.h
	$(CC) $(CFLAGS) -c lowMap.c

mt-rand.o : mt-rand.c mt-rand.h
	$(CC) $(CFLAGS) -c mt-rand.c

ThisRobot.o : ThisRobot.h basic.h
	$(CC) $(CFLAGS) -c ThisRobot.cpp

basic.o : basic.h
	$(CC) $(CFLAGS) -c basic.c

map.o : laser.h map.h map.c
	$(CC) $(CFLAGS) -c map.c


