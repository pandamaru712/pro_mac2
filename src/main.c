#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include "nodeInfo.h"
#include "setting.h"
#include "initialization.h"
#include "success.h"
#include "idle.h"
#include "result.h"
#include "macro.h"
#include "probability.h"
#include "success.h"
#include "perModel.h"

#include "engine.h"
#include "matrix.h"
#include "tmwtypes.h"

double gElapsedTime;
std11 gStd;
simSpec gSpec;
FILE *gFileSta;
FILE *gFp;

Engine *gEp;
double r[(NUM_STA+1)*(NUM_STA+1)] = {};//{-1, -4, -45, -51, -29, -42, -16, -1, -25, -39, -24, -35, -3, -23, -1, -56, -78, -10, -11, -34, -22, -1, -7, -67, -45, -23, -65, -55, -1, -76, -12, -6, -95, -67, -52, -1};
double pro[NUM_STA+1][NUM_STA+1];
double dummyA[NUM_STA*2][(NUM_STA+1)*(NUM_STA+1)];
double A[NUM_STA*2][(NUM_STA+1)*(NUM_STA+1)];
double u[NUM_STA*2];
double dummyAeq[2][(NUM_STA+1)*(NUM_STA+1)];
double Aeq[2][(NUM_STA+1)*(NUM_STA+1)];
double beq[2] = {100, 0};
double lb[(NUM_STA+1)*(NUM_STA+1)] = {};

//void simSetting(int, char**);

int main(int argc, char *argv[]){
	//Check option values from command line.
	//checkOption(argc, argv);
	//Apply option values to simulation settings.
	printf("Start simulation.\n");

	simSetting(argc,argv);
	printf("Set simulation parameters.\n");
	/*if((gFileSta=fopen("sta's buffer.txt", "w"))==NULL){
		printf("File cannot open! 3");
		exit(33);
	}*/

	staInfo *sta;
	sta = (staInfo *)malloc(sizeof(staInfo)*gSpec.numSta);

	//sta = (staInfo*)malloc(sizeof(staInfo)*gSpec.numSta);
	apInfo ap;
	resultInfo result;
	//Intialize result information.
	initializeResult(&result);

	int numTx = 0;
	int trialID;
	bool fEmpty = false;
	double lastBeacon = 0;

	gFp = fopen("topology.txt", "r");
	if(gFp==NULL){
		printf("Can not open topology file.\n");
		exit(92);
	}

	if(!(gEp = engOpen(""))){
		fprintf(stderr, "\nCan't start MATLAB engine.\n");
		return EXIT_FAILURE;
	}
	printf("Open MATLAB engine.\n");

	for (trialID=0; trialID<gSpec.numTrial; trialID++){
		printf("\n***** No. %d *****\n", trialID);
		srand(trialID);
		numTx = 0;
		fEmpty = false;
		lastBeacon = 0;
		initializeNodeInfo(sta, &ap);
		gElapsedTime = gStd.difs;
		initializeMatrix();
		printf("Initialization NodeInfo and Matrix.\n");
		calculateProbability(sta, &ap, 1);

		for( ;gElapsedTime<gSpec.simTime*1000000; ){
			transmission(sta, &ap);

			if(lastBeacon+100000<gElapsedTime){
				if(gSpec.proMode==1 || gSpec.proMode==2){
					calculateProbability(sta, &ap, 1);
				}
				lastBeacon = gElapsedTime;
			}
		}

		simulationResult(sta, &ap, &result, trialID);
	}

	if(gSpec.fOutput==true){
		fclose(gSpec.output);
	}
	//fclose(gFileSta);
	free(sta);
	printf("Free memory space.\n");

	engEvalString(gEp, "close;");
	engClose(gEp);
	fclose(gFp);
	printf("Close MATLAB.\nFinish.\n");
	return 0;
}
