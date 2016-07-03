#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "probability.h"
#include "perModel.h"
#include "engine.h"
#include "matrix.h"
#include "tmwtypes.h"
#include "macro.h"
#include "setting.h"
#include "nodeInfo.h"
#include "Initialization.h"

extern Engine *gEp;
extern double r[(NUM_STA+1)*(NUM_STA+1)];
extern double pro[NUM_STA+1][NUM_STA+1];
extern double dummyA[NUM_STA*2][(NUM_STA+1)*(NUM_STA+1)];
extern double A[NUM_STA*2][(NUM_STA+1)*(NUM_STA+1)];
extern double u[NUM_STA*2];
extern double dummyAeq[2][(NUM_STA+1)*(NUM_STA+1)];
extern double Aeq[2][(NUM_STA+1)*(NUM_STA+1)];
extern double beq[2];
extern double lb[(NUM_STA+1)*(NUM_STA+1)];

void solveLP(){
	int i, j;
	int tate = NUM_STA * 2;
	int yoko = pow(NUM_STA+1, 2);
	//char buffer[EP_BUFFER_SIZE] = {'\0'};

	printf("Setting matrixes.\n");

	mxArray *mx_p = NULL;
	mxArray *mx_fval = NULL;
	mxArray *mx_r = NULL;
	mxArray *mx_A = NULL;
	mxArray *mx_u = NULL;
	mxArray *mx_Aeq = NULL;
	mxArray *mx_beq = NULL;
	mxArray *mx_lb = NULL;
	mx_p = mxCreateDoubleMatrix(1, yoko, mxREAL);
	mx_fval = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *p, *fval;

	mx_r = mxCreateDoubleMatrix(1, yoko, mxREAL);
	memcpy((void *)mxGetPr(mx_r), (void *)r, sizeof(r));
	mx_A = mxCreateDoubleMatrix(tate, yoko, mxREAL);
	memcpy((void *)mxGetPr(mx_A), (void *)A, sizeof(A));
	mx_u = mxCreateDoubleMatrix(1, tate, mxREAL);
	memcpy((void *)mxGetPr(mx_u), (void *)u, sizeof(u));
	mx_Aeq = mxCreateDoubleMatrix(2, yoko, mxREAL);
	memcpy((void *)mxGetPr(mx_Aeq), (void *)Aeq, sizeof(Aeq));
	mx_beq = mxCreateDoubleMatrix(1, 2, mxREAL);
	memcpy((void *)mxGetPr(mx_beq), (void *)beq, sizeof(beq));
	mx_lb = mxCreateDoubleMatrix(1, yoko, mxREAL);
	memcpy((void *)mxGetPr(mx_lb), (void *)lb, sizeof(lb));

	engPutVariable(gEp, "mx_r", mx_r);
	engPutVariable(gEp, "mx_A", mx_A);
	engPutVariable(gEp, "mx_u", mx_u);
	engPutVariable(gEp, "mx_Aeq", mx_Aeq);
	engPutVariable(gEp, "mx_beq", mx_beq);
	engPutVariable(gEp, "mx_lb", mx_lb);

	//engOutputBuffer(gEp, buffer, EP_BUFFER_SIZE);
	//engEvalString(gEp, "mx_r");
	//printf("%s", buffer);

	printf("Optimization starts.\n");

	engEvalString(gEp, "[p, fval] = linprog(mx_r, mx_A, mx_u, mx_Aeq, mx_beq, mx_lb, []);");
	//printf("%s", buffer);
	engEvalString(gEp, "p = p ./ 100;");
	engEvalString(gEp, "fval = fval / (-100);");
	//printf("%s", buffer);
	mx_p = engGetVariable(gEp, "p");
	p = mxGetPr(mx_p);
	mx_fval = engGetVariable(gEp, "fval");
	fval = mxGetPr(mx_fval);

	for(i=0; i<yoko; i++){
		if(p[i]>=0.000001){
			pro[i/(NUM_STA+1)][i%(NUM_STA+1)] = p[i];
			//printf("%f,", p[i]);
		}else{
			pro[i/(NUM_STA+1)][i%(NUM_STA+1)] = 0;
		}
		//printf("%f, ", p[i]);
	}
	//printf("\n\n");
	printf("Optimization terminated.\n");
	/*printf("***** Probability *****\n");
	for(i=0; i<=NUM_STA; i++){
		for(j=0; j<=NUM_STA; j++){
			printf("%f,", pro[i][j]);
		}
		printf("\n");
	}

	for(i=0; i<yoko; i++){
		if(p[i]>0.00001){
			//printf("\n   p[%d] = %f\n", i, p[i]);
		}
	}
	printf("   fval = %f\n", *fval);
	printf("***** Probability *****\n\n ");*/

	mxDestroyArray(mx_r);
	mxDestroyArray(mx_A);
	mxDestroyArray(mx_u);
	mxDestroyArray(mx_Aeq);
	mxDestroyArray(mx_beq);
	mxDestroyArray(mx_lb);
	mxDestroyArray(mx_p);
	mxDestroyArray(mx_fval);
	//engEvalString(gEp, "close;");
}

void calculateProbability(staInfo sta[], apInfo *ap, int mode){
	int nodeID;
	int i;
	double delay[NUM_STA+1] = {};
	if(mode == 0){   //random
		nodeID = rand() % gSpec.numSta;
		for(i=0; i<gSpec.numSta; i++){
			if(i==nodeID){
				sta[i].fTx = true;
			}else{
				sta[i].fTx = false;
			}
		}
	}else{   //probability
		//calculateDelay
		if(gSpec.proMode==1||gSpec.proMode==2){
			calculateDelay(ap, sta, delay);
		}
		calculateRSSI(ap, sta, delay);
		solveLP();
		//selectNode(sta, fUpColl, fNoUplink);
	}
}

void initializeMatrix(){
	int tate = NUM_STA * 2;
	int yoko = pow((NUM_STA+1), 2);
	int i, j, no;

	for(i=0; i<NUM_STA; i++){
		for(j=0; j<yoko; j++){
			if((j/(NUM_STA+1)==i+1)&&(j%(NUM_STA+1)!=i+1)){
				dummyA[i][j] = -1;
			}else{
				dummyA[i][j] = 0;
			}
			//printf("%f ", dummyA[i][j]);
		}
		//printf("\n");
	}
	for(i=NUM_STA; i<NUM_STA*2; i++){
		for(j=0; j<yoko; j++){
			if((j%(NUM_STA+1)==(i-NUM_STA+1))&&(j/(NUM_STA+1)!=(i-NUM_STA+1))){
				dummyA[i][j] = -1;
			}else{
				dummyA[i][j] = 0;
			}
			//printf("%f ", dummyA[i][j]);
		}
		//printf("\n");
	}
	for(j=0; j<yoko; j++){
		for(i=0; i<tate; i++){
			no = tate * j + i + 1;
			if(no%yoko!=0){
				A[no/yoko][no%yoko-1] = dummyA[i][j];
			}else{
				A[no/yoko-1][yoko-1] = dummyA[i][j];
			}
		}
	}

	tate = 2;

	for(j=0; j<yoko; j++){
		if(j/(NUM_STA+1)==j%(NUM_STA+1)){
			dummyAeq[1][j] = 1;
		}else{
			dummyAeq[0][j] = 1;
		}
	}
	for(j=0; j<yoko; j++){
		for(i=0; i<tate; i++){
			no = tate * j + i + 1;
			if(no%yoko!=0){
				Aeq[no/yoko][no%yoko-1] = dummyAeq[i][j];
			}else{
				Aeq[no/yoko-1][yoko-1] = dummyAeq[i][j];
			}
		}
	}
	for(i=0; i<NUM_STA*2; i++){
		u[i] = -100/(2*NUM_STA);
	}
}

int selectNode(staInfo sta[], bool *fUpColl, bool *fNoUplink, bool *fNoDownlink, int *upNode, int *downNode){
	double *proDown;
	proDown = (double*)malloc(sizeof(double)*(NUM_STA+1));// = {};
	double *proUp;
	proUp = (double*)malloc(sizeof(double)*(NUM_STA+1));// = {};
	double *proTempDown;
	proTempDown = (double*)malloc(sizeof(double)*(NUM_STA+1));// = {};   //確率判定のため
	//int upNode, downNode;   //0がなし．1--NUM_STAまで．配列とずれてるから注意．
	double downRand;
	int i, j;
	int numTx = 0;
	int minBackoff = INT_MAX;
	int dummyNode = INT_MAX;

	//配列の初期化
	initializeDoubleArray(proDown, NUM_STA+1, 0);
	initializeDoubleArray(proUp, NUM_STA+1, 0);
	initializeDoubleArray(proTempDown, NUM_STA+1, 0);

	//下り通信を受信する端末の決定
	//printf("***** Probability that each node is selected as a destination node of AP. *****\n");
	for(i=0; i<NUM_STA+1; i++){
		if(i!=0){
			proTempDown[i] += proTempDown[i-1];
		}
		for(j=0; j<NUM_STA+1; j++){
			proDown[i] += pro[i][j];
		}
		proTempDown[i] += proDown[i];
		//printf("p_d[%d] is %f.\n", i, proTempDown[i]);
	}

	if(proTempDown[NUM_STA]<=0.999 || 1.001<=proTempDown[NUM_STA]){
		printf("Probability is wrong.%f\n", proTempDown[NUM_STA]);
	}

	downRand = (double)rand() / RAND_MAX;
	//printf("downRand is %f\n", downRand);
	for(i=0; i<=NUM_STA; i++){
		if(i==0){
			if(downRand<=proTempDown[i]){
				*downNode = i;
				//printf("Dummy STA is selected as a destination node.\n");
				*fNoDownlink = true;
				break;
			}
		}else if(proTempDown[i-1]<downRand && downRand<=proTempDown[i]){
			*downNode = i;
			sta[i-1].fRx = true;
			//printf("STA %d is selected as a destination node.\n", i-1);
			break;
		}
		if(i==NUM_STA){
			printf("Error. downRand = %f\n", downRand);
			*downNode = 0;
			*fNoDownlink = true;
		}
	}

	//上り通信端末の選択
	//printf("***** Probabiliy that each node is selected as a source node of AP. *****\n");
	for(j=0; j<NUM_STA+1; j++){
		if(*downNode==j){
			proUp[j] = 0;
		}else{
			proUp[j] = pro[*downNode][j]/proDown[*downNode];
			//if(proUp[j]<0.000001){
				//printf("%f, %f ", pro[*downNode][j], proDown[*downNode]);
			//}*/
		}
		//printf("p_u[%d] is %f\n", j, proUp[j]);
	}
	//printf("\n\n");
	//int temp = 0;
	for(i=0; i<NUM_STA+1; i++){
		if(proUp[i]!=0){
			if(i==0){
				dummyNode = rand() % ((int)(1/proUp[i])+1);
			}else{
				sta[i-1].cw = (int)(1/proUp[i]);
				sta[i-1].backoffCount = rand() % (sta[i-1].cw+1);
				//printf("%f, %d ", proUp[i], sta[i-1].backoffCount);
				//temp++;
			}
		}
		if(i!=0){
			//printf("%d, ", sta[i-1].cw);
		}
	}
	//printf("%d, ", temp);
	//printf("\n\n");

	for(i=0; i<gSpec.numSta; i++){
		if(proUp[i+1]!=0){
			if((minBackoff>sta[i].backoffCount)&&(sta[i].buffer[0].lengthMsdu!=0)){
				minBackoff = sta[i].backoffCount;
			}
		}
	}
	//printf("%d, ", minBackoff);
	if(minBackoff==INT_MAX){
		//printf("All STAs don't have a frame.\n");   //フレームが無いときだけじゃないかも ダミーが選ばれる場合も
	}
	if(dummyNode<minBackoff){
		minBackoff = dummyNode;
		*fNoUplink = true;
		*upNode = 0;
	}else{
		if(*fNoUplink==false){
			for(i=0; i<gSpec.numSta; i++){
				if(proUp[i+1]!=0 && minBackoff==sta[i].backoffCount && sta[i].fRx==false){
					sta[i].fTx = true;
					//sta[i].backoffCount = rand() % (sta[i].cw + 1);
					numTx++;
					*upNode = i+1;
					//printf("STA %d has minimum backoff count.\n", i);
				}else{
					//sta[i].backoffCount -= minBackoff;
					sta[i].fTx = false;
				}
			}
		}else{
			sta[i].fTx = false;
		}
	}
	//printf("%d, ", numTx);
	if(numTx==0 && *fNoUplink==false){
		printf("undefined\n");
	}else if(numTx==1){
		*fUpColl = false;
	}else{
		*fUpColl = true;
	}
	//printf("(%d, %d),", *downNode, *upNode);

	free(proUp);
	free(proDown);
	free(proTempDown);

	return minBackoff;
}
