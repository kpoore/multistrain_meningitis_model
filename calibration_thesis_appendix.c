//Libraries required for the program
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>

//Define constant values
int POP_SIZE = 1.7e6;
int ARRAY_SIZE = 2200000;
int SIM_LENGTH = 12*500;
double CONTACT_MATRIX[101][101];
double TRANS_PROB_A[101][101], TRANS_PROB_B [101][101];
double TRANS_PROB_C [101][101], TRANS_PROB_NL [101][101];
double DEATH_RATE [100];
int SEROORDER [7];
int BIRTHS = 1300;
int NUM_PROGRAMS = 0;

//Undefined constant values
char file_name [101], fullPath [150];
double FINAL_CALIBRATION_TIME, RHO_MONO, RHO_QUAD;
double EPSILON_A_1, EPSILON_A_2, EPSILON_A_3, EPSILON_A_4, EPSILON_A_5;
double EPSILON_B_1, EPSILON_B_2, EPSILON_B_3, EPSILON_B_4, EPSILON_B_5;
double EPSILON_C_1, EPSILON_C_2, EPSILON_C_3, EPSILON_C_4, EPSILON_C_5;
double BETA_A_1, BETA_A_2, BETA_A_3, BETA_A_4, BETA_A_5;
double BETA_B_1, BETA_B_2, BETA_B_3, BETA_B_4, BETA_B_5;
double BETA_C_1, BETA_C_2, BETA_C_3, BETA_C_4, BETA_C_5;
double BETA_NL_1, BETA_NL_2, BETA_NL_3, BETA_NL_4, BETA_NL_5;
double AVG_DUR_INFECT_A, AVG_DUR_INFECT_B, AVG_DUR_INFECT_C, AVG_DUR_INFECT_NL;
double AVG_DUR_IMMUNE_A, AVG_DUR_IMMUNE_B, AVG_DUR_IMMUNE_C, AVG_DUR_IMMUNE_NL;
double CROSS_A, CROSS_B, CROSS_C, CROSS_NL;
double ACCEPT_PARAMS[400][75];


//Structure for each individual
typedef struct pop{//create agent.
    int status; //0 susc, 1 infectious, 2 recovery, 3 IMD. 
	int vaccine;//binary 1 for vaccineated, 0 for not-vaccinated
    char vaccineType;
    float vaccineEfficacy;
    float infectTime;
    float recoverTime;
    float vaccineTime;
    float age;
    char serogroup;
    char immuneSero;
    float immuneEfficacy;
}pop;

//Structure for each vaccine program
typedef struct vaccineProgram{
    char serogroupProtection;
    float coverage;
    float efficacy;
    float sTime;
    float fTime;
    float lengthVaccinated;
    int ageApplied;
}vaccineProgram;

//list of functions
int SIR(int flagOpt ,double *targets,  char *filename);
float phaseTimer (double average, double sigma);
char probSerogroup(pop *population);
double coinFlip ();
int ranCheck (double prob, double ranNum);
float gammaDist (double a, double b);
void infection(pop *population, float curTime);
void applyVaccine(pop *agent, vaccineProgram *vaccine, double currTime);
void orderSero (int numA, int numB, int numC, int numNl, int numInfect);
int compare(const double* a, const double* b);
double runningAVG (int size, double *a);
void shuffle(int *index);

//main program
int
main(int argc, char *argv[])
{
    int i, numAccept = 0, optFlag, itemp;
    double flag1_params [147], flag2_params [139], flag3_params[21];
    char *filenames;
    FILE *flag1_beta, *flag1_epsilon, *flag3_accepted;
    optFlag = 1;
    i = 0;
	filenames = argv[1];
	//parameters for calibration
    flag1_params[0] = 312000;//Length of time for Calibration
    flag1_params[1] = 0.0000007;//stepsize
    flag1_params[2] = 0.1; //tolerance A
    flag1_params[3] = 0.1; //tolerance B
    flag1_params[4] = 0.1; //tolerance C
    flag1_params[5] = 0.1; //tolerance W
    flag1_params[6] = 0.1; //tolerance X
    flag1_params[7] = 0.1; //tolerance Y
    flag1_params[8] = 0.29697; //IMD target, A
    flag1_params[9] = 0.057986;
    flag1_params[10] = 0.011992;
    flag1_params[11] = 0.159819;
    flag1_params[12] = 0.092761;
    flag1_params[13] = 1.712876; //IMD target, B
    flag1_params[14] = 0.422106;
    flag1_params[15] = 0.095002;
    flag1_params[16] = 2.652981;
    flag1_params[17] = 0.253702;
    flag1_params[18] = 2.051521; //IMD target, C
    flag1_params[19] = 0.194890;
    flag1_params[20] = 0.221264;
    flag1_params[21] = 1.534278;
    flag1_params[22] = 0.396670;
    flag1_params[23] = 0.016756; //IMD target, W
    flag1_params[24] = 0.075799;
    flag1_params[25] = 0.129617;
    flag1_params[26] = 0.368358;
    flag1_params[27] = 0.059707;
    flag1_params[28] = 0.0000; //IMD target, X
    flag1_params[29] = 0.0000;
    flag1_params[30] = 0.0000;
    flag1_params[31] = 0.0000;
    flag1_params[32] = 0.0000;
    flag1_params[33] = 0.352212; //IMD target, Y
    flag1_params[34] = 0.002644;
    flag1_params[35] = 0.11314;
    flag1_params[36] = 0.315598;
    flag1_params[37] = 0.44419;
	
    flag1_params[38] = 0.1; //tolerance, A carriage
    flag1_params[39] = 0.1; //tolerance, B carriage
    flag1_params[40] = 0.1; //tolerance, C carriage
    flag1_params[41] = 0.1; //tolerance, W carriage
    flag1_params[42] = 0.1; //tolerance, X carriage
    flag1_params[43] = 0.1; //tolerance, Y carriage
    flag1_params[44] = 0.1; //tolerance, Nl carriage
    flag1_params[45] = 0.027142; //prevalence target, A carriage
    flag1_params[46] = 0.03393;
    flag1_params[47] = 0.040725;
    flag1_params[48] = 0.0723964;
    flag1_params[49] = 0.047964;
    flag1_params[50] = 0.01530; //prevalence target, B carriage
    flag1_params[51] = 0.01913;
    flag1_params[52] = 0.02295;
    flag1_params[53] = 0.04080;
    flag1_params[54] = 0.02703;
    flag1_params[55] = 0.00862; //prevalence target, C carriage
    flag1_params[56] = 0.010775;
    flag1_params[57] = 0.012935;
    flag1_params[58] = 0.022995;
    flag1_params[59] = 0.015235;
    flag1_params[60] = 0.11444; //prevalence target, Nl carriage
    flag1_params[61] = 0.1100;
    flag1_params[62] = 0.075000;
    flag1_params[63] = 0.0400;
    flag1_params[64] = 0.021660;
    flag1_params[65] = 0.000;//not used anywhere
    flag1_params[66] = 0.00001; //prev_stepsize betas are on the order of 1e-5 -- 1e-6
    flag1_params[67] = 1.0;//tolerance Prev A_1
    flag1_params[68] = 1.0;//tolerance Prev A_2
    flag1_params[69] = 1.0;//tolerance Prev A_3
    flag1_params[70] = 0.9;//tolerance Prev A_4
    flag1_params[71] = 0.20;//tolerance Prev A_5
    flag1_params[72] = 0.55;//tolerance Prev B_1
    flag1_params[73] = 0.5;//tolerance Prev B_2
    flag1_params[74] = 0.5;//tolerance Prev B_3
    flag1_params[75] = 0.33;//tolerance Prev B_4
    flag1_params[76] = 0.1;//tolerance Prev B_5
    flag1_params[77] = 0.999;//tolerance Prev C_1
    flag1_params[78] = 0.999;//tolerance Prev C_2
    flag1_params[79] = 0.999;//tolerance Prev C_3
    flag1_params[80] = 0.999;//tolerance Prev C_4
    flag1_params[81] = 0.999;//tolerance Prev C_5
    flag1_params[82] = 0.2;//tolerance Prev Nl_1
    flag1_params[83] = 0.33;//tolerance Prev Nl_2
    flag1_params[84] = 0.3;//tolerance Prev Nl_3
    flag1_params[85] = 0.33;//tolerance Prev NL_4
    flag1_params[86] = 0.5;//tolerance Prev NL_5
    flag1_params[90] = 0.50;//tolerance IMD A_1
    flag1_params[91] = 2.0;//tolerance IMD A_2
    flag1_params[92] = 1.8;//tolerance IMD A_3
    flag1_params[93] = 1.0;//tolerance IMD A_4
    flag1_params[94] = 1.5;//tolerance IMD A_5
    flag1_params[95] = 0.5;//tolerance IMD B_1
    flag1_params[96] = 0.8;//tolerance IMD B_2
    flag1_params[97] = 1.2;//tolerance IMD B_3
    flag1_params[98] = 0.5;//tolerance IMD B_4
    flag1_params[99] = 1.0;//tolerance IMD B_5
    flag1_params[100] = 0.5;//tolerance IMD C_1
    flag1_params[101] = 0.8;//tolerance IMD C_2
    flag1_params[102] = 0.8;//tolerance IMD C_3
    flag1_params[103] = 0.5;//tolerance IMD C_4
    flag1_params[104] = 0.8;//tolerance IMD C_5
	
    FINAL_CALIBRATION_TIME = flag1_params[0];
	
    if(optFlag == 1){
		//Variables for calibration
        SIM_LENGTH = flag1_params[0];
        
        BETA_A_1 = 0.001000000000; BETA_A_2 = 0.000700000000;
		BETA_A_3 = 0.000295000000; BETA_A_4 = 0.001050000000;
		BETA_A_5 = 0.001162500000;
        BETA_B_1 = 0.001079892600; BETA_B_2 = 0.000515991200;
		BETA_B_3 = 0.000643461300; BETA_B_4 = 0.000795039000;
		BETA_B_5 = 0.001084655800;
        BETA_C_1 = 0.001065819200; BETA_C_2 = 0.000628453100;
		BETA_C_3 = 0.000657328100; BETA_C_4 = 0.000935981500; BETA_C_5 = 0.001029868200;
        BETA_NL_1= 0.039575000000; BETA_NL_2= 0.185629819500;
		BETA_NL_3= 0.018825000000; BETA_NL_4= 0.005337500000; BETA_NL_5= 0.002015645800;
		
        CROSS_A = 0.2; CROSS_B = 0.2;
		CROSS_C = 0.2; CROSS_NL = 0.2;
        AVG_DUR_IMMUNE_A = 24.0;	AVG_DUR_IMMUNE_B = 24.0;
		AVG_DUR_IMMUNE_C = 24.0;	AVG_DUR_IMMUNE_NL = 56.4;
        AVG_DUR_INFECT_A = 13.0;	AVG_DUR_INFECT_B = 13.0;
		AVG_DUR_INFECT_C = 13.0;	AVG_DUR_INFECT_NL = 4.8;
		
        itemp = SIR(optFlag, flag1_params, filenames);
	}
    return 1;
}


int
SIR(int flagOpt ,double *target, char *filenames){
	//Declaration of variables
    int numSusc = 0, numInfect= 0, numIMD = 0, numNl = 0, numRecover = 0;
    int numSeroA = 0, numSeroB = 0, numSeroC = 0;
	int vaccineStart = 0, totalPasses = 0, numVac = 0;
    int initInfect = 9000;
    int infectByAge [101] ={0};
	int condenseAge[5], condenseSero[4][5];
	int numAgentSusc[101][4], numAgentImmune[101][4];
    int infectSeroA [101]={0}, infectSeroB [101]={0};
	int infectSeroC [101]={0}, infectSeroNl [101]={0};
    int IMDPerYear = 0;
    int i, j, k, itemp, agetemp;
	int numAgentPop [101] ={0};
	int *randInfectCtr;
	int firstVac, vaccineIndex[4200]={-1}, signPrev[20]={1}, signIMD[15]={1};
    double pS2IA[101], pS2IB[101], pS2IC[101], pS2INl[101];
    double prevAVG[4][5], imdAVG[3][5], stepsizePrev[4][5], stepsizeIMD[3][5]
	double avgDataPrev[4][5][3901], avgDataIMD[3][5][3901];
    double avgDataPrev3[4][5][3901], avgDataIMD3[3][5][3901];
    double IMD100KPop [3][5], imdCtr[3][5];
    double randNumInfect;
    double condensePrev[4][5], condenseIMD[3][5];
    double infectBySero [4] = {0};
    double pR2IA[3][101], pR2IB[3][101], pR2IC[3][101], pR2INl[3][101];
    double t = 0.0, dt = 1.0, dtemp = 0.0, vaccineLoss = 0.00481;
    unsigned int iseed;
    char serogroup = ' ', ctemp =' ';
    FILE *fout, *fage, *fdemograph, *fsero, *fAVGsero, *fcontact;
	FILE *fdeathRate, *fvaccineProg, *flag1_beta, *flag1_stepsize;
	FILE *flag1_epsilon, *fconditions;
    FILE *fAVGIMD, *fConDemo;
    FILE *fprevA, *fprevB, *fprevC, *fprevNl;
    FILE *Beta_B, *Epsilon_B, *infectionB, *prevalenceB;
    vaccineProgram *monoval, *monovalTmp;
    pop *population;
    
    randInfectCtr = (int *)calloc(ARRAY_SIZE, sizeof(int));
    SIM_LENGTH = target[0];
	
	//Input files containing contact rates, death rates, and vaccine programs.
    fcontact = fopen("post_UK_contact.dat","r");
    fdeathRate = fopen("death_rate.dat","r");
    fvaccineProg = fopen("vaccine_program_6.dat","r");//Vaccine Program file
	
	//Output files 
	sprintf(fullPath, "%s_SIR_OPTFLAG_1.dat",filenames);
    fout = fopen(fullPath,"w");
    sprintf(fullPath, "%s_conditions_OPTFLAG_1.dat",filenames);
    fconditions = fopen(fullPath,"w");
    sprintf(fullPath, "%s_SEROGROUP_AVG_OPTFLAG_1.dat", filenames);
    fAVGsero = fopen(fullPath,"w");
    sprintf(fullPath, "%s_INFECTION_AGE_OPTFLAG_1.dat",filenames);
    fage = fopen(fullPath,"w");
    sprintf(fullPath, "%s_SEROGROUP_OPTFLAG_1.dat",filenames);
    fsero = fopen(fullPath,"w");
    sprintf(fullPath, "%s_BETA_OPTFLAG_1.dat",filenames);
    flag1_beta = fopen(fullPath, "w");
    sprintf(fullPath, "%s_EPSILON_OPTFLAG_1.dat",filenames);
    flag1_epsilon = fopen(fullPath, "w");
    sprintf(fullPath, "%s_AVG_IMD_1.dat",filenames);
    fAVGIMD = fopen(fullPath,"w");
    sprintf(fullPath, "%s_Condense_demography_1.dat",filenames);
    fConDemo = fopen(fullPath,"w");
    sprintf(fullPath, "%s_Demography_1.dat",filenames);
    fdemograph = fopen(fullPath,"w");
    flag1_stepsize = fopen("stepsize_3_strain_test.dat", "r");
	
	//Create individuals
    population = (pop *)malloc(sizeof(pop)*ARRAY_SIZE);
	
	//Stepsizes to increase/decrease BETA values during calibration
    if(flagOpt == 1){
        for(i=0;i<5;i++){
            for(j=0;j<4;j++){
                //stepsizePrev[j][i] = target[81];
                fscanf(flag1_stepsize,"%lf", &stepsizePrev[j][i]);
                printf("%.8lf, ", stepsizePrev[j][i]);
            }
        }
    }
	//Initialize counters
    for (i=0;i<101;i++){
        numAgentPop[i] = 0;
        infectByAge[i] = 0;
        randInfectCtr[i] = i;
    }
	
	//Obtain contact matric from file
    for(i=0; i<101;i++){
        for(j=0;j<101;j++){
            fscanf(fcontact, "%lf\t", &CONTACT_MATRIX[i][j]);
        }
    }

	//Obtain death rates from file
    for(i=0;i<100;i++){//obtain all-mortality death rates
        fscanf(fdeathRate, "%lf", &DEATH_RATE[i]);
    }
	
	//create all susceptible individuals
    for(i=0;i < (POP_SIZE); i++){
        population[i].status = 0;
        population[i].vaccine = 0;
        population[i].age = (float)(i%90);
        population[i].serogroup = serogroup;
        population[i].vaccineType = ' ';
        population[i].immuneSero = ' ';
        agetemp = (int)(floorf(population[i].age));
        numAgentPop[agetemp]++;
        //printf("%d\n", i);
    }
	
	//Create "dead" individuals
    for(i=POP_SIZE;i<ARRAY_SIZE; i++){
        population[i].status = -1;
        population[i].vaccine = 0;
        population[i].age = -1.0;
        population[i].serogroup = serogroup;
        population[i].vaccineType = serogroup;
        population[i].immuneSero = ' ';
    }
	
    t = 0;
    ////////////////////////////////////////
	//150 year burn-in to obtain proper demographics
    while(t < 150*12){
		
        for(i=0;i<ARRAY_SIZE;i++){//increase age of all indviduals
            if(population[i].status != -1){
                population[i].age = population[i].age + (1.0/12);
                if(fmod(population[i].age, 1.0) > 0.985){
                    population[i].age = ceil(population[i].age);
                }
                else if(fmod(population[i].age, 1.0) < 0.0190){
                    population[i].age = floor(population[i].age);
                }
            }
        }
		
        for(i=0;i<ARRAY_SIZE;i++){//death by old age
            if(population[i].age >= 100.0 && population[i].status > -1){
                population[i].status = -1;
                //POP_SIZE--;
                population[i].serogroup = ' ';
                population[i].immuneSero = ' ';
                population[i].infectTime = 0.0;
                population[i].recoverTime = 0.0;
                population[i].immuneSero = ' ';
                population[i].age = -1.0;
            }
        }
		
        for(i=0;i<ARRAY_SIZE;i++){
            if(population[i].status != -1){
                iseed = (unsigned int)rand()*time(NULL);
                srand(iseed);
                dtemp = (double)(1.0*rand()/RAND_MAX);
                agetemp = (int)(floorl(population[i].age));
                if(agetemp > 90){
                    agetemp = 91;
                }
                if(dtemp < DEATH_RATE[agetemp]) {
                    //POP_SIZE--;
                    population[i].status = -1;
                    population[i].serogroup = ' ';
                    population[i].immuneSero = ' ';
                    population[i].infectTime = 0.0;
                    population[i].recoverTime = 0.0;
                    population[i].immuneSero = ' ';
                    population[i].age = -1.0;
                }
            }
        }
		
		
        //Births of new omni-susceptible individuals
        itemp = BIRTHS;
        for(i=0;i<ARRAY_SIZE;i++){
            if(population[i].status == -1 && itemp > 0){
                population[i].status = 0;
                population[i].age = 0.0;
                population[i].vaccine = 0;
                population[i].vaccineEfficacy = 0.0;
                population[i].vaccineType = 0;
                population[i].infectTime = 0.0;
                population[i].recoverTime = 0.0;
                population[i].serogroup = ' ';
                population[i].immuneSero = ' ';
                population[i].immuneEfficacy = 0.0;
                //POP_SIZE++;
                itemp--;
				
            }
        }
		
        for(i=0;i<100;i++){
            numAgentPop[i] = 0;
        }
        POP_SIZE = 0;
        for(i=0;i<ARRAY_SIZE;i++){
            agetemp = (int)(floor(population[i].age));
            if(population[i].status > -1 && (population[i].age >=0.0 && population[i].age < 100.0)){
                numAgentPop[agetemp]++;
                POP_SIZE++;
            }
        }
        t += dt;
    }
	
	
	//Initial infection
    i = 0;
    while(initInfect > 0){
        infection(&population[i], 0.00);
        population[i].serogroup = probSerogroup(&population[i]);
        population[i].vaccine = 0;
        population[i].vaccineType = serogroup;
        population[i].immuneSero = ' ';
        agetemp = (int)(floorf(population[i].age));
        infectByAge[agetemp]++;
        numAgentPop[agetemp]++;
        i++;
        initInfect --;
    }
	
    for(i=0;i<POP_SIZE;i++){
        agetemp = (int)(floor(population[i].age));
        //printf("agetemp is %d\n", agetemp);
    }
	
    for(i=0;i<ARRAY_SIZE;i++){
        agetemp = (int)(floor(population[i].age));
        if (population[i].status == 1 || population[i].status == 2){
            infectByAge[agetemp]++;
            if(population[i].serogroup == 'b'){
                infectSeroB[agetemp]++;
            }
            else if(population[i].serogroup == 'a'){
                infectSeroA[agetemp]++;
            }
            else if(population[i].serogroup == 'c'){
                infectSeroC[agetemp]++;
            }
            else if(population[i].serogroup == '2'){
                infectSeroNl[agetemp]++;
            }
			
			
        }
        if(population[i].status == 0){
            if(population[i].vaccineType != 'a' && population[i].vaccineType != 'q'){
                numAgentSusc[agetemp][0]++;
            }
            if(population[i].vaccineType != 'b'){
                numAgentSusc[agetemp][1]++;
            }
            if(population[i].vaccineType != 'c' && population[i].vaccineType != 'q'){
                numAgentSusc[agetemp][2] += 1;
            }
            if(population[i].vaccineType != '2'){
                numAgentSusc[agetemp][3] += 1;
            }
        }
        if(population[i].status >= 0 && population[i].status <= 4){
            numAgentPop[agetemp]++;
        }
		
        if (population[i].immuneSero == 'a' && population[i].vaccineType != 'a' && population[i].vaccineType != 'q') {
            numAgentImmune[agetemp][0]++;
        }
        else if (population[i].immuneSero == 'b' && population[i].vaccineType != 'b') {
            numAgentImmune[agetemp][1]++;
        }
        else if (population[i].immuneSero == 'c' && population[i].vaccineType != 'c' && population[i].vaccineType != 'q') {
            numAgentImmune[agetemp][2]++;
        }
        else if (population[i].immuneSero == '2' && population[i].vaccineType != '2') {
            numAgentImmune[agetemp][3]++;
        }
		
    }
	
    for(i=0;i<101;i++){
        numSeroA += infectSeroA[i];
        numSeroB += infectSeroB[i];
        numSeroC += infectSeroC[i];
        numNl += infectSeroNl[i];
    }
    if(flagOpt == 3){
        t = target[1];
    }
    else{
        t = 0.0;
    }

    numInfect = initInfect;
    i = 0;

	
	
    ////////////////////////////////////////////////////////////////////////////
	//Start infection simulation
    while (t <= SIM_LENGTH){
		//If vaccine is being applied, determine who receives the vaccine and apply vaccine to individuals
        if(flagOpt == 4){
            //////////////////////////////////////////////////////////
            //			Vaccination
            for(j=0;j<NUM_PROGRAMS;j++){
                if(monoval[j].sTime <= t  && t <= monoval[j].fTime){
					//When vaccine is available, find people to vaccinate
                    itemp = (int)(floorf(monoval[j].ageApplied));
                    vaccineStart = 1;
                    numVac = (int)(round((double)(numAgentPop[itemp]) * monoval[j].coverage * monoval[j].efficacy));
                    i =0;
                    while(i < ARRAY_SIZE && numVac > 0)
                    {
                        agetemp = (int)(floorf(population[i].age));
                        if((population[i].status == 0 || (population[i].status == 3 && population[i].immuneSero != monoval[j].serogroupProtection)) && population[i].vaccine == 0){
                            if(agetemp == itemp)
                            {
                                applyVaccine(&population[i], &monoval[j], t);
                                k = 0;
                                while(vaccineIndex[k] != -1 && k < 4200){
                                    k++;
                                }
                                vaccineIndex[k] = i;
                                numVac--;
                            }
                        }
                        if (numVac == 0) {
                            firstVac = i;
                        }
                        i++;
                    }
                }
            }
        }
        ///////////////////////////////////////////////////////////////
        //	Calculate the probability of becoming Infectious from Susceptible
		
        for(j=0;j<101;j++){
            pS2IA[j] = 1;
            pS2IB[j] = 1;
            pS2IC[j] = 1;
            pS2INl[j] = 1;
        }
		
        for(i=0;i<101;i++){
            for(j=0;j<101;j++){
                if(numAgentPop[j] <= 0){
                    numAgentPop[j] = 1;
                }
                if(i < 5){
                    TRANS_PROB_A[i][j] = BETA_A_1*CONTACT_MATRIX[i][j]*(double)(infectSeroA[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_B[i][j] = BETA_B_1*CONTACT_MATRIX[i][j]*(double)(infectSeroB[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_C[i][j] = BETA_C_1*CONTACT_MATRIX[i][j]*(double)(infectSeroC[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_NL[i][j] = BETA_NL_1*CONTACT_MATRIX[i][j]*(double)(infectSeroNl[j])/(double)(numAgentPop[j]);
                }
                else if(i >=5 && i < 10){
                    TRANS_PROB_A[i][j] = BETA_A_2*CONTACT_MATRIX[i][j]*(double)(infectSeroA[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_B[i][j] = BETA_B_2*CONTACT_MATRIX[i][j]*(double)(infectSeroB[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_C[i][j] = BETA_C_2*CONTACT_MATRIX[i][j]*(double)(infectSeroC[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_NL[i][j] = BETA_NL_2*CONTACT_MATRIX[i][j]*(double)(infectSeroNl[j])/(double)(numAgentPop[j]);
                }
                else if(i >=10 && i < 15){
                    TRANS_PROB_A[i][j] = BETA_A_3*CONTACT_MATRIX[i][j]*(double)(infectSeroA[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_B[i][j] = BETA_B_3*CONTACT_MATRIX[i][j]*(double)(infectSeroB[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_C[i][j] = BETA_C_3*CONTACT_MATRIX[i][j]*(double)(infectSeroC[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_NL[i][j] = BETA_NL_3*CONTACT_MATRIX[i][j]*(double)(infectSeroNl[j])/(double)(numAgentPop[j]);
                }
                else if(i >=15 && i < 20){
                    TRANS_PROB_A[i][j] = BETA_A_4*CONTACT_MATRIX[i][j]*(double)(infectSeroA[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_B[i][j] = BETA_B_4*CONTACT_MATRIX[i][j]*(double)(infectSeroB[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_C[i][j] = BETA_C_4*CONTACT_MATRIX[i][j]*(double)(infectSeroC[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_NL[i][j] = BETA_NL_4*CONTACT_MATRIX[i][j]*(double)(infectSeroNl[j])/(double)(numAgentPop[j]);
                }
                else if(i >= 20){
                    TRANS_PROB_A[i][j] = BETA_A_5*CONTACT_MATRIX[i][j]*(double)(infectSeroA[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_B[i][j] = BETA_B_5*CONTACT_MATRIX[i][j]*(double)(infectSeroB[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_C[i][j] = BETA_C_5*CONTACT_MATRIX[i][j]*(double)(infectSeroC[j])/(double)(numAgentPop[j]);
                    TRANS_PROB_NL[i][j]= BETA_NL_5*CONTACT_MATRIX[i][j]*(double)(infectSeroNl[j])/(double)(numAgentPop[j]);
                }
				
            }
        }
		
        for(i=0;i<101;i++)
        {
            for(j=0;j<101;j++){
                pS2IA[i] *= 1.0 - TRANS_PROB_A[i][j];
                pS2IB[i] *= 1.0 - TRANS_PROB_B[i][j];
                pS2IC[i] *= 1.0 - TRANS_PROB_C[i][j];
                pS2INl[i] *= 1.0 - TRANS_PROB_NL[i][j];
				
            }
        }
		
        for(i=0;i<101;i++){
            pS2IA[i] = 1.0 - pS2IA[i];
            pS2IB[i] = 1.0 - pS2IB[i];
            pS2IC[i] = 1.0 - pS2IC[i];
            pS2INl[i] = 1.0 - pS2INl[i];
        }
		
        ////////////////////////////////////////////////////////////////////////////////////////////
        //		End of Probabilities and start of transitions
        orderSero(infectBySero[0], infectBySero[1],infectBySero[2],infectBySero[3], numInfect+numIMD);
        for(i=0;i<4;i++){
            infectBySero[i] = 0;
        }
		
		for(i=0;i<ARRAY_SIZE;i++){
            randInfectCtr[i] = i;
        }
        
		for(k=0;k<4;k++){
            shuffle(randInfectCtr);
            switch (SEROORDER[k]){
                case 0:
                    for(i=0;i<ARRAY_SIZE;i++){
                        itemp = randInfectCtr[i];
                        agetemp = (int)(floorf(population[itemp].age));
                        if(population[itemp].status == 0 && (population[itemp].vaccineType != 'a' && population[itemp].vaccineType != 'q'))
                        {
                            randNumInfect = coinFlip();
                            if(ranCheck(pS2IA[agetemp], randNumInfect)==1){
                                population[itemp].serogroup ='a';
                                infection(&population[itemp], t);
                                if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_A)){
                                    population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_A);
                                }
                            }
                        }
                    }
					
                    break;
					
                case 1:
                    for(i=0;i<ARRAY_SIZE;i++){
                        itemp = randInfectCtr[i];
                        agetemp = (int)(floorf(population[itemp].age));
                        if(population[itemp].status == 0 && (population[itemp].vaccineType != 'b'))
                        {
                            randNumInfect = coinFlip();
                            if(ranCheck(pS2IB[agetemp], randNumInfect)==1){
                                population[itemp].serogroup ='b';
                                infection(&population[itemp], t);
                                if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_B)){
                                    population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_B);
                                }
                            }
                        }
                    }
                    break;
					
                case 2:
                    for(i=0;i<ARRAY_SIZE;i++){
                        itemp = randInfectCtr[i];
                        agetemp = (int)(floorf(population[itemp].age));
                        if(population[itemp].status == 0 && (population[itemp].vaccineType != 'c' && population[itemp].vaccineType !='q'))
                        {
                            randNumInfect = coinFlip();
                            if(ranCheck(pS2IC[agetemp], randNumInfect)==1){
                                population[itemp].serogroup ='c';
                                infection(&population[itemp], t);
                                if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_C)){
                                    population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_C);
                                }
                            }
                        }
                    }
                    break;
					
                case 3:
                    for(i=0;i<ARRAY_SIZE;i++){
                        itemp = randInfectCtr[i];
                        agetemp = (int)(floorf(population[itemp].age));
                        if(population[itemp].status == 0 && (population[itemp].vaccineType != '2'))
                        {
                            randNumInfect = coinFlip();
                            if(ranCheck(pS2INl[agetemp], randNumInfect)==1){
                                population[itemp].serogroup ='2';
                                infection(&population[itemp], t);
                                if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_NL)){
                                    population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_NL);
                                }
                            }
                        }
                    }
                    break;
					
            }
        }
		
        ///////////////////////////////////////////////////////////////////
        //	Calculate the probability of becoming Infectious from Removed
        for(j=0;j<101;j++){
            for(i=0;i<6;i++){
                pR2IA[i][j] = 1;
                pR2IB[i][j] = 1;
                pR2IC[i][j] = 1;
                pR2INl[i][j] = 1;
            }
        }
		
        for(i=0;i<101;i++){
            for(j=0;j<101;j++){
                pR2IA[0][i] *= 1.0 - TRANS_PROB_A[i][j] *(1.0-CROSS_B);
                pR2IA[1][i] *= 1.0 - TRANS_PROB_A[i][j] *(1.0-CROSS_C);
                pR2IA[2][i] *= 1.0 - TRANS_PROB_A[i][j] *(1.0-CROSS_NL);
				
                pR2IB[0][i] *= 1.0 - TRANS_PROB_B[i][j] *(1.0-CROSS_A);
                pR2IB[1][i] *= 1.0 - TRANS_PROB_B[i][j] *(1.0-CROSS_C);
                pR2IB[2][i] *= 1.0 - TRANS_PROB_B[i][j] *(1.0-CROSS_NL);
				
                pR2IC[0][i] *= 1.0 - TRANS_PROB_C[i][j] *(1.0-CROSS_A);
                pR2IC[1][i] *= 1.0 - TRANS_PROB_C[i][j] *(1.0-CROSS_B);
                pR2IC[2][i] *= 1.0 - TRANS_PROB_C[i][j] *(1.0-CROSS_NL);
				
                pR2INl[0][i] *= 1.0 - TRANS_PROB_NL[i][j] *(1.0-CROSS_A);
                pR2INl[1][i] *= 1.0 - TRANS_PROB_NL[i][j] *(1.0-CROSS_B);
                pR2INl[2][i] *= 1.0 - TRANS_PROB_NL[i][j] *(1.0-CROSS_C);
				
            }
        }
		
        for(i=0;i<101;i++){
            for(j=0;j<3;j++){
                pR2IA[j][i] = 1.0 - pR2IA[j][i];
                pR2IB[j][i] = 1.0 - pR2IB[j][i];
                pR2IC[j][i] = 1.0 - pR2IC[j][i];
                pR2INl[j][i] = 1.0 - pR2INl[j][i];
            }
        }
		////////////////////////////////////////////////////////////////////////////////////////////
        //		End of Probabilities and start of transitions	
        for(k=0;k<4;k++){
            shuffle(randInfectCtr);
            switch (SEROORDER[k]){
                case 0:
                    for(j=0;j<3;j++){
                        for(i=0;i<ARRAY_SIZE;i++){
                            itemp = randInfectCtr[i];
							agetemp = (int)(floorf(population[itemp].age));
							if(population[itemp].status == 3 && (population[itemp].vaccineType != 'a' && population[itemp].vaccineType != 'q' && population[itemp].immuneSero != 'a'))
							{
								randNumInfect = coinFlip();
								if(j==0 && population[itemp].immuneSero != 'b' && ranCheck(pR2IA[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='a';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_A)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_A);
									}
								}
								if(j==1 && population[itemp].immuneSero != 'c' && ranCheck(pR2IA[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='a';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_A)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_A);
									}
								}
								if(j==2 && population[itemp].immuneSero != '2' && ranCheck(pR2IA[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='a';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_A)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_A);
									}
								}
							}
						}
					}
                    break;
					
                case 1:
                    for(j=0;j<3;j++){
                        for(i=0;i<ARRAY_SIZE;i++){
							
							itemp = randInfectCtr[i];
							agetemp = (int)(floorf(population[itemp].age));
							if(population[itemp].status == 3 && (population[itemp].vaccineType != 'b' && population[itemp].immuneSero != 'b'))
							{
								randNumInfect = coinFlip();
								if(j==0 && population[itemp].immuneSero != 'a' && ranCheck(pR2IB[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='b';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_B)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_B);
									}
								}
								if(j==1 && population[itemp].immuneSero != 'c' && ranCheck(pR2IB[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='b';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_B)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_B);
									}
								}
								if(j==2 && population[itemp].immuneSero != '2' && ranCheck(pR2IB[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='b';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_B)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_B);
									}
								}
							}
						}
					}
					break;
					
                case 2:
                    for(j=0;j<3;j++){
                        for(i=0;i<ARRAY_SIZE;i++){
							
							itemp = randInfectCtr[i];
							agetemp = (int)(floorf(population[itemp].age));
							if(population[itemp].status == 3 && (population[itemp].vaccineType != 'c' && population[itemp].vaccineType != 'q' && population[itemp].immuneSero != 'c'))
							{
								randNumInfect = coinFlip();
								if(j==0 && population[itemp].immuneSero != 'a' && ranCheck(pR2IC[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='c';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_C)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_C);
									}
								}
								if(j==1 && population[itemp].immuneSero != 'b' && ranCheck(pR2IC[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='c';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_C)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_C);
									}
								}
								if(j==2 && population[itemp].immuneSero != '2' && ranCheck(pR2IC[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='c';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_C)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_C);
									}
								}
							}
						}
					}
					break;
					
                case 3:
                    for(j=0;j<3;j++){
                        for(i=0;i<ARRAY_SIZE;i++){
							itemp = randInfectCtr[i];
							agetemp = (int)(floorf(population[itemp].age));
							if(population[itemp].status == 3 && (population[itemp].vaccineType != '2' && population[itemp].immuneSero != '2'))
							{
								randNumInfect = coinFlip();
								if(j==0 && population[itemp].immuneSero != 'a' && ranCheck(pR2INl[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='2';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_NL)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_NL);
									}
								}
								if(j==1 && population[itemp].immuneSero != 'b' && ranCheck(pR2INl[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='2';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_NL)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_NL);
									}
								}
								if(j==2 && population[itemp].immuneSero != 'c' && ranCheck(pR2INl[j][agetemp], randNumInfect)==1){
									population[itemp].serogroup ='2';
									infection(&population[itemp], t);
									if(population[itemp].infectTime > t+(1.1*AVG_DUR_INFECT_NL)){
										population[itemp].infectTime = t+(1.1*AVG_DUR_INFECT_NL);
									}
								}
							}
						}
					}
					break;
            }
        }
		
        ///////////////////////////////////////////////////////////////////
        //	Check timers to see if individuals change from Infectious to Removed
        for(i=0;i<ARRAY_SIZE;i++){
			if((population[i].status == 1 || population[i].status == 2) && population[i].infectTime <= t){
                population[i].status = 3;
                population[i].infectTime = 0;
                if(population[i].serogroup == 'a'){
                    population[i].recoverTime = phaseTimer(AVG_DUR_IMMUNE_A, 1.0) + t;
                }
                if(population[i].serogroup == 'b'){
                    population[i].recoverTime = phaseTimer(AVG_DUR_IMMUNE_B, 1.0) + t;
                }
                if(population[i].serogroup == 'c'){
                    population[i].recoverTime = phaseTimer(AVG_DUR_IMMUNE_C, 1.0) + t;
                }
                if(population[i].serogroup == '2'){
                    population[i].recoverTime = phaseTimer(AVG_DUR_IMMUNE_NL, 1.50) + t;
                }
                population[i].infectTime = 0.0;
                population[i].immuneSero = population[i].serogroup;
                population[i].serogroup = ' ';
                population[i].immuneEfficacy = 1.00;
				
            }
			
            /////////////////////////////////////////////////////////
            //		Check timers to see if individuals change from Removed to Susceptible
            if(population[i].status == 3 && population[i].recoverTime <= t)
            {
                population[i].status = 0;
                population[i].vaccine = 0;
                population[i].infectTime = 0.0;
                population[i].recoverTime = 0.0;
                population[i].serogroup = ' ' ;
                population[i].immuneSero = ' ';
            }
        }
		
        //////////////////////////////////////////////////////
        //	Update demographics (age, births, deaths due to old age or from all-mortality death rates)
		
        for(i=0;i<ARRAY_SIZE;i++){//increase age of all indviduals
            if(population[i].status != -1){
                population[i].age = population[i].age + (1.0/12);//increase by one week if they are not dead
                if(fmod(population[i].age, 1.0) > 0.985){
                    population[i].age = ceil(population[i].age);
                }
                else if(fmod(population[i].age, 1.0) < 0.0190){
                    population[i].age = floor(population[i].age);
                }
            }
        }
		
        itemp = BIRTHS;
        for(i=0;i<ARRAY_SIZE;i++){
            if(population[i].status == -1 && itemp > 0){
                population[i].status = 0;
                population[i].age = 0.0;
                population[i].vaccine = 0;
                population[i].vaccineEfficacy = 0.0;
                population[i].vaccineType = 0;
                population[i].infectTime = 0.0;
                population[i].recoverTime = 0.0;
                population[i].serogroup = ' ';
                population[i].immuneSero = ' ';
                population[i].immuneEfficacy = 0.0;
                itemp--;
				
            }
        }
		
        for(i=0;i<ARRAY_SIZE;i++){//death by old age
            if(population[i].age >= 100.0 && population[i].status > -1){
                population[i].status = -1;
                population[i].serogroup = ' ';
                population[i].immuneSero = ' ';
                population[i].infectTime = 0.0;
                population[i].recoverTime = 0.0;
                population[i].age = -1.0;
            }
        }
		
        for(i=0;i<ARRAY_SIZE;i++){
            if(population[i].status != -1){
                iseed = (unsigned int)rand()*time(NULL);
                srand(iseed);
                dtemp = (double)(1.0*rand()/RAND_MAX);
                agetemp = (int)(floorf(population[i].age));
                if(agetemp > 90){
                    agetemp = 91;
                }
                if(dtemp < DEATH_RATE[agetemp]) {
                    population[i].status = -1;
                    population[i].serogroup = ' ';
                    population[i].immuneSero = ' ';
                    population[i].infectTime = 0.0;
                    population[i].recoverTime = 0.0;
                    population[i].immuneSero = ' ';
                    population[i].age = -1.0;
                }
            }
        }
		
		
        /////////////////////////////
        // Counter updates
		
		
        numSusc = 0; numInfect =0; numIMD = 0; numRecover = 0;
        //// Status tracking update
        for(i=0;i<ARRAY_SIZE;i++){
            agetemp = (int)(floorf(population[i].age));
            if(population[i].status == 0 && agetemp < 100)
            {numSusc++;}
            else if(population[i].status == 1 && agetemp < 100)
            {numInfect++;}
            else if(population[i].status == 2 && agetemp < 100)
            {numIMD++;}
            else if(population[i].status == 3 && agetemp < 100)
            {numRecover++;}
			
        }

        POP_SIZE = numSusc+numInfect+numIMD+numRecover;

		
        //Infect by Age  and number of agents status update
        for(i=0;i< 101;i++){
            infectByAge[i] = 0;
            numAgentPop[i] = 0;
            numAgentSusc[i][0] = 0;
            numAgentSusc[i][1] = 0;
            numAgentSusc[i][2] = 0;
            numAgentSusc[i][3] = 0;
            numAgentImmune[i][0] = 0;
            numAgentImmune[i][1] = 0;
            numAgentImmune[i][2] = 0;
            numAgentImmune[i][3] = 0;
            infectSeroA[i] = 0;
            infectSeroB[i] = 0;
            infectSeroC[i] = 0;
            infectSeroNl[i] = 0;
        }
        for(i=0;i<5;i++){
            for(j=0;j<4;j++){
                condensePrev[j][i] = 0.0;
                condenseSero[j][i] = 0;
                infectBySero[j] = 0;
            }
            if(flagOpt == 2){
                for(j=0;j<3;j++){
                    condenseIMD[j][i] = imdCtr[j][i];
                    IMDPerYear += (int)(imdCtr[j][i]);
                }
            }
            condenseAge[i] = 0;
        }
		
        for(i=0;i<ARRAY_SIZE;i++){
            if(population[i].status == 1){
                if(population[i].serogroup == 'a'){
                    infectBySero[0] = infectBySero[0] + 1.0;
                }
                else if(population[i].serogroup == 'b'){
                    infectBySero[1] = infectBySero[1] + 1.0;
                }
                else if(population[i].serogroup == 'c'){
                    infectBySero[2] = infectBySero[2] + 1.0;
                }
                else{
                    infectBySero[3] = infectBySero[3] + 1.0;
                }
            }
        }
		
		
        for(i=0;i<ARRAY_SIZE;i++){
            agetemp = (int)(floorf(population[i].age));
            if(agetemp >= 0 && agetemp < 5)
            {condenseAge[0]++;}
            else if(agetemp >= 5 && agetemp < 10)
            {condenseAge[1]++;}
            else if(agetemp >= 10 && agetemp < 15)
            {condenseAge[2]++;}
            else if(agetemp >= 15 && agetemp < 20)
            {condenseAge[3]++;}
            else if(agetemp >= 20)
            {condenseAge[4]++;}
			
			
            if ((population[i].status == 1 || population[i].status == 2) && agetemp >= 0){
                infectByAge[agetemp]++;
                if(population[i].serogroup == 'a'){
                    infectSeroA[agetemp]++;
                    if(agetemp >= 0 && agetemp < 5)
                    {
                        condenseSero[0][0]++;
                        if(population[i].status==1){
                            condensePrev[0][0]++;
                        }
                    }
                    else if(agetemp >= 5 && agetemp < 10)
                    {
                        condenseSero[0][1]++;
                        if(population[i].status==1){
                            condensePrev[0][1]++;
                        }
                    }
                    else if(agetemp >= 10 && agetemp < 15)
                    {
                        condenseSero[0][2]++;
                        if(population[i].status==1){
                            condensePrev[0][2]++;
                        }
                    }
                    else if(agetemp >= 15 && agetemp < 20)
                    {
                        condenseSero[0][3]++;
                        if(population[i].status==1){
                            condensePrev[0][3]++;
                        }
                    }
                    else if(agetemp >= 20)
                    {
                        condenseSero[0][4]++;
                        if(population[i].status==1){
                            condensePrev[0][4]++;
                        }
                    }
                }
                else if(population[i].serogroup == 'b'){
                    infectSeroB[agetemp]++;
                    if(agetemp >= 0 && agetemp < 5)
                    {
                        condenseSero[1][0]++;
                        if(population[i].status==1){
                            condensePrev[1][0]++;
                        }
                    }
                    else if(agetemp >= 5 && agetemp < 10)
                    {
                        condenseSero[1][1]++;
                        if(population[i].status==1){
                            condensePrev[1][1]++;
                        }
                    }
                    else if(agetemp >= 10 && agetemp < 15)
                    {
                        condenseSero[1][2]++;
                        if(population[i].status==1){
                            condensePrev[1][2]++;
                        }
                    }
                    else if(agetemp >= 15 && agetemp < 20)
                    {
                        condenseSero[1][3]++;
                        if(population[i].status==1){
                            condensePrev[1][3]++;
                        }
                    }
                    else if(agetemp >= 20)
                    {
                        condenseSero[1][4]++;
                        if(population[i].status==1){
                            condensePrev[1][4]++;
                        }
                    }
                }
				
                else if(population[i].serogroup == 'c'){
                    infectSeroC[agetemp]++;
                    if(agetemp >= 0 && agetemp < 5)
                    {
                        condenseSero[2][0]++;
                        if(population[i].status==1){
                            condensePrev[2][0]++;
                        }
                    }
                    else if(agetemp >= 5 && agetemp < 10)
                    {
                        condenseSero[2][1]++;
                        if(population[i].status==1){
                            condensePrev[2][1]++;
                        }
                    }
                    else if(agetemp >= 10 && agetemp < 15)
                    {
                        condenseSero[2][2]++;
                        if(population[i].status==1){
                            condensePrev[2][2]++;
                        }
                    }
                    else if(agetemp >= 15 && agetemp < 20)
                    {
                        condenseSero[2][3]++;
                        if(population[i].status==1){
                            condensePrev[2][3]++;
                        }
                    }
                    else if(agetemp >= 20)
                    {
                        condenseSero[2][4]++;
                        if(population[i].status==1){
                            condensePrev[2][4]++;
                        }
                    }
                }
				
                else if(population[i].serogroup == '2'){
                    infectSeroNl[agetemp]++;
                    if(agetemp >= 0 && agetemp < 5)
                    {
                        condenseSero[3][0]++;
                        if(population[i].status==1){
                            condensePrev[3][0]++;
                        }
                    }
                    else if(agetemp >= 5 && agetemp < 10)
                    {
                        condenseSero[3][1]++;
                        if(population[i].status==1){
                            condensePrev[3][1]++;
                        }
                    }
                    else if(agetemp >= 10 && agetemp < 15)
                    {
                        condenseSero[3][2]++;
                        if(population[i].status==1){
                            condensePrev[3][2]++;
                        }
						
                    }
                    else if(agetemp >= 15 && agetemp < 20)
                    {
                        condenseSero[3][3]++;
                        if(population[i].status==1){
                            condensePrev[3][3]++;
                        }
                    }
                    else if(agetemp >= 20)
                    {
                        condenseSero[3][4]++;
                        if(population[i].status==1){
                            condensePrev[3][4]++;
                        }
                    }
                }
				
            }
            if(population[i].status == 0){
                if(population[i].vaccineType != 'a' && population[i].vaccineType != 'q'){
                    numAgentSusc[agetemp][0]++;
                }
                if(population[i].vaccineType != 'b'){
                    numAgentSusc[agetemp][1]++;
                }
                if(population[i].vaccineType != 'c' && population[i].vaccineType != 'q'){
                    numAgentSusc[agetemp][2] += 1;
                }
                if(population[i].vaccineType != '2'){
                    numAgentSusc[agetemp][3] += 1;
                }
            }
            if(population[i].status >= 0 && population[i].status <= 4){
                numAgentPop[agetemp]++;
            }
			
            if (population[i].immuneSero == 'a' && population[i].vaccineType != 'a' && population[i].vaccineType != 'q') {
                numAgentImmune[agetemp][0]++;
            }
            else if (population[i].immuneSero == 'b' && population[i].vaccineType != 'b') {
                numAgentImmune[agetemp][1]++;
            }
            else if (population[i].immuneSero == 'c' && population[i].vaccineType != 'c' && population[i].vaccineType != 'q') {
                numAgentImmune[agetemp][2]++;
            }
            else if (population[i].immuneSero == '2' && population[i].vaccineType != '2') {
                numAgentImmune[agetemp][3]++;
            }
        }
        
        totalPasses++;
		///////////////////////////////////////////
		//Find the running average of prevalence of carriage
        if((int)(t)%1200==0){
            for(j=0;j<5;j++){
                for(i=0;i<4;i++){
                    for(k=0;k<=600;k++){
                        avgDataPrev[i][j][k] = 0;
                    }
				}
                for(i=0;i<3;i++){
                    for(k=0;k<=600;k++){
                        avgDataIMD[i][j][k] = 0;
                    }
                }
            }
        }
		
        if(t >= 0.0 && (flagOpt == 4)){
            if((int)(t) >= 120){
                itemp = (int)(t) - 120;
                for(j=0;j<5;j++){
                    for(i=0;i<4;i++){
                        avgDataPrev3[i][j][itemp] = (double)(condensePrev[i][j])/(double)(condenseAge[j]);
                        prevAVG[i][j] = runningAVG(itemp+1, &avgDataPrev3[i][j]);
                    }
                    for(i=0;i<3;i++){
                        avgDataIMD3[i][j][itemp] = IMD100KPop[i][j];
                        imdAVG[i][j] = runningAVG(itemp+1, &avgDataIMD3[i][j]);
                    }
                }
            }
        }
        //produce a running average
        else{
            if(((int)(t)%1200 >= 300) && t > 0.0){
                itemp = ((int)(t)%1200) - 300;
                for(j=0;j<5;j++){
                    for(i=0;i<4;i++){
                        avgDataPrev[i][j][itemp] = (double)(condensePrev[i][j])/(double)(condenseAge[j]);
                        prevAVG[i][j] = runningAVG(itemp+1, &avgDataPrev[i][j]);
                    }
                    for(i=0;i<3;i++){
                        avgDataIMD[i][j][itemp] = IMD100KPop[i][j];
                        imdAVG[i][j] = runningAVG(itemp+1, &avgDataIMD[i][j]);
                    }
                }
            }
        }
		
		
		
        itemp = 0;
		
        ///////////////////////////////////////////////////////////////////////////////////////////
        //									Opt Flag 1											 //
        ///////////////////////////////////////////////////////////////////////////////////////////
		// Compare running average to targets. If relative error is too large, adjustments are needed to BETA values
        if(flagOpt <= 1){
            //itemp = (int)(t/52.0) % 150;
            if ((int)(t)%1200 == 1199 && t != 0.0){
                if(fabs((target[45] - prevAVG[0][0])/target[45]) < target[67])//A_1
				{itemp++;}
                if(fabs((target[46] - prevAVG[0][1])/target[46]) < target[68])//A_2
				{itemp++;}
                if(fabs((target[47] - prevAVG[0][2])/target[47]) < target[69])//A_3
				{itemp++;}
                if(fabs((target[48] - prevAVG[0][3])/target[48]) < target[70])//A_4
				{itemp++;}
                if(fabs((target[49] - prevAVG[0][4])/target[49]) < target[71])//A_5
				{itemp++;}
                if(fabs((target[50] - prevAVG[1][0])/target[50]) < target[72])//B_1
				{itemp++;}
                if(fabs((target[51] - prevAVG[1][1])/target[51]) < target[73])//B_2
				{itemp++;}
                if(fabs((target[52] - prevAVG[1][2])/target[52]) < target[74])//B_3
				{itemp++;}
                if(fabs((target[53] - prevAVG[1][3])/target[53]) < target[75])//B_4
				{itemp++;}
                if(fabs((target[54] - prevAVG[1][4])/target[54]) < target[76])//B_5
				{itemp++;}
                if(fabs((target[55] - prevAVG[2][0])/target[55]) < target[77])//C_1
				{itemp++;}
                if(fabs((target[56] - prevAVG[2][1])/target[56]) < target[78])//C_2
				{itemp++;}
                if(fabs((target[57] - prevAVG[2][2])/target[57]) < target[79])//C_3
				{itemp++;}
                if(fabs((target[58] - prevAVG[2][3])/target[58]) < target[80])//C_4
				{itemp++;}
                if(fabs((target[59] - prevAVG[2][4])/target[59]) < target[81])//C_5
				{itemp++;}
                if(fabs((target[60] - prevAVG[3][0])/target[60]) < target[82])//NL_1
				{itemp++;}
                if(fabs((target[61] - prevAVG[3][1])/target[61]) < target[83])//NL_2
				{itemp++;}
                if(fabs((target[62] - prevAVG[3][2])/target[62]) < target[84])//NL_3
				{itemp++;}
                if(fabs((target[63] - prevAVG[3][3])/target[63]) < target[85])//NL_4
				{itemp++;}
                if(fabs((target[64] - prevAVG[3][4])/target[64]) < target[86])//NL_5
				{itemp++;}
				
				///////////////////////////////////////////////////////////////
				// If calibration is successful start vaccine scenarios
                if(itemp >=20 && condensePrev[2][0] > 0.0 && condensePrev[2][1] > 0.0 && condensePrev[2][2] > 0.0 && condensePrev[2][3] > 0.0 && condensePrev [2][4] > 0.0){
					flagOpt = 4;
					//print initial conditions of each individual
                    for(i=0;i<ARRAY_SIZE;i++){
                        fprintf(fconditions,"%d, %f, %d, %c, %f, %c, %f, %f, %c\n", population[i].status, population[i].age, population[i].vaccine, population[i].vaccineType, population[i].vaccineEfficacy, population[i].serogroup, population[i].infectTime, population[i].recoverTime, population[i].immuneSero);
                    }
                    
					t= -600.0;
					SIM_LENGTH = 300*12;
					
					//Close old output files and open new output files
                    fclose(fout);
                    fclose(fconditions);
                    fclose(fAVGsero);
                    fclose(fAVGIMD);
                    fclose(fage);
                    fclose(fsero);
                    fclose(flag1_epsilon);
                    fclose(flag1_beta);
                    fclose(flag1_stepsize);
                    fclose(fConDemo);
                    fclose(fdemograph);
					
                    sprintf(fullPath, "%s_SIR_OPTFLAG_4.dat",filenames);
                    fout = fopen(fullPath,"w");
		
                    sprintf(fullPath, "%s_SEROGROUP_AVG_OPTFLAG_4.dat", filenames);
                    fAVGsero = fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_INFECTION_AGE_OPTFLAG_4.dat",filenames);
                    fage = fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_SEROGROUP_OPTFLAG_4.dat",filenames);
                    fsero = fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_Condense_demography_4.dat",filenames);
                    fConDemo = fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_Demography_4.dat",filenames);
                    fdemograph = fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_age_strat_A.dat",filenames);
                    fprevA = fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_age_strat_B.dat",filenames);
                    fprevB= fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_age_strat_C.dat",filenames);
                    fprevC= fopen(fullPath,"w");
					
                    sprintf(fullPath, "%s_age_strat_NL.dat",filenames);
                    fprevNl= fopen(fullPath,"w");
					
					// Load vaccine programs
					monoval = (vaccineProgram *)malloc(sizeof(vaccineProgram));
					i = 0;
					while(fscanf(fvaccineProg, "%d", &itemp) != EOF){
						monovalTmp = (vaccineProgram *)realloc(monoval,(i+1)*sizeof(vaccineProgram));
						monoval = monovalTmp;
						
						if(itemp == 1){
							monoval[i].serogroupProtection = 'c';	
						}
						else if(itemp == 2){
							monoval[i].serogroupProtection = 'q';	
						}
						else if(itemp == 3){
							monoval[i].serogroupProtection = 'b';	
						}
						fscanf(fvaccineProg, "%f", &monoval[i].coverage);
						fscanf(fvaccineProg, "%f", &monoval[i].efficacy);
						fscanf(fvaccineProg, "%f", &monoval[i].sTime);
						fscanf(fvaccineProg, "%f", &monoval[i].fTime);
						fscanf(fvaccineProg, "%d", &monoval[i].ageApplied);
						fscanf(fvaccineProg, "%f", &monoval[i].lengthVaccinated);
						i++; NUM_PROGRAMS++;
					}
					for (i=0; i< NUM_PROGRAMS;i++) {
						printf("Vaccine: %c %f %f %.1f %d\n", monoval[i].serogroupProtection, monoval[i].sTime, monoval[i].fTime, monoval[i].ageApplied, i);
					}
                }
				
				///////////////////////////////////////////////////////////////
				// If calibration is unsuccessful change parameters
				if(flagOpt == 1){
					// Sero A Prev
					i = 0; j = 0;
					if(fabs((target[45] - prevAVG[0][0])/target[45]) > target[67] && target[45] != 0.0){
						if (((target[46] - prevAVG[0][1])/target[46]) > 0.0 && signPrev[1] == 1){
							if(fabs((target[45] - prevAVG[0][0])/target[45]) > (2*target[67])){
								BETA_A_1 += 2*stepsizePrev[0][0];
							}
							else{
								BETA_A_1 += stepsizePrev[0][0];
								
							}
							signPrev[0] = 1;
							if (prevAVG[0][0] < 0.001 && target[45] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						
						else if(((target[45] - prevAVG[0][0])/target[45]) > 0.0 && signPrev[0] == -1 ){
							stepsizePrev[0][0] = stepsizePrev[0][0]/2.0;
							signPrev[0] = 1;
							BETA_A_1 += stepsizePrev[0][0];
							printf("A_1 test 2 has passed BETA A_1: %lf\n", BETA_A_1);
							if (prevAVG[0][0] < 0.001 && target[45] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[45] - prevAVG[0][0])/target[45]) < 0.0 && signPrev[0] == -1){
							if(fabs((target[45] - prevAVG[0][0])/target[45]) > (2*target[67])){
								BETA_A_1 -= 2*stepsizePrev[0][0];
							}
							else{
								BETA_A_1 -= stepsizePrev[0][0];
								
							}
							
							signPrev[0] = -1;
						}
						
						else if(((target[45] - prevAVG[0][0])/target[45]) < 0.0 && signPrev[0] == 1){
							stepsizePrev[0][0] = stepsizePrev[0][0]/2.0;
							BETA_A_1 -= stepsizePrev[0][0];
							printf("A_1 test 4 has passed BETA A_1: %lf\n", BETA_A_1);
							signPrev[0] = -1;
						}
					}
					else if(target[45] != 0.0 && BETA_A_1 == 0.0){
						BETA_A_1 += stepsizePrev[0][0];
					}
					else if(target[45] != 0.0 && prevAVG[0][0] == 0.0){
						BETA_A_1 += stepsizePrev[0][0];
					}
					else if(target[45] == 0.0 && prevAVG[0][0] != 0.0){
						BETA_A_1 = 0.0;
						if (BETA_A_1 < 0.0){
							BETA_A_1 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[46] - prevAVG[0][1])/target[46]) > target[68] && target[46] != 0.0){
						if (((target[46] - prevAVG[0][1])/target[46]) > 0.0 && signPrev[1] == 1){
							if(fabs((target[46] - prevAVG[0][1])/target[46]) > (2*target[68])){
								BETA_A_2 += 2*stepsizePrev[0][1];
							}
							else{
								BETA_A_2 += stepsizePrev[0][1];
							}
							signPrev[1] = 1;
							if (prevAVG[0][1] < 0.001 && target[46] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[46] - prevAVG[0][1])/target[46]) > 0.0 && signPrev[1] == -1 ){
							stepsizePrev[0][1] = stepsizePrev[0][1]/2.0;//Halve the stepsize
							BETA_A_2 += stepsizePrev[0][1];
							printf("A_2 test 2 has passed BETA A_2: %lf\n", BETA_A_2);
							signPrev[1] = 1;
							if (prevAVG[0][1] < 0.001 && target[46] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[46] - prevAVG[0][1])/target[46]) < 0.0 && signPrev[1] == -1){
							if(fabs((target[46] - prevAVG[0][1])/target[46]) > (2*target[68])){
								BETA_A_2 -= 2*stepsizePrev[0][1];
							}
							else{
								BETA_A_2 -= stepsizePrev[0][1];
							}
							signPrev[1] = -1;
						}
						
						else if(((target[46] - prevAVG[0][1])/target[46]) < 0.0 && signPrev[1] == 1){
							stepsizePrev[0][1] = stepsizePrev[0][1]/2.0;//halve the stepsize
							BETA_A_2 -= stepsizePrev[0][1];
							printf("A_2 test 4 has passed BETA A_2: %lf\n", BETA_A_2);
							signPrev[1] = -1;
						}
					}
					else if(target[46] != 0.0 && BETA_A_2 == 0.0){
						BETA_A_2 += stepsizePrev[0][1];
					}
					else if(target[46] != 0.0 && prevAVG[0][1] == 0.0){
						BETA_A_2 += stepsizePrev[0][1];
					}
					else if(target[46] == 0.0 && prevAVG[0][1] != 0.0){
						BETA_A_2 = 0.0;
						if (BETA_A_2 < 0.0){
							BETA_A_2 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[47] - prevAVG[0][2])/target[47]) > target[69] && target[47] != 0.0){
						if (((target[47] - prevAVG[0][2])/target[47]) > 0.0 && signPrev[2] == 1){
							if(fabs((target[47] - prevAVG[0][2])/target[47]) > (2*target[69])){
								BETA_A_3 += 2*stepsizePrev[0][2];
							}
							else{
								BETA_A_3 += stepsizePrev[0][1];
							}
							signPrev[2] = 1;
							if (prevAVG[0][2] < 0.001 && target[47] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[47] - prevAVG[0][2])/target[47]) > 0.0 && signPrev[2] == -1 ){
							stepsizePrev[0][2] = stepsizePrev[0][2]/2.0;//Halve the stepsize
							BETA_A_3 += stepsizePrev[0][2];
							printf("A_3 test 2 has passed BETA A_3: %lf\n", BETA_A_3);
							signPrev[2] = 1;
							if (prevAVG[0][2] < 0.001 && target[47] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[47] - prevAVG[0][2])/target[47]) < 0.0 && signPrev[2] == -1){
							if(fabs((target[47] - prevAVG[0][2])/target[47]) > (2*target[69])){
								BETA_A_3 -= 2*stepsizePrev[0][2];
							}
							else{
								BETA_A_3 -= stepsizePrev[0][2];
							}
							signPrev[2] = -1;
						}
						
						else if(((target[47] - prevAVG[0][2])/target[47]) < 0.0 && signPrev[2] == 1){
							stepsizePrev[0][2] = stepsizePrev[0][2]/2.0;//halve the stepsize
							BETA_A_3 -= stepsizePrev[0][2];
							printf("A_3 test 4 has passed BETA A_3: %lf\n", BETA_A_3);
							signPrev[2] = -1;
						}
					}
					else if(target[47] != 0.0 && BETA_A_3 == 0.0){
						BETA_A_3 += stepsizePrev[0][2];
					}
					else if(target[47] != 0.0 && prevAVG[0][2] == 0.0){
						BETA_A_3 += stepsizePrev[0][2];
					}
					else if(target[47] == 0.0 && prevAVG[0][2] != 0.0){
						BETA_A_3 = 0.0;
						if (BETA_A_3 < 0.0){
							BETA_A_3 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[48] - prevAVG[0][3])/target[48]) > target[70] && target[48] != 0.0){
						if (((target[48] - prevAVG[0][3])/target[48]) > 0.0 && signPrev[3] == 1){
							if(fabs((target[48] - prevAVG[0][3])/target[48]) > (2*target[70])){
								BETA_A_4 += 2*stepsizePrev[0][3];
							}
							else{
								BETA_A_4 += stepsizePrev[0][3];
							}						signPrev[3] = 1;
							if (prevAVG[0][3] < 0.001 && target[48] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[48] - prevAVG[0][3])/target[48]) > 0.0 && signPrev[3] == -1 ){
							stepsizePrev[0][3] = stepsizePrev[0][3]/2.0;//Halve the stepsize
							BETA_A_4 += stepsizePrev[0][3];
							signPrev[3] = 1;
							if (prevAVG[0][3] < 0.001 && target[48] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[48] - prevAVG[0][3])/target[48]) < 0.0 && signPrev[3] == -1){
							if(fabs((target[48] - prevAVG[0][3])/target[48]) > (2*target[70])){
								BETA_A_4 -= 2*stepsizePrev[0][3];
							}
							else{
								BETA_A_4 -= stepsizePrev[0][3];
							}						signPrev[3] = -1;
						}
						
						else if(((target[48] - prevAVG[0][3])/target[48]) < 0.0 && signPrev[3] == 1){
							stepsizePrev[0][3] = stepsizePrev[0][3]/2.0;//halve the stepsize
							BETA_A_4 -= stepsizePrev[0][3];
							signPrev[3] = -1;
						}
					}
					else if(target[48] != 0.0 && BETA_A_4 == 0.0){
						BETA_A_4 += stepsizePrev[0][3];
					}
					else if(target[48] != 0.0 && prevAVG[0][3] == 0.0){
						BETA_A_4 += stepsizePrev[0][3];
					}
					else if(target[48] == 0.0 && prevAVG[0][3] != 0.0){
						BETA_A_4 = 0.0;
						if (BETA_A_4 < 0.0){
							BETA_A_4 = 0.0;
						}
						
					}

					i = 0; j = 0;
					if(fabs((target[49] - prevAVG[0][4])/target[49]) > target[71] && target[49] != 0.0){
						if (((target[49] - prevAVG[0][4])/target[49]) > 0.0 && signPrev[4] == 1){
							if(fabs((target[49] - prevAVG[0][4])/target[49]) > (2*target[71])){
								BETA_A_5 += 2*stepsizePrev[0][4];
							}
							else{
								BETA_A_5 += stepsizePrev[0][4];
							}
							signPrev[4] = 1;
							if (prevAVG[0][4] < 0.001 && target[49] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[49] - prevAVG[0][4])/target[49]) > 0.0 && signPrev[4] == -1 ){
							stepsizePrev[0][4] = stepsizePrev[0][4]/2.0;//Halve the stepsize
							BETA_A_5 += stepsizePrev[0][4];
							signPrev[4] = 1;
							if (prevAVG[0][4] < 0.001 && target[49] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = 'a';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[49] - prevAVG[0][4])/target[49]) < 0.0 && signPrev[4] == -1){
							if(fabs((target[49] - prevAVG[0][4])/target[49]) > (2*target[71])){
								BETA_A_5 -= 2*stepsizePrev[0][4];
							}
							else{
								BETA_A_5 -= stepsizePrev[0][4];
							}
							
							signPrev[4] = -1;
						}
						
						else if(((target[49] - prevAVG[0][4])/target[49]) < 0.0 && signPrev[4] == 1){
							stepsizePrev[0][4] = stepsizePrev[0][4]/2.0;//halve the stepsize
							BETA_A_5 -= stepsizePrev[0][4];
							signPrev[4] = -1;
						}
					}
					else if(target[49] != 0.0 && BETA_A_5 == 0.0){
						BETA_A_5 += stepsizePrev[0][4];
					}
					else if(target[49] != 0.0 && prevAVG[0][4] == 0.0){
						BETA_A_5 += stepsizePrev[0][4];
					}
					else if(target[49] == 0.0 && prevAVG[0][4] != 0.0){
						BETA_A_5 = 0.0;
						if (BETA_A_5 < 0.0){
							BETA_A_5 = 0.0;
						}
						
					}
					
					//////////////////////////////////////////////////
					// Sero B Prev
					i = 0; j = 0;
					if(fabs((target[50] - prevAVG[1][0])/target[50]) > target[72] && target[50] != 0.0){//check to see if Prev B for 0-4 y.o. are outside target
						if (((target[50] - prevAVG[1][0])/target[50]) > 0.0 && signPrev[5] == 1){
							if(fabs((target[50] - prevAVG[1][0])/target[50]) > (2*target[72])){
								BETA_B_1 += 2*stepsizePrev[1][0];
							}
							else{
								BETA_B_1 += stepsizePrev[1][0];
							}
							
							signPrev[5] = 1;
							if (prevAVG[1][0] < 0.001 && target[50] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[50] - prevAVG[1][0])/target[50]) > 0.0 && signPrev[5] == -1 ){
							stepsizePrev[1][0] = stepsizePrev[1][0]/2.0;
							signPrev[5] = 1;
							BETA_B_1 += stepsizePrev[1][0];
							if (prevAVG[1][0] < 0.001 && target[50] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[50] - prevAVG[1][0])/target[50]) < 0.0 && signPrev[5] == -1){
							if(fabs((target[50] - prevAVG[1][0])/target[50]) > (2*target[72])){
								BETA_B_1 -= 2*stepsizePrev[1][0];
							}
							else{
								BETA_B_1 -= stepsizePrev[1][0];
							}
							
							signPrev[5] = -1;
						}
						
						else if(((target[50] - prevAVG[1][0])/target[50]) < 0.0 && signPrev[5] == 1){
							stepsizePrev[1][0] = stepsizePrev[1][0]/2.0;
							BETA_B_1 -= stepsizePrev[1][0];
							signPrev[5] = -1;
						}
					}
					else if(target[50] != 0.0 && BETA_B_1 == 0.0){
						BETA_B_1 += stepsizePrev[1][0];
					}
					else if(target[50] != 0.0 && prevAVG[1][0] == 0.0){
						BETA_B_1 += stepsizePrev[1][0];
					}
					else if(target[50] == 0.0 && prevAVG[1][0] != 0.0){
						BETA_B_1 -= stepsizePrev[1][0];
						if (BETA_B_1 < 0.0){
							BETA_B_1 = 0.0;
						}
						
					}
					
					////
					
					i = 0; j = 0;
					if(fabs((target[51] - prevAVG[1][1])/target[51]) > target[73] && target[51] != 0.0){
						if (((target[51] - prevAVG[1][1])/target[51]) > 0.0 && signPrev[6] == 1){
							BETA_B_2 += stepsizePrev[1][1];
							signPrev[6] = 1;
							if (prevAVG[1][1] < 0.001 && target[51] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[51] - prevAVG[1][1])/target[51]) > 0.0 && signPrev[6] == -1 ){
							stepsizePrev[1][1] = stepsizePrev[1][1]/2.0;//Halve the stepsize
							BETA_B_2 += stepsizePrev[1][1];
							signPrev[6] = 1;
							if (prevAVG[1][1] < 0.001 && target[51] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[51] - prevAVG[1][1])/target[51]) < 0.0 && signPrev[6] == -1){
							BETA_B_2 -= stepsizePrev[1][1];
							signPrev[6] = -1;
						}
						
						else if(((target[51] - prevAVG[1][1])/target[51]) < 0.0 && signPrev[6] == 1){
							stepsizePrev[1][1] = stepsizePrev[1][1]/2.0;//halve the stepsize
							BETA_B_2 -= stepsizePrev[1][1];
							signPrev[6] = -1;
						}
					}
					else if(target[51] != 0.0 && BETA_B_2 == 0.0){
						BETA_B_2 += stepsizePrev[1][1];
					}
					else if(target[51] != 0.0 && prevAVG[1][1] == 0.0){
						BETA_B_2 += stepsizePrev[1][1];
					}
					else if(target[51] == 0.0 && prevAVG[1][1] != 0.0){
						BETA_B_3 -= stepsizePrev[1][1];
						if (BETA_B_2 < 0.0){
							BETA_B_2 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[52] - prevAVG[1][2])/target[52]) > target[74] && target[52] != 0.0){
						if (((target[52] - prevAVG[1][2])/target[52]) > 0.0 && signPrev[7] == 1){
							BETA_B_3 += stepsizePrev[1][2];
							signPrev[7] = 1;
							if (prevAVG[1][2] < 0.001 && target[52] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[52] - prevAVG[1][2])/target[52]) > 0.0 && signPrev[7] == -1 ){
							stepsizePrev[1][2] = stepsizePrev[1][2]/2.0;//Halve the stepsize
							BETA_B_3 += stepsizePrev[1][2];
							signPrev[7] = 1;
							if (prevAVG[1][2] < 0.001 && target[52] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[52] - prevAVG[1][2])/target[52]) < 0.0 && signPrev[7] == -1){
							BETA_B_3 -= stepsizePrev[1][2];
							signPrev[7] = -1;
						}
						
						else if(((target[52] - prevAVG[1][2])/target[52]) < 0.0 && signPrev[7] == 1){
							stepsizePrev[1][2] = stepsizePrev[1][2]/2.0;//halve the stepsize
							BETA_B_3 -= stepsizePrev[1][2];
							signPrev[7] = -1;
						}
					}
					else if(target[52] != 0.0 && BETA_B_3 == 0.0){
						BETA_B_3 += stepsizePrev[1][2];
					}
					else if(target[52] != 0.0 && prevAVG[1][2] == 0.0){
						BETA_B_3 += stepsizePrev[1][2];
					}
					else if(target[52] == 0.0 && prevAVG[1][2] != 0.0){
						BETA_B_3 -= stepsizePrev[1][2];
						if (BETA_B_3 < 0.0){
							BETA_B_3 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[53] - prevAVG[1][3])/target[53]) > target[75] && target[53] != 0.0){
						if (((target[53] - prevAVG[1][3])/target[53]) > 0.0 && signPrev[8] == 1){
							BETA_B_4 += stepsizePrev[1][3];
							signPrev[8] = 1;
							if (prevAVG[1][3] < 0.001 && target[53] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[53] - prevAVG[1][3])/target[53]) > 0.0 && signPrev[8] == -1 ){
							stepsizePrev[1][3] = stepsizePrev[1][3]/2.0;//Halve the stepsize
							BETA_B_4 += stepsizePrev[1][3];
							signPrev[8] = 1;
							if (prevAVG[1][3] < 0.001 && target[53] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[53] - prevAVG[1][3])/target[53]) < 0.0 && signPrev[8] == -1){
							BETA_B_4 -= stepsizePrev[1][3];
							signPrev[8] = -1;
						}
						
						else if(((target[53] - prevAVG[1][3])/target[53]) < 0.0 && signPrev[8] == 1){
							stepsizePrev[1][3] = stepsizePrev[1][3]/2.0;//halve the stepsize
							BETA_B_4 -= stepsizePrev[1][3];
							signPrev[8] = -1;
						}
					}
					else if(target[53] != 0.0 && BETA_B_4 == 0.0){
						BETA_B_4 += stepsizePrev[1][3];
					}
					else if(target[53] != 0.0 && prevAVG[1][3] == 0.0){
						BETA_B_4 += stepsizePrev[1][3];
					}
					else if(target[53] == 0.0 && prevAVG[1][3] != 0.0){
						BETA_B_4 -= stepsizePrev[1][3];
						if (BETA_B_4 < 0.0){
							BETA_B_4 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[54] - prevAVG[1][4])/target[54]) > target[76] && target[54] != 0.0){
						if (((target[54] - prevAVG[1][4])/target[54]) > 0.0 && signPrev[9] == 1){
							BETA_B_5 += stepsizePrev[1][4];
							signPrev[9] = 1;
							if (prevAVG[1][4] < 0.001 && target[54] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[54] - prevAVG[1][4])/target[54]) > 0.0 && signPrev[9] == -1 ){
							stepsizePrev[1][4] = stepsizePrev[1][4]/2.0;//Halve the stepsize
							BETA_B_5 += stepsizePrev[1][4];
							signPrev[9] = 1;
							if (prevAVG[1][4] < 0.001 && target[54] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = 'b';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[54] - prevAVG[1][4])/target[54]) < 0.0 && signPrev[9] == -1){
							BETA_B_5 -= stepsizePrev[1][4];
							signPrev[9] = -1;
						}
						
						else if(((target[54] - prevAVG[1][4])/target[54]) < 0.0 && signPrev[9] == 1){
							stepsizePrev[1][4] = stepsizePrev[1][4]/2.0;//halve the stepsize
							BETA_B_5 -= stepsizePrev[1][4];
							signPrev[9] = -1;
						}
					}
					else if(target[54] != 0.0 && BETA_B_5 == 0.0){
						BETA_B_5 += stepsizePrev[1][4];
					}
					else if(target[54] != 0.0 && prevAVG[1][4] == 0.0){
						BETA_B_5 += stepsizePrev[1][4];
					}
					else if(target[54] == 0.0 && prevAVG[1][4] != 0.0){
						BETA_B_5 -= stepsizePrev[1][4];
						if (BETA_B_5 < 0.0){
							BETA_B_5 = 0.0;
						}
						
					}
					
					
					/////////////////////////////////////////////
					// Sero C Prev
					
					i = 0; j = 0;
					if(fabs((target[55] - prevAVG[2][0])/target[55]) > target[77] && target[55] != 0.0){//check to see if Prev C for 0-4 y.o. are outside target
						if (((target[55] - prevAVG[2][0])/target[55]) > 0.0 && signPrev[10] == 1){
							if(fabs((target[55] - prevAVG[2][0])/target[55]) > (2*target[77])){
								BETA_C_1 += 2*stepsizePrev[2][0];
							}
							else{
								BETA_C_1 += stepsizePrev[2][0];
							}
							signPrev[10] = 1;
							if (prevAVG[2][0] < 0.001 && target[55] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[55] - prevAVG[2][0])/target[55]) > 0.0 && signPrev[10] == -1 ){
							stepsizePrev[2][0] = stepsizePrev[2][0]/2.0;
							signPrev[10] = 1;
							BETA_C_1 += stepsizePrev[2][0];
							if (prevAVG[2][0] < 0.001 && target[55] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[55] - prevAVG[2][0])/target[55]) < 0.0 && signPrev[10] == -1){
							if(fabs((target[55] - prevAVG[2][0])/target[55]) > (2*target[77])){
								BETA_C_1 -= 2*stepsizePrev[2][0];
							}
							else{
								BETA_C_1 -= stepsizePrev[2][0];
							}
							signPrev[10] = -1;
						}
						
						else if(((target[55] - prevAVG[2][0])/target[55]) < 0.0 && signPrev[10] == 1){
							stepsizePrev[2][0] = stepsizePrev[2][0]/2.0;
							BETA_C_1 -= stepsizePrev[2][0];
							signPrev[10] = -1;
						}
					}
					else if(target[55] != 0.0 && BETA_C_1 == 0.0){
						BETA_C_1 += stepsizePrev[2][0];
					}
					else if(target[55] != 0.0 && prevAVG[2][0] == 0.0){
						BETA_C_1 += stepsizePrev[2][0];
					}
					else if(target[55] == 0.0 && prevAVG[2][0] != 0.0){
						BETA_C_1 -= stepsizePrev[2][0];
						if (BETA_C_1 < 0.0){
							BETA_C_1 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[56] - prevAVG[2][1])/target[56]) > target[78] && target[56] != 0.0){
						if (((target[56] - prevAVG[2][1])/target[56]) > 0.0 && signPrev[11] == 1){
							BETA_C_2 += stepsizePrev[2][1];
							signPrev[11] = 1;
							if (prevAVG[2][1] < 0.001 && target[56] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[56] - prevAVG[2][1])/target[56]) > 0.0 && signPrev[11] == -1 ){
							stepsizePrev[2][1] = stepsizePrev[2][1]/2.0;//Halve the stepsize
							BETA_C_2 += stepsizePrev[2][1];
							signPrev[11] = 1;
							if (prevAVG[2][1] < 0.001 && target[56] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[56] - prevAVG[2][1])/target[56]) < 0.0 && signPrev[11] == -1){
							BETA_C_2 -= stepsizePrev[2][1];
							//AVG_DUR_INFECT_C -= 1.0;
							signPrev[11] = -1;
						}
						
						else if(((target[56] - prevAVG[2][1])/target[56]) < 0.0 && signPrev[11] == 1){
							stepsizePrev[2][1] = stepsizePrev[2][1]/2.0;//halve the stepsize
							BETA_C_2 -= stepsizePrev[2][1];
							//AVG_DUR_INFECT_C -= 1.0;
							signPrev[11] = -1;
						}
					}
					else if(target[56] != 0.0 && BETA_C_2 == 0.0){
						BETA_C_2 += stepsizePrev[2][1];
					}
					else if(target[56] != 0.0 && prevAVG[2][1] == 0.0){
						BETA_C_2 += stepsizePrev[2][1];
					}
					else if(target[56] == 0.0 && prevAVG[2][1] != 0.0){
						BETA_C_2 -= stepsizePrev[2][1];
						if (BETA_C_2 < 0.0){
							BETA_C_2 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[57] - prevAVG[2][2])/target[57]) > target[79] && target[57] != 0.0){
						if (((target[57] - prevAVG[2][2])/target[57]) > 0.0 && signPrev[12] == 1){
							BETA_C_3 += stepsizePrev[2][2];
							//AVG_DUR_INFECT_C += 1.0;
							signPrev[12] = 1;
							if (prevAVG[2][2] < 0.001 && target[57] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[57] - prevAVG[2][2])/target[57]) > 0.0 && signPrev[12] == -1 ){
							stepsizePrev[2][2] = stepsizePrev[2][2]/2.0;//Halve the stepsize
							BETA_C_3 += stepsizePrev[2][2];
							//AVG_DUR_INFECT_C += 1.0;
							signPrev[12] = 1;
							if (prevAVG[2][2] < 0.001 && target[57] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[57] - prevAVG[2][2])/target[57]) < 0.0 && signPrev[12] == -1){
							BETA_C_3 -= stepsizePrev[2][2];
							//AVG_DUR_INFECT_C -= 1.0;
							signPrev[12] = -1;
						}
						
						else if(((target[57] - prevAVG[2][2])/target[57]) < 0.0 && signPrev[12] == 1){
							stepsizePrev[2][2] = stepsizePrev[2][2]/2.0;//halve the stepsize
							BETA_C_3 -= stepsizePrev[2][2];
							//AVG_DUR_INFECT_C -= 1.0;
							signPrev[12] = -1;
						}
					}
					else if(target[57] != 0.0 && BETA_C_3 == 0.0){
						BETA_C_3 += stepsizePrev[2][2];
					}
					else if(target[57] != 0.0 && prevAVG[2][2] == 0.0){
						BETA_C_3 += stepsizePrev[2][2];
					}
					else if(target[57] == 0.0 && prevAVG[2][2] != 0.0){
						BETA_C_3 -= stepsizePrev[2][2];
						if (BETA_C_3 < 0.0){
							BETA_C_3 = 0.0;
						}
						
					}
					
					////
					
					i = 0; j = 0;
					if(fabs((target[58] - prevAVG[2][3])/target[58]) > target[80] && target[58] != 0.0){
						if (((target[58] - prevAVG[2][3])/target[58]) > 0.0 && signPrev[13] == 1){
							BETA_C_4 += stepsizePrev[2][3];
							signPrev[13] = 1;
							if (prevAVG[2][3] < 0.001 && target[58] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[58] - prevAVG[2][3])/target[58]) > 0.0 && signPrev[13] == -1 ){
							stepsizePrev[2][3] = stepsizePrev[2][3]/2.0;//Halve the stepsize
							BETA_C_4 += stepsizePrev[2][3];
							signPrev[13] = 1;
							if (prevAVG[2][3] < 0.001 && target[58] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[58] - prevAVG[2][3])/target[58]) < 0.0 && signPrev[13] == -1){
							BETA_C_4 -= stepsizePrev[2][3];
							signPrev[13] = -1;
						}
						
						else if(((target[58] - prevAVG[2][3])/target[58]) < 0.0 && signPrev[13] == 1){
							stepsizePrev[2][3] = stepsizePrev[2][3]/2.0;//halve the stepsize
							BETA_C_4 -= stepsizePrev[2][3];
							signPrev[13] = -1;
						}
					}
					else if(target[58] != 0.0 && BETA_C_4 == 0.0){
						BETA_C_4 += stepsizePrev[2][3];
					}
					else if(target[58] != 0.0 && prevAVG[2][3] == 0.0){
						BETA_C_4 += stepsizePrev[2][3];
					}
					else if(target[58] == 0.0 && prevAVG[2][3] != 0.0){
						BETA_C_4 -= stepsizePrev[2][3];
						if (BETA_C_4 < 0.0){
							BETA_C_4 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[59] - prevAVG[2][4])/target[59]) > target[81] && target[59] != 0.0){
						if (((target[59] - prevAVG[2][4])/target[59]) > 0.0 && signPrev[14] == 1){
							BETA_C_5 += stepsizePrev[2][4];
							signPrev[14] = 1;
							if (prevAVG[2][4] < 0.001 && target[59] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[59] - prevAVG[2][4])/target[59]) > 0.0 && signPrev[14] == -1 ){
							stepsizePrev[2][4] = stepsizePrev[2][4]/2.0;//Halve the stepsize
							BETA_C_5 += stepsizePrev[2][4];
							signPrev[14] = 1;
							if (prevAVG[2][4] < 0.001 && target[59] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = 'c';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[59] - prevAVG[2][4])/target[59]) < 0.0 && signPrev[14] == -1){
							BETA_C_5 -= stepsizePrev[2][4];
							signPrev[14] = -1;
						}
						
						else if(((target[59] - prevAVG[2][4])/target[59]) < 0.0 && signPrev[14] == 1){
							stepsizePrev[2][4] = stepsizePrev[2][4]/2.0;//halve the stepsize
							BETA_C_5 -= stepsizePrev[2][4];
							signPrev[14] = -1;
						}
					}
					else if(target[59] != 0.0 && BETA_C_5 == 0.0){
						BETA_C_5 += stepsizePrev[2][4];
					}
					else if(target[59] != 0.0 && prevAVG[2][4] == 0.0){
						BETA_C_5 += stepsizePrev[2][4];
					}
					else if(target[59] == 0.0 && prevAVG[2][4] != 0.0){
						BETA_C_5 -= stepsizePrev[2][4];
						if (BETA_C_5 < 0.0){
							BETA_C_5 = 0.0;
						}
						
					}
					
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//Sero Nl Prev
					i = 0; j = 0;
					if(fabs((target[60] - prevAVG[3][0])/target[60]) > target[82] && target[60] != 0.0){//check to see if Prev NL for 0-4 y.o. are outside target
						if (((target[60] - prevAVG[3][0])/target[60]) > 0.0 && signPrev[15] == 1){
							if(fabs((target[60] - prevAVG[3][0])/target[60]) > (2*target[82])){
								BETA_NL_1 += 2*stepsizePrev[3][0];
							}
							else{
								BETA_NL_1 += stepsizePrev[3][0];
							}
							signPrev[30] = 1;
							if (prevAVG[3][0] < 0.001 && target[60] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[60] - prevAVG[3][0])/target[60]) > 0.0 && signPrev[15] == -1 ){
							stepsizePrev[3][0] = stepsizePrev[3][0]/2.0;
							signPrev[15] = 1;
							BETA_NL_1 += stepsizePrev[3][0];
							if (prevAVG[3][0] < 0.001 && target[60] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=0.00 && population[i].age < 5.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[60] - prevAVG[3][0])/target[60]) < 0.0 && signPrev[15] == -1){
							if(fabs((target[60] - prevAVG[3][0])/target[60]) > (2*target[82])){
								BETA_NL_1 -= 2*stepsizePrev[3][0];
							}
							else{
								BETA_NL_1 -= stepsizePrev[3][0];
							}
							
							signPrev[15] = -1;
						}
						
						else if(((target[60] - prevAVG[3][0])/target[60]) < 0.0 && signPrev[15] == 1){
							stepsizePrev[3][0] = stepsizePrev[3][0]/2.0;
							BETA_NL_1 -= stepsizePrev[3][0];
							signPrev[15] = -1;
						}
					}
					else if(target[60] != 0.0 && BETA_NL_1 == 0.0){
						BETA_NL_1 += stepsizePrev[3][0];
					}
					else if(target[60] != 0.0 && prevAVG[3][0] == 0.0){
						BETA_NL_1 += stepsizePrev[3][0];
					}
					else if(target[60] == 0.0 && prevAVG[3][0] != 0.0){
						BETA_NL_1 -= stepsizePrev[3][0];
						if (BETA_NL_1 < 0.0){
							BETA_NL_1 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[61] - prevAVG[3][1])/target[61]) > target[83] && target[61] != 0.0){
						if (((target[61] - prevAVG[3][1])/target[61]) > 0.0 && signPrev[16] == 1){
							BETA_NL_2 += stepsizePrev[3][1];
							signPrev[16] = 1;
							if (prevAVG[3][1] < 0.001 && target[61] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[61] - prevAVG[3][1])/target[61]) > 0.0 && signPrev[16] == -1 ){
							stepsizePrev[3][1] = stepsizePrev[3][1]/2.0;//Halve the stepsize
							BETA_NL_2 += stepsizePrev[3][1];
							signPrev[16] = 1;
							if (prevAVG[3][1] < 0.001 && target[61] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=5.00 && population[i].age < 10.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[61] - prevAVG[3][1])/target[61]) < 0.0 && signPrev[16] == -1){
							BETA_NL_2 -= stepsizePrev[3][1];
							signPrev[16] = -1;
						}
						
						else if(((target[61] - prevAVG[6][1])/target[61]) < 0.0 && signPrev[16] == 1){
							stepsizePrev[3][1] = stepsizePrev[3][1]/2.0;//halve the stepsize
							BETA_NL_2 -= stepsizePrev[6][1];
							signPrev[16] = -1;
						}
					}
					else if(target[61] != 0.0 && BETA_NL_2 == 0.0){
						BETA_NL_2 += stepsizePrev[3][1];
					}
					else if(target[61] != 0.0 && prevAVG[6][1] == 0.0){
						BETA_NL_2 += stepsizePrev[3][1];
					}
					else if(target[61] == 0.0 && prevAVG[6][1] != 0.0){
						BETA_NL_2 -= stepsizePrev[3][1];
						if (BETA_NL_2 < 0.0){
							BETA_NL_2 = 0.0;
						}
						
					}
					
					i = 0; j = 0;
					if(fabs((target[62] - prevAVG[3][2])/target[62]) > target[84] && target[62] != 0.0){
						if (((target[62] - prevAVG[3][2])/target[62]) > 0.0 && signPrev[17] == 1){
							BETA_NL_3 += stepsizePrev[3][2];
							signPrev[17] = 1;
							if (prevAVG[3][2] < 0.001 && target[62] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[62] - prevAVG[3][2])/target[62]) > 0.0 && signPrev[32] == -1 ){
							stepsizePrev[6][2] = stepsizePrev[3][2]/2.0;//Halve the stepsize
							BETA_NL_3 += stepsizePrev[3][2];
							signPrev[17] = 1;
							if (prevAVG[3][2] < 0.001 && target[62] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=10.00 && population[i].age < 15.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[62] - prevAVG[3][2])/target[62]) < 0.0 && signPrev[17] == -1){
							BETA_NL_3 -= stepsizePrev[3][2];
							signPrev[17] = -1;
						}
						
						else if(((target[62] - prevAVG[3][2])/target[62]) < 0.0 && signPrev[17] == 1){
							stepsizePrev[3][2] = stepsizePrev[3][2]/2.0;//halve the stepsize
							BETA_NL_3 -= stepsizePrev[3][2];
							signPrev[17] = -1;
						}
					}
					else if(target[62] != 0.0 && BETA_NL_3 == 0.0){
						BETA_NL_3 += stepsizePrev[3][2];
					}
					else if(target[62] != 0.0 && prevAVG[3][2] == 0.0){
						BETA_NL_3 += stepsizePrev[3][2];
					}
					else if(target[62] == 0.0 && prevAVG[3][2] != 0.0){
						BETA_NL_3 -= stepsizePrev[3][2];
						if (BETA_NL_3 < 0.0){
							BETA_NL_3 = 0.0;
						}
						
					}

					i = 0; j = 0;
					if(fabs((target[63] - prevAVG[3][3])/target[63]) > target[85] && target[63] != 0.0){
						if (((target[63] - prevAVG[3][3])/target[63]) > 0.0 && signPrev[18] == 1){
							BETA_NL_4 += stepsizePrev[3][3];
							signPrev[18] = 1;
							if (prevAVG[3][3] < 0.001 && target[63] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[63] - prevAVG[3][3])/target[63]) > 0.0 && signPrev[18] == -1 ){
							stepsizePrev[3][3] = stepsizePrev[3][3]/2.0;//Halve the stepsize
							BETA_NL_4 += stepsizePrev[3][3];
							signPrev[18] = 1;
							if (prevAVG[3][3] < 0.001 && target[63] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >=15.00 && population[i].age < 20.0)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[63] - prevAVG[3][3])/target[63]) < 0.0 && signPrev[18] == -1){
							BETA_NL_4 -= stepsizePrev[3][3];
							signPrev[18] = -1;
						}
						
						else if(((target[63] - prevAVG[3][3])/target[63]) < 0.0 && signPrev[18] == 1){
							stepsizePrev[3][3] = stepsizePrev[3][3]/2.0;//halve the stepsize
							BETA_NL_4 -= stepsizePrev[3][3];
							signPrev[18] = -1;
						}
					}
					else if(target[63] != 0.0 && BETA_NL_4 == 0.0){
						BETA_NL_4 += stepsizePrev[3][3];
					}
					else if(target[63] != 0.0 && prevAVG[3][3] == 0.0){
						BETA_NL_4 += stepsizePrev[3][3];
					}
					else if(target[63] == 0.0 && prevAVG[3][3] != 0.0){
						BETA_NL_4 -= stepsizePrev[3][3];
						if (BETA_NL_4 < 0.0){
							BETA_NL_4 = 0.0;
						}
						
					}
					//
					
					i = 0; j = 0;
					if(fabs((target[64] - prevAVG[3][4])/target[64]) > target[86] && target[64] != 0.0){
						if (((target[64] - prevAVG[3][4])/target[64]) > 0.0 && signPrev[19] == 1){
							BETA_NL_5 += stepsizePrev[3][4];
							signPrev[19] = 1;
							if (prevAVG[3][4] < 0.001 && target[64] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[64] - prevAVG[6][4])/target[64]) > 0.0 && signPrev[19] == -1 ){
							stepsizePrev[3][4] = stepsizePrev[6][4]/2.0;//Halve the stepsize
							BETA_NL_5 += stepsizePrev[6][4];
							signPrev[19] = 1;
							if (prevAVG[3][4] < 0.001 && target[64] != 0.0) {
								i = 0;	j = 15;
								while (i < ARRAY_SIZE && j >0){
									if(population[i].status == 0 && (population[i].age >= 20.00)){
										infection(&population[i], t);
										population[i].serogroup = '2';
										j--;
									}
									i++;
								}
							}
						}
						else if(((target[64] - prevAVG[3][4])/target[64]) < 0.0 && signPrev[19] == -1){
							BETA_NL_5 -= stepsizePrev[6][4];
							signPrev[19] = -1;
						}
						
						else if(((target[64] - prevAVG[3][4])/target[64]) < 0.0 && signPrev[19] == 1){
							stepsizePrev[3][4] = stepsizePrev[3][4]/2.0;//halve the stepsize
							BETA_NL_5 -= stepsizePrev[3][4];
							signPrev[19] = -1;
						}
					}
					else if(target[64] != 0.0 && BETA_NL_5 == 0.0){
						BETA_NL_5 += stepsizePrev[3][4];
					}
					else if(target[64] != 0.0 && prevAVG[3][4] == 0.0){
						BETA_NL_5 += stepsizePrev[3][4];
					}
					else if(target[64] == 0.0 && prevAVG[3][4] != 0.0){
						BETA_NL_5 -= stepsizePrev[3][4];
						if (BETA_NL_5 < 0.0){
							BETA_NL_5 = 0.0;
						}
					}
				}
			}
			//Prevents BETA from being negative
			if(BETA_A_1 < 0.0)
            {BETA_A_1 = 0.0;}
			if(BETA_A_2 < 0.0)
            {BETA_A_2 = 0.0;}
			if(BETA_A_3 < 0.0)
            {BETA_A_3 = 0.0;}
			if(BETA_A_4 < 0.0)
            {BETA_A_4 = 0.0;}
			if(BETA_A_5 < 0.0)
            {BETA_A_5 = 0.0;}
			if(BETA_B_1 < 0.0)
            {BETA_B_1 = 0.0;}
			if(BETA_B_2 < 0.0)
            {BETA_B_2 = 0.0;}
			if(BETA_B_3 < 0.0)
            {BETA_B_3 = 0.0;}
			if(BETA_B_4 < 0.0)
            {BETA_B_4 = 0.0;}
			if(BETA_B_5 < 0.0)
            {BETA_B_5 = 0.0;}
			if(BETA_C_1 < 0.0)
            {BETA_C_1 = 0.0;}
			if(BETA_C_2 < 0.0)
            {BETA_C_2 = 0.0;}
			if(BETA_C_3 < 0.0)
            {BETA_C_3 = 0.0;}
			if(BETA_C_4 < 0.0)
            {BETA_C_4 = 0.0;}
			if(BETA_C_5 < 0.0)
            {BETA_C_5 = 0.0;}
			if(BETA_NL_1 < 0.0)
            {BETA_NL_1 = 0.0;}
			if(BETA_NL_2 < 0.0)
            {BETA_NL_2 = 0.0;}
			if(BETA_NL_3 < 0.0)
            {BETA_NL_3 = 0.0;}
			if(BETA_NL_4 < 0.0)
            {BETA_NL_4 = 0.0;}
			if(BETA_NL_5 < 0.0)
            {BETA_NL_5 = 0.0;}
		}
		/////////////////////////////////////////////
		// Display update to user
        if(flagOpt == 1){
            if(((int)(t)%1200 == 1199)){
				
                printf("Time %.1lf \t Number of parameters within target: %d\t Population size: %d\n", t/12.0, itemp, POP_SIZE);
                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
                printf("|\t Prevalence Targets \t\t\t\t\t\t Sim run \t\t\t\t\t\t|\n");
                printf("|\t A \t B \t C \t Nl\t \t A \t B \t C \t Nl\t|\n");
                printf("|\t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t \t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t|\n", target[45], target[50], target[55], target[60], prevAVG[0][0], prevAVG[1][0], prevAVG[2][0], prevAVG[3][0]);
                printf("|\t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t \t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t|\n", target[46], target[51], target[56], target[61], prevAVG[0][1], prevAVG[1][1], prevAVG[2][1], prevAVG[3][1]);
                printf("|\t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t \t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t|\n", target[47], target[52], target[57], target[62], prevAVG[0][2], prevAVG[1][2], prevAVG[2][2], prevAVG[3][2]);
                printf("|\t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t \t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t|\n", target[48], target[53], target[58], target[63], prevAVG[0][3], prevAVG[1][3], prevAVG[2][3], prevAVG[3][3]);
                printf("|\t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t \t %0.4lf\t %0.4lf\t %0.4lf\t %0.4lf\t|\n", target[49], target[54], target[59], target[64], prevAVG[0][4], prevAVG[1][4], prevAVG[2][4], prevAVG[3][4]);
                printf("|\t BETA VALUES\t\t\t\t\t\t\t\t\t\t\t\t\t\t|\n");
                printf("|\t A \t\t B \t\t C \t\t Nl\t\t  | \n");
                printf("|\t %0.10lf\t %0.10lf\t %0.10lf\t %0.10lf\t \t |\n", BETA_A_1, BETA_B_1, BETA_C_1, BETA_NL_1);
                printf("|\t %0.10lf\t %0.10lf\t %0.10lf\t %0.10lf\t \t |\n", BETA_A_2, BETA_B_2, BETA_C_2, BETA_NL_2);
                printf("|\t %0.10lf\t %0.10lf\t %0.10lf\t %0.10lf\t \t |\n", BETA_A_3, BETA_B_3, BETA_C_3, BETA_NL_3);
                printf("|\t %0.10lf\t %0.10lf\t %0.10lf\t %0.10lf\t \t |\n", BETA_A_4, BETA_B_4, BETA_C_4, BETA_NL_4);
                printf("|\t %0.10lf\t %0.10lf\t %0.10lf\t %0.10lf\t \t |\n", BETA_A_5, BETA_B_5, BETA_C_5, BETA_NL_5);
                printf("|\t Stepsize\t\t\t\t\t\t\t\t\t\t\t\t\t\t |\n");
                printf("|\t A \t\t B \t\t C \t\t Nl\t\t  | \n");
                printf("|\t %0.8lf\t %0.8lf\t %0.8lf\t %0.8lf\t \t |\n", stepsizePrev[0][0], stepsizePrev[1][0], stepsizePrev[2][0], stepsizePrev[3][0]);
                printf("|\t %0.8lf\t %0.8lf\t %0.8lf\t %0.8lf\t \t |\n", stepsizePrev[0][1], stepsizePrev[1][1], stepsizePrev[2][1], stepsizePrev[3][1]);
                printf("|\t %0.8lf\t %0.8lf\t %0.8lf\t %0.8lf\t \t |\n", stepsizePrev[0][2], stepsizePrev[1][2], stepsizePrev[2][2], stepsizePrev[3][2]);
                printf("|\t %0.8lf\t %0.8lf\t %0.8lf\t %0.8lf\t \t |\n", stepsizePrev[0][3], stepsizePrev[1][3], stepsizePrev[2][3], stepsizePrev[3][3]);
                printf("|\t %0.8lf\t %0.8lf\t %0.8lf\t %0.8lf\t \t |\n", stepsizePrev[0][4], stepsizePrev[1][4], stepsizePrev[2][4], stepsizePrev[3][4]);
                printf("|\t %lf\t %lf\t %lf\t %lf\t \t\n", AVG_DUR_INFECT_A, AVG_DUR_INFECT_B, AVG_DUR_INFECT_C, AVG_DUR_INFECT_NL);
                printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
				
                fprintf(flag1_beta, "%.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf, %.12lf\n", BETA_A_1, BETA_A_2, BETA_A_3, BETA_A_4, BETA_A_5,BETA_B_1, BETA_B_2, BETA_B_3, BETA_B_4, BETA_B_5, BETA_C_1, BETA_C_2, BETA_C_3, BETA_C_4, BETA_C_5, BETA_NL_1, BETA_NL_2, BETA_NL_3, BETA_NL_4, BETA_NL_5);		
            }
        }

		////////////////////////////
		// Write data to file
        if(t >= 0.0 && flagOpt <= 4){
            fprintf(fout,"%.2f,%d,%f,%f,%f,%f \n", t/12.0, POP_SIZE/POP_SIZE, (double)(numSusc)/(double)(POP_SIZE), (double)(numInfect)/(double)(POP_SIZE), (double)(IMDPerYear)/(double)(POP_SIZE), (double)(numRecover)/(double)(POP_SIZE));
			
            /////////////////////////
            //output infections by age
            fprintf(fage,"%.2f,",t/12.0);
            for(j=0;j<101;j++){
                fprintf(fage," %f, ", (double)(infectByAge[j])/(double)(POP_SIZE));
            }
			
            if((int)(round(t))%12 == 0){
                IMDPerYear = 0;
            }
            //printf("\n\n");
            fprintf(fage,"\n");
			
            /////////////////////////////////////////
            // Demographic Update
            fprintf(fdemograph, "%.2lf", t/12.0);
			for(i=0;i<100;i++){
				fprintf(fdemograph, ",%d", numAgentPop[i]);
			}
			fprintf(fdemograph, "\n");
            
			//Condensed demographics
            fprintf(fConDemo, "%.2lf,%d,%d,%d,%d,%d\n", t/12.0, condenseAge[0], condenseAge[1], condenseAge[2], condenseAge[3], condenseAge[4]);
			
            /////////////////////////////////////////////////
            //Write to serogroup file  'other, b, c, nl
			
            if(flagOpt == 4){
                fprintf(fprevA, "%.2lf,%lf,%lf,%lf,%lf,%lf\n",t/12.0, condensePrev[0][0]/(double)(condenseAge[0]), condensePrev[0][1]/(double)(condenseAge[1]), condensePrev[0][2]/(double)(condenseAge[2]), condensePrev[0][3]/(double)(condenseAge[3]), condensePrev[0][4]/(double)(condenseAge[4]));
                fprintf(fprevB, "%.2lf,%lf,%lf,%lf,%lf,%lf\n",t/12.0, condensePrev[1][0]/(double)(condenseAge[0]), condensePrev[1][1]/(double)(condenseAge[1]), condensePrev[1][2]/(double)(condenseAge[2]), condensePrev[1][3]/(double)(condenseAge[3]), condensePrev[1][4]/(double)(condenseAge[4]));
                fprintf(fprevC, "%.2lf,%lf,%lf,%lf,%lf,%lf\n",t/12.0, condensePrev[2][0]/(double)(condenseAge[0]), condensePrev[2][1]/(double)(condenseAge[1]), condensePrev[2][2]/(double)(condenseAge[2]), condensePrev[2][3]/(double)(condenseAge[3]), condensePrev[2][4]/(double)(condenseAge[4]));
                fprintf(fprevNl,"%.2lf,%lf,%lf,%lf,%lf,%lf\n",t/12.0, condensePrev[3][0]/(double)(condenseAge[0]), condensePrev[3][1]/(double)(condenseAge[1]), condensePrev[3][2]/(double)(condenseAge[2]), condensePrev[3][3]/(double)(condenseAge[3]), condensePrev[3][4]/(double)(condenseAge[4]));
            }
			
			fprintf(fsero,"%.2f,%lf,%lf,%lf,%lf\n", t/12.0, infectBySero[0]/POP_SIZE, infectBySero[1]/POP_SIZE, infectBySero[2]/POP_SIZE, infectBySero[3]/POP_SIZE);
			
            infectBySero[0] = prevAVG[0][0]*condenseAge[0]+prevAVG[0][1]*condenseAge[1]+prevAVG[0][2]*condenseAge[2]+prevAVG[0][3]*condenseAge[3]+prevAVG[0][4]*condenseAge[4];
            infectBySero[1] = prevAVG[1][0]*condenseAge[0]+prevAVG[1][1]*condenseAge[1]+prevAVG[1][2]*condenseAge[2]+prevAVG[1][3]*condenseAge[3]+prevAVG[1][4]*condenseAge[4];
            infectBySero[2] = prevAVG[2][0]*condenseAge[0]+prevAVG[2][1]*condenseAge[1]+prevAVG[2][2]*condenseAge[2]+prevAVG[2][3]*condenseAge[3]+prevAVG[2][4]*condenseAge[4];
            infectBySero[3] = prevAVG[3][0]*condenseAge[0]+prevAVG[3][1]*condenseAge[1]+prevAVG[3][2]*condenseAge[2]+prevAVG[3][3]*condenseAge[3]+prevAVG[3][4]*condenseAge[4];
			
            fprintf(fAVGsero,"%.2f,%lf,%lf,%lf,%lf\n",  t/12.0, infectBySero[0]/POP_SIZE, infectBySero[1]/POP_SIZE, infectBySero[2]/POP_SIZE, infectBySero[3]/POP_SIZE);
        }
		
		//Update counters
        for(i=0;i<ARRAY_SIZE;i++){
            if(population[i].status == 1){
                if(population[i].serogroup == 'a'){
                    infectBySero[0] = infectBySero[0] + 1.0;
                }
                else if(population[i].serogroup == 'b'){
                    infectBySero[1] = infectBySero[1] + 1.0;
                }
                else if(population[i].serogroup == 'c'){
                    infectBySero[2] = infectBySero[2] + 1.0;
                }
                else{
                    infectBySero[3] = infectBySero[3] + 1.0;
                }
            }
        }
		
        t = t + dt;
		
    }
	//close all open files
    fclose(fdemograph);
	fclose(fout);
	fclose(fage);
	fclose(fsero);
	free(monoval);
    fclose(flag1_beta);
	
    return 0;
}

// A function to produce a random time accorinding to a gamma distribution
float
phaseTimer (double average, double sigma)
{
    float phaseTime;
    phaseTime = gammaDist(average, sigma);
    return phaseTime;
}

// Determine initially infectious individuals' serogroup
char
probSerogroup(pop *population){
    //use a case switch to perform the vaccine policy. AKA if vaccineSerogroup == b, decrease chance of b
    double probA, probB, probC, dtemp = 0.0;
    char sero = population->serogroup;
    unsigned int iseed;
    probA = 0.25; probB = 0.25; probC = 0.25;
    while (sero == population->serogroup){
        iseed = (unsigned int)rand()*time(NULL);
        srand(iseed);
        dtemp = 1.0*rand()/RAND_MAX;
		
        if(population->vaccine == 0){
            if(dtemp >= 0 && dtemp < probA){
                sero = 'a';
            }
            else if(dtemp >= probA && dtemp < (probA+probB)){
                sero = 'b';
            }
            else if(dtemp >= (probA+probB) && dtemp < (probA+probB+probC)){
                sero = 'c';
            }
            else{
                sero = '2';
            }
        }
    }
    return sero;
}

// A function to produce a random number to compare to the probability of becoming infectious 
double
coinFlip (){
    unsigned int iseed = (unsigned int)time(NULL)*(unsigned int)(rand());
    srand (iseed);
	
    return ((double)(rand())/(double)(RAND_MAX));
}

// A funtion to compare 2 numbers to determine if individual becomes infectious
int
ranCheck (double prob, double ranNum){
	
    if(ranNum <= prob)
        return 1;
    else
        return 0;
}

//Gamma distribution random number generator to determine duration of infection and duration of immunity
float
gammaDist (double a, double b){
    const gsl_rng_type * T;
    gsl_rng * r;
	
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    float k = gsl_ran_gamma (r, a, b);
	
    gsl_rng_free (r);
    return k;
}

//A function to change the status of an individual to become infectious
void
infection(pop *population, float curTime){
    population->status = 1;
    if(population->serogroup == 'a'){
        population->infectTime = phaseTimer(AVG_DUR_INFECT_A, 1.30)+ curTime;
    }
    else if(population->serogroup == 'b'){
        population->infectTime = phaseTimer(AVG_DUR_INFECT_B, 1.30)+ curTime;
    }
    else if(population->serogroup == 'c'){
        population->infectTime = phaseTimer(AVG_DUR_INFECT_C, 1.30)+ curTime;
    }
    else{
        population->infectTime = phaseTimer(AVG_DUR_INFECT_NL, 1.0) + curTime;
    }
    population->recoverTime = 0.0;
}

// A function to apply vaccine properties to an individual
void
applyVaccine(pop *agent, vaccineProgram *vaccine, double currTime)
{
    agent->vaccine = 1;
    agent->vaccineEfficacy = 1.0;
    agent->vaccineTime = vaccine->lengthVaccinated+currTime;
    agent->vaccineType = vaccine->serogroupProtection;
}

// A function to determine which serogroup infects individuals in what order.
// The order is determined by number of currently infectious individuals with one type of serogroup.
// Most infectious serogroup infects individuals first.  
void
orderSero (int numA, int numB, int numC, int numNl, int numInfect)
{
    double tempA, tempB, tempC, tempNl;
    int i;
    tempA = (double)(numA)/(double)(numInfect); 	tempB = (double)(numB)/(double)(numInfect);		tempC = (double)(numC)/(double)(numInfect);
    tempNl = (double)(numNl)/(double)(numInfect);
    double array[7] = {tempA, tempB, tempC, tempNl};
    gsl_heapsort(array, 7, sizeof(double), compare);
    //printf("Flag 3b passed\n");
    for(i=0;i<7;i++){
        if(array[i] == tempA && SEROORDER[i-1] != 0)
        {SEROORDER[i] = 0;}
        else if(array[i] == tempB && SEROORDER[i-1] != 1)
        {SEROORDER[i] = 1;}
        else if(array[i] == tempC && SEROORDER[i-1] != 2)
        {SEROORDER[i] = 2;}
        else if(array[i] == tempNl && SEROORDER[i-1] != 3)
        {SEROORDER[i] = 3;}
    }
}

// A comparison funtion used in orderSero
int
compare(const double* a, const double* b)
{
	
    if (*a > *b)
        return -1;
    else if (*a < *b)
        return 1;
    else
        return 0;
}

// a function to calculate the running average
double
runningAVG (int size, double *a)
{
    double sum = 0.0, output = 0.0;
    int count;
    //printf("Previous data: ");
    for(count=0;count <= size;count++)
    {
        //printf("%0.5lf, ", a[count]);
        sum+=a[count];
    }
	
    output=sum/(double)(size);
    //printf("\t output: %0.5lf", output);
    return output;
}

// Randomize the order of individuals to be checked for infection
void
shuffle(int *index){
    unsigned long int iseed;
    const gsl_rng_type * T;
    gsl_rng * r;
    iseed = (unsigned long int)rand()*time(NULL);
	
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
	
    gsl_rng_set(r, iseed);
	
    gsl_ran_shuffle (r, index, ARRAY_SIZE, sizeof(int));
	
    gsl_rng_free (r);
}

