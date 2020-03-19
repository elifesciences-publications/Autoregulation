/* 
The following code is to find the fold change of the expression of a negatively
autoregulated gene and its target in presence of competing decoy binding sites.

Decoy = # of competing sites
Nm1 (Nm2) = # of mRNA for TF (target)
Np1 (Np2) = # of protein monomer for TF (target)
Nd = # of TF dimer
comp1, comp2, CDecoy = # of TF bound TF promoter, TF bound target promoter, TF bond competing sites
ku1, ku2, kud = TF unbinding rate for TF promoter, target promoter and competing sites
kb1, kb2, kbd = TF binding rate for TF promoter, target promoter and competing sites
rp1, rp2 = degradation rate for TF monomner, dimer and target protein
gm1, gm2 = transcription rate for TF and target
gp1, gp2 = translation rate for TF and target
rm1, rm2 = mRNA degradation rate for TF and target
Xst, Yst = Constituive expression form TF and target gene
*/ 



#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "time.h"

/* Memory allocation*/

double* dvmalloc(int n);
double* dvmalloc(int n) {
  double *vector;
  vector = (double *)malloc(sizeof(double)*n);
  for (int i=0;i<n;i++) {
    vector[i] = 0;
  }
  return vector;
}

int* ivmalloc(int n) {
  int *vector;
  vector = (int *)malloc(sizeof(int)*n);
  for (int i=0;i<n;i++) {
    vector[i] = 0;
  }
  return vector;
}

/****************************/
/* Find cumulative sum */

double* cumsum(double *A, int num);
double* cumsum(double *A, int num)
{
	double *CSum;
	CSum = dvmalloc(num);
	CSum[0] = A[0];
	for (int i=1;i<num;i++){
		CSum[i] = A[i] + CSum[i-1];
//		printf("%f %d\n",CSum[i],i);
	}
	return CSum;
}
/***********************/
/* Find index */

int FindInd(double *A, double flag);
int FindInd(double *A, double flag)
{
	int check = 0;
	int n=0, Ind;
	while(check==0){
		if (A[n]>=flag){
			check = 1;
			Ind = n;
//		printf("Psum %f Index %d\n",A[n],Ind);
		}
		n++;
		
	}
	return Ind;
}

/* Find a random floating number between min and max */

double unifrnd(double min, double max);
double unifrnd(double min, double max)
{
	double x1;
	/* x1 will be an element of [min,max] */
	x1=((double)rand()/RAND_MAX)*(max-min) + min;
	return x1;
}


/* FInd a random number from exponential distribution with mean K*/
/* x = -ln(1-y)/L, where y is uniform random number between 0 and 1*/
/* K = L */

double exprnd(double K);
double exprnd(double K)
{
	double x,y;
	y = unifrnd(0,0.999);
	x = -log(1-y)/K;
	return x;
}

/* Find a random floating number from a discrete gaussian distribution */

int gaussrnd(int mu, int sig);
int gaussrnd(int mu, int sig)
{
	double x,y, ss = 0;
	double CDF[2*mu+1];
	y = unifrnd(0,0.9999);

	for (int i=0;i<=2*mu;i++){
		ss = ss + exp(-pow(mu-i,2)/(2*sig*sig)) / sqrt(2*M_PI*sig*sig);
		CDF[i] = ss; // Cumulative probability
	}

	int check = 0;
	int n=0, Ind;
//	while(check==0){
	for (int i=0;i<=2*mu;i++){
		if (CDF[n] < y){
			check = 1;
			Ind = n+1;
		}
		n++;	
	}
//	printf("Done step1 = %d\n",Ind);
	return Ind;	
}


/****************************/

int* Gillespie(int Decoy, double tend, int Nm1, int Np1, int comp1, int Nm2, int Np2, int comp2, int CDecoy, int Nd, double ku1, double ku2, double kud, double kb, double rp1, double rp2, double gm1, double gm2, double gp1, double gp2, double rm1, double Xst, double Yst);
int* Gillespie(int Decoy, double tend, int Nm1, int Np1, int comp1, int Nm2, int Np2, int comp2, int CDecoy, int Nd, double ku1, double ku2, double kud, double kb, double rp1, double rp2, double gm1, double gm2, double gp1, double gp2, double rm1, double Xst, double Yst)
//int* Gillespie(int gene1, int gene2, int Decoy, double tend)
{	
// Rates
	int gene1 = 1, gene2 = 1; // gene copy number, kept one throughout
	double kb2 = kb, kb1 = kb, kbd = kb; // Binding
	//double kud = 0.0024; // unbinding from decoy site	
	//double gm1 = 0.44; // mRNA production TF gene
	//double gm2 = 0.44; // mRNA production target gene
	//double gp1 = 0.01; // Protein production
	//double gp2 = 0.01; // Protein production
	//double rp1 = 0.00005; //Degradation of protein
	//double rm1 = 0.0055; // Degradation of RNA 0.011
	
	double rd = 1.38, rm = 0.001/60; // Rate for dimerization/ monomerization 1.38 and 0.001/60

	double t = 0;


// Initiate the rates, To start all the rates are zero except the initiation of the first available site
	int NR = 20;
	double Rates[NR]; // Rate matrix
	Rates[0] = kb1*(gene1-comp1)*Nd;
	Rates[1] = ku1*comp1;

	Rates[2] = kbd*(Decoy-CDecoy)*Nd;
	Rates[3] = kud*CDecoy;

	Rates[4] = gm1*(gene1-comp1); // mRNA production

	Rates[5] = gp1*Nm1; // protein production

	Rates[6] = rm1*Nm1; // mRNA degradation

	Rates[7] = rp1*Nd; // Dimer degradation
	
	Rates[8] = rp1*comp1; // Degradation of TFs bound to gene

	Rates[9] = rp1*CDecoy; // Degradation of TFs bound to decoy

	Rates[10] = 0.5*rd*Np1*(Np1-1);; // Dimerization
	Rates[11] = Nd*rm; // Monomerization
	Rates[12] = Np1*rp1; // Monomer degradation

	// Gene 2
	Rates[13] = kb2*(gene2-comp2)*Nd;
	Rates[14] = ku2*comp2;
	Rates[15] = gm2*(gene2-comp2);
	Rates[16] = gp2*Nm2; // protein production
	Rates[17] = rm1*Nm2; // mRNA degradation
	Rates[18] = Np2*rp2; // Monomer degradation
	Rates[19] = rp1*comp2; // Degradation of TFs bound to gene
		
//Start Gillespie

	/*FILE *filet;
	char str[200], buf[4];
	strcpy(str,"TimeSeries_Decoy-");
	sprintf(buf,"%d",Decoy);
	strcat(str,buf);
        strcat(str,".txt");	
	filet = fopen(str,"w");*/

	double *PSum, K, dt, flag;
	int Ind, Reaction;

	int *NOut;
	NOut = ivmalloc(8);

	while (t < tend){

		NOut[0] = Nm1;
		NOut[1] = Np1;
		NOut[2] = comp1;
		NOut[3] = Nm2;
		NOut[4] = Np2;
		NOut[5] = comp2;
		NOut[6] = CDecoy;
		NOut[7] = Nd;

		PSum = cumsum(Rates,NR);
		K = PSum[NR-1];
		dt = exprnd(K);
//		printf("%f %f\n",K,dt);

		t = t+dt; // update time

/********* chose reaction **************/
		flag = unifrnd(0.00000001,K);
		Ind = FindInd(PSum,flag); // index of the reaction to execute
		free(PSum);
		Reaction = Ind; 
		switch (Reaction){
			case 0:
				comp1 = comp1 + 1; //Binding of protein1 to gene1
				Nd = Nd - 1;
				break;
				
			case 1:
				comp1 = comp1 - 1; // unbinding
				Nd = Nd + 1;
				break;
				
			case 2:
				CDecoy = CDecoy + 1; // Binding of protein to Decoy
				Nd = Nd - 1;
				break;

			case 3:
				CDecoy = CDecoy - 1; // Unbinding of protein to Decoy
				Nd = Nd + 1;
				break;
				
			case 4:
				Nm1 = Nm1 + 1; // mrna production
				break;
				
			case 5:
				Np1 = Np1 + 1;
				break;
				
			case 6:
				Nm1 = Nm1 - 1; // mrna degradation
				break;
				
			case 7:
				Nd = Nd - 1;
				break;
				
			case 8:
				comp1 = comp1 - 1;
				break;
				
			case 9:
				CDecoy = CDecoy - 1;
				break;

			case 10:
				Nd = Nd + 1;
				Np1 = Np1-2;
				break;
				
			case 11:
				Nd = Nd - 1;
				Np1 = Np1 + 2;
				break;

			case 12:
				Np1 = Np1 - 1;
				break;

			case 13:
				comp2 = comp2 + 1;
				Nd = Nd - 1;
				break;

			case 14:
				comp2 = comp2 - 1;
				Nd = Nd + 1;
				break;

			case 15:
				Nm2 = Nm2 + 1;
				break;

			case 16:
				Np2 = Np2 + 1;
				break;

			case 17:
				Nm2 = Nm2 - 1;
				break;

			case 18:
				Np2 = Np2 - 1;
				break;

			default:
				comp2 = comp2 - 1;

		}

/*************** Recalculate rates ************************************/
		Rates[0] = kb1*(gene1-comp1)*Nd;
		Rates[1] = ku1*comp1;
		Rates[2] = kbd*(Decoy-CDecoy)*Nd;
		Rates[3] = kud*CDecoy;
		Rates[4] = gm1*(gene1-comp1);
		Rates[5] = gp1*Nm1;
		Rates[6] = rm1*Nm1;
		Rates[7] = rp1*Nd; 
		Rates[8] = rp1*comp1;
		Rates[9] = rp1*CDecoy;  
		Rates[10] = 0.5*rd*Np1*(Np1-1);
		Rates[11] = Nd*rm;
		Rates[12] = Np1*rp1;
		Rates[13] = kb2*(gene2-comp2)*Nd;
		Rates[14] = ku2*comp2;
		Rates[15] = gm2*(gene2-comp2);
		Rates[16] = gp2*Nm2;
		Rates[17] = rm1*Nm2;
		Rates[18] = Np2*rp2;
		Rates[19] = rp1*comp2;

		//fprintf(filet,"%f %d %d\n",t,Np1+2*(Nd+comp1+comp2+CDecoy),Np2);
	}
	//fclose(filet);
	 
	return NOut;
}

/********* section to run code for a given input and compute fold change ***********/
double* Run_code(int Decoy, double ku1, double ku2, double kud, double kb, double rp, double rp2, double gp1, double gm1, double rm);
double* Run_code(int Decoy, double ku1, double ku2, double kud, double kb, double rp, double rp2, double gp1, double gm1, double rm)
{
	srand(time(NULL));
	time_t tstart, tfinish;
	tstart = time (NULL);
	FILE *file;
	char str[200], buf4[6], buf1[6], buf2[6], buf3[6], buf5[6], buf6[6], buf7[6];

	//int Decoy = 300;

	sprintf(buf1,"%1.5f",ku1);
	sprintf(buf2,"%1.5f",ku2);
	sprintf(buf3,"%1.5f",kud);
	sprintf(buf4,"%1.5f",kb);
	sprintf(buf5,"%1.5f",rp);
	sprintf(buf6,"%1.2f",gp1);
	sprintf(buf7,"%d",Decoy);



	/*strcpy(str,"DATA//FINAL/Oct_2019/RAW/Ku1-");
	strcat(str,buf1);
	strcat(str,"_Ku2-");
	strcat(str,buf2);
	strcat(str,"_Kud-");
	strcat(str,buf3);
	strcat(str,"_Kb-");
	strcat(str,buf4);
	strcat(str,"_rp-");
	strcat(str,buf5);
	strcat(str,"_gp-");
	strcat(str,buf6);
	strcat(str,"_Decoy-");
	strcat(str,buf7);
        strcat(str,".txt");
	file = fopen(str,"w");*/

	double gp2 = gp1, gm2 = gm1;

	double Xst = gp1*gm1/(rp*rm); // Constitutive epxression
	double Yst = gp2*gm2/(rp2*rm);

//	double ku1 = 0.0024, ku2 = 0.0024,kud = 0.0024;

	int Nit = 5*pow(10,3); // Number of iteration for averaging
	int *NOut, *Ntemp;
	int Nin[8];
	double tend = 1*pow(10,5); // sampling time for averaging 
	double sum = 0, sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
	
	Ntemp = Gillespie(Decoy,pow(10,6), 0, 0, 0, 0, 0, 0, 0, 0, ku1, ku2, kud,kb, rp, rp2,gm1,gm2,gp1,gp2,rm,Xst,Yst);
	Nin[0] = Ntemp[0], Nin[1] = Ntemp[1], Nin[2] = Ntemp[2], Nin[3] = Ntemp[3];
	Nin[4] = Ntemp[4], Nin[5] = Ntemp[5], Nin[6] = Ntemp[6], Nin[7] = Ntemp[7];
	free(Ntemp);
	//printf("%d %d %d %d\n",Nin[0],Nin[1],Nin[2],Nin[3]);
	
	for (int ni=0;ni<Nit;ni++){
		NOut = Gillespie(Decoy,tend,Nin[0],Nin[1],Nin[2],Nin[3],Nin[4],Nin[5],Nin[6],Nin[7],ku1,ku2,kud,kb,rp,rp2,gm1,gm2,gp1,gp2,rm,Xst,Yst);
		Nin[0] = NOut[0], Nin[1] = NOut[1], Nin[2] = NOut[2], Nin[3] = NOut[3];
		Nin[4] = NOut[4], Nin[5] = NOut[5], Nin[6] = NOut[6], Nin[7] = NOut[7];
		sum = sum + Nin[1]+2*(Nin[7]+Nin[2]+Nin[5]+Nin[6]);
		sum1 = sum1 + Nin[4];
		sum2 = sum2 + Nin[7];
		sum3 = sum3 + Nin[2];
		sum4 = sum4 + Nin[5];
		//fprintf(file,"%d %f %f\n",Decoy,(Nin[1]+2*(Nin[7]+Nin[2]+Nin[5]+Nin[6]))/Xst,Nin[4]/Yst);
		free(NOut);
	}
	//fclose(file);
	tfinish = time (NULL);
	long comp_time = (tfinish - tstart);
	printf("Computation time = %ld %f\n", comp_time, Xst);
	//printf("%f\n",sum2/Nit);
	double *Exp;
	Exp = dvmalloc(5);
	Exp[0] = sum/(Xst*Nit); // Fold change TF 
	Exp[1] = sum1/(Yst*Nit); // Fold change target
	Exp[2] = sum2/Nit; // Average TF dimers
	Exp[3] = sum3/Nit;	// average TF gene occupancy
	Exp[4] = sum4/Nit;	// average target gene occupancy
	return Exp;
}

/* Main code to for providing inputs and storing final results */
int main()
{	
	time_t tstart, tfinish;
	tstart = time (NULL);

	FILE *file;

	int Decoy,i,n;
	double *Exp;

	int T12 = 25; // Cell Division time in minutes
	double rp = log(2)/(T12*60), rp2 = log(2)/(T12*60); // rp==TF, rp2 = Target
	double EO1 = 15.3, EO2 = 13., EOid = 17; // Binding energy of TF binding sites in kT units

	double kb = 0.0015;
	int XX = 2600; // Constitutive protein expression
	double gm = 0.1;
	double rm = 1./30;
	double gp = XX*rp*rm/gm;
//	double gm = XX*rp*rm/gp;
//	double ku1 =  5*pow(10,6)*kb/exp(EO1)-rp, ku2 = 5*pow(10,6)*kb/exp(EO1)-rp, kud = 5*pow(10,6)*kb/exp(EOid)-rp;
	double ku1 =  0.00149, ku2 = ku1, kud = 0.000415;
//	printf("%f %f %f\n",ku1,ku2,kud);
	char str[200], buf1[6], buf2[6], buf3[6], buf4[6], buf5[6], buf6[6], buf7[6],buf8[6],buf9[6];
	sprintf(buf1,"%1.5f",ku1);
	sprintf(buf2,"%1.5f",ku2);
	sprintf(buf3,"%1.5f",kud);
	sprintf(buf4,"%1.5f",kb);
	sprintf(buf5,"%d",T12);
	sprintf(buf6,"%1.2f",gp);
	sprintf(buf7,"%1.3f",gm);
	sprintf(buf8,"%1.3f",rm);
	sprintf(buf9,"%d",XX);

	strcpy(str,"Asymmetry_Ku-");
	strcat(str,buf1);
	strcat(str,"_Kud-");
	strcat(str,buf3);
	strcat(str,"_kb-");
	strcat(str,buf4);
	strcat(str,"_rp-");
	strcat(str,buf5);
	strcat(str,"min_XX-");
	strcat(str,buf9);
	strcat(str,"_gm-");
	strcat(str,buf7);
	strcat(str,"_rm-");
	strcat(str,buf8);
	strcat(str,"-3.txt");
	
	file = fopen(str,"w");

	for (i=0;i<1;i++){
		Decoy = 10*i; // vary number of decoy sites
		Exp = Run_code(Decoy,ku1,ku2,kud,kb,rp,rp2,gp,gm,rm);
		fprintf(file,"%d %f %f %f\n",Decoy,Exp[0],Exp[1],Exp[2]); // print fold change for TF, target and average TF dimers
		printf("%f %f %d %f %f %f\n",ku1,rp,Decoy,Exp[0],Exp[1],Exp[2]);
		free(Exp);
	}
	fclose(file);
	tfinish = time (NULL);
	long comp_time = (tfinish - tstart);
	printf("Computation time in min = %ld\n", comp_time);
}
