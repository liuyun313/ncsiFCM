#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<cstring>
#include<fstream>
#include<memory.h>
#include<math.h>
#include "FCM.h"
//#include"class.h"
#include<time.h>
//#include"FCM.cpp"



using namespace ker;
using namespace std;

//******************************end****************************//

void kmer(char *s,int *p,int bp)//kmer calculation****n reads with l bp
{
    int k=4,l=bp;
    int *g = new int[l];
    memset(g,0,sizeof(int)*l);
    int fVector[256]={0};
	//********convert the "ATGC" to "0123"
	for(int i=0;i<l;i++)
	{
		if(*s=='A')g[i]=0;
		if(*s=='T')g[i]=1;
		if(*s=='G')g[i]=2;
		if(*s=='C')g[i]=3;
		s++;
	}
	//*****kmer calculate*****
	for(int i=0;i<l-3;i++)
	{
		fVector[64*g[i]+16*g[i+1]+4*g[i+2]+g[i+3]]++;		
	}
	delete g;
	//return fVector;
	for(int i=0;i<256;i++)
	{	
		*(p+i)=fVector[i];
	}
}





void normalization(int **X1,double **X,int num)
{
	double max,min;
	for(int i=0;i<256;i++)
	{
		max=X1[0][i];
		for(int j=1;j<num;j++)
		{
			if(X1[j][i]>max)max=X1[j][i];
		}
		for(int j=0;j<num;j++)
		{
			if(max!=0)X[j][i]=X1[j][i]/max;
		}
	}
}



int main(int argc, char *argv[])
{
	time_t t1,t2;
	t1=time(NULL);
	char f[100];
	strcpy(f, argv[1]);
	int k_=atoi(argv[2]);
//***********************file read*******************//
	//calculate the number of symbols in the input file
	ifstream in(f);//("1a.fna");
	in.seekg(0, ios::end);      //set the file point to the end of the file stream
	streampos ps = in.tellg();  //read the position the poing,while ps will be the symbol number
	in.close();                 //关闭文件流

	FILE *fd;
	int len=int(ps);
	char *str=new char[len];
	fd=fopen(f,"r");//("1a.fna","r");
	fread(str,1,len,fd);//read the file into str
	fclose(fd);
//******************************end*******************//

//************separate the title and the DNA read according to '\n'********//
	int num=0;
	int length=2e4;//number of reads contains in the dataset
	for(int i=0;i<len;i++)
	{
		if(str[i]=='>')num++;
	}
	char title[num][length],read[num][length];
	int bp[num];//length of each read
	num=0;
	int i=0; 
	while(i<len)
	{
		if(str[i]=='>')
		{
			int j=0;
			while(!(str[i]=='\n'))
			{
				title[num][j]=str[i];
				i++;
				j++;
			}
			i++;
			j=0;
			while(!(str[i]=='>')&&i<len)
			{
			   if(!(str[i]=='\n'))
			   {
			        read[num][j]=str[i];
				j++;
			   }	
				i++;
			}
			bp[num]=j-1;
			i--;
			num++;
		}
		i++;
	} 

	delete str;
//*********************end**********************//
       
	
//***************calculate the k-mer frequency*************//
	//int x[256]={0},*p;
	int **X1=new int*[num];//feature vector with dimension num*256
	for(i=0;i<num;i++)
	{
		X1[i]=new int[256];
		memset(X1[i],0,sizeof(int)*256);
	}
	//p=x;
	cout<<"Counting kmer frequencies......"<<endl;
	for(i=0;i<num;i++)
	{
		kmer(read[i],X1[i],bp[i]);
		
	}
//********************end******************//


//**************feature vector normalization******//
	

       double **X = new double*[num];
	for (int i = 0; i < num; i++)
	{
		X[i] = new double[256];
	}

	normalization(X1,X,num);
	
//////////////////////////////////////////////////////
	for(i=0;i<num;i++)delete X1[i];
	delete[] X1;
//for(i=0;i<256;i++)cout<<X[0][i]<<' ';
//******************end*********************//



//******c-means clustering for the initials of model-based method****//

	int *label;//[n][num];
	label=new int[num];
	memset(label,0,sizeof(int)*num);

	//cmeans(p_X,p_label,num,kmax,kmin);
//************************end**************//

//*****************FCM clustering method********************//
	
	cout<<"ncsiFCM based binning......"<<endl;
	FCM *fcm=new FCM(num,256,k_);
	fcm->Cluster(X,label);
	
	delete fcm;


//*******************end********************//	


//output the binning result
	cout<<"Writing binning result into file......"<<endl;
	ofstream out("result.txt");
	if(out.is_open())
	{
		for(int j=0;j<k_;j++)
		{
			out<<"Clusters "<<j<<endl;
			for(i=0;i<num;i++)
			{
				if(label[i]==j)out<<title[i]<<endl;								
			}
		}
	}

	t2=time(NULL);
	cout<<"The running time is "<<difftime(t2,t1)<<'s'<<endl;



//****************end*****************//
}

