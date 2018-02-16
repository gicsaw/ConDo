/*
    Program name: ConDo(Contact based Domain boundary prediction)
    File name: feature.c 
    This program generate input features for machine learning
    This program was developed by Seung Hwan Hong 
    (shhong@kias.re.kr)
    Ref: Protein Domain boundary prediction using Co-evolutionary information
    icc feature.c -o feature -qopenmp -O2     or 
    gcc feature.c -o feature -lm -fopenmp -O2
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

int maxline=2000 ;
char seqcode1[]="-ARNDCQEGHILKMFPSTWYVX" ;
int Ntype=21;

char *cstring(int n1)
{
	char *arr;
	arr=(char *) malloc(n1*sizeof(char));
	return arr;
}

int *iarray1(int n1)
{
	int *arr;
	arr=(int *) malloc(n1*sizeof(int));
	return arr;
}

double *darray1(int n1)
{
	double *arr;
	arr=(double *) malloc(n1*sizeof(double));
	return arr;
}

int **iarray2(int n1, int n2)
{
	int i;
	int **arr;

	arr=(int **) malloc(n1 *sizeof(int*));
	arr[0]=(int *) malloc(n1*n2*sizeof(int));
	for(i=1;i<n1;i++)
    {
        arr[i]=arr[i-1]+n2;
    }
	return arr;
}
double **darray2(int n1, int n2)
{
	int i;
	double **arr;

	arr=(double **) malloc(n1 *sizeof(double*));
	arr[0]=(double *) malloc(n1*n2*sizeof(double));
	for(i=1;i<n1;i++)
    {
        arr[i]=arr[i-1]+n2;
    }
	return arr;
}

int ****iarray4(int n1, int n2, int n3, int n4)
{
    int i,j,k ;
    int ****arr;

    arr=(int ****) malloc(n1*sizeof(int***));

    arr[0]=(int ***) malloc(n1*n2*sizeof(int**));
    arr[0][0]=(int **) malloc(n1*n2*n3*sizeof(int*));
    arr[0][0][0]=(int *) malloc(n1*n2*n3*n4*sizeof(int));

    for(i=0;i<n1;i++)
    {
        if(i>0)
        {
            arr[i]=arr[i-1] + n2 ;
            arr[i][0]=arr[i-1][0]+n2*n3;
            arr[i][0][0]=arr[i-1][0][0] + n2*n3*n4 ;
        }
        for(j=0;j<n2;j++)
        {
            if(j>0)
            {
                arr[i][j]=arr[i][j-1] + n3 ;
                arr[i][j][0]=arr[i][j-1][0] + n3*n4 ;
            }
            for(k=0;k<=n3;k++)
            {
                if(k>0)
                {
                    arr[i][j][k]=arr[i][j][k-1] + n4 ;
                }
            }
        }
    }
    return arr;
}



void free_cstring(char *arr)
{
    free(arr);
}

void free_iarray1(int *arr)
{
    free(arr);
}
void free_darray1(double *arr)
{
    free(arr);
}

void free_iarray2(int **arr)
{
    free(arr[0]);
    free(arr);
}

void free_darray2(double **arr)
{
    free(arr[0]);
    free(arr);
}

void free_iarray4(int ****arr)
{
    free(arr[0][0][0]);
    free(arr[0][0]);
    free(arr[0]);
    free(arr);
}



void print_time(time_t time_0)
{
    time_t time_now, running_time ;
    time(&time_now);
    printf("%s",ctime(&time_now));
    running_time = time_now - time_0   ;   
    printf("h : %ld m : %ld s : %ld \n" , running_time/3600, (running_time/60)%60, running_time%60);}

char *read_seq(char *target, int *Nres )
{
    FILE *fp;
    char filename[100];
    char line[maxline];
    char line2[maxline];
    char *seq;
    int i,j,Nres_temp;
    int Nchar;
    char *eof;
    strcpy(line,"");
    strcpy(line2,"");

    sprintf(filename,"%s.fasta",target);
    fp=fopen(filename, "r");
    while(fgets(line,maxline,fp)!=NULL)
    {
        if(line[0]!='>')
        {
            Nchar=strlen(line)-1;
            strncat(line2,line,Nchar);
        }
    }
    fclose(fp) ;

    Nres_temp=strlen(line2) ;
    *Nres=Nres_temp;

    seq=cstring(Nres_temp+1);
    strncpy(seq,line2,Nres_temp);
    seq[Nres_temp]='\0';

    return seq ;
}


int *mod_seq(char *seq, int Nres)
{
    int i,j;
    int *seqn;

    seqn=iarray1(Nres);
    for(i=0;i<Nres;i++)
    {
        for(j=0;j<=Ntype;j++)
        {
            if(seq[i]==seqcode1[j] || j==Ntype)
            {
                seqn[i]=j;
                break ;
            }
        }
    }

    return seqn ;
}

double **read_ss2(char *target, int Nres )
{
    FILE *fp;
    char filename[100];
    char line[maxline];
    int i,j;
    char *eof;
    char amino,ss;
    double **ss2;

    ss2=darray2(Nres,3);

    sprintf(filename,"%s.ss2",target);
    fp=fopen(filename, "r");
    i=-1;
    while(fgets(line,maxline,fp)!=NULL)
    {
        if(strlen(line)<=10)
        {
            continue;
        }
        if(line[0]=='#')
        {
            continue;
        }
        i+=1;
        sscanf(line,"%d %c %c %lf %lf %lf",&j,&amino,&ss,&ss2[i][0],&ss2[i][1],&ss2[i][2]);

    }

    fclose(fp) ;

    return ss2 ;
}

double **read_sa2(char *target, int Nres )
{
    FILE *fp;
    char filename[100];
    char line[maxline];
    int i,j;
    char *eof;
    char amino,sa;
    double rsa, tsa, sab,sai,sae,con;
    double **sa2;

    sa2=darray2(Nres,3);
    sprintf(filename,"%s.a22",target);
    fp=fopen(filename, "r");

    i=-1;
    while(fgets(line,maxline,fp)!=NULL)
    {
        if(strlen(line)<=10)
        {
            continue;
        }
        if(line[0]=='#')
        {
            continue;
        }
        i+=1;
        sscanf(line,"%d %c %c %lf %lf %lf",&j,&amino,&sa,&sa2[i][1],&sa2[i][2],&con);
        if(con<=20.0)
        {
            sa2[i][3]=con/20.0;
        }
        else
        {
            sa2[i][3]=1.0;
        }

    }

    fclose(fp) ;

    sprintf(filename,"%s.a3",target);
    fp=fopen(filename, "r");

    i=-1;
    while(fgets(line,maxline,fp)!=NULL)
    {
        if(strlen(line)<=10)
        {
            continue;
        }
        if(line[0]=='#')
        {
            continue;
        }
        i+=1;
        sscanf(line,"%d %c %c %lf %lf %lf %lf %lf %lf",&j,&amino,&sa,&sab,&sai,&sae,&rsa,&tsa,&con);
        if(rsa<=1.0)
        {
            sa2[i][0]=rsa;
        }
        else
        {
            sa2[i][0]=1.0;
        }
    }

    fclose(fp) ;

    return sa2 ;
}

void find_inifin(char line[], int *ini, int *fin)
{
    char cini[10], cfin[10];
    char chr='-';
    int i, k, len ;
    len=strlen(line);
    for(i=0;i<len;i++)
    {
        if (line[i]==chr)
        {
            k=i;
            break;
        }
    }

    for(i=0;i<k;i++)
    {
        cini[i]=line[i];
    }
    cini[k]='\0';
    *ini=atoi(cini);
    for(i=0;i<len-k-1;i++)
    {
        cfin[i]=line[i+k+1];
    }
    cfin[len-k-1]='\0';
    len=strlen(cfin);

    *fin=atoi(cfin);

}

int find_Nmsa(char *target)
{
    FILE *fp;
    char filename[100];
    char line[maxline];
    int k;
    sprintf(filename,"%s.msa",target);
    fp=fopen(filename, "r");
    k=0;
    while(fgets(line,maxline,fp)!=NULL)
    {
        if(line[0]=='>')
        {
            k+=1;
        }
    }
    fclose(fp);
    return k;
}

void read_msa(char *target, int **msa, int **msa2, int *pasinfo, int Nres, int Nmsa)
{
    FILE *fp;
    char filename[100];
    char line[maxline];
    int i,j,k, l;
    int ini, fin;
    char *ptr;
    char delim[]="/,";
    char inifin[10];
    
    sprintf(filename,"%s.msa",target);
    fp=fopen(filename, "r");
    k=-1;
    while(fgets(line,maxline,fp)!=NULL)
    {
        if(line[0]=='>')
        {
            k+=1;
            for(i=0;i<Nres;i++)
            {
                msa2[k][i]=0;
            }
            ptr=strtok(line,delim);
            while(ptr=strtok(NULL,delim))
            {
                strcpy(inifin,ptr);

                find_inifin(inifin,&ini,&fin); 
                for(i=ini-1;i<=fin-1;i++)
                {
                    msa2[k][i]=1;
                }
            }
            if((ini>20)||(fin<=Nres-20))
            {
                pasinfo[k]=1;
            }
            else
            {
                pasinfo[k]=0;
            }
            continue;
        }
        for(i=0;i<Nres;i++)
        {
            for(j=0;j<=Ntype;j++)
            {
                if(line[i]==seqcode1[j] || j==Ntype)
                {
                    msa[k][i]=j;
                    break ;
                }
            }
        }

    }
    fclose(fp);

} 


double **gen_profile(int **msa, int Nres, int Nmsa)
{
    int i,j,k;

    double **profile;
    int **frequency;
    int type;


    profile=darray2(Nres,Ntype);
    frequency=iarray2(Nres,Ntype);

    #pragma omp parallel for private(i,j)
    for(i=0;i<Nres;i++)
    {
        for(j=0;j<Ntype;j++)
        {
            frequency[i][j]=0;
        }
    }

    #pragma omp parallel for private(i,k,type)
    for(i=0;i<Nres;i++)
    { 
        for(k=0;k<Nmsa;k++)
        {
            type=msa[k][i];
            if(type<Ntype)
            {
                frequency[i][type]+=1;
            }
        }
    }

    #pragma omp parallel for private(i,j)
    for(i=0;i<Nres;i++)
    {
        for(j=0;j<Ntype;j++)
        {
            profile[i][j]=(double) frequency[i][j] /(double) Nmsa;
        }
    }

    free_iarray2(frequency);
    return profile;
}

double **read_profile(char *target, int Nres)
{

    FILE *fp;
    char filename[100];
    char line[maxline];
    char *ptr;

    double **profile;
    int i,j,k, m;
    char delim[]=" ";
    int Nres2;
    char seq2[maxline];

    profile=darray2(Nres,Ntype);
   
    sprintf(filename,"%s.ck2",target);
    fp=fopen(filename, "r");

    fgets(line,maxline,fp);
    Nres2=atoi(line);
    fgets(line,maxline,fp);
    strncat(seq2,line,Nres2);

    k=-1;

    while(fgets(line,maxline,fp)!=NULL)
    {
        k+=1;
        profile[k][0]=0.0;
        m=0;

        ptr=strtok(line,delim);
        m+=1;
        profile[k][m]=atof(ptr);
//        printf("%s ",ptr);
        while(ptr=strtok(NULL,delim))
        {
            m+=1;
            profile[k][m]=atof(ptr);
//            printf("%s ",ptr);
        }
//        printf("\n");
    }
//    for(i=0;i<Nres;i++)
//    {
//        for(m=0;m<Ntype;m++)
//        {
//            printf("%6.4f ", profile[i][m]);
//        }
//        printf("\n");
//    }

    return profile;
}



void get_n_signal(int **msa2, double **msig5, double **psig5, int **Nsite, int Nres, int Nmsa)
{

    int i,j,k,l;
    int ini, fin;
    int **Nsig, **Nsig5;
    int maxsig5[3];
    double sum[3], sum2[3];
    double ave[3],ave2[3],sigma[3];
    int ggap;

    Nsig=iarray2(Nres,3);
    Nsig5=iarray2(Nres,3);

    #pragma omp parallel for private(i,j)
    for(i=0;i<Nres;i++)
    {
        for(j=0;j<=2;j++)
        {
            Nsig[i][j]=0;
            Nsig5[i][j]=0;
            Nsite[i][j]=0;
            msig5[i][j]=0.0;
            psig5[i][j]=0.0;
        }
    }

//Nsig[i][0]: left_gap &right_gap
//Nsig[i][1]: left_gap 
//Nsig[i][2]: right_gap

    #pragma omp parallel for private(i,k)
    for(i=1;i<Nres-1;i++)
    {
        for(k=0;k<Nmsa;k++)
        {
            if((msa2[k][i-1]==0)&&(msa2[k][i]==1))
            {
                Nsig[i][0]+=1;
                Nsig[i][1]+=1;
            }
            else if((msa2[k][i]==1)&&(msa2[k][i+1]==0))
            {
                Nsig[i][0]+=1;
                Nsig[i][2]+=1;
            }
        }
    }

    #pragma omp parallel for private(i,j,ini,fin,l)
    for(i=0;i<Nres;i++)
    {
        ini=-5;
        fin=5;
        if(i<5)
        {
            ini=-i;
            Nsite[i][0]+=(5-i);
            Nsite[i][1]+=(5-i);
        }

        if(i>Nres-6)
        {
            fin=Nres-i-1;
            Nsite[i][0]+=(i-(Nres-6));
            Nsite[i][2]+=(i-(Nres-6));

        }
//            printf("i:%d ini:%d fin:%d \n",i,ini,fin);


        for(l=ini;l<=fin;l++)
        {
            for(j=0;j<=2;j++)
            {
                Nsig5[i][j]+=Nsig[i+l][j];
                if(Nsig[i+l][j]>0)
                {

                    Nsite[i][j]+=1;
                }

            }
        }
    }

    #pragma omp parallel for private(i)
    for(i=0;i<Nres;i++)
    {
        Nsite[i][0]= Nsite[i][1]+ Nsite[i][2];
//        printf("%d %d %d %d",i,Nsig[i][0],Nsig[i][1],Nsig[i][2]);
//        printf("%d %d %d %d\n",i,Nsite[i][0],Nsite[i][1],Nsite[i][2]);
    }

    free_iarray2(Nsig);


    ggap=20*Nres/100;
    if(ggap<60)
    {
        ggap=60;
    }
//    ggap=0;
    ini=ggap;
    fin=Nres-ggap;
//    printf("ggap %d \n",ggap);
    for(j=0;j<=2;j++)
    {
        maxsig5[j]=1;
        for(i=ini;i<fin;i++)
        {
            if(maxsig5[j]<Nsig5[i][j])
            {
                maxsig5[j]=Nsig5[i][j];
            }
        }
//        printf("maxsig= %d\n", maxsig5[j]);
        for(i=0;i<Nres;i++)
        {
            if(maxsig5[j]<Nsig5[i][j])
            {
                Nsig5[i][j]=maxsig5[j];
            }
        }
    }

    for(j=0;j<=2;j++)
    {
        sum[j]=0.0;
        sum2[j]=0.0;
        for(i=0;i<Nres;i++)
        {
            if(j>0)
            {
                psig5[i][j]=((double) Nsig5[i][j]/ ((double) Nmsa)) ;
                msig5[i][j]=((double) Nsig5[i][j]/ ((double) maxsig5[j])) ;
            }  
            else
            {
                psig5[i][j]=((double) Nsig5[i][j]/ ((double) Nmsa * 2.0)) ;
                msig5[i][j]=((double) Nsig5[i][j]/ ((double) maxsig5[j] * 2.0)) ;
            }
        }
    }

    free_iarray2(Nsig5);

    return ;
}

void make_PAS(double **PAS, double **PAS2, double** PAS3, int **msa2, int *pasinfo, int Nres, int Nmsa,int Npasinfo)
{
    int ****fc ;
    int **fs ;

    int i,j,k, m,n;

    fc=iarray4(Nres,Ntype+1,Nres,Ntype+1);
    fs=iarray2(Nres,Ntype+1) ;

    #pragma omp parallel for private(i,j,m,n)
    for(i=0;i<Nres;i++)
    {
        for(j=0;j<Nres;j++)
        {
            PAS[i][j]=0.0;
            PAS2[i][j]=0.0;
            PAS3[i][j]=0.0;
        }
        for(m=0;m<Ntype;m++)
        {
            fs[i][m]=0;
            for(j=0;j<Nres;j++)
            {
                for(n=0;n<Ntype;n++)
                {
                    fc[i][m][j][n]=0;
                }
            }
        }
    }


    for(k=0;k<Nmsa;k++)
    {
        if(pasinfo[k]==0)
        {
            continue ;
        }
        #pragma omp parallel for private(i,j,m,n)
        for(i=0;i<Nres;i++)
        {
            m=msa2[k][i] ;
            fs[i][m]+=1 ;
            for(j=i;j<Nres;j++)
            {
                n=msa2[k][j];
                fc[i][m][j][n]+=1;
            }
        }
    }

    #pragma omp parallel for private(i,m,j,n)
    for(i=0;i<Nres;i++)
    {
        for(j=i;j<Nres;j++)
        {
            m=0;
            n=0;
            PAS[i][j] += (double) fc[i][m][j][n]/(double) Npasinfo;
            PAS2[i][j] += (double) fc[i][m][j][n]/(double) Npasinfo;
//            PAS[i][j] += (double) fc[i][m][j][n]/(double) Nmsa;
//            PAS2[i][j] += (double) fc[i][m][j][n]/(double) Nmsa;
            m=1;
            n=1;
            PAS[i][j] += (double) fc[i][m][j][n]/(double) Npasinfo;
            PAS3[i][j] += (double) fc[i][m][j][n]/(double) Npasinfo;
//            PAS[i][j] += (double) fc[i][m][j][n]/(double) Nmsa;
//            PAS3[i][j] += (double) fc[i][m][j][n]/(double) Nmsa;    
        }
    }

    #pragma omp parallel for private(i,j)
    for(i=0;i<Nres;i++)
    {
        for(j=i+1;j<Nres;j++)
        {
            PAS[j][i]=PAS[i][j];
            PAS2[j][i]=PAS2[i][j];
            PAS3[j][i]=PAS3[i][j];    
        }
    }
    free_iarray4(fc);
    free_iarray2(fs);
}


void sum_PAS(double **PASaveN, double **PASaveC, double **PAS3, int Nres, int win)
{
    int i,j,k, ii ;
    int ini, fin;

    double PASsumC, PASsumN;

    #pragma omp parallel for private(i,j,k,ii,ini,fin,PASsumC,PASsumN)
    for(i=0;i<Nres;i++)
    {
        for(k=-win;k<=win;k++)
        {
            PASaveN[i][k+win]=0.0;
            PASaveC[i][k+win]=0.0;
        }

        ini=-win;
        fin=+win;
        if(i+ini<0)
        {
            ini=-i;
        }
        if(i+fin>Nres-1)
        {
            fin=Nres-1-i;
        }
        for(k=ini;k<=fin;k++)
        {
            ii=i+k;
            PASsumC=0.0;
            PASsumN=0.0;

            for(j=i;j<Nres;j++)
            {
                PASsumC+=PAS3[ii][j];
            }
            PASaveC[i][k+win]=PASsumC/((double) (Nres-i));
            for(j=0;j<=i;j++)
            {
                PASsumN+=PAS3[ii][j];
            }
            PASaveN[i][k+win]=PASsumN/((double) (i+1));

        }
        for(k=-win;k<ini;k++)
        {
            PASaveN[i][k+win]=PASaveN[i][ini+win];
            PASaveC[i][k+win]=PASaveC[i][ini+win];
        }
        for(k=fin+1;k<=win;k++)
        {
            PASaveN[i][k+win]=PASaveN[i][fin+win];
            PASaveC[i][k+win]=PASaveC[i][fin+win];
        }
    }
}

void interpolation(double *PASaveN2, double *PASaveC2, double *PASaveN5, double *PASaveC5, int Nwin, int Nres)
{
    int j,j1,j2;
    double step, step2;
    double x, x1, x2 ;
    step=1.0/( (double) Nwin);
    step2=1.0/( (double) Nres);

    for(j=0;j<Nwin;j++) 
    {
        x=(j+1)*step;
        j1=(int) (x*Nres)-1;
        j2=j1+1;
        x1=(j1+1)*step2;
        x2=(j2+1)*step2;
        if(j1<0)
        {
            PASaveN2[j]=PASaveN5[j2] ;
            PASaveC2[j]=PASaveC5[j2] ;
        }
        else
        {
            PASaveN2[j]=(x-x1)*(PASaveN5[j2]-PASaveN5[j1])/step2+PASaveN5[j1];
            PASaveC2[j]=(x-x1)*(PASaveC5[j2]-PASaveC5[j1])/step2+PASaveC5[j1];
        }
    }

}

void sum_PAS2(double **PASaveN2, double **PASaveC2, double *PASmaxN, double *PASmaxC, double **PAS3, int Nres, int Nwin)
{
    int i,j,k, ii ;
    int ini, fin;

    double *PASsumN, *PASsumC;
    double *PASaveN5, *PASaveC5;
    double maxN, maxC;


    PASsumN= darray1(Nres);
    PASsumC= darray1(Nres);
    PASaveN5= darray1(Nres);
    PASaveC5= darray1(Nres);

    #pragma omp parallel for private(j,k)
    for(j=0;j<Nres;j++)
    {
        PASsumC[j]=0.0;
        PASsumN[j]=0.0;

        for(k=0;k<=0;k++)
        {
            PASsumN[j]+=PAS3[k][j];
        }

        for(k=0;k<Nres;k++)
        {
            PASsumC[j]+=PAS3[k][j];
        }

        PASaveN5[j]=PASsumN[j]/((double) (1));
        PASaveC5[j]=PASsumC[j]/((double) (Nres));
    }
    maxN=0;
    maxC=0;
    for(j=0;j<Nres;j++)
    {
        if(maxN<PASaveN5[j])
        {
            maxN=PASaveN5[j];
        }
        if(maxC<PASaveC5[j])
        {
            maxC=PASaveC5[j];
        }
    }
    PASmaxN[0]=maxN;
    PASmaxC[0]=maxC;

//    interpolation
    interpolation(PASaveN2[0],PASaveC2[0],PASaveN5,PASaveC5,Nwin,Nres);

    for(i=1;i<Nres;i++)
    {
        #pragma omp parallel for private(j)
        for(j=0;j<Nres;j++)
        {
            PASsumN[j]+=PAS3[i][j];
            PASsumC[j]-=PAS3[i][j];
            PASaveN5[j]=PASsumN[j]/((double) (i+1));
            PASaveC5[j]=PASsumC[j]/((double) (Nres-i));
        }
        maxN=0;
        maxC=0;
        for(j=0;j<Nres;j++)
        {
            if(maxN<PASaveN5[j])
            {
                maxN=PASaveN5[j];
            }
            if(maxC<PASaveC5[j])
            {
                maxC=PASaveC5[j];
            }
        }
        PASmaxN[i]=maxN;
        PASmaxC[i]=maxC;

        interpolation(PASaveN2[i],PASaveC2[i],PASaveN5,PASaveC5,Nwin,Nres);
    }

    free_darray1(PASsumN);
    free_darray1(PASsumC);
    free_darray1(PASaveN5);
    free_darray1(PASaveC5);

}

void write_PAS(double **PAS, char *out_PAS, int Nres)
{
    FILE *fp_PAS;
    int i,j;

    fp_PAS=fopen(out_PAS,"w");

    for(i=0;i<Nres;i++)
    {
        for(j=0;j<Nres;j++)
        {
            fprintf(fp_PAS,"%4d %4d %8.6f\n", i+1,j+1,PAS[i][j]);
        }
    }

    fclose(fp_PAS) ;
}

double **read_ccmpred(char *target, int Nres)
{

    FILE *fp;
    char filename[100];
    char line[maxline];

    double **ccmpred;
    int i,j,k, m;

    ccmpred=darray2(Nres,Nres);
   
    sprintf(filename,"%s.ccmpred",target);
    fp=fopen(filename, "r");

    for(i=0;i<Nres;i++)
    {
        for(j=0;j<Nres;j++)
        {
            fscanf(fp,"%lf",&ccmpred[i][j]);
//            fscanf(fp,"%lf\t",&ccmpred[i][j]);
        }
//        fscanf(fp,"\n");
    }
    fclose(fp);

    for(i=0;i<Nres;i++)
    {
        for(j=i;j<i+12;j++)
        {
            if(j<Nres)
            {
                ccmpred[i][j]=0.0;
                ccmpred[j][i]=0.0;
            }
        }
    }


    return ccmpred;
}

double cal_ccm_cut(double **ccmpred, int Nres)
{

    int i, j ;
    double sum, sum2, Nedge, ave,ave2, dd, stddev;
    double ccm_cut;

    sum=0;
    sum2=0;
    Nedge=0;
    for(i=0;i<Nres;i++)
    {
        for(j=i+12;j<Nres;j++)
        {
            Nedge+=1 ;
            sum+=ccmpred[i][j] ;
            sum2+=pow(ccmpred[i][j],2) ;
        } 
    }

    ave=sum/Nedge;
    ave2=sum2/Nedge;
    dd=ave2-pow(ave,2);
    stddev=sqrt(dd);

    ccm_cut=ave+2*stddev;
//    printf("Nedge: %f ave: %f stddev: %f, ccm_cut: %f \n",Nedge,ave,stddev,ccm_cut);

    return ccm_cut;
}


double *cal_ccm_community(double **ccmpred,double ccm_cut, int Nres ) 
{

    int i, j,k ;
    int Nedge1, Nedge2;
    double sum_con, M;
    double *community;
    double community1, community2;

    double ddd1,ddd2,dm1,dm2,dd1,dd2,dd0;

    community=darray1(Nres);

    sum_con=0;

    for(i=0;i<Nres;i++)
    {
        for(j=i+1;j<Nres;j++)
        {
            if(ccmpred[i][j]>ccm_cut)
            {
               sum_con+=ccmpred[i][j] ;
            }
        } 
    }

    M=sum_con ;

    ddd1=0;
    ddd2=M;
    dm1=0;
    dm2=2*M;


    for(k=0;k<Nres;k++)
    {
        dd1=0;
        dd2=0;
        dd0=0;
        for(i=0;i<k;i++)
        {
            if(ccmpred[k][i]>ccm_cut)
            {
                dd1+=ccmpred[k][i];
                dd0+=ccmpred[k][i];
            }
        }

        for(i=k+1;i<Nres;i++)
        {
            if(ccmpred[k][i]>ccm_cut)
            {
                dd2+=ccmpred[k][i];
                dd0+=ccmpred[k][i];
            }
        }

        ddd1+=dd1;
        ddd2-=dd2;
        dm1+=dd0;
        dm2-=dd0;
        community1=ddd1/M+ddd2/M;
        community2=pow(dm1/(2*M),2)+pow(dm2/(2*M),2);
        community[k]=community1-community2;
    }

    return community;
}

double *interpolation_community(double *ccm_community, int Nwin, int Nres)
{
    int j,j1,j2;
    double step, step2;
    double x, x1, x2 ;

    double *ccm_comm;

    step=1.0/( (double) Nwin);
    step2=1.0/( (double) Nres);

    ccm_comm=darray1(Nwin);

    for(j=0;j<Nwin;j++) 
    {
        x=(j+0.5)*step;
        j1=(int) (x*Nres)-1;
        j2=j1+1;
        x1=(j1+1.5)*step2;
        x2=(j2+1.5)*step2;
        if(j1<0)
        {
            ccm_comm[j]=ccm_community[j2];
        }
        else
        {
            ccm_comm[j]=(x-x1)*(ccm_community[j2]-ccm_community[j1])/step2+ccm_community[j1];
        }
    }

    return ccm_comm;
}

void find_ccm_comm_min_max(double *ccm_community, double *comm_min, double *comm_max, int Nres)
{
    int i;
    double max,min;

    max=-10;
    min=10;
    for(i=0;i<Nres;i++)
    {
        if(max<ccm_community[i])
        {
            max=ccm_community[i];
        }
        if(min>ccm_community[i])
        {
            min=ccm_community[i];
        }
    }
    *comm_max=max;
    *comm_min=min;
}


void write_ccm(double **ccmpred, double ccm_cut, char *out_ccm, int Nres)
{
    FILE *fp_ccm;
    int i,j;

    fp_ccm=fopen(out_ccm,"w");

    for(i=0;i<Nres;i++)
    {
        for(j=0;j<Nres;j++)
        {
            if (ccmpred[i][j]>ccm_cut)
            {
                fprintf(fp_ccm,"%4d %4d %8.6f\n", i+1,j+1,ccmpred[i][j]);
            }
            else
            {
                fprintf(fp_ccm,"%4d %4d %8.6f\n", i+1,j+1,0.0);
            }
        }
    }

    fclose(fp_ccm) ;
}

void write_community(double *ccm_community, char *out_community, int Nres)
{
    FILE *fp_com;
    int i;

    fp_com=fopen(out_community,"w");

    for(i=0;i<Nres;i++)
    {
        fprintf(fp_com,"%4d %8.6f\n", i+1,ccm_community[i]);
    }

    fclose(fp_com) ;
   
}


int main(int argc, char *argv[])
{
    FILE *fp_feature ;
    FILE *fp_feature2 ;

    FILE *fp_PAS, *fp_PAS2, *fp_PAS3 ;

    char out_feature[100];
    char out_feature2[100];

    char *target ;
    char *seq, line[maxline] ;
    int *seqn;
    int **msa, **msa2;
    double **ss2;
    double **sa2;
    int **Nsig;
    int **Nsite;
    double **psig5, **msig5;

    int ****fc ;
    int **fs ;
    int **fz ;
    double **profile;

    int *pasinfo;
    double **PAS, **PAS2, **PAS3;
    double **PASaveN, **PASaveC ;
    double **PASaveN2, **PASaveC2 ;
    double *PASmaxN, *PASmaxC;

    double **ccmpred ;
    double *ccm_community;
    double *ccm_comm;
    double comm_max, comm_min;
    double ccm_cut ;

    int Npasinfo;

    char out_PAS[100], out_PAS2[100], out_PAS3[100] ;
    char out_ccm[100], out_community[100];
   
    double ccc,ccc2 ;
    double Nter, Cter, NCter;
    int Nres, Nmsa ;
    double dmsa, dpas, rpas;

    int win_prof, win_ss, win_sa;
    int win_PAS=50;
    int win_PAS2=100;
    int win_ccm=50;
    int win_ccm2=100;

    int Ncpu ;

    int i,j,k, l, m,n ; 
    int fsim,fsjn;
    time_t time_0, time_now, running_time;

    if(argc==1)
    {
        printf("mi target\n");
        abort();
    }
    target=argv[1];

    printf("target: %s\n",target);

    Ncpu=1 ;
    if(argc>=3)
    {
        Ncpu=atoi(argv[2]);
    }
    omp_set_num_threads(Ncpu);
    printf("Ncpu: %d\n",Ncpu);

    time(&time_0);
    printf("%s",ctime(&time_0));

    seq=read_seq(target,&Nres);
    seqn=mod_seq(seq,Nres);

    printf("Nres: %d\n",Nres) ;

    ss2=read_ss2(target,Nres);
//    printf("read_ss2\n");

//    for(i=0;i<Nres;i++)
//    {
//        printf("%4d %5.3f %5.3f %5.3f\n",i, ss2[i][0],ss2[i][1],ss2[i][2]);
//    }

    sa2=read_sa2(target,Nres);
//    printf("read_sa2\n");

//    for(i=0;i<Nres;i++)
//    {
//        printf("%4d %5.3f %5.3f %5.3f %5.3f\n",i, sa2[i][0],sa2[i][1],sa2[i][2],sa2[i][3]);
//    }

    Nmsa=find_Nmsa(target);
    msa=iarray2(Nmsa+1,Nres+1);
    msa2=iarray2(Nmsa+1,Nres+1);
    pasinfo=iarray1(Nmsa+1);
    read_msa(target,msa,msa2,pasinfo,Nres,Nmsa);

//    for(k=0;k<Nmsa;k++)
//    {
//        printf("%4d ",k);
//        for(i=0;i<Nres;i++)
//        {
//            printf("%d",msa2[k][i]);
//        }
//        printf("\n");
//    }

    Npasinfo=0;
    for(k=0;k<Nmsa;k++)
    {
        if(pasinfo[k]!=0)
        {
            Npasinfo+=1;
        }
    }

    if(Nmsa<10)
    {
        dmsa=0.0;
    }
    else if (Nmsa<10000)
    {
        dmsa= log10(Nmsa)/5.0;
    }
    else
    {
        dmsa=1.0;
    }

    if(Npasinfo<10)
    {
        dpas=0.0;
    }
    else if (Npasinfo<10000)
    {
        dpas=log10(Npasinfo)/5.0;
    }
    else
    {
        dpas=1.0;
    }

    if(Npasinfo==0)
    {
        rpas=0.0;
    }
    else
    {
        rpas=(double) Npasinfo/(double)Nmsa;
    }

    printf("Nmsa: %d dmsa: %f \n", Nmsa,dmsa);
    printf("Npas: %d dpas: %f rpas: %f \n", Npasinfo,dpas,rpas);


//    profile=gen_profile(msa,Nres,Nmsa);
    profile=read_profile(target,Nres);

//    for(i=0;i<Nres;i++)
//    {
//        printf("%4d %5.3f %5.3f %5.3f\n",i, profile[i][0],profile[i][1],profile[i][2]);
//    }

//    printf("initialize\n");

//    print_time(time_0);
    msig5=darray2(Nres,3);
    psig5=darray2(Nres,3);
    Nsite=iarray2(Nres,3);

    get_n_signal(msa2,msig5,psig5,Nsite, Nres,Nmsa);
//    for(i=0;i<Nres;i++)
//    {
//        printf("%d %d %d %d\n",i,Nsite[i][0],Nsite[i][1],Nsite[i][2]);
//    }

//    for(i=0;i<Nres;i++)
//    {
//        printf("%d %f %f %f\n",i,psig5[i][0],psig5[i][1],psig5[i][2]);
//    }

    PAS=darray2(Nres,Nres);
    PAS2=darray2(Nres,Nres);
    PAS3=darray2(Nres,Nres);
    PASaveN=darray2(Nres,2*win_PAS+1);
    PASaveC=darray2(Nres,2*win_PAS+1);
    PASaveN2=darray2(Nres,win_PAS2);
    PASaveC2=darray2(Nres,win_PAS2);
    PASmaxN=darray1(Nres);
    PASmaxC=darray1(Nres);

    if(Npasinfo>0)
    {
        make_PAS(PAS,PAS2,PAS3,msa2,pasinfo, Nres,Nmsa,Npasinfo) ;
        sum_PAS(PASaveN,PASaveC,PAS3,Nres, win_PAS) ;
        sum_PAS2(PASaveN2,PASaveC2,PASmaxN,PASmaxC,PAS3,Nres, win_PAS2) ;
    }

    else
    {
        for(i=0;i<Nres;i++)
        {
            PASmaxN[i]=0;
            PASmaxC[i]=0;

            for(j=0;j<Nres;j++)
            {
                PAS[i][j]=0.0;
                PAS2[i][j]=0.0;
                PAS3[i][j]=0.0;

            }
            for(k=0;k<=2*win_PAS;k++)
            {
                PASaveN[i][k]=0.0;
                PASaveC[i][k]=0.0;
            }
            for(j=0;j<win_PAS2;j++)
            {
                PASaveN2[i][j]=0.0;
                PASaveC2[i][j]=0.0;
            }
        }
    }

//    for(j=0;j<win_PAS2;j++)
//    {
//        printf("%d %f %f\n",j, PASaveN2[1][j], PASaveC2[1][j]);
//    }

//  for(i=0;i<Nres;i++)
//  {
//      printf("%d %f %f\n",i, PASmaxN[i], PASmaxC[i]);
//  }

    ccm_cut=0.0;
    if(Nmsa>5)
    {
        ccmpred=read_ccmpred(target,Nres);
        ccm_cut=cal_ccm_cut(ccmpred,Nres);
        printf("ccm_cutoff: %f\n",ccm_cut);
        ccm_community=cal_ccm_community(ccmpred,ccm_cut, Nres);
    }
    else
    {
        ccmpred=darray2(Nres,Nres);
        ccm_community=darray1(Nres);
        for(i=0;i<Nres;i++)
        {
            ccm_community[i]=0.0;
            for(j=0;j<Nres;j++)
            {
                ccmpred[i][j]=0.0;
            }
        }
    }


//    for(i=0;i<Nres;i++)
//    {       
//        printf("%d %f\n",i, ccmpred[i][Nres-1]);
//    }

    ccm_comm=interpolation_community(ccm_community,win_ccm2, Nres);
    find_ccm_comm_min_max(ccm_community,&comm_min,&comm_max,Nres);

//    printf("max: %f, min: %f \n",comm_max, comm_min);

//    for(i=0;i<Nres;i++)
//    {       
//        printf("%f %f\n",(i+1.5), ccm_community[i]);
//    }
//    for(i=0;i<win_ccm;i++)
//    {
//        printf("%f %f\n",Nres*(i+0.5)/100.0, ccm_comm[i]);
//    }

//    for(j=0;j<win_ccm2;j++)
//    {
//        printf("%d %6.4f\n",j, ccm_comm[j]);
//    }

    sprintf(out_feature,"%s_feature.txt",target);
    sprintf(out_feature2,"%s_feature2.txt",target);
    fp_feature=fopen(out_feature,"w");
    fp_feature2=fopen(out_feature2,"w");

//    fprintf(fp_feature,"#res exist(100) ccm_community[100] comm_max comm_min PAS_N(k,101) PAS_C(k,101) PASmaxN PASmaxC position profile(k,m) ss(41,3) sa2(41,4) dmsa dpas rpas Nter Cter Len Nsite msig5 psig5 ");
//    fprintf(fp_feature2,"#res exist(40) ccm_comm[100] PAS_N(k,100) PAS_C(k,100) position profile(k,m) ss(41,3) sa2(41,4) dmsa dpas rpas Nter Cter Len Nsite msig5 psig5 ");

    win_prof=10;
    win_ss=20;
    win_sa=20;


    for(i=0;i<Nres;i++)
    {

        fprintf(fp_feature,"%4d", i+1);
        fprintf(fp_feature2,"%4d", i+1);
//exist
        for(j=-win_PAS;j<=win_PAS;j++)
        {
            if(j==0)
            {
                continue;
            }
            k=i+j;
            if((k>=0)&&(k<Nres))
            {
                fprintf(fp_feature," 1");
            }
            else
            {
                fprintf(fp_feature," 0");
            }
        }

        for(j=-win_ss;j<=win_ss;j++)
        {
            if(j==0)
            {
                continue;
            }
            k=i+j;
            if((k>=0)&&(k<Nres))
            {
                fprintf(fp_feature2," 1");
            }
            else
            {
                fprintf(fp_feature2," 0");
            }
        }

// ccm_community
        for(j=-win_ccm;j<win_ccm;j++)
        {
            k=i+j;
            if((k>=0)&&(k<Nres))
            {
                if (ccm_community[k]>0)
                {
                    fprintf(fp_feature," %6.4f", ccm_community[k]/comm_max);
                }
                else if (ccm_community[k]<0)
                {
                    fprintf(fp_feature," %6.4f", -ccm_community[k]/comm_min);
                }
				else
                {
                    fprintf(fp_feature," %6.4f", 0.0);
                }
//                fprintf(fp_feature," %6.4f", ccm_community[k]);
            }
            else 
            {
                fprintf(fp_feature," %6.4f", 0.0);
            }
        }

        for(j=0;j<win_ccm2;j++)
        {
//            fprintf(fp_feature2," %6.4f", ccm_comm[j]);
            if(ccm_comm[j]>0)
            {
                fprintf(fp_feature2," %6.4f", ccm_comm[j]/comm_max);
            }
            else if(ccm_comm[j]<0)
            {
                fprintf(fp_feature2," %6.4f", -ccm_comm[j]/comm_min);
            }
			else 
            {
                fprintf(fp_feature2," %6.4f", 0.0);
            }

        }

//community max min
        fprintf(fp_feature," %6.4f", comm_max);
        fprintf(fp_feature," %6.4f", comm_min);

//PAS
        for(j=-win_PAS;j<=win_PAS;j++)
        {
            k=j+win_PAS;
            fprintf(fp_feature," %6.4f", PASaveN[i][k]);
        }

        for(j=-win_PAS;j<=win_PAS;j++)
        {
            k=j+win_PAS;
            fprintf(fp_feature," %6.4f", PASaveC[i][k]);
        }
        for(j=0;j<win_PAS2;j++)
        {
            fprintf(fp_feature2," %6.4f", PASaveN2[i][j]);
        }

        for(j=0;j<win_PAS2;j++)
        {
            fprintf(fp_feature2," %6.4f", PASaveC2[i][j]);
        }

//PASmax
        fprintf(fp_feature," %6.4f", PASmaxN[i]);
        fprintf(fp_feature," %6.4f", PASmaxC[i]);

//position
        fprintf(fp_feature," %6.4f", (i+1)/(double) Nres);
        fprintf(fp_feature2," %6.4f", (i+1)/(double) Nres);

        for(j=-win_prof;j<=win_prof;j++)
        {
            k=i+j;
            if((k>=0)&&(k<Nres))
            {
                for(m=1;m<Ntype;m++)
                {
                    fprintf(fp_feature," %6.4f", profile[k][m]);
                    fprintf(fp_feature2," %6.4f", profile[k][m]);
                }
            }
            else
            {
                for(m=1;m<Ntype;m++)
                {
                    fprintf(fp_feature," %6.4f", 0.0);
                    fprintf(fp_feature2," %6.4f", 0.0);
                }
            }
        }

        for(j=-win_ss;j<=win_ss;j++)
        {
            k=i+j;
            if((k>=0)&&(k<Nres))
            {
                fprintf(fp_feature," %6.4f %6.4f %6.4f", ss2[k][0],ss2[k][1],ss2[k][2]);
                fprintf(fp_feature2," %6.4f %6.4f %6.4f", ss2[k][0],ss2[k][1],ss2[k][2]);
            }
            else
            {
                fprintf(fp_feature," %6.4f %6.4f %6.4f", 0.0,0.0,0.0);
                fprintf(fp_feature2," %6.4f %6.4f %6.4f", 0.0,0.0,0.0);
            }
        }
        for(j=-win_sa;j<=win_sa;j++)
        {
            k=i+j;
            if((k>=0)&&(k<Nres))
            {
                fprintf(fp_feature," %6.4f %6.4f %6.4f %6.4f",sa2[k][0],sa2[k][1],sa2[k][2],sa2[k][3]);
                fprintf(fp_feature2," %6.4f %6.4f %6.4f %6.4f",sa2[k][0],sa2[k][1],sa2[k][2],sa2[k][3]);
            }
            else
            {
                fprintf(fp_feature," %6.4f %6.4f %6.4f %6.4f",0.0,0.0,0.0,0.0);
                fprintf(fp_feature2," %6.4f %6.4f %6.4f %6.4f",0.0,0.0,0.0,0.0);
            }
        }


        Nter=((double) i) /1000.0;
        if(Nter>1.0)
        {
            Nter=1.0;
        }
        Cter=((double) (Nres-i-1)) /1000.0;
        if(Cter>1.0)
        {
            Cter=1.0;
        }
        NCter=((double) Nres) / 1000.0;
        if(NCter>1.0)
        {
            NCter=1.0;
        }

        fprintf(fp_feature," %6.4f %6.4f %6.4f",dmsa,dpas,rpas);
        fprintf(fp_feature," %6.4f %6.4f %6.4f",Nter,Cter,NCter);
        fprintf(fp_feature," %6.4f %6.4f %6.4f",Nsite[i][0]/22.0,Nsite[i][1]/11.0,Nsite[i][2]/11.0);

        fprintf(fp_feature," %6.4f %6.4f %6.4f",msig5[i][0],msig5[i][1],msig5[i][2]);
        fprintf(fp_feature," %6.4f %6.4f %6.4f",psig5[i][0],psig5[i][1],psig5[i][2]);

        fprintf(fp_feature,"\n");

        fprintf(fp_feature2," %6.4f %6.4f %6.4f",dmsa,dpas,rpas);
        fprintf(fp_feature2," %6.4f %6.4f %6.4f",Nter,Cter,NCter);
        fprintf(fp_feature2," %6.4f %6.4f %6.4f",Nsite[i][0]/22.0,Nsite[i][1]/11.0,Nsite[i][2]/11.0);

        fprintf(fp_feature2," %6.4f %6.4f %6.4f",msig5[i][0],msig5[i][1],msig5[i][2]);
        fprintf(fp_feature2," %6.4f %6.4f %6.4f",psig5[i][0],psig5[i][1],psig5[i][2]);

        fprintf(fp_feature2,"\n");

    }

    fclose(fp_feature) ;
    fclose(fp_feature2) ;


    sprintf(out_PAS,"%s_PAS.txt",target);
    write_PAS(PAS,out_PAS,Nres) ;

    sprintf(out_PAS2,"%s_PAS2.txt",target);
    write_PAS(PAS2,out_PAS2,Nres) ;

    sprintf(out_PAS3,"%s_PAS3.txt",target);
    write_PAS(PAS3,out_PAS3,Nres) ;


    sprintf(out_ccm,"result_ccm2.txt");
    write_ccm(ccmpred,ccm_cut,out_ccm,Nres);
    sprintf(out_community,"community_ccm2.txt");
    write_community(ccm_community,out_community,Nres);

    print_time(time_0);

    free_cstring(seq);
    
    free_iarray1(seqn);
    free_darray1(ccm_community);
    free_darray1(ccm_comm);
    free_darray1(PASmaxN);
    free_darray1(PASmaxC);

    free_iarray2(msa);
    free_iarray2(msa2);
    free_iarray2(Nsite);

    free_darray2(ss2);
    free_darray2(sa2);
    free_darray2(profile);
    free_darray2(ccmpred);

    free_darray2(msig5);
    free_darray2(psig5);

    free_darray2(PAS);
    free_darray2(PAS2);
    free_darray2(PAS3);

    free_darray2(PASaveN);
    free_darray2(PASaveC);
    free_darray2(PASaveN2);
    free_darray2(PASaveC2);


}



