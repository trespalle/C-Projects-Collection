#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define dim 2
#define z 4


unsigned int larueda[256], number = 0;
unsigned char ind_ran, i1, i2, i3;


void ini_ran(int seed);
double Rand(void);

void read_params(int *L_dir, double *T_dir, int *n_mcs_dir, int *n_meas_dir, int *seed_dir, char *name_file);
int the_size_of(char *name_file);
void read_time_series(double *eners, double *magnets, int N, char* namefile);
int closest_smaller_power_of_two(int n); ///devuelve la parte entera del logaritmo en base 2 de n
void group_by_pairs(double *data, int N, double *binned_data);
void med_var(double *data, int N, double *med_dir, double *var_dir);
void binning_analysis(double *binned_data, double *auxiliar, int n_bins_ini, double **stats, int n_binnings);

int main()
{

    double T, factor;
    int i,j, L, N, n_mcs, n_meas, seed, N_medidas, N_data, n_binnings, bin_size;

    ///we read the relevant parameters: L, T, n_mcs, n_meas, seed
    read_params(&L, &T, &n_mcs, &n_meas, &seed, "input.txt");
    N_medidas = n_mcs/n_meas;
    N = L*L;

    ///we read the Monte Carlo data
    N_data = the_size_of("evolution.txt");
    if(N_data != N_medidas) printf("Warning: length of the timeseries differs from N_medidas.\n");
    n_binnings = closest_smaller_power_of_two(N_data);


    double *eners = (double*) malloc(N_data*sizeof(double));
    double *magnets = (double*) malloc(N_data*sizeof(double));
    double *auxiliar = (double*) malloc(N_data*sizeof(double));

    read_time_series(eners, magnets, N_data, "evolution.txt");
    for(i=0; i<N_data; i++) magnets[i] = fabs(magnets[i]);



    ///we erase thermalization data
    int N_term = 2*L*L; ///we assume that the thermalization time grows with the lattice size
    //printf("nterm = %d\n", N_term);

    for(i=0; i<N_data; i++) auxiliar[i] = eners[i]; ///we readjust energy vector size
    free(eners);
    eners = (double*) malloc((N_data-N_term)*sizeof(double));
    for(i=N_term; i<N_data; i++) eners[i-N_term] = auxiliar[i];

    for(i=0; i<N_data; i++) auxiliar[i] = magnets[i]; ///we readjust magnetizaction vector size
    free(magnets);
    magnets = (double*) malloc((N_data-N_term)*sizeof(double));
    for(i=N_term; i<N_data; i++) magnets[i-N_term] = auxiliar[i];

    free(auxiliar);
    auxiliar = (double*) malloc((N_data-N_term)*sizeof(double)); ///we readjust auxiliar vector size

    N_data = N_data-N_term;
    printf("ndata = %d\n", N_data);

    //for(i=0; i<N_data; i++) printf("%lf %lf\n", eners[i], magnets[i]);




    ///BINNING ANALYSIS

    double *binned_data = (double*) malloc(N_data*sizeof(double));
    double **ener_stats = (double**) malloc(n_binnings*sizeof(double*));
    for(i=0; i<n_binnings; i++) ener_stats[i] = (double*) malloc(2*sizeof(double));
    double **magnet_stats = (double**) malloc(n_binnings*sizeof(double*));
    for(i=0; i<n_binnings; i++) magnet_stats[i] = (double*) malloc(2*sizeof(double));

    ///energy:
    for(i=0; i<N_data; i++) binned_data[i] = eners[i];
    binning_analysis(binned_data, auxiliar, N_data, ener_stats, n_binnings);

    ///magnetization:
    for(i=0; i<N_data; i++) binned_data[i] = magnets[i];
    binning_analysis(binned_data, auxiliar, N_data, magnet_stats, n_binnings);


    ///we save the output
    FILE *f = fopen("error_vs_binning_size.txt", "wt");
    if(f==NULL) exit(5);
    bin_size = 1;
    for(i=0; i<n_binnings; i++)
    {
        fprintf(f, "%d %lf %lf %lf %lf\n", bin_size, ener_stats[i][0], ener_stats[i][1], magnet_stats[i][0], magnet_stats[i][1]);
        bin_size*=2;
    }
    fclose(f);




    free(eners);
    free(magnets);
    free(auxiliar);
    free(binned_data);
    for(i=0; i<n_binnings; i++) free(ener_stats[i]);
    for(i=0; i<n_binnings; i++) free(magnet_stats[i]);
    free(ener_stats);
    free(magnet_stats);




    return 0;
}


void ini_ran(int seed)
{
    int INI, FACTOR, SUM, i;
    INI = seed;
    FACTOR = 67397;
    SUM = 7364893;

    for(i=0; i<256; i++)
    {
        INI = (INI*FACTOR + SUM);
        larueda[i] = INI;
    }
    ind_ran = i1 = i2 = i3 = 0;
}

double Rand(void)
{
    double r;
    i1 = ind_ran-24;
    i2 = ind_ran-55;
    i3 = ind_ran-61;
    larueda[ind_ran] = larueda[i1] + larueda[i2];
    number = (larueda[ind_ran]^larueda[i3]);
    ind_ran++;
    r = number*NormRANu;
    return r;
}



void read_params(int *L_dir, double *T_dir, int *n_mcs_dir, int *n_meas_dir, int *seed_dir, char *name_file)
{
    FILE* f;
    f = fopen(name_file, "rt");
    if(f==NULL) exit(1);
    fscanf(f, "%d\n%lf\n%d\n%d\n%d\n", L_dir, T_dir, n_mcs_dir, n_meas_dir, seed_dir);
    fclose(f);

}

int the_size_of(char *name_file)
{
    int dum_int, size = 0;
    double dummy;
    FILE* f = fopen(name_file, "rt");
    if(f==NULL) exit(2);
    while(feof(f)==0)
    {
        fscanf(f,"%d %lf %lf\n", &dum_int, &dummy, &dummy);
        size++;
    }
    fclose(f);
    return size;
}

void read_time_series(double *eners, double *magnets, int N, char* namefile)
{
    int i, dum_int;
    FILE *f = fopen(namefile,"rt");
    if(f==NULL) exit(3);
    for(i=0; i<N; i++) fscanf(f, "%d %lf %lf\n", &dum_int, eners + i, magnets + i);
    fclose(f);
}

int closest_smaller_power_of_two(int n) ///devuelve la parte entera del logaritmo en base 2 de n
{
    int thegrefg, power;
    power = 2;
    thegrefg = 0;
    if(n<=power)
    {
        printf("Tu time series tiene menos de tres datos.");
        exit(3);
    }
    while(power<n)
    {
        power*=2;
        thegrefg++;

    }
    return thegrefg;
}

void group_by_pairs(double *data, int N, double *binned_data)
{
    int i, j;
    j=0;
    for(i=0; i<(N-1); i+=2)
    {
        binned_data[j] = 0.5*(data[i]+data[i+1]);
        j++;
    }
}

void med_var(double *data, int N, double *med_dir, double *var_dir)
{
    int i;
    double med, suma2, factor;

    factor = N;
    factor/=(N-1);

    med = suma2 = 0.0;

    for(i=0; i<N; i++)
    {
        med+=data[i];
        suma2+=data[i]*data[i];
    }
    med/=N;
    suma2/=N;

    *med_dir = med;
    *var_dir = factor*(suma2-med*med);

}

void binning_analysis(double *binned_data, double *auxiliar, int n_bins_ini, double **stats, int n_binnings)
{
    int i,j, n_bins = n_bins_ini;
    for(i=0; i<n_binnings; i++)
    {
        med_var(binned_data, n_bins, &stats[i][0], &stats[i][1]);
        stats[i][1] = stats[i][1]/n_bins; ///esta es la varianza del promedio

        for(j=0; j<n_bins; j++) auxiliar[j] = binned_data[j];
        group_by_pairs(auxiliar, n_bins, binned_data);
        n_bins/=2;
    }
}













