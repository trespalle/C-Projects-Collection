#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NormRANu (2.3283063671E-10F)
#define normNorm (4.656612873E-10F) ///1/2^31, so that the value of r1279 is never 1
#define NBITM1 31
#define dim 2
#define z 4


unsigned int larueda[256], number = 0;
unsigned char ind_ran, i1, i2, i3;

int secuencia_rand[2048], index1[2048], index2[2048], iaux;


void ini_ran(int seed);
double Rand(void);

void crea_lista_vecinos(int** vecinos, int L, int N);
void muestra_lista(int **vecinos, int N);
void genera_output(double* eners, double* magnets, int n_meas, int N_medidas, char *name_file);
void lee_input(int *L_dir, double *T_dir, int *n_mcs_dir, int *n_meas_dir, int *seed_dir, char *name_file);
void measure_observables(char *spins, int N, int **vecinos, double *eners, double *magnets, int index);
void lee_configuracion(char *spins, int N, char *name_file);
void mc_step(char *spins, int N, int **vecinos, double *probs);
void ini_r1279();
double r1279(void);



int main()
{

    double T, factor;
    int i,j, L, N, n_mcs, n_meas, seed, N_medidas;

    ///we read the input: L, T, n_mcs, n_meas, seed
    lee_input(&L, &T, &n_mcs, &n_meas, &seed, "input.txt");
    N_medidas = n_mcs/n_meas;
    N = L*L;

    ///we initialize the RNG using the seed
    ini_ran(seed);
    ini_r1279();

    ///we create the array of spins, the neigbours list (vecinos), y and the vectors that store the measurements of m and E
    char* spins = (char*) malloc(N*sizeof(char));
    int **vecinos = (int**) malloc(N*sizeof(int*));
    for(i=0; i<N; i++)  vecinos[i] = (int*) malloc(z*sizeof(int));
    double *eners = (double*) malloc(N_medidas*sizeof(double));
    double *magnets = (double*) malloc(N_medidas*sizeof(double));


    ///we initialize the lattice geometry:
    crea_lista_vecinos(vecinos, L, N);
    //muestra_lista(vecinos, N);
    //lee_configuracion(spins, N, "mys3.txt");
    for(i=0; i<N; i++) spins[i] = 2*floor(2.0*Rand())-1; ///random initial configuration
    //for(i=0; i<N; i++) printf("%d\n", spins[i]);


    ///we create the probability vector
    double* probs = (double*) malloc((2*z+1)*sizeof(double));
    factor = exp(-2.0/T);
    probs[z] = 1.0;
    for(i=1; i<=z; i++)
    {
        probs[z-i] = probs[z-i+1]/factor;
        probs[z+i] = probs[z+i-1]*factor;
    }
    //for(i=0; i<(2*z+1); i++) printf("%lf\n", probs[i]);



    ///we run the simulation
    j=0;
    for(i=0; i<n_mcs; i++)
    {
        mc_step(spins, N, vecinos, probs);
        if((i%n_meas) == 0)
        {
            measure_observables(spins, N, vecinos, eners, magnets, j);
            j++;
        }
    }

    for(i=0; i<N; i++) printf("%d\n", spins[i]);



    ///we generate the output:
    genera_output(eners, magnets, n_meas, N_medidas, "evolution.txt");




    for(i=0; i<N; i++) free(vecinos[i]);
    free(vecinos);
    free(spins);
    free(eners);
    free(magnets);
    free(probs);



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

void crea_lista_vecinos(int** vecinos, int L, int N)
{
    int x, y, index, aux;
    for(x=0; x<L; x++) for(y=0; y<L; y++)
    {
        index = x + L*y;
        vecinos[index][0] = index + L; ///vecino de arriba
        vecinos[index][1] = index + 1; ///vecino de la derecha
        vecinos[index][2] = index - L; ///vecino de abajo
        vecinos[index][3] = index -1; ///vecino de la izquierda
    }

    for(index=0; index<N; index++)
    {
        printf("%d\n",vecinos[index][2]);
        printf("%lf\n", floor(((double)vecinos[index][2])/N)*N);
        vecinos[index][0] -= floor(((double)vecinos[index][0])/N)*N;
        vecinos[index][2] = vecinos[index][2] - floor(((double)vecinos[index][2])/N)*N;
        printf("%d\n",vecinos[index][2]);


        aux = vecinos[index][1]-floor((double)index/L)*L;
        aux-= floor((double)aux/L)*L;
        vecinos[index][1] = floor((double)index/L)*L + aux;

        aux = vecinos[index][3]-floor((double)index/L)*L;
        aux-= floor((double)aux/L)*L;
        vecinos[index][3] = floor((double)index/L)*L + aux;

    }
}

void muestra_lista(int **vecinos, int N)
{
    int i,j;
    for(i=0; i<N; i++)
    {
        printf("%d   ", i);
        for(j=0; j<z; j++) printf("%d ", vecinos[i][j]);
        printf("\n");
    }
}

void genera_output(double* eners, double* magnets, int n_meas, int N_medidas, char *name_file)
{
    int i;
    FILE *f = fopen(name_file, "wt");
    if(f==NULL) exit(2);
    for(i=0; i<N_medidas; i++) fprintf(f, "%d %lf %lf\n", i*n_meas, eners[i], magnets[i]);
    fclose(f);
}

void lee_input(int *L_dir, double *T_dir, int *n_mcs_dir, int *n_meas_dir, int *seed_dir, char *name_file)
{
    FILE* f;
    f = fopen(name_file, "rt");
    if(f==NULL) exit(1);
    fscanf(f, "%d\n%lf\n%d\n%d\n%d\n", L_dir, T_dir, n_mcs_dir, n_meas_dir, seed_dir);
    fclose(f);

}

void measure_observables(char *spins, int N, int **vecinos, double *eners, double *magnets, int index)
{
    int i, e;
    double energia, imanacion;
    energia = imanacion = 0.0;
    for(i=0; i<N; i++)
    {
        imanacion += spins[i];
        //eneriga-=spins[i]*(spins[vecinos[i][0]]+spins[vecinos[i][1]]); ///optimizado
        for(e=0; e<z; e++) energia-=0.5*spins[i]*spins[vecinos[i][e]]; ///generalizado
    }
    eners[index] = energia;
    magnets[index] = imanacion;
}

void lee_configuracion(char *spins, int N, char *name_file)
{
    int dummy, i;

    FILE *f = fopen(name_file, "rt");
    if(f==NULL) exit(3);
    i=0;
    while(feof(f)==0)
    {
        fscanf(f, "%d %d\n", &dummy, &spins[i]);
        i++;
    }

    fclose(f);

    if(i!=N)
    {
        printf("El numero de espines del archivo no coincide con N.\nError fatal.");
        exit(4);
    }


}

void mc_step(char *spins, int N, int **vecinos, double *probs)
{
    int i,j, index, ind2, local_field;

    for(i=0; i<N; i++) ///proponemos N updates
    {
        do
        {
            index = N*r1279();
        } while(index == N);

        local_field = 0;
        for(j=0; j<z; j++) local_field += spins[vecinos[index][j]];
        ind2 = spins[index]*local_field + z;
        if(probs[ind2]>1.0) spins[index] = -spins[index];
        else if(r1279()<probs[ind2]) spins[index] = -spins[index];

    }
}

void ini_r1279()
{
    int i,j, one_bit;

    iaux = 0;
    for(i=0; i<2048; i++)
    {
        index1[i] = (i-1279)&2047;
        index2[i] = (i-418)&2047;

        secuencia_rand[i] = 0;
        for(j=0; j<=NBITM1; j++)
        {
            one_bit = 0;
            if(Rand()>0.5) one_bit = 1;

            one_bit = one_bit << j; ///we shift j places to the left
            secuencia_rand[i] = secuencia_rand[i]|one_bit; ///hacemos un or bit a bit
        }
        secuencia_rand[i] = 2*secuencia_rand[i] + 1;
    }
}

double r1279(void)
{
    int numerin;
    iaux = (iaux+1)&2047;
    secuencia_rand[iaux] = secuencia_rand[index1[iaux]]*secuencia_rand[index2[iaux]];
    numerin = (secuencia_rand[iaux]>>1)&2147483647;
    return ((double)numerin)*normNorm;
}












