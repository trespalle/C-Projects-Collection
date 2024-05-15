#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define K0 50
#define N0 10000
#define NormRANu (2.3283063671E-10F)
#define pi 3.141592653F
#define N_part 32///number of monomers
#define m 1.0F ///mass of each monomer
#define T 1.0F ///temperature
#define b 1.0F ///equilibrium length of each bond
#define gamma 0.2F ///the initial deviation of each bond's length from the equilibrium length is gamma*b
///(gamma is a fraction between 0 and 1)
#define tolerancia 0.0001F
#define k_e 1000.0F ///spring constant of each bond


unsigned int larueda[256], number = 0;
unsigned char ind_ran, i1, i2, i3;


void ini_ran(int seed);
double Rand(void);

double rand_gauss(void);
void med_var(double *serie, int N, double *mediadir, double *vardir);
double maximo (double *numeros, int N);
double minimo (double *numeros, int N);
void histograma (double *secuencia, double *histo, int N, int K, double *maxdir, double *mindir, double *deltadir);
void guarda_histograma(double *histo, int K, double min, double delta, char *name_file);
double change_in_x (double p);
double change_in_p(double x, double x_next, double p, double dist, double eta, double *fuerza_elastica_dir);
void genera_condiciones_iniciales_random(double**x, double**p);
void pon_cambios_a_cero(double **delta_x_1, double **delta_x_2, double **delta_p_1, double **delta_p_2, double *T_dir, double *V_dir, double *r_giro_dir, double *r_ee_dir);
double distancia(double *x, double *y);
double energia_cinetica(double **p); ///average kinetic energy per monomer
void calcula_cdm(double *r_cdm, double **x);
double radio_de_giro(double *r_cdm, double **x);


int main()
{

    ini_ran(123457989); /// we initialize the RNG

    int i,j, i_sim, n_pasos;
    double eta, h, periodo, max, min, delta, ep_med, ek_med, r_giro, r_ee;
    double tiempo_fisico, fuerza_elastica, dist, dummy[3], next_dummy[3], r_cdm[3];
    FILE *f;

    double **x = (double**) malloc(N_part*sizeof(double*));
    for(i=0; i<N_part; i++) x[i] = (double*) malloc(3*sizeof(double)); ///3 dimensiones
    double **p = (double**) malloc(N_part*sizeof(double*));
    for(i=0; i<N_part; i++) p[i] = (double*) malloc(3*sizeof(double)); ///3 dimensiones
    double **delta_p_1 = (double**) malloc(N_part*sizeof(double*));
    for(i=0; i<N_part; i++) delta_p_1[i] = (double*) malloc(3*sizeof(double)); ///3 dimensiones
    double **delta_x_1 = (double**) malloc(N_part*sizeof(double*));
    for(i=0; i<N_part; i++) delta_x_1[i] = (double*) malloc(3*sizeof(double)); ///3 dimensiones
    double **delta_x_2 = (double**) malloc(N_part*sizeof(double*));
    for(i=0; i<N_part; i++) delta_x_2[i] = (double*) malloc(3*sizeof(double)); ///3 dimensiones
    double **delta_p_2 = (double**) malloc(N_part*sizeof(double*));
    for(i=0; i<N_part; i++) delta_p_2[i] = (double*) malloc(3*sizeof(double)); ///3 dimensiones






    eta = 0.1; ///damping parameter
    h = 0.001; ///time step
    tiempo_fisico = 500.0; ///"physical" time duration of the simulation
    n_pasos = tiempo_fisico/h; ///number of steps of the simulation


    ///we set initial conditions:
    genera_condiciones_iniciales_random(x,p);

    f = fopen("polymer_evolution.txt", "wt");
    if(f==NULL) exit(3);

    for(i_sim=0; i_sim<n_pasos; i_sim++)
    {
        pon_cambios_a_cero(delta_x_1, delta_x_2, delta_p_1, delta_p_2, &ek_med, &ep_med, &r_giro, &r_ee);
        //printf("%d %lf\n", i_sim, ek_med);
        for(i=0; i<N_part; i++) ///we compute "delta_1" (Runge-Kutta 4 algorithm)
        {
            if(i<(N_part-1))
            {
                dist = distancia(x[i], x[i+1]);
                if(dist<tolerancia) dist = tolerancia;
                ep_med+= (0.5*k_e*(dist-b)*(dist-b)); ///we compute potential energy at the same time
            }

            for(j=0; j<3; j++)
            {
                delta_x_1[i][j] = change_in_x( p[i][j] + sqrt(2.0*eta*T*h)*rand_gauss());
                if(i<(N_part-1))
                {
                    delta_p_1[i][j] += change_in_p(x[i][j], x[i+1][j], p[i][j] + sqrt(2.0*eta*T*h)*rand_gauss(), dist, eta, &fuerza_elastica);
                    delta_p_1[i+1][j]-=fuerza_elastica;
                }
                else delta_p_1[i][j] += change_in_p(0.0, 0.0, p[i][j] + sqrt(2.0*eta*T*h)*rand_gauss(), b, eta, &fuerza_elastica);
            }

        }

        for(i=0; i<N_part; i++) ///we compute "delta_2" (Runge-Kutta 4 algorithm)
        {
            if(i<(N_part-1))
            {
                for(j=0; j<3; j++)
                {
                    dummy[j] = x[i][j] + h*delta_x_1[i][j];
                    next_dummy[j] = x[i+1][j] +h*delta_x_1[i+1][j];
                }
                dist = distancia(dummy, next_dummy);
                if(dist<tolerancia) dist = tolerancia;
            }

            for(j=0; j<3; j++)
            {
                delta_x_2[i][j] = change_in_x( p[i][j] + h*delta_p_1[i][j]);
                if(i<(N_part-1))
                {
                    delta_p_2[i][j] += change_in_p(dummy[j], next_dummy[j], p[i][j] + h*delta_p_1[i][j], dist, eta, &fuerza_elastica);
                    delta_p_2[i+1][j]-=fuerza_elastica;
                }
                else delta_p_2[i][j] += change_in_p(0.0, 0.0, p[i][j] + h*delta_p_1[i][j], b, eta, &fuerza_elastica);
            }

        }

        ///We compute the observables. We compute them before updating the system so we don't have to compute distances again.
        ep_med /=N_part;  ///potential energy per monomer
        ek_med = energia_cinetica(p); ///kinetic energy per monomer
        calcula_cdm(r_cdm, x); ///we compute the center of mass of the polymer
        r_giro = radio_de_giro(r_cdm, x); ///radius of gyration of the polymer
        fprintf(f, "%lf %lf %lf %lf %lf\n", i_sim*h, distancia(x[0], x[1]), ep_med, ek_med, r_giro); ///we print observables in output file

        ///after computing changes and observables, we update the system:

        for(i=0; i<N_part; i++) for(j=0; j<3; j++)
        {
            x[i][j]+=0.5*h*(delta_x_1[i][j] + delta_x_2[i][j]);
            p[i][j]+= (0.5*h*(delta_p_1[i][j] + delta_p_2[i][j]) + sqrt(2.0*eta*T*h)*rand_gauss());

        }

    }

    fclose(f);

    f = fopen("final_picture.txt", "wt");
    if(f==NULL)
    {
        free(x);
        free(p);
        free(delta_p_1);
        free(delta_p_2);
        free(delta_x_1);
        free(delta_x_2);
        exit(2);
    }

    for(i=0; i<N_part; i++) fprintf(f, "%lf %lf %lf\n", x[i][0],x[i][1],x[i][2]);
    fclose(f);



    free(x);
    free(p);
    free(delta_p_1);
    free(delta_p_2);
    free(delta_x_1);
    free(delta_x_2);




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

double rand_gauss(void)
{
    double ro = sqrt(-2.0*log(1.0-Rand()));
    return ro*cos(2.0*pi*Rand());
}

void med_var(double *serie, int N, double *mediadir, double *vardir)
{
    double suma, sumacuadrados;
    suma = sumacuadrados = 0.0;
    for(int i=0; i<N; i++)
    {
        suma += serie[i];
        sumacuadrados += serie[i]*serie[i];
    }

    *mediadir = suma/N;
    *vardir = (sumacuadrados - suma*suma/N)/(N-1); /// ver formula del estimador de la varianza en los apuntes
}

double maximo (double *numeros, int N)
{
    double min = numeros[0];
    for(int i=0; i<N; i++) if(numeros[i]>min) min = numeros[i];
    return min;
}


double minimo (double *numeros, int N)
{
    double min = numeros[0];
    for(int i=0; i<N; i++) if(numeros[i]<min) min = numeros[i];
    return min;
}

void histograma (double *secuencia, double *histo, int N, int K, double *maxdir, double *mindir, double *deltadir)
{
    double max, min, delta;
    int i, j;
    max = maximo(secuencia, N);
    min = minimo(secuencia, N);
    delta = (max - min)/K;

    for(j=0;j<K;j++) histo[j] = 0; ///inicializamos el hsitograma a 0

    for(i=0;i<N;i++)
        for(j=0;j<K;j++) if(((min+j*delta)<=secuencia[i])&&(secuencia[i]<=min+(j+1)*delta)) histo[j]+=1.0; ///contamos cuántos números de la secuencia están en cada casilla

    for(i=0;i<K;i++) histo[i] = histo[i]/N; ///normalizamos el histograma

    *maxdir = max;
    *mindir = min;
    *deltadir = delta;

}
void guarda_histograma(double *histo, int K, double min, double delta, char *name_file)
{
    FILE *f;
    f = fopen(name_file, "w");
    for(int i=0; i<K; i++) fprintf(f, "%lf %lf\n", min + 0.5*delta + i*delta, histo[i]); /// el +0.5*delta es para que las barras del histograma queden centradas
    fclose(f);
}

double change_in_x (double p)
{
    return p/m;
}
double change_in_p(double x, double x_next, double p, double dist, double eta, double *fuerza_elastica_dir)
{
    double fuerza_elastica = -k_e*(1.0 - (b/dist))*(x - x_next);
    double change = -eta*p/m + fuerza_elastica;
    (*fuerza_elastica_dir) = fuerza_elastica;
    return change;
}

void genera_condiciones_iniciales_random(double**x, double**p)
{
    int i, j;
    double phi, omega, rho, theta;
    FILE *f;
    f = fopen("puntos_generados.txt", "wt");
    if(f==NULL)
    {
        free(x);
        free(p);
        exit(1);
    }

    for(i=0; i<N_part; i++)
    {
        phi = 2.0*pi*Rand();
        omega = Rand();
        theta = acos(1.0-2.0*omega);
        omega = Rand();
        rho = 2.0*gamma*(3.0+gamma*gamma)*omega + (1.0-gamma)*(1.0-gamma)*(1.0-gamma); ///rho generado uniformemente en ((1-gamma)b, (1+gamma)b)
        rho = pow(rho, 1.0/3.0)*b;

        x[i][0] = rho*sin(theta)*cos(phi);
        x[i][1] = rho*sin(theta)*sin(phi);
        x[i][2] = rho*cos(theta);

        for(j=0; j<3; j++) p[i][j] = 0.0; ///todas las particulas parten del reposo

        if(i>0) for(j=0; j<3; j++) x[i][j] += x[i-1][j]; ///recordemos que las particulas forman una cadena



        fprintf(f, "%lf %lf %lf\n", x[i][0], x[i][1], x[i][2]);
    }
    fclose(f);
}

void pon_cambios_a_cero(double **delta_x_1, double **delta_x_2, double **delta_p_1, double **delta_p_2, double *T_dir, double *V_dir, double *r_giro_dir, double *r_ee_dir)
{
    int i, j;
    for(i=0; i<N_part; i++) for(j=0; j<3; j++)
        delta_p_1[i][j] = delta_p_2[i][j] = delta_x_1[i][j] = delta_x_2[i][j] = 0.0;
    (*T_dir) = (*V_dir) = (*r_giro_dir) = (*r_ee_dir) = 0.0;

}

double distancia(double *x, double *y)
{
    double dist = 0.0;
    int i;
    for(i=0; i<3; i++) dist+= (x[i]-y[i])*(x[i]-y[i]);
    dist = sqrt(dist);
    return dist;
}

double energia_cinetica(double **p) ///energía cinética por partícula, en realidad
{
    int i, j;
    double suma = 0.0;
    for(i=0; i<N_part; i++) for(j=0; j<3; j++) suma+=(p[i][j]*p[i][j]/2.0/m);
    return suma/=N_part;

}

void calcula_cdm(double *r_cdm, double **x)
{
    int i, j;
    for(i=0; i<3; i++) r_cdm[i] = 0.0;
    for(i=0; i<N_part; i++) for(j=0; j<3; j++) r_cdm[j] += x[i][j];
    for(i=0; i<3; i++) r_cdm[i] /=N_part;
}

double radio_de_giro(double *r_cdm, double **x)
{
    int i;
    double suma, dist;
    suma = 0.0;
    for(i=0; i<N_part; i++)
    {
        dist = distancia(r_cdm, x[i]);
        suma+=dist*dist;
    }
    suma/=N_part;
    return sqrt(suma);
}



