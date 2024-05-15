#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define n 4
#define NormRANu (2.3283063671E-10F)
#define TINY 1.0e-20 ///a small number
#define beta 0.2F
#define p 0.5F
#define r 0.5F
#define k1 4.0F
#define k2 4.0F
#define tolerancia_en_F 0.000001F
#define tolerancia_en_ro_punto 0.000001F
#define n_intentos 100
#define impulso
#define n_explorers 100
#define epsilon 0.05F
#define tolerancia_en_dist 0.01F

unsigned int larueda[256], number = 0;
unsigned char ind_ran, i1, i2, i3;
double alpha;

///these are functions used in the Newton-Raphson algorithm:
void calcula_las_f_s(double *ro, double *f); ///hh ch hc cc
void calcula_vector_F (double *ro, double *f, double *vector_F); ///hh hc ch cc
void calcula_jacobiano(double *ro, double **jacob, double *f); ///this function computes the jacobian of the system evaluated at a particular point
double LU_decomposition(double**matriz, int *permutaciones);
void resuelve_sistema(double **matriz, double *be, int *permutaciones);

void muestra_matriz(double**matriz, int n_filas, int n_columnas); ///this shows on screen the content of a matrix
void muestra_vector(double*vector); ///this displays the content of a vector
void ini_ran(int seed); ///this function initalizes the RNG
double Rand(void); ///uniform RNG between 0 and 1
int encuentra_punto_fijo( double *ro, double *f, double *vector_F, double **jacob, int *permutaciones); ///this applies Newton-Raphson to a particicular explorer, so it ends up at a fixed point
void impulso_random(double **explorers);
void multiplica(double**matriz, double *vector, double *result, int dim); ///this multiplies two matrices
void normaliza(double *vector, int dim); ///this normalizes a vector
double distancia(double *x, double *y, int dim); ///distance function
int estabilidad_del_punto(double *rho, double **jacob); ///this computes the stability of a fixed point
void explorers_listos_para_salir(double **explorers); ///this initializes the positions of the explorers uniformly over a 4D tetrahedron
void construye_lista_descubrimientos(int *size_of_list_dir, double **explorers, double** *descubrimientos_dir); ///this creates the "discovered fixed points" list from the final positions of the explorers
void crea_old_findings(double **explorers, double *f, double *vector_F, double **jacob, int *permutaciones, double** *old_findings_dir, int *old_size_dir);
void lo_nuevo_es_lo_viejo( double**descub, double** *old_list_dir, int size_of_list, int old_size); ///this resets the main matrices at the end of each iteration

int main()
{
    ini_ran(123456789);

    int i, j, e, n_tries, dentro;
    double suma, delta_alpha = 0.001;

    double*ro = (double*) malloc(n*sizeof(double));
    double*delta_ro = (double*) malloc(n*sizeof(double));
    double*f = (double*) malloc(n*sizeof(double));
    double*vector_F = (double*) malloc(n*sizeof(double));
    int*permutaciones = (int*)malloc(n*sizeof(int));
    double**jacob = (double**) malloc(n*sizeof(double*));
    for(i=0; i<n; i++) jacob[i] = (double*) malloc(n*sizeof(double));
    double **explorers = (double**) malloc(n_explorers*sizeof(double*));
    for(i=0; i<n_explorers; i++) explorers[i] = (double*) malloc(4*sizeof(double));
    double **old_findings = (double**) malloc(sizeof(double*));
    old_findings[0] = (double*) malloc(6*sizeof(double));
    double **descubrimientos = (double**) malloc(sizeof(double*));
    descubrimientos[0] = (double*) malloc(6*sizeof(double));

    int old_size, size_of_list;
    old_size = size_of_list = 1;


    FILE*g;
    g = fopen("Raphson_lines.txt", "wt"); ///we create the output file
    fclose(g);


    for(alpha = 0.0; alpha<=1.0; alpha+=delta_alpha)
    {

        //printf("Now, the value of alpha is %lf\n", alpha);
        crea_old_findings(explorers, f, vector_F, jacob, permutaciones, &descubrimientos, &size_of_list); ///we apply Newton-Raphson to each explorer and construct the "discovered fixed points" list


        g = fopen("Raphson_lines.txt", "at"); ///we save the data:
        for(i=0; i<size_of_list; i++)
        {
            dentro = 0;
            for(j=0; j<4; j++) if(fabs(descubrimientos[i][j]-0.5)>3.0) dentro ++; ///we only print in the output file those fixed points that lie within our region of interest.
            ///in this case, we have chosen as our "region of interest" a cube of side length 6 centered at (0.5, 0.5, 0.5, 0.5)
            if(dentro==0)
            {
                fprintf(g, "%lf ", alpha);
                for(j=0; j<5; j++) fprintf(g, "%lf ", descubrimientos[i][j]);
                fprintf(g, "\n");
            }

        }
        fclose(g);

        lo_nuevo_es_lo_viejo(descubrimientos, &old_findings, size_of_list, old_size); ///we reset our matrices
        old_size = size_of_list;
        size_of_list =1;


    }





    free(f);
    free(ro);
    free(delta_ro);
    free(vector_F);
    free(permutaciones);
    for(i=0; i<n; i++) free(jacob[i]);
    free(jacob);
    for(i=0; i<n_explorers; i++) free(explorers[i]);
    free(explorers);
    for(i=0; i<old_size; i++) free(old_findings[i]);
    free(old_findings);
    for(i=0; i<size_of_list; i++) free(descubrimientos[i]);
    free(descubrimientos);


    return 0;
}

void calcula_las_f_s(double *ro, double *f)
{
    int j;
    ///calculamos las f's
    f[0] = alpha*(ro[1]+ro[3]); ///f_1_alpha
    f[1] = alpha*(ro[2]+ro[3]); ///f_2_alpha
    f[2] = beta*(ro[0]+ro[2]); ///f_1_beta
    f[3] = beta*(ro[0]+ro[1]); ///f_2_beta

    for(j=0; j<4; j+=2)
    {
        f[j] = 1.0 - f[j];
        f[j] = pow(f[j],k1);
        f[j]= 1.0 - f[j];
    }

    for(j=1; j<4; j+=2)
    {
        f[j] = 1.0 - f[j];
        f[j] = pow(f[j],k2);
        f[j]= 1.0 - f[j];
    }

}

void calcula_vector_F (double *ro, double *f, double *vector_F)
{
    vector_F[0] = -(r*(1.0 -ro[0] -ro[1] -ro[2] -ro[3]) - (p*f[0] + (1.0-p)*f[1])*ro[0]);
    vector_F[2] = -((1.0-p)*f[1]*ro[0] - (p*f[0] + (1.0-p)*f[3])*ro[2]);
    vector_F[1] = -(p*f[0]*ro[0] - (p*f[2] + (1.0-p)*f[1])*ro[1]);
    vector_F[3] = -(p*f[0]*ro[2] + (1.0-p)*f[1]*ro[1] - (p*f[2] + (1.0-p)*f[3])*ro[3]);
    ///le meto un signo menos a cada bicho porque estoy llamando F a lo que el libro llama -F
}

void calcula_jacobiano(double *ro, double **jacob, double *f)
{

    jacob[0][0] = p*(-f[0]) - r - (-f[1])*(p - 1.0);
    jacob[0][1] = -r - alpha*ro[0]*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]),k1 - 1.0);
    jacob[0][2] = alpha*ro[0]*k2*pow(1.0 - alpha*( ro[3] +  ro[2]),k2 - 1.0)*(p - 1.0) - r;
    jacob[0][3] = - r -ro[0]*(alpha*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]), k1 - 1.0) - alpha*k2*pow(1.0 - alpha*( ro[3] +  ro[2]), k2 - 1.0)*(p - 1.0));

    jacob[2][0] = (-f[1])*(p - 1.0) + beta*ro[2]*k2*pow(1.0 - beta*( ro[1] +  ro[0]), k2 - 1.0)*(p - 1.0);
    jacob[2][1] = -ro[2]*(alpha*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]), k1 - 1.0) - beta*k2*pow(1.0 - beta*( ro[1] +  ro[0]), k2 - 1.0)*(p - 1.0));
    jacob[2][2] = p*(-f[0]) - (-f[3])*(p - 1.0) - alpha*ro[0]*k2*pow(1.0 - alpha*( ro[3] +  ro[2]), k2 - 1.0)*(p - 1.0);
    jacob[2][3] = -alpha*ro[2]*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]), k1 - 1.0) - alpha*ro[0]*k2*pow(1.0 - alpha*( ro[3] +  ro[2]), k2 - 1.0)*(p - 1.0);

    jacob[1][0] = -p*(-f[0]) - beta*ro[1]*k1*p*pow(1.0 - beta*( ro[2] +  ro[0]), k1 - 1.0);
    jacob[1][1] = p*(-f[2]) - (-f[1])*(p - 1.0) + alpha*ro[0]*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]), k1 - 1.0);
    jacob[1][2] = -ro[1]*(beta*k1*p*pow((1.0 - beta*( ro[2] +  ro[0])), (k1 - 1.0)) - alpha*k2*pow(1.0 - alpha*( ro[3] +  ro[2]), k2 - 1.0)*(p - 1.0));
    jacob[1][3] = alpha*ro[0]*k1*p*pow((1.0 - alpha*(ro[3] +  ro[1])), (k1 - 1.0)) + alpha* ro[1]*k2*pow(1.0 - alpha*( ro[3] +  ro[2]), k2 - 1.0)*(p - 1.0);

    jacob[3][0] = -ro[3]*(beta*k1*p*pow(1.0 - beta*( ro[2] +  ro[0]), k1 - 1.0) - beta*k2*pow(1.0 - beta*( ro[1] +  ro[0]), k2 - 1.0)*(p - 1.0));
    jacob[3][1] = (-f[1])*(p - 1.0) + alpha* ro[2]*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]), k1 - 1.0) + beta* ro[3]*k2*pow(1.0 - beta*( ro[1] +  ro[0]), k2 - 1.0)*(p - 1.0);
    jacob[3][2] = -p*(-f[0]) - beta*ro[3]*k1*p*pow(1.0 - beta*(ro[2] + ro[0]), k1 - 1.0) - alpha*ro[1]*k2*pow(1.0 - alpha*(ro[3] + ro[2]), k2 - 1.0)*(p - 1.0);
    jacob[3][3] = p*(-f[2]) - (-f[3])*(p - 1.0) + alpha*ro[2]*k1*p*pow(1.0 - alpha*( ro[3] +  ro[1]), k1 - 1.0) - alpha* ro[1]*k2*pow(1.0 - alpha*( ro[3] +  ro[2]), k2 - 1.0)*(p - 1.0);

}

void muestra_matriz(double**matriz, int n_filas, int n_columnas)
{
    int i,j;
    for(i=0; i<n_filas; i++)
    {
        for(j=0; j<n_columnas; j++) printf("%lf ", matriz[i][j]);
        printf("\n");
    }
    printf("\n");
}

void muestra_vector(double*vector)
{
    int i;
    for(i=0; i<n; i++) printf("%lf ", vector[i]);
    printf("\n\n");
}

double LU_decomposition(double**matriz, int *permutaciones)
{
    int i,j,e, imax;
    double maximo, dummy, suma, paridad;
    double*escalado = (double*)malloc(n*sizeof(double));

    paridad = 1.0;

    ///We create the "scaling" vector:

    for(i=0; i<n; i++)
    {
        maximo = 0.0;
        for(j=0; j<n; j++)
            if((dummy = fabs(matriz[i][j]))> maximo) maximo = dummy;

        if(fabs(maximo)<TINY)
        {
            printf("Error: singular matrix.\n");
            exit(1);
        }
        escalado[i]=1.0/maximo;
    }

    ///Here Crout's method begins:
    for(j=0; j<n; j++)
    {
        ///we start with cases i<j
        for(i=0; i<j; i++)
        {
            suma = matriz[i][j];
            for(e=0; e<i; e++) suma-=matriz[i][e]*matriz[e][j];
            matriz[i][j] = suma;
        }

        maximo = 0.0;
        for(i=j; i<n; i++)
        {
            suma = matriz[i][j];
            for(e=0; e<j; e++) suma-=matriz[i][e]*matriz[e][j];
            matriz[i][j] = suma;
            if((dummy = escalado[i]*fabs(suma))>=maximo)
            {
                maximo = dummy;
                imax = i;
            }
        }

        if(j!=imax)
        {
            for(e=0; e<n; e++)
            {
                dummy = matriz[imax][e];
                matriz[imax][e] = matriz[j][e];
                matriz[j][e] = dummy;
            }
            paridad = -paridad;
            dummy = escalado[imax];
            escalado[imax] = escalado[j];
            escalado[j] = dummy;
        }

        permutaciones[j] = imax;
        if(matriz[j][j]==0.0) matriz[j][j] = TINY;
        if(j<(n-1))
        {
            dummy = 1.0/matriz[j][j];
            for(i=j+1; i<n; i++) matriz[i][j]*=dummy;
        }

    }

    free(escalado);
    return paridad;

}
void resuelve_sistema(double **matriz, double *be, int *permutaciones)
{
    int i, ii=-1, ip, j;
    double suma;

    for(i=0; i<n; i++)
    {
        ip = permutaciones[i];
        suma = be[ip];
        be[ip] = be[i];
        if(ii>(-1)) for(j=ii; j<i; j++) suma-=matriz[i][j]*be[j];
        else if(suma) ii = i;
        be[i] = suma;
    }

    for(i=n-1; i>=0; i--)
    {
        suma = be[i];
        for(j=i+1; j<n; j++) suma-=matriz[i][j]*be[j];
        be[i] = suma/matriz[i][i];
    }
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
    double r_0;
    i1 = ind_ran-24;
    i2 = ind_ran-55;
    i3 = ind_ran-61;
    larueda[ind_ran] = larueda[i1] + larueda[i2];
    number = (larueda[ind_ran]^larueda[i3]);
    ind_ran++;
    r_0 = number*NormRANu;
    return r_0;
}



int encuentra_punto_fijo( double *ro, double *f, double *vector_F, double **jacob, int *permutaciones)
{

    double suma1, suma2, aux, modulo, modulo_test, ro_test[4], F_test[4], control;
    int i, j, vigilante, thegrefg;
    vigilante = 0;

    for(i=0; (i<n_intentos)&&(vigilante<1); i++)
    {
        calcula_las_f_s(ro, f);
        calcula_vector_F(ro, f, vector_F);

        suma1 = 0.0;
        for(j=0; j<n; j++) suma1+=fabs(vector_F[j]);
        if(suma1<tolerancia_en_F) vigilante++;
        //printf("%d\n", vigilante);

        modulo = 0.0;
        for(j=0; j<4; j++) modulo += vector_F[j]*vector_F[j];
        modulo = sqrt(modulo);

        calcula_jacobiano(ro, jacob, f);
        aux = LU_decomposition(jacob, permutaciones);
        resuelve_sistema(jacob, vector_F, permutaciones);

        suma2 = 0.0;
        for(j=0; j<n; j++) suma2+=fabs(vector_F[j]);
        if(suma2<tolerancia_en_ro_punto) vigilante++;

        ///control
        control = 1.0;
        thegrefg = 0;
        do
        {
            for(j=0; j<4; j++) ro_test[j] = ro[j] + control*vector_F[j];
            calcula_las_f_s(ro_test, f);
            calcula_vector_F(ro_test, f, F_test);
            modulo_test = 0.0;
            for(j=0; j<4; j++) modulo_test += F_test[j]*F_test[j];
            modulo_test = sqrt(modulo_test);
            control/=2.0;
            thegrefg++;
        } while((modulo_test>modulo)&&(thegrefg<100));

        if(thegrefg>999) printf("thegrefg\n");


        for(j=0; j<n; j++) ro[j] = ro_test[j];

        //muestra_vector(ro);
        //printf("%d\n", vigilante);

    }

    return i;
}

void impulso_random(double **explorers)
{
    int i, j;
    for(i=0; i<n_explorers; i++) for(j=0; j<4; j++) explorers[i][j]+= epsilon*(Rand()-0.5);

}

void multiplica(double**matriz, double *vector, double *result, int dim)
{
    int i, j;
    double suma;
    for(i=0; i<dim; i++)
    {
        suma = 0.0;
        for(j=0; j<dim; j++) suma+= matriz[i][j]*vector[j];
        result[i] = suma;
    }
}

void normaliza(double *vector, int dim)
{
    int i;
    double norma = 0.0;
    for(i=0; i<dim; i++) norma+= vector[i]*vector[i];
    norma = sqrt(norma);

    for(i=0; i<dim; i++) vector[i]/=norma;
}

double distancia(double *x, double *y, int dim)
{
    int i;
    double suma = 0.0;
    for(i=0; i<dim; i++) suma+= (x[i]-y[i])*(x[i]-y[i]);
    suma = sqrt(suma);
    if(suma<tolerancia_en_dist) suma = tolerancia_en_dist/2.0;
    return suma;
}

int estabilidad_del_punto(double *rho, double **jacob) ///the jacobian has to be already computed
{
    int i, contador, ind_max;
    double dummy_vect[4], new_dummy[4], dist, autovalor;

    for(i=0; i<4; i++) dummy_vect[i] = rho[i]; ///any arbitrary vector could be taken as the initial one
    contador = 0;
    normaliza(dummy_vect,4);

    do
    {
        multiplica(jacob, dummy_vect, new_dummy, 4);
        normaliza(new_dummy, 4);
        dist = distancia(dummy_vect, new_dummy, 4);
        contador++;
        for(i=0; i<4; i++) dummy_vect[i] = new_dummy[i];
    } while((contador<10000)&&(dist>0.0001));

    if(contador>10000) return 0;

    multiplica(jacob, dummy_vect, new_dummy,4);
    ind_max = 0;
    for(i=0; i<4; i++) if(fabs(new_dummy[i])>fabs(new_dummy[ind_max])) ind_max = i;

    autovalor = new_dummy[ind_max]/dummy_vect[ind_max];
    //printf("%lf\n", autovalor);

    if(autovalor<0.0) return -1;
    else return 1;

}

void explorers_listos_para_salir(double **explorers)
{
    int i;
    double omega;

    for(i=0; i<n_explorers; i++)
    {
        omega = Rand();
        explorers[i][0] = 1.0 - pow(1.0-omega, 1.0/4.0);
        omega = Rand();
        explorers[i][1] = (1.0-explorers[i][0])*(1.0 - pow(1.0-omega, 1.0/3.0));
        omega = Rand();
        explorers[i][2] = (1.0-explorers[i][0]-explorers[i][1])*(1.0-sqrt(1.0-omega));
        explorers[i][3] = (1.0-explorers[i][0]-explorers[i][1]-explorers[i][2])*Rand();
    }
}

void construye_lista_descubrimientos(int *size_of_list_dir, double **explorers, double** *descubrimientos_dir)
{
    int i, j, size_of_list = 1;
    double dist;

    (*descubrimientos_dir) = (double**) realloc((*descubrimientos_dir), sizeof(double*));
    (*descubrimientos_dir)[0] = (double*) malloc(6*sizeof(double));

    //printf("pepe\n");

    for(i=0; i<4; i++) (*descubrimientos_dir)[0][i] = explorers[0][i];

    for(i=0; i<n_explorers; i++)
    {
        j=0;
        do
        {
            dist = distancia((*descubrimientos_dir)[j], explorers[i], 4);
            j++;
        } while((dist>tolerancia_en_dist)&&(j<size_of_list));

        if((j==size_of_list)&&(dist>tolerancia_en_dist)&&(explorers[i][0]<950.0))
        {
            size_of_list++;
            (*descubrimientos_dir) = (double**) realloc((*descubrimientos_dir), size_of_list*sizeof(double*));
            (*descubrimientos_dir)[size_of_list-1] = (double*) malloc(6*sizeof(double));

            for(j=0; j<4; j++) (*descubrimientos_dir)[size_of_list-1][j] = explorers[i][j];
            ///muestra_matriz((*descubrimientos_dir), size_of_list, 6);
        }
    }

    (*size_of_list_dir) = size_of_list;
}

void crea_old_findings(double **explorers, double *f, double *vector_F, double **jacob, int *permutaciones, double** *old_findings_dir, int *old_size_dir)
{
    int i,j, n_tries;

    //explorers_listos_para_salir(explorers);
    for(i=0; i<n_explorers; i++) for(j=0; j<4; j++) explorers[i][j] = Rand();

    for(i=0; i<n_explorers; i++)
    {
        n_tries = encuentra_punto_fijo(explorers[i], f, vector_F, jacob, permutaciones);
        if(n_tries==n_intentos) explorers[i][0] = 1000.0; ///these are the "failed" explorers that got too far from the region of interest

    }

    construye_lista_descubrimientos(old_size_dir, explorers, old_findings_dir);
    for(i=0; i<(*old_size_dir); i++)
    {
        calcula_las_f_s((*old_findings_dir)[i], f);
        calcula_jacobiano((*old_findings_dir)[i], jacob, f);
        (*old_findings_dir)[i][4] = (double) estabilidad_del_punto((*old_findings_dir)[i], jacob);
        (*old_findings_dir)[i][5] = (double)i;
    }

}



void lo_nuevo_es_lo_viejo( double**descub, double** *old_list_dir, int size_of_list, int old_size)
{
    int i,j;
    for(i=0; i<old_size; i++) free((*old_list_dir)[i]);
    free((*old_list_dir));

    (*old_list_dir) = (double**) malloc(size_of_list*sizeof(double*));
    for(i=0; i<size_of_list; i++) (*old_list_dir)[i] = malloc(6*sizeof(double));

    for(i=0; i<size_of_list; i++) for(j=0; j<6; j++) (*old_list_dir)[i][j] = descub[i][j];

}






