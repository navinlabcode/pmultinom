#include <gsl/gsl_sf_gamma.h>
#include <math.h>

/* Poisson P(X >= x) */
double pgreater(double lambda, int x)
{
    return gsl_sf_gamma_inc_P((double) x, lambda);
}

/* Poisson P(X=x) */
double pequal(double lambda, int x)
{
    return exp(-lambda + (double) x * log(lambda)
               - gsl_sf_lngamma((double) x + 1));
}

/* One iteration of Newton root-finding */
void newtonUpdate(
    /* To define the update function, need for each clone probabilities,
    cutoffs and multiplicities. Also the target. */
    double npop, double ps[], int cs[], int ms[], double target,
    /* The resulting function updates a pair, containing the current value of
    the parameter, and how much it changed in this iteration, so you can
    decide when to terminate */
    double *lambda, double *dl)
{
    double funcValue = 0;
    double funcDeriv = 0;
    int i;
    double p;
    int c;
    int m;
    double ratio;

    for (i=0; i < npop; i++)
    {
        p = ps[i];
        c = cs[i];
        m = ms[i];

        funcValue += (double) m * log(pgreater(*lambda*p, c));
        funcDeriv += (double) m * pequal(*lambda*p, c-1) /
                        pgreater(*lambda*p, c);
    }
    funcValue -= log(target);

    ratio = funcValue / funcDeriv;

    *lambda = *lambda - ratio;
    *dl = ratio;
}

/* Newton root-finding, to find the number of individuals which must be
 * sampled to have a probability of success equal to some target */
void partialNewtonRoot(
    /* Defining the root-finding function requires the target probability, and
    the desired precision */
    double target, double eps,
    /* The root-finding function is made to be folded. It takes in a list of
    information on each subpopulation, and a Poisson rate which works for
    those subpopulations. Then, as its second argument, it takes information
    about an additional subpopulation. From that information, it calculates a
    new Poisson rate. */
    double ps[], int cs[], int ms[], int k, double *lambda)
{
    double dl;
    dl = 1000 * eps;
    while (fabs(dl) > eps)
    {
        newtonUpdate(k, ps, cs, ms, target, lambda, &dl);
    }
}

/* A function which sequentially does Newton root-finding. Given a list of
 * populations and conditions, sorted from hardest to easiest to satisfy, it
 * gives a Poisson rate that will satisfy all of them simultaneously. */
double sequentialNewtonRoot(
    /* The information about each population: probability, cutoff, and
    multiplicity. */
        int npop, double ps[], int cs[], int ms[],
    /* The lowest Poisson rate, which we need as a starting value */
        double lambdaStart,
    /* The target probability and required precision */
        double target, double eps)
{
    int i;
    double lambda;
    lambda = lambdaStart;
    /* i will be interpreted as the number of populations in the calculation.
     * When the number of populations is 1, we already know the answer. So it
     * starts at 2. */
    for (i=2; i <= npop; i++)
    {
        partialNewtonRoot(target, eps, ps, cs, ms, i, &lambda);
    }
    return lambda;
}

/* The interface to R for sequentialNewtonRoot */
void sequentialNewtonRootR(
        int *npop, double ps[], int cs[], int ms[],
        double *lambdaStart,
        double *target, double *eps,
        double *output)
{   
    *output =  sequentialNewtonRoot(*npop, ps, cs, ms, *lambdaStart, *target, *eps); 
}
