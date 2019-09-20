#include <math.h>

double Power(double a, int b)
{
    double current = 1;
    int i;
    for (i=0; i < b; i++)
    {
        current = current * a;
    }
    return current;
}

void calculate_correction(int *cumulant_ptr, int *c_ptr, double *ld_ptr, double *x_ptr, double *correction_ptr)
{
    int cumulant = *cumulant_ptr;
    int c = *c_ptr;
    double ld = *ld_ptr;
    double x = *x_ptr;
    double correction;

    if (cumulant == 1)
    {
        correction = x;
    }
    
    if (cumulant == 2)
    {
        correction = c*x - ld*x - Power(x,2);
    }
    
    if (cumulant == 3)
    {
        correction =
            Power(c,2)*x - ld*x - 2*c*ld*x + Power(ld,2)*x - 
            3*c*Power(x,2) + 3*ld*Power(x,2) + 2*Power(x,3);
    }
    
    if (cumulant == 4)
    {
        correction = 
          Power(c,3)*x - ld*x - 3*c*ld*x - 3*Power(c,2)*ld*x + 
          3*Power(ld,2)*x + 3*c*Power(ld,2)*x - Power(ld,3)*x - 
          7*Power(c,2)*Power(x,2) + 4*ld*Power(x,2) + 
          14*c*ld*Power(x,2) - 7*Power(ld,2)*Power(x,2) + 
          12*c*Power(x,3) - 12*ld*Power(x,3) - 6*Power(x,4);
    }
    
    if (cumulant == 5)
    {
        correction = 
         Power(c,4)*x - ld*x - 4*c*ld*x - 6*Power(c,2)*ld*x - 
           4*Power(c,3)*ld*x + 7*Power(ld,2)*x + 12*c*Power(ld,2)*x + 
           6*Power(c,2)*Power(ld,2)*x - 6*Power(ld,3)*x - 
           4*c*Power(ld,3)*x + Power(ld,4)*x - 
           15*Power(c,3)*Power(x,2) + 5*ld*Power(x,2) + 
           25*c*ld*Power(x,2) + 45*Power(c,2)*ld*Power(x,2) - 
           25*Power(ld,2)*Power(x,2) - 45*c*Power(ld,2)*Power(x,2) + 
           15*Power(ld,3)*Power(x,2) + 50*Power(c,2)*Power(x,3) - 
           20*ld*Power(x,3) - 100*c*ld*Power(x,3) + 
           50*Power(ld,2)*Power(x,3) - 60*c*Power(x,4) + 
           60*ld*Power(x,4) + 24*Power(x,5);
    }
    
    if (cumulant == 6)
    {
        correction = 
         Power(c,5)*x - ld*x - 5*c*ld*x - 10*Power(c,2)*ld*x - 
           10*Power(c,3)*ld*x - 5*Power(c,4)*ld*x + 15*Power(ld,2)*x + 
           35*c*Power(ld,2)*x + 30*Power(c,2)*Power(ld,2)*x + 
           10*Power(c,3)*Power(ld,2)*x - 25*Power(ld,3)*x - 
           30*c*Power(ld,3)*x - 10*Power(c,2)*Power(ld,3)*x + 
           10*Power(ld,4)*x + 5*c*Power(ld,4)*x - Power(ld,5)*x - 
           31*Power(c,4)*Power(x,2) + 6*ld*Power(x,2) + 
           39*c*ld*Power(x,2) + 101*Power(c,2)*ld*Power(x,2) + 
           124*Power(c,3)*ld*Power(x,2) - 67*Power(ld,2)*Power(x,2) - 
           202*c*Power(ld,2)*Power(x,2) - 
           186*Power(c,2)*Power(ld,2)*Power(x,2) + 
           101*Power(ld,3)*Power(x,2) + 124*c*Power(ld,3)*Power(x,2) - 
           31*Power(ld,4)*Power(x,2) + 180*Power(c,3)*Power(x,3) - 
           30*ld*Power(x,3) - 210*c*ld*Power(x,3) - 
           540*Power(c,2)*ld*Power(x,3) + 210*Power(ld,2)*Power(x,3) + 
           540*c*Power(ld,2)*Power(x,3) - 180*Power(ld,3)*Power(x,3) - 
           390*Power(c,2)*Power(x,4) + 120*ld*Power(x,4) + 
           780*c*ld*Power(x,4) - 390*Power(ld,2)*Power(x,4) + 
           360*c*Power(x,5) - 360*ld*Power(x,5) - 120*Power(x,6);
    }
    
    if (cumulant == 7)
    {
        correction = 
          Power(c,6)*x - ld*x - 6*c*ld*x - 15*Power(c,2)*ld*x - 
            20*Power(c,3)*ld*x - 15*Power(c,4)*ld*x - 6*Power(c,5)*ld*x + 
            31*Power(ld,2)*x + 90*c*Power(ld,2)*x + 
            105*Power(c,2)*Power(ld,2)*x + 60*Power(c,3)*Power(ld,2)*x + 
            15*Power(c,4)*Power(ld,2)*x - 90*Power(ld,3)*x - 
            150*c*Power(ld,3)*x - 90*Power(c,2)*Power(ld,3)*x - 
            20*Power(c,3)*Power(ld,3)*x + 65*Power(ld,4)*x + 
            60*c*Power(ld,4)*x + 15*Power(c,2)*Power(ld,4)*x - 
            15*Power(ld,5)*x - 6*c*Power(ld,5)*x + Power(ld,6)*x - 
            63*Power(c,5)*Power(x,2) + 7*ld*Power(x,2) + 
            56*c*ld*Power(x,2) + 189*Power(c,2)*ld*Power(x,2) + 
            336*Power(c,3)*ld*Power(x,2) + 315*Power(c,4)*ld*Power(x,2) - 
            161*Power(ld,2)*Power(x,2) - 651*c*Power(ld,2)*Power(x,2) - 
            1008*Power(c,2)*Power(ld,2)*Power(x,2) - 
            630*Power(c,3)*Power(ld,2)*Power(x,2) + 
            462*Power(ld,3)*Power(x,2) + 1008*c*Power(ld,3)*Power(x,2) + 
            630*Power(c,2)*Power(ld,3)*Power(x,2) - 
            336*Power(ld,4)*Power(x,2) - 315*c*Power(ld,4)*Power(x,2) + 
            63*Power(ld,5)*Power(x,2) + 602*Power(c,4)*Power(x,3) - 
            42*ld*Power(x,3) - 378*c*ld*Power(x,3) - 
            1372*Power(c,2)*ld*Power(x,3) - 2408*Power(c,3)*ld*Power(x,3) + 
            644*Power(ld,2)*Power(x,3) + 2744*c*Power(ld,2)*Power(x,3) + 
            3612*Power(c,2)*Power(ld,2)*Power(x,3) - 
            1372*Power(ld,3)*Power(x,3) - 2408*c*Power(ld,3)*Power(x,3) + 
            602*Power(ld,4)*Power(x,3) - 2100*Power(c,3)*Power(x,4) + 
            210*ld*Power(x,4) + 1890*c*ld*Power(x,4) + 
            6300*Power(c,2)*ld*Power(x,4) - 1890*Power(ld,2)*Power(x,4) - 
            6300*c*Power(ld,2)*Power(x,4) + 2100*Power(ld,3)*Power(x,4) + 
            3360*Power(c,2)*Power(x,5) - 840*ld*Power(x,5) - 
            6720*c*ld*Power(x,5) + 3360*Power(ld,2)*Power(x,5) - 
            2520*c*Power(x,6) + 2520*ld*Power(x,6) + 720*Power(x,7);
    }
    
    if (cumulant == 8)
    {
        correction =
     Power(c,7)*x - ld*x - 7*c*ld*x - 21*Power(c,2)*ld*x - 
       35*Power(c,3)*ld*x - 35*Power(c,4)*ld*x - 
       21*Power(c,5)*ld*x - 7*Power(c,6)*ld*x + 63*Power(ld,2)*x + 
       217*c*Power(ld,2)*x + 315*Power(c,2)*Power(ld,2)*x + 
       245*Power(c,3)*Power(ld,2)*x + 105*Power(c,4)*Power(ld,2)*x + 
       21*Power(c,5)*Power(ld,2)*x - 301*Power(ld,3)*x - 
       630*c*Power(ld,3)*x - 525*Power(c,2)*Power(ld,3)*x - 
       210*Power(c,3)*Power(ld,3)*x - 35*Power(c,4)*Power(ld,3)*x + 
       350*Power(ld,4)*x + 455*c*Power(ld,4)*x + 
       210*Power(c,2)*Power(ld,4)*x + 35*Power(c,3)*Power(ld,4)*x - 
       140*Power(ld,5)*x - 105*c*Power(ld,5)*x - 
       21*Power(c,2)*Power(ld,5)*x + 21*Power(ld,6)*x + 
       7*c*Power(ld,6)*x - Power(ld,7)*x - 
       127*Power(c,6)*Power(x,2) + 8*ld*Power(x,2) + 
       76*c*ld*Power(x,2) + 316*Power(c,2)*ld*Power(x,2) + 
       734*Power(c,3)*ld*Power(x,2) + 1002*Power(c,4)*ld*Power(x,2) + 
       762*Power(c,5)*ld*Power(x,2) - 367*Power(ld,2)*Power(x,2) - 
       1826*c*Power(ld,2)*Power(x,2) - 
       3801*Power(c,2)*Power(ld,2)*Power(x,2) - 
       4008*Power(c,3)*Power(ld,2)*Power(x,2) - 
       1905*Power(c,4)*Power(ld,2)*Power(x,2) + 
       1798*Power(ld,3)*Power(x,2) + 5400*c*Power(ld,3)*Power(x,2) + 
       6012*Power(c,2)*Power(ld,3)*Power(x,2) + 
       2540*Power(c,3)*Power(ld,3)*Power(x,2) - 
       2333*Power(ld,4)*Power(x,2) - 4008*c*Power(ld,4)*Power(x,2) - 
       1905*Power(c,2)*Power(ld,4)*Power(x,2) + 
       1002*Power(ld,5)*Power(x,2) + 762*c*Power(ld,5)*Power(x,2) - 
       127*Power(ld,6)*Power(x,2) + 1932*Power(c,5)*Power(x,3) - 
       56*ld*Power(x,3) - 616*c*ld*Power(x,3) - 
       2884*Power(c,2)*ld*Power(x,3) - 7196*Power(c,3)*ld*Power(x,3) - 
       9660*Power(c,4)*ld*Power(x,3) + 1736*Power(ld,2)*Power(x,3) + 
       9856*c*Power(ld,2)*Power(x,3) + 
       21588*Power(c,2)*Power(ld,2)*Power(x,3) + 
       19320*Power(c,3)*Power(ld,2)*Power(x,3) - 
       6972*Power(ld,3)*Power(x,3) - 21588*c*Power(ld,3)*Power(x,3) - 
       19320*Power(c,2)*Power(ld,3)*Power(x,3) + 
       7196*Power(ld,4)*Power(x,3) + 9660*c*Power(ld,4)*Power(x,3) - 
       1932*Power(ld,5)*Power(x,3) - 10206*Power(c,4)*Power(x,4) + 
       336*ld*Power(x,4) + 3864*c*ld*Power(x,4) + 
       17976*Power(c,2)*ld*Power(x,4) + 40824*Power(c,3)*ld*Power(x,4) - 
       6552*Power(ld,2)*Power(x,4) - 35952*c*Power(ld,2)*Power(x,4) - 
       61236*Power(c,2)*Power(ld,2)*Power(x,4) + 
       17976*Power(ld,3)*Power(x,4) + 40824*c*Power(ld,3)*Power(x,4) - 
       10206*Power(ld,4)*Power(x,4) + 25200*Power(c,3)*Power(x,5) - 
       1680*ld*Power(x,5) - 18480*c*ld*Power(x,5) - 
       75600*Power(c,2)*ld*Power(x,5) + 18480*Power(ld,2)*Power(x,5) + 
       75600*c*Power(ld,2)*Power(x,5) - 25200*Power(ld,3)*Power(x,5) - 
       31920*Power(c,2)*Power(x,6) + 6720*ld*Power(x,6) + 
       63840*c*ld*Power(x,6) - 31920*Power(ld,2)*Power(x,6) + 
       20160*c*Power(x,7) - 20160*ld*Power(x,7) - 5040*Power(x,8);
    }
    
    if (cumulant == 9)
    {
        correction = 
     Power(c,8)*x - ld*x - 8*c*ld*x - 28*Power(c,2)*ld*x - 
       56*Power(c,3)*ld*x - 70*Power(c,4)*ld*x - 
       56*Power(c,5)*ld*x - 28*Power(c,6)*ld*x - 8*Power(c,7)*ld*x + 
       127*Power(ld,2)*x + 504*c*Power(ld,2)*x + 
       868*Power(c,2)*Power(ld,2)*x + 840*Power(c,3)*Power(ld,2)*x + 
       490*Power(c,4)*Power(ld,2)*x + 168*Power(c,5)*Power(ld,2)*x + 
       28*Power(c,6)*Power(ld,2)*x - 966*Power(ld,3)*x - 
       2408*c*Power(ld,3)*x - 2520*Power(c,2)*Power(ld,3)*x - 
       1400*Power(c,3)*Power(ld,3)*x - 420*Power(c,4)*Power(ld,3)*x - 
       56*Power(c,5)*Power(ld,3)*x + 1701*Power(ld,4)*x + 
       2800*c*Power(ld,4)*x + 1820*Power(c,2)*Power(ld,4)*x + 
       560*Power(c,3)*Power(ld,4)*x + 70*Power(c,4)*Power(ld,4)*x - 
       1050*Power(ld,5)*x - 1120*c*Power(ld,5)*x - 
       420*Power(c,2)*Power(ld,5)*x - 56*Power(c,3)*Power(ld,5)*x + 
       266*Power(ld,6)*x + 168*c*Power(ld,6)*x + 
       28*Power(c,2)*Power(ld,6)*x - 28*Power(ld,7)*x - 
       8*c*Power(ld,7)*x + Power(ld,8)*x - 
       255*Power(c,7)*Power(x,2) + 9*ld*Power(x,2) + 
       99*c*ld*Power(x,2) + 489*Power(c,2)*ld*Power(x,2) + 
       1401*Power(c,3)*ld*Power(x,2) + 2505*Power(c,4)*ld*Power(x,2) + 
       2787*Power(c,5)*ld*Power(x,2) + 1785*Power(c,6)*ld*Power(x,2) - 
       813*Power(ld,2)*Power(x,2) - 4755*c*Power(ld,2)*Power(x,2) - 
       12201*Power(c,2)*Power(ld,2)*Power(x,2) - 
       17331*Power(c,3)*Power(ld,2)*Power(x,2) - 
       13935*Power(c,4)*Power(ld,2)*Power(x,2) - 
       5355*Power(c,5)*Power(ld,2)*Power(x,2) + 
       6429*Power(ld,3)*Power(x,2) + 24078*c*Power(ld,3)*Power(x,2) + 
       36963*Power(c,2)*Power(ld,3)*Power(x,2) + 
       27870*Power(c,3)*Power(ld,3)*Power(x,2) + 
       8925*Power(c,4)*Power(ld,3)*Power(x,2) - 
       13278*Power(ld,4)*Power(x,2) - 31953*c*Power(ld,4)*Power(x,2) - 
       27870*Power(c,2)*Power(ld,4)*Power(x,2) - 
       8925*Power(c,3)*Power(ld,4)*Power(x,2) + 
       9816*Power(ld,5)*Power(x,2) + 13935*c*Power(ld,5)*Power(x,2) + 
       5355*Power(c,2)*Power(ld,5)*Power(x,2) - 
       2787*Power(ld,6)*Power(x,2) - 1785*c*Power(ld,6)*Power(x,2) + 
       255*Power(ld,7)*Power(x,2) + 6050*Power(c,6)*Power(x,3) - 
       72*ld*Power(x,3) - 936*c*ld*Power(x,3) - 
       5364*Power(c,2)*ld*Power(x,3) - 17316*Power(c,3)*ld*Power(x,3) - 
       33252*Power(c,4)*ld*Power(x,3) - 36300*Power(c,5)*ld*Power(x,3) + 
       4374*Power(ld,2)*Power(x,3) + 30420*c*Power(ld,2)*Power(x,3) + 
       88998*Power(c,2)*Power(ld,2)*Power(x,3) + 
       133008*Power(c,3)*Power(ld,2)*Power(x,3) + 
       90750*Power(c,4)*Power(ld,2)*Power(x,3) - 
       29720*Power(ld,3)*Power(x,3) - 126048*c*Power(ld,3)*Power(x,3) - 
       199512*Power(c,2)*Power(ld,3)*Power(x,3) - 
       121000*Power(c,3)*Power(ld,3)*Power(x,3) + 
       54366*Power(ld,4)*Power(x,3) + 133008*c*Power(ld,4)*Power(x,3) + 
       90750*Power(c,2)*Power(ld,4)*Power(x,3) - 
       33252*Power(ld,5)*Power(x,3) - 36300*c*Power(ld,5)*Power(x,3) + 
       6050*Power(ld,6)*Power(x,3) - 46620*Power(c,5)*Power(x,4) + 
       504*ld*Power(x,4) + 7056*c*ld*Power(x,4) + 
       42084*Power(c,2)*ld*Power(x,4) + 
       134316*Power(c,3)*ld*Power(x,4) + 
       233100*Power(c,4)*ld*Power(x,4) - 19656*Power(ld,2)*Power(x,4) - 
       143136*c*Power(ld,2)*Power(x,4) - 
       402948*Power(c,2)*Power(ld,2)*Power(x,4) - 
       466200*Power(c,3)*Power(ld,2)*Power(x,4) + 
       101052*Power(ld,3)*Power(x,4) + 402948*c*Power(ld,3)*Power(x,4) + 
       466200*Power(c,2)*Power(ld,3)*Power(x,4) - 
       134316*Power(ld,4)*Power(x,4) - 233100*c*Power(ld,4)*Power(x,4) + 
       46620*Power(ld,5)*Power(x,4) + 166824*Power(c,4)*Power(x,5) - 
       3024*ld*Power(x,5) - 42336*c*ld*Power(x,5) - 
       239904*Power(c,2)*ld*Power(x,5) - 
       667296*Power(c,3)*ld*Power(x,5) + 71568*Power(ld,2)*Power(x,5) + 
       479808*c*Power(ld,2)*Power(x,5) + 
       1000944*Power(c,2)*Power(ld,2)*Power(x,5) - 
       239904*Power(ld,3)*Power(x,5) - 667296*c*Power(ld,3)*Power(x,5) + 
       166824*Power(ld,4)*Power(x,5) - 317520*Power(c,3)*Power(x,6) + 
       15120*ld*Power(x,6) + 196560*c*ld*Power(x,6) + 
       952560*Power(c,2)*ld*Power(x,6) - 196560*Power(ld,2)*Power(x,6) - 
       952560*c*Power(ld,2)*Power(x,6) + 317520*Power(ld,3)*Power(x,6) + 
       332640*Power(c,2)*Power(x,7) - 60480*ld*Power(x,7) - 
       665280*c*ld*Power(x,7) + 332640*Power(ld,2)*Power(x,7) - 
       181440*c*Power(x,8) + 181440*ld*Power(x,8) + 40320*Power(x,9);
    }
    
    if (cumulant == 10)
    {
        correction =
     Power(c,9)*x - ld*x - 9*c*ld*x - 36*Power(c,2)*ld*x - 
       84*Power(c,3)*ld*x - 126*Power(c,4)*ld*x - 
       126*Power(c,5)*ld*x - 84*Power(c,6)*ld*x - 
       36*Power(c,7)*ld*x - 9*Power(c,8)*ld*x + 255*Power(ld,2)*x + 
       1143*c*Power(ld,2)*x + 2268*Power(c,2)*Power(ld,2)*x + 
       2604*Power(c,3)*Power(ld,2)*x + 1890*Power(c,4)*Power(ld,2)*x + 
       882*Power(c,5)*Power(ld,2)*x + 252*Power(c,6)*Power(ld,2)*x + 
       36*Power(c,7)*Power(ld,2)*x - 3025*Power(ld,3)*x - 
       8694*c*Power(ld,3)*x - 10836*Power(c,2)*Power(ld,3)*x - 
       7560*Power(c,3)*Power(ld,3)*x - 3150*Power(c,4)*Power(ld,3)*x - 
       756*Power(c,5)*Power(ld,3)*x - 84*Power(c,6)*Power(ld,3)*x + 
       7770*Power(ld,4)*x + 15309*c*Power(ld,4)*x + 
       12600*Power(c,2)*Power(ld,4)*x + 5460*Power(c,3)*Power(ld,4)*x + 
       1260*Power(c,4)*Power(ld,4)*x + 126*Power(c,5)*Power(ld,4)*x - 
       6951*Power(ld,5)*x - 9450*c*Power(ld,5)*x - 
       5040*Power(c,2)*Power(ld,5)*x - 1260*Power(c,3)*Power(ld,5)*x - 
       126*Power(c,4)*Power(ld,5)*x + 2646*Power(ld,6)*x + 
       2394*c*Power(ld,6)*x + 756*Power(c,2)*Power(ld,6)*x + 
       84*Power(c,3)*Power(ld,6)*x - 462*Power(ld,7)*x - 
       252*c*Power(ld,7)*x - 36*Power(c,2)*Power(ld,7)*x + 
       36*Power(ld,8)*x + 9*c*Power(ld,8)*x - Power(ld,9)*x - 
       511*Power(c,8)*Power(x,2) + 10*ld*Power(x,2) + 
       125*c*ld*Power(x,2) + 715*Power(c,2)*ld*Power(x,2) + 
       2435*Power(c,3)*ld*Power(x,2) + 5377*Power(c,4)*ld*Power(x,2) + 
       7853*Power(c,5)*ld*Power(x,2) + 7387*Power(c,6)*ld*Power(x,2) + 
       4088*Power(c,7)*ld*Power(x,2) - 1771*Power(ld,2)*Power(x,2) - 
       11838*c*Power(ld,2)*Power(x,2) - 
       35758*Power(c,2)*Power(ld,2)*Power(x,2) - 
       62706*Power(c,3)*Power(ld,2)*Power(x,2) - 
       68032*Power(c,4)*Power(ld,2)*Power(x,2) - 
       44322*Power(c,5)*Power(ld,2)*Power(x,2) - 
       14308*Power(c,6)*Power(ld,2)*Power(x,2) + 
       21879*Power(ld,3)*Power(x,2) + 97010*c*Power(ld,3)*Power(x,2) + 
       185967*Power(c,2)*Power(ld,3)*Power(x,2) + 
       193598*Power(c,3)*Power(ld,3)*Power(x,2) + 
       110805*Power(c,4)*Power(ld,3)*Power(x,2) + 
       28616*Power(c,5)*Power(ld,3)*Power(x,2) - 
       67671*Power(ld,4)*Power(x,2) - 205324*c*Power(ld,4)*Power(x,2) - 
       251132*Power(c,2)*Power(ld,4)*Power(x,2) - 
       147740*Power(c,3)*Power(ld,4)*Power(x,2) - 
       35770*Power(c,4)*Power(ld,4)*Power(x,2) + 
       76686*Power(ld,5)*Power(x,2) + 154333*c*Power(ld,5)*Power(x,2) + 
       110805*Power(c,2)*Power(ld,5)*Power(x,2) + 
       28616*Power(c,3)*Power(ld,5)*Power(x,2) - 
       36620*Power(ld,6)*Power(x,2) - 44322*c*Power(ld,6)*Power(x,2) - 
       14308*Power(c,2)*Power(ld,6)*Power(x,2) + 
       7387*Power(ld,7)*Power(x,2) + 4088*c*Power(ld,7)*Power(x,2) - 
       511*Power(ld,8)*Power(x,2) + 18660*Power(c,7)*Power(x,3) - 
       90*ld*Power(x,3) - 1350*c*ld*Power(x,3) - 
       9150*Power(c,2)*ld*Power(x,3) - 36210*Power(c,3)*ld*Power(x,3) - 
       90210*Power(c,4)*ld*Power(x,3) - 
       141630*Power(c,5)*ld*Power(x,3) - 
       130620*Power(c,6)*ld*Power(x,3) + 10590*Power(ld,2)*Power(x,3) + 
       86280*c*Power(ld,2)*Power(x,3) + 
       309750*Power(c,2)*Power(ld,2)*Power(x,3) + 
       619620*Power(c,3)*Power(ld,2)*Power(x,3) + 
       708150*Power(c,4)*Power(ld,2)*Power(x,3) + 
       391860*Power(c,5)*Power(ld,2)*Power(x,3) - 
       115140*Power(ld,3)*Power(x,3) - 606720*c*Power(ld,3)*Power(x,3) - 
       1317600*Power(c,2)*Power(ld,3)*Power(x,3) - 
       1416300*Power(c,3)*Power(ld,3)*Power(x,3) - 
       653100*Power(c,4)*Power(ld,3)*Power(x,3) + 
       333180*Power(ld,4)*Power(x,3) + 
       1137180*c*Power(ld,4)*Power(x,3) + 
       1416300*Power(c,2)*Power(ld,4)*Power(x,3) + 
       653100*Power(c,3)*Power(ld,4)*Power(x,3) - 
       348990*Power(ld,5)*Power(x,3) - 708150*c*Power(ld,5)*Power(x,3) - 
       391860*Power(c,2)*Power(ld,5)*Power(x,3) + 
       141630*Power(ld,6)*Power(x,3) + 130620*c*Power(ld,6)*Power(x,3) - 
       18660*Power(ld,7)*Power(x,3) - 204630*Power(c,6)*Power(x,4) + 
       720*ld*Power(x,4) + 11880*c*ld*Power(x,4) + 
       86400*Power(c,2)*ld*Power(x,4) + 
       354600*Power(c,3)*ld*Power(x,4) + 
       870120*Power(c,4)*ld*Power(x,4) + 
       1227780*Power(c,5)*ld*Power(x,4) - 54450*Power(ld,2)*Power(x,4) - 
       484380*c*Power(ld,2)*Power(x,4) - 
       1813770*Power(c,2)*Power(ld,2)*Power(x,4) - 
       3480480*Power(c,3)*Power(ld,2)*Power(x,4) - 
       3069450*Power(c,4)*Power(ld,2)*Power(x,4) + 
       470940*Power(ld,3)*Power(x,4) + 
       2563740*c*Power(ld,3)*Power(x,4) + 
       5220720*Power(c,2)*Power(ld,3)*Power(x,4) + 
       4092600*Power(c,3)*Power(ld,3)*Power(x,4) - 
       1104570*Power(ld,4)*Power(x,4) - 
       3480480*c*Power(ld,4)*Power(x,4) - 
       3069450*Power(c,2)*Power(ld,4)*Power(x,4) + 
       870120*Power(ld,5)*Power(x,4) + 
       1227780*c*Power(ld,5)*Power(x,4) - 
       204630*Power(ld,6)*Power(x,4) + 1020600*Power(c,5)*Power(x,5) - 
       5040*ld*Power(x,5) - 85680*c*ld*Power(x,5) - 
       619920*Power(c,2)*ld*Power(x,5) - 
       2404080*Power(c,3)*ld*Power(x,5) - 
       5103000*Power(c,4)*ld*Power(x,5) + 
       236880*Power(ld,2)*Power(x,5) + 
       2101680*c*Power(ld,2)*Power(x,5) + 
       7212240*Power(c,2)*Power(ld,2)*Power(x,5) + 
       10206000*Power(c,3)*Power(ld,2)*Power(x,5) - 
       1481760*Power(ld,3)*Power(x,5) - 
       7212240*c*Power(ld,3)*Power(x,5) - 
       10206000*Power(c,2)*Power(ld,3)*Power(x,5) + 
       2404080*Power(ld,4)*Power(x,5) + 
       5103000*c*Power(ld,4)*Power(x,5) - 
       1020600*Power(ld,5)*Power(x,5) - 2739240*Power(c,4)*Power(x,6) + 
       30240*ld*Power(x,6) + 498960*c*ld*Power(x,6) + 
       3331440*Power(c,2)*ld*Power(x,6) + 
       10956960*Power(c,3)*ld*Power(x,6) - 
       841680*Power(ld,2)*Power(x,6) - 
       6662880*c*Power(ld,2)*Power(x,6) - 
       16435440*Power(c,2)*Power(ld,2)*Power(x,6) + 
       3331440*Power(ld,3)*Power(x,6) + 
       10956960*c*Power(ld,3)*Power(x,6) - 
       2739240*Power(ld,4)*Power(x,6) + 4233600*Power(c,3)*Power(x,7) - 
       151200*ld*Power(x,7) - 2268000*c*ld*Power(x,7) - 
       12700800*Power(c,2)*ld*Power(x,7) + 
       2268000*Power(ld,2)*Power(x,7) + 
       12700800*c*Power(ld,2)*Power(x,7) - 
       4233600*Power(ld,3)*Power(x,7) - 3780000*Power(c,2)*Power(x,8) + 
       604800*ld*Power(x,8) + 7560000*c*ld*Power(x,8) - 
       3780000*Power(ld,2)*Power(x,8) + 1814400*c*Power(x,9) - 
       1814400*ld*Power(x,9) - 362880*Power(x,10);
    }
    
    if (cumulant == 11)
    {
        correction = 
     Power(c,10)*x - ld*x - 10*c*ld*x - 45*Power(c,2)*ld*x - 
       120*Power(c,3)*ld*x - 210*Power(c,4)*ld*x - 
       252*Power(c,5)*ld*x - 210*Power(c,6)*ld*x - 
       120*Power(c,7)*ld*x - 45*Power(c,8)*ld*x - 
       10*Power(c,9)*ld*x + 511*Power(ld,2)*x + 
       2550*c*Power(ld,2)*x + 5715*Power(c,2)*Power(ld,2)*x + 
       7560*Power(c,3)*Power(ld,2)*x + 6510*Power(c,4)*Power(ld,2)*x + 
       3780*Power(c,5)*Power(ld,2)*x + 1470*Power(c,6)*Power(ld,2)*x + 
       360*Power(c,7)*Power(ld,2)*x + 45*Power(c,8)*Power(ld,2)*x - 
       9330*Power(ld,3)*x - 30250*c*Power(ld,3)*x - 
       43470*Power(c,2)*Power(ld,3)*x - 36120*Power(c,3)*Power(ld,3)*x - 
       18900*Power(c,4)*Power(ld,3)*x - 6300*Power(c,5)*Power(ld,3)*x - 
       1260*Power(c,6)*Power(ld,3)*x - 120*Power(c,7)*Power(ld,3)*x + 
       34105*Power(ld,4)*x + 77700*c*Power(ld,4)*x + 
       76545*Power(c,2)*Power(ld,4)*x + 42000*Power(c,3)*Power(ld,4)*x + 
       13650*Power(c,4)*Power(ld,4)*x + 2520*Power(c,5)*Power(ld,4)*x + 
       210*Power(c,6)*Power(ld,4)*x - 42525*Power(ld,5)*x - 
       69510*c*Power(ld,5)*x - 47250*Power(c,2)*Power(ld,5)*x - 
       16800*Power(c,3)*Power(ld,5)*x - 3150*Power(c,4)*Power(ld,5)*x - 
       252*Power(c,5)*Power(ld,5)*x + 22827*Power(ld,6)*x + 
       26460*c*Power(ld,6)*x + 11970*Power(c,2)*Power(ld,6)*x + 
       2520*Power(c,3)*Power(ld,6)*x + 210*Power(c,4)*Power(ld,6)*x - 
       5880*Power(ld,7)*x - 4620*c*Power(ld,7)*x - 
       1260*Power(c,2)*Power(ld,7)*x - 120*Power(c,3)*Power(ld,7)*x + 
       750*Power(ld,8)*x + 360*c*Power(ld,8)*x + 
       45*Power(c,2)*Power(ld,8)*x - 45*Power(ld,9)*x - 
       10*c*Power(ld,9)*x + Power(ld,10)*x - 
       1023*Power(c,9)*Power(x,2) + 11*ld*Power(x,2) + 
       154*c*ld*Power(x,2) + 1001*Power(c,2)*ld*Power(x,2) + 
       3949*Power(c,3)*ld*Power(x,2) + 10373*Power(c,4)*ld*Power(x,2) + 
       18733*Power(c,5)*ld*Power(x,2) + 23177*Power(c,6)*ld*Power(x,2) + 
       18898*Power(c,7)*ld*Power(x,2) + 9207*Power(c,8)*ld*Power(x,2) - 
       3817*Power(ld,2)*Power(x,2) - 28611*c*Power(ld,2)*Power(x,2) - 
       98890*Power(c,2)*Power(ld,2)*Power(x,2) - 
       204402*Power(c,3)*Power(ld,2)*Power(x,2) - 
       274120*Power(c,4)*Power(ld,2)*Power(x,2) - 
       241296*Power(c,5)*Power(ld,2)*Power(x,2) - 
       132286*Power(c,6)*Power(ld,2)*Power(x,2) - 
       36828*Power(c,7)*Power(ld,2)*Power(x,2) + 
       72204*Power(ld,3)*Power(x,2) + 367158*c*Power(ld,3)*Power(x,2) + 
       834273*Power(c,2)*Power(ld,3)*Power(x,2) + 
       1085700*Power(c,3)*Power(ld,3)*Power(x,2) + 
       858825*Power(c,4)*Power(ld,3)*Power(x,2) + 
       396858*Power(c,5)*Power(ld,3)*Power(x,2) + 
       85932*Power(c,6)*Power(ld,3)*Power(x,2) - 
       322212*Power(ld,4)*Power(x,2) - 
       1165967*c*Power(ld,4)*Power(x,2) - 
       1799710*Power(c,2)*Power(ld,4)*Power(x,2) - 
       1485880*Power(c,3)*Power(ld,4)*Power(x,2) - 
       661430*Power(c,4)*Power(ld,4)*Power(x,2) - 
       128898*Power(c,5)*Power(ld,4)*Power(x,2) + 
       525723*Power(ld,5)*Power(x,2) + 
       1345135*c*Power(ld,5)*Power(x,2) + 
       1369995*Power(c,2)*Power(ld,5)*Power(x,2) + 
       661430*Power(c,3)*Power(ld,5)*Power(x,2) + 
       128898*Power(c,4)*Power(ld,5)*Power(x,2) - 
       375738*Power(ld,6)*Power(x,2) - 650232*c*Power(ld,6)*Power(x,2) - 
       396858*Power(c,2)*Power(ld,6)*Power(x,2) - 
       85932*Power(c,3)*Power(ld,6)*Power(x,2) + 
       125411*Power(ld,7)*Power(x,2) + 132286*c*Power(ld,7)*Power(x,2) + 
       36828*Power(c,2)*Power(ld,7)*Power(x,2) - 
       18898*Power(ld,8)*Power(x,2) - 9207*c*Power(ld,8)*Power(x,2) + 
       1023*Power(ld,9)*Power(x,2) + 57002*Power(c,8)*Power(x,3) - 
       110*ld*Power(x,3) - 1870*c*ld*Power(x,3) - 
       14630*Power(c,2)*ld*Power(x,3) - 68530*Power(c,3)*ld*Power(x,3) - 
       209594*Power(c,4)*ld*Power(x,3) - 
       427966*Power(c,5)*ld*Power(x,3) - 
       570284*Power(c,6)*ld*Power(x,3) - 
       456016*Power(c,7)*ld*Power(x,3) + 24992*Power(ld,2)*Power(x,3) + 
       232056*c*Power(ld,2)*Power(x,3) + 
       977306*Power(c,2)*Power(ld,2)*Power(x,3) + 
       2402532*Power(c,3)*Power(ld,2)*Power(x,3) + 
       3681854*Power(c,4)*Power(ld,2)*Power(x,3) + 
       3421704*Power(c,5)*Power(ld,2)*Power(x,3) + 
       1596056*Power(c,6)*Power(ld,2)*Power(x,3) - 
       420948*Power(ld,3)*Power(x,3) - 
       2618440*c*Power(ld,3)*Power(x,3) - 
       7074144*Power(c,2)*Power(ld,3)*Power(x,3) - 
       10447756*Power(c,3)*Power(ld,3)*Power(x,3) - 
       8554260*Power(c,4)*Power(ld,3)*Power(x,3) - 
       3192112*Power(c,5)*Power(ld,3)*Power(x,3) + 
       1813482*Power(ld,4)*Power(x,3) + 
       7779068*c*Power(ld,4)*Power(x,3) + 
       13531804*Power(c,2)*Power(ld,4)*Power(x,3) + 
       11405680*Power(c,3)*Power(ld,4)*Power(x,3) + 
       3990140*Power(c,4)*Power(ld,4)*Power(x,3) - 
       2897862*Power(ld,5)*Power(x,3) - 
       8307926*c*Power(ld,5)*Power(x,3) - 
       8554260*Power(c,2)*Power(ld,5)*Power(x,3) - 
       3192112*Power(c,3)*Power(ld,5)*Power(x,3) + 
       1969990*Power(ld,6)*Power(x,3) + 
       3421704*c*Power(ld,6)*Power(x,3) + 
       1596056*Power(c,2)*Power(ld,6)*Power(x,3) - 
       570284*Power(ld,7)*Power(x,3) - 456016*c*Power(ld,7)*Power(x,3) + 
       57002*Power(ld,8)*Power(x,3) - 874500*Power(c,7)*Power(x,4) + 
       990*ld*Power(x,4) + 18810*c*ld*Power(x,4) + 
       161370*Power(c,2)*ld*Power(x,4) + 
       808830*Power(c,3)*ld*Power(x,4) + 
       2559150*Power(c,4)*ld*Power(x,4) + 
       5133150*Power(c,5)*ld*Power(x,4) + 
       6121500*Power(c,6)*ld*Power(x,4) - 
       143550*Power(ld,2)*Power(x,4) - 
       1492920*c*Power(ld,2)*Power(x,4) - 
       6839910*Power(c,2)*Power(ld,2)*Power(x,4) - 
       17493300*Power(c,3)*Power(ld,2)*Power(x,4) - 
       25665750*Power(c,4)*Power(ld,2)*Power(x,4) - 
       18364500*Power(c,5)*Power(ld,2)*Power(x,4) + 
       1976040*Power(ld,3)*Power(x,4) + 
       13332660*c*Power(ld,3)*Power(x,4) + 
       37125000*Power(c,2)*Power(ld,3)*Power(x,4) + 
       51331500*Power(c,3)*Power(ld,3)*Power(x,4) + 
       30607500*Power(c,4)*Power(ld,3)*Power(x,4) - 
       7301580*Power(ld,4)*Power(x,4) - 
       32006700*c*Power(ld,4)*Power(x,4) - 
       51331500*Power(c,2)*Power(ld,4)*Power(x,4) - 
       30607500*Power(c,3)*Power(ld,4)*Power(x,4) + 
       9815850*Power(ld,5)*Power(x,4) + 
       25665750*c*Power(ld,5)*Power(x,4) + 
       18364500*Power(c,2)*Power(ld,5)*Power(x,4) - 
       5133150*Power(ld,6)*Power(x,4) - 
       6121500*c*Power(ld,6)*Power(x,4) + 
       874500*Power(ld,7)*Power(x,4) + 5921520*Power(c,6)*Power(x,5) - 
       7920*ld*Power(x,5) - 158400*c*ld*Power(x,5) - 
       1393920*Power(c,2)*ld*Power(x,5) - 
       6922080*Power(c,3)*ld*Power(x,5) - 
       20603880*Power(c,4)*ld*Power(x,5) - 
       35529120*Power(c,5)*ld*Power(x,5) + 
       716760*Power(ld,2)*Power(x,5) + 
       7753680*c*Power(ld,2)*Power(x,5) + 
       35287560*Power(c,2)*Power(ld,2)*Power(x,5) + 
       82415520*Power(c,3)*Power(ld,2)*Power(x,5) + 
       88822800*Power(c,4)*Power(ld,2)*Power(x,5) - 
       7513440*Power(ld,3)*Power(x,5) - 
       49808880*c*Power(ld,3)*Power(x,5) - 
       123623280*Power(c,2)*Power(ld,3)*Power(x,5) - 
       118430400*Power(c,3)*Power(ld,3)*Power(x,5) + 
       21443400*Power(ld,4)*Power(x,5) + 
       82415520*c*Power(ld,4)*Power(x,5) + 
       88822800*Power(c,2)*Power(ld,4)*Power(x,5) - 
       20603880*Power(ld,5)*Power(x,5) - 
       35529120*c*Power(ld,5)*Power(x,5) + 
       5921520*Power(ld,6)*Power(x,5) - 21538440*Power(c,5)*Power(x,6) + 
       55440*ld*Power(x,6) + 1108800*c*ld*Power(x,6) + 
       9424800*Power(c,2)*ld*Power(x,6) + 
       42966000*Power(c,3)*ld*Power(x,6) + 
       107692200*Power(c,4)*ld*Power(x,6) - 
       3049200*Power(ld,2)*Power(x,6) - 
       31878000*c*Power(ld,2)*Power(x,6) - 
       128898000*Power(c,2)*Power(ld,2)*Power(x,6) - 
       215384400*Power(c,3)*Power(ld,2)*Power(x,6) + 
       22453200*Power(ld,3)*Power(x,6) + 
       128898000*c*Power(ld,3)*Power(x,6) + 
       215384400*Power(c,2)*Power(ld,3)*Power(x,6) - 
       42966000*Power(ld,4)*Power(x,6) - 
       107692200*c*Power(ld,4)*Power(x,6) + 
       21538440*Power(ld,5)*Power(x,6) + 
       46070640*Power(c,4)*Power(x,7) - 332640*ld*Power(x,7) - 
       6320160*c*ld*Power(x,7) - 48565440*Power(c,2)*ld*Power(x,7) - 
       184282560*Power(c,3)*ld*Power(x,7) + 
       10644480*Power(ld,2)*Power(x,7) + 
       97130880*c*Power(ld,2)*Power(x,7) + 
       276423840*Power(c,2)*Power(ld,2)*Power(x,7) - 
       48565440*Power(ld,3)*Power(x,7) - 
       184282560*c*Power(ld,3)*Power(x,7) + 
       46070640*Power(ld,4)*Power(x,7) - 
       59875200*Power(c,3)*Power(x,8) + 1663200*ld*Power(x,8) + 
       28274400*c*ld*Power(x,8) + 179625600*Power(c,2)*ld*Power(x,8) - 
       28274400*Power(ld,2)*Power(x,8) - 
       179625600*c*Power(ld,2)*Power(x,8) + 
       59875200*Power(ld,3)*Power(x,8) + 
       46569600*Power(c,2)*Power(x,9) - 6652800*ld*Power(x,9) - 
       93139200*c*ld*Power(x,9) + 46569600*Power(ld,2)*Power(x,9) - 
       19958400*c*Power(x,10) + 19958400*ld*Power(x,10) + 
       3628800*Power(x,11);
    }
    *correction_ptr = correction;
}
