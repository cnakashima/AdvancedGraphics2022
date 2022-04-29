#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include <math.h>

double x[50][1000], y[50][1000];
double z[1000] = {0};
int numobj;

int graph_function(double (*fx)(double x), double (*fy)(double y),
 	    double upbound,
	    double lowbound,
	    double move[4][4], double color[3], int onum, double step){
  
  G_rgb(color[0],color[1],color[2]);
  int index = 0;
  for (double t = lowbound; t < upbound; t += step){
    x[onum][index] = fx(t);
    y[onum][index] = fy(t);
    M3d_mat_mult_points(x[onum], y[onum], z, move, x[onum], y[onum], z, index+1);
    G_point(x[onum][index],y[onum][index]);
    printf("%f, %f\n", x[onum][index], y[onum][index]);
    index ++;
  }        

}


int main(void){


int onum = 0;

// must do this before you do 'almost' any other graphical tasks 
double swidth = 800 ;  double sheight = 800 ;
G_init_graphics (swidth,sheight) ;  // interactive graphics

   
// clear the screen in a given color
G_rgb (0.0, 0.0, 0.0) ; // dark gray
G_clear () ;

int Tn = 0;
int Ttypelist[100]; double Tvlist[100];
double v[4][4], vi[4][4];
double color[3];

// circle

  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;



M3d_make_movement_sequence_matrix(v, vi, Tn, Ttypelist, Tvlist);
color[0]=1, color[1] = 0, color [2] = 0;
graph_function(cos, sin, 1.5*M_PI, .25*M_PI, v, color, onum, .005*M_PI);

onum++;

G_wait_key();

//screen = 800x800

}
//also 3 cricles SX20, SY20, TX100, TY30
// x = 1-sin(t), y = 1-cos(t) (0,2pi)
//made from a unit circle rolling follow bottom point

// sphere parametric equation 
// x = cosv* cosu y = sinv z = cosv*sinu u =[0,2pi] v = [-pi/2,pi/2]
// better: x=sqrt(1-v^2)cos(u) y = v z = sqrt(1-v^2) sinu u = [0,2pi]v=[-1,1]