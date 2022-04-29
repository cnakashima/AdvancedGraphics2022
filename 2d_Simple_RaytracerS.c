#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
int    num_objects ;



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////





void Draw_ellipsoid (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = i*2*M_PI/n ;
    xyz[0] = cos(t) ;
    xyz[1] = sin(t) ;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }

}




void Draw_the_scene()
{
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    Draw_ellipsoid(onum) ;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int quad (double A, double B, double C, double t[2]){
	double g = B*B - 4*A*C;
	if (g < 0){
		return 0;
	}
	else if (g == 0){
		t[0] = -B / (2* A);
		return 1;
	}
	else{
		t[0] = (-B + sqrt(g))/ (2* A);
		t[1] = (-B - sqrt(g))/ (2* A);
		return 2;
	}
}


int normal (int n, double* pt, double* npt){
	// for a circle
	double dx = 2*pt[0];
	double dy = 2*pt[1];
	//	printf("object number = %d\n",n) ;
	//	double v[4][4] ;
	//      M3d_mat_mult(v, obinv[n], obmat[n]) ;
	//      M3d_print_mat(obinv[n]) ;
	
	npt[0] = obinv[n][0][0] * dx + obinv[n][1][0] * dy;
	npt[1] = obinv[n][0][1] * dx + obinv[n][1][1] * dy;
	
}


int ray (double* Rsource, double * Rtip, double* argb){
	double close[3];
	double oclose[3];
	double closet = 1000000000000;
	double Otip[3];
	double Osource[3];
	double t[2];
	int pt_found = -1;
	double npt[2];
	for (int i = 0; i < num_objects; i++){
		M3d_mat_mult_pt(Otip, obinv[i], Rtip);
		M3d_mat_mult_pt(Osource, obinv[i], Rsource);
		double dx = Otip[0] - Osource[0];
		double dy = Otip[1] - Osource[1];
		double A = dx*dx + dy*dy;
		double B = 2* Osource[0] * dx + 2* Osource[1] * dy;
		double C = Osource[0] * Osource[0] + Osource[1]*Osource[1] -1;
		int n = quad(A, B, C, t);
		for (int j = 0; j < n; j++){	
			if (t[j] > 0 && t[j] < closet){
				closet = t[j];
				oclose[0] = Osource[0] + closet * dx;
				oclose[1] = Osource[1] + closet * dy;
				oclose[2] = 0 ;
				M3d_mat_mult_pt(close, obmat[i], oclose);
				argb[0] = color[i][0]; argb[1] = color[i][1];argb[2] = color[i][2];
				pt_found = i;
			}
		}
	}
	if (pt_found > -1){
		G_rgb(argb[0], argb[1], argb[2]);
		G_line(Rsource[0], Rsource[1], close[0], close[1]);
		normal(pt_found, oclose, npt);
		double mag = sqrt(((npt[0])*(npt[0])) + ((npt[1])*(npt[1]))); 
		double xcomp = (npt[0])/mag;
		double ycomp = (npt[1])/mag;
		G_rgb(1,1,1);
		//		G_line((npt[0])+close[0], (npt[1])+close[1], close[0], close[1]);
		G_line(50*xcomp+close[0], 50*ycomp+close[1], close[0], close[1]);		
		
	}
}


int test01()
{
  double vm[4][4], vi[4][4];
  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  double Rsource[3];
  double Rtip[3];
  double argb[3] ;

    //////////////////////////////////////////////////////////////////////
    M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
    //////////////////////////////////////////////////////////////////////

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.0 ;
    color[num_objects][1] = 0.8 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   25   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  200   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  180   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   40   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  400   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  550   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.3 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   75   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   35   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  150   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  360   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  500   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  130   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   30   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -15   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  700   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////

    

    G_rgb(0,0,0) ;
    G_clear() ;

    Draw_the_scene() ;
    
    Rsource[0] =  20 ;  Rsource[1] =  400 ;  Rsource[2] = 0 ;    
    G_rgb(1,0,1) ; G_fill_circle(Rsource[0], Rsource[1], 3) ;
    G_rgb(1,0,1) ; G_line(100,200,  100,600) ;
    
    G_wait_key() ;
    
    double ytip ;
    for (ytip = 200 ; ytip <= 600 ; ytip++) {
      Rtip[0]    = 100 ;  Rtip[1]    = ytip ;  Rtip[2]   = 0  ;    

      G_rgb(1,1,0) ; G_line(Rsource[0],Rsource[1],  Rtip[0],Rtip[1]) ;
      ray (Rsource, Rtip, argb) ; 

      Draw_the_scene() ;
      G_wait_key() ;
    }

    G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
    while (G_wait_key() != 'q') ;
    G_save_image_to_file("2d_Simple_Raytracer.xwd") ;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(800,800);
  test01() ;
}

