#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"
#include <math.h>

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define HITHER .001
#define HALFANGLE 45
double ORIGIN[3] =  {0,0,0};

double obmat[100000][4][4] ;
double obinv[100000][4][4] ;
double color[100000][3] ;
double reflection[100000]; // [0.0,1.0] 1.0 complete reflects colors, 0.0 reflects no other colors
double transparent[100000]; // [0.0,1.0] 1.0 complete refracts colors, 0.0 refracts no other colors
int    obtype[100000]; // type of shape: 0 = sphere, 1 = plane, 2 = hyperboloid
int    num_objects ;
int texture[100000]; // tells what texture to put on object
double triangle_coords[100000][3][3];

int d1[2];
int idA ;
int d2[2];
int idB ;
int d3[2];
int idC ;

double vm[4][4], vi[4][4];




  int destinationid;


////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////


int ray (double* Rsource, double * Rtip, double* argb, int depth, int flag, int from);

// To support the light model (taken from view_test05MbS.c by Jeff) :
double light_in_eye_space[3] = {0.0, 0, 0.0};
double AMBIENT      = 0.3 ;
double MAX_DIFFUSE  = 0.7 ;
double SPECPOW      = 50 ;



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye(0,0,0))
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; // unit vector with direction of n
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ; // Eu, eye vector
  E[0] = s[0] - p[0] ;
  E[1] = s[1] - p[1] ;
  E[2] = s[2] - p[2] ;
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ; // Light vector
  L[0] = light_in_eye_space[0] - p[0] ;
  L[1] = light_in_eye_space[1] - p[1] ;
  L[2] = light_in_eye_space[2] - p[2] ;
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;

  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below

  // shadow

  //vector from close to light





   double lsource[3];
   double ltip[3];
   lsource[0] = p[0] + .001*L[0];
   lsource[1] = p[1] + .001*L[1];
   lsource[2] = p[2] + .001*L[2];

   ltip[0] = p[0] + L[0];
   ltip[1] = p[1] + L[1];
   ltip[2] = p[2] + L[2];

   double dummy[3];



  double intensity ;

  int object_found = ray(lsource, ltip, dummy, 0, 0, -1);
  if (object_found && transparent[object_found - 1] == 0){
  // test for shadow issues: correct light vector make sure it points toward light
      intensity = AMBIENT;
      goto LLL;
      //argb[0] = 0; argb[1] = 0; argb[2] = 0;
      //return 0;
  }

  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ;
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ;
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) /
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}


///////////////////////////////////////////////////////////////////////////////////

// set color base on texture
int Texture(double color[3], int texture_code, double u, double v, double lowu, double upu, double lowv, double upv, double texture[3]){

// (((int)((v+lowv)/box_height))%2) (((int)((u+lowu)/box_width))%2)

	if (texture_code == 3){
		double box_height = (upv - lowv)/10.0;
		double box_width = (upu - lowu)/10.0;
		if ((((int)((v-lowv)/box_height))%2) ==  (((int)((u-lowu)/box_width))%2)){
			texture[0] = 1.0;
			texture[1] = 1.0;
			texture[2] = 1.0;
		}
		else{
		texture[0] = color[0];
		texture[1]=color[1];
		texture[2]= color[2];
		}
	}
  if (texture_code == 2){
    double height = upv-lowv;
    double width = upu-lowu;
    // printf("%d %d\n",d[0],d[1]) ;
    double wratio = d1[0]/width;
    double hratio = d1[1]/height;


    double x = wratio * (u-lowu);
    double y = hratio * (v-lowv);

    // need to map u,v coordinates (double) to x,y coordinates (int) within dimensions of file in d
    int e = get_xwd_map_color(idA, x, y,texture) ;

    // returns -1 on error, 1 if ok
        // if (e == -1) { continue ; }
  }
  if (texture_code == 4){
    double height = upv-lowv;
    double width = upu-lowu;
    // printf("%d %d\n",d[0],d[1]) ;
    double wratio = d2[0]/width;
    double hratio = d2[1]/height;


    double x = wratio * (u-lowu);
    double y = hratio * (v-lowv);

    // need to map u,v coordinates (double) to x,y coordinates (int) within dimensions of file in d
    int e = get_xwd_map_color(idB, x, y,texture) ;

    // returns -1 on error, 1 if ok
        // if (e == -1) { continue ; }
  }if (texture_code == 1){
    double height = upv-lowv;
    double width = upu-lowu;
    // printf("%d %d\n",d[0],d[1]) ;
    double wratio = d3[0]/width;
    double hratio = d3[1]/height;


    double x = wratio * (u-lowu);
    double y = hratio * (v-lowv);

    // need to map u,v coordinates (double) to x,y coordinates (int) within dimensions of file in d
    int e = get_xwd_map_color(idC, x, y,texture) ;

    // returns -1 on error, 1 if ok
        // if (e == -1) { continue ; }
  }


  if (texture_code == 0){
     texture[0] = color[0];
     texture[1] = color[1];
     texture[2] = color[2];
  }
	return 0;

}

////////////////////////////////////////////////////////////////////////

int normal (int n, double* pt, double* npt){
	// for a circle
	double dx;
	double dy;
	double dz;
	if (obtype[n] == 0){
		dx = 2*pt[0];
		dy = 2*pt[1];
		dz = 2*pt[2];
	}
	else if (obtype[n] == 1){
		dx = 0;
		dy = 1;
		dz = 0;
	}
	else if (obtype[n] == 2){
		dx = 2*pt[0];
		dy = -2*pt[1];
		dz = 2*pt[2];
	}
	else if (obtype[n] == 3){
		dx = 2*pt[0];
		dy = 0;
		dz = 2*pt[2];
	}
  else if (obtype[n] == 4){
    double eyeA[3];
    double eyeB[3];
    double eyeC[3];
    M3d_mat_mult_pt(eyeA, vm, triangle_coords[n][0]);
    M3d_mat_mult_pt(eyeB, vm, triangle_coords[n][1]);
    M3d_mat_mult_pt(eyeC, vm, triangle_coords[n][2]);

    double AB[3] = {eyeA[0] - eyeB[0], eyeA[1] - eyeB[1], eyeA[2] - eyeB[2]};
    double AC[3] = {eyeA[0] - eyeC[0], eyeA[1] - eyeC[1], eyeA[2] - eyeC[2]};
    M3d_x_product(npt, AB, AC);
    return 0;
  }
	npt[0] = obinv[n][0][0] * dx + obinv[n][1][0] * dy + obinv[n][2][0] * dz;
	npt[1] = obinv[n][0][1] * dx + obinv[n][1][1] * dy + obinv[n][2][1] * dz;
	npt[2] = obinv[n][0][2] * dx + obinv[n][1][2] * dy + obinv[n][2][2] * dz;



}


int sphere (double A, double B, double C, double t[2]){
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

int plane (double* source, double dx, double dy, double dz, double* t){
	if (dy == 0){
		return 0;
	}
	double k = -source[1] / dy;
	double x = source[0] + k * dx;
	double z = source[2] + k * dz;
	if ((x <= 1 && x >= -1) && (z <= 1 && z >= -1)){
		t[0] = k;
		return 1;
	}
	return 0;
}

int hyperboloid (double A, double B, double C, double t[2]){
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

int cylinder (double A, double B, double C, double t[2]){
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

int triangle(double S[3], double E[3], int trinum, double* t){
  // solve for t based on system of equations S-A = u(B-A) + v(C-A) + t(S-E)
  double num[3][3];
  double denom[3][3];

  double Bx_minus_Ax = triangle_coords[trinum][1][0] - triangle_coords[trinum][0][0];//B-A
  double By_minus_Ay = triangle_coords[trinum][1][1] - triangle_coords[trinum][0][1];
  double Bz_minus_Az = triangle_coords[trinum][1][2] - triangle_coords[trinum][0][2];

  double Cx_minus_Ax = triangle_coords[trinum][2][0] - triangle_coords[trinum][0][0];//C-A
  double Cy_minus_Ay = triangle_coords[trinum][2][1] - triangle_coords[trinum][0][1];
  double Cz_minus_Az = triangle_coords[trinum][2][2] - triangle_coords[trinum][0][2];

  double Sx_minus_Ex = S[0]-E[0];
  double Sy_minus_Ey = S[1]-E[1];
  double Sz_minus_Ez = S[2]-E[2];

  double Sx_minus_Ax = S[0] - triangle_coords[trinum][0][0];//S-A
  double Sy_minus_Ay = S[1] - triangle_coords[trinum][0][1];
  double Sz_minus_Az = S[2] - triangle_coords[trinum][0][2];

  num[0][0] = Bx_minus_Ax; num[0][1] = Cx_minus_Ax; num[0][2] = Sx_minus_Ax;
  num[1][0] = By_minus_Ay; num[1][1] = Cy_minus_Ay; num[1][2] = Sy_minus_Ay;
  num[2][0] = Bz_minus_Az; num[2][1] = Cz_minus_Az; num[2][2] = Sz_minus_Az;

  denom[0][0] = Bx_minus_Ax; denom[0][1] = Cx_minus_Ax; denom[0][2] = Sx_minus_Ex;
  denom[1][0] = By_minus_Ay; denom[1][1] = Cy_minus_Ay; denom[1][2] = Sy_minus_Ey;
  denom[2][0] = Bz_minus_Az; denom[2][1] = Cz_minus_Az; denom[2][2] = Sz_minus_Ez;

  if (det(denom) == 0){
    return 0;
  }

  double parameter_t = det(num) / det(denom);

  num[0][0] = Sx_minus_Ax; num[0][1] = Cx_minus_Ax; num[0][2] = Sx_minus_Ex;
  num[1][0] = Sy_minus_Ay; num[1][1] = Cy_minus_Ay; num[1][2] = Sy_minus_Ey;
  num[2][0] = Sz_minus_Az; num[2][1] = Cz_minus_Az; num[2][2] = Sz_minus_Ez;

  double u = det(num)/det(denom);

  num[0][0] = Bx_minus_Ax; num[0][1] = Sx_minus_Ax; num[0][2] = Sx_minus_Ex;
  num[1][0] = By_minus_Ay; num[1][1] = Sy_minus_Ay; num[1][2] = Sy_minus_Ey;
  num[2][0] = Bz_minus_Az; num[2][1] = Sz_minus_Az; num[2][2] = Sz_minus_Ez;

  double v = det(num)/det(denom);


  if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && u+v <= 1 && t >= 0){
    t[0] = parameter_t;
    return 1;
  }
  return 0;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

double ripple (double x, double z, int frame){
  double step = .2;
  double progress = -2 + (frame-1)*step;
  if (frame < 1){
    return 0;
  }
  if (sqrt(x*x + z*z) > (frame-1)*step/27){
    return 0;
  }
  double dim = frame/13.5;
  if (dim < 1){
    dim = 1;
  }
  double amplitude = 1.5/((dim)*(10+80*sqrt(x*x + z*z)));
  double wave = cos(30*sqrt(x*x + z*z) - progress);
  double y = amplitude * wave;
  return y;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int reflect(int depth, int pt_found, double* close, double* nu, double* vu, double* argb, double* reflect_rgb){
	double source[3], tip[3], lrgb[3];
	double h = nu[0] * vu[0] + nu[1] * vu[1] + nu[2] * vu[2];
	double ru[3];
	ru[0] = 2*h*nu[0] - vu[0];
	ru[1] = 2*h*nu[1] - vu[1];
	ru[2] = 2*h*nu[2] - vu[2];

  // printf("%i: %f,%f,%f\n", obtype[pt_found], nu[0],nu[1],nu[2]);

	// set new source aprintfnd tip
	source[0] = close[0] + .001 * ru[0];
	source[1] = close[1] + .001 * ru[1];
	source[2] = close[2] + .001 * ru[2];

	tip[0] = close[0] + ru[0];
	tip[1] = close[1] + ru[1];
	tip[2] = close[2] + ru[2];
  int ob = ray(source, tip, lrgb, depth - 1, 1, pt_found);
	if (!ob){
		lrgb[0] = 0;
		lrgb[1] = 0;
		lrgb[2] = 0;
	}
  // reflect_rgb[0] = (1 - reflection[pt_found]) * argb[0] + reflection[pt_found] * lrgb[0];
	// reflect_rgb[1] = (1 - reflection[pt_found]) * argb[1] + reflection[pt_found] * lrgb[1];
	// reflect_rgb[2] = (1 - reflection[pt_found]) * argb[2] + reflection[pt_found] * lrgb[2];


  reflect_rgb[0] = lrgb[0];
	reflect_rgb[1] = lrgb[1];
	reflect_rgb[2] = lrgb[2];

}

int refract (int from, int depth, int pt_found, double* close, double* nu, double* vu, double* argb, double* refract_rgb){
// vu needs to be negated
	double source[3], tip[3], lrgb[3];
	vu[0] *= -1; vu[1] *= -1; vu[2] *= -1;
	double mu[3];
	mu[0] = -nu[0]; mu[1] = -nu[1]; mu[2] = -nu[2];
	double k = 1/1.33; // refraction index of glass for
  int sign = 1;
	if (pt_found == from){
		k = 1/k;
    sign = -1;
	}
	double cosl = mu[0]*vu[0] + mu[1]*vu[1] + mu[2]*vu[2];
	// note 1 - (k*k) * (1 - (cosl*cosl)) will be negative if total internal reflection for dense to less dense material have to work on this...
	double ru[3];

	if (1 - (k*k) * (1 - (cosl*cosl)) < 0) {
		vu[0] *= -1; vu[1] *= -1; vu[2] *= -1;
		reflect(depth - 1, pt_found, close, nu, vu, argb, refract_rgb);
	}
	else{
		ru[0] = k * vu[0] + sign * (sqrt(1 - (k*k) * (1 - (cosl*cosl))) - k*cosl)*mu[0];
		ru[1] = k * vu[1] + sign * (sqrt(1 - (k*k) * (1 - (cosl*cosl))) - k*cosl)*mu[1];
		ru[2] = k * vu[2] + sign * (sqrt(1 - (k*k) * (1 - (cosl*cosl))) - k*cosl)*mu[2];
		// set new source and tip
		source[0] = close[0] + .001 * ru[0];
		source[1] = close[1] + .001 * ru[1];
		source[2] = close[2] + .001 * ru[2];

		tip[0] = close[0] + ru[0];
		tip[1] = close[1] + ru[1];
		tip[2] = close[2] + ru[2];

		if (!ray(source, tip, lrgb, depth, 1, pt_found)){
			lrgb[0] = 0;
			lrgb[1] = 0;
			lrgb[2] = 0;
		}

    // refract_rgb[0] = (1 - transparent[pt_found]) * argb[0] + transparent[pt_found] * lrgb[0];
		// refract_rgb[1] = (1 - transparent[pt_found]) * argb[1] + transparent[pt_found] * lrgb[1];
		// refract_rgb[2] = (1 - transparent[pt_found]) * argb[2] + transparent[pt_found] * lrgb[2];

    refract_rgb[0] = lrgb[0];
		refract_rgb[1] =  lrgb[1];
		refract_rgb[2] = lrgb[2];
	}
}


int ray (double* Rsource, double * Rtip, double* argb, int depth, int flag, int from){
	double close[4];
	double oclose[4];
	double closet = 1000000000000;
	double Otip[4];
	double tip[4];
	double Osource[4];
	double source[4];
	double t[2];
	int pt_found = -1;
	double npt[3];
	double irgb[3];
	double lrgb[3];

	source[0] = Rsource[0];
	source[1] = Rsource[1];
	source[2] = Rsource[2];
	source[3] = Rsource[3];

	tip[0] = Rtip[0];
	tip[1] = Rtip[1];
	tip[2] = Rtip[2];
	tip[3] = Rtip[3];


	for (int i = 0; i < num_objects; i++){

		M3d_mat_mult_pt(Otip, obinv[i], tip);
		M3d_mat_mult_pt(Osource, obinv[i], source);

		double dx = Otip[0] - Osource[0];
		double dy = Otip[1] - Osource[1];
		double dz = Otip[2] - Osource[2];
		double A, B, C;
		int n;
		if (obtype[i] == 0){//sphere
			A = dx*dx + dy*dy + dz*dz;
			B = 2* Osource[0] * dx + 2* Osource[1] * dy + 2*Osource[2] * dz;
			C = Osource[0] * Osource[0] + Osource[1]*Osource[1] + Osource[2]*Osource[2] -1;
			n = sphere(A, B, C, t);
		}
		else if (obtype[i] == 1){//plane
			n = plane(Osource, dx, dy, dz, t);
		}
		else if (obtype[i] == 2){//hyperboloid
			A = dx*dx - dy*dy + dz*dz;
			B = 2* Osource[0] * dx - 2* Osource[1] * dy + 2*Osource[2] * dz;
			C = Osource[0] * Osource[0] - Osource[1]*Osource[1] + Osource[2]*Osource[2] -1;
			n = hyperboloid(A, B, C, t);
		}
		else if (obtype[i] == 3){//cylinder
			A = dx*dx + dz*dz;
			B = 2* Osource[0] * dx + 2*Osource[2] * dz;
			C = Osource[0] * Osource[0]  + Osource[2]*Osource[2] -1;
			n = cylinder(A, B, C, t);
		}
    else if (obtype[i] == 4){// 2d triangle
      n = triangle(Osource, Otip, i, t);

    }
		for (int j = 0; j < n; j++){
			if (t[j] > 0 && t[j] < closet){
				closet = t[j];
				double x = Osource[0] + closet * dx;
				double y = Osource[1] + closet * dy;
				double z = Osource[2] + closet * dz;
				if (obtype[i] == 2 && (x > cosh(1) || x < -cosh(1) || y < sinh(-1) || y > sinh(1))){ // skip if it is a hyperbola and we are out of range
					continue;
				}
				if (obtype[i] == 3 && (y > 1 || y < -1 )){ // skip if it is a cylinder and we are out of range
					continue;
				}
				oclose[0] = x;
				oclose[1] = y;
				oclose[2] = z ;
				oclose[3] = 0;
				M3d_mat_mult_pt(close, obmat[i], oclose);
				M3d_mat_mult_pt(source, obmat[i], Osource);
				pt_found = i;
			}
		}
		}
		if (pt_found > -1){
			if (flag == 0){
			   double dl[3];
			   dl [0] = source[0] - light_in_eye_space[0];
			   dl [1] = source[1] - light_in_eye_space[1];
			   dl [2] = source[2] - light_in_eye_space[2];
			   double pl[3];
			   pl [0] = source[0] - close[0];
			   pl [1] = source[1] - close[1];
			   pl [2] = source[2] - close[2];
			   double distl = sqrt(dl[0]*dl[0] + dl[1]*dl[1] + dl[2]*dl[2]);
			   double distp = sqrt(pl[0]*pl[0] + pl[1]*pl[1] + pl[2]*pl[2]);
			   if (distl <= distp){ return 0;}
			   return pt_found + 1;
			}

			// convert x, y, z to u, v
			// for a sphere
			if (obtype[pt_found] == 0){
			   double u = atan2(oclose[2],oclose[0]);
			   double v = acos(oclose[1]);

			   irgb[0] = color[pt_found][0]; irgb[1] = color[pt_found][1]; irgb[2] = color[pt_found][2];
			   Texture(color[pt_found], texture[pt_found], u, v, -M_PI, M_PI, M_PI, 0, irgb);
			}
      else if (obtype[pt_found] == 1){
          double u = oclose[0];
          double v = oclose[2];

           irgb[0] = color[pt_found][0]; irgb[1] = color[pt_found][1]; irgb[2] = color[pt_found][2];
           Texture(color[pt_found], texture[pt_found], u, v, -1, 1, -1, 1, irgb);
      }
			else{
			   irgb[0] = color[pt_found][0]; irgb[1] = color[pt_found][1]; irgb[2] = color[pt_found][2];
			}

			normal(pt_found, oclose, npt);


			// Get the light model color to combine with reflection if applicable
			Light_Model(irgb, source, close, npt, argb);

			// calculat vu and nu for reflection and refraction
			double mag = sqrt(((npt[0])*(npt[0])) + ((npt[1])*(npt[1])) + ((npt[2])*(npt[2])));
				double nu[3];
				nu[0] = (npt[0])/mag;
				nu[1] = (npt[1])/mag;
				nu[2] = (npt[2])/mag;

				// vu = unit form of incoming ray
				double vu[3];
				vu[0] = (source[0] - close[0]) / sqrt(((close[0]-source[0])*(close[0]-source[0])) + ((close[1]-source[1])*(close[1]-source[1]))+ ((close[2]-source[2])*(close[2]-source[2])));
				vu[1] = (source[1] - close[1]) / sqrt(((close[0]-source[0])*(close[0]-source[0])) + ((close[1]-source[1])*(close[1]-source[1]))+ ((close[2]-source[2])*(close[2]-source[2])));
				vu[2] = (source[2] - close[2]) / sqrt(((close[0]-source[0])*(close[0]-source[0])) + ((close[1]-source[1])*(close[1]-source[1]))+ ((close[2]-source[2])*(close[2]-source[2])));
				// if normal vector and income vector are on opposite side negate normal vector
				if (nu[0]*vu[0] + nu[1]*vu[1] + nu[2]*vu[2] < 0){
					nu[0] *= -1;
					nu[1] *= -1;
					nu[2] *= -1;
				}

        double reflect_rgb[3] = {0, 0, 0};
        double refract_rgb[3] = {0, 0, 0};

			// if this object is reflective and depth greater than 0 reflect recursively to get actual color
      if (reflection[pt_found] > 0.0 && depth > 0){
				reflect(depth, pt_found, close, nu, vu, argb, reflect_rgb);
			}

			if (transparent[pt_found] > 0){
				refract(from, depth, pt_found, close, nu, vu, argb, refract_rgb);
			}


      double trgb[3];
      trgb[0] = refract_rgb[0];
      trgb[1] = refract_rgb[1];
      trgb[2] = refract_rgb[2];

      argb[0] = (1 - (reflection[pt_found] + transparent[pt_found]))*argb[0] + reflection[pt_found] * reflect_rgb[0] + transparent[pt_found]*refract_rgb[0];
      argb[1] = (1 - (reflection[pt_found] + transparent[pt_found]))*argb[1] + reflection[pt_found] * reflect_rgb[1] + transparent[pt_found]*refract_rgb[1];
      argb[2] = (1 - (reflection[pt_found] + transparent[pt_found]))*argb[2] + reflection[pt_found] * reflect_rgb[2] + transparent[pt_found]*refract_rgb[2];


			return pt_found + 1;


			// todo: other shapes and shadows

			//closet = 1000000000000;

			//pt_found = -1;
		}
		else{ return 0;}

}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int test01()
{
  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  double Rsource[4];
  double Rtip[4];
  double argb[3] ;
  int depth = 6;
  double H = tan(HALFANGLE * M_PI/180);
  double eye[3], coi[3], up[3] ;
  double theta;
  int fnum = 0;
  double t = .01;

  int e;
  char nameA[] = "starry_night.xwd";
  idA = init_xwd_map_from_file (nameA) ;// returns -1 on error, 1 if ok
  if (idA == -1) { printf("failure\n") ;  exit(0) ; }
  e = get_xwd_map_dimensions(idA, d1) ;
  if (e == -1) { printf("failure\n") ;  exit(0) ; }
  char nameB[] = "marbletexture.xwd";
  idB = init_xwd_map_from_file (nameB) ;// returns -1 on error, 1 if ok
  if (idB == -1) { printf("failure\n") ;  exit(0) ; }
  e = get_xwd_map_dimensions(idB, d2) ;
  if (e == -1) { printf("failure\n") ;  exit(0) ; }
  char nameC[] = "mosaic_skin.xwd";
  idC = init_xwd_map_from_file (nameC) ;// returns -1 on error, 1 if ok
  if (idC == -1) { printf("failure\n") ;  exit(0) ; }
  e = get_xwd_map_dimensions(idC, d3) ;
  if (e == -1) { printf("failure\n") ;  exit(0) ; }

    destinationid = create_new_xwd_map (SCREEN_WIDTH,SCREEN_HEIGHT) ;// returns -1 on error, 1 if ok

    //////////////////////////////////////////////////////////////////////
    M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
    //////////////////////////////////////////////////////////////////////

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // floor of pool
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.2 ;
    color[num_objects][1] = 0.2 ;
    color[num_objects][2] = 0.2 ;
    texture[num_objects] = 1;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -3   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // front of pool
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RX ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  1   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -6   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // back of pool
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RX ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -1   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -6   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // right of pool
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  1   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -6   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // left of pool
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  -1   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -6   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this


    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // far floor
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  3   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2   ; Tn++ ;


    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // close floor
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -3   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // right floor
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = SX ; Tvlist[Tn] =  2    ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  3    ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // left floor
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = SX ; Tvlist[Tn] =  2    ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  -3   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // far wall
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RX ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  5   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  1   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // close wall
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RX ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -5   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  1   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // right wall
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  1   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    // left wall
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 4;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  4   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  -6   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  1   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // ceiling
    obtype[num_objects] = 1;
    color[num_objects][0] = 0.8 ;
    color[num_objects][1] = 0.8 ;
    color[num_objects][2] = 0.8;
    texture[num_objects] = 2;
    reflection[num_objects] = 0.0;
    transparent[num_objects] = 0.0;

    Tn = 0 ;

    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  4   ; Tn++ ;

    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;
    num_objects++ ; // don't forget to do this
    //
    // // test water
    // obtype[num_objects] = 1;
    // color[num_objects][0] = 0.0 ;
    // color[num_objects][1] = 0.0 ;
    // color[num_objects][2] = 1.0;
    // texture[num_objects] = 0;
    // reflection[num_objects] = .1;
    // transparent[num_objects] = 0.9;
    //
    // Tn = 0 ;
    //
    // // Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2.25   ; Tn++ ;
    // Ttypelist[Tn] = TX ; Tvlist[Tn] =  -1  ; Tn++ ;
    // Ttypelist[Tn] = TY ; Tvlist[Tn] =  -1.9   ; Tn++ ;
    //
    // M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    // M3d_mat_mult(obmat[num_objects], vm, m) ;
    // M3d_mat_mult(obinv[num_objects], mi, vi) ;
    //
    // num_objects++ ; // don't forget to do this


    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    int triangle_ind = num_objects;

    double step = .02;
    for (double i = -1; i < 1; i += step){
      for (double j = -1+step; j <= 1; j+=step){
        // for every i. j pair make two triangles to comprise a square
        obtype[num_objects] = 4;
        color[num_objects][0] = 0.0 ;
        color[num_objects][1] = 0.0 ;
        color[num_objects][2] = 1.0;
        texture[num_objects] = 0;
        reflection[num_objects] = 0.2;
        transparent[num_objects] = 0.8;
        triangle_coords[num_objects][0][0] = i;triangle_coords[num_objects][0][1] = ripple(i, j, 0);triangle_coords[num_objects][0][2] = j;
        triangle_coords[num_objects][1][0] = i;triangle_coords[num_objects][1][1] = ripple(i, j-step, 0);triangle_coords[num_objects][1][2] = j-step;
        triangle_coords[num_objects][2][0] = i+step;triangle_coords[num_objects][2][1] = ripple(i+step, j-step,0);triangle_coords[num_objects][2][2] = j-step;

        Tn = 0 ;

        // Ttypelist[Tn] = RX ; Tvlist[Tn] =  20   ; Tn++ ;
        Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2.25   ; Tn++ ;
        // Ttypelist[Tn] = TY ; Tvlist[Tn] =  -1.9   ; Tn++ ;
        //
        // Ttypelist[Tn] = TX ; Tvlist[Tn] = 1   ; Tn++ ;


        M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
        M3d_mat_mult(obmat[num_objects], vm, m) ;
        M3d_mat_mult(obinv[num_objects], mi, vi) ;

        num_objects++ ; // don't forget to do this

        obtype[num_objects] = 4;
        color[num_objects][0] = 0.0 ;
        color[num_objects][1] = 0.0 ;
        color[num_objects][2] = 1.0;
        texture[num_objects] = 0;
        reflection[num_objects] = 0.2;
        transparent[num_objects] = 0.8;
        triangle_coords[num_objects][0][0] = i;triangle_coords[num_objects][0][1] = ripple(i, j,0);triangle_coords[num_objects][0][2] = j;
        triangle_coords[num_objects][1][0] = i+step;triangle_coords[num_objects][1][1] = ripple(i+step, j,0);triangle_coords[num_objects][1][2] = j;
        triangle_coords[num_objects][2][0] = i+step;triangle_coords[num_objects][2][1] = ripple(i+step, j-step,0);triangle_coords[num_objects][2][2] = j-step;

        Tn = 0 ;
        // Ttypelist[Tn] = RX ; Tvlist[Tn] =  20   ; Tn++ ;
        Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2.25   ; Tn++ ;
        // Ttypelist[Tn] = TY ; Tvlist[Tn] =  -1.9   ; Tn++ ;
        //
        // Ttypelist[Tn] = TX ; Tvlist[Tn] =  1   ; Tn++ ;


        M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
        M3d_mat_mult(obmat[num_objects], vm, m) ;
        M3d_mat_mult(obinv[num_objects], mi, vi) ;

        num_objects++ ; // don't forget to do this
      }
    }


int key;





  char name[100] ;

      theta = t*fnum;

      eye[0] = 0;
      eye[1] =  0;
      eye[2] =  -2.1;

      coi[0] =  0.0 ;
      coi[1] =  -1.3 ;
      coi[2] =  0.0 ;

      up[0]  = eye[0] ;
      up[1]  = eye[1] + 1 ;
      up[2]  = eye[2] ;

      M3d_view (vm, vi,  eye,coi,up) ;
      for (int i = 0; i < num_objects; i++){
      M3d_mat_mult(obmat[i], vm, obmat[i]) ;
      M3d_mat_mult(obinv[i], obinv[i], vi) ;
      }



      Rsource[0] =  ORIGIN[0];  Rsource[1] =  ORIGIN[1] ;  Rsource[2] = ORIGIN[2] ;


// to generate many frames
  for (int frame = 0; frame <= 46; frame ++){

    for (int i = triangle_ind; i < num_objects; i++){
      triangle_coords[i][0][1] = ripple(triangle_coords[i][0][0], triangle_coords[i][0][2],frame);
      triangle_coords[i][1][1] = ripple(triangle_coords[i][1][0], triangle_coords[i][1][2],frame);
      triangle_coords[i][2][1] = ripple(triangle_coords[i][2][0], triangle_coords[i][2][2],frame);
    }





      for (int ypix = 0 ; ypix <= SCREEN_HEIGHT ; ypix++) {
         for (int xpix = 0; xpix <= SCREEN_WIDTH ; xpix ++){

             Rtip[0]    = (H/(SCREEN_WIDTH/2)) * (xpix - SCREEN_WIDTH/2);  Rtip[1]    = (H/(SCREEN_HEIGHT/2)) * (ypix - SCREEN_HEIGHT/2) ;  Rtip[2]   =   1;

             if (ray (Rsource, Rtip, argb, depth, 1, -1)){
               int e = set_xwd_map_color(destinationid, xpix,ypix, argb[0],argb[1], argb[2]) ;
             }
             else{
               int e = set_xwd_map_color(destinationid, xpix,ypix, 0,0, 0) ;
             }
         }
      }

  sprintf(name, "wave/wave%04d.xwd", frame);
  xwd_map_to_named_xwd_file(destinationid, name) ;
    //
    // G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
    // while (G_wait_key() != 'q') ;
  }
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  test01() ;
}
