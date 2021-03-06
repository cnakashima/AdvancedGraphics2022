#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"
#include <math.h>

#define SCREEN_SIZE 800
#define HITHER .001
double ORIGIN[3] =  {0,0,0};
char nameA[100] ;
int d[2];
int idA ;


double zbuff[SCREEN_SIZE][SCREEN_SIZE];


Bool inbounds(double result[3],  int map[2]){
  // closer than hither the out
  // if tan(angle) > tan(half angle) y/z > tan(half angle) pr x/z > tan(half angle)
  // or if x/z < -tan(half angle) y/z < -tan(half angle)
  // make sure map points are between 0 and 800
  
  if (result[2] < HITHER){
   return False;
   }
  if (result[0]/result[2] > SCREEN_SIZE/2) {
  return False;
  }
  if (result[1]/result[2] > SCREEN_SIZE/2) {
  return False;
  }
  if (result[1]/result[2] < -SCREEN_SIZE/2) { 
  return False;
  }
  if (result[0]/result[2] < -SCREEN_SIZE/2) {
  return False;
  }
  if (map[0] < 0 && map[0] >= 800) {
  return False;
  }
  if (map[1] < 0 && map[1] >= 800) {
  return False;
  }
  return True;
  
}

// To support the light model (taken from view_test05MbS.c by Jeff) :
double light_in_eye_space[3] = {0.0, 0.0, 0.0};
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
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




  double intensity ;
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




// set color base on texture
int Texture(double color[3], int texture_code, double u, double v, double lowu, double upu, double lowv, double upv, double texture[3]){

// (((int)((v+lowv)/box_height))%2) (((int)((u+lowu)/box_width))%2)

	if (texture_code == 1){
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
    double wratio = d[0]/width;
    double hratio = d[1]/height;
    int nl ;
    int tlist[100] ;
    double plist[100] ;
    double move[4][4];
  double movei[4][4];

    nl = 0 ;
    tlist[nl] = TX; plist[nl] = -lowu ; nl++ ;
    tlist[nl] = TY; plist[nl] = -lowv ; nl++ ;
    tlist[nl] = SX ; plist[nl] = wratio ; nl++ ;
    tlist[nl] = SY; plist[nl] = fabs(hratio) ; nl++ ;
    

    M3d_make_movement_sequence_matrix (move,movei,  nl,tlist,plist) ; 

    double coords[3];
    coords[0] = u;
    coords[1] = v;
    coords[2] = 0;

    M3d_mat_mult_pt(coords, move, coords);

    // need to map u,v coordinates (double) to x,y coordinates (int) within dimensions of file in d
    int e = get_xwd_map_color(idA, (int)coords[0],(int)coords[1],texture) ; 
    
    // returns -1 on error, 1 if ok
        // if (e == -1) { continue ; }
  }
	
	return 0;

}



// take in parametric fuctions for x, y , and z, compute each point in 3d space and then map to 2d plane and plot
int graph_function_3d(double (*fx)(double u, double v), double (*fy)(double u, double v), double (*fz)(double u, double v),
 	    double lowboundu,
	    double upboundu,
      double lowboundv,
      double upboundv,
	    double move[4][4], double color[3], double step, int texturekey){
  // printf("check\n");
  G_rgb(color[0],color[1],color[2]);
  double result[3];// store the original 3d point
  double temp[3]; // mess with 3d point to get it 2d
  int map_point[2]; // map onto 2d with integers for zbuff
  double helperq[3]; // points to help create the normal vector
  double helperr[3];
  double a[3];// vectors to get normal vector
  double b[3]; 
  double n[3];// normal vector for light model
  double h = .001;
  double rgb[3];
  
  for (double u = lowboundu; u < upboundu; u += step){
    for (double v = lowboundv;  v < upboundv; v+=step){
    // printf("%f, %f\n", u, v);
    result[0] = fx(u, v); // generate coordinates for point in object space
    result[1] = fy(u, v);
    result[2] = fz(u, v);
    
    helperq[0] = fx(u + h, v); // generate coordinates for helper points in object space
    helperq[1] = fy(u + h, v);
    helperq[2] = fz(u + h, v);
    
    helperr[0] = fx(u, v+h);
    helperr[1] = fy(u, v+h);
    helperr[2] = fz(u, v+h);
    
    //set texture in object space
    double texture[3];
    Texture(color, texturekey, u, v, lowboundu, upboundu, lowboundv, upboundv, texture);

    // printf("1\n");
    M3d_mat_mult_pt(result, move, result); // move all points to eye space
    M3d_mat_mult_pt(helperq, move, helperq); 
    M3d_mat_mult_pt(helperr, move, helperr); 


    // from here map to 2d
    // printf("2\n");
    temp[0] = (result[0]/result[2]) * SCREEN_SIZE/2;
    temp[1] = (result[1]/result[2]) * SCREEN_SIZE/2; 
    temp[2] = 0;
    // currently centered at 0,0 so move to center at the center of the screen
    double translator[4][4];
    M3d_make_translation(translator, SCREEN_SIZE/2, SCREEN_SIZE/2, 0);
    M3d_mat_mult_pt(temp, translator, temp);
    map_point[0] = temp[0];
    map_point[1] = temp[1];

    // printf("%i, %i\n", map_point[0], map_point[1]);
    // M3d_print_mat(move);
    // printf("%f, %f\n", result[0],result[1]);

    // printf("%lf\n", zbuff[map_point[0]][map_point[1]]);

    // check if point is closer than a current point on grid, if so plot it
    if (zbuff[map_point[0]][map_point[1]] > result[2] && inbounds(result, map_point)){
    
      a[0] = helperq[0] - result[0]; a[1] = helperq[1] - result[1]; a[2] = helperq[2] - result[2];
      b[0] = helperr[0] - result[0]; b[1] = helperr[1] - result[1]; b[2] = helperr[2] - result[2];
    
    M3d_x_product (n, a, b);
    
    Light_Model (texture, ORIGIN, result,n,rgb);// make color dependent on u,v -> make textures
    
      G_rgb(rgb[0], rgb[1], rgb[2]);
      
      
      G_point(temp[0],temp[1]);
      zbuff[map_point[0]][map_point[1]] = result[2];
    }
  }
  }        
}

double cylinderx(double u, double v){
  return cos(u);
}
double cylindery(double u, double v){
  return v;
}
double cylinderz(double u, double v){
  return sin(u);
}

double spherex(double u, double v){
  return sin(u)*cos(v);
}
double spherey(double u, double v){
  return sin(u)*sin(v);
}
double spherez(double u, double v){
  return cos(u);
}

int init_zbuff(){
  for (int i = 0; i <SCREEN_SIZE; i++){
    for (int j = 0; j < SCREEN_SIZE; j++){
      zbuff[i][j] = 1e50;
    }
  }
}

int main(void){


  

  int e ;
  double rgb[3] ;

  printf("enter name of xwd file\n") ;
  scanf("%s",nameA) ;
  idA = init_xwd_map_from_file (nameA) ;// returns -1 on error, 1 if ok
  if (idA == -1) { printf("failure\n") ;  exit(0) ; }
  e = get_xwd_map_dimensions(idA, d) ;
  if (e == -1) { printf("failure\n") ;  exit(0) ; }



  double V[4][4];
  double Vi[4][4];
  double additional_movement[4][4];
  double additional_movementi[4][4];
  double move[4][4];
  double movei[4][4];
  char q = ' ';

  double eye[3], coi[3], up[3] ;

  int nl ;
  int tlist[100] ;
  double plist[100] ;
  
  int fnum ;
  double t ;

  fnum = 0 ;

  G_init_graphics (SCREEN_SIZE, SCREEN_SIZE) ;


  double color[3];  


  while (q != 'q') {
    init_zbuff();
    t = 0.01*fnum ;

    eye[0] = 15*cos(2*M_PI*t) ; 
    eye[1] =  6*t ; 
    eye[2] =  7*sin(2*M_PI*t) ; 

    // printf("t = %lf   eye = %lf %lf %lf\n",t, eye[0],eye[1],eye[2]) ;

    coi[0] =  1.0 ;
    coi[1] =  2.0 ; 
    coi[2] =  0.5 ;

    up[0]  = eye[0] ; 
    up[1]  = eye[1] + 1 ;
    up[2]  = eye[2] ; 

    M3d_view (V, Vi,  eye,coi,up) ;

  G_rgb(0,0,0) ; 
  G_clear() ;





  //sphere

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 4 ; nl++ ;
  tlist[nl] = SY; plist[nl] = 4 ; nl++ ;
  tlist[nl] = SZ; plist[nl] = 4 ; nl++ ;

  M3d_make_movement_sequence_matrix (additional_movement,additional_movementi,  nl,tlist,plist) ; 

  M3d_mat_mult(move, V, additional_movement);


  color[0] = 1.0; color[1] = 1.0; color[2] = 0.0;
  graph_function_3d(spherex, spherey, spherez, 0, 2*M_PI, -.5*M_PI, .5*M_PI, move, color, .001*M_PI, 2 );


  q = G_wait_key() ;


    fnum++ ;
  } // end while (1)

}
