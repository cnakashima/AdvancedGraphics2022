#include "FPToolkit.c"
#include <unistd.h>


int main (int argc, char **argv)
{
  if (argc != 8){
    printf("pgm usage: pgm_name screen_width screen_height prefix_name back_and_forth(0 if no else 1) start_index end_index wait_time\n");
    exit(0);
  }

  int width,height,back_and_forth, start,end;
  long wait_time;
  char prefix[100] ;
  width = atoi(argv[1]);
  height = atoi(argv[2]);
  back_and_forth = atoi(argv[4]);
  start = atoi(argv[5]);
  end = atoi(argv[6]);
  wait_time = atol(argv[7]);
  sprintf(prefix, "%s", argv[3]);


  G_init_graphics(width,height) ;
 

  char filename[200];
  char q;
  int i = start;
  int direction = 1;
  while (q != 'q'){
    sprintf(filename, "%s%04d.xwd", prefix, i);

    G_get_image_from_file(filename,0,0) ;

    G_display_image();

    q = G_no_wait_key() ;

    usleep(wait_time);
    i += direction;
    if (i == end + 1){
      if (back_and_forth == 1){
        direction *= -1;
        i = end - 1;
      }
      else{
        i = start;
      }
    }
    if (i == start - 1){
      direction *= -1;
      i = start + 1;
    }
  }


}