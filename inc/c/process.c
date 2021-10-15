#include <stdio.h>

void PrintProcBar(int cnt, int total)
{

  int i, j;

  int percent=cnt*100/total;
  int barnum=cnt*50/total;

  char ch[4] = {'|', '/', '-','\\'};
  j = cnt%4;

  printf("\033[0G  [ ");

  for (i = 0; i < barnum; i++)
  {
    printf("#");
  }

  if (cnt != total)
    printf("\033[55G ] %d%%  %c", percent, ch[j]);
  else
    printf("\033[55G ] %d%%  %s\n", percent, "Done!");
  fflush(stdout);
}
