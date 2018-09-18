#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "cifparse.c"

void cif2pdb(char *file, char *fpdb);

int main(int argc, char **argv)
{
    char inpfile[512], outfile[512],command[512];
    
    strcpy(outfile, "CIF2PDB.pdb");
    if (argc==3)strcpy(outfile, argv[2]);
    sprintf(command,"rm -f %s", outfile);
    system(command);
    cif2pdb(argv[1], outfile);
    printf("The output pdb =%s\n",outfile);
}

