// pointsuite 0.5.7 
// author Cathy Lawson
// RCSB-PDB
// JUNE 2007

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXMAT 500 //needed for mat44lib.c
#define MAXCHAR 80 //needed for readcif.c
#define MAXCHAIN 80
#include "mat44lib.c"
#include "readcif.c"
#define PR8 0.00000001
#define PR5 0.00001

int main(int argc, char * argv[]){
  int maxmat;
  char idc[MAXMAT][MAXCHAR];
  double mmat[MAXMAT][4][4];
  int i=0;
  int ij, j, k, l;
  int ierr;
  char cifstring[80], spacegroup[80], tmpstring[80];
  int plist[MAXMAT];
  int srlist[96],stlist[96];
  int snum, pnum, p2xnum, pcount, ncsnum, chainum;
  FILE * file, * outxframefile, *outciffile, *outncsfile, *outmatfile;
  FILE * sgciffile;
  char  buffer[200], tmpstr[80];
  char * t;
  char  pdbid[5], mattypestring[30];   
  char matframestring[30];
  double transtest;
  double p2i[4][4], i2p[4][4], p2x[4][4], x2p[4][4];
  double i2x[4][4],  x2i[4][4], i2xt[4][4], x2it[4][4];
  double p2xall[8][4][4];
  double pmat[MAXMAT][4][4]; 
  double ptemp[4][4], ptemp2[4][4]; 
  double psmat[96][MAXMAT][4][4];
  double smfixed[96][4][4];  //symmetry matrices read in
  double smat[96][4][4];    //manipulated symmetry matrices 
  double M_PRECISION, T_PRECISION, diff;
  double unitcell[6];
  double orthmat[4][4];
  double fracmat[4][4];
  double det;
  double ar, arot, arise;
  double copyrot, copyrise;
  int divisor, numbio, xsym;
  double nmat[4][4][4];
  int exp[4];
  int ident, stype, mtype, zsym, xznum, rewind;
  stype = 'U';
  zsym = 1;
 
  //precision to accept for matrix and translation comparison errors:
  //really loose helps for a few remediation entries, doesn't appear to hurt others
    // M_PRECISION = 0.0029;
    //T_PRECISION = 0.03;
    //REALLY loose: 
      M_PRECISION = 0.03;
      T_PRECISION = 0.08;

   printf ("\n\nPROGRAM findxsym\n\n");

  if(argc != 3) {
    printf("USAGE: findxsym file.cif\n");
    return 0;
  }

 
  //defaults
      unitcell[0] = 1.; unitcell[1] = 1.; unitcell[2] = 1.;
      unitcell[3] = 90.; unitcell[4] = 90.; unitcell[5] = 90.;
      smfixed[0][0][0]=1.; smfixed[0][0][1]=0.; smfixed[0][0][2]=0.; smfixed[0][0][3]=0.;
      smfixed[0][1][0]=0.; smfixed[0][1][1]=1.; smfixed[0][1][2]=0.; smfixed[0][1][3]=0.;
      smfixed[0][2][0]=0.; smfixed[0][2][1]=0.; smfixed[0][2][2]=1.; smfixed[0][2][3]=0.;
      smfixed[0][3][0]=0.; smfixed[0][3][1]=0.; smfixed[0][3][2]=0.; smfixed[0][3][3]=1.;
      snum=1;

// read cell, asym ids, symmetry info from cif
 ierr=0;
 ierr=ierr+ readOneCif(argv[1],"_cell.entry_id", pdbid);
 ierr=ierr+ readOneCif(argv[1],"_cell.length_a", cifstring); unitcell[0] = atof(cifstring);
 ierr=ierr+ readOneCif(argv[1],"_cell.length_b", cifstring); unitcell[1] = atof(cifstring);
 ierr=ierr+ readOneCif(argv[1],"_cell.length_c", cifstring); unitcell[2] = atof(cifstring);
 ierr=ierr+ readOneCif(argv[1],"_cell.angle_alpha", cifstring); unitcell[3] = atof(cifstring);
 ierr=ierr+ readOneCif(argv[1],"_cell.angle_beta", cifstring); unitcell[4] = atof(cifstring);
 ierr=ierr+ readOneCif(argv[1],"_cell.angle_gamma", cifstring); unitcell[5] = atof(cifstring);
 ierr=ierr+ readOneCif(argv[1],"_symmetry.space_group_name_H-M", spacegroup); 
 if (!ierr) {
    for (i=0;i<3;i++) {if (unitcell[i]==0 || unitcell[i]<0)unitcell[i]=1.; }
    for (i=3;i<6;i++) {if (unitcell[i]==0 || unitcell[i]<0)unitcell[i]=90.; }
    printf ("pdb id: %s\n", pdbid);
    printf ("a:      %f\n", unitcell[0]);
    printf ("b:      %f\n", unitcell[1]);
    printf ("c:      %f\n", unitcell[2]);
    printf ("alpha:  %f\n", unitcell[3]);
    printf ("beta:   %f\n", unitcell[4]);
    printf ("gamma:  %f\n", unitcell[5]);
    printf ("s.g.:   %s\n", spacegroup);   } 
    else {printf ("error: problem with  cell/s.g. info in file1\n"); return 1;}

/*
  char chain[MAXCHAIN][MAXCHAR];
  char entity[MAXCHAIN][MAXCHAR];
  int maxch;

ierr = readAsymCif(argv[1], chain, entity, &maxch);
if (ierr==1) {printf("\nError: Problem reading cif \n"); return 1;}
chainum = maxch;
printf ("number of chains read: %i\n", maxch);
  for (i=0;i<maxch;i++) printf("chain asym id: %s   entity id:  %s\n", chain[i], entity[i]);

  printf("full chain list:   ");
  for (i=0;i<maxch-1;i++) printf("%s,", chain[i]);
  printf("%s\n", chain[maxch-1]);
*/
  //calculate orthmat and fracmat
  getOrthMat(unitcell,orthmat); det=invMat(orthmat,fracmat); 

  int smtry[96][16];
  double numerator;
  char * tmppath;
    tmppath = getenv("PTSUITE");
    if (tmppath !=NULL){
    strcpy(tmpstring,tmppath);
    strcat(tmpstring,"/data/space_group.cif");
    printf("location of space group data: %s\n", tmpstring);
    sgciffile = fopen(tmpstring, "r");
  i=0;
  while (fgets(buffer, 200, sgciffile)) {
  if(strncmp(buffer, spacegroup, strlen(spacegroup))!=0) continue;
   printf ("%s", (buffer + strlen(spacegroup)));
  sscanf((buffer + strlen(spacegroup)), "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
    &smtry[i][0], &smtry[i][1], &smtry[i][2], 
    &smtry[i][3], &smtry[i][4], &smtry[i][5],
    &smtry[i][6], &smtry[i][7], &smtry[i][8], 
    &smtry[i][9], &smtry[i][10], &smtry[i][11],
    &smtry[i][12], &smtry[i][13],&smtry[i][14], &smtry[i][15]);

  //transfer to 4x4:
   smfixed[i][0][0] = smtry[i][1];
   smfixed[i][0][1] = smtry[i][2];
   smfixed[i][0][2] = smtry[i][3];
   smfixed[i][1][0] = smtry[i][4];
   smfixed[i][1][1] = smtry[i][5];
   smfixed[i][1][2] = smtry[i][6];
   smfixed[i][2][0] = smtry[i][7];
   smfixed[i][2][1] = smtry[i][8];
   smfixed[i][2][2] = smtry[i][9];
   smfixed[i][3][0] = 0;
   smfixed[i][3][1] = 0;
   smfixed[i][3][2] = 0;
   smfixed[i][3][3] = 1;
   numerator = smtry[i][10]; divisor = smtry[i][11];
   if (smtry[i][11]==0) {smfixed[i][0][3] = numerator; }
   else {smfixed[i][0][3]=numerator/divisor;}
   numerator = smtry[i][12]; divisor = smtry[i][13];
   if (smtry[i][13]==0) {smfixed[i][1][3] = numerator; }
   else {smfixed[i][1][3]=numerator/divisor;}
   numerator = smtry[i][14]; divisor = smtry[i][15];
   if (smtry[i][15]==0) {smfixed[i][2][3] = numerator; }
   else {smfixed[i][2][3]=numerator/divisor;}

 //calc SMTRY equiv
 multMat(smfixed[i],fracmat,ptemp);
 multMat(orthmat,ptemp,smfixed[i]);

   i++;
}
  fclose (sgciffile);
  snum = i;
} //end of getting info from PTSUITE/data/space_group.cif
 
printf ("number of symops read for s.g. %s:    %i\n", spacegroup, snum);
for (i=0;i<snum;i++) {
printf("symop %i\n", i+1);
for (j=0;j<4;j++) {
printf("%f    %f     %f     %f\n", smfixed[i][j][0], smfixed[i][j][1], smfixed[i][j][2], smfixed[i][j][3]);  }}

  if (snum == 1 ) { printf ("warning: only one symmetry matrix\n"); }

  //copy to smat matrix set
     for (i=0; i<snum; i++) {
          for (j=0;j<4;j++) for (k=0;k<4;k++) { smat[i][j][k] = smfixed[i][j][k];}
          }     

 /*   
  //read matrix and symmetry info from "findframe.cif"
  //default symmetry class is point
  mtype = 'P';
  stype = ' ';
  ierr=readOneCif(argv[2],"_pdbx_point_symmetry.Schoenflies_symbol",cifstring);
  if (!ierr) stype = cifstring[0];

  if (!(stype == 'C' || stype =='D' || stype =='T' || stype =='I' || stype =='O')) {
     //check for helical symmetry
     ierr=readOneCif(argv[2],"_pdbx_helical_symmetry.rotation_per_n_subunits", cifstring);
     if (!ierr) {stype = 'H';}
     else {return 1;}}

  if (stype == 'H') {
     arot = atof (cifstring);
     ierr=readOneCif(argv[2],"_pdbx_helical_symmetry.rise_per_n_subunits", cifstring);
     arise = atof (cifstring);
     ierr=ierr+readOneCif(argv[2],"_pdbx_helical_symmetry.n_subunits_divisor", cifstring);
     divisor = atof (cifstring);
     ierr=ierr+readOneCif(argv[2],"_pdbx_helical_symmetry.number_of_operations", cifstring);
     numbio = atoi (cifstring);
     ierr=ierr+readOneCif(argv[2],"_pdbx_helical_symmetry.circular_symmetry", cifstring);
     zsym = atoi (cifstring);
     ierr=ierr+readOneCif(argv[2],"_pdbx_helical_symmetry.dyad_axis", cifstring);
     if (ierr > 0) {printf("incomplete helical parameter definition\n"); return 1;}
     if (strncmp(cifstring,"yes",3)==0) xsym = 2;
     if (strncmp(cifstring,"no",2)==0) xsym = 1;
     if (numbio > MAXMAT) numbio=MAXMAT;
     if (divisor == 0) divisor=1;
     copyrot = arot;
     copyrise = arise;
     arot = arot / divisor;
     arise = arise / divisor;
     if (xsym < 1 || xsym > 2) xsym = 1;
     if (zsym < 1) zsym = 1;
     xznum = abs(zsym*xsym);    //how many mats per x/z symmetry
     pnum = (numbio/xznum)*xznum;  //pnum must be multiple of xznum
     mtype = stype;
     printf ("\n");
     if (arot*arise < 0)printf("%i Matrices requested for a left-handed helix\n",numbio);
     if (arot*arise == 0)printf("%i Matrices requested for an untwisted fiber or disc\n",numbio);
     if (arot*arise > 0)printf("%i Matrices requested for a right-handed helix\n",numbio);
     printf("   %i matrices will be written\n", pnum);
     rewind =   pnum/(2*xznum); //helix will be wound back by rewind arot/arise units
     ident = rewind*xznum+1;
     printf("   %i =identity matrix index\n", ident);
     printf("   %10.3f degree rotation per subunit around z\n", arot);
     printf("   %10.3f Angstrom rise per subunit along z\n", arise);
     printf("   %i-fold about the x axis\n", xsym);
     printf("   %i-fold about the z axis\n", zsym);
     }
  else if (stype == 'C' || stype == 'D'){ 
       ierr=readOneCif(argv[2],"_pdbx_point_symmetry.circular_symmetry",cifstring);
           zsym = atoi (cifstring); if (zsym < 1) zsym = 1;
           printf(" symmetry info (from  %s): %s %i\n", argv[2], &stype, zsym); }
     else  printf(" symmetry info (from  %s): %s\n", argv[2], &stype);

//default for p2i
  p2i[0][0]=1.0; p2i[0][1]=0.0; p2i[0][2]=0.0; p2i[0][3]=0.0;
  p2i[1][0]=0.0; p2i[1][1]=1.0; p2i[1][2]=0.0; p2i[1][3]=0.0;
  p2i[2][0]=0.0; p2i[2][1]=0.0; p2i[2][2]=1.0; p2i[2][3]=0.0;
  p2i[3][0]=0.0; p2i[3][1]=0.0; p2i[3][2]=0.0; p2i[3][3]=1.0;

  ierr = readAssemblyCif(argv[2], idc, mmat,&maxmat);
  if (!ierr) {for (i=0;i<maxmat;i++)
    {if ((!strncmp("P",idc[i],MAXCHAR)) ||(!strncmp("H",idc[i],MAXCHAR)) ) break;}
     if (i<maxmat) { for (j=0;j<4;j++) for (k=0;k<4;k++){ p2i[j][k]=mmat[i][j][k]; }}}


    printf (" transform to point or helical frame matrix: \n");
    for (i=0;i<3;i++) { 
    printf ("      %10.6f%10.6f%10.6f%15.5f\n",p2i[i][0], p2i[i][1], p2i[i][2], p2i[i][3]);
    }

  //calculate the inverse (i2p) matrix
    det=invMat (p2i, i2p);

  //get pdb2xtal matrices
  p2xnum=0;
  if (!ierr) {for (ij=0;ij<maxmat;ij++) {
      buffer[0]='X'; buffer[1]='\0'; itoa(ij,tmpstr); strcat(buffer,tmpstr);
      for (i=0;i<maxmat;i++) { if ((!strncmp(buffer,idc[i],MAXCHAR))  ) break;}
      if (i<maxmat) {for (j=0;j<4;j++) for (k=0;k<4;k++){ p2xall[p2xnum][j][k]=mmat[i][j][k]; }
            p2xnum++;}   } }

//for now cannot do ncs for helical structures 
 if (mtype == 'H') p2xnum = 0;
 printf ("%i transform to crystal frame matrices (from %s)\n", p2xnum, argv[2]);
    for (j=0;j<p2xnum;j++){
    for (i=0;i<3;i++) { 
    printf ("      %10.6f%10.6f%10.6f%15.5f\n",p2xall[j][i][0], 
        p2xall[j][i][1], p2xall[j][i][2], p2xall[j][i][3]);
    } printf ("\n");}

//END READING INPUT FROM FILES GIVEN IN COMMAND LINE

//BEGIN CIF OUT SECTION 1
  printf ("\n\n writing assembly cif records to *assembly.cif*\n");
   outciffile = fopen("assembly.cif","w");

if (mtype == 'P'){
fprintf(outciffile,"_pdbx_point_symmetry.entry_id               %s\n",pdbid);
fprintf(outciffile,"_pdbx_point_symmetry.Schoenflies_symbol        %s \n", &stype);
if (stype == 'C' || stype == 'D') {
fprintf(outciffile,"_pdbx_point_symmetry.circular_symmetry           %i\n\n",zsym);}}

if (mtype=='H'){
   fprintf(outciffile,"_pdbx_helical_symmetry.entry_id                    %s\n",pdbid);
   fprintf(outciffile,"_pdbx_helical_symmetry.number_of_operations        %i\n",pnum );
   fprintf(outciffile,"_pdbx_helical_symmetry.rotation_per_n_subunits     %f\n",copyrot);
   fprintf(outciffile,"_pdbx_helical_symmetry.rise_per_n_subunits         %f\n",copyrise);
   fprintf(outciffile,"_pdbx_helical_symmetry.n_subunits_divisor          %i\n",divisor);
   if(xsym==1)fprintf(outciffile,"_pdbx_helical_symmetry.dyad_axis                   no\n");
   if(xsym==2)fprintf(outciffile,"_pdbx_helical_symmetry.dyad_axis                   yes\n");
   fprintf(outciffile,"_pdbx_helical_symmetry.circular_symmetry           %i\n",zsym);
           }

  //WRITE matrix tokens to cif
  fprintf(outciffile,"#\nloop_\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.id\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.type \n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[1][1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[1][2]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[1][3]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.vector[1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[2][1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[2][2]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[2][3] \n");
  fprintf(outciffile,"_pdbx_struct_oper_list.vector[2]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[3][1]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[3][2] \n");
  fprintf(outciffile,"_pdbx_struct_oper_list.matrix[3][3]\n");
  fprintf(outciffile,"_pdbx_struct_oper_list.vector[3]\n");


  //WRITE OUT pdb2pt mat to cif
   if (mtype == 'P') strcpy (matframestring,"\'transform to point frame\'");
   if (mtype == 'H') strcpy (matframestring,"\'transform to helical frame\'");
           fprintf(outciffile,"   %s     %30s\n", &mtype, matframestring);
           for (j=0;j<3;j++){
           fprintf(outciffile,"    %14.8f%14.8f%14.8f%20.5f\n",
              checkZero(p2i[j][0],PR8), 
              checkZero(p2i[j][1],PR8), 
              checkZero(p2i[j][2],PR8), 
              checkZero(p2i[j][3],PR5)); }
                                                                                                                             
  //WRITE OUT pdb2xtal mats to cif
   strcpy (matframestring,"\'transform to crystal frame\'");
    for (i=0; i<p2xnum; i++) {
           fprintf(outciffile,"   X%-4i %30s\n",i, matframestring);
           for (j=0;j<3;j++){
           fprintf(outciffile,"    %14.8f%14.8f%14.8f%20.5f\n",
              checkZero(p2xall[i][j][0],PR8), 
              checkZero(p2xall[i][j][1],PR8), 
              checkZero(p2xall[i][j][2],PR8), 
              checkZero(p2xall[i][j][3],PR5));}
           }

    //get the std point set
    if (mtype == 'P') defPointMats(exp,nmat,&stype,&zsym);
    if (mtype == 'H') defHelixMats(exp,nmat,arot,arise,(numbio/xznum-1),xsym,zsym);
    //total mats
    pnum = (exp[3]+1) * (exp[2]+1) * (exp[1]+1) * (exp[0]+1);
    calcSymmMats(exp,nmat,pmat);


   //"rewind" if helical
   if (mtype == 'H') {
   rewind = -1.0 * rewind;
   ar = arot * (M_PI/180.0);
   ptemp[0][0] = cos(ar*rewind); ptemp[0][1] = -1.0*sin(ar*rewind);    ptemp[0][2]=0.0;  ptemp[0][3]=0.0;
   ptemp[1][0] = sin(ar*rewind); ptemp[1][1] = cos(ar*rewind);         ptemp[1][2]=0.0;  ptemp[1][3]=0.0;
   ptemp[2][0] = 0.0;            ptemp[2][1] = 0.0;                    ptemp[2][2]=1.0;  ptemp[2][3]=arise*rewind;
   ptemp[3][0] = 0.0;            ptemp[3][1] = 0.0;                    ptemp[3][2]=0.0;  ptemp[3][3]=1.0;

   //rewind the matrix set
   for (i=0;i<pnum;i++) {
       multMat(ptemp,pmat[i],ptemp2);
          for (j=0;j<4;j++) for (k=0;k<4;k++) {pmat[i][j][k]=ptemp2[j][k]; }}
                    }

                                                                            
   // move std mats from point frame to the deposited frame,
   // [i2p][stdmat][p2i], 
    for (i=0;i<pnum; i++) {
       multMat(i2p,pmat[i],ptemp);
       multMat(ptemp,p2i,pmat[i]);
       }

   //print out standard point mats in deposited frame
   if (mtype == 'P') strcpy (matframestring, "\'point symmetry operation\'");
   if (mtype == 'H') strcpy (matframestring, "\'helical symmetry operation\'");
  j =0;
    for (i=0; i<pnum; i++) {
       j++;
           fprintf(outciffile,"%4i %30s\n",j, matframestring);
           for(k=0;k<3;k++){
           fprintf(outciffile,"    %14.8f%14.8f%14.8f%20.5f\n",
              checkZero(pmat[i][k][0],PR8),
              checkZero(pmat[i][k][1],PR8),
              checkZero(pmat[i][k][2],PR8),
              checkZero(pmat[i][k][3],PR5));}
           }

fprintf(outciffile,"\n##assembly information\n");
fprintf(outciffile,"loop_\n");
fprintf(outciffile,"_pdbx_struct_assembly.id\n");
fprintf(outciffile,"_pdbx_struct_assembly.details\n");
if (mtype == 'P') {
if (stype == 'I') {
fprintf(outciffile,"      1           \'complete icosahedral assembly\'\n");
fprintf(outciffile,"      2           \'icosahedral asymmetric unit\'\n");
fprintf(outciffile,"      3           \'icosahedral pentamer\'\n");
fprintf(outciffile,"      4           \'icosahedral 23 hexamer\'\n");
fprintf(outciffile,"      PAU         \'icosahedral asymmetric unit, std point frame\'\n");}
else {
fprintf(outciffile,"      1           \'complete point assembly\'\n");
fprintf(outciffile,"      2           \'point asymmetric unit\'\n");
fprintf(outciffile,"      PAU         \'point asymmetric unit, std point frame\'\n");}}
if (mtype == 'H') {
fprintf(outciffile,"      1           \'representative helical assembly\'\n");
fprintf(outciffile,"      2           \'helical asymmetric unit\'\n");
fprintf(outciffile,"      HAU         \'helical asymmetric  unit, std helical frame\'\n"); }


if (p2xnum>0){
fprintf(outciffile,"      XAU         \'crystal asymmetric unit, crystal frame\'\n");}


fprintf(outciffile,"\nloop_\n");
fprintf(outciffile,"_pdbx_struct_assembly_gen.assembly_id\n");
fprintf(outciffile,"_pdbx_struct_assembly_gen.oper_expression\n");
fprintf(outciffile,"_pdbx_struct_assembly_gen.asym_id_list\n");

if (mtype == 'P') {
  fprintf (outciffile,"     1            (1-%i)            ", pnum);
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);
  fprintf (outciffile,  "     2              1                 ");
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);


 if (stype == 'I') { 
  fprintf (outciffile,"     3            (1-5)             ");
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);
  fprintf (outciffile,"     4            (1,2,6,10,23,24)  ");
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);}


  fprintf (outciffile,  "     PAU            P                 ");
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);
                 }                                                                                                            

if (mtype == 'H') {
  fprintf (outciffile,  "     1            (1-%i)            ",pnum);
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);
  fprintf (outciffile,  "     2              %i               ",ident);
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]);
  fprintf (outciffile,  "     HAU            %c              ",mtype);
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]); fprintf(outciffile,"%1s\n",chain[chainum-1]); }


//END CIF OUT SECTION 1

  //print out BIOMT to biomats file
  printf(" writing assembly BIOMT records to *assembly.biomt*\n");
  printf("  number of assembly matrices: %i\n\n\n", pnum);
  outmatfile = fopen("assembly.biomt","w");
    for (i=0; i<pnum; i++) {
       for (k=0;k<3;k++) {
           fprintf(outmatfile,"REMARK 350   BIOMT%1i%4i%10.6f%10.6f%10.6f%15.5f\n",
           k+1, i+1, 
            checkZero (pmat[i][k][0],PR8), 
            checkZero (pmat[i][k][1],PR8),
            checkZero (pmat[i][k][2],PR8),
            checkZero (pmat[i][k][3],PR5));}}
    fclose(outmatfile);



//BEGIN PDB/CIF outputs for crystal structure
if (p2xnum>0) {
   printf(" %2i independent particle(s) in crystal unit cell\n",p2xnum);

  //initialize ncs cif file
  outncsfile = fopen("assembly.ncs","w");
  printf(" writing ncs cif records to *assembly.ncs*\n");
  fprintf(outncsfile,"## ncs operations to be applied to coordinates in the crystal frame\n");
  fprintf(outncsfile,"loop_\n");
  fprintf(outncsfile,"_struct_ncs_oper.id\n");
  fprintf(outncsfile,"_struct_ncs_oper.code \n");
  fprintf(outncsfile,"_struct_ncs_oper.details \n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[1][1]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[1][2]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[1][3]\n");
  fprintf(outncsfile,"_struct_ncs_oper.vector[1]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[2][1]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[2][2]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[2][3] \n");
  fprintf(outncsfile,"_struct_ncs_oper.vector[2]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[3][1]\n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[3][2] \n");
  fprintf(outncsfile,"_struct_ncs_oper.matrix[3][3]\n");
  fprintf(outncsfile,"_struct_ncs_oper.vector[3]\n");

  //OPEN & READ select PDB LINES TO OUTPUT XTAL FRAME FILE
   outxframefile = fopen("assembly_xframe.pdb", "w");
   printf(" writing crystal-frame ncs/coord file to *assembly_xframe.pdb*\n\n");
   strncpy(cifstring,spacegroup+1,(strlen(spacegroup)-2));
   itoa(snum,tmpstring);
   if (strlen(tmpstring)<1 || strlen(tmpstring)>3) return 1;
   if (strlen(tmpstring)==1) {for (i=(strlen(spacegroup)-2);i<14;i++) cifstring[i]=' '; cifstring[14]='\0'; strcat(cifstring,tmpstring);  }
   if (strlen(tmpstring)==2) {for (i=(strlen(spacegroup)-2);i<13;i++) cifstring[i]=' '; cifstring[13]='\0'; strcat(cifstring,tmpstring);  }
   if (strlen(tmpstring)==3) {for (i=(strlen(spacegroup)-2);i<12;i++) cifstring[i]=' '; cifstring[12]='\0'; strcat(cifstring,tmpstring);  }
   cifstring[15]='\0';
   fprintf(outxframefile,"CRYST1 %8.3f %8.3f %8.3f %6.2f %6.2f %6.2f %s\n", unitcell[0], unitcell[1],unitcell[2], unitcell[3],
            unitcell[4],unitcell[5],cifstring);
for(i=0;i<3;i++)
   fprintf(outxframefile,"SCALE%i     %9.6f %9.6f %9.6f %14.5f\n", i+1, fracmat[i][0], fracmat[i][1], fracmat[i][2], fracmat[i][3]);


//need to write out valid CRYST1 and SCALE matrix from cif tokens


}  //end if p2xnum > 0



//**************************************** 
//begin big loop over each pdb2xtal matrix
//**************************************** 

    ncsnum = 0;  //initialize counter for ncs matrices

 for (pcount=0;pcount<p2xnum;pcount++){

    printf ("*************************************************\n");
    printf (" PARTICLE %2i \n", (pcount+1),(pcount+1));

    //get the std symmetry matrix set
    calcSymmMats(exp,nmat,pmat);

  //assign p2x mat for this round
        for (j=0;j<4;j++) for (k=0;k<4;k++) { p2x[j][k] = p2xall[pcount][j][k];}

  //calculate the x2p matrix
    det=invMat (p2x, x2p);

  //calculate the i2x matrix
    multMat (p2x, i2p, i2x);

    printf (" deposited frame to crystal frame transformation:\n" );
    for (i=0;i<3;i++) { 
    printf (" %1i%10.6f%10.6f%10.6f%15.5f\n",
          (pcount+1),p2x[i][0],p2x[i][1],p2x[i][2],p2x[i][3]); }


  //calculate the x2i matrix, for 1st x2p only
  if (pcount==0) { multMat (p2i, x2p, x2i);}


   // move std mats from point frame into the crystal frame,
   // [i2x][stdmat][x2i], 
    for (i=0;i<pnum; i++) {
       multMat(i2x,pmat[i],ptemp);
       multMat(ptemp,x2i,pmat[i]);
       }

  //make translation-only versions of i2x and x2i with fractional translations
    for (j=0;j<3;j++){
        i2xt[j][3]= fracmat[j][0] *i2x[0][3] + fracmat[j][1] *i2x[1][3] + fracmat [j][2] *i2x[2][3];}
    i2xt[0][0]=1; i2xt[0][1]=0; i2xt[0][2]=0; 
    i2xt[1][0]=0; i2xt[1][1]=1; i2xt[1][2]=0; 
    i2xt[2][0]=0; i2xt[2][1]=0; i2xt[2][2]=1; 
    i2xt[3][0]=0; i2xt[3][1]=0; i2xt[3][2]=0; i2xt[3][3]=1;
    det=invMat(i2xt,x2it);

    //printf (" translation from origin in fractional coordinates:\n");
    //printf (" %12.6f%12.6f%12.6f\n", i2xt[0][3],i2xt[1][3],i2xt[2][3]);

  //initialize lists that keep track of operator status
       for (i=0;i<pnum;i++){ plist[i] = 1; }
       for (i=0;i<snum;i++){ stlist[i] = 0; }  //first pass list:  symops w/o translations
       for (i=0;i<snum;i++){ srlist[i] = 0; }  //second pass list: symops giving equivalent positions
    //always include identity mat
       stlist[0] = 1;
       srlist[0] = 1;

 
     for (i=0; i<snum; i++) {   //begin loop over symops

        //transform cartesion SYMOPS to their fractional form
       // cartesian SYMOPS are really [orth][symop-frac][frac]
       // to get at [symop-frac], we do [frac][SYMOPcart][orth]=[frac][orth][symop-frac][frac][orth]

         multMat(fracmat,smat[i],ptemp);
         multMat(ptemp,orthmat,smat[i]);

        //transform symop by inverse translation of particle
        //this  moves crystal symop to frame where the particle is at the origin
        //e.g., for C 2 2 21, x,y,z and -x,y,-z+1/2 are related by two-fold on x=0, z=1/4.
        //by moving the coordinate system such that z=1/4 --> z=0, the operation pair
        //becomes  x,y,z and -x,y,-z.

         multMat(x2it,smat[i],ptemp);
         multMat(ptemp,i2xt,smat[i]);

        //if fractional translation is close to an integer,  reset it to zero.
        for (j=0;j<3;j++){
          transtest = fabs(smat[i][j][3]);
          if ((transtest-floor(transtest))<T_PRECISION || (ceil(transtest)-transtest)<T_PRECISION) smat[i][j][3]=0;
           }

          //printf ("  mod %2i frac: %12.8f%12.8f%12.8f\n", i+1,smat[i][0][3],smat[i][1][3],smat[i][2][3]);

         //test for symops with remaining translations
          transtest = fabs(smat[i][0][3]) + fabs(smat[i][1][3]) + fabs(smat[i][2][3]) ;

          if (transtest < T_PRECISION) { stlist[i]= 1; }
         
         //reset symop to original cartesian transform
          for (j=0;j<4;j++) for (k=0;k<4;k++) { smat[i][j][k] = smfixed[i][j][k];}


           }  //end of loop over symops



   //apply symmetry mats to point mats, save in psmat
      for (i=0; i<snum; i++){
          if (stlist[i]==1) {
             for (j=0 ; j< pnum ; j++) {
                multMat (smat[i], pmat[j], psmat[i][j]);
        }}}


   //eliminate redundant mats
   //report symops that yield redundant mats
      for (i=0; i<snum; i++){
          if (stlist[i]==1) {
             for (j=0 ; j<pnum ; j++) 
                 for (k=(j+1); k<pnum; k++) {
                    diff = equiMat(psmat[i][k], pmat[j]);
                    if (diff <= M_PRECISION){
                         plist[k] = 0 ;
                         srlist[i] = 1 ; }}}
        }

  //print remark with selected symmetry indices
   j =0;
  printf (" crystal symops contributing to point symmetry:");
    for (i=0; i<snum; i++) {
       if (srlist[i] == 1) {
          if (j%15 == 0.0) printf ("\n "); 
          printf("%3i", (i+1));
          j++;
          }}
  printf ("  of %3i crystal symops total\n", snum);

 
  //print remark with selected pmat indices
  printf (" point operations corresponding to ncs ops:");
   j =0;
    for (i=0; i<pnum; i++) {
       if (plist[i] == 1) {
          if (j%15 == 0.0) printf ("\n "); 
          printf("%3i", (i+1));
          j++;
          }}
  printf ("\n %3i of %3i point ops \n", j, pnum );


  fprintf (outciffile,"     XAU          (X%i)",pcount);
  fprintf (outciffile,"(");
    for (i=0; i<pnum; i++) {
          if (plist[i]==1) {
              if (i == 0) fprintf(outciffile,"%i", (i+1)); // "1"
              if (i == 1 && plist[i-1]==0) fprintf(outciffile,"%i", (i+1)); //"2"
              if (i == 1 && plist[i-1]==1 && plist[i+1]==0) fprintf(outciffile,",%i", (i+1)); // "1,2"
              if (i == 1 && plist[i-1]==1 && plist[i+1]==1) fprintf(outciffile,"-"); //"1-"
              if (i > 1 && i < pnum-1)  {
                 if (plist[i-1]==0) fprintf(outciffile,",%i", (i+1)); //singleton or begin new run 
                 if (plist[i-2] == 0 && plist[i-1]==1 && plist[i+1]==1) fprintf(outciffile,"-", (i+1)); //run identified
                 if (plist[i-2] == 1 && plist[i-1]==1 && plist[i+1]==0) fprintf(outciffile,"%i", (i+1));//end run identified
                                        }
             if (i == pnum-1 && plist[i-2]==1 && plist[i-1]==1) fprintf(outciffile,"%i",(i+1)); //end last run
             if (i == pnum-1 && plist[i-2]==0 && plist[i-1]==1) fprintf(outciffile,",%i",(i+1)); //last singleton
             if (i == pnum-1 && plist[i-2]==1 && plist[i-1]==0) fprintf(outciffile,",%i",(i+1)); //last singleton
          }}
 
  fprintf (outciffile,")    ");
  for (i=0;i<chainum-1;i++) fprintf(outciffile,"%1s,",chain[i]);
  fprintf(outciffile,"%1s\n",chain[chainum-1]);


  //printf(  "REMARK PNCS  APPLY THE FOLLOWING TO COORDINATES TRANSFORMED BY 2XTAL MAT 1\n");

  //print out the selected ncs matrices as MTRIX cards to PDB and CIF records to separate file
    for (i=0; i<pnum; i++) {
       if (plist[i] == 1) {
       ncsnum++;
       fprintf(outxframefile, "MTRIX1%4i%10.6f%10.6f%10.6f   %12.5f\nMTRIX2%4i%10.6f%10.6f%10.6f   %12.5f\nMTRIX3%4i%10.6f%10.6f%10.6f   %12.5f\n",
       ncsnum, pmat[i][0][0], pmat[i][0][1], pmat[i][0][2], pmat[i][0][3],
       ncsnum, pmat[i][1][0], pmat[i][1][1], pmat[i][1][2], pmat[i][1][3],
       ncsnum, pmat[i][2][0], pmat[i][2][1], pmat[i][2][2], pmat[i][2][3]);

           if (ncsnum == 1) fprintf(outncsfile,"%4i    given ?\n",ncsnum);
           if (ncsnum >  1) fprintf(outncsfile,"%4i generate ?\n",ncsnum);
           for (j=0;j<3;j++){
           fprintf(outncsfile,"    %14.8f%14.8f%14.8f%20.5f\n",
              checkZero(pmat[i][j][0],PR8),
              checkZero(pmat[i][j][1],PR8),
              checkZero(pmat[i][j][2],PR8),
              checkZero(pmat[i][j][3],PR5));}}
       }}

//*******************************
//end big loop over p2x matrices
//*******************************
if (p2xnum>0) {

/*
    fprintf(outxframefile,"%30s%8.3f%8.3f%8.3f%26s",linestart,xyz[0],xyz[1],xyz[2],lineend);
    i++;
*/
/*
  fclose(outxframefile);
  fclose(outncsfile);
}
   printf("closing assembly.cif\n");
   fclose(outciffile);
*/
  return 0;
}
