// pointsuite 0.5.8 
// author Huanwang Yang 
// RCSB-PDB
// JUNE 2007
//now with double precision, 20 Jun 2011

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cifparse.c"

int strncmp_case (const char *s1, const char *s2, int n);
void reformat_file(char *inpfile);

void write_matrix(FILE *fpdb, int ncs, double rt[]);
void write_remark350(FILE *fpdb, int ncs, double rt[]);
void write_remark350_cif(FILE *fcif, int ncs, double rt[]);



int main(int argc, char * argv[]){

    if(argc==1){
        printf("This program convert the matrix files to the standard format.\nIt supports 10 different matrix formats.\n");
        printf("\nUsage:  importmats  matrix_file\n\n");
        exit(0);
    }
    
        
    printf ("\n******converting matrix file (%s) to CIF/BIOMT/MATRIX*******\n", argv[1]);
    reformat_file(argv[1]);
    printf("\nThe output files = import.biomt & import.cif & import.matrix\n");
    exit(0);
}

void reformat_file(char *inpfile)
/*  get the correct format ready for matrix extraction */
{
    int ncs=0,i,j, n=0, nc1=0, nc2=0, nc3=0, nc4=0, nn;
    double rt[14];
    char str[256], tmp[256];
    char **m11,**m12,**m13,**t1,**m21,**m22,**m23,**t2, **m31,**m32,**m33,**t3,**mtype, **matr;
    
    FILE *fp,*fpdb,*fcif,*fmtr;

    remove("import.biomt");
    remove("import.cif");
    remove("import.matrix");
    
    fpdb=fopen("import.biomt","w");
    fcif=fopen("import.cif","w");
    fmtr=fopen("import.matrix","w");
    
    if((fp=fopen(inpfile, "r"))==NULL){
        printf("The file=%s can not be opened!\n", inpfile);
        exit(0);
    }

    if (is_cif(inpfile) && cifparse(inpfile, "_pdbx_struct_oper_list.")){
        mtype=parse_values("_pdbx_struct_oper_list.type", &nn);
        m11=parse_values("_pdbx_struct_oper_list.matrix[1][1]", &nn);
        m12=parse_values("_pdbx_struct_oper_list.matrix[1][2]", &nn);
        m13=parse_values("_pdbx_struct_oper_list.matrix[1][3]", &nn);
        t1 =parse_values("_pdbx_struct_oper_list.vector[1]", &nn);
        m21=parse_values("_pdbx_struct_oper_list.matrix[2][1]", &nn);
        m22=parse_values("_pdbx_struct_oper_list.matrix[2][2]", &nn);
        m23=parse_values("_pdbx_struct_oper_list.matrix[2][3]", &nn);
        t2 =parse_values("_pdbx_struct_oper_list.vector[2]", &nn);
        m31=parse_values("_pdbx_struct_oper_list.matrix[3][1]", &nn);
        m32=parse_values("_pdbx_struct_oper_list.matrix[3][2]", &nn);
        m33=parse_values("_pdbx_struct_oper_list.matrix[3][3]", &nn);
        t3 =parse_values("_pdbx_struct_oper_list.vector[3]", &nn);
        strcpy(tmp,"REMARK 350   BIOMT");
        j=0;
        for (i=0;  i < nn; i++) {
/*            
            if (!strstr(mtype[i], "point symmetry operation") &&
                !strstr(mtype[i], "helical symmetry operation") ) continue;
*/       
        
            j++;
            rt[0]=atof(m11[i]);
            rt[1]=atof(m12[i]);
            rt[2]=atof(m13[i]);
            
            rt[3]=atof(m21[i]);
            rt[4]=atof(m22[i]);
            rt[5]=atof(m23[i]);
            
            rt[6]=atof(m31[i]);
            rt[7]=atof(m32[i]);
            rt[8]=atof(m33[i]);
            
            rt[9] =atof(t1[i]);
            rt[10]=atof(t2[i]);
            rt[11]=atof(t3[i]);
            
            write_matrix(fmtr, j, rt);
            write_remark350(fpdb, j, rt);
            write_remark350_cif(fcif, j, rt);
        }
        
            
    } else if (is_cif(inpfile) && cifparse(inpfile, "_pdbx_struct_assembly_gen_depositor_info.")){
        matr=parse_values("_pdbx_struct_assembly_gen_depositor_info.full_matrices", &nn);
        j=0;
        for (i=0;  i < nn; i++) {
            j++;
            
            nc1=sscanf( matr[i],"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",  &rt[0], &rt[1], &rt[2], &rt[9], &rt[3], &rt[4], &rt[5], &rt[10],   &rt[6], &rt[7], &rt[8], &rt[11]);
            
            write_matrix(fmtr, j, rt);
            write_remark350(fpdb, j, rt);
            write_remark350_cif(fcif, j, rt);
/*            printf("%d  %s\n", nc1, matr[i]); */
        }
        
    }else if (is_cif(inpfile)){
        printf("Error, input cif file has no pdbx_struct_oper_list or pdbx_struct_assembly_gen_depositor_info.\n");
    }

    if (is_cif(inpfile)){
        fclose(fpdb);
        fclose(fcif);
        fclose(fp);
        return;
    }
    
    
    

    while(fgets(str, sizeof str,fp)!=NULL){
        rid_of_front_space(str);
        for(i=0; i<14; i++)rt[i]=0;
        nc1=sscanf(str,"%lf%lf%lf%lf", &rt[0], &rt[1], &rt[2], &rt[9]);

/* PDBSET */
        if(!strncmp_case(str, "transform ",10) && strstr(str, " -")){
            ncs++;
            
            for(i=0; i<14; i++)rt[i]=0;
            strcpy(tmp,strstr(str,"form")+4);
            
            sscanf(tmp,"%lf%lf%lf", &rt[0], &rt[1], &rt[2]);

            fgets(str, sizeof str,fp);
            sscanf(str,"%lf%lf%lf", &rt[3], &rt[4], &rt[5]);
            
            fgets(str, sizeof str,fp);
            sscanf(str,"%lf%lf%lf", &rt[6], &rt[7], &rt[8]);
            
            fgets(str, sizeof str,fp);
            sscanf(str,"%lf%lf%lf", &rt[9], &rt[10], &rt[11]);
            write_matrix(fmtr, ncs, rt);
            write_remark350(fpdb, ncs, rt);
            write_remark350_cif(fcif, ncs, rt);

        }else if(strstr(str,"ROTATION MATRIX:")){
            for(i=0; i<14; i++)rt[i]=0;
            ncs++;
            fgets(str, sizeof str,fp);
            sscanf(str,"%lf%lf%lf", &rt[0], &rt[1], &rt[2]);
            printf("str1=%s\n", str);
            
            
            fgets(str, sizeof str,fp);
            sscanf(str,"%lf%lf%lf", &rt[3], &rt[4], &rt[5]);
            printf("str2=%s\n", str);
            
            fgets(str, sizeof str,fp);
            sscanf(str,"%lf%lf%lf", &rt[6], &rt[7], &rt[8]);
            printf("str3=%s\n", str);
            
        }else if(strstr(str,"TRANSLATION VECTOR IN AS")){
            
            sscanf(str+27,"%lf%lf%lf", &rt[9], &rt[10], &rt[11]);
            printf("str4=%s\n", str);
            write_matrix(fmtr, ncs, rt);
            write_remark350(fpdb, ncs, rt);
            write_remark350_cif(fcif, ncs, rt);
            
            
            
/* CNS matrix */                
        }else if(strstr(str,"things below this line do not normally need to be changed")){
            break;
                     
        }else if(!strncmp(str,"{===>}", 6) && strstr(str,"ncs_op_") && strstr(str,"true")){
            ncs++;
            n=0;
            for(i=0; i<14; i++)rt[i]=0;
            while(fgets(str, sizeof str,fp)!=NULL){
                rid_of_front_space(str);
                if(strlen(str)<3 || !strncmp(str,"{*", 2) ) continue;
                
                if(!strncmp(str,"{===>}", 6) && strstr(str,"ncs_vector_") && strchr(str, '(')){
                    if(strchr(str,'(')) {
                        strcpy(tmp,strchr(str,'(') +1);
                        sscanf(tmp,"%lf%lf%lf", &rt[9], &rt[10], &rt[11]);
/*                        printf("The trans=%2d  %s\n",ncs, tmp); */
                    }
/*                    for(i=0; i<12; i++)printf("%.2f \n",rt[i] );*/
                    
                    write_matrix(fmtr, ncs, rt);
                    write_remark350(fpdb, ncs, rt);
                    write_remark350_cif(fcif, ncs, rt); 
                    break;
                }
                if(!strncmp(str,"{===>}", 6) && strstr(str,"ncs_matrix_") && strchr(str, '(')){
                    if(strchr(str,'(')) {
                        strcpy(tmp,strchr(str,'(') +1);
                        sscanf(tmp,"%lf%lf%lf", &rt[0], &rt[1], &rt[2]);
                    }
                    
                    fgets(str, sizeof str,fp);
                    if(strchr(str,'(')) {
                        strcpy(tmp,strchr(str,'(') +1);
                        sscanf(tmp,"%lf%lf%lf", &rt[3], &rt[4], &rt[5]);
                    }
                    
                    fgets(str, sizeof str,fp);
                    if(strchr(str,'(')) {
                        strcpy(tmp,strchr(str,'(') +1);
                        sscanf(tmp,"%lf%lf%lf", &rt[6], &rt[7], &rt[8]);
                    }
                }
            }
            
/* 60 matrix, no identifer :
   3 columns, 4 rows for each matrix. last row for tran!
   matrix rt[9], [10], [11] for translation
*/            
        }else if(strstr(str, "xncsrel")){
            while(fgets(str, sizeof str,fp)!=NULL){
                
                if(strstr(str,"translation")){
                    nc4=sscanf(strchr(str,'(')+1,"%lf%lf%lf",  &rt[9], &rt[10], &rt[11]);
                    if(nc1==3 && nc2==3 && nc3==3 && nc4==3){
                        ncs++;
                        write_matrix(fmtr, ncs, rt);
                        write_remark350(fpdb, ncs, rt);
                        write_remark350_cif(fcif, ncs, rt);
                        
                        nc1=0;
                        nc2=0;
                        nc3=0;
                        nc4=0;
                        
                        break;
                    }
                }else if(strstr(str,"matrix")){
                    if (!strchr(str,'(') && !strchr(str,')'))  fgets(str, sizeof str,fp);
                    nc1=sscanf(strchr(str,'(')+1,"%lf%lf%lf", &rt[0], &rt[1], &rt[2]);

                    fgets(str, sizeof str,fp);
                    nc2=sscanf(strchr(str,'(')+1,"%lf%lf%lf", &rt[3], &rt[4], &rt[5]);
                    
                    fgets(str, sizeof str,fp);
                    nc3=sscanf(strchr(str,'(')+1,"%lf%lf%lf",  &rt[6], &rt[7], &rt[8]);
                    
                }
            }

           
        }else if(nc1==3){
            fgets(str, sizeof str,fp);
            nc2=sscanf(str,"%lf%lf%lf", &rt[3], &rt[4], &rt[5]);
            
            fgets(str, sizeof str,fp);
            nc3=sscanf(str,"%lf%lf%lf",  &rt[6], &rt[7], &rt[8]);
            
            if(fgets(str, sizeof str,fp)!=NULL){
                nc4=sscanf(str,"%lf%lf%lf",  &rt[9], &rt[10], &rt[11]);
//                printf("nc1=%d %d %d %d \n", nc1,nc2, nc3,nc4  );
            }
            
            if(nc1==3 && nc2==3 && nc3==3 && nc4==3){
                ncs++;
                write_matrix(fmtr, ncs, rt);
                write_remark350(fpdb, ncs, rt);
                write_remark350_cif(fcif, ncs, rt);
            }else if(nc1==3 && nc2==3 && nc3==3 && nc4!=3){ /* for rotation only */
                rt[9]=0;
                rt[10]=0;
                rt[11]=0;
                ncs++;
                write_matrix(fmtr, ncs, rt);
                write_remark350(fpdb, ncs, rt);
                write_remark350_cif(fcif, ncs, rt);
            }
            
            
            nc1=0;
            nc2=0;
            nc3=0;
            nc4=0;

/* 60 matrix, no identifer :
   4 columns, 3 rows for each matrix. !
   matrix rt[9], [10], [11] for translation
*/            
            
        }else if(nc1==4){
            fgets(str, sizeof str,fp);
            nc2=sscanf(str,"%lf%lf%lf%lf", &rt[3], &rt[4], &rt[5], &rt[10]);
            
            fgets(str, sizeof str,fp);
            nc3=sscanf(str,"%lf%lf%lf%lf",  &rt[6], &rt[7], &rt[8], &rt[11]);
            
            if(nc1==4 && nc2==4 && nc3==4){
                ncs++;
                write_matrix(fmtr, ncs, rt);
                write_remark350(fpdb, ncs, rt);
                write_remark350_cif(fcif, ncs, rt);
            }else{
                
                printf("Warnning! matrix is not 3(row)X4(column), please check output!\n");
            }
            
            
            nc1=0;
            nc2=0;
            nc3=0;
        
/* only  MTRIX1  in first column */           
        }else if(!strncmp(str,"MTRIX1", 6)){
            nc1=0;
            nc2=0;
            nc3=0;
            
            strcpy(tmp,strstr(str,"MTRIX1") +6);
            nc1=sscanf(tmp,"%*s%lf%lf%lf%lf", &rt[0], &rt[1], &rt[2], &rt[9]);
            
            fgets(str, sizeof str,fp);
            if(strstr(str,"MTRIX2")){
                strcpy(tmp,strstr(str,"MTRIX2") +6);
                nc2=sscanf(tmp,"%*s%lf%lf%lf%lf", &rt[3], &rt[4], &rt[5], &rt[10]);
            }
            
            fgets(str, sizeof str,fp);
            if(strstr(str,"MTRIX3")){
                strcpy(tmp,strstr(str,"MTRIX3") +6);
                nc3=sscanf(tmp,"%*s%lf%lf%lf%lf",&rt[6], &rt[7], &rt[8],  &rt[11]);
            }
            if(nc1==4 && nc2==4 && nc3==4){
                ncs++;
                write_matrix(fmtr, ncs, rt);
                write_remark350(fpdb, ncs, rt);
                write_remark350_cif(fcif, ncs, rt);
            }
            
/* only  BIOMT1  in first column */           
        }else if(!strncmp(str,"BIOMT1", 6) ||
                 (!strncmp(str,"REMARK", 6) && strstr(str,"BIOMT1"))){
            nc1=0;
            nc2=0;
            nc3=0;
            
            strcpy(tmp,strstr(str,"BIOMT1") +6);
            nc1=sscanf(tmp,"%*s%lf%lf%lf%lf", &rt[0], &rt[1], &rt[2], &rt[9]);
            
            fgets(str, sizeof str,fp);
            if(strstr(str,"BIOMT2")){
                strcpy(tmp,strstr(str,"BIOMT2") +6);
                nc2=sscanf(tmp,"%*s%lf%lf%lf%lf", &rt[3], &rt[4], &rt[5], &rt[10]);
            }
            
            fgets(str, sizeof str,fp);
            if(strstr(str,"BIOMT3")){
                strcpy(tmp,strstr(str,"BIOMT3") +6);
                nc3=sscanf(tmp,"%*s%lf%lf%lf%lf",&rt[6], &rt[7], &rt[8],  &rt[11]);
            }
            if(nc1==4 && nc2==4 && nc3==4){
                ncs++;
                write_matrix(fmtr, ncs, rt);
                write_remark350(fpdb, ncs, rt);
                write_remark350_cif(fcif, ncs, rt);
            }
            
        }else if(strstr(str,"REMARK Icosahedral Matrix") && strstr(str,"applied to")){
                // by VIPERdb: Oligomer generation
            nc1=0; 
            nc2=0;
            nc3=0;
            nc4=0;
            
            fgets(str, sizeof str,fp);
            nc1=sscanf(str,"%*s%lf%lf%lf", &rt[0], &rt[1], &rt[2]);
            fgets(str, sizeof str,fp);
            nc2=sscanf(str,"%*s%lf%lf%lf", &rt[3], &rt[4], &rt[5]);
            fgets(str, sizeof str,fp);
            nc3=sscanf(str,"%*s%lf%lf%lf", &rt[6], &rt[7], &rt[8]);
            fgets(str, sizeof str,fp);
            fgets(str, sizeof str,fp);
            nc4=sscanf(str,"%*s%lf%lf%lf", &rt[9], &rt[10], &rt[11]);
            
            if(nc1==3 && nc2==3 && nc3==3 && nc4==3){
                ncs++;
                write_matrix(fmtr, ncs, rt);
                write_remark350(fpdb, ncs, rt);
                write_remark350_cif(fcif, ncs, rt);
            }
            
        }
    }
    fclose(fpdb);
    fclose(fcif);
    fclose(fp);
}


void write_remark350(FILE *fpdb, int ncs, double rt[])
{
    fprintf(fpdb,"REMARK 350   BIOMT1 %3d%10.6lf%10.6lf%10.6lf%15.5lf\n", ncs, rt[0], rt[1], rt[2],rt[9]);
    fprintf(fpdb,"REMARK 350   BIOMT2 %3d%10.6lf%10.6lf%10.6lf%15.5lf\n", ncs, rt[3], rt[4], rt[5],rt[10]);
    fprintf(fpdb,"REMARK 350   BIOMT3 %3d%10.6lf%10.6lf%10.6lf%15.5lf\n", ncs, rt[6], rt[7], rt[8],rt[11]);
}

void write_matrix(FILE *fpdb, int ncs, double rt[])
{
    fprintf(fpdb,"MTRIX1 %3d%10.6lf%10.6lf%10.6lf%15.5lf\n", ncs, rt[0], rt[1], rt[2],rt[9]);
    fprintf(fpdb,"MTRIX2 %3d%10.6lf%10.6lf%10.6lf%15.5lf\n", ncs, rt[3], rt[4], rt[5],rt[10]);
    fprintf(fpdb,"MTRIX3 %3d%10.6lf%10.6lf%10.6lf%15.5lf\n", ncs, rt[6], rt[7], rt[8],rt[11]);
}

void write_remark350_cif(FILE *fcif, int ncs, double rt[])
{
    if(ncs==1){
        fprintf(fcif,"loop_\n");
        fprintf(fcif,"_pdbx_struct_oper_list.id\n");              
        fprintf(fcif,"_pdbx_struct_oper_list.type\n");          
        fprintf(fcif,"_pdbx_struct_oper_list.name\n");          
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[1][1]\n"); 
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[1][2]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[1][3]\n"); 
        fprintf(fcif,"_pdbx_struct_oper_list.vector[1]\n");   
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[2][1]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[2][2]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[2][3]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.vector[2]\n");  
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[3][1]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[3][2]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.matrix[3][3]\n");
        fprintf(fcif,"_pdbx_struct_oper_list.vector[3]\n");  
    }
    fprintf(fcif,"%i \"general operation\" . \n",ncs);
    fprintf(fcif,"%14.8lf%14.8lf%14.8lf%14.5lf\n", rt[0], rt[1], rt[2],rt[9]);
    fprintf(fcif,"%14.8lf%14.8lf%14.8lf%14.5lf\n", rt[3], rt[4], rt[5],rt[10]);
    fprintf(fcif,"%14.8lf%14.8lf%14.8lf%14.5lf\n", rt[6], rt[7], rt[8],rt[11]);

}



int strncmp_case (const char *s1, const char *s2, int n){

  int same =0;
  char a1, a2;
  int i=0;

  i=0;
  while (s1[0]!='\0' && s2[0]!='\0' && !same && i<n) {
    a1 = tolower(s1[0]);
    a2 = tolower(s2[0]);
    if (a1!=a2)
    {
      if (a1<a2)
	same = 1;
      else
	same = -1;
    }
  s1++; s2++;
  i++;
  }
  if ((s1[0]!='\0'|| s2[0]!='\0')&& same==0)
    if (i<n) {
      if (s1[0]!='\0')
	same = 1;
      else
	same = -1;
    }
  return same;
}

