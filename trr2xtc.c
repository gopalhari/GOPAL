/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id: trr2xtc.c,v 1.3 2009/05/18 09:06:38 spoel Exp $
 *
 * Copyright (c) Erik Lindahl, David van der Spoel 2003,2004.
 * Coordinate compression (c) by Frans van Hoesel. 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
#include <signal.h>
#define PI 3.141592653589793
/* This program tests reading and writing to XDR files */

static void _die(char *msg, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die(msg) _die(msg,__LINE__,__FILE__)

static void _die_r(char *msg, int result, int line, char *file) {
	fprintf(stderr, "Fatal error: %s\n", msg);
	fprintf(stderr,"result = %d\n", result);
	fprintf(stderr, "Death occurred at %s, line %d\n", file, line);
	exit(1);
}
#define die_r(msg,res) _die_r(msg,res,__LINE__,__FILE__)


void ReadWrite(char *rfile, char *wfile, int in_xtcBool, int out_xtcBool, int in_trrBool, int out_trrBool, double romi, double rohi, double angi,  double rohj, double angj )
{
  XDRFILE *xd_read, *xd_write;
  int result_xtc, result_trr;
  int natoms_xtc, natoms_trr;
  int step_xtc, step_trr;
  float time_xtc, time_trr;
  matrix box_xtc, box_trr;
  rvec *x_xtc, *x_trr, *v_trr, *f_trr;
  float prec_xtc = 1000.0;
  float lambda_trr = 0.0;
  
      
  xd_read = xdrfile_open(rfile, "r");
  if (NULL == xd_read)
    die("Opening xdrfile for reading");
  
  /* Test whether output file exists */
  if ((xd_write = xdrfile_open(wfile,"r")) != NULL) {
    xdrfile_close(xd_write);
    die("Output file exists.");
  }
  
  /* Output file does not exist. Now we can open it for writing */
  xd_write = xdrfile_open(wfile, "w");
  if (NULL == xd_write)
    die("Opening xdrfile for writing");
  
  
  /* .xtc -> .xtc */
  if(in_xtcBool && out_xtcBool)
    {
      result_xtc = read_xtc_natoms(rfile, &natoms_xtc);
      if (exdrOK != result_xtc)
	die_r("read_xtc_natoms",result_xtc);
      
      x_xtc = calloc(natoms_xtc, sizeof(x_xtc[0]));

      while(1)
	{
	  result_xtc = read_xtc(xd_read, natoms_xtc, &step_xtc, &time_xtc,
				box_xtc, x_xtc, &prec_xtc);
	    
	  if (result_xtc == 0) // if not reach the end of file, write it to the output.xtc file
	    {
	      if (exdrOK != result_xtc)
		die_r("Reading xtc file", result_xtc);
	      
	      result_xtc = write_xtc(xd_write, natoms_xtc, step_xtc, time_xtc,
				     box_xtc, x_xtc, prec_xtc);
	      
	      if (result_xtc != 0)
		die_r("Writing xtc file", result_xtc);
	    }
	  else
	    break;
	}
    }
  

  /* .xtc -> .trr */
  if(in_xtcBool && out_trrBool)
    {
      result_xtc = read_xtc_natoms(rfile, &natoms_xtc);
      if (exdrOK != result_xtc)
	die_r("read_xtc_natoms",result_xtc);
      
      x_xtc = calloc(natoms_xtc, sizeof(x_xtc[0]));
      
      while(1)
	{
	  result_xtc = read_xtc(xd_read, natoms_xtc, &step_xtc, &time_xtc,
				box_xtc, x_xtc, &prec_xtc);
	  
	  if (result_xtc == 0) // if not reach the end of file, write it to the output.trr file
	    {
	      if (exdrOK != result_xtc)
		die_r("Reading xtc file", result_xtc);
	      
	      result_trr = write_trr(xd_write, natoms_xtc, step_xtc, time_xtc, lambda_trr,
				     box_xtc, x_xtc, NULL, NULL);
	      
	      if (0 != result_trr)
		die_r("Writing trr file",result_trr);
	      
	    }
	  else
	    break;
	}
    }
  
  
  /* .trr -> .trr */
  if(in_trrBool && out_trrBool)
    {
      result_trr = read_trr_natoms(rfile, &natoms_trr);
      
      if (exdrOK != result_trr)
	die_r("read_trr_natoms",result_trr);
      
      x_trr = calloc(natoms_trr, sizeof(x_trr[0]));
      v_trr = calloc(natoms_trr, sizeof(v_trr[0]));
      f_trr = calloc(natoms_trr, sizeof(f_trr[0]));
          
	  double m1=15.9994; /*atomic mass of oxygen*/
          double m2=1.008; /*atomic mass of hydrogen*/
          double m3=1.008; /*atomic mass of hydrogen*/
          double M ;
          M = m1+m2+m3;
          /*double romi,rohi,angi,rohj,angj; */
          double OH1 = 0.09572;
          double OH0 = 0.1;
          double HH1 = 0.151390;
	  double HH0 = 0.163300;
	  /*double OHi = OH0*(1.0-lami)+OH1*lami; 
	  double OHj = OH0*(1.0-lamj)+OH1*lamj;*/
	  double OHi = rohi; 
	  double OHj = rohj;
	  double OHo = OH1;
         
	  /*double HHi = HH0*(1.0-lami)+HH1*lami; 
	  double HHj = HH0*(1.0-lamj)+HH1*lamj;*/
          double sinthetaby2i = sin(angi/2.0*PI/180.0);
          double costhetaby2i = cos(angi/2.0*PI/180.0);
          double sinthetaby2j = sin(angj/2.0*PI/180.0);
          double costhetaby2j = cos(angj/2.0*PI/180.0);
          double rhhi = rohi*sinthetaby2i*2.0;
          double rhhj = rohj*sinthetaby2j*2.0; 
	  double HHi = rhhi; 
	  double HHj = rhhj;
	  double HHo = HH1;

	  double ai = romi/2/rohi/costhetaby2i;
          double bi = ai;
          /*ai = 0.0;
          bi = 0.0;*/ 
          double Kio = OHi/OHo*(sqrt((1.0-(HHi/2.0/OHi)*(HHi/2.0/OHi))/(1.0-(HHo/2.0/OHo)*(HHo/2.0/OHo))));
          double Koj = OHo/OHj*(sqrt((1.0-(HHo/2.0/OHo)*(HHo/2.0/OHo))/(1.0-(HHj/2.0/OHj)*(HHj/2.0/OHj))));
          
          double bc = -m1/M;
        
          double l11oj = m1/M+(1.0+bc)*Koj;           
          double l12oj = m2/M-(1.0+bc)*Koj/2.0;           
          double l13oj = m3/M-(1.0+bc)*Koj/2.0;           
          double l21oj = m1/M+bc*Koj;           
          double l22oj = m2/M-bc*Koj/2.0+0.5*HHo/HHj;           
          double l23oj = m3/M-bc*Koj/2.0-0.5*HHo/HHj;           
          double l31oj = m1/M+bc*Koj;           
          double l32oj = m2/M-bc*Koj/2.0-0.5*HHo/HHj;           
          double l33oj = m3/M-bc*Koj/2.0+0.5*HHo/HHj;           
 

          double l11io = m1/M+(1.0+bc)*Kio;           
          double l12io = m2/M-(1.0+bc)*Kio/2.0;           
          double l13io = m3/M-(1.0+bc)*Kio/2.0;           
          double l21io = m1/M+bc*Kio;           
          double l22io = m2/M-bc*Kio/2.0+0.5*HHi/HHo;           
          double l23io = m3/M-bc*Kio/2.0-0.5*HHi/HHo;           
          double l31io = m1/M+bc*Kio;           
          double l32io = m2/M-bc*Kio/2.0-0.5*HHi/HHo;           
          double l33io = m3/M-bc*Kio/2.0+0.5*HHi/HHo;


          double l11 = l11io*l11oj + l12io*l21oj + l13io*l31oj;
          double l12 = l11io*l12oj + l12io*l22oj + l13io*l32oj;
          double l13 = l11io*l13oj + l12io*l23oj + l13io*l33oj;
          double l21 = l21io*l11oj + l22io*l21oj + l23io*l31oj;
          double l22 = l21io*l12oj + l22io*l22oj + l23io*l32oj;
          double l23 = l21io*l13oj + l22io*l23oj + l23io*l33oj;
          double l31 = l31io*l11oj + l32io*l21oj + l33io*l31oj;
          double l32 = l31io*l12oj + l32io*l22oj + l33io*l32oj;
          double l33 = l31io*l13oj + l32io*l23oj + l33io*l33oj;
          
           
			/* debugging
          printf("%2.15f \n",l11); 
          printf("%2.15f \n",l12); 
          printf("%2.15f \n",l13); 
          printf("%2.15f \n",l21); 
          printf("%2.15f \n",l22); 
          printf("%2.15f \n",l23); 
          printf("%2.15f \n",l31); 
          printf("%2.15f \n",l32); 
          printf("%2.15f \n",l33); 
			 debugging*/
          /*raise(SIGINT);*/
      while (1)
	{
	  result_trr = read_trr(xd_read, natoms_trr, &step_trr, &time_trr, &lambda_trr,
				box_trr, x_trr, v_trr, f_trr);
	      
	  int ii_trr, jj_trr, x_ck=0, v_ck=0, f_ck=0;
	  int x_ck_bool=0, v_ck_bool=0, f_ck_bool=0;
          
 

            /*SPCE to TIP3P sart */
            for (ii_trr = 0; ii_trr < natoms_trr; ii_trr=ii_trr+4)
	    {

               for(jj_trr = 0; jj_trr < DIM; jj_trr++ )  
		{            
			/* debugging
		        if (ii_trr == 0 && time_trr ==10.0)
                        	{       printf("ii_trr %2d ,jj_trr %2d \n",ii_trr,jj_trr);
					printf("old %2.15f %2.15f %2.15f %2d\n",x_trr[ii_trr][jj_trr],x_trr[ii_trr+1][jj_trr],x_trr[ii_trr+2][jj_trr],jj_trr);
					printf("old %2.15f %2.15f %2.15f %2d\n",x_trr[ii_trr][jj_trr]*l11,x_trr[ii_trr+1][jj_trr]*l12,x_trr[ii_trr+2][jj_trr]*l13,jj_trr);
					printf("old %2.15f %2.15f %2.15f %2d\n",x_trr[ii_trr][jj_trr]*l21,x_trr[ii_trr+1][jj_trr]*l22,x_trr[ii_trr+2][jj_trr]*l23,jj_trr);
					printf("old %2.15f %2.15f %2.15f %2d\n",x_trr[ii_trr][jj_trr]*l31,x_trr[ii_trr+1][jj_trr]*l32,x_trr[ii_trr+2][jj_trr]*l33,jj_trr);
				}	
			*/
			double temp1 = x_trr[ii_trr][jj_trr];
                        double temp2 = x_trr[ii_trr+1][jj_trr];
			double temp3 = x_trr[ii_trr+2][jj_trr];

                        x_trr[ii_trr][jj_trr] = l11*temp1 + l12*temp2 + l13*temp3;
			x_trr[ii_trr+1][jj_trr] = l21*temp1 + l22*temp2 + l23*temp3;
			x_trr[ii_trr+2][jj_trr] = l31*temp1 + l32*temp2 + l33*temp3;
			x_trr[ii_trr+3][jj_trr] = (1.0-ai-bi)*x_trr[ii_trr][jj_trr] + ai*x_trr[ii_trr+1][jj_trr] + bi*x_trr[ii_trr+2][jj_trr];
			
			/* debugging
		        if (ii_trr == 0 && time_trr==10.0)
                        	{       printf("ii_trr %2d ,jj_trr %2d \n",ii_trr,jj_trr);
					printf("new %2.15f %2.15f %2.15f %2d\n",x_trr[ii_trr][jj_trr],x_trr[ii_trr+1][jj_trr],x_trr[ii_trr+2][jj_trr],jj_trr);
				}	
                        
			 debugging*/
		}	
            }	
	    
            /*printf(" %2.5f  %Le \n",time_trr,J);*/

	    for (ii_trr = 0; ii_trr < natoms_trr; ii_trr++)
	    {
	      for(jj_trr = 0; jj_trr < DIM; jj_trr++)
		{
		  if (x_trr[ii_trr][jj_trr] == 0)
		    x_ck++;
		  if (v_trr[ii_trr][jj_trr] == 0)
		    v_ck++;
		  if (f_trr[ii_trr][jj_trr] == 0)
		    f_ck++;
		}
	    }
	  
	  if (x_ck == natoms_trr*DIM)
	    x_ck_bool = 1;
	  if (v_ck == natoms_trr*DIM)
	    v_ck_bool = 1;
	  if (f_ck == natoms_trr*DIM)
	    f_ck_bool = 1;
	      
	  if (result_trr == 0) // if not reach the end of file, write it to the output.trr file
	    {
	      if (exdrOK != result_trr)
		die_r("Reading trr file",result_trr);
	      
	      if(v_ck_bool && f_ck_bool)
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, NULL, NULL);
	      else if(v_ck_bool)
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, NULL, f_trr);
	      else if(f_ck_bool)
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, v_trr, NULL);
	      else
		result_trr = write_trr(xd_write, natoms_trr, step_trr, time_trr, lambda_trr,
				       box_trr, x_trr, v_trr, f_trr);
	      
	      if (0 != result_trr)
		die_r("Writing trr file",result_trr);
	      
	    }
	  else
	    break;
	}
    }
  
  
  /* .trr -> .xtc */
  if(in_trrBool && out_xtcBool)
    {
      result_trr = read_trr_natoms(rfile, &natoms_trr);
      
      if (exdrOK != result_trr)
	die_r("read_trr_natoms",result_trr);
      
      x_trr = calloc(natoms_trr, sizeof(x_trr[0]));
      v_trr = calloc(natoms_trr, sizeof(v_trr[0]));
      f_trr = calloc(natoms_trr, sizeof(f_trr[0]));
      
      while(1)
	{
	  result_trr = read_trr(xd_read, natoms_trr, &step_trr, &time_trr, &lambda_trr,
				box_trr, x_trr, v_trr, f_trr);
	  
	  if (result_trr == 0) // if not reach the end of file, write it to the output.trr file
	    {
	      if (exdrOK != result_trr)
		die_r("Reading trr file", result_trr);
	      
	      result_xtc = write_xtc(xd_write, natoms_trr, step_trr, time_trr,
				     box_trr, x_trr, prec_xtc);
	      
	      if (result_xtc != 0)
		die_r("Writing xtc file", result_xtc);
		}
	  else
	    break;
	}
    }
  
  xdrfile_close(xd_read);
  xdrfile_close(xd_write);
  
}



int main(int argc, char *argv[])
{
  int inFileBool = 0;
  int outFileBool = 0;
  int in_xtcBool = 0;
  int out_xtcBool = 0;
  int in_trrBool = 0;
  int out_trrBool = 0;
  char *rfile=NULL, *wfile=NULL;
  double romi,rohi,angi,rohj,angj;
  

  if(argc != 15)
    {
      fprintf(stderr,"Usage: %s -i inFile -o outFile -romi romi_val -rohi rohi_val -angi angi_val -rohj -rohj_val -angj angj_val \n",argv[0]);
      exit(1);
    }
  else
    {
      int ii = 1;

      while(ii < argc)
	{
	  if(strcmp(argv[ii], "-i") == 0)         // if (argv[ii] == "-i")
	    {
	      ii++;
	      inFileBool = 1;
	 
	      if(strstr(argv[ii], ".xtc") != NULL)
		{
		  in_xtcBool = 1;
		  in_trrBool = 0;
		}
	      if(strstr(argv[ii], ".trr") != NULL)
		{
		  in_trrBool = 1;
		  in_xtcBool = 0;
		}

	      rfile = argv[ii];
	    }

	  if(strcmp(argv[ii], "-o") == 0)
	    {
	      ii++;
	      outFileBool = 1;
	 
	      if(strstr(argv[ii], ".xtc") != NULL)
		{		
		  out_xtcBool = 1;
		  out_trrBool = 0;
		}
	      if(strstr(argv[ii], ".trr") != NULL)
		{
		  out_trrBool = 1;
		  out_xtcBool = 0;
		}

	      wfile = argv[ii];
	    }

	  if(strcmp(argv[ii], "-romi") == 0)
	    {
	      ii++;
	      romi =  atof(argv[ii]);
              printf("romi = %f \n",romi);    
	    }

	  if(strcmp(argv[ii], "-rohi") == 0)
	    {
	      ii++;
	      rohi =  atof(argv[ii]);
              printf("rohi = %f \n",rohi);    
	    }
	  if(strcmp(argv[ii], "-angi") == 0)
	    {
	      ii++;
	      angi =  atof(argv[ii]);
              printf("angi = %f \n",angi);    
	    }
	  if(strcmp(argv[ii], "-rohj") == 0)
	    {
	      ii++;
	      rohj =  atof(argv[ii]);
              printf("rohj = %f \n",rohj);    
	    }

	  if(strcmp(argv[ii], "-angj") == 0)
	    {
	      ii++;
	      angj =  atof(argv[ii]);
              printf("angj = %f \n",angj);    
	    }
	  ii++;
	}
    }

  if(!inFileBool || !outFileBool)
    {
      perror("Usage : ./ReadWrite -i inFile -o outFile -romi romi_val -rohi rohi_val -angi angi_val -rohj -rohj_val -angj angj_val");
      exit(1);
    }

  ReadWrite(rfile, wfile, in_xtcBool, out_xtcBool, in_trrBool, out_trrBool,romi,rohi,angi,rohj,angj);
  

  return 0;
}
