#include <string.h>
#include "syspara.h"

void vm_data(FILE *fp7, double time)
{

	int i,j;

	fprintf(fp7,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp7,"%lf ",var.u[0][i][j]);
		}
	}
	fprintf(fp7,"\n");

} 

void ina_data(FILE *fp8, double time)
{

	int i,j;
	
	fprintf(fp8,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp8,"%lf ",ina.total[i][j]);
		}
	}
	fprintf(fp8,"\n");
	
}

void ical_data(FILE *fp9, double time)
{

	int i,j;
	
	fprintf(fp9,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp9,"%lf ",ical.ica[i][j]);
		}
	}
	fprintf(fp9,"\n");
	 	
}

void cleft(FILE *fp10, double time)
{

	int i;

	fprintf(fp10,"%lf ",time);
	for(i=0;i<MEDIA_SITE-1;i++){
		fprintf(fp10,"%lf ",var.cleft_potential[i]);
	}
	fprintf(fp10,"\n");
	 	
}

void intra_v_data(FILE *fp11, double time)
{

	int i,j;

	fprintf(fp11,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp11,"%lf ",var.inner_v[i][j]);
		}
	}
	fprintf(fp11,"\n");

}

void nai_data(FILE *fp12, double time)
{

	int i,j;

	fprintf(fp12,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp12,"%e ",var.u[32][i][j]);
		}
	}
	fprintf(fp12,"\n");

}

void ki_data(FILE *fp13, double time)
{

	int i,j;

	fprintf(fp13,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp13,"%e ",var.u[34][i][j]);
		}
	}
	fprintf(fp13,"\n");

}

