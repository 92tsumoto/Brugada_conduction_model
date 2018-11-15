#include "syspara.h"

typedef double Number;
typedef long long Lint;

void initial_mem()
{
	int i,k,m,l,z;

// initialized memory for state variables and arguments for Eular method 
        var.u=(double***)malloc(sizeof(double**)*NN);
        var.k1=(double***)malloc(sizeof(double**)*NN);
        if(var.u==NULL || var.k1==NULL) exit(1);
        for(i=0;i<NN;i++){
                var.u[i] = (double**)malloc(sizeof(double*)*MEDIA_SITE);
                var.k1[i] = (double**)malloc(sizeof(double*)*MEDIA_SITE);
                if(var.u[i]==NULL || var.k1[i]==NULL) exit(1);
                for(k=0;k<MEDIA_SITE;k++){
                        var.u[i][k] = (double*)malloc(sizeof(double)*MEDIA_PATCH);
                        var.k1[i][k] = (double*)malloc(sizeof(double)*MEDIA_PATCH);
                        if(var.u[i][k]==NULL || var.k1[i][k]==NULL) exit(1);
                }
        }

// initialized tablization memorys for Exp functions

	// ina_fast
	ina.Tmss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttaum=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauh=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tjss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauj=(Number *)calloc(VNMAX,sizeof(Number));
	ina.ThCaMKss=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmss==NULL || ina.Ttaum==NULL 
		|| ina.Thss==NULL || ina.Ttauh==NULL 
		|| ina.Tjss==NULL || ina.Ttauj==NULL 
		|| ina.ThCaMKss==NULL ) exit(1);
	// ina_late
	ina.Tmlss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauml=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thlss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.ThlCaMKss=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmlss==NULL || ina.Ttauml==NULL
		|| ina.Thlss==NULL || ina.ThlCaMKss==NULL ) exit(1);

	// ito
	ito.Tass=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaua=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Tiss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaui_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaui_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TAi_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TaCaMKss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TdeltaCaMK_dev=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TdeltaCaMK_rec=(Number *)calloc(VNMAX,sizeof(Number));
	if(var.celltype == 3){
		ito.Tdepi=(Number *)calloc(VNMAX,sizeof(Number));
		if(ito.Tdepi == NULL) exit(1);
	}
	if( ito.Tass==NULL || ito.Ttaua==NULL 
		|| ito.Tiss==NULL || ito.Ttaui_fast==NULL || ito.Ttaui_slow==NULL
		|| ito.TAi_fast==NULL || ito.TaCaMKss==NULL
		|| ito.TdeltaCaMK_dev==NULL || ito.TdeltaCaMK_rec==NULL ) exit(1);
	
	// ical
	ical.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaufca_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaufca_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ical.TAfca_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Texp_Ca=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Texp_Na=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Texp_K=(Number *)calloc(VNMAX,sizeof(Number));
	if( ical.Tdss==NULL || ical.Ttaud==NULL 
		|| ical.Tfss==NULL || ical.Ttauf_fast==NULL || ical.Ttauf_slow==NULL
		|| ical.Ttaufca_fast == NULL || ical.Ttaufca_slow==NULL || ical.TAfca_fast==NULL 
		|| ical.Texp_Ca==NULL || ical.Texp_Na==NULL || ical.Texp_K == NULL ) exit(1);

	// ikr
	ikr.Txrss=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.TAxr_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Trkr=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikr.Txrss==NULL || ikr.Ttauxr_fast==NULL || ikr.Ttauxr_slow==NULL 
		|| ikr.TAxr_fast==NULL || ikr.Trkr==NULL ) exit(1);

	// iks
	iks.Txs1ss=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs1=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs2=(Number *)calloc(VNMAX,sizeof(Number));
	if( iks.Txs1ss==NULL || iks.Ttauxs1==NULL || iks.Ttauxs2==NULL ) exit(1);

	// ik1
	ik1.Tk1ss=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Ttauk1=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Trk1=(Number *)calloc(VNMAX,sizeof(Number));
	if( ik1.Tk1ss == NULL || ik1.Ttauk1==NULL || ik1.Trk1==NULL ) exit(1);

	// inaca
	var.Thca=(Number *)calloc(VNMAX,sizeof(Number));
	var.Thna=(Number *)calloc(VNMAX,sizeof(Number));
	if( var.Thca==NULL || var.Thna==NULL ) exit(1);

	// inak
	inak.Tknai=(Number *)calloc(VNMAX,sizeof(Number));
	inak.Tknao=(Number *)calloc(VNMAX,sizeof(Number));
	if( inak.Tknai==NULL || inak.Tknao==NULL ) exit(1);

	// ikb	
	ikb.Txkb=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikb.Txkb==NULL ) exit(1);

	// icab
	icab.Texp=(Number *)calloc(VNMAX,sizeof(Number));
	if( icab.Texp==NULL ) exit(1);

	// inab
	inab.Texp=(Number *)calloc(VNMAX,sizeof(Number));
	if( inab.Texp==NULL ) exit(1);

}
		
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
void close_mem()
{

	int i,j;

		for(i=0;i<NN;i++){
			for(j=0;j<MEDIA_SITE;j++){
				free(var.u[i][j]);
				free(var.k1[i][j]);
			}
			free(var.u[i]);
			free(var.k1[i]);
		}
		free(var.u); free(var.k1);

	// ina_fast
		free(ina.Tmss); free(ina.Ttaum); 
		free(ina.Thss); free(ina.Ttauh); 
		free(ina.Tjss); free(ina.Ttauj);
		free(ina.ThCaMKss);

	// ina_late
		free(ina.Tmlss); free(ina.Ttauml); free(ina.Thlss); free(ina.ThlCaMKss);

	// ito
		free(ito.Tass); free(ito.Ttaua); free(ito.Tiss); free(ito.Ttaui_fast); free(ito.Ttaui_slow);
		free(ito.TAi_fast); free(ito.TaCaMKss); if(var.celltype==3){free(ito.Tdepi);}
		free(ito.TdeltaCaMK_dev); free(ito.TdeltaCaMK_rec);

	// ical
		free(ical.Tdss); free(ical.Ttaud); free(ical.Tfss); free(ical.Ttauf_fast); free(ical.Ttauf_slow);
		free(ical.Ttaufca_fast); free(ical.Ttaufca_slow); free(ical.TAfca_fast);
		free(ical.Texp_Ca); free(ical.Texp_Na); free(ical.Texp_K);

	// ikr
		free(ikr.Txrss); free(ikr.Ttauxr_fast); free(ikr.Ttauxr_slow); 
		free(ikr.TAxr_fast); free(ikr.Trkr);

	// iks
		free(iks.Txs1ss); free(iks.Ttauxs1); free(iks.Ttauxs2);
	
	// ik1
		free(ik1.Tk1ss); free(ik1.Ttauk1); free(ik1.Trk1);

	// inaca
		free(var.Thca); free(var.Thna);

	// inak
		free(inak.Tknai); free(inak.Tknao);

	// ikb	
		free(ikb.Txkb);

	// icab	
		free(icab.Texp);

	// inab	
		free(inab.Texp);

}

