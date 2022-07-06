parameter {
	// idset_hru is a ";"" separator list, e.g., "435;441;476" will select HRUs with IDs 435, 441, and 476
	// leave blank (i.e., "") to use all HRU IDs
	//idSet_hru ""
    idSet_hru "13;17;42;44;18;33;43;45;65;63;60;54;64;62;59;52" 	
	// idset_reach uses actual (negative) reach IDs, i.e., "-87;-91;-95" will select reaches with IDs -87, -91, and -95
	// leave blank (i.e., "") to use all stream reach IDs
	idSet_reach ""
	
	outFile_hru                "-"
	attrSet_hru                "ID;precip;tmean;rhum;wind;rain;snow;solRad;netRad;netRainOut;netSnowOut;irrigation_amount;interception;intercStorage;throughfall;snowDepth;snowTotSWE;snowMelt;potentialET;potET;actET;soilSat;soilSatMPS;soilSatLPS;soilActMPS;soilActLPS;percolation;inRD1;inRD2;inRG1;inRG2;actRG1;actRG2;outRD1;outRD2;outRG1;outRG2";
	attrSet_hru_w              "rain;snow;potentialET;potET;actET;netRainOut;netSnowOut;throughfall;irrigation_amount;percolation;interception;intercStorage;snowTotSWE;snowMelt;soilActMPS;soilActLPS;inRD1;inRD2;inRG1;inRG2;actRG1;actRG2;outRD1;outRD2;outRG1;outRG2"
	
	outFile_hru_crop           "-"
	attrSet_hru_crop           "ID;zrootd;BioAct;BioNAct;FPHUact;PHUact;nstrs;wstrs;tstrs;nmin;NYield;BioYield;Addresidue_pool;Addresidue_pooln;cropid;nfert"
	attrSet_hru_crop_w         "-"
	
	outFile_hru_layer          "HRU_Layered.csv"
	attrSet_hru_layer          "ID;swc_h[0];swc_h[0];swc_h[1];swc_h[2];swc_h[3];swc_h[4];swc_h[5];swc_h[6];swc_h[7];swc_h[8];swc_h[9];swc_h[10];swc_h[11];swc_h[12];swc_h[13];swc_h[14];swc_h[15];swc_h[16];swc_h[17];swc_h[18];swc_h[19];outRD2_h;aEP_h;NO3_N;InterflowN"
	attrSet_hru_layer_w        "outRD2_h;aEP_h"
	
	outFile_hru_n_mb           "-"
	attrSet_hru_n_mb           "ID;SurfaceN_in;SurfaceNabs;sinterflowN_in;sinterflowNabs;PercoNabs;Denit_trans;Volati_trans;actnup;SurfaceN;sinterflowN;sum_Ninput;PercoN"
	attrSet_hru_n_mb_w         "-"
	
	outFile_hru_n_pool         "-"
	attrSet_hru_n_pool         "ID;sNO3_Pool;sNH4_Pool;sN_activ_pool;sN_stable_pool;sNResiduePool"
	attrSet_hru_n_pool_w       "-"
	
	outFile_hru_sediment       "-"
	attrSet_hru_sediment       "ID;precip;tmean;outRD1;insed;gensed;outsed;sedpool"
	attrSet_hru_sediment_w     "outRD1;insed;gensed;outsed;sedpool"
	
	outFile_hru_crop_upgm      "-"
	attrSet_hru_crop_upgm      "ID;Addresidue_pool;Addresidue_pooln;BioYield;LAI;PHUact;BioAct;BioagAct;BioOpt_delta;deltabiomass;CanHeightAct;zrootd;HarvIndex;FPHUact;wstrs;tstrs;standleaf;standstem;standstore;flatleaf;flatstem;flatstore;grainf;root"
	attrSet_hru_crop_upgm_w    "-"
	
	outFile_catch              "-"
	attrSet_catch              "precip;tmean;rhum;wind;rain;snow;solRad;netRad;netRainOut;netSnowOut;irrigation_amount;interception;intercStorage;throughfall;snowDepth;snowTotSWE;snowMelt;potentialET;potET;actET;soilSat;soilSatMPS;soilSatLPS;soilActMPS;soilActLPS;percolation;inRD1;inRD2;inRG1;inRG2;actRG1;actRG2;outRD1;outRD2;outRG1;outRG2"
	attrSet_catch_w            "precip;tmean;rhum;wind;rain;snow;solRad;netRad;netRainOut;netSnowOut;irrigation_amount;interception;intercStorage;throughfall;snowDepth;snowTotSWE;snowMelt;potentialET;potET;actET;soilSatMPS;soilSatLPS;soilActMPS;soilActLPS;percolation;inRD1;inRD2;inRG1;inRG2;actRG1;actRG2;outRD1;outRD2;outRG1;outRG2"
	
	outFile_catch_crop         "-"
	attrSet_catch_crop         "zrootd;BioAct;BioNAct;FPHUact;PHUact;nstrs;wstrs;tstrs;nmin;NYield;BioYield;Addresidue_pool;Addresidue_pooln;nfert"
	attrSet_catch_crop_w       "zrootd;BioAct;FPHUact;PHUact;nmin;NYield;BioYield;nfert"
	
	outFile_catch_n_mb         "-"
	attrSet_catch_n_mb         "SurfaceN_in;SurfaceNabs;sinterflowN_in;sinterflowNabs;PercoNabs;Denit_trans;Volati_trans;actnup;SurfaceN;sinterflowN;sum_Ninput;PercoN"
	attrSet_catch_n_mb_w       "-"
	
	outFile_catch_n_pool       "-"
	attrSet_catch_n_pool       "sNO3_Pool;sNH4_Pool;sN_activ_pool;sN_stable_pool;sNResiduePool"
	attrSet_catch_n_pool_w     "-"
	
	outFile_outlet             "-"
	attrSet_outlet             "catchmentSimRunoff;catchmentSimRunoffN;catchmentSimRunoff_NO3;catchmentSimRunoff_NO3_N;catchmentNO3_N_Load_kg_ha;catchmentRD1_w;catchmentRD2_w;catchmentRG1_w;catchmentRG2_w;catchmentNRD1_w;catchmentNRD2_w;catchmentNRG1_w;catchmentNRG2_w;catchmentSed;deepsinkW;deepsinkN"
	
	outFile_catch_sediment     "-"
	attrSet_catch_sediment     "precip;tmean;outRD1;insed;gensed;outsed;sedpool"
	attrSet_catch_sediment_w   "precip;tmean;outRD1;insed;gensed;outsed;sedpool"
	
	outFile_catch_crop_upgm    "-"
	attrSet_catch_crop_upgm    "BioNAct"
	attrSet_catch_crop_upgm_w  "-"
	
	outFile_reach              "-"
	attrSet_reach              "ID;simRunoff;actRD1;actRD2;actRG1;actRG2;outRD1;outRD2;outRG1;outRG2"
}
