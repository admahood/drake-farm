calibration_start "2003-04-01"         // Calibration start date
    start_month_of_year 1
    rounds 1                               // calibration rounds, default 1
    summary_file "drake_luca_data.txt"
    trace_file "trace_of.txt"  
   // model_stdout false
    // step definitions
   
     step()
        {  
        
        max_exec(1000)
        
        parameter {

           soilPolRed(lower:1, upper:2, calib_strategy:INDIVIDUAL)
           BetaW(lower:8, upper:25, calib_strategy:INDIVIDUAL)
		soilMaxDPS(lower:0, upper:2, calib_strategy:INDIVIDUAL)
           soilMaxInfSummer(lower:0, upper:90, calib_strategy:INDIVIDUAL)
           soilMaxInfWinter(lower:0, upper:80, calib_strategy:INDIVIDUAL)
           soilMaxInfSnow(lower:0, upper:80, calib_strategy:INDIVIDUAL)
           soilDistMPSLPS(lower:2, upper:4, calib_strategy:INDIVIDUAL)	
		soilDiffMPSLPS(lower:3, upper:10, calib_strategy:INDIVIDUAL)
           soilOutLPS(lower:3, upper:10, calib_strategy:INDIVIDUAL)
           kdiff_layer(lower:0.05, upper:0.2, calib_strategy:INDIVIDUAL)
		soilLatVertLPS(lower:0.1, upper:3, calib_strategy:INDIVIDUAL)
		kf(lower:0, upper:200, calib_strategy:MEAN, filter_param:"soils\$SID", subset:"10101-10120,10201-10220,10401-10420,10501-10520,10601-10620,10801-10820,10901-10920,11001-11020,11201-11220,11301-11320,11501-11520,11701-11720,11801-11820,12001-12020,12101-12120,12201-12220,12301-12320,12401-12420,12501-12520,12601-12620,12701-12720,12801-12820,12901-12920,13001-13020,13101-13120,13201-13220,13301-13320,13401-13420,13501-13520,13601-13620,13701-13720,13801-13820,13901-13920,14001-14020,14101-14120,14201-14220,14301-14320,14401-14420,14501-14520,14601-14620,14801-14820,14901-14920,15001-15020,15101-15120,15201-15220,15301-15320,15401-15420,15501-15520,15601-15620,15701-15720,15801-15820,15901-15920,16001-16020,16101-16120,16201-16220,16301-16320,16401-16420,16501-16520")
           
                  }
        
        //30cm
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_13.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_a1.csv", table:"A1",column:"swc[0]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_17.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_a2.csv", table:"A2",column:"swc[0]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_42.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_a3.csv", table:"A3",column:"swc[0]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_44.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_a4.csv", table:"A4",column:"swc[0]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_65.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_c1.csv", table:"C1",column:"swc[0]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_63.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_c2.csv", table:"C2",column:"swc[0]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_60.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_c3.csv", table:"C3",column:"swc[0]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_54.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_c4.csv", table:"C4",column:"swc[0]")
                                                }
 
      objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_18.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_b1.csv", table:"B1",column:"swc[0]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_33.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_b2.csv", table:"B2",column:"swc[0]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_43.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_b3.csv", table:"B3",column:"swc[0]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_45.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_b4.csv", table:"B4",column:"swc[0]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_64.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_d1.csv", table:"D1",column:"swc[0]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_62.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_d2.csv", table:"D2",column:"swc[0]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_59.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_d3.csv", table:"D3",column:"swc[0]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_52.csv", table:"output", column:"swc_h[3]")
          observed(file:"$data/oswc_d4.csv", table:"D4",column:"swc[0]")
                                                }
        //60cm
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_13.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_a1.csv", table:"A1",column:"swc[1]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_17.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_a2.csv", table:"A2",column:"swc[1]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_42.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_a3.csv", table:"A3",column:"swc[1]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_44.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_a4.csv", table:"A4",column:"swc[1]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_65.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_c1.csv", table:"C1",column:"swc[1]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_63.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_c2.csv", table:"C2",column:"swc[1]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_60.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_c3.csv", table:"C3",column:"swc[1]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_54.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_c4.csv", table:"C4",column:"swc[1]")
                                                }
 
      objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_18.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_b1.csv", table:"B1",column:"swc[1]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_33.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_b2.csv", table:"B2",column:"swc[1]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_43.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_b3.csv", table:"B3",column:"swc[1]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_45.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_b4.csv", table:"B4",column:"swc[1]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_64.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_d1.csv", table:"D1",column:"swc[1]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_62.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_d2.csv", table:"D2",column:"swc[1]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_59.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_d3.csv", table:"D3",column:"swc[1]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_52.csv", table:"output", column:"swc_h[6]")
          observed(file:"$data/oswc_d4.csv", table:"D4",column:"swc[1]")
                                                }
        //90cm
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_13.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_a1.csv", table:"A1",column:"swc[2]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_17.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_a2.csv", table:"A2",column:"swc[2]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_42.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_a3.csv", table:"A3",column:"swc[2]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_44.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_a4.csv", table:"A4",column:"swc[2]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_65.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_c1.csv", table:"C1",column:"swc[2]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_63.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_c2.csv", table:"C2",column:"swc[2]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_60.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_c3.csv", table:"C3",column:"swc[2]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_54.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_c4.csv", table:"C4",column:"swc[2]")
                                                }
 
      objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_18.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_b1.csv", table:"B1",column:"swc[2]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_33.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_b2.csv", table:"B2",column:"swc[2]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_43.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_b3.csv", table:"B3",column:"swc[2]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_45.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_b4.csv", table:"B4",column:"swc[2]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_64.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_d1.csv", table:"D1",column:"swc[2]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_62.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_d2.csv", table:"D2",column:"swc[2]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_59.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_d3.csv", table:"D3",column:"swc[2]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_52.csv", table:"output", column:"swc_h[9]")
          observed(file:"$data/oswc_d4.csv", table:"D4",column:"swc[2]")
                                                }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
        //120cm or 150cm
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_13.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_a1.csv", table:"A1",column:"swc[3]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_17.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_a2.csv", table:"A2",column:"swc[3]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_42.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_a3.csv", table:"A3",column:"swc[3]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_44.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_a4.csv", table:"A4",column:"swc[3]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_65.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_c1.csv", table:"C1",column:"swc[3]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_63.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_c2.csv", table:"C2",column:"swc[3]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_60.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_c3.csv", table:"C3",column:"swc[3]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_54.csv", table:"output", column:"swc_h[12]")
          observed(file:"$data/oswc_c4.csv", table:"C4",column:"swc[3]")
                                                }
 
      objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_18.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_b1.csv", table:"B1",column:"swc[3]")
                                                }    
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_33.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_b2.csv", table:"B2",column:"swc[3]")
                                                }                                              
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_43.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_b3.csv", table:"B3",column:"swc[3]")
                                                }
        
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_45.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_b4.csv", table:"B4",column:"swc[3]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_64.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_d1.csv", table:"D1",column:"swc[3]")
                                                }
                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_62.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_d2.csv", table:"D2",column:"swc[3]")
                                                }                                        
         
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_59.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_d3.csv", table:"D3",column:"swc[3]")
                                                }                                                                                                                        
                                                                                                                                
        objfunc(method:KGE, timestep:DAILY) {
          simulated(file:"HRU_Layered_52.csv", table:"output", column:"swc_h[15]")
          observed(file:"$data/oswc_d4.csv", table:"D4",column:"swc[3]")
                                                }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
        }


        
