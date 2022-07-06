parameter (file:"$data/drake_watershed_58hru.csv")
parameter (file:"$data/soils_hor_swc.csv")

parameter {
	params RAWPARAM

	// HRU and reach files
	hruFilePath          "$data/hrus.csv"    // HRU parameter file
	landuseFilePath      "$data/landuse.csv" // landuse parameter file
	hydroGeologyFilePath "$data/hgeo.csv"    // hydrogeology parameter file
	reachFilePath        "$data/reach.csv"   // stream reach parameter file

	// Cropping and management files
	fertilizerFilePath "$data/fert.csv"       // fertilizer operations file
	tillageFilePath    "$data/till.csv"       // tillage operations file
	cropFilePath       "$data/crop.csv"       // crop parameter file
	managementFilePath "$data/management.csv" // management operations file
	soilFilePath       "$data/soils_hor_swc.csv"  // soil properties parameter file
	routingFilePath    "$data/topo.csv"    // watershed topology (routing) file

	// Irrigation files
	//irrigationFilePath "$data/irrigation.csv" // irrigation file
	//managementIrrigationFilePath "$data/man_irri.csv" // irrigation ID (IID) and management ID (MID) linkage file

	// Crop rotation files
	rotationFilePath     "$data/croprotation.csv" // crop rotation parameter file
	hruRotationFilePath  "$data/hrus_rot.csv"     // HRU rotation parameter file

	// Climate files
	dataFilePrecip "$data/precip.csv"          // precipitation climate file
	dataFileTmin   "$data/tmin.csv"            // minimum temperature climate file
	dataFileTmax   "$data/tmax.csv"            // maximum temperature climate file
	dataFileHum    "$data/rhum.csv"            // humidity (absolute or relative) climate file
	dataFileSol    "$data/srad.csv"            // radiation (solar or sunlight hours) climate file
	dataFileWind   "$data/wind.csv"            // wind speed climate file
	dataFileCO2    "$data/co2.csv"             // CO2 parameter file
	dataFileReg    "$data/regionalization.csv" // regionalization parameter file

	humidity   "rel"  // absolute ("abs") or relative ("rel") humidity
	solar      "srad" // sunlight hours ("sunh") or solar radiation ("srad")

	// Additional parameter files
	// hruAdditionalFilePath      "$data/hrus_additional.csv"     // HRU additional parameter file
	// reachAdditionalFilePath    "$data/reach_additional.csv"    // Reach additional parameter file
	// temporalAdditionalFilePath "$data/temporal_additional.csv" // Temporal additional parameter file

	// HRU parameter override file
	//hruOverrideFilePath "$data/hru_override.csv" // HRU parameter over-ride file

	saveStateFilePath "$data/state.bin"
	loadStateFilePath "$data/state.bin"
}
