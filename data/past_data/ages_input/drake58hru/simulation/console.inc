analysis(title:"drake58hru") {
	def f = "%last/Catchment.csv"
	def fn = "%last/Outlet.csv"
	def fnp = "%last/Catchment Nitrogen Pool.csv"
	def fcg = "%last/Catchment Crop Growth.csv"
	def fsed = "%last/Catchment Sediment.csv"

	timeseries(title:"Runoff") {
		x(file:fn, column:"time")
		y(file:fn, column:"catchmentSimRunoff")
	}

	timeseries(title:"Precip", view:STACKED) {
		x(file:f, column:"time")
		y(file:f, column:"precip")
		y(file:f, column:"netRainOut")
		y(file:f, column:"netSnowOut")
	}

	timeseries(title:"Snow") {
		x(file:f, column:"time")
		y(file:f, column:"snowTotSWE")
		y(file:f, column:"snowMelt")
		y(file:f, column:"snowDepth")
	}

	timeseries(title:"ET", view:STACKED) {
		x(file:f, column:"time")
		y(file:f, column:"potET")
		y(file:f, column:"actET")
		y(file:f, column:"rhum")
		y(file:f, column:"tmean")
	}

	timeseries(title:"Interception") {
		x(file:f, column:"time")
		y(file:f, column:"interception")
		y(file:f, column:"throughfall")
		y(file:f, column:"intercStorage")
	}

	timeseries(title:"Rad", view:STACKED) {
		x(file:f, column:"time")
		y(file:f, column:"solRad")
		y(file:f, column:"netRad")
	}

	timeseries(title:"Soil" , view:STACKED) {
		x(file:f, column:"time")
		y(file:f, column:"soilSatMPS")
		y(file:f, column:"soilSatLPS")
		y(file:f, column:"soilActMPS")
		y(file:f, column:"soilActLPS")
		y(file:f, column:"soilSat")
		y(file:f, column:"precip")
	}

	timeseries(title:"RD/RG", view:STACKED) {
		x(file:fn, column:"time")
		y(file:fn, column:"catchmentRD1_w")
		y(file:fn, column:"catchmentRD2_w")
		y(file:fn, column:"catchmentRG1_w")
		y(file:fn, column:"catchmentRG2_w")
		y(file:f, column:"percolation")
	}

	timeseries(title:"act RD/RG", view:STACKED) {
		x(file:f, column:"time")
		y(file:f, column:"actRG1")
		y(file:f, column:"actRG2")
	}

	timeseries(title:"N Pools" , view:STACKED) {
		x(file:fnp, column:"time")
		y(file:fnp, column:"sNO3_Pool")
		y(file:fnp, column:"sNH4_Pool")
		y(file:fnp, column:"sNResiduePool")
	}

	timeseries(title:"Biomass" , view:STACKED) {
		x(file:fcg, column:"time")
		y(file:fcg, column:"BioAct")
		y(file:fcg, column:"BioNAct")
		y(file:fcg, column:"FPHUact")
	}

	timeseries(title:"Stress" , view:STACKED) {
		x(file:fcg, column:"time")
		y(file:fcg, column:"nstrs")
		y(file:fcg, column:"wstrs")
		y(file:fcg, column:"tstrs")
	}

	timeseries(title:"Erosion" , view:STACKED) {
		x(file:fsed, column:"time")
		y(file:fsed, column:"gensed")
		y(file:fsed, column:"outsed")
		y(file:fn, column:"catchmentSed")
		y(file:fsed, column:"precip")
		y(file:fsed, column:"outRD1")
	}
}
