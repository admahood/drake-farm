resource "$oms_prj/dist"   // dll/so files (JNA)

// define output strategy: output base dir and
// the strategy NUMBERED|SIMPLE|DATE
outputstrategy(dir:"$data_prj/output", scheme:SIMPLE)

//  define the simulation model
model(classname:"ages.AgES") {
	#include "$sim/input.inc"
	#include "$sim/program.inc"
	#include "$sim/flags.inc"
	#include "$sim/output.inc"
	#include "$sim/simulation.inc"
}
#include "$sim/console.inc"
