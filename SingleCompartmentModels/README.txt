//////////////////////  README for Schild 1994 Modeling Project  /////////////////////////

Author: David Catherall, Nicole Pelot
Last Updated: 2/22/2020

///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////   Overview    ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
This file contains a description of the file structure and instructions on how to run and modify
single compartment voltage clamp and current clamp simulations of the Schild 1994 and 1997 models.
the C-fiber Modeling project. Note that the mod files containing the ion channel mechanisms must 
be compiled using mknrndll after migrating the files to a different machine or after making any edits.
The resulting dll file must be stored in the directory in which the hoc files are stored. 

The membrane mechanisms are:

		Mechanism							Suffix
		
		Transient Calcium					cat
		Long-Lasting Calcium				can
		TTX-Sensitive Sodium				naf -- Schild 1994 dynamics
											naf97mean -- Schild 1997 dynamics (defaults to Schild 1994 max conductance)
		TTX-Insensitive Sodium				nas -- Schild 1994 dynamics
											nas97mean -- Schild 1997 dynamics (defaults to Schild 1994 max conductance)
		Delayer Rectifier Potassium			kd -- Note this is IK in Schild 1994, i.e. Kd in the publication
		Transient Potassium					ka -- Note this is IA in Schild 1994, one component of Ktrans in the publication
		Delay Potassium						kds -- Note this is ID in Schild 1994, one component of Ktrans in the publication
		Calcium Activated Potassium			kca
		Background Ca and Na				leak
		Na-Ca Exchanger						NaCaPump
		Calcium Pump						CaPump
		NaK Pump							NaKPump
		Internal Calcium Handling			caint -----------
		External Calcium Handling			caext			| See notes below on Ca handling
		Scaled Internal Calcium Handling	caintscale		| & on geometry
		Scaled External Calcium Handling	caextscale-------
		
///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Implementation Notes   /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

///////////// Form of the SS gating parameter equations ////////////////
The steady-state equation for each gating parameter is of the form:

	z_inf = 1.0/(1.0 * exp((V_half-V)/S_half)) = 1.0/(1.0 * exp((V-V_half)/(-S_half)))
	
However, in Schild 1997, while they defined the equation in this standard way (Equ. 1), 
the provided values for S_half (e.g., in Table 1 and in the caption for Fig 8) have the 
incorrect signs. S_half should have the same sign as the slope of z_inf vs. Vm.

///////////// Q10 factor and temperature ////////////////
We simulated the model at room and body temperatures. The text in the original publication
indicates that simulations were run at 22oC and 37oC Schild 1994, but in the table of constants
(Table 4), the temperature is listed as 296 K, which is 22.85oC. Therefore, we ran the room
temperature simulations at 22.85oC and used 22.85oC for the Schild 1994 Q10 reference temperature,
although the 0.85oC difference in room temperature values did not noticeably affect the results.
For the Schild 1997 sodium ion channels, we used a Q10 reference temperature of 22oC, taking the
mean of their stated range from 21 to 23oC.

Q10 factors provide temperature scaling for simulation parameters, specifically a multiplicative
factor for each 10oC change in temperature. Q10 values are provided in Schild without explanation 
of implementation methods. We assumed that they were used in the standard form with the following
equation:

τz(T)=τz(T0) * Q10z^(((T0-T)/10))

where z is a gating mechanism, T is the temperature in degrees Celsius, and T0 is the temperature 
at which τz(T0) was quantified. While conductance values are also temperature-dependent, Table 4 
of Schild 1994 only lists the gating variables as the targets for the Q10 scale factors. 
As discussed above, T0 was 296 K = 22.85oC for Schild 1994 and 22oC for Schild 1997. Note that
the Q10 value for the j gating variable for Naf is not listed in the table of constants in Schild 
1994 (Table 4); therefore, we did not apply Q10 temperature scaling to the tau equation for Naf,j. 
For the NaK pump, NaCa exchanger, and CaP pump, the provided Q10 values were used to scale
the I_NaK_bar, KNaCa, and I_CaP_bar values, respectively, rather than a time constant equation, as
indicated in Table 4 of Schild 1994. 

///////////// Faraday’s constant and gas constant ////////////////
Although F and R constants are built into Neuron, Schild 1994 used rounded values, especially for F,
F = 96500 C/mol and R = 8.314 J/(kg*mol*K), which had a significant effect on the output of the 
model. Therefore, we hard-coded the values for Faraday's constant and the gas constant in the Matlab
implementation (Schild_Int_Test.m, BackgroundVoltageClamps.m, CaWithHandling.m), in the Brian
implementation (SchildHolding.py & SchildRelease.py), in each HOC wrapper file 
(Schild_1994_A-C-Fiber_Model_Threshold.hoc, VclampModel.hoc), and in multiple mod files, when 
those constants are used in that mod file. To use these values for F and R, the ion_style function is 
used in the HOC scripts to let the reversal potentials be calculated manually. This was done for ENa 
and EK, but ECa is calculated in individual mod files, as the Schild equation is very different from 
the typical Nernst Equation. Therefore, amongst the mod files.... 
- F is hard-coded in all the calcium mechanisms except the Ca pump (caext, caextscale, caint,
  caintscale, can, cat, leak, and NaCapump)
- R is hard-coded in all the calcium current mechanisms except the Ca pump (can, cat, leak, and
  NaCaPump)

///////////// Calcium handling and model geometry ////////////////
The model assumes constant concentrations for the sodium and potassium ions, as well as equal
concentrations of each ion between the perineural and extracellular bath compartments; conversely,
the model includes calcium ion accumulation dynamics for the intracellular and perineural
compartments with constant extracellular concentration, as well as intracellular calcium buffering
and diffusion of calcium ions between the perineural space and extracellular bath. 

Note that in the NEURON implementation of the model, the perineural calcium concentration [Ca]s is 
called cao, and the extracellular/bath calcium concentration [Ca]o is a mechanism-specific parameter
called cabath. 

The perineural calcium concentrations generally at ~2 mM at rest. However, the intracellular calcium 
concentration at rest depends on the model tempereature and on the model dynamics & maximum 
conductances (Schild 1994 vs. 1997). See the publication supplemental material for specific values.

Regarding the diffusion of calcium ions between the perineural and extracellular spaces:
We examined the stability of the model in the absence of any inputs. We simulated the single 
compartment model in NEURON, and only inserted the calcium mechanisms. When the model was simulated
for 100 s without any stimulation input, the cell started to quickly depolarize at ~40 s. This was
due to an error in the original paper that caused the equation for the derivative of the perineural
calcium concentration ([Ca]s) to be unstable (Figure 14). The first term of the equation:

	d[Ca]s/dt = ([Ca]s-[Ca]o)/tau_Ca + ...
	
drives [Ca]s away from [Ca]o, instead of equilibrating the two. The equation should read:

	d[Ca]s/dt = ([Ca]o-[Ca]s)/tau_Ca + ...
	
This correction fixed the instability in the perineural calcium accumulation mechanism. This issue
was not evident in typical simulations given the long time constant of tau_Ca = 4511 ms.

Error in value for thickness of perineural space:
The text in (Schild et al., 1994) incorrectly states the perineural space thickness as 1.0 µm, but
the perineural volume (Vols) in Table 4 of the paper corresponds to a spherical shell with thickness
of 0.5 µm.

Single compartment geometry:
While the original model defined a spherical compartment (30 µm diameter), NEURON operates with
cylinders; therefore, we used a cylinder that was 20 µm in diameter and 45 µm in length to match the
volumes and surface areas, noting that the circular end-caps of a cylindrical compartment in NEURON
are not counted towards the surface area.

We originally coded the caint and caext mechanisms (i.e., intracellular and perineural calcium 
buffering/diffusion mechanisms, respectively) to replicate the Schild 1994 single compartment 
equations:

	surface area 	= pi * d * L 				
					= pi * (20 um) * (45 um)
					= 2.8274e-5 cm^2
	vol_intra    	= 0.9 * (pi*(d^2)/4) * L 	
					= 0.9 * (pi*((20 um)^2)/4) * (45 um) 	
					= 1.2723e-8 cm^3
	vol_peri     	= (4*pi/3) * (d_out^3 - d_in^3) 
					= (4*pi/3) * ((15.5 um)^3 - (15 um)^3) 
					= 1.4614e-9 cm^2
	
where the intracellular volume vol_intra) was 90% of the complete sphere (Schild 1994 paper) 
or cylinder (NEURON implementation) to account for the organelles taking up 10% of the 
intracellular volume. We calculated the volume of the perineural space (vol_peri) based on 
the spherical Schild geometry.

However, given that the surface area and volumes were hard-coded, they did not scale with model 
geometry. We modified the mechanisms, named caintscale and caextscale, to scale with geometry.
Note that in these versions, we excluded the 0.9 multiplicative factor for vol_intra given that
there are not many organelles inside the axon; this change did not affect the thresholds of the 
single compartment model. In these scaled mechanisms, we used a 0.5 um thickness for the 
cylindrical perineural space/shell. Note that this results in a different perineural volume as 
compared to the original single compartment model.

///////////////////////////////////////////////////////////////////////////////////	
///////////////////////////  Simulation Wrappers   ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
There are two HOC files which complete different simulation tasks. These are:

	Schild_1994_A-C-Fiber_Model.hoc
	VclampModel.hoc


///////////// Schild_1994_A-C-Fiber_Model.hoc ////////////////
	
The main .hoc file which puts all of the pieces together to run the single compartment
C-Fiber Model described in Schild 1994 is called Schild_1994_A-C-Fiber_Model.hoc. This file
is set up to run a current clamp test in order to observe action potential properties of
the model. The model can be run at varying temperatures, varying initial voltages, and 
also includes the ability to run the model as an A-Fiber as per the parameters for the 
A-Fiber model provided in Schild 1994. Stimulation can be initiated from resting voltage
or from a predefined holding voltage. In the latter case, a voltage clamp holds the cell at the 
holding voltage for the initialization time before the cell is stimulated. The current clamp  
can be delivered at a user-defined amplitude, or the thresholds can be found using a binary 
search. 

The parameters that should be edited to tailor simulations to a specific task are, in order of 
appearance:

	isafiber				If set to 1, the conductances and shift values for the Schild 1994
							A-type cell are inserted into the model. v_init is also set to -59 mV, 
							corresponding to the resting potential of an A-Fiber cell.
							
	insert97na				If set to 1, the Schild 1997 sodium channel dynamics are used (i.e., 
							naf97mean.mod and nas97mean.mod, instead of naf.mod and nas.mod)
							
	conductances97			If set to 1, the Schild 1997 maximum conductances are used.
	
	secondorder				NEURON variable which defines the integration method. Usually secondorder = 0,
							which corresponds to Backward Euler, works fine.

	stimfromrest 			Determines whether or not the cell is stimulated from rest or from 
							some holding voltage. This functionality was added because in Schild
							1994, many current clamps were performed after the cell was held at
							a hyperpolarized voltage.
							
	holdingv				Voltage which cell is held at if cell is being stimulated from some 
							other voltage besides resting potential (if stimfromrest == 0).
					
	simtime 				The duration of the simulation from t=0, after attaining steady-state. 
	
	inittime  				Length of time cell is initialized during t<0.
							If stimfromrest == 1, then the cell is initialized at v_init, then
							allowed to drift to rest.
							If stimfromrest == 0, then the cell is held at holdingv for the 
							duration of inittime.

	dt_initSS				Timestep of simulation during inittime. This is usually	much larger than 
							the simulation dt, in order to speed up the simulation. However, if it is 
							too large, then the numerical integration can diverge. 
							1 us suggested; 10 us can be too large.
							
	dt						Timestep of the simulation during t>0. This should be sufficiently small 
							such that a smaller dt would not change the value of the outcome measure 
							(e.g., threshold) within a chosen tolerance.
							
	celsius					Temperature at which simulation is run

	find_thresh				If find_thresh == 0, then a current pulse at CurrentClampAmplitude is delivered.
							If find_thresh == 1, then the binary search in FindThresh.hoc is used to find threshold.
							
	CurrentClampDuration	Duration of current clamp
	
	CurrentClampDelay		Delay from t = 0 for current clamp to start. Current clamp will start at 
							t = CurrentClampDelay.
	
	CurrentClampAmplitude	Amplitude of current clamp; only used if find_thresh == 0.
	
	stimamp_bottom_init 	Only used if find_thresh == 1.
							Defines the lower stimulation current bound for the binary search algorithm.
							If this is too high, such that it produces an action potential, an error
							will be thrown. Usually, 0 nA works fine for this parameter.
							
	stimpamp_top_init		Only used if find_thresh == 1.
							Defines the upper stimulation current bound for the binary search algorithm.
							If this is too low, such that it does not produce an action potential, an 
							error will be thrown. 

	thresh_resoln			Only used if find_thresh == 1.
							Defines the exit condition of the binary search, expressed as a fraction, such
							that when abs((stimamp_bottom - stimamp_top) / stimamp_top) < thresh_resoln, 
							the binary search is deemed to have found threshold.
							
	ap_thresh				The rising edge must pass ap_thresh to identify an action potential.
							
	
	filename				The filename for the file that is created to save simulation data. 
							If only a name is given, the file will be created in the same directory
							as the HOC code. Other folders can be created in the parent
							directory and the filename will then be: filename = "NewFolder/examplefile.dat"
							NEURON will not create NewFolder, it must already exist in the file
							system for the new file to be created. Files will be overwritten without
							warning if the filename refers to a file which already exists.						
							The user may want to modify what variables are recorded and saved to the file.
							
	The other HOC simulation structure is designed for the multicompartment simulations. However, to 
	modify this HOC code for that purpose...
		- diam = 1 (for a 1 um fiber).
		- L = 5000 (for a 5000 um long fiber).
		- nseg = 600 (to achieve 8.33 um segments with L = 5 mm; NOTE that this is equivalent to setting 
		L = 5000/600, and setting nseg = 1, and connecting the individual sections, rather than using a 
		single multi-segmented section).
		- Define and connect passive end compartments, if desired.
		- Ra = 100 (intracellular resistivity in ohm-cm; if not specified, NEURON will use its default value; 
		Ra isn't used for a single compartment simulations).
		- Make adjustments as needed for stimulus location, recording variables, and APCount.
		- Make sure to use caintscale/caextscale for multicompartment simulations, rather than caint/caext.
	
///////////  VClampModel.hoc  //////////////
							
VClamp_Model.hoc is a wrapper which is set up to run voltage clamps. The purpose of this setup
is to enable easy voltage clamp simulations for individual ionic current mechanisms. Usually, only
one mechanism at a time is inserted into the cell. Additionally, calcium buffering and diffusion 
mechanisms can be inserted in conjunction with a calcium current to see the interplay between those
two mechanisms. Many of the parameters, mainly related to the time of simulation, are similar to the 
parameters in Schild_1994_A-C-Fiber_Model.hoc. However, other relevant parameters which can be edited are:

	inittime				In this wrapper, inittime refers to the length of time the cell is held
							at the conditioning voltage. In many cases, this can be set to 0 if 
							v_init is set to the conditioning voltage level, as the state variables
							are initialized at their steady state value at this voltage, so holding
							the voltage at this level would only waste resources. The exception to 
							this is when calcium dynamics are inserted into the cell. In this case,
							calcium concentrations are initialized at some pre-defined value that is
							not the steady state value of these variables. It may be pertinent to then
							hold the cell at this conditioning voltage for some time before the test
							voltage is applied.
							
	clampinit				Defines the conditioning level of the voltage clamp. Default value is
							whatever v_init is set to, so this can be edited by changing v_init.
							However, there might be a case where the conditioning voltage needs to
							be some other value
							
	ClampAmp				This is the amplitude of the test level for the voltage clamp
	
	ClampDur				This is the duration of the test level for the voltage clamp
	
	cadynamic				This variable defines whether or not the calcium is treated as constant
							or dynamic. If cadynamic is set to 1, caint and caext should be inserted
							into the cell.
							
	napresent				This should be set to 1 if a sodium mechanism is inserted into the cell.
							If a sodium mechanism is inserted, the sodium concentrations are inserted
							into the cell. This is included because if there is no sodium mechanism
							inserted, and you try to define a sodium concentration, NEURON will throw an error.
							
	kpresent				Analogous to previous, but when a potassium mechanism is inserted
	
	capresent				Analogous to previous, but when a calcium mechanism is inserted
	

There is a separate folder called \mod_Vclamp, which are mod files which can be used to reproduce
voltage clamp plots found in Schild 1994 in combination with VClamp_Model.hoc. These mod files have
different conductances than are found in the C-Fiber or A-Fiber models outlined in Schild 1994. These 
files should not be used to run any simulations in current clamp. To use these files, they must be
compiled separately from the files in \mod, and the nrnmech.dll that is created must be moved into
the parent directory.
	
	Note: Voltage Clamp mod files in \mod_Vclamp have the same suffixes as the files in \mod, except 
	with an "a" prepended, e.g. anaf, akd, and acat

At the end of all wrapper files, one or two different files are created to store data. If there are two 
files, by default, one of them is commented out. These files save variables which were useful for
verification purposes when the model was being created, however they can be changed in order to save any
data that is relevant to the simulation. These lines of code should be altered to reflect whatever data
needs to be saved.
	
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Matlab Scripts  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

Also included in the file tree are two Matlab scripts which run Voltage Clamps for the Schild
1994 model. They have no functionality in the NEURON model, but are useful for verification
of model function.

	Schild_Int_Test.m				Calculate Vclamp data for primary currents analytically and numerically
									with constant calcium concentrations; run Vclamps for C-type cell and for 
									reproductions of the figures in Schild 1994 (just comment/uncomment conductance and shift values)
	CaWithHandling.m				Calculate Vclamp data for calcium currents including dynamic calcium concentrations

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Brian Scripts  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

Additionally, the Brian model of Schild 1994 used to validate the NEURON model is also 
included in the files. There are two files, Schild1994HOLD.py and Schild1994Release.py. 
The HOLD script holds the cell at a particular voltage and saves the states at the end 
so that the Release script can use those states as its initialization. To run the Brian 
model, you will need an installation of Python3 on your machine. Check the Brian 
documentation for installation recommendations. (http://brian2.readthedocs.io/en/stable/index.html)