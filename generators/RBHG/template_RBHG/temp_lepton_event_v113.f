c
c	To compile and run, in the directory just above this,  "python mac_RBHG.py”
c	Go to "output" directory to find the 3-vector file
c
c	version 6 of code weights the distribution of events by xs= acos( costheta ).  This is the most efficient version
c	of the code to run.  The nuclear and atomic form factors are included. 
c	version 8 allows for normal event weighting, or non-standard event weighting, theta=x**2, switchable with ‘standard_phase_space’. 
c		Non-standard is the most efficient. 
c	version 9: histograms J_T phi angle, uses simplified version of W_pol. 
c	version 10: histograms cross section arrays
c	version 11: histogram integrated cross section as a function of x and phi1. 
c	version 12: introduced logicals to do electron, muon, and muon with pion mass hypothesis histogramming. 
c	version 13: added some histogram options and controls
c	version 14: added some more controls, got rid of things in the code that weren’t useful 
c	version 15: brought over some features from the pi+pi- code, fixed some bugs.   Can control shape of the proton form factor. 
c	version 16: consistently use transverse momentum transfer squared everywhere, changed name to lepton_event
c	Uses rejection sampling: 
c	version 17: added feature to run one kinematic point
c	version 17_1: use weighting dcos theta/dx = (1-cos theta)**N , got rid of the option for bending the FF
c	version 17_2: use weighting theta = x**phase_space, with phase_space an integer >1 .   Set phase_space = 0 for dcos theta/dx =1. 
c		      makes a guess for the largest cross section, which seems to be correct
c	version 17_3: add option for log t plot, control proton rms radius independent of dipole FF parameter
c	version 17_4: reads Brem. file, independent histogram output control, reworked hist. binning, 
c		      cut on w distribution. This is the preferred version. 
c	version 17_5: added many options for plot output, added control for selecting transverse q2, or 3-momentum q2 dependence in cross section
c	version 17_6x1: allow for energy conservation, energy is transferred to the recoiling proton, got rid of "analysis" subroutine 
c	version 17_6x2: Energy is conserved. Define x as E1-m_part=x*(E1+E2-2*m_part), where E1 and E2 are the energies of the leptons 
c	version 17_6x3: Standard formulation of the program: 
c			Went back to old definition of x, E1=x*E0, but using energy conservation for second lepton.  Added option for doing Heitler 
c			formulation of cross section, with or w/O polarization. 
c			Eliminated options for non-energy conservation, and don't use q2_T in the denominator of the Bakmaev formula. 
c			For Bakmaev require that transverse momenta of the leptons should be at least 10*m_lepton. I've only 
c			tested this for e+e-, this will probably have to be changed for mu+mu-. The stated requirement is pT**2 >> m_lepton**2.
c			Removed option for W cut, not clear why we would want to do this. 
c			Added feature for para or perp polarization 
c	version 100:	Added internal radiation, added histograms for elasticity and missing mass squared
c	version 101: 	Randomized lepton 1 and lepton 2 on the output 
C	version 102:	Setup logicals for Bakmaev, Heitler, or Berlin options. Allowed phi_JT histogramming to be for either Bakmaev or Berlin
C	version 103:	Seem to have fixed rare occurrences of square-roots of negative numbers, producing NaN in the output file. 
C	version 104:	Variable t-bin width according to Andrew's prescription
C	version 105:	Allow for polarization at 45 deg or 135 deg
C	version 106:	Use logical flag GlueX=.true. for GlueX configuration, GlueX=.false. for CPP configuration. Use Fourier transform of hydrogen electronic
C			charge density for the atomic FF. The formula used in v105 of the code and earlier, from Bakmaev, is wrong. Cleaned up the code in many places,
C			and reverified the 208Pb FF.  Simplified the t-distribution output, less confusing. Increased the # of output digits in the txt event file. 
C	version 107:	Use q**2 in place of Q2 consistently in the Heitler and Berlin forms of the cross section, where q is the 3-momentum transfer. Output 3-momentum
C			vector for Brem. All form factors	(nuclear and atomic) are now a function of q. 
C	version 108:	Cleaned up the code in several places: (1) the emitted Brem photon is now emitted in the direction of either lepton #1 or #2. (2) use
C			energy conservation to find the energy of the recoil nucleus, and momentum conservation to find the direction vector of the recoil nucleus, 
C			(3) fixed a bug for rotating nucleus and Brem. momentum vectors into the 45/135 degree plane, (4) check energy and momentum conservation for 
C			the event. Tau
c	version 110:	Allow both leptons to radiate.  Went back to the original formula from Mo and Tsai for calculating the internal radiation length.
C			Modify proton dipole FF parameter so it is consistent with PDG value for proton RMS charge radius
C	version 111:	Introduced logical switch that allows control over one or two leptons radiating
C	version 112:	Calculate the energy/momentum for a ghost particle that when summed with the energy/momentum of the incident photon, two leptons and recoil
C			target gives energy/momentum conservation.   Output energy/momentum for the ghost particle in the ascii file, no longer output momentum 
C			for the two Brem. gamma rays.  
C	version 113:	Output 4-vector for leptons and recoil nucleus, got rid of ghost particle
C			
c
	implicit none
	real*8 E0,theta_min,theta_max
	real*8 pi,cos_max,cos_min,cross_sum,cross_max,cross_max_old,delta_cross_max,
     &	cross_test,x_min,x_max,theta1,phi1,costheta1,theta2,phi2,
     &	costheta2,cross_section,total_xscn,total_xscn_old,delta_total_xscn,
     &	xsctn,k1(3),k2(3),m_e,m_part,m_muon,pol,x1,x2,Q2,q,failure,			
     &	w_mumu,xs1,xs2,xs_max,xs_min,
     &	Atgt(100),tgtlen(100),radlen(100),c_den(100),a_den(100),rho0,hbarc,W,jacobian,
     &	m_pi,phi_JT,delta_w,delta_x,delta_phi,mtgt,ktgt(3),Etgt,kghost(3),Eghost,Ptgt,Ebrem1(3),Ebrem2(3),
     &	missing_mass,t,delta_t,x_value,delta_theta,
     &	delta_elasticity,delta_mm,
     &	Rexp,xmax,theta1_max,theta2_max,phi1_max,phi2_max,temp,frac_delta_t,delta_Egamma,
     &	data_array_w(1:200),data_array_x(1:200),data_array_t(1:200),data_array_phi_JT(1:200),data_array_t_variable(1:200),
     &	data_array_Egamma(1:200),data_array_theta(1:200),data_array_q2_ratio(1:200),data_array_elasticity(200),
     &	data_array_mm(200),delta_k(3),Econs, 
     &	error_array(1:200),Egamma_max,E_hi,E_lo,Brem,E_coherent,Egamma,bin_width,
     &	x_value_lo,x_value_hi,delta_x_value,delta_q2_ratio,eta_rad,b_rad,Eloss1,Eloss2,p1,p2,E1,E2,JT(2),
     &	test,E_rad,ks(3),elasticity,phiR,p1_tran,p2_tran,c1,c2,vec(3),sqrt2,count,delta
	integer*4 iseed
	integer*8 i,itest(4),nevent,j,nfail,bad_max,i_array,j_array, ztgt,phase_space
	logical hist_w,hist_x,hist_t,hist_phi_JT,hist_Egamma,hist_theta,output_event,hist_t_variable,
     &	nuc_FF,muon,electron,brem_init,cobrems,cobrems_varbin,integral_xsctn,verbose_output,
     &	bakmaev,heitler,berlin,para,radiation,single_radiation,hist_elasticity,hist_mm,
     &	zero_ninety,GlueX,hypgeom_radiation
c
c	Standard CPP configuration ***** These assignments are now made through a conditional statement later on
c	data ztgt,E0,E_coherent,pol,theta_min,theta_max /82,11.,5.7,1.0,0.75, 5.3 /   !Standard CPP configuration, min angle to TOF and max angle to TOF
c
c	Standard GlueX configuration
C 	data ztgt,E0,E_coherent,pol,theta_min,theta_max /1,11.0,8.8,1.0,1.50,13.12/	!standard GlueX, use these params for B.H. analysis, production angle min, max angle (deg) to TOF
c
c	*****Set tagging interval***** These assignments are now made through a conditional statement later on
c	data E_lo,E_hi /5.2,5.8/	!CPP
C	data E_lo,E_hi /8.2,8.8/	!Gluex
c
	character(len=200) :: genDir, vectors, hists, outputDirTop, oFile
	character(len=200) :: Whist, Xhist, Thist,Tvarhist, Elashist 
	character(len=200) :: MMhist, Thetahist, JThist, Egammahist
	character(len=20) :: arg_seed, arg_jobnum, arg_nevents
	integer :: num_args, arg_seed_int, arg_jobnum_int, arg_nevents_int

	data itest(1),itest(2),itest(3),itest(4),nevent 
     &	/100000,1000000,10000000,100000000,RBHGNEVENTS/  !the variable nevent is where you set the number of output events.
	data m_e,m_muon,m_pi,hbarc /0.00051099,0.105658,0.139570,0.197326/
c	
c  Histogram parameters
	data delta_w,delta_t,delta_x,delta_phi,frac_delta_t,delta_Egamma,delta_theta,delta_q2_ratio,delta_elasticity,delta_mm 
     &	/0.01,.001,0.02,2.0,0.2,0.05,0.02,.01,.01,.02/	
c	units of GeV, GeV^2, dimensionless, degrees, dimensionless, GeV, dimensionless
c
c  Target information
	data Atgt(1), tgtlen(1), radlen(1) /1., .0338,63.04/	!tgtlen = # Rad lengths, radlen= rad length of material in g/cm^2
c								this target is based on 30 cm LH2
	data Atgt(6), tgtlen(6), radlen(6), c_den(6), a_den(6)  /12.,  .05, 42.70, 2.45,  .524/	!RL is in g/cm^2
	data Atgt(14),tgtlen(14),radlen(14),c_den(14),a_den(14) /28., .05, 21.82, 3.14,  .537/ !units of c_den and a_den are fm
	data Atgt(20),tgtlen(20),radlen(20),c_den(20),a_den(20) /40., .05, 16.12, 3.51,  .563/ !target length is in units of RL
	data Atgt(26),tgtlen(26),radlen(26),c_den(26),a_den(26) /56., .05, 13.84, 3.971, .5935/ !RL is g/cm^2
	data Atgt(50),tgtlen(50),radlen(50),c_den(50),a_den(50) /116.,.05, 8.82,  5.416, .552/
	data Atgt(82),tgtlen(82),radlen(82) /208.,.05,6.37/  
c
	common /charge_density/ rho0,c_den,a_den		!rho0 is the overall normalization factor for nuclear charge densities
c
	double precision ZBQLUAB,zlo,zhi
	data zlo,zhi /0.,1./
c
c      BLOCK DATA ZBQLBD01
*
*       Initializes seed array etc. for random number generator.
*       The values below have themselves been generated using the
*       NAG generator.
*
      DOUBLE PRECISION ZBQLIX(43),B,C
      COMMON /ZBQL0001/ ZBQLIX,B,C
c      INTEGER I
      DATA (ZBQLIX(I),I=1,43) /8.001441D7,5.5321801D8,
     +1.69570999D8,2.88589940D8,2.91581871D8,1.03842493D8,
     +7.9952507D7,3.81202335D8,3.11575334D8,4.02878631D8,
     +2.49757109D8,1.15192595D8,2.10629619D8,3.99952890D8,
     +4.12280521D8,1.33873288D8,7.1345525D7,2.23467704D8,
     +2.82934796D8,9.9756750D7,1.68564303D8,2.86817366D8,
     +1.14310713D8,3.47045253D8,9.3762426D7 ,1.09670477D8,
     +3.20029657D8,3.26369301D8,9.441177D6,3.53244738D8,
     +2.44771580D8,1.59804337D8,2.07319904D8,3.37342907D8,
     +3.75423178D8,7.0893571D7 ,4.26059785D8,3.95854390D8,
     +2.0081010D7,5.9250059D7,1.62176640D8,3.20429173D8,
     +2.63576576D8/
      DATA B / 4.294967291D9 /
      DATA C / 0.0D0 /
c
	! Note: iseed will be set from command-line args or compiled default
	! RNG initialization moved to after argument parsing (see below)
c
c  Initializations
c
	temp=-1.
	pi = dacos(temp)
	sqrt2=sqrt(2.)
c

c	Initialize histogram arrays, parameters
	do i=1,200
		data_array_w(i)		=0.
		data_array_x(i)		=0.
		data_array_t(i)		=0.
		data_array_elasticity(i)	=0.
		data_array_mm(i)		=0.
		data_array_phi_JT(i)	=0.
		data_array_t_variable(i)=0.
		data_array_Egamma(i)	=0.
		data_array_theta(i)	=0.
		error_array(i)		=0.
	end do

c
c	Command-line argument support for job parameters
c	=====================================================
c	This executable accepts 3 optional command-line arguments:
c	  Argument 1: Random seed (integer)
c	  Argument 2: Job number (integer, used in output filenames)
c	  Argument 3: Number of events to generate (integer)
c
c	The compiled-in values RBHGSEEDVALUE, RBHGJOBNUMBER, RBHGNEVENTS
c	represent the FIRST job's parameters (job 0). When running multiple
c	jobs with the same executable, pass different arguments to each job.
c	
c	Example usage:
c	  Job 0: ./exe 12345 0 5000     (uses args, matches compiled values)
c	  Job 1: ./exe 67890 1 5000     (uses args, overrides compiled values)
c	  Job 2: ./exe 24680 2 5000     (uses args, overrides compiled values)
c	
c	The complete job parameter history is saved in rbhg_config.json
c	in the job_parameters array for full audit trail.
c	=====================================================
c
c	Parse command-line arguments (seed, job_number, nevents)
c	Usage: ./lepton_event_v113.exe <seed> <job_num> <nevents>
	num_args = command_argument_count()
	if (num_args >= 3) then
		call get_command_argument(1, arg_seed)
		call get_command_argument(2, arg_jobnum)
		call get_command_argument(3, arg_nevents)
		
		! Convert strings to integers
		read(arg_seed, *) arg_seed_int
		read(arg_jobnum, *) arg_jobnum_int
		read(arg_nevents, *) arg_nevents_int
		
		! Override compiled-in values with command-line args
		iseed = arg_seed_int
		nevent = arg_nevents_int
		
		print *, 'Command-line args: seed=', iseed, 
     &           ' job=', arg_jobnum_int, ' nevents=', nevent
	else
		! No args: use system clock for random seed (iseed=0 triggers clock-based seeding)
		iseed = 0
		arg_jobnum_int = RBHGJOBNUMBER
		print *, 'Using compiled defaults: seed=CLOCK', 
     &           ' job=', arg_jobnum_int, ' nevents=', nevent
	endif
	
	! Initialize random number generator with the seed (from args or clock if iseed=0)
	call ZBQLINI(ISEED)


c
c  Start logical assignments
c
c******* Histogram output control
	hist_w        =RBHGHIST_W
	hist_x        =RBHGHIST_X
	hist_t        =RBHGHIST_T
	hist_elasticity =RBHGHIST_ELASTICITY
	hist_mm         =RBHGHIST_MM
	hist_theta    =RBHGHIST_THETA
	hist_phi_JT    =RBHGHIST_PHI_JT
	hist_t_variable    =RBHGHIST_T_VARIABLE
	hist_Egamma    =RBHGHIST_EGAMMA
	output_event    =RBHGOUTPUT_EVENT
c*******Output control
	integral_xsctn    =RBHGINTEGRAL_XSCTN
	verbose_output    =RBHGVERBOSE_OUTPUT
c****** The B.H. cross section: only one of the following 3 logicals should be true
	bakmaev        =RBHGBAKMAEV
	heitler        =RBHGHEITLER
	berlin        =RBHGBERLIN
c****** Setting run conditions: polarization direction, Brem type, and GlueX or CPP
	para        =RBHGPARA
	zero_ninety    =RBHGZERO_NINETY
	cobrems        =RBHGCOBREMS
	cobrems_varbin =RBHGCOBREMS_VARBIN
	if (.not.cobrems) pol=0. !fixes a possible problem, incoherent Brem is necessarily unpolarized
	GlueX        =RBHGGLUEX
c*******e+e- or mu+mu-
	electron        =RBHGELECTRON
	muon=.not.electron	!true if electron=.false.
c*******	Radiation control
	radiation    =RBHGRADIATION
	if (muon) radiation =.false. !if muon then turn off internal radiation
	single_radiation	=RBHGSINGLE_RAD	!if true then only one lepton radiates, if false then both leptons radiate
	hypgeom_radiation	=RBHGHYPGEOM	!future: hypergeometric radiation for muons
c******* Form factor choice
	nuc_FF        =RBHGNUC_FF	!set true for nuclear form factor, set false for FF=1
C       
c       End logical assignments
c       
	genDir = 'RBHGHOMEDIRECTORY'
	outputDirTop = 
     &  'RBHGOUTPUTDIRTOP'
	vectors = trim(outputDirTop) // 'RBHGVECTORSPATH'
	hists = trim(outputDirTop) // 'RBHGHISTSPATH'
	
	! Build output filenames using job number (from cmdline or compiled-in)
	write(arg_jobnum, '(I0)') arg_jobnum_int
	oFile = trim(vectors) // '/vectors_' // trim(arg_jobnum) // '.txt'
	Whist = trim(hists) // '/array_w_' // trim(arg_jobnum) // '.txt'
	Xhist = trim(hists) // '/array_x_' // trim(arg_jobnum) // '.txt'
	Thist = trim(hists) // '/array_t_' // trim(arg_jobnum) // '.txt'
	Tvarhist = trim(hists) // '/array_tvar_' // trim(arg_jobnum) 
     &          // '.txt'
	Elashist = trim(hists) // '/array_elasticity_' 
     &           // trim(arg_jobnum) // '.txt'
	MMhist = trim(hists) // '/array_mm_' // trim(arg_jobnum) // '.txt'
	Thetahist = trim(hists) // '/array_theta_' // trim(arg_jobnum) 
     &            // '.txt'
	JThist = trim(hists) // '/array_phiJT_' // trim(arg_jobnum) 
     &         // '.txt'
	Egammahist = trim(hists) // '/array_Egamma_' // trim(arg_jobnum) 
     &           // '.txt'
c

c*******Set run conditions to either GlueX or CPP
	if (GlueX) then
	   E_lo=RBHG_ELO_GLUEX	!tagging interval
	   E_hi=RBHG_EHI_GLUEX
	   ztgt=RBHG_ZTGT_GLUEX
	   E0=RBHG_E0_GLUEX
	   E_coherent=RBHG_ECOHERENT_GLUEX
	   pol=RBHG_POL_GLUEX
	   theta_min=RBHG_THETAMIN_GLUEX
	   theta_max=RBHG_THETAMAX_GLUEX
c       
	else			!then CPP
	   E_lo=RBHG_ELO_CPP
	   E_hi=RBHG_EHI_CPP
	   ztgt=RBHG_ZTGT_CPP
	   E0=RBHG_E0_CPP
	   E_coherent=RBHG_ECOHERENT_CPP
	   pol=RBHG_POL_CPP
	   theta_min=RBHG_THETAMIN_CPP
	   theta_max=RBHG_THETAMAX_CPP
	end if
C       
c       .  Set target mass in GeV
	mtgt=Atgt(ztgt)*.931494
	if (ztgt.eq.1) mtgt=0.93827208
c       
c       Set lepton mass: 
	if (muon) 	m_part=m_muon
	if (electron) 	m_part=m_e
C       
c       .  Set radiation constants, see details in the Walker PhD
	eta_rad=log(1194./ztgt**(2./3.))/log(184.15/ztgt**(1./3.))
	b_rad=4./3.*( 1.+1./9.*(ztgt+1.)/(ztgt+eta_rad)/log(184.15/ztgt**(1./3.)) )
c       
	if (output_event) open(unit=4,status='new',file=oFile)
c       
c       .  Initialize Brem. distribution: select 1/Egamma or coherent Brems. file
	if (cobrems) then	! read coherent Brem. file
	   brem_init=.true.
	   temp=Brem(brem_init,cobrems,cobrems_varbin,E0,Egamma,GlueX)	!read coherent brems file, then set brem_init = .false.
	endif
c       
c       Phase space assignment:
	phase_space= 8		!theta=xs**phase_space, with integer phase_space >1 . Phase_space=8 is the fastest for e+e- and also mu+mu-
c       !Set phase_space=0 for standard dcos theta/dx =1
	Rexp=real(phase_space)
c       
c       
c       Initialize the nuclear form factor for the 2-parameter Fermi model if needed:
	if (nuc_FF.and.ztgt.ne.1.and.ztgt.ne.82) then !use 2-parameter Fermi model	
	   c_den(ztgt)=c_den(ztgt)/hbarc
	   a_den(ztgt)=a_den(ztgt)/hbarc
	   call density_init(ztgt) !for initializing the nuclear form factor 
	end if
c       
c***********************************
c	Evaluate the cross section at one kinematic point at coherent peak
	x1=0.4
	theta1=  theta_min*(pi/180.)	
	theta2=  theta_min*(pi/180.)
	phi1=   90.*(pi/180.)
	phi2=   270.*(pi/180.)	
	cross_section=xsctn(E_coherent,ztgt,x1,x2,theta1,phi1,theta2,phi2,pol,m_part,mtgt,t,
     &	E1,p1,k1,E2,p2,k2,Etgt,ktgt,nuc_FF,bakmaev,heitler,berlin,para)
	if(verbose_output) print *, ' cross section nb/sr^2 = ', cross_section	
c*************************************
	theta_min=theta_min*pi/180.	!switch to radians
	theta_max=theta_max*pi/180.
	cos_max=cos(theta_min)
	cos_min=cos(theta_max)
c Limits on xs
	if (phase_space.eq.0) then 	 !dcos theta/dx =1
		xs_max=cos_max	
		xs_min=cos_min
	else
		xs_max=theta_max**(1./Rexp)
		xs_min=theta_min**(1./Rexp)
	end if
c
c***************************************************************************
	print *, 'Finding maximum cross section (this may take a while)...'
c	cross_max=0.
	do j=1,4	!loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the maximum 
			!cross section*Brem converges
		print *, 'Starting sampling iteration', j, 'with', itest(j), 'test events'
		cross_max_old=cross_max
	do i=1,itest(j)	!find maximum cross section in allowed phase space
		if (i == 1) print *, 'Entered inner loop, calling xsctn...'
25		continue
		Egamma=ZBQLUAB(E_lo,E_hi)	!get tagged photon energy
		x1=0.5		!x_min+(x_max-x_min)*ZBQLUAB(zlo,zhi)	!make a guess for the energy fraction
		if (para) then 
			phi1=90.*pi/180.	!2.*pi*ZBQLUAB(zlo,zhi)		!make a guess for phi1
			phi2=270.*pi/180.	!2.*pi*ZBQLUAB(zlo,zhi)		!make a guess for phi2
		else
			phi1=0.*pi/180.	!2.*pi*ZBQLUAB(zlo,zhi)		!make a guess for phi1
			phi2=180.*pi/180.	!2.*pi*ZBQLUAB(zlo,zhi)		!make a guess for phi2
		end if		
		xs1=ZBQLUAB(xs_min,xs_max)	 
		xs2=ZBQLUAB(xs_min,xs_max)	 
		if (phase_space.eq.0) then	! dcos theta/dx = 1
			theta1=acos(xs1)
			theta2=acos(xs2)
			jacobian=1.
		else
			theta1=xs1**phase_space
			theta2=xs2**phase_space
			jacobian=(Rexp*xs1**(phase_space-1)*sin(xs1**phase_space))*(Rexp*xs2**(phase_space-1)*sin(xs2**phase_space))
		end if
		cross_section=xsctn(Egamma,ztgt,x1,x2,theta1,phi1,theta2,phi2,pol,m_part,mtgt,t,
     &		E1,p1,k1,E2,p2,k2,Etgt,ktgt,nuc_FF,bakmaev,heitler,berlin,para)*
     &		jacobian*Brem(brem_init,cobrems,cobrems_varbin,E0,Egamma,GlueX)
c
c		print *, 'phi1', phi1, 'phi2', phi2
c
		if (cross_section.gt.cross_max)  then 
			cross_max=cross_section
			Egamma_max=Egamma
			xmax=x1
			theta1_max=theta1*180./pi
			theta2_max=theta2*180./pi
			phi1_max=phi1*180./pi
			phi2_max=phi2*180./pi
		end if
	end do
		delta_cross_max=(cross_max-cross_max_old)/cross_max*100.
		if(verbose_output) print *, 'test events', itest(j), 'maximum xsctn*Brem' , cross_max,' percent change ' , 
     &		delta_cross_max 
		if(verbose_output) print *, 'Egamma', Egamma_max, 'x', xmax, 'theta1', theta1_max, 'theta2', theta2_max,
     &		'phi1', phi1_max, 'phi2', phi2_max 
	end do
c
c**********************************************************************************************************
c Loop over 4 samplings of the phase space at coherent peak, each a factor of x10 larger, to see if the integrated cross section converges
c Set limits on energy fraction
	if (integral_xsctn) then
c
	x_max=(E_coherent-m_part)/E_coherent
	x_min=m_part/E_coherent
c
	total_xscn=0.
	do j=1,4	!loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the integrated
			!cross section at the coherent peak converges
		total_xscn_old=total_xscn
		cross_sum=0.
	do i=1,itest(j)	
50		continue
		x1=ZBQLUAB(x_min,x_max)	!energy fraction
		phi1=2.*pi*ZBQLUAB(zlo,zhi)
		phi2=2.*pi*ZBQLUAB(zlo,zhi)
		xs1=ZBQLUAB(xs_min,xs_max)	
		xs2=ZBQLUAB(xs_min,xs_max)	 
		if (phase_space.eq.0) then	! dcos theta/dx = 1
			theta1=acos(xs1)
			theta2=acos(xs2)
			jacobian=1.
		else
			theta1=xs1**phase_space
			theta2=xs2**phase_space
			jacobian=(Rexp*xs1**(phase_space-1)*sin(xs1**phase_space))*(Rexp*xs2**(phase_space-1)*sin(xs2**phase_space))
		end if
		cross_section=xsctn(E_coherent,ztgt,x1,x2,theta1,phi1,theta2,phi2,pol,m_part,mtgt,t,
     &		E1,p1,k1,E2,p2,k2,Etgt,ktgt,nuc_FF,bakmaev,heitler,berlin,para)*
     &		jacobian
c
		cross_sum=cross_sum+cross_section
	end do
		total_xscn=cross_sum/float(itest(j))*(xs_max-xs_min)**2*(2.*pi)**2*(x_max-x_min)
		delta_total_xscn=(total_xscn-total_xscn_old)/total_xscn*100.
		print *, 'test events', itest(j), 'Egamma', E_coherent, 'total cross section nb', total_xscn, 
     &		' percent change ', delta_total_xscn
	end do
c
	end if
c*************************************************************************
c
c	Start event generation
c
c	Use the widest possible range in x by using the maximum accepted tagged photon energy, then test it
	x_max=(E_hi-m_part)/E_hi	!largest possible x
	x_min=m_part/E_hi		!smallest possible x
c
	nfail=0
	bad_max=0
c
	print *, 'Starting event generation...'
	do i=1,nevent
		if (mod(i,1000) == 0) print *, 'Generated', i, 'of', nevent, 'events'
100		continue			!try again
		Egamma=ZBQLUAB(E_lo,E_hi)	!get tagged photon energy
		x1=ZBQLUAB(x_min,x_max)		!energy fraction	
c	Test x1 to make sure it's within the allowed range for the photon energy Egamma
		if (x1.ge.((Egamma-m_part)/Egamma).or.x1.le.(m_part/Egamma)) go to 100 ! x1 is out of range, try again
c
		phi1=2.*pi*ZBQLUAB(zlo,zhi)
		phi2=2.*pi*ZBQLUAB(zlo,zhi)
		xs1=ZBQLUAB(xs_min,xs_max)
		xs2=ZBQLUAB(xs_min,xs_max)
c
		if (phase_space.eq.0) then	! dcos theta/dx = 1
			theta1=acos(xs1)
			theta2=acos(xs2)
			jacobian=1.
		else 
			theta1=xs1**phase_space
			theta2=xs2**phase_space
			jacobian=(Rexp*xs1**(phase_space-1)*sin(xs1**phase_space))*(Rexp*xs2**(phase_space-1)*sin(xs2**phase_space))
		end if
c
		cross_section=xsctn(Egamma,ztgt,x1,x2,theta1,phi1,theta2,phi2,pol,m_part,mtgt,t,
     &		E1,p1,k1,E2,p2,k2,Etgt,ktgt,nuc_FF,bakmaev,heitler,berlin,para)*
     &		jacobian*Brem(brem_init,cobrems,cobrems_varbin,E0,Egamma,GlueX)
c
		if (cross_section.eq.0.) go to 100
c
		if (cross_section.gt.cross_max) then	!bad cross max
			bad_max=bad_max+1	!cross section is larger than cross_max, not supposed to happen
			if(verbose_output) print *, ' bad max cross section= ', cross_section, ' max cross section update '
			cross_max=cross_section
		end if
		cross_test=cross_max*ZBQLUAB(zlo,zhi)
		if (cross_test.gt.cross_section) then	!selection fails
			nfail=nfail+1
			go to 100
		end if
c
c	******************************Event selection succeeds!!! Do analysis and internal radiation if required ********************************* 
c
c
c	Do radiation if called for: 
c
		Eloss1=0.
		Eloss2=0.
		do j=1,3
			Ebrem1(j)=0.
			Ebrem2(j)=0.
		end do
c
		if (radiation) then	!allow internal radiation to change E1 and E2. 
C					 Use unradiated t to calculate internal rad length. 
			if (single_radiation) then
				test=ZBQLUAB(zlo,zhi)				
				if (test.lt.0.5) then ! radiate lepton 1		
					E1=E_rad(b_rad,m_part,t,E1,Eloss1)	!radiated E1, recalculate 3-vector for lepton 1
					x1=E1/Egamma
					temp=E1**2-m_part**2
					if (temp.lt.0.) temp=0.			
					p1=sqrt(temp)
					k1(1)=p1*sin(theta1)*cos(phi1)
					k1(2)=p1*sin(theta1)*sin(phi1)
					k1(3)=p1*cos(theta1)
					Ebrem1(1)=Eloss1*sin(theta1)*cos(phi1)
					Ebrem1(2)=Eloss1*sin(theta1)*sin(phi1)
					Ebrem1(3)=Eloss1*cos(theta1)		
				else			! radiate lepton 2	
					E2=E_rad(b_rad,m_part,t,E2,Eloss2)	!radiated E2, recalculate 3-vector for lepton 2
					x2=E2/Egamma
					temp=E2**2-m_part**2
					if (temp.lt.0.) temp=0.			
					p2=sqrt(temp)
					k2(1)=p2*sin(theta2)*cos(phi2)
					k2(2)=p2*sin(theta2)*sin(phi2)
					k2(3)=p2*cos(theta2)
					Ebrem2(1)=Eloss2*sin(theta2)*cos(phi2)
					Ebrem2(2)=Eloss2*sin(theta2)*sin(phi2)
					Ebrem2(3)=Eloss2*cos(theta2)		
				end if							
			else	!then allow both leptons to radiate
c			! radiate lepton 1
				E1=E_rad(b_rad,m_part,t,E1,Eloss1)	!radiated E1, recalculate 3-vector for lepton 1
				x1=E1/Egamma
				temp=E1**2-m_part**2
				if (temp.lt.0.) temp=0.			
				p1=sqrt(temp)
				k1(1)=p1*sin(theta1)*cos(phi1)
				k1(2)=p1*sin(theta1)*sin(phi1)
				k1(3)=p1*cos(theta1)
				Ebrem1(1)=Eloss1*sin(theta1)*cos(phi1)
				Ebrem1(2)=Eloss1*sin(theta1)*sin(phi1)
				Ebrem1(3)=Eloss1*cos(theta1)		
c			! radiate lepton 2	
				E2=E_rad(b_rad,m_part,t,E2,Eloss2)	!radiated E2, recalculate 3-vector for lepton 2
				x2=E2/Egamma
				temp=E2**2-m_part**2
				if (temp.lt.0.) temp=0.			
				p2=sqrt(temp)
				k2(1)=p2*sin(theta2)*cos(phi2)
				k2(2)=p2*sin(theta2)*sin(phi2)
				k2(3)=p2*cos(theta2)
				Ebrem2(1)=Eloss2*sin(theta2)*cos(phi2)
				Ebrem2(2)=Eloss2*sin(theta2)*sin(phi2)
				Ebrem2(3)=Eloss2*cos(theta2)		
			end if	! end of 'single_radiation' test
		end if	! end of 'radiation' test
C
c
C	Create ghost particle that gives energy/momentum conservation when summed with target and lepton energy/momentum 
C	
		kghost(1)=-k1(1)-k2(1)-ktgt(1)	!momenta of ghost particle required by momentum conservation
		kghost(2)=-k1(2)-k2(2)-ktgt(2)	
		kghost(3)=-k1(3)-k2(3)-ktgt(3)+Egamma
		Eghost	= Egamma+mtgt-E1-E2-Etgt	!energy of ghost particle required by momentum conservation
c
c
		do j=1,3
			ks(j)=k1(j)+k2(j)	!lepton summed momentum
		end do
		t=ks(1)**2+ks(2)**2+(Egamma-ks(3))**2-(Egamma-E1-E2)**2  !(radiated) 4-momentum transfer squared to nucleus, this is positive
		q=sqrt(ks(1)**2+ks(2)**2+(Egamma-ks(3))**2)	!(radiated) 3-momentum transfer to nucleus		
		w_mumu=sqrt((E1+E2)**2-ks(1)**2-ks(2)**2-ks(3)**2)	!(radiated) invariant mass
		elasticity=(E1+E2)/Egamma				!(radiated) elasticity
		missing_mass=sqrt( (Egamma+mtgt-E1-E2)**2-ks(1)**2-ks(2)**2-(Egamma-ks(3))**2 ) !(radiated) missing mass
c
		if (.not.zero_ninety) then !rotate coordinate system by -45 deg, so the former "para" direction is now at 45 deg, and
c					  the former "perp" direction is at 135 deg. Find 3-vectors in this coordinate system
			vec(1)=k1(1)
			vec(2)=k1(2)
			k1(1)=vec(1)/sqrt2-vec(2)/sqrt2
			k1(2)=vec(1)/sqrt2+vec(2)/sqrt2
c			
			vec(1)=k2(1)
			vec(2)=k2(2)
			k2(1)=vec(1)/sqrt2-vec(2)/sqrt2
			k2(2)=vec(1)/sqrt2+vec(2)/sqrt2
C
			vec(1)=Ebrem1(1)
			vec(2)=Ebrem1(2)
			Ebrem1(1)=vec(1)/sqrt2-vec(2)/sqrt2
			Ebrem1(2)=vec(1)/sqrt2+vec(2)/sqrt2	
c
			vec(1)=Ebrem2(1)
			vec(2)=Ebrem2(2)
			Ebrem2(1)=vec(1)/sqrt2-vec(2)/sqrt2
			Ebrem2(2)=vec(1)/sqrt2+vec(2)/sqrt2	
c
			vec(1)=ktgt(1)
			vec(2)=ktgt(2)
			ktgt(1)=vec(1)/sqrt2-vec(2)/sqrt2
			ktgt(2)=vec(1)/sqrt2+vec(2)/sqrt2	
		end if		
c
		if (berlin) then	!calculate polarized Berlin JT angle
			do j=1,2
				JT(j)=   2.*E2/(E1-k1(3))*k1(j)   +2.*E1/(E2-k2(3))*k2(j)	!vector current, units of GeV^-1
			end do
			phi_JT=acos(JT(1)/sqrt(JT(1)**2+JT(2)**2))	!phi angle of JT wrt to x axis, radians
			if (JT(2).lt.0.) phi_JT=2.*pi-phi_JT
c
		else if (bakmaev) then
			p1_tran=sqrt(k1(1)**2+k1(2)**2)
			p2_tran=sqrt(k2(1)**2+k2(2)**2)
			c1=p1_tran**2+m_part**2
			c2=p2_tran**2+m_part**2
			do j=1,2
				JT(j)=k1(j)/c1+k2(j)/c2	!vector current, units of GeV^-1
			end do
			phi_JT=acos(JT(1)/sqrt(JT(1)**2+JT(2)**2))	!phi angle of JT wrt to x axis, radians
			if (JT(2).lt.0.) phi_JT=2.*pi-phi_JT
		end if
c							
c	********Do the histograming you need to do:
c
		if (hist_w) then
			i_array=int(w_mumu/delta_w)+1			!w distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_w(i_array)=data_array_w(i_array)+1.
		end if
c
		if (hist_x) then
			i_array=int(x1/delta_x)+1		!	x distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_x(i_array)=data_array_x(i_array)+1.
		end if
c
		if (hist_t) then
			i_array=int(t/delta_t)+1		!	t distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_t(i_array)=data_array_t(i_array)+1.
		end if
c						
		if (hist_elasticity) then
			i_array=int(elasticity/delta_elasticity)+1	!  elasticity distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_elasticity(i_array)=data_array_elasticity(i_array)+1.
		end if
c
		if (hist_mm) then
			i_array=int(missing_mass/delta_mm)+1		!missing mass distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_mm(i_array)=data_array_mm(i_array)+1.
		end if
c
		if (hist_phi_JT) then
			i_array=int(phi_JT*180./pi/delta_phi)+1 !	J_T phi distribution in degrees
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_phi_JT(i_array)=data_array_phi_JT(i_array)+1.
		end if
c
		if (hist_t_variable) then
			if ( 0.0   .le.t.and. t.lt.0.0002)  	i_array=1
			if ( 0.0002.le.t.and. t.lt.0.0006) 	i_array=2
			if ( 0.0006.le.t.and. t.lt.0.0013) 	i_array=3
			if ( 0.0013.le.t.and. t.lt.0.002) 	i_array=4
			if ( 0.002. le.t.and. t.lt.0.003)		i_array=5
			if ( 0.003. le.t.and. t.lt.0.0042)	i_array=6
			if ( 0.0042.le.t.and. t.lt.0.0055)	i_array=7
			if ( 0.0055.le.t.and. t.lt.0.007)		i_array=8
			if ( 0.007. le.t.and. t.lt.0.009)		i_array=9
			if (t.ge.0.009) i_array=int((t-0.009)/.002)+10
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_t_variable(i_array)=data_array_t_variable(i_array)+1.
		end if
c
		if (hist_Egamma) then
			i_array=int((Egamma-E_lo)/delta_Egamma)+1	!photon energy distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_Egamma(i_array)=data_array_Egamma(i_array)+1.
		end if
c
		if (hist_theta) then
			i_array=int(theta1*180./pi/delta_theta)+1	!theta1 distribution
			if (i_array.lt.1) i_array=1
			if (i_array.gt.200) i_array=200		
			data_array_theta(i_array)=data_array_theta(i_array)+1.
		end if
c
c********************* Output 3-vectors ***************************
c

c	Output 4-vectors vectors
		if (output_event) then
c
			write(4,200) 
     &			Egamma,E1,k1(1),k1(2),k1(3),E2,k2(1),k2(2),k2(3),Etgt,ktgt(1),ktgt(2),ktgt(3)

		end if
200		format(2x,f18.15,12(','e23.16))
c
c************************************************************
		if ((verbose_output).and.(mod(i,100).eq.0)) print *,' event # ', i
	end do
c
c   Event generation ends
c	
	failure=float(nfail)/float(nevent)
	if(verbose_output) print *, 'phase space parameter = ', phase_space, 'failures/event = ', 
     &		failure , ' Events exceeding max xsctn = ', bad_max
c
c   Write out the histograms 
c
	if (hist_w) then 
		open(unit=2,status='new',file= Whist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_w
			error_array(i)=sqrt(data_array_w(i))
			write(2,20) x_value, data_array_w(i), error_array(i)		
		end do
		close(2)
	end if	
c
	if (hist_x) then
		open(unit=2,status='new',file= Xhist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_x
			if (x_value.gt.1.0) go to 40
			error_array(i)=sqrt(data_array_x(i))
			write(2,20) x_value, data_array_x(i), error_array(i)	
		end do
40		close(2)
	end if
c
	if (hist_t) then 
		open(unit=2,status='new',file= Thist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_t
			error_array(i)=sqrt(data_array_t(i))
			write(2,20) x_value, data_array_t(i), error_array(i)		
		end do
		close(2)
	end if
c
	if (hist_elasticity) then 
		open(unit=2,status='new',file= Elashist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_elasticity
			error_array(i)=sqrt(data_array_elasticity(i))
			write(2,20) x_value, data_array_elasticity(i), error_array(i)		
		end do
		close(2)
	end if
c
	if (hist_mm) then 
		open(unit=2,status='new',file= MMhist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_mm
			error_array(i)=sqrt(data_array_mm(i))
			write(2,20) x_value, data_array_mm(i), error_array(i)		
		end do
		close(2)
	end if
c
	if (hist_phi_JT) then 
		open(unit=2,status='new',file= JThist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_phi
			if (x_value.gt.360.) go to 41
			error_array(i)=sqrt(data_array_phi_JT(i))
			write(2,20) x_value, data_array_phi_JT(i), error_array(i)		
		end do
41		close(2)
	end if
c
	if (hist_t_variable) then
		open(unit=2,status='new',file= Tvarhist)
		do i=1,200		!write out the arrays			
			if (i.eq.1) then
				x_value=0.0001
				bin_width=0.0002
			else if (i.eq.2) then
				x_value=0.0004
				bin_width=0.0004
			else if (i.eq.3) then
				x_value=0.00095
				bin_width=0.0007
			else if (i.eq.4) then
				x_value=0.0018
				bin_width=0.0007
			else if (i.eq.5) then
				x_value=0.0025
				bin_width=0.001
			else if (i.eq.6) then
				x_value=0.0036
				bin_width=0.0012
			else if (i.eq.7) then 
				x_value=0.0049
				bin_width=0.0013
			else if (i.eq.8) then
				x_value=0.0063
				bin_width=0.0015
			else if (i.eq.9) then
				x_value=0.008
				bin_width=0.0018
			else if (i.ge.10) then 
				x_value=float(i-9)*0.002+0.008
				bin_width=0.002
			end if
			delta=bin_width/2.
			write(2,30) x_value,delta,data_array_t_variable(i)
		end do
		close(2)
	end if
c
	if (hist_Egamma) then
		open(unit=2,status='new',file= Egammahist)
		do i=1,200		!write out the arrays
			x_value=E_lo+(float(i)-.5)*delta_Egamma
			error_array(i)=sqrt(data_array_Egamma(i))
			write(2,20) x_value,data_array_Egamma(i),error_array(i)
		end do
		close(2)
	end if
c
	if (hist_theta) then
		open(unit=2,status='new',file= Thetahist)
		do i=1,200		!write out the arrays
			x_value=(float(i)-.5)*delta_theta
			error_array(i)=sqrt(data_array_theta(i))
			write(2,20) x_value,data_array_theta(i),error_array(i)
		end do
		close(2)
	end if
c
c
20	format(2x,f10.5,2x,f10.0,2x,f10.1)
30	format(2x,2(e10.4,2x),f12.0)
c
	stop
	end
c
c*******************************************************************
c
	real*8 function Brem(brem_init,cobrems,cobrems_varbin,E0,Egamma,GlueX)
	implicit none
	real*8 E0,Egamma,Eg(500),Br(500)
	real*8 Eg_low(1000),Eg_high(1000),Br_var(1000)
	integer*8 i,imax,imax_var,ipoint
	logical brem_init,cobrems,cobrems_varbin,GlueX
	character(len=200) :: genDir2, templateDir, cobremsGlueX, cobremsCPP, cobremsGlueX_var, cobremsCPP_var
	character(len=200) :: line_buffer
	common /Brem_spect/ Eg,Br,Eg_low,Eg_high,Br_var,imax,imax_var		!this common block is required to hold over stored values for Eg and Br to 
C					!next function call to Brem! Added varbin arrays too.
c	save Eg_low,Eg_high,Br_var,imax_var	!save variable bin arrays between calls  
c	
	Brem=0.
c
	! Note: selection moved into initialization to avoid repeated calls/prints
	! Determine cobrems filenames based on configuration (will be selected at init)
c	! call select_cobrems_files(GlueX, cobrems_varbin, 
c     !      cobremsGlueX, cobremsCPP, cobremsGlueX_var, 
c     !      cobremsCPP_var, templateDir)

	if (brem_init) then  !open and read coherent Brem file
		! Set up paths for cobrems files
		genDir2 = 'RBHGHOMEDIRECTORY'
		templateDir = trim(genDir2) // '/template_RBHG/'
		
		! Select file names once at initialization to avoid repeating selection and prints
		call select_cobrems_files(GlueX, cobrems_varbin, 
     &       cobremsGlueX, cobremsCPP, cobremsGlueX_var, 
     &       cobremsCPP_var, templateDir)
		if (cobrems_varbin) then
			! Read variable bin file with validation
			print *, 'Opening variable-width cobrems file...'
			if (GlueX) open(unit=2,file=cobremsGlueX_var)
			if (.not.GlueX) open(unit=2,file=cobremsCPP_var)
			
			! First, read metadata and validate energy limits
			! call validate_cobrems_limits(2, E_lo, E_hi) ! TODO: Fix variable references
			
			! Skip header lines (lines starting with #)
5			continue
			read(2,'(A)',end=20) line_buffer
			if (line_buffer(1:1) == '#') go to 5
			! First non-comment line is data - parse it
			read(line_buffer,*) Eg_low(1), Eg_high(1), Br_var(1)
			
			! Now read the rest of the data
			i=2
10			continue
			read(2,*,end=20) Eg_low(i), Eg_high(i), Br_var(i)
			i=i+1
			go to 10
20			continue
			imax_var=i-1
			close(2)
			print *, 'Variable-width cobrems: loaded', imax_var, 'bins'
			print *, 'Energy range:', Eg_low(1), 'to', Eg_high(imax_var)
		else
			! Read fixed bin file (2 columns: energy counts) - original method
			if (GlueX) open(unit=2,file=cobremsGlueX)
			if (.not.GlueX) open(unit=2,file=cobremsCPP)
			i=1
30			continue
			read(2,*,end=40) Eg(i),Br(i) 
			i=i+1
			go to 30
40			continue
			imax=i-1
			close(2)
		endif
		brem_init=.false.	!done with initialization
		return
	end if
c
	if(.not.cobrems) then
		Brem=E0/Egamma
		return
	end if
c
	if (cobrems) then !return coherent brems distribution
		if (cobrems_varbin) then
			! Variable bin lookup using binary search
			call find_bin(Egamma, Eg_low, Eg_high, imax_var, ipoint)
			if (ipoint > 0 .and. ipoint <= imax_var) then
				Brem = Br_var(ipoint)
			else
				Brem = 0.d0  ! Outside energy range
			endif
		else
			! Fixed bin lookup (original method)
			if(GlueX) ipoint=nint((Egamma +.02)/.04)	!GlueX Brem. File
			if(.not.GlueX) ipoint=nint((EGamma + .12)/.12)	!CPP brem. file
			if (ipoint >= 1 .and. ipoint <= imax) then
				Brem=Br(ipoint)
			else
				Brem=0.d0  ! Outside range
			endif
		endif
		return
	end if
	end
c	
c******************************************************************
c The reference for this is Bakmaev et al., Physics Letters B 660 (2008) 494-500, eqn. 23
c Had to correct a mistake in Eqn. 23 of their paper.   The cos(phi1 + phi2) term should be 
c multiplied by 2.  You can see this by comparing Wp in Eqn. 22 with the vector current part of Eqn. 23
c
	real*8 function xsctn(E0,ztgt,x1,x2,theta1,phi1,theta2,phi2,pol,m_part,mtgt,t,
     &	E1,p1,k1,E2,p2,k2,Etgt,ktgt,nuc_FF,bakmaev,heitler,berlin,para)	!units of nb/sr^2
	implicit none
	real*8 E0,Z,x1,x2,theta1,phi1,theta2,phi2,k1(3),k2(3),pol,W_unpol,W_pol,q2_T,xsctn_point
	real*8 m_part,mtgt,alpha,hbarc
	real*8 pi,E1,E2,k1_mag,k2_mag,p1,p2,p1_tran,p2_tran,Q2,q,c1,c2,JS,JT(2),KT(2),pT(2),		
     &	FF2,FF_nuc,FF_TFM,phi_JT,phi_KT,t,ks(3),A,B,C,Ap,Bp,Cp,E2_tst,D2,w_mumu,missing_mass,
     &	Etgt,ktgt(3),p_mult,W,
     &	x_tmp,E_tmp,p_tmp,k_tmp(3),theta_tmp,phi_tmp
	integer*8 i,ztgt
	logical nuc_FF,bakmaev,heitler,berlin,para
	double precision ZBQLUAB,zlo,zhi
	data zlo,zhi /0.,1./
c
	data alpha,hbarc /7.297352e-3,0.197326/
c
	pi=acos(-1.)
	Z=float(ztgt)
c	
	E1=E0*x1	!use standard expression for x1
c
c   Set cos(2*phi) multiplier
	if (para) then 
		p_mult=1.
	else
		p_mult=-1.
	endif
c
	A=E1**2-m_part**2
	if (A.lt.0.) A=0.
	p1=sqrt(A)
	k1(1)=p1*sin(theta1)*cos(phi1)
	k1(2)=p1*sin(theta1)*sin(phi1)
	k1(3)=p1*cos(theta1)
c
c	Work through the exact equations of energy conservation for the second lepton. Need to be careful with floating point precision
	Ap=mtgt*(E0-E1)-E0*(E1-k1(3))+m_part**2
	Bp=E0+mtgt-E1
	Cp=k1(1)*sin(theta2)*cos(phi2)+k1(2)*sin(theta2)*sin(phi2)+(k1(3)-E0)*cos(theta2)
	A=Cp**2-Bp**2
	B=2.*Ap*Bp
	C=-(Cp**2*m_part**2+Ap**2)
	D2=B**2-4.*A*C
	if (D2.lt.0.) D2=0.
	E2=(-B-sqrt(D2))/(2.*A)	! the negative root seems to give the correct solution
c
c	print *, ' x ', x, 'E2 non cons', E2_tst, ' D2', D2 ,' E2 cons', E2 
	if (E2.lt.m_part) E2=m_part
c
	x2=E2/E0			!standard expression for energy fraction of second lepton
c	
	A=E2**2-m_part**2
	if (A.lt.0.) A=0.	
	p2=sqrt(A)
	k2(1)=p2*sin(theta2)*cos(phi2)
	k2(2)=p2*sin(theta2)*sin(phi2)
	k2(3)=p2*cos(theta2)
c
	do i=1,3
		ks(i)=k1(i)+k2(i)	!lepton summed momentum
	end do
c
	t=ks(1)**2+ks(2)**2+(E0-ks(3))**2-(E0-E1-E2)**2	!4-momentum transfer squared to nucleus, positive quantity
	if (t.lt.0.) t=0.
	Q2=t
	q=sqrt(ks(1)**2+ks(2)**2+(E0-ks(3))**2)	!magnitude 3-momentum transfer to nucleus	! xxx***xxx
C	Target kinematics:
	Etgt=E0+mtgt-E1-E2			!energy of recoil target
	ktgt(1)=-ks(1)
	ktgt(2)=-ks(2)
	ktgt(3)=E0-ks(3)
c
c  Start evaluating cross sections
c
	if (bakmaev) then	!use Bakmaev equations
c
c	Apply the restriction on transverse momenta for Bakmaev:
		p1_tran=sqrt(k1(1)**2+k1(2)**2)
		p2_tran=sqrt(k2(1)**2+k2(2)**2)
		if (p1_tran.lt.10.*m_part.or.p2_tran.lt.10.*m_part) then
			xsctn=0.
			return
		end if		
c
		c1=p1_tran**2+m_part**2
		c2=p2_tran**2+m_part**2
		JS=1./c1-1./c2	!units of 1/GeV^2	!scalar current, units of GeV^-2
		do i=1,2
			JT(i)=k1(i)/c1+k2(i)/c2	!vector current, units of GeV^-1
		end do
c
		phi_JT=acos(JT(1)/sqrt(JT(1)**2+JT(2)**2))	!phi angle of JT wrt to x axis, radians
		if (JT(2).lt.0.) phi_JT=2.*pi-phi_JT
c
     		W_unpol = m_part**2*JS**2+(x1**2+(1.-x1)**2)*(JT(1)**2+JT(2)**2)
		W_pol = -2.*x1*(1.-x1)*(JT(1)**2+JT(2)**2)	!this is my reduction of the Bakmaev equations
c
		xsctn_point=2.*alpha**3*Z**2*E0**4*x1**2*(1.-x1)**2/(pi**2*q2**2)*(W_unpol+p_mult*pol*cos(2.*phi_JT)*W_pol) !contains the cos(2phi_JT) term
     &		*hbarc**2/100.*1.e9 				!units of nb/sr^2   
c
	else if (heitler) then	!  Use Heitler expression for cross section, differential in dx, domega1, domega2, units of nb/sr^2. 
c
		xsctn_point=Z**2*alpha**3/(2.*pi)**2*p1*p2/E0**2/q**4*(		! xxx***xxx
     &		   (p1*sin(theta1)/(E1-p1*cos(theta1)))**2*(q**2-4.*E2**2) 		! xxx***xxx
     &		  +(p2*sin(theta2)/(E2-p2*cos(theta2)))**2*(q**2-4.*E1**2) 		! xxx***xxx
     &		  -2.*p1*p2*sin(theta1)*sin(theta2)*cos(phi2-phi1)/( (E1-p1*cos(theta1))*(E2-p2*cos(theta2)) )*(4.*E1*E2+q**2) 	! xxx***xxx
     &		  +2.*E0**2*( (p1*sin(theta1))**2+(p2*sin(theta2))**2+2.*p1*p2*sin(theta1)*sin(theta2)*cos(phi1-phi2) )
     &			/( (E1-p1*cos(theta1))*(E2-p2*cos(theta2)) )
     &		  )*hbarc**2/100.*1.e9 
c
	else if (berlin)	then! then use vector formulation of the cross section that allows for polarization, derived from Berlin and Madansky
c
			do i=1,2
				JT(i)=   2.*E2/(E1-p1*cos(theta1))*k1(i)   +2.*E1/(E2-p2*cos(theta2))*k2(i)	!vector current, units of GeV^-1
				KT(i)=q/(E1-p1*cos(theta1))*k1(i)-q/(E2-p2*cos(theta2))*k2(i)		! xxx***xxx
				pT(i)=k1(i)+k2(i)
			end do
			phi_JT=acos(JT(1)/sqrt(JT(1)**2+JT(2)**2))	!phi angle of JT wrt to x axis, radians
			if (JT(2).lt.0.) phi_JT=2.*pi-phi_JT
			phi_KT=acos(KT(1)/sqrt(KT(1)**2+KT(2)**2))	!phi angle of KT wrt to x axis, radians
			if (KT(2).lt.0.) phi_KT=2.*pi-phi_KT
c
			xsctn_point=Z**2*alpha**3/(2.*pi)**2*p1*p2/E0**2/q**4*(		! xxx***xxx
     &			-(JT(1)**2+JT(2)**2)+(KT(1)**2+KT(2)**2)+2.*E0**2*(pT(1)**2+pT(2)**2)/
     &			(E1-p1*cos(theta1))/(E2-p2*cos(theta2)) 
     &			+p_mult*pol*( -cos(2.*phi_JT)*(JT(1)**2+JT(2)**2)+cos(2.*phi_KT)*(KT(1)**2+KT(2)**2) )
     &			)*hbarc**2/100.*1.e9 
	end if				
c
	xsctn=xsctn_point*FF2(q,ztgt,nuc_FF) 	!units of nb/sr^2  	! xxx***xxx
c
c    Randomize lepton 1 and lepton 2
c
	W=ZBQLUAB(zlo,zhi)
	if (W.le.0.5) then !flip all the lepton indices
		x_tmp=x2		!assign particle #2 to temp
		E_tmp=E2
		p_tmp=p2
		do i=1,3
			k_tmp(i)=k2(i)
		end do
		theta_tmp=theta2
		phi_tmp=phi2
c				assign particle #1 to particle #2
		x2=x1
		E2=E1
		p2=p1
		do i=1,3
			k2(i)=k1(i)
		end do
		theta2=theta1
		phi2=phi1
C				assign temp to particle #1
		x1=x_tmp
		E1=E_tmp
		p1=p_tmp
		do i=1,3
			k1(i)=k_tmp(i)
		end do
		theta1=theta_tmp
		phi1=phi_tmp
	end if
c
	return
	end
c
c****************************************************************************************
c
	real*8 function FF2(q,ztgt,nuc_FF)	!evaluate nuclear and atomic FFs, and take square of difference	! xxx***xxx
	implicit none
	real*8 q,qA,FF_nuc,FF_atomic,FF,aBohr,a(4),b(4),c,hbarc,pi		! xxx***xxx
	integer*8 ztgt,i
	logical nuc_FF
	data aBohr /268176./	!Hydrogen Bohr radius in GeV^-1, 6 digit precision
c  Sum of Gaussian parameterization of Pb atomic form factor
	data a(1),a(2),a(3),a(4) /31.0617, 13.0637, 18.442, 5.9696/	!dimensionless
	data b(1),b(2),b(3),b(4) /0.6902, 2.3576, 8.618, 47.2579/	!units of Angstrom squared
	data c	/13.4630/	!constant term, dimensionless
c
	data hbarc,pi /.197326,3.14159/
	
c
	if (nuc_FF) then 
		FF_nuc=FF(q,ztgt)		! xxx***xxx
	else
		FF_nuc=1.
	end if
c
	if (ztgt.eq.1) then !use Fourier transform of hydrogen ground state charge density
		FF_atomic=1./( 1.+(aBohr/2.)**2*q**2 )**2	!with FF_atomic(0)=1.
c
	else if (ztgt.eq.82) then !use sum of Gaussian form for atomic FF
C		Reference: http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors.php
C
		qA=q/hbarc*100000.	!qA in Angstrom^-1		! xxx***xxx
		FF_atomic=0.
		do i=1,4
			FF_atomic=FF_atomic+a(i)*exp(-b(i)*(qA/(4.*pi))**2)	! xxx***xxx
		end do
		if (qA.lt.25.) then ! add simple constant term		! xxx***xxx
			FF_atomic=FF_atomic+c	
		else 	! If qA greater than 25 Angstrom^-1, make ad hoc assumption that the constant term 
C			decreases with increasing q^2 at the same rate as the a(1) term.
			FF_atomic=FF_atomic +c*exp(-b(1)*(qA/(4.*pi))**2)*exp(+b(1)*(25./(4.*pi))**2)	! xxx***xxx
		end if
		FF_atomic=FF_atomic/float(ztgt)	!normalize to FF_atomic(0)=1.
	else !all other cases 
		FF_atomic=0.
	end if
c
	FF2=(FF_nuc-FF_atomic)**2
c
	return
	end		
c
c******************************************************************
c
	subroutine density_init(ztgt)
	implicit none
	real*8 pi,c_den(100),a_den(100),rho0,w
	integer*8 i,ztgt
	common /charge_density/ rho0,c_den,a_den
c
	pi=acos(-1.)
c
	rho0=0.
c	These equations have to do with Fermi distribution, reference? , worth rechecking this before using
	w=4.*pi*c_den(ztgt)/3.*( (pi*a_den(ztgt))**2+c_den(ztgt)**2 )
	do i=1,10
		w=w+8.*pi*a_den(ztgt)**3*(-1.)**(i-1)*exp(-float(i)*c_den(ztgt)/a_den(ztgt))/float(i)**3
	end do
	rho0=1./w
	return
	end
c
c******************************************************************

	real*8 function FF(q,ztgt)		! xxx***xxx
	implicit none
	real*8 q02,hbarc,q,qF,gamma,r(12),A(12),rho0,c_den(100),a_den(100),pi,norm,proton_rms	! xxx***xxx
	integer*8 i,ztgt
c
c	data q02 /0.71/    !standard nucleon dipole form factor parameter, GeV^2
	data q02 /0.6608/  !dipole form factor parameter, GeV^2, assuming PDG RMS proton charge radius of .8409 fm
	data hbarc /.197326/
	common /charge_density/ rho0,c_den,a_den
c
c	Density parameters for 208Pb from Frois et al. 
	data R(1),A(1)   /0.1,0.003845/	! R is in units of fm, the A terms are the same as the Qi terms in the Frois et al. 208Pb paper
	data R(2),A(2)   /0.7,0.009724/	! and sum to 1.0
	data R(3),A(3)   /1.6,0.033093/
	data R(4),A(4)   /2.1,0.000120/
	data R(5),A(5)   /2.7,0.083107/
	data R(6),A(6)   /3.5,0.080869/
	data R(7),A(7)   /4.2,0.139957/
	data R(8),A(8)   /5.1,0.260892/
	data R(9),A(9)   /6.0,0.336013/
	data R(10),A(10) /6.6,0.033637/
	data R(11),A(11) /7.6,0.018729/
	data R(12),A(12) /8.7,0.000020/
	data gamma /1.388/	!units of fm
c
c  Select the FF
c
	pi=acos(-1.)
c
	if (ztgt.eq.1) then	!proton, with FF(0)=1	
		FF=1./(1.+q**2/q02)**2	! Standard dipole form factor		! xxx***xxx
c
	else if (ztgt.eq.82) then	!lead, with FF(0)=1
		qF=q/hbarc		!units of fm^-1				! xxx***xxx		
		FF=0.
		do i=1,12
			FF=FF+A(i)*(gamma**2*cos(qF*R(i))+2.*R(i)/qF*sin(qF*R(i)))/(gamma**2+2.*R(i)**2)*		! xxx***xxx
     &			exp(-qF**2*gamma**2/4.)									! xxx***xxx
		end do
c
	else	!for everything else use 2-parameter fermi model, reference Maximon and Schrack, "The Form Factor of the Fermi Model Spatial Distribution",
C		Journal of Research of the National Bureau of Standards - B. Mathematics and Mathematical Physics, Vol. 70B, No. 1, January-March 1966, pg. 85. 
C		The equations here need to be rechecked
c
		FF=4.*pi**2*rho0*a_den(ztgt)**3/( (qF*a_den(ztgt))**2*(sinh(pi*qF*a_den(ztgt)))**2 )*		! xxx***xxx
     &		( pi*qF*a_den(ztgt)*cosh(pi*qF*a_den(ztgt))*sin(qF*c_den(ztgt))-qF*c_den(ztgt)*cos(qF*c_den(ztgt))*sinh(pi*qF*a_den(ztgt)) )	! xxx***xxx
		do i=1,10
			FF=FF+8.*pi*rho0*a_den(ztgt)**3*(-1.)**(i-1)*float(i)*exp(-float(i)*c_den(ztgt)/a_den(ztgt))/
     &			(float(i)**2+(qF*a_den(ztgt))**2)**2							! xxx***xxx
		end do
	end if
c
	return
	end
c
c*********************************************************************
c Calculate internal radiation
c
	real*8 function E_rad(b_rad,m_part,t,E,Eloss)
	implicit none
	real*8 b_rad,t_rad,m_part,t,E,alpha,pi,delta,bt,const,W,W0,W1,u0,u1,Eloss,dudW,ln2
	integer*8 loop_count
	double precision ZBQLUAB,zlo,zhi
	data zlo,zhi /0.,1./
	data alpha,pi,ln2 /7.297352e-3,3.141593,0.693147/
	data delta /.001/	!test for precision on W, the integrated probablity
c
	
c	t_rad=3./4.*alpha/pi*( log( t/m_part**2+2. )-ln2 ) ! Modified the Mo and Tsai/Walker equation so that t_rad=0 at t=0
	t_rad=3./4.*alpha/pi*( log( t/m_part**2+2. )-1. ) ! this is the length as given by Mo and Tsai and Walker PhD
	bt=b_rad*t_rad
c	print *, ' Q2 ', t, ' t_rad ', t_rad,' bt ', bt
c	if (bt.le.0.) then
c		E_rad=E
c		Eloss=0.
c		return
c	end if
	const=1./bt-1./(bt+1.)+3./4./(bt+2.)	! this is the normalization constant for the PDF, see my notes
c
c	***Pick UDRV for integrated probability
10	W=ZBQLUAB(zlo,zhi)
c	***Find approx. solution for radiated energy, see my notes on this
	u0=W**(1./bt)	
c	***Find W0 at u0
	W0=1./const*( u0**bt/bt-u0**(bt+1.)/(bt+1.)+3./4.*u0**(bt+2.)/(bt+2.) )
c	print *, ' W ', W, ' W0 ', W0, ' u0 ', u0 
	loop_count=0	!count how many iterations needed to reach a solution
20	continue
	loop_count=loop_count+1
c	***Test W0:
	if (abs(W0-W).lt.delta.or.u0.le.0.0) then ! reached required precision
		E_rad=(1.-u0)*E
		if (E_rad.le.m_part) go to 10		!too much radiation! try again
		Eloss=u0*E
		return
	end if
c	***Find derivative du/dW at u0, need to exclude the case u0=0.
	dudW=const/( u0**(bt-1.)*(1.-u0+3./4.*u0**2) )
c	***Find new value for u1:
	u1=u0+dudW*(W-W0)	
c	***Find W1 at u1
	W1=1./const*( u1**bt/bt-u1**(bt+1.)/(bt+1.)+3./4.*u1**(bt+2.)/(bt+2.) )
c	print *,' t ', t, ' t_rad ', t_rad,  ' const ', const,' W ', W, ' W1 ', W1, ' u1 ', u1 , ' dudW ', dudW
	u0=u1
	W0=W1
	if (loop_count.ge.30) go to 10 ! Average number of iterations = 3, appears the solution isn't converging, try to radiate again
c 
	go to 20	! Iterate for solution
	end
c

*******************************************************************
********	FILE: randgen.f				***********
********	AUTHORS: Richard Chandler		***********
********		 (richard@stats.ucl.ac.uk)	***********
********		 Paul Northrop 			***********
********		 (northrop@stats.ox.ac.uk)	***********
********	LAST MODIFIED: 26/8/03			***********
********	See file randgen.txt for details	***********
*******************************************************************


******************************************************************
******************************************************************
      SUBROUTINE ZBQLINI(SEED)
******************************************************************
*       To initialize the random number generator - either
*       repeatably or nonrepeatably. Need double precision
*       variables because integer storage can't handle the
*       numbers involved
******************************************************************
*	ARGUMENTS
*	=========
*	SEED	(integer, input). User-input number which generates
*		elements of the array ZBQLIX, which is subsequently used 
*		in the random number generation algorithm. If SEED=0,
*		the array is seeded using the system clock if the 
*		FORTRAN implementation allows it.
******************************************************************
*	PARAMETERS
*	==========
*	LFLNO	(integer). Number of lowest file handle to try when
*		opening a temporary file to copy the system clock into.
*		Default is 80 to keep out of the way of any existing
*		open files (although the program keeps searching till
*		it finds an available handle). If this causes problems,
*               (which will only happen if handles 80 through 99 are 
*               already in use), decrease the default value.
******************************************************************
      INTEGER LFLNO
      PARAMETER (LFLNO=80)
******************************************************************
*	VARIABLES
*	=========
*	SEED	See above
*	ZBQLIX	Seed array for the random number generator. Defined
*		in ZBQLBD01
*	B,C	Used in congruential initialisation of ZBQLIX
*	SS,MM,}	System clock secs, mins, hours and days
*	HH,DD }
*	FILNO	File handle used for temporary file
*	INIT	Indicates whether generator has already been initialised
*
      INTEGER SEED,SS,MM,HH,DD,FILNO,I
      INTEGER INIT
      DOUBLE PRECISION ZBQLIX(43),B,C
      DOUBLE PRECISION TMPVAR1,DSS,DMM,DHH,DDD

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE INIT

*
*	Ensure we don't call this more than once in a program
*
      IF (INIT.GE.1) THEN
       IF(INIT.EQ.1) THEN
        WRITE(*,1)
        INIT = 2
       ENDIF
       RETURN
      ELSE
       INIT = 1
      ENDIF
*
*       If SEED = 0, cat the contents of the clock into a file
*       and transform to obtain ZQBLIX(1), then use a congr.
*       algorithm to set remaining elements. Otherwise take
*       specified value of SEED.
*
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>	NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
*>>>>>>>	(NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
*>>>>>>>	THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
*>>>>>>>	OF THE FOLLOWING IF BLOCK SHOULD BE	 >>>>>>>
*>>>>>>>	COMMENTED OUT.				 >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF (SEED.EQ.0) THEN
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
*>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>
*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       CALL SYSTEM(' date +%S%M%H%j > zbql1234.tmp')
*
*       Try all file numbers for LFLNO to 999 
*
       FILNO = LFLNO
 10    OPEN(FILNO,FILE='zbql1234.tmp',ERR=11)
       GOTO 12
 11    FILNO = FILNO + 1
       IF (FILNO.GT.999) THEN
        WRITE(*,2)
        RETURN
       ENDIF
       GOTO 10
 12    READ(FILNO,'(3(I2),I3)') SS,MM,HH,DD
       CLOSE(FILNO)
       CALL SYSTEM('rm zbql1234.tmp')
       DSS = DINT((DBLE(SS)/6.0D1) * B)
       DMM = DINT((DBLE(MM)/6.0D1) * B)
       DHH = DINT((DBLE(HH)/2.4D1) * B)
       DDD = DINT((DBLE(DD)/3.65D2) * B)
       TMPVAR1 = DMOD(DSS+DMM+DHH+DDD,B)
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<
*<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ELSE
       TMPVAR1 = DMOD(DBLE(SEED),B)
      ENDIF
      ZBQLIX(1) = TMPVAR1
      DO 100 I = 2,43
       TMPVAR1 = ZBQLIX(I-1)*3.0269D4
       TMPVAR1 = DMOD(TMPVAR1,B)       
       ZBQLIX(I) = TMPVAR1
 100  CONTINUE

 1    FORMAT(//5X,'****WARNING**** You have called routine ZBQLINI ',
     +'more than',/5X,'once. I''m ignoring any subsequent calls.',//)
 2    FORMAT(//5X,'**** ERROR **** In routine ZBQLINI, I couldn''t',
     +' find an',/5X,
     +'available file number. To rectify the problem, decrease the ',
     +'value of',/5X,
     +'the parameter LFLNO at the start of this routine (in file ',
     +'randgen.f)',/5X,
     +'and recompile. Any number less than 100 should work.')
      END
******************************************************************
      FUNCTION ZBQLU01(DUMMY)
*
*       Returns a uniform random number between 0 & 1, using
*       a Marsaglia-Zaman type subtract-with-borrow generator.
*       Uses double precision, rather than integer, arithmetic 
*       throughout because MZ's integer constants overflow
*       32-bit integer storage (which goes from -2^31 to 2^31).
*       Ideally, we would explicitly truncate all integer 
*       quantities at each stage to ensure that the double
*       precision representations do not accumulate approximation
*       error; however, on some machines the use of DNINT to
*       accomplish this is *seriously* slow (run-time increased
*       by a factor of about 3). This double precision version 
*       has been tested against an integer implementation that
*       uses long integers (non-standard and, again, slow) -
*       the output was identical up to the 16th decimal place
*       after 10^10 calls, so we're probably OK ...
*
      DOUBLE PRECISION ZBQLU01,DUMMY,B,C,ZBQLIX(43),X,B2,BINV
      INTEGER CURPOS,ID22,ID43

      COMMON /ZBQL0001/ ZBQLIX,B,C
      SAVE /ZBQL0001/
      SAVE CURPOS,ID22,ID43
      DATA CURPOS,ID22,ID43 /1,22,43/

      B2 = B
      BINV = 1.0D0/B
 5    X = ZBQLIX(ID22) - ZBQLIX(ID43) - C
      IF (X.LT.0.0D0) THEN
       X = X + B
       C = 1.0D0
      ELSE
       C = 0.0D0
      ENDIF
      ZBQLIX(ID43) = X
*
*     Update array pointers. Do explicit check for bounds of each to
*     avoid expense of modular arithmetic. If one of them is 0 the others
*     won't be
*
      CURPOS = CURPOS - 1
      ID22 = ID22 - 1
      ID43 = ID43 - 1
      IF (CURPOS.EQ.0) THEN
       CURPOS=43
      ELSEIF (ID22.EQ.0) THEN
       ID22 = 43
      ELSEIF (ID43.EQ.0) THEN
       ID43 = 43
      ENDIF
*
*     The integer arithmetic there can yield X=0, which can cause 
*     problems in subsequent routines (e.g. ZBQLEXP). The problem
*     is simply that X is discrete whereas U is supposed to 
*     be continuous - hence if X is 0, go back and generate another
*     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
*
      IF (X.LT.BINV) THEN
       B2 = B2*B
       GOTO 5
      ENDIF

      ZBQLU01 = X/B2

      END
******************************************************************
      FUNCTION ZBQLUAB(A,B)
*
*       Returns a random number uniformly distributed on (A,B)
*
      DOUBLE PRECISION A,B,ZBQLU01,ZBQLUAB
      
*
*       Even if A > B, this will work as B-A will then be -ve
*
      IF (A.NE.B) THEN
       ZBQLUAB = A + ( (B-A)*ZBQLU01(0.0D0) )
      ELSE
       ZBQLUAB = A
       WRITE(*,1)
      ENDIF
 1    FORMAT(/5X,'****WARNING**** (function ZBQLUAB) Upper and ',
     +'lower limits on uniform',/5X,'distribution are identical',/)
      END
******************************************************************
      FUNCTION ZBQLEXP(MU)
*
*       Returns a random number exponentially distributed with
*       mean MU
*
      DOUBLE PRECISION MU,ZBQLEXP,ZBQLU01

      ZBQLEXP = 0.0D0

      IF (MU.LT.0.0D0) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      ZBQLEXP = -DLOG(ZBQLU01(0.0D0))*MU

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLEXP',/)

      END
******************************************************************
      FUNCTION ZBQLNOR(MU,SIGMA)
*
*       Returns a random number Normally distributed with mean
*       MU and standard deviation |SIGMA|, using the Box-Muller
*       algorithm
*
      DOUBLE PRECISION THETA,R,ZBQLNOR,ZBQLU01,PI,MU,SIGMA
      DOUBLE PRECISION SPARE
      INTEGER STATUS
      SAVE STATUS,SPARE,PI
      DATA STATUS /-1/

      IF (STATUS.EQ.-1) PI = 4.0D0*DATAN(1.0D0)

      IF (STATUS.LE.0) THEN
       THETA = 2.0D0*PI*ZBQLU01(0.0D0)
       R = DSQRT( -2.0D0*DLOG(ZBQLU01(0.0D0)) )
       ZBQLNOR = (R*DCOS(THETA))
       SPARE = (R*DSIN(THETA))
       STATUS = 1
      ELSE
       ZBQLNOR = SPARE
       STATUS = 0
      ENDIF
      
      ZBQLNOR = MU + (SIGMA*ZBQLNOR)

      END
******************************************************************
      FUNCTION ZBQLBIN(N,P)
*
*       Returns a random number binomially distributed (N,P)
*
      DOUBLE PRECISION P,ZBQLBET1
      DOUBLE PRECISION PP,PPP,G,Y,TINY
      INTEGER N,ZBQLBIN,ZBQLGEO,IZ,NN

      TINY = 1.0D-8
      ZBQLBIN = 0

      IF (.NOT.( (P.GE.0.0D0).AND.(P.LE.1.0D0) )) THEN
       WRITE(*,1)
       RETURN
      ELSEIF(N.LE.0) THEN
       WRITE(*,1)
       RETURN
      ENDIF
*
*	First step: if NP > 10, say, things will be expensive, and 
*	we can get into the right ballpark by guessing a value for
*	ZBQLBIN (IZ, say), and simulating Y from a Beta distribution 
*	with parameters IZ and NN-IZ+1 (NN starts off equal to N).
*	If Y is less than PP (which starts off as P) then the IZth order 
*	statistic from NN U(0,1) variates is less than P, and we know 
*	that there are at least IZ successes. In this case we focus on 
*	the remaining (NN-IZ) order statistics and count how many are
*	less than PP, which is binomial (NN-IZ,(PP-Y)/(1-Y)). 
*	Otherwise, if Y is greater than PP there must be less 
*	than IZ successes, so we can count the number of order statistics
*	under PP, which is binomial (IZ-1,P/Y). When we've got NN*PP
*	small enough, we go to the next stage of the algorithm and 
*	generate the final bits directly.
*
      NN = N
      PP = P
 10   IZ = INT(DBLE(NN)*PP) + 1
      IF ( (IZ.GT.10).AND.(IZ.LT.NN-10) ) THEN
       Y = ZBQLBET1(DBLE(IZ),DBLE(NN-IZ+1))
       IF (Y.LT.PP) THEN
        ZBQLBIN = ZBQLBIN + IZ
        NN = NN - IZ
        PP = (PP-Y) / (1.0D0-Y)
       ELSE
        NN = IZ-1
        PP = PP/Y
       ENDIF
       GOTO 10
      ENDIF
*
*	PP is the probability of the binomial we're currently
*	simulating from. For the final part, we simulate either number 
*	of failures or number of success, depending which is cheaper.
*      
 20   IF (PP.GT.0.5) THEN
       PPP = 1.0D0-PP
      ELSE
       PPP = PP
      ENDIF

      G = 0
      IZ = 0
*
*     ZBQLGEO falls over for miniscule values of PPP, so ignore these
*     (tiny probability of any successes in this case, anyway)
* 
      IF (PPP.GT.TINY) THEN
 30    G = G + ZBQLGEO(PPP)
       IF (G.LE.NN) THEN
        IZ = IZ + 1
        GOTO 30
       ENDIF
      ENDIF

      IF (PP.GT.0.5) IZ = NN - IZ
      ZBQLBIN = ZBQLBIN + IZ

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLBIN',/)
      END
******************************************************************
      FUNCTION ZBQLGEO(P)
*
*       Returns a random number geometrically distributed with 
*       parameter P ie. mean 1/P
* 
  
      DOUBLE PRECISION P,ZBQLU01,U,TINY
      INTEGER ZBQLGEO

      TINY = 1.0D-12
      ZBQLGEO = 0
 
      IF (.NOT.( (P.GE.0.0D0).AND.(P.LE.1.0D0) )) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      IF (P.GT.0.9D0) THEN
 10    ZBQLGEO = ZBQLGEO + 1 
       U = ZBQLU01(0.0D0)
       IF (U.GT.P) GOTO 10
      ELSE
       U = ZBQLU01(0.0D0)
*
*	For tiny P, 1-p will be stored inaccurately and log(1-p) may
*	be zero. In this case approximate log(1-p) by -p
*
       IF (P.GT.TINY) THEN
        ZBQLGEO = 1 + INT( DLOG(U)/DLOG(1.0D0-P) )
       ELSE
        ZBQLGEO = 1 + INT(-DLOG(U)/P)
       ENDIF
      ENDIF

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLGEO',/)
      END
******************************************************************
      FUNCTION ZBQLPOI(MU)
*
*       Returns a random number Poisson distributed with mean MU
*
      
      DOUBLE PRECISION ZBQLU01,X,Y,MU,PI
      DOUBLE PRECISION ZBQLLG,ZBQLGAM,MU1,TMP1,TMP2,T
      INTEGER ZBQLPOI,ZBQLBIN,K,INIT
      SAVE INIT,PI
      DATA INIT /0/

      IF (INIT.EQ.0) THEN
       PI = 4.0D0*DATAN(1.0D0)
       INIT = 1
      ENDIF

      ZBQLPOI = 0

      IF (MU.LT.0.0D0) THEN
       WRITE(*,1)
       RETURN
      ENDIF
*
*      For small MU, generate exponentials till their sum exceeds 1
*      (equivalently, uniforms till their product falls below e^-MU)
*
      IF (MU.LE.1.0D3) THEN
       MU1 = MU
*
*     For values of MU less than 1000, use order statistics - the Kth
*     event in a Poisson process of rate MU has a Gamma distribution
*     with parameters (MU,K); if it's greater than 1 we know that there 
*     are less than K events in (0,1) (and the exact number is binomial)
*     and otherwise the remaining number is another Poisson. Choose K so
*     that we'll get pretty close to 1 in the first go but are unlikely
*     to overshoot it.
*
 19    IF (MU1.GT.1.0D1) THEN        
        K = INT(MU1-DSQRT(MU1))
        Y = ZBQLGAM(DBLE(K),MU1)
        IF (Y.GT.1.0D0) THEN
         ZBQLPOI = ZBQLPOI + ZBQLBIN(K-1,(1.0D0/Y))
         RETURN
        ENDIF
        ZBQLPOI = ZBQLPOI + K
        MU1 = MU1  * (1.0D0-Y)
        GOTO 19
       ENDIF
       Y = DEXP(-MU1)
       X = 1.0D0
 20    X = X*ZBQLU01(0.0D0)
       IF (X.GT.Y) THEN
        ZBQLPOI = ZBQLPOI + 1
        GOTO 20
       ENDIF
*
*     For really huge values of MU, use rejection sampling as in 
*     Press et al (1992) - large numbers mean some accuracy may be
*     lost, but it caps the execution time.
*
      ELSE
       TMP1 = DSQRT(2.0D0*MU)
       TMP2 = ZBQLLG(MU+1.0D0)-(MU*DLOG(MU))
 30    Y = DTAN(PI*ZBQLU01(0.0D0))
       ZBQLPOI = INT(MU + (TMP1*Y))
       IF (ZBQLPOI.LT.0) GOTO 30
       X = DBLE(ZBQLPOI)
       T = (X*DLOG(MU)-ZBQLLG(X+1.0D0)) + TMP2
       IF (DABS(T).LT.1.0D2) THEN
        T = 0.9D0*(1.0D0+(Y*Y))*DEXP(T)
        IF (ZBQLU01(0.0D0).GT.T) GOTO 30
       ELSE
        T = DLOG(0.9D0*(1.0D0+(Y*Y))) + T
        IF (DLOG(ZBQLU01(0.0D0)).GT.T) GOTO 30
       ENDIF
      ENDIF 

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLPOI',/)
      END
******************************************************************
      FUNCTION ZBQLGAM(G,H)
*
*       Returns a random number with a gamma distribution with mean
*       G/H and variance G/(H^2). (ie. shape parameter G & scale
*       parameter H)
*
      DOUBLE PRECISION C,D,R,ZBQLGAM,ZBQLU01,G,H,A,z1,z2,B1,B2,M
      DOUBLE PRECISION U1,U2,U,V,TEST,X
      double precision c1,c2,c3,c4,c5,w

      ZBQLGAM = 0.0D0

      IF ( (G.LE.0.0D0).OR.(H.LT.0.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      IF (G.LT.1.0D0) THEN
889    u=ZBQLU01(0.0d0)
       v=ZBQLU01(0.0d0)
       if (u.gt.exp(1.0d0)/(g+exp(1.0d0))) goto 891
       ZBQLGAM=((g+exp(1.0d0))*u/exp(1.0d0))**(1.0d0/g)
       if (v.gt.exp(-ZBQLGAM)) then
	goto 889
       else
	goto 892
       endif
891    ZBQLGAM=-log((g+exp(1.0d0))*(1.0d0-u)/(g*exp(1.0d0)))
       if (v.gt.ZBQLGAM**(g-1.0)) goto 889
892    ZBQLGAM=ZBQLGAM/h
       RETURN
      ELSEIF (G.LT.2.0D0) THEN
       M = 0.0D0
      elseif (g.gt.10.0d0) then
       c1=g-1.0d0
       c2=(g-1.0d0/(6.0d0*g))/c1
       c3=2.0d0/c1
       c4=c3+2.0d0
       c5=1.0d0/sqrt(g)
777    u=ZBQLU01(0.0d0)
       v=ZBQLU01(0.0d0)
       if (g.gt.2.50d0) then
	u=v+c5*(1.0d0-1.860d0*u)
       endif 
       if (u.le.0.0d0.or.u.ge.1.0d0) goto 777 
       w=c2*v/u 
       if (c3*u+w+1.0d0/w.le.c4) goto 778 
       if (c3*log(u)-log(w)+w.ge.1.0d0) goto 777
778    ZBQLGAM=c1*w/h 
       return
      ELSE
       M = -(G-2.0D0) 
      ENDIF
      R = 0.50D0
      a = ((g-1.0d0)/exp(1.0d0))**((g-1.0d0)/(r+1.0d0))
      C = (R*(M+G)+1.0D0)/(2.0D0*R)
      D = M*(R+1.0D0)/R
      z1 = C-DSQRT(C*C-D)
*
*     On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
*     is not exactly zero - this needs trapping if negative.
*
      IF ((Z1-M.LT.0.0D0).AND.(Z1-M.GT.-1.0D-12)) Z1 = M
      z2 = C+DSQRT(C*C-D)
      B1=(z1*(z1-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z1-M)/(R+1.0D0))
      B2=(z2*(z2-M)**(R*(G-1.0D0)/(R+1.0D0)))*DEXP(-R*(z2-M)/(R+1.0D0))
50    U1=ZBQLU01(0.0D0)
      U2=ZBQLU01(0.0D0)
      U=A*U1
      V=B1+(B2-B1)*U2
      X=V/(U**R)
      IF (X.LE.M) GOTO 50
      TEST = ((X-M)**((G-1)/(R+1)))*EXP(-(X-M)/(R+1.0D0))
      IF (U.LE.TEST) THEN
       ZBQLGAM = (X-M)/H
      ELSE
       GOTO 50
      ENDIF
 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLGAM',/5X, '(both parameters must be positive)',/)

      END
***************************************************************
      FUNCTION ZBQLBET1(NU1,NU2)
*
*       Returns a random number, beta distributed with degrees
*       of freedom NU1 and NU2.
*
      DOUBLE PRECISION NU1,NU2,ZBQLGAM,ZBQLBET1,ZBQLU01,X1,X2

      ZBQLBET1 = 0.0D0

      IF ( (NU1.LE.0.0).OR.(NU2.LE.0.0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF
*      
*       If parameters are too small, gamma subroutine tends to return zero
*       as all the probability goes to the origin and we get rounding
*       errors, even with double precision. In this case, we use Johnk's
*       method, suitably scaled to avoid rounding errors as much as possible.
*
      
      IF ( (NU1.LT.0.9D0).AND.(NU2.LT.0.9D0) ) THEN
 10    X1 = ZBQLU01(0.0D0)
       X2 = ZBQLU01(0.0D0)
       IF ( (X1**(1.0D0/NU1))+(X2**(1.0D0/NU2)).GT.1.0D0) GOTO 10    
       X1 = (DLOG(X2)/NU2) - (DLOG(X1)/NU1)
       ZBQLBET1 = (1.0D0 + DEXP(X1))**(-1)
       IF (ZBQLBET1.GT.1.0D0) GOTO 10    
      ELSE
       X1 = ZBQLGAM(NU1,1.0D0)
       X2 = ZBQLGAM(NU2,1.0D0)
       ZBQLBET1 = X1/(X1+X2)
      ENDIF
       
      RETURN

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLBET1',/5X,
     +'(both degrees of freedom must be positive)',/)

      END
***************************************************************
      FUNCTION ZBQLWEI(A,B)
*
*       Returns a random number, Weibull distributed with shape parameter
*       A and location parameter B, i.e. density is
*	f(x) = ( A/(B**A) ) * x**(A-1) * EXP( -(x/B)**A )
*
      DOUBLE PRECISION A,B,ZBQLU01,ZBQLWEI,U

      ZBQLWEI = 0.0D0

      IF ( (A.LE.0.0).OR.(B.LE.0.0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF
 
      U = ZBQLU01(0.0D0)
      ZBQLWEI = B * ( (-DLOG(U))**(1.0D0/A) )

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLWEI',/5X,
     +'(both parameters must be positive)',/)
      END
***************************************************************
      FUNCTION ZBQLNB(R,P)
*
*       Returns a pseudo-random number according to a Negative
*	Binomial distribution with parameters (R,P). NB these are
*	both DOUBLE - it copes with non-integer R as well. The
*       form of the distribution is *not* the no. of trials to 
*       the Rth success - see documentation for full spec.
*
      DOUBLE PRECISION R,P,ZBQLGAM,Y
      INTEGER ZBQLNB,ZBQLPOI

      ZBQLNB = 0

      IF ( (R.LE.0.0D0).OR.(P.LE.0.0D0).OR.(P.GE.1.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF

      Y = ZBQLGAM(R,1.0D0)
      Y = Y*P/(1.0D0-P)
      ZBQLNB = ZBQLPOI(Y)

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLNB')
      END
***************************************************************
      FUNCTION ZBQLPAR(A,B)
*
*     Returns a random number, Pareto distributed with parameters
*     A and B. The density is A*(B**A) / (B+X)**(A+1) for X > 0.
*     (this is slightly nonstandard - see documentation in 
*     randgen.txt). The algorithm is straightforward - it uses the
*     inverse CDF method.
*
      DOUBLE PRECISION A,B,ZBQLPAR,ZBQLU01,U

      ZBQLPAR = 0.0D0

      IF ( (A.LE.0.0D0).OR.(B.LE.0.0D0) ) THEN
       WRITE(*,1)
       RETURN
      ENDIF
 
      U = ZBQLU01(0.0D0)
      ZBQLPAR = B * (U**(-1.0D0/A)-1.0D0)

 1    FORMAT(/5X,'****ERROR**** Illegal parameter value in ',
     +' ZBQLPAR',/5X,
     +'(both parameters must be positive)',/)
      END
***************************************************************
      FUNCTION ZBQLLG(X)
*
*     Returns log(G(X)) where G is the Gamma function. The algorithm is
*     that given in Press et al (1992), Section 6.1, although this
*     version also allows for arguments less than 1.
*
      DOUBLE PRECISION X,Z,Z2,ZBQLLG,PI,RLN2P,C(0:6),TMP,SUM
      INTEGER INIT,I
      SAVE INIT,C,RLN2P,PI
      DATA INIT /0/
      DATA (C(I),I=0,6) / 
     +              1.000000000190015D0,76.18009172947146D0,
     +              -86.50532032941677D0,24.01409824083091D0,
     +              -1.231739572450155D0,0.1208650973866179D-2,
     +              -0.5395239384953D-5/

      IF (INIT.EQ.0) THEN
        PI = 4.0D0*DATAN(1.0D0)
        RLN2P = 0.5D0*DLOG(2.0D0*PI)
        INIT = 1
      ENDIF
*
*     Compute for x > 1, then use transformation if necessary. Z is
*     our working argument.
*
      IF (X.GE.1.0D0) THEN
       Z = X
      ELSE 
       Z = 2.0D0-X
       Z2 = 1.0D0-X
      ENDIF

      IF (DABS(Z-1.0D0).LT.1.0D-12) THEN
       ZBQLLG = 0.0D0
       RETURN
      ENDIF

      TMP = Z + 4.5D0
      TMP = ( (Z-0.5D0)*DLOG(TMP) ) - TMP + RLN2P

      SUM = C(0)
      DO 50 I=1,6
       SUM = SUM + (C(I)/(Z+DBLE(I-1)))
 50   CONTINUE
      ZBQLLG = TMP + DLOG(SUM)
*
*     Transformation required if X<1
*
      IF (X.LT.1.0D0) THEN
       TMP = PI*Z2
       ZBQLLG = DLOG(TMP/DSIN(TMP)) - ZBQLLG
      ENDIF


      END

c*******************************************************************
c	Binary search subroutine for variable-width energy bins
c	Finds which bin contains the given energy E
c*******************************************************************
	subroutine find_bin(E, E_low, E_high, n, ibin)
	implicit none
	integer*8 n, ibin, left, right, mid
	real*8 E, E_low(n), E_high(n)
	
	! Initialize binary search
	left = 1
	right = n
	ibin = 0
	
	! Binary search loop
	do while (left <= right)
		mid = (left + right) / 2
		if (E >= E_low(mid) .and. E < E_high(mid)) then
			ibin = mid
			return
		else if (E < E_low(mid)) then
			right = mid - 1
		else
			left = mid + 1
		endif
	enddo
	
	! Check last bin (inclusive upper edge)
	if (E >= E_low(n) .and. E <= E_high(n)) then
		ibin = n
	endif
	
	return
	end

c*******************************************************************
c	Validate user energy limits against bremsstrahlung file limits
c	Reads metadata from varbin cobrems file and checks compatibility
c*******************************************************************
	subroutine validate_cobrems_limits(unit_num, E_lo_user, E_hi_user)
	implicit none
	integer unit_num
	real*8 E_lo_user, E_hi_user, E_min_avail, E_max_avail
	character*200 line
	logical found_min, found_max
	
	found_min = .false.
	found_max = .false.
	
	! Rewind file to start
	rewind(unit_num)
	
	! Read through comments to find energy limits
	do while (.true.)
		read(unit_num, '(A)', end=100) line
		
		! Skip non-comment lines (start reading data)
		if (line(1:1) /= '#') then
			backspace(unit_num)
			exit
		endif
		
		! Look for energy limit metadata
		if (index(line, 'E_MIN_AVAILABLE') > 0) then
			read(line(index(line, 'E_MIN_AVAILABLE')+15:), *) E_min_avail
			found_min = .true.
		endif
		
		if (index(line, 'E_MAX_AVAILABLE') > 0) then
			read(line(index(line, 'E_MAX_AVAILABLE')+15:), *) E_max_avail  
			found_max = .true.
		endif
	enddo
	
100	continue
	
	! Perform validation if limits were found
	if (found_min .and. found_max) then
		print *, '=================================================='
		print *, 'BREMSSTRAHLUNG FILE VALIDATION:'
		print '(A,F8.4,A,F8.4,A)', ' File energy range: ', E_min_avail, ' -', E_max_avail, ' GeV'
		print '(A,F8.4,A,F8.4,A)', ' User energy range: ', E_lo_user, ' -', E_hi_user, ' GeV'
		
		if (E_lo_user < E_min_avail) then
			print *, 'ERROR: User E_LO is below file minimum!'
			print '(A,F8.4,A,F8.4)', ' E_LO = ', E_lo_user, ' < E_MIN_AVAILABLE = ', E_min_avail
			print *, 'Set E_LO >= ', E_min_avail, ' GeV or use different cobrems file'
			stop 'ENERGY RANGE ERROR'
		endif
		
		if (E_hi_user > E_max_avail) then
			print *, 'ERROR: User E_HI is above file maximum!'
			print '(A,F8.4,A,F8.4)', ' E_HI = ', E_hi_user, ' > E_MAX_AVAILABLE = ', E_max_avail
			print *, 'Set E_HI <= ', E_max_avail, ' GeV or use different cobrems file'
			stop 'ENERGY RANGE ERROR'
		endif
		
		print *, 'Energy ranges compatible - proceeding with generation'
		print *, '=================================================='
	else
		print *, 'WARNING: Could not read energy limits from cobrems file'
		print *, 'Proceeding without validation (may cause runtime errors)'
	endif
	
	return
	end

c*******************************************************************
c	Smart bremsstrahlung file selection with future flexibility
c	Handles defaults, overrides, and future polarization-specific files
c*******************************************************************
	subroutine select_cobrems_files(GlueX, cobrems_varbin, 
     &    cobremsGlueX, cobremsCPP, cobremsGlueX_var, cobremsCPP_var, templateDir)
	implicit none
	logical GlueX, cobrems_varbin
	character*200 cobremsGlueX, cobremsCPP, cobremsGlueX_var, cobremsCPP_var
	character*200 templateDir, override_file
	character*200 pol_suffix, exp_name
	
	! Get user overrides from regex replacements
	! This will be a relative path like ./CobremsDistribution_12345_0deg_varbin.dat
	override_file = 'RBHGCOBREMS_FILE_OVERRIDE'
	
	! Determine experiment name and polarization suffix for future use
	if (GlueX) then
		exp_name = 'GlueX'
	else
		exp_name = 'CPP'
	endif
	
	! Future: Add polarization-specific suffix logic here
	! Example: pol_suffix = '_0DEG', '_45DEG', '_90DEG', '_135DEG', '_AMO'
	pol_suffix = ''  ! For now, no polarization-specific files
	
	! Select appropriate files
	! Check if override_file is not empty and not the original placeholder
	if (len_trim(override_file) > 0 .and. 
     &      index(override_file, 'OVERRIDE') == 0 .and.
     &      trim(override_file) /= '#') then
		! User specified custom file (copied to templateDir/CobremsDistributions/)
		! Prepend templateDir/CobremsDistributions/ to the filename
		cobremsGlueX = trim(templateDir) // '/CobremsDistributions/' // trim(override_file)
		cobremsCPP = trim(templateDir) // '/CobremsDistributions/' // trim(override_file)
		cobremsGlueX_var = trim(templateDir) // '/CobremsDistributions/' // trim(override_file)
		cobremsCPP_var = trim(templateDir) // '/CobremsDistributions/' // trim(override_file)
		print *, 'Using custom cobrems: ', trim(override_file)
	else
		! Use smart defaults
		if (GlueX) then
			if (cobrems_varbin) then
				! GlueX variable-width (future: add polarization suffix)
				cobremsGlueX_var = trim(templateDir) // '/CobremsDistribution_varbin' // trim(pol_suffix) // '.dat'
				print *, 'Using GlueX variable-width cobrems (future: polarization-aware)'
			else
				! GlueX fixed-width (legacy compatibility)
				cobremsGlueX = trim(templateDir) // '/CobremsDistribution.dat'
				print *, 'Using legacy GlueX fixed-width cobrems'
			endif
		else
			if (cobrems_varbin) then
				! CPP variable-width default (future: add polarization suffix)
				cobremsCPP_var = trim(templateDir) // '/CPP_Cobrems_varbin' // trim(pol_suffix) // '.dat'
				print *, 'Using CPP variable-width cobrems (future: polarization-aware)'
			else
				! CPP fixed-width (legacy)
				cobremsCPP = trim(templateDir) // '/CPP_Cobrems.dat'
				print *, 'Using legacy CPP fixed-width cobrems'
			endif
		endif
	endif
	
	return
	end



	



				


