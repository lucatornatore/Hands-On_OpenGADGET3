#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  enable/Disable compile-time options as needed #
##################################################

#-------------------------------------- Some General settings
VERBOSE_LEVEL=2
EXPERT_LEVEL=42
#DEBUG
#DO_NOT_SPLIT_OVERMASSIVE_GAS

#-------------------------------------- Suggested settings for special behaviour
KD_LIMIT_DIVV=30
#SUBFIND_DO_NOT_COUNT_IN_PARALLEL
KD_FAST_SINGLE_GUESS
NOTEST_FOR_IDUNIQUENESS
KD_RESTRICT_NEIGHBOURS
KD_BUFFER_MANAGEMENT=0.3
KD_RESTRICT_MAXVEL=5000.0
WRITE_INFO_BLOCK               # Enables writing the INFO block
#WRITE_KEY_FILES=1
SYNCRONIZ_OUTPUT
NOSTOP_WHEN_BELOW_MINTIMESTEP

#--------------------------------------- Special settings for this setup
#DS_SHIFT_BOX
#LONGIDS
PERIODIC
MULTIPLEDOMAINS=8
TOPNODEFACTOR=3
DOUBLEPRECISION=3

#--------------------------------------- Basic operation mode of code
COOLING
SFR

#--------------------------------------- SPH Settings
WENDLAND_C6_KERNEL        # Switch to Wendland C6 kernel

#--------------------------------------- SPH viscosity options
#TIME_DEP_ART_VISC=2           # Enables time dependend viscosity
#ARTIFICIAL_CONDUCTIVITY    # enables Price-Monaghan artificial conductivity
#TIME_DEP_ART_COND=2

#--------------------------------------- TreePM Options
PMGRID=256
KD_PM_NOSENDPOS
KD_PM_FASTGREENS
ASMTH=1.25
RCUT=4.0
#PLACEHIGHRESREGION=55
#ENLARGEREGION=1.1
#HIGHRES_PMGRID=256
EVALPOTENTIAL           # computes gravitational potential            

#--------------------------------------- SFR/feedback model
STELLARAGE
WINDS

#------------------------------------ # [options from LT and DF]
LT_STELLAREVOLUTION         # enable stellar evolution+

LT_METAL_COOLING_WAL

LT_NMet=16                  # number of species
LT_PM_LIFETIMES             # use Padovani&matteucci 1993 lifetimes

LT_MAX_TEMP_FEEDBACK=5.0e8
LT_USE_TRUENGB
LT_STARS_GUESSHSML
LT_WIND_VELOCITY=350.0     # set the winds' velocity in km/s
ALL_EXCLUDE_WINDS_FROM_COOLING
ALL_EXCLUDE_WINDS_FROM_FEEDBACK

#-------------------------------------- AGN stuff
#BLACK_HOLES             # enables Black-Holes (master switch)
#AD_DYNFRIC

#--------------------------------------- Thermal conduction
#CONDUCTION
#CONDUCTION_SATURATION

#---------------------------------------- On the fly FOF groupfinder
#FOF                                # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32   # 2^type for the types linked to nearest primaries                                                             
#FOF_GROUP_MIN_LEN=32              # default is 32
#SUBFIND                            # enables substructure finder
#DENSITY_SPLIT_BY_TYPE=1+2+16+32    # 2^type for whch the densities should be calculated seperately
#KD_ALTERNATIVE_GROUP_SORT
#MAX_NGB_CHECK=3                      		# Max number of neighbours for saddle-point detection (default = 2)
#SUBFIND_SAVE_PARTICLELISTS          		# saves also phase-space and type variables parallel to IDs
#SO_VEL_DISPERSIONS                  		# computes velocity dispersions for part of FOF SO-properties
#KEEP_HSML_AS_GUESS                   		# keeps using hsml for gas particles in subfind_density
#LINKLENGTH=0.16                      		# Linkinglength for FoF (default=0.2)
#NO_GAS_CLOUDS                        		# does not accept pure gaseous sub-structures
#WRITE_SUB_IN_SNAP_FORMAT             		# saves subfind results in snap format
#WRITE_SUB_COLLECTIVELY               		# combines several subfind outputs from MPI tasks
#LT_ADD_GAL_TO_SUB=12                 		# adds optical luminosities in 6 bands to subhalos
#INTERP_OBSERVER_FRAME=144
#DUSTATT=11                           		# includes dust attenuation into the luminosity calculation (using 11 radial bins)
#OBSERVER_FRAME                       		# if defined, use CB07 Observer Frame Luminosities, otherwise CB07 Rest Frame Luminosities
#SO_BAR_INFO                          		# adds temperature, Lx, bfrac, etc to Groups
#SUBFIND_COUNT_BIG_HALOS=1e4          		# adds extra blocks for halos with M_TopHat > SUBFIND_COUNT_BIG_HALOS
#NO_SAVE_OF_HSML
#NO_FOF_OUTPUT                        		# does not output the group tabs in separate files
#NEW_FILTERS
#KD_MAIN_HALO
#KD_ID_LIST_OPTIMIZATION
#SUBFIND_ID_OVERBOOK_FACTOR=20
#KD_EXTRA_TREEALLOC_FACTOR_SUBFIND=2.0

