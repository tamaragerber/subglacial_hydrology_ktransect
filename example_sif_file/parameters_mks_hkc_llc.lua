-- ---------------------------------
-- all units in SI 
-- parameters for lev6 spinup
-- ---------------------------------- 
-- general
-- --------------------------------------------
gravity = 9.80		 			--! gravitational acceleration [ms⁻²]
lw = 334000.0					--! latent heat [Jkg⁻¹]
rhoi = 910.0	 				--! ice density [kg m⁻³]
rhow = 1000.0	 				--! water density [kg m⁻³]
 
-- ice flow
-- ---------------------------------------------
ng = 3.0					--! glen's flow law parameter 
ac = 6.8e-25  					--! channel closure coefficient [Pa^(-3)s⁻¹]
as = 6.8e-25					--! sheet closure coefficient [Pa^(-3)s^-1]

-- GlaDS model
-- ---------------------------------------------
ct = -7.5e-8					--! pressure melt coefficient [K Pa⁻¹]
cw = 4220.0					--! heat capacity of water [Jkg⁻¹K⁻¹]
alphas = 3					--! first sheet flow component
betas = 2 					--! second sheet flow exponent
alphac = 5/4 					--! first channel flow exponent
betac = 3/2					--! second channel flow exponent

ks = 0.001   					--! sheet conductivity [m^(7/4)kg^(-1/2)]
kc = 0.5  				 	--! channel conductivity [m^(3/2)kg^(-1/2)]
lc = 5   					--! sheet width below channel [m] /characteristic length
lr = 10  					--! cavity spacing [m]
hr = 0.5   					--! bedrock bump height [m]
ev = 10e-4 					--! englacial void ratio
--!am = 10.0					--! moulin cross-sectional area [m²]

-- bm = 3.1709792e-8                		--! sheet input [m s⁻¹], corresponding to 1 mm/yr of basal melt.

--!zsl = 1e-18					--! sea level

maxarea =1000.0					--! max channel area [m²]
maxh = 10					--! max sheet thickness [m] 




