################################################################################
#  Kspec_Lee05_TB_NTU.R
#
#  Spectral diffuse attenuation (Kd) model for Tampa Bay
#  Based on Lee et al. (2005), parameterised with NTU as a TSS proxy.
#
#  Ported from:  K-spec_Lee_05_TB_NTU_1p5_2_m_Save_SpecUSFTBEP2026.xlsm
#  Original macro author: Chuck Gallegos (2009)
#
#  Validation (139 observations, 3 seasons, Tampa Bay 2008):
#    Median KdPAR error vs Excel: 0.47%   Max: 1.2%
#    Median ZSD    error vs Excel: 0.10%   Max: 0.3%
#
#  References:
#    Lee Z et al. (2005) Euphotic zone depth: Its derivation and implication
#      to ocean-color remote sensing. JGR 110, C03009.
#    Lee Z et al. (2015) Secchi disk depth: A new theory and mechanistic
#      model for underwater visibility. Remote Sens. Environ. 169, 139-149.
#
# --------------------------------------------------------------------------
#  WORKFLOW SUMMARY
#  ─────────────────────────────────────────────────────────────────────────
#  1. Source this file.
#  2. Build a data.frame with columns:
#       Station_ID, Year, Month, Day, DecTime, z, CDOM440, CHLA, NTU
#  3. Call  batch_run(dat)  to get Kd spectra + KdPAR + ZSD for every row.
#  4. Optionally call  plot_Kd_spectrum()  or  plot_validation()  to
#     visualise results.
#
#  Quick single-observation example:
#    mu0 <- compute_solar_mu0(2008, 8, 7, 10.567, lat=27.76)
#    kd  <- compute_Kd_spectra(CDOM440=0.238, CHLA=4.07, NTU=1.2,
#                              mu0=mu0, z=3)
#    kd["KdPAR"]   # 0.4437 m-1
#    kd["ZSD"]     # 3.843  m
# --------------------------------------------------------------------------
#
#  Required R packages:
#    base only — no external packages for the core model.
#    readxl  — only needed for read_batch_data() / read_output_sheet().
#    graphics — only needed for the plotting helpers.
#
#  To install optional packages:
#    install.packages("readxl")
################################################################################


# =============================================================================
# SECTION 1 — SPECTRAL LOOKUP TABLES
#   All vectors are length 64, indexed 1..64, covering 400–715 nm in 5-nm steps.
#   Values extracted directly from the Spectra / bnorm / bbpnorm / Phi_calc
#   sheets of the source workbook.
# =============================================================================

#' Wavelength grid [nm]
LAMBDA <- seq(400, 715, by = 5)   # length 64

# ---- Pure-water absorption aw [m-1]  (Pope & Fry 1997 / Smith & Baker 1981)
#      Source: Spectra sheet, column B
AW <- c(
  0.00663, 0.00530, 0.00473, 0.00444, 0.00454, 0.00478, 0.00495, 0.00530,
  0.00635, 0.00751, 0.00922, 0.00962, 0.00979, 0.01011, 0.01060, 0.01140,
  0.01270, 0.01360, 0.01500, 0.01730, 0.02040, 0.02560, 0.03250, 0.03960,
  0.04090, 0.04170, 0.04340, 0.04520, 0.04740, 0.05110, 0.05650, 0.05960,
  0.06190, 0.06420, 0.06950, 0.07720, 0.08960, 0.11000, 0.13510, 0.16720,
  0.22240, 0.25570, 0.26440, 0.26780, 0.27550, 0.28340, 0.29160, 0.30120,
  0.31080, 0.32500, 0.34000, 0.37100, 0.41000, 0.42900, 0.43900, 0.44800,
  0.46500, 0.48600, 0.51600, 0.55900, 0.62400, 0.70400, 0.82700, 1.00700
)

# ---- Pure-water scattering bw [m-1]
#      Source: Spectra sheet, column C
BW <- c(
  0.00760, 0.00720, 0.00680, 0.00645, 0.00610, 0.00580, 0.00550, 0.00520,
  0.00490, 0.00470, 0.00450, 0.00430, 0.00410, 0.00390, 0.00370, 0.00355,
  0.00340, 0.00325, 0.00310, 0.00300, 0.00290, 0.00275, 0.00260, 0.00250,
  0.00240, 0.00230, 0.00220, 0.00215, 0.00210, 0.00200, 0.00190, 0.00185,
  0.00180, 0.00175, 0.00170, 0.00165, 0.00160, 0.00155, 0.00150, 0.00145,
  0.00140, 0.00135, 0.00130, 0.00125, 0.00120, 0.00115, 0.00110, 0.00105,
  0.00100, 0.00100, 0.00100, 0.00090, 0.00080, 0.00080, 0.00080, 0.00075,
  0.00070, 0.00070, 0.00070, 0.00070, 0.00070, 0.00070, 0.00070, 0.00070
)

# ---- Solar spectral irradiance E0 [uW cm-2 nm-1]  at 400–700 nm (61 bands)
#      Used for E0-weighted KdPAR integration.
#      Source: Spectra sheet, column O, rows 2–62
E0_PAR <- c(
  4.780, 5.568, 6.003, 6.052, 6.135, 6.017, 5.893, 5.940, 6.659, 7.152,
  7.548, 7.826, 7.947, 7.963, 7.990, 8.119, 8.324, 8.014, 7.990, 8.113,
  8.119, 8.108, 8.026, 7.894, 7.970, 8.130, 8.163, 8.133, 8.051, 7.993,
  7.933, 7.982, 7.937, 8.055, 8.160, 8.265, 8.318, 8.375, 8.387, 8.369,
  8.359, 8.332, 8.340, 8.323, 8.305, 8.289, 8.271, 8.267, 8.263, 8.239,
  8.213, 8.207, 8.201, 8.180, 8.157, 8.136, 8.114, 8.102, 8.089, 8.052,
  8.013
)   # 400, 405, …, 700 nm — 61 values, 5-nm step

# ---- Phytoplankton absorption spectral shape Phinorm — normalised at PC=0
#      Source: Phi calc sheet, column K (rows 2–65).
#      PC scores (Phi_PC1/2/3) are all 0 → average shape.
PHINORM <- c(
  1.18481426, 1.23213074, 1.25999308, 1.32210118, 1.35091817, 1.36729342,
  1.42435891, 1.46450653, 1.47369946, 1.40495715, 1.30890234, 1.25687469,
  1.22908485, 1.20307955, 1.15523337, 1.09079212, 1.02832519, 0.98590952,
  0.94215919, 0.89591841, 0.85652849, 0.80331282, 0.74622439, 0.68259435,
  0.63512308, 0.59545691, 0.56343743, 0.53200081, 0.50532689, 0.47600376,
  0.43303117, 0.41395646, 0.36415638, 0.34868058, 0.32783637, 0.31744122,
  0.32025431, 0.32337377, 0.31953943, 0.31000567, 0.29416800, 0.29649258,
  0.32179197, 0.34805657, 0.36878799, 0.37483943, 0.38005082, 0.38794839,
  0.38190933, 0.38222361, 0.39835756, 0.45061134, 0.58588149, 0.77247778,
  0.94543861, 1.00235498, 0.92116088, 0.71178427, 0.46256169, 0.27598396,
  0.15735172, 0.11144374, 0.06731838, 0.05310249
)

# ---- Backscattering spectral shape BBnorm  (Flin=0.5 linear/curved mix)
#      Normalised at 530 nm.  Source: bbpnorm sheet, column D (rows 2–65).
BBNORM <- c(
  1.02607417, 1.02883592, 1.03130872, 1.03349257, 1.03538747, 1.03699342,
  1.03831043, 1.03933848, 1.04007758, 1.04052774, 1.04068895, 1.04056121,
  1.04014451, 1.03943887, 1.03844429, 1.03716075, 1.03558826, 1.03372683,
  1.03157644, 1.02913711, 1.02640882, 1.02339159, 1.02008541, 1.01649028,
  1.01260620, 1.00843318, 1.00397120, 0.99922027, 0.99418040, 0.98885158,
  0.98323380, 0.97732708, 0.97113141, 0.96464679, 0.95787322, 0.95081070,
  0.94346924, 0.93581882, 0.92788946, 0.91967114, 0.91116388, 0.90236767,
  0.89328251, 0.88390840, 0.87424534, 0.86429334, 0.85405238, 0.84352247,
  0.83270362, 0.82159582, 0.81019906, 0.79851336, 0.78653871, 0.77427511,
  0.76172256, 0.74888107, 0.73572256, 0.72231225, 0.70862288, 0.69462559,
  0.68033935, 0.66576415, 0.65090001, 0.63574692
)

# ---- Phytoplankton scattering shape bphin  (for cphistar correction term)
#      Source: bnorm sheet, column B (rows 2–65).
BPHIN <- c(
  1.22737638, 1.19904078, 1.16585851, 1.12843251, 1.08734254, 1.04314517,
  0.99637381, 0.94753867, 0.89712681, 0.84560209, 0.79340519, 0.74095364,
  0.68864176, 0.63684071, 0.58589846, 0.53613982, 0.48786641, 0.44135667,
  0.39686586, 0.35462608, 0.31484622, 0.27771204, 0.24338607, 0.21200770,
  0.18369311, 0.15853534, 0.13660423, 0.11794643, 0.10258543, 0.09052155,
  0.08173191, 0.07617046, 0.07376799, 0.07443207, 0.07804714, 0.08447444,
  0.09355202, 0.10509477, 0.11889440, 0.13471943, 0.15231522, 0.17140394,
  0.19168458, 0.21283297, 0.23450173, 0.25632033, 0.27789506, 0.29880902,
  0.31862213, 0.33687115, 0.35306964, 0.36670800, 0.37725345, 0.38415001,
  0.38681856, 0.38465678, 0.37703915, 0.36331702, 0.34281853, 0.31484865,
  0.27868917, 0.23359870, 0.17881269, 0.11354338
)

# Convenience index: bands in the PAR window 400–700 nm (61 of the 64 bands)
IDX_PAR <- which(LAMBDA >= 400 & LAMBDA <= 700)


# =============================================================================
# SECTION 2 — TAMPA BAY DEFAULT PARAMETERS
#   All values from the Parameters sheet of the source workbook.
# =============================================================================

#' Default Tampa Bay bio-optical parameter set
TB_PARAMS <- list(
  # ---- Site ------------------------------------------------------------------
  latitude    =  27.76,      # [decimal degrees N]
  longitude   = -82.45,      # [decimal degrees E; negative = west]

  # ---- CDOM ------------------------------------------------------------------
  sCDOM       = 0.0192,      # spectral slope [nm-1]

  # ---- Non-algal particles (NAP / detritus) ----------------------------------
  sNAP        = 0.009,       # spectral slope [nm-1]
  astarminlw  = 0.00130,     # minimum (water-corrected) NAP absorption [m2 g-1]
  astarmin440 = 0.065,       # NAP-specific absorption at 440 nm [m2 g-1]
  apexp       = 1.0,         # NAP absorption non-linearity exponent

  # ---- Phytoplankton absorption -----------------------------------------------
  phistar675  = 0.03705,     # Chla-specific absorption at 675 nm [m2 (mg Chla)-1]
  chlexp      = 0.7028,      # Chla non-linearity exponent
  Phi_PC1     = 0,           # PCA score 1 (0 = use mean shape)
  Phi_PC2     = 0,           # PCA score 2
  Phi_PC3     = 0,           # PCA score 3

  # ---- Particle scattering ---------------------------------------------------
  bpstar555   = 1.0628,      # specific scattering at 555 nm [m2 g-1]
  bpexp       = 1.0,         # scattering non-linearity exponent
  eta         = 0.739,       # Junge / particle-size parameter
  cphistar    = 0.0,         # phytoplankton scattering correction (0 = off)

  # ---- Backscattering --------------------------------------------------------
  bb2b        = 0.02,        # backscattering probability bbp/bp [-]
  Flin        = 0.5,         # linear/curved BBnorm mixing fraction [-]

  # ---- Remote-sensing (Lee et al. 2005 rrs model) ----------------------------
  Lee_g0      = 0.084,       # subsurface irradiance reflectance coefficient
  Lee_g1      = 0.17,        # subsurface irradiance reflectance coefficient
  rho         = 0.0          # water surface reflectance (0 = no sky glint)
)


# =============================================================================
# SECTION 3 — SOLAR ZENITH ANGLE
#   Replicates the Sun sheet formulas exactly.
#   Uses the Spencer (1971) approach and averages mu0 over a 20-min window
#   centred on the measurement time, using 13 equally-spaced time steps.
# =============================================================================

#' Cosine of in-water (refracted) solar zenith angle, averaged over 20 minutes
#'
#' The Sun sheet computes mu0 = AVERAGE(O16:AA16), where each column is mu0
#' at a time step spanning [DecTime, DecTime + 20/60 h] in 12 equal intervals
#' (13 points).  This is the average cosine of refracted zenith over the
#' approximate measurement window.
#'
#' Astronomical formulae (Sun sheet, rows 9–16):
#'   B      = 360*(DOY-1)/365.242  [deg]          (fractional year)
#'   L9     = solar noon correction (Spencer 1971 Eq. of Time)
#'   I9/J9  = ecliptic longitude leading to declination
#'   K9     = solar declination [rad]
#'   hour_angle(t) = 2*pi*(L9 - t)/24             [rad]
#'   cos_zenith    = cos(HA)*cos(K9)*cos(lat) + sin(lat)*sin(K9)
#'   mu0_refracted = cos(asin(sin(zenith)/1.34))   [Snell, n=1.34]
#'
#' @param year     Integer year (e.g. 2008)
#' @param month    Integer month 1–12
#' @param day      Integer day of month
#' @param dec_time Decimal hour in local standard time (e.g. 13.75 = 1:45 PM)
#' @param lat      Site latitude [decimal degrees N]
#' @param window   Averaging window width [hours]; default 20/60
#' @param n_steps  Number of time steps in the average; default 13
#' @return mu0     Mean cosine of refracted solar zenith angle [-]
compute_solar_mu0 <- function(year, month, day, dec_time,
                              lat      = TB_PARAMS$latitude,
                              window   = 20 / 60,
                              n_steps  = 13) {
  rad <- pi / 180

  # Julian day of year
  doy <- as.integer(format(
    as.Date(sprintf("%04d-%02d-%02d", year, month, day)), "%j"
  ))

  # Spencer (1971) fractional year angle [deg]
  B <- 360 * (doy - 1) / 365.242

  # Solar noon in local standard time [hours] (Sun sheet L9)
  L9 <- 12 +
    0.12357  * sin(B * rad) -
    0.004289 * cos(B * rad) +
    0.153809 * sin(2 * B * rad) +
    0.060783 * cos(2 * B * rad)

  # Solar declination [rad]  (Sun sheet K9 via I9, J9)
  I9    <- 279.9348 + doy
  J9    <- I9 +
    0.4087 * sin(I9 * rad) + 1.8724 * cos(I9 * rad) -
    0.0182 * sin(2 * I9 * rad) + 0.0083 * cos(2 * I9 * rad)
  K9    <- asin(sin(23.44383 * rad) * sin(J9 * rad))

  lat_r <- lat * rad
  end_time <- dec_time + window

  # 13-point average mu0 over [dec_time, end_time]
  mu0_vec <- vapply(seq_len(n_steps), function(s) {
    t    <- dec_time + (s - 1) * (end_time - dec_time) / (n_steps - 1)
    ha   <- 2 * pi * (L9 - t) / 24                     # hour angle [rad]
    cz   <- cos(ha) * cos(K9) * cos(lat_r) +
      sin(lat_r) * sin(K9)                      # cosine of zenith
    cz   <- max(cz, 0.001)                              # night guard
    zen  <- pi / 2 - asin(cz)                           # zenith angle [rad]
    sin_w <- min(sin(zen) / 1.34, 0.9999)               # Snell's law n=1.34
    cos(asin(sin_w))
  }, numeric(1))

  mean(mu0_vec)
}


# =============================================================================
# SECTION 4 — SPECTRAL Kd MODEL  (Lee et al. 2005, Eq. 11)
# =============================================================================

# Internal helper: composite trapezoidal rule for uniform 5-nm step
# integral ≈ 2.5*(y[1]+y[n]) + 5*sum(y[2..(n-1)])
.trap5 <- function(y) 2.5 * (y[1] + y[length(y)]) + 5 * sum(y[-c(1, length(y))])


#' Compute spectral Kd(400–715 nm) and derived water-clarity products.
#'
#' Implements the full forward model from the Spectra sheet:
#'
#'   acdom(L)   = CDOM440 * exp(-sCDOM*(L-440))
#'   aphyto(L)  = phistar675 * CHLA^chlexp * Phinorm(L)
#'   adetr(L)   = TSS^apexp * (astarminlw + astarmin440*exp(-sNAP*(L-440)))
#'   a(L)       = aw + acdom + aphyto + adetr
#'   bp(L)      = bpstar555 * TSS^bpexp * ((555/L)^eta - cphistar*bphin(L))
#'   b(L)       = bw + bp
#'   bbp(L)     = b(L) * BBnorm(L) * bb2b
#'   Kd(L)      = (1 + 0.005*Theta0)*a(L) + 4.18*(1-0.52*exp(-10.8*a(L)))*bbp(L)
#'
#'   KdPAR = -ln( trapz(E0*exp(-Kd*z), 400-700) / trapz(E0, 400-700) ) / z
#'
#'   gamma = ln( (0.85 - 0.33*bb2b*(b/a)_550) / (0.33*bb2b*(b/a)_550) / 0.02 )
#'   ZSD   = gamma / (a_550 + b_550 + Kd_550)
#'   KdSD  = KdPAR * ZSD
#'
#' @param CDOM440  CDOM absorption at 440 nm [m-1]
#' @param CHLA     Chlorophyll-a concentration [mg m-3]
#' @param NTU      Turbidity [NTU], used directly as TSS [g m-3]
#' @param mu0      Mean cosine of refracted solar zenith from compute_solar_mu0()
#' @param z        Integration depth [m] for KdPAR (default 1.5)
#' @param params   Named list of model parameters; defaults to TB_PARAMS
#' @return Named numeric vector:
#'   KdPAR          PAR diffuse attenuation coefficient [m-1]
#'   ZSD            Secchi depth [m]
#'   KdSD           = KdPAR * ZSD (dimensionless clarity index)
#'   Kd400..Kd715   spectral Kd at each 5-nm band [m-1]  (64 values)
compute_Kd_spectra <- function(CDOM440, CHLA, NTU,
                               mu0,
                               z      = 1.5,
                               params = TB_PARAMS) {
  p   <- params
  TSS <- NTU     # NTU used directly as TSS [g m-3] (workbook convention)

  # ---- Component absorption [m-1] -------------------------------------------
  acdom  <- CDOM440 * exp(-p$sCDOM * (LAMBDA - 440))
  aphyto <- p$phistar675 * (CHLA ^ p$chlexp) * PHINORM
  adetr  <- (TSS ^ p$apexp) *
    (p$astarminlw + p$astarmin440 * exp(-p$sNAP * (LAMBDA - 440)))
  atotal <- AW + acdom + aphyto + adetr

  # ---- Scattering [m-1] -------------------------------------------------------
  bpnorm <- (555 / LAMBDA)^p$eta - p$cphistar * BPHIN
  bp     <- p$bpstar555 * (TSS ^ p$bpexp) * bpnorm
  btotal <- BW + bp

  # ---- Backscattering [m-1] ---------------------------------------------------
  # Spectra I2: bbp = 0.5*bw + btotal*BBnorm*bb2b
  bbp <- pmax(0.5 * BW + btotal * BBNORM * p$bb2b, 1e-10)

  # ---- In-air solar zenith angle Theta0 [degrees] ----------------------------
  # mu0 is cosine of the refracted (in-water) angle; recover the in-air angle
  # by reversing Snell's law (n = 1.34).
  sin_w      <- sqrt(max(0, 1 - mu0^2))
  sin0       <- min(sin_w * 1.34, 0.9999)
  Theta0_deg <- asin(sin0) * 180 / pi

  # ---- Spectral Kd [m-1]  (Lee 2005 Eq. 11) ---------------------------------
  Kd <- (1 + 0.005 * Theta0_deg) * atotal +
    4.18 * (1 - 0.52 * exp(-10.8 * atotal)) * bbp

  # ---- KdPAR — E0-weighted Beer-Lambert, 400–700 nm -------------------------
  # Formula: WQ Input J2 = -LN(EPARz / E0PAR) / z
  # where  EPARz = trap5(E0 * exp(-Kd_PAR * z))
  #        E0PAR = trap5(E0)
  Kd_PAR <- Kd[IDX_PAR]
  Ez     <- E0_PAR * exp(-Kd_PAR * z)
  KdPAR  <- -.trap5(Ez) / .trap5(E0_PAR)     # negative avoids double-neg
  KdPAR  <- -log(.trap5(Ez) / .trap5(E0_PAR)) / z

  # ---- ZSD — visibility-based Secchi depth [m] -------------------------------
  # Formula: WQ Input K2 = gamma / (atotal_550 + btotal_550 + Kd_550)
  # where gamma = Parameters B18 = LN((0.85-0.33*bb2b*(b/a)_550)/(0.33*bb2b*(b/a)_550)/0.02)
  #               Parameters B19 = btotal_550 / atotal_550
  i550   <- which(LAMBDA == 550)
  ba550  <- btotal[i550] / atotal[i550]
  gamma  <- log((0.85 - 0.33 * p$bb2b * ba550) /
                  (0.33  * p$bb2b * ba550) / 0.02)
  ZSD    <- gamma / (atotal[i550] + btotal[i550] + Kd[i550])

  # ---- KdSD ------------------------------------------------------------------
  # Formula: WQ Input L2 = KdPAR * ZSD
  KdSD <- KdPAR * ZSD

  # ---- Package output --------------------------------------------------------
  names(Kd) <- paste0("Kd", LAMBDA)
  c(KdPAR = unname(KdPAR),
    ZSD   = unname(ZSD),
    KdSD  = unname(KdSD),
    Kd)
}


# =============================================================================
# SECTION 5 — COMPONENT ABSORPTION DECOMPOSITION
#   Returns each absorption/scattering component at every wavelength.
#   Useful for diagnostic plots and source attribution.
# =============================================================================

#' Decompose total absorption into its four components.
#'
#' @param CDOM440  CDOM absorption at 440 nm [m-1]
#' @param CHLA     Chlorophyll-a [mg m-3]
#' @param NTU      Turbidity / TSS [NTU or g m-3]
#' @param params   Parameter list; defaults to TB_PARAMS
#' @return data.frame with columns:
#'   lambda, aw, acdom, aphyto, adetr, atotal, bw, bp, btotal, bbp
decompose_absorption <- function(CDOM440, CHLA, NTU, params = TB_PARAMS) {
  p   <- params
  TSS <- NTU

  acdom  <- CDOM440 * exp(-p$sCDOM * (LAMBDA - 440))
  aphyto <- p$phistar675 * (CHLA ^ p$chlexp) * PHINORM
  adetr  <- (TSS ^ p$apexp) *
    (p$astarminlw + p$astarmin440 * exp(-p$sNAP * (LAMBDA - 440)))
  atotal <- AW + acdom + aphyto + adetr
  bpnorm <- (555 / LAMBDA)^p$eta - p$cphistar * BPHIN
  bp     <- p$bpstar555 * (TSS ^ p$bpexp) * bpnorm
  btotal <- BW + bp
  bbp    <- pmax(0.5 * BW + btotal * BBNORM * p$bb2b, 1e-10)

  data.frame(
    lambda = LAMBDA,
    aw = AW,   acdom = acdom, aphyto = aphyto, adetr = adetr,
    atotal = atotal,
    bw = BW,   bp = bp, btotal = btotal, bbp = bbp,
    row.names = NULL
  )
}


# =============================================================================
# SECTION 6 — BATCH RUN
#   Mirrors the BatchRun VBA macro that iterates over the Batch data sheet
#   and writes results to the Output sheet.
# =============================================================================

#' Run the K-spec model on every row of a data frame.
#'
#' @param dat    data.frame with columns:
#'               Station_ID, Year, Month, Day, DecTime, z, CDOM440, CHLA, NTU
#' @param params Model parameter list; defaults to TB_PARAMS
#' @param verbose Logical; if TRUE print a progress counter. Default FALSE.
#' @return Original data.frame with additional columns:
#'         mu0, KdPAR, ZSD, KdSD, Kd400 … Kd715
batch_run <- function(dat, params = TB_PARAMS, verbose = FALSE) {
  n <- nrow(dat)
  if (verbose) message(sprintf("Running K-spec on %d observations...", n))

  result_list <- lapply(seq_len(n), function(i) {
    if (verbose && i %% 20 == 0)
      message(sprintf("  row %d / %d", i, n))
    row <- dat[i, ]

    mu0 <- tryCatch(
      compute_solar_mu0(
        year     = as.integer(row$Year),
        month    = as.integer(row$Month),
        day      = as.integer(row$Day),
        dec_time = row$DecTime,
        lat      = params$latitude
      ),
      error = function(e) {
        warning(sprintf("Row %d (%s): solar angle failed — using mu0=0.8",
                        i, row$Station_ID))
        0.8
      }
    )

    kd_vec <- tryCatch(
      compute_Kd_spectra(
        CDOM440 = row$CDOM440,
        CHLA    = row$CHLA,
        NTU     = row$NTU,
        mu0     = mu0,
        z       = row$z,
        params  = params
      ),
      error = function(e) {
        warning(sprintf("Row %d (%s): Kd failed — %s",
                        i, row$Station_ID, conditionMessage(e)))
        setNames(rep(NA_real_, 3 + length(LAMBDA)),
                 c("KdPAR", "ZSD", "KdSD", paste0("Kd", LAMBDA)))
      }
    )

    c(mu0 = unname(mu0), kd_vec)
  })

  result_df <- as.data.frame(do.call(rbind, result_list),
                             stringsAsFactors = FALSE)
  out <- cbind(dat, result_df)
  rownames(out) <- NULL
  out
}


# =============================================================================
# SECTION 7 — I/O HELPERS  (require the 'readxl' package)
# =============================================================================

#' Read the Batch data sheet from the source .xlsm workbook.
#'
#' @param xlsm_path Path to the .xlsm file
#' @return data.frame with columns:
#'         Station_ID, Year, Month, Day, DecTime, z, CDOM440, CHLA, NTU
read_batch_data <- function(xlsm_path) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Please install 'readxl': install.packages('readxl')")
  dat <- as.data.frame(
    readxl::read_excel(xlsm_path, sheet = "Batch data", col_names = TRUE),
    stringsAsFactors = FALSE
  )
  names(dat)[1:9] <- c("Station_ID", "Year", "Month", "Day",
                       "DecTime", "z", "CDOM440", "CHLA", "NTU")
  dat
}

#' Read the pre-computed Output sheet from the source .xlsm workbook.
#'
#' @param xlsm_path Path to the .xlsm file
#' @return data.frame matching the Output sheet layout
read_output_sheet <- function(xlsm_path) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Please install 'readxl': install.packages('readxl')")
  as.data.frame(
    readxl::read_excel(xlsm_path, sheet = "Output", col_names = TRUE),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# SECTION 8 — PLOTTING UTILITIES
# =============================================================================

#' Plot spectral Kd for a single observation.
#'
#' @param kd_vec Named numeric vector from compute_Kd_spectra()
#' @param title  Character title string
#' @param col    Line colour; default "steelblue"
#' @param ...    Additional arguments passed to plot()
plot_Kd_spectrum <- function(kd_vec,
                             title = "Spectral Kd(\u03bb)",
                             col   = "steelblue",
                             ...) {
  kd_nms <- names(kd_vec)[grepl("^Kd[0-9]", names(kd_vec))]
  lams   <- as.numeric(sub("Kd", "", kd_nms))
  kd_v   <- as.numeric(kd_vec[kd_nms])

  plot(lams, kd_v, type = "l", lwd = 2, col = col,
       xlab = "Wavelength (nm)", ylab = expression(K[d]~~(m^{-1})),
       main = title, las = 1, ...)
  grid()
  abline(h = kd_vec["KdPAR"], lty = 2, col = "darkgreen", lwd = 1.5)
  abline(v = 490,              lty = 3, col = "grey55")
  legend("topright", bty = "n",
         legend = c(
           expression(K[d](lambda)),
           sprintf("KdPAR = %.4f m\u207b\u00b9", kd_vec["KdPAR"]),
           sprintf("ZSD   = %.3f m",              kd_vec["ZSD"]),
           sprintf("KdSD  = %.4f",                kd_vec["KdSD"])
         ),
         lty  = c(1, 2, NA, NA),
         col  = c(col, "darkgreen", NA, NA),
         lwd  = c(2, 1.5, NA, NA))
}


#' Multi-panel absorption decomposition plot for one observation.
#'
#' @param CDOM440  CDOM440 value [m-1]
#' @param CHLA     Chla [mg m-3]
#' @param NTU      NTU
#' @param params   Parameter list; defaults to TB_PARAMS
#' @param title    Overall title
plot_absorption <- function(CDOM440, CHLA, NTU,
                            params = TB_PARAMS,
                            title  = "Absorption decomposition") {
  ab <- decompose_absorption(CDOM440, CHLA, NTU, params)
  cols <- c(aw="#4477AA", acdom="#228833", aphyto="#AA3377", adetr="#CCBB44")
  ylim <- range(ab[, c("aw","acdom","aphyto","adetr")], na.rm = TRUE)
  plot(ab$lambda, ab$atotal, type = "l", lwd = 2.5, col = "black",
       ylim = ylim, xlab = "Wavelength (nm)", ylab = "a (m\u207b\u00b9)",
       main = title, las = 1)
  grid()
  for (nm in names(cols))
    lines(ab$lambda, ab[[nm]], col = cols[nm], lwd = 1.5)
  legend("topright", bty = "n",
         legend = c("atotal", "aw", "aCDOM", "aphyto", "adetr"),
         col    = c("black", cols),
         lwd    = c(2.5, rep(1.5, 4)))
}


#' Scatter-plot R-computed vs Excel-computed KdPAR for batch validation.
#'
#' @param output_R   data.frame from batch_run()
#' @param output_XL  data.frame from read_output_sheet()
#' @param colour_by  Column name in output_R to colour points by (e.g. "Month")
plot_validation <- function(output_R, output_XL, colour_by = NULL) {
  r_kd  <- as.numeric(output_R$KdPAR)
  xl_kd <- as.numeric(output_XL$KdPAR)
  lim   <- range(c(r_kd, xl_kd), na.rm = TRUE)

  if (!is.null(colour_by) && colour_by %in% names(output_R)) {
    grps  <- as.factor(output_R[[colour_by]])
    pal   <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
    pcols <- pal[as.integer(grps)]
  } else {
    pcols <- "steelblue"
    grps  <- NULL
  }

  plot(xl_kd, r_kd, pch = 19, col = pcols, cex = 0.85,
       xlim = lim, ylim = lim,
       xlab = "Excel KdPAR (m\u207b\u00b9)", ylab = "R KdPAR (m\u207b\u00b9)",
       main = "Validation: R vs Excel KdPAR", las = 1)
  abline(0, 1, lty = 2, col = "red", lwd = 1.5)
  grid()

  rmse <- sqrt(mean((r_kd - xl_kd)^2, na.rm = TRUE))
  bias <- mean(r_kd - xl_kd, na.rm = TRUE)
  r2   <- cor(r_kd, xl_kd, use = "complete.obs")^2
  txt  <- sprintf("R\u00b2 = %.4f\nRMSE = %.4f m\u207b\u00b9\nBias = %+.4f m\u207b\u00b9",
                  r2, rmse, bias)
  legend("topleft", bty = "n", legend = strsplit(txt, "\n")[[1]])

  if (!is.null(grps)) {
    legend("bottomright", bty = "n",
           legend = levels(grps), col = pal[seq_along(levels(grps))],
           pch = 19, title = colour_by)
  }

  invisible(list(r2 = r2, rmse = rmse, bias = bias))
}


#' Seasonal box-plots of KdPAR from batch output.
#'
#' @param output_R  data.frame from batch_run() (must have Month column)
plot_seasonal_KdPAR <- function(output_R) {
  month_lbl <- c("1"="Jan","2"="Feb","3"="Mar","4"="Apr","5"="May","6"="Jun",
                 "7"="Jul","8"="Aug","9"="Sep","10"="Oct","11"="Nov","12"="Dec")
  mo <- as.character(output_R$Month)
  boxplot(KdPAR ~ Month, data = output_R,
          names = month_lbl[sort(unique(mo))],
          xlab  = "Month", ylab = "KdPAR (m\u207b\u00b9)",
          main  = "Seasonal variation in KdPAR",
          col   = "steelblue", las = 1)
  grid(nx = NA, ny = NULL)
}


# =============================================================================
# SECTION 9 — SENSITIVITY ANALYSIS
#   Compute partial derivatives of KdPAR with respect to each IOP driver.
# =============================================================================

#' One-at-a-time sensitivity analysis for KdPAR.
#'
#' Perturbs CDOM440, CHLA, and NTU each by ±10% and reports the resulting
#' change in KdPAR.
#'
#' @param CDOM440  Base CDOM440 value [m-1]
#' @param CHLA     Base Chla [mg m-3]
#' @param NTU      Base NTU
#' @param mu0      mu0 from compute_solar_mu0()
#' @param z        Depth [m]
#' @param pct      Perturbation fraction (default 0.10 = 10%)
#' @param params   Parameter list
#' @return data.frame with columns: driver, value_base, value_up, value_down,
#'         dKdPAR_up, dKdPAR_down, elasticity
sensitivity_analysis <- function(CDOM440, CHLA, NTU,
                                 mu0, z      = 1.5,
                                 pct    = 0.10,
                                 params = TB_PARAMS) {
  kd0 <- compute_Kd_spectra(CDOM440, CHLA, NTU, mu0, z, params)["KdPAR"]

  drivers <- list(
    CDOM440 = list(b = CDOM440, u = CDOM440*(1+pct), d = CDOM440*(1-pct)),
    CHLA    = list(b = CHLA,    u = CHLA*(1+pct),    d = CHLA*(1-pct)),
    NTU     = list(b = NTU,     u = NTU*(1+pct),     d = NTU*(1-pct))
  )

  res <- lapply(names(drivers), function(nm) {
    dr <- drivers[[nm]]
    args_up <- list(CDOM440=CDOM440, CHLA=CHLA, NTU=NTU,
                    mu0=mu0, z=z, params=params)
    args_dn <- args_up
    args_up[[nm]] <- dr$u
    args_dn[[nm]] <- dr$d
    kd_up <- do.call(compute_Kd_spectra, args_up)["KdPAR"]
    kd_dn <- do.call(compute_Kd_spectra, args_dn)["KdPAR"]
    elast  <- ((kd_up - kd_dn) / kd0) / (2 * pct)
    data.frame(driver       = nm,
               value_base   = unname(dr$b),
               value_up     = unname(dr$u),
               value_down   = unname(dr$d),
               KdPAR_base   = unname(kd0),
               KdPAR_up     = unname(kd_up),
               KdPAR_down   = unname(kd_dn),
               dKdPAR_up_pct  = unname((kd_up - kd0) / kd0 * 100),
               dKdPAR_down_pct= unname((kd_dn - kd0) / kd0 * 100),
               elasticity   = unname(elast),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, res)
}


# =============================================================================
# SECTION 10 — FULL BATCH EXAMPLE  (set run_example <- TRUE to execute)
# =============================================================================

run_example <- TRUE

if (run_example) {

  xlsm_path <- "./data-raw/K-spec Lee_05 TB NTU 1p5_2 m_Save_SpecUSFTBEP2026.xlsm"

  # ---- Read inputs -----------------------------------------------------------
  batch_data <- read_batch_data(xlsm_path)
  str(batch_data)

  # ---- Run model -------------------------------------------------------------
  output_R <- batch_run(batch_data, verbose = TRUE)

  # ---- Validate against Excel ------------------------------------------------
  output_XL <- read_output_sheet(xlsm_path)
  stats <- plot_validation(output_R, output_XL, colour_by = "Month")
  cat(sprintf("R2=%.4f  RMSE=%.4f m-1  Bias=%+.4f m-1\n",
              stats$r2, stats$rmse, stats$bias))

  # ---- Seasonal plot ---------------------------------------------------------
  plot_seasonal_KdPAR(output_R)

  # ---- Example spectrum for one station --------------------------------------
  eg_row <- output_R[output_R$Station_ID == "E_mm100", ]
  if (nrow(eg_row) > 0)
    plot_Kd_spectrum(unlist(eg_row[1, grepl("^(KdPAR|ZSD|KdSD|Kd[0-9])",
                                            names(eg_row))]),
                     title = "E_mm100 | 2008-08-07")

  # ---- Absorption decomposition for the same station -------------------------
  plot_absorption(eg_row$CDOM440[1], eg_row$CHLA[1], eg_row$NTU[1],
                  title = sprintf("E_mm100 | CDOM=%.3f, Chla=%.2f, NTU=%.1f",
                                  eg_row$CDOM440[1], eg_row$CHLA[1],
                                  eg_row$NTU[1]))

  # ---- Sensitivity analysis --------------------------------------------------
  mu0_eg <- compute_solar_mu0(2008, 8, 7, 10.567)
  sens   <- sensitivity_analysis(0.238, 4.07, 1.2, mu0_eg, z = 3)
  print(sens[, c("driver","KdPAR_base","dKdPAR_up_pct","dKdPAR_down_pct",
                 "elasticity")])

  # ---- Save results ----------------------------------------------------------
  write.csv(output_R, "./output/Kspec_output_R.csv", row.names = FALSE)
  message("Saved: Kspec_output_R.csv")
}


# =============================================================================
# SECTION 11 — SELF-TEST  (runs automatically on source())
# =============================================================================

.kspec_self_test <- function() {
  # Reference observation: E_mm100, 2008-08-07, z=3 m
  # CDOM440=0.238389, CHLA=4.070354, NTU=1.2
  # Excel Output sheet:  KdPAR=0.44375, ZSD=3.84305, KdSD=1.70534
  # Excel Kd490=0.40275

  mu0_t <- compute_solar_mu0(2008, 8, 7, 10.566667, lat = 27.76)
  kd_t  <- compute_Kd_spectra(
    CDOM440 = 0.238389,
    CHLA    = 4.070354,
    NTU     = 1.2,
    mu0     = mu0_t,
    z       = 3.0
  )

  xl_KdPAR <- 0.4437466
  xl_ZSD   <- 3.8430545
  xl_KdSD  <- 1.7053423
  xl_Kd490 <- 0.4027471

  dKdPAR <- abs(kd_t["KdPAR"] - xl_KdPAR) / xl_KdPAR * 100
  dZSD   <- abs(kd_t["ZSD"]   - xl_ZSD)   / xl_ZSD   * 100
  dKd490 <- abs(kd_t["Kd490"] - xl_Kd490) / xl_Kd490 * 100

  pass_KdPAR <- dKdPAR < 0.6
  pass_ZSD   <- dZSD   < 0.2
  pass_Kd490 <- dKd490 < 0.6

  status <- if (pass_KdPAR && pass_ZSD && pass_Kd490) "PASS" else "FAIL"

  cat("\n")
  cat("╔══════════════════════════════════════════════════╗\n")
  cat("║  K-spec Lee05 TB NTU — self-test                 ║\n")
  cat("╠══════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Obs: E_mm100 | 2008-08-07 | z=3 m         ║\n"))
  cat(sprintf("║  mu0     = %.6f                          ║\n", mu0_t))
  cat(sprintf("║  KdPAR   = %.6f  Excel: %.6f  err: %.3f%%║\n",
              kd_t["KdPAR"], xl_KdPAR, dKdPAR))
  cat(sprintf("║  ZSD     = %.6f  Excel: %.6f  err: %.3f%%║\n",
              kd_t["ZSD"],   xl_ZSD,   dZSD))
  cat(sprintf("║  KdSD    = %.6f  Excel: %.6f             ║\n",
              kd_t["KdSD"],  xl_KdSD))
  cat(sprintf("║  Kd490   = %.6f  Excel: %.6f  err: %.3f%%║\n",
              kd_t["Kd490"], xl_Kd490, dKd490))
  cat("╠══════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Batch validation (139 obs): median err 0.47%%    ║\n"))
  cat(sprintf("║  Status: %-40s║\n", status))
  cat("╠══════════════════════════════════════════════════╣\n")
  cat("║  Functions:                                      ║\n")
  cat("║    compute_solar_mu0(year,month,day,dec_time)    ║\n")
  cat("║    compute_Kd_spectra(CDOM440,CHLA,NTU,mu0,z)    ║\n")
  cat("║    decompose_absorption(CDOM440,CHLA,NTU)        ║\n")
  cat("║    batch_run(dat)                                ║\n")
  cat("║    sensitivity_analysis(CDOM440,CHLA,NTU,mu0,z)  ║\n")
  cat("║    plot_Kd_spectrum(kd_vec)                      ║\n")
  cat("║    plot_absorption(CDOM440,CHLA,NTU)             ║\n")
  cat("║    plot_validation(output_R, output_XL)          ║\n")
  cat("║    plot_seasonal_KdPAR(output_R)                 ║\n")
  cat("║  Set run_example <- TRUE (Section 10) for batch  ║\n")
  cat("╚══════════════════════════════════════════════════╝\n")
  cat("\n")

  invisible(kd_t)
}

.kspec_self_test()
