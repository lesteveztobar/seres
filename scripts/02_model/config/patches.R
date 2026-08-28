
# ecmwfr >= 2.0 added a `job_name` argument that the CDS API endpoint rejects.
# This wraps the original and silently drops job_name before forwarding the call.
if (!exists("original_wf_request")) {
  original_wf_request <- ecmwfr::wf_request
  assignInNamespace("wf_request",
                    function(request, user = "ecmwfr", transfer = TRUE,
                             path = tempdir(), time_out = 7200, retry = 240,
                             job_name, verbose = TRUE) {
                      original_wf_request(
                        request = request, user = user, transfer = transfer,
                        path = path, time_out = time_out, retry = retry,
                        verbose = verbose
                      )
                    },
                    ns = "ecmwfr"
  )
}

# microclimdata::reflectance_calc uses `wgt` in the iterative update but the
# variable is named `bwgt` — causing it to fail. This replaces the function
# with the corrected version.
assignInNamespace("reflectance_calc",
                  function(alb, lai, x, plotprogress = TRUE, maxiter = 50, tol = 0.001, bwgt = 0.5) {
                    e1 <- terra::intersect(terra::ext(lai), terra::ext(alb))
                    e  <- terra::intersect(e1, terra::ext(x))
                    lai <- terra::crop(lai, e)
                    alb <- terra::crop(alb, e)
                    x   <- terra::crop(x, e)
                    all_same <- terra::compareGeom(lai, alb, x)
                    if (all_same) {
                      tst  <- exp(-mean(as.vector(lai), na.rm = TRUE))
                      # guard: if lai is all-NA the mean is NaN, making tst NaN and breaking if()
                      if (is.nan(tst) || is.na(tst)) tst <- 0
                      lref <- (x * 0 + 0.5) * (1 - bwgt) + bwgt * alb  # fix: wgt -> bwgt
                      gref <- x * 0 + 0.15
                      mxdif <- tol * 10
                      paim  <- as.matrix(lai, wide = TRUE)
                      xm    <- as.matrix(x,   wide = TRUE)
                      albm  <- as.matrix(alb, wide = TRUE)
                      itr   <- 1
                      while (mxdif > tol) {
                        if (tst < 0.5) {
                          lref2 <- microclimdata:::.rast(microclimdata:::find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
                          lref2 <- microclimdata:::.fillna(lref2, x, zerotoNA = FALSE)
                          gref2 <- microclimdata:::.rast(microclimdata:::find_gref(as.matrix(lref2, wide = TRUE), paim, xm, albm), x)
                          gref2 <- microclimdata:::.fillna(gref2, x, zerotoNA = FALSE)
                        } else {
                          gref2 <- microclimdata:::.rast(microclimdata:::find_gref(as.matrix(lref, wide = TRUE), paim, xm, albm), x)
                          gref2 <- microclimdata:::.fillna(gref2, x, zerotoNA = FALSE)
                          lref2 <- microclimdata:::.rast(microclimdata:::find_lref(paim, as.matrix(gref, wide = TRUE), xm, albm), x)
                          lref2 <- microclimdata:::.fillna(lref2, x, zerotoNA = FALSE)
                        }
                        gref  <- bwgt * gref + (1 - bwgt) * gref2
                        lref  <- bwgt * lref + (1 - bwgt) * lref2
                        mxdif1 <- mean(abs(as.vector(gref) - as.vector(gref2)), na.rm = TRUE)
                        mxdif2 <- mean(abs(as.vector(lref) - as.vector(lref2)), na.rm = TRUE)
                        mxdif  <- max(mxdif1, mxdif2)
                        # if all inputs are NA (e.g. from a degenerate LAI raster), mxdif is NaN
                        if (is.nan(mxdif) || is.na(mxdif))
                          stop("reflectance_calc: all-NA result — check that lai_agg and alb_agg cover the site extent")
                        itr    <- itr + 1
                        if (itr > maxiter) mxdif <- 0
                      }
                    } else {
                      stop("Geometries of input rasters do not match")
                    }
                    return(list(gref = gref, lref = lref))
                  },
                  ns = "microclimdata"
)

# microclimdata::lai_mosaic (reso=500) calls rast(fi)[[2]] and rast(fi)[[4]] to
# open MOD15A2H HDF4 files. Newer terra/GDAL requires the subdataset to be named
# explicitly. Subdataset 2 = Lai_500m, subdataset 4 = FparExtra_QC (QC mask).
# The reso=10 branch uses .tif files and is unchanged.
assignInNamespace("lai_mosaic",
                  function(r, pathin, reso = 10, msk = TRUE) {
                    isres <- reso %in% c(10, 500)
                    if (isres == FALSE) stop("reso must be one of 10 or 500")
                    lst <- list.files(pathin)
                    if (reso == 10) {
                      # reso=10 uses HRVPP .tif files — no HDF4 issue, keep original logic
                      ns        <- getNamespace("microclimdata")
                      .chcktif  <- get(".chcktif",  envir = ns)
                      .cleanlai <- get(".cleanlai", envir = ns)
                      wlai <- rep(0, length(lst))
                      for (i in 1:length(lst)) {
                        n <- nchar(lst[[i]]); ed <- n - 4; st <- n - 6
                        ss <- substr(lst[[i]], st, ed)
                        if (ss == "LAI") wlai[i] <- 1
                      }
                      s <- which(wlai == 1)
                      lfile <- lst[s]
                      qfile <- paste0(substr(lfile, 1, 40), "QFLAG2.tif")
                      ro <- list(); tile <- 0; n <- length(lfile)
                      pb <- utils::txtProgressBar(min = 0, max = n * 1.5, style = 3)
                      for (i in 1:n) {
                        utils::setTxtProgressBar(pb, i)
                        fi <- paste0(pathin, lfile[i])
                        ri <- terra::rast(fi)
                        e  <- terra::ext(ri); re <- terra::rast(e); terra::crs(re) <- terra::crs(ri)
                        if (terra::crs(r) != terra::crs(re)) re <- terra::project(re, terra::crs(r))
                        oip <- terra::intersect(terra::ext(r), terra::ext(re))
                        if (!is.null(oip)) {
                          xx <- suppressWarnings(.chcktif(fi))
                          if (!is.na(xx)) {
                            fiq <- paste0(pathin, qfile[i])
                            if (file.exists(fiq)) {
                              qu <- terra::rast(fiq)
                              ri <- .cleanlai(ri, qu)
                              if (terra::crs(ri) != terra::crs(r)) ri <- terra::project(ri, r)
                              ri <- terra::resample(ri, r, method = "near")
                              if (msk) ri <- terra::mask(ri, r)
                              v <- as.vector(r); v <- v[!is.na(v)]
                              if (length(v) > 0) { tile <- tile + 1; ro[[tile]] <- ri }
                            }
                          }
                        }
                      }
                      rma <- ro[[1]]; nn <- length(ro)
                      if (length(ro) > 1) {
                        for (i in 2:length(ro)) {
                          utils::setTxtProgressBar(pb, ((i * 0.5 * n) / nn) + n)
                          rma <- terra::mosaic(rma, ro[[i]])
                        }
                      }
                    } else {
                      # reso=500: MOD15A2H HDF4 — open subdatasets by name, not by index
                      .open_hdf <- function(fi, layer) {
                        # pathin already ends with "/" so use paste0 to avoid double-slash in path
                        sds <- paste0('HDF4_EOS:EOS_GRID:"', fi, '":MOD_Grid_MOD15A2H:', layer)
                        terra::rast(sds)
                      }
                      # Project r to sinusoidal once so we can work in metric space throughout.
                      # All tiles are in sinusoidal — mosaic there, then reproject the combined
                      # result to WGS84 at the end. Reprojecting each tile separately produces
                      # slightly different resolutions in WGS84, which breaks terra::mosaic().
                      r_sinu <- terra::project(terra::rast(terra::ext(r), crs = terra::crs(r)),
                                               "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")
                      # single shared 500m template in sinusoidal — all tiles resample to this
                      # so they share an identical grid before mosaicking
                      tmpl_sinu <- terra::rast(terra::ext(r_sinu), res = 500,
                                               crs = terra::crs(r_sinu))
                      ro <- list(); tile <- 0; n <- length(lst)
                      pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
                      for (i in 1:n) {
                        utils::setTxtProgressBar(pb, i)
                        fi   <- paste0(pathin, lst[i])
                        ri   <- .open_hdf(fi, "Lai_500m")
                        mskr <- .open_hdf(fi, "FparExtra_QC")
                        ri   <- terra::mask(ri, mskr)
                        oip  <- terra::intersect(terra::ext(r_sinu), terra::ext(ri))
                        if (!is.null(oip)) {
                          tile <- tile + 1
                          ro[[tile]] <- terra::resample(ri, tmpl_sinu)
                        }
                      }
                      if (length(ro) == 0) stop("No LAI tiles intersect the site extent — check downloaded tiles")
                      # mosaic in sinusoidal (all tiles share same 500m resolution)
                      rma <- ro[[1]]
                      if (length(ro) > 1) {
                        for (i in 2:length(ro)) rma <- terra::mosaic(rma, ro[[i]])
                      }
                      # reproject combined mosaic to WGS84 once
                      rma <- terra::project(rma, terra::crs(r))
                      if (msk) {
                        rmsk <- terra::resample(r, rma)
                        rma  <- terra::mask(rma, rmsk)
                      }
                    }
                    return(rma)
                  },
                  ns = "microclimdata"
)

# microclimdata::albedo_process calls rast(fi)[[30]] to open the HDF4 file and
# select the WSA shortwave band. Newer terra/GDAL requires the subdataset to be
# specified explicitly via the HDF4_EOS:EOS_GRID: prefix — opening the raw .hdf
# path fails with "not recognized as a supported file format". This replaces only
# the file-open line; all other logic is identical to the original.
assignInNamespace("albedo_process",
                  function(r, pathin) {
                    lst <- list.files(pathin)
                    if (length(lst) == 0) stop("No files to process!")
                    ns      <- getNamespace("microclimdata")
                    apply3D <- get("apply3D", envir = ns)
                    .rast   <- get(".rast",   envir = ns)
                    
                    .open_wsa <- function(fi) {
                      # pathin already ends with "/" so use paste0, not file.path, to avoid
                      # a double slash that breaks GDAL's HDF4_EOS: subdataset path parsing
                      sds <- paste0('HDF4_EOS:EOS_GRID:"', fi, '":MOD_Grid_BRDF:Albedo_WSA_shortwave')
                      terra::rast(sds)
                    }
                    
                    pb   <- utils::txtProgressBar(min = 0, max = length(lst) + 1, style = 3)
                    fi   <- paste0(pathin, lst[1])
                    modr <- .open_wsa(fi)
                    bbx  <- terra::rast(terra::ext(r)); terra::crs(bbx) <- terra::crs(r)
                    bbx  <- terra::project(bbx, terra::crs(modr))
                    e    <- terra::ext(bbx)
                    e$xmin <- e$xmin - 1000; e$xmax <- e$xmax + 1000
                    e$ymin <- e$ymin - 1000; e$ymax <- e$ymax + 1000
                    modr <- terra::extend(modr, e); modr <- terra::crop(modr, e)
                    utils::setTxtProgressBar(pb, 1)
                    
                    for (i in 2:length(lst)) {
                      fi <- paste0(pathin, lst[i])
                      mr <- .open_wsa(fi)
                      mr <- terra::extend(mr, e); mr <- terra::crop(mr, e)
                      modr <- c(modr, mr)
                      utils::setTxtProgressBar(pb, i)
                    }
                    
                    m      <- apply3D(as.array(modr))
                    albedo <- .rast(m, modr)
                    albedo <- terra::project(albedo, terra::crs(r))
                    albedo <- terra::crop(albedo, terra::ext(r))
                    utils::setTxtProgressBar(pb, length(lst) + 1)
                    return(albedo)
                  },
                  ns = "microclimdata"
)
# microclimf::checkinputs computes the minimum allowed pressure (mnp) using
# mnelev (lowest DTM point) instead of mxelev (highest point). At high-elevation
# sites like Maquipucuna, ERA5 surface pressure (~76 kPa at 2300m) falls below
# the incorrectly-computed floor (~87 kPa near sea-level). Fix: swap mnelev→mxelev
# for mnp so the floor reflects pressure at the highest terrain, not the lowest.
local({
  ns  <- getNamespace("microclimf")
  orig <- get("checkinputs", envir = ns)
  body_txt <- deparse(body(orig))
  # Fix the one wrong line: mnp uses mnelev but should use mxelev
  body_txt <- gsub(
    "mnp <- 87 \\* \\(\\(293 - 0\\.0065 \\* mnelev\\)/293\\)\\^5\\.26",
    "mnp <- 87 * ((293 - 0.0065 * mxelev)/293)^5.26",
    body_txt
  )
  fixed_body <- parse(text = paste(body_txt, collapse = "\n"))[[1]]
  fixed_fn   <- orig
  body(fixed_fn) <- fixed_body
  assignInNamespace("checkinputs", fixed_fn, ns = "microclimf")
})

# microclimf::runmicro takes a `method` argument ("R" | "Cpp") and forwards
# it to .runmicrosnow1()/.runmicrosnow2() (snow=TRUE path -- neither of
# those two actually has a method parameter, so nothing to fix there) and
# .runmicronosnow() (snow=FALSE, the path this project always uses) -- but
# the call to .runmicronosnow() never actually includes `method` in its
# argument list, so it silently falls back to .runmicronosnow()'s own
# method="Cpp" default regardless of what runmicro(method=...) was called
# with. Discovered 2026-08-10 while investigating why an earlier method="R"
# vs method="Cpp" timing/correctness comparison found them near-identical
# in both speed and output (max|diff|=0) -- turns out both runs silently
# took the exact same Cpp-branch code path the whole time. Doesn't change
# this project's own behavior (method="Cpp" is what run_microclimate_site.R
# already requests, and "Cpp" was always the effective default here either
# way), but the parameter should genuinely reach where it's supposed to.
local({
  ns   <- getNamespace("microclimf")
  orig <- get("runmicro", envir = ns)
  body_txt <- paste(deparse(body(orig)), collapse = "\n")
  fixed_txt <- gsub(
    "mout\\s*<-\\s*\\.runmicronosnow\\(micropoint,\\s*reqhgt,\\s*vegp,\\s*soilc,\\s*dtm,\\s*dtmc,\\s*altcorrect,\\s*runchecks,\\s*pai_a,\\s*tfact,\\s*out,\\s*slr,\\s*apr,\\s*hor,\\s*twi,\\s*wsa,\\s*svf\\)",
    "mout <- .runmicronosnow(micropoint, reqhgt, vegp, soilc, dtm, dtmc, altcorrect, runchecks, pai_a, tfact, out, slr, apr, hor, twi, wsa, svf, method)",
    body_txt, perl = TRUE
  )
  # Fail loudly rather than silently patching nothing if a future package
  # update reformats this call and the pattern stops matching.
  if (identical(fixed_txt, body_txt)) {
    stop("patches.R: runmicro() method-forwarding patch pattern no longer matches -- microclimf package source may have changed, update the patch.")
  }
  fixed_body <- parse(text = fixed_txt)[[1]]
  fixed_fn   <- orig
  body(fixed_fn) <- fixed_body
  assignInNamespace("runmicro", fixed_fn, ns = "microclimf")
})

# microclimdata::vegheight_download was written for Windows and appends
# "python.exe" to pathtopython. On macOS/Linux the binary is just "python".
assignInNamespace("vegheight_download",
                  function(r, GoogleDrivefolder, pathtopython, projectname = NA, silent = FALSE) {
                    reticulate::use_python(paste0(pathtopython, "python"), required = TRUE)
                    if (!is.na(projectname)) rgee::ee$Initialize(project = projectname)
                    e  <- terra::ext(r)
                    r2 <- terra::rast(e); terra::crs(r2) <- terra::crs(r)
                    r2 <- terra::project(r2, "EPSG:4326")
                    e  <- terra::ext(r2)
                    proj_string <- terra::crs(r, describe = TRUE)
                    epsg_code   <- paste0("EPSG:", proj_string$code)
                    aoi         <- rgee::ee$Geometry$Rectangle(c(e$xmin, e$ymin, e$xmax, e$ymax))
                    aoi_bounds  <- aoi$bounds()$getInfo()
                    aoi_coords  <- aoi_bounds$coordinates[[1]]
                    canopy_height <- rgee::ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1")
                    task <- rgee::ee$batch$Export$image$toDrive(
                      image          = canopy_height,
                      description    = "canopy_height_export",
                      folder         = GoogleDrivefolder,
                      fileNamePrefix = "canopy_height_2020",
                      region         = aoi_coords,
                      scale          = 10,
                      crs            = epsg_code
                    )
                    task$start()
                    if (!silent) microclimdata:::.monitor_task(task$id)
                  },
                  ns = "microclimdata"
)