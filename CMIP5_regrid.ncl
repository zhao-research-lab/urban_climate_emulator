; Regridding from GCM's original grid to CESM 1 degree grid
; Writing data to a NetCDF file using the ESMF regridding framework
;
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
    PLOT_RESULTS  = True
    WRITE_RESULTS = True
    INTERP_METHOD = "patch"           ;ESMF methods: bilinear, patch, conserve, neareststod

;----------------------------------------------------------------------
; Regridding section
;----------------------------------------------------------------------

;**** Define const variables
    RCP         = (/"RCP45", "RCP85"/)
    CMIP5List   = asciiread("[YOUR CMIP5 NAME LIST FILE]", -1, "string")  ; load CMIP5 model name list
    varList     = (/"ta", "uas", "vas", "huss", "ps", "pr", "rsds", "rlds"/)
    dataDir     = "/[YOUR CMIP5 DATA DIR]/"    
;**** Read data from input files containing source grid and destination grid
    do i=0,1                                  ;loop through RCPs
    	rcp         = RCP(i)
        print("-----Processing "+rcp+" ......")
        CESM_TSA_U_Filename  = "b40.rcp4_5.1deg.001.clm2.h0.TSA_U.200501-210012.nc"
        dfile       = addfile(dataDir + "RCP45" + "/CESM/" + CESM_TSA_U_Filename ,"r")   ; time * 192 * 288
        lat1d       = dfile->lat
        lon1d       = dfile->lon            
	do j=0,24                              ;loop through ESMs
    	    ESM      = CMIP5List(j)
            print("-----Processing ESM: "+ESM)
            srcDir      = dataDir + rcp + "/" + ESM + "/"
            prefix      = "_Amon_" 
    	    if (rcp .eq. "RCP45") then
    		suffix = "_rcp45_r1i1p1_200601-210012.nc"
    	    else
    		suffix = "_rcp85_r1i1p1_200601-210012.nc"
    	    end if
            do v=0,7                          ;loop through variables
    	    	var         = varList(v)     
    	    	srcFileName = srcDir + var + prefix + ESM + suffix
            	print("Reading CMIP5 ESM: " + ESM + " from Directory: " + srcDir)
            	sfile       = addfile(srcFileName, "r")     ;reading source data file
    	    	data        = sfile->$var$        ; variable dimension: (time, lat, lon)
    	    	printVarSummary(data)

;---Set up options for regridding
    	    	Opt                   = True
    	    	Opt@SrcFileName       = srcDir + ESM + "_SCRIP.nc"   ; Output intermediate files
            	Opt@DstFileName       = srcDir + "CESM_SCRIP.nc"      ; Output intermediate files
    	    	Opt@WgtFileName       = srcDir + ESM + "_2_CESM_patch_wgt.nc"   ; weights file name, default is "weights_file.nc"
    	    	Opt@SrcInputFileName  = srcFileName        ; optional, but good idea
    	    	Opt@DstInputFileName  = dataDir + rcp + "/CESM/" + CESM_TSA_U_Filename
    	    	;Opt@SrcTitle          = "Original Grid"    ; the ESM original grid
    	    	;Opt@SrcMask2D         = where(ismissing(data(0,:,:)),0,1)
    	    	Opt@ForceOverwrite    = True
    
    		Opt@DstGridLat        = lat1d
    		Opt@DstGridLon        = lon1d
    		Opt@InterpMethod      = INTERP_METHOD
    		;;Needed if you want to regrid without generating these files again.
    		Opt@SkipSrcGrid       = False
    		Opt@SkipDstGrid       = False
    		Opt@SkipWgtGen        = False
    		Opt@PrintTimings     = True
    		Opt@Debug            = True

;----------------------------------------------------------------------
; Regrid the data. "data" contains 1D coordinate arrays,
; so these are used automatically as the "source" grid.
;----------------------------------------------------------------------
    		if (v .eq. 0) then
		  data_regrid = ESMF_regrid(data,Opt)         ;this initial call to ESMF_regrid will generate weights_file for this pair of src-dst grids
    		  printVarSummary(data_regrid)
		else
		  data_regrid = ESMF_regrid_with_weights(data,srcDir + ESM + "_2_CESM_patch_wgt.nc",False)     ;regrid using existing weights file generated from the initial call of ESMF_regrid
		  printVarSummary(data_regrid)
                end if
                print("************regridding done!*******************")
;----------------------------------------------------------------------
; Write the regridded variable to a file using the "inefficient" method.
;----------------------------------------------------------------------
    		if(WRITE_RESULTS) then
      		  print("****************writing regridded data to nc file.....")
      		  rgrdFileName = var + "_" + ESM + "_regrided.nc"
      		  system("rm -f " + srcDir + rgrdFileName)
                  rgrd_file = addfile(srcDir + rgrdFileName,"c")

;---Create variable to hold global file attributes
      		  global = True
      		  copy_VarAtts(sfile, global)

      		  if (isatt(sfile,"title")) then
        	    global@TITLE = "REMAPPED: " + sfile@title
      		  end if
      	  	  global@remap         = "NCL: ESMF_regrid_with_weights (NCL version '" + \
                                         get_ncl_version() + "')"
      		  global@remap_method  = INTERP_METHOD
         	  global@creation_date = systemfunc("date")
    
      		  fileattdef( rgrd_file, global )        ; copy global file attributes
      		  filedimdef(rgrd_file,"time",-1,True)   ; force an unlimited dimension
;
; Write variables to file. Coordinate arrays will be written automatically
;
      		  rgrd_file->$var$     = data_regrid
    		end if

;----------------------------------------------------------------------
; Plot the original and regridded data.
;----------------------------------------------------------------------
    		if(PLOT_RESULTS) then  
      		  print("*******************plot original data and regrided data.......")
      		  wks = gsn_open_wks("ps",srcDir + var + "_" + ESM +"_regrid")     ; send graphics to PNG file

      		  res                     = True              ; Plot mods desired.
      		  res@gsnDraw             = False             ; We will panel later.
      		  res@gsnFrame            = False

      		  res@gsnMaximize         = True              ; Maximize plot  
      		  res@cnFillOn            = True              ; color plot desired
      		  res@cnFillPalette       = "gui_default"     ; set color map
      		  res@cnLinesOn           = False             ; turn off contour lines
      		  res@cnLineLabelsOn      = False             ; turn off contour lines
      		  res@cnLevelSelectionMode = "AutomaticLevels"
      		  res@lbLabelBarOn        = False             ; Labelbar will be in panel

      		  title_orig   = ESM + " original grid"
      		  title_regrid = "Regridded to CESM grid"

;---Resources for paneling
      		  pres                  = True
      		  pres@gsnMaximize      = True
      		  pres@gsnPanelLabelBar = True

;---Make a copy of each resource list
      		  res_orig                = res
      		  res_orig@gsnAddCyclic   = False
      		  res_orig@tiMainString   = title_orig

      		  res_regrid              = res
      		  res_regrid@gsnAddCyclic = True
      		  res_regrid@tiMainString = title_regrid

;---Create plot of original variable
      		  plot_orig = gsn_csm_contour_map(wks,data(6,:,:),res_orig)  ;6 is a random number: checking a snapshot 

;---Create plot of regridded variable
      		  plot_regrid = gsn_csm_contour_map(wks,data_regrid(6,:,:),res_regrid)

;---Draw both plots on one frame
      		  gsn_panel(wks,(/plot_orig,plot_regrid/),(/1,2/),pres)
    		end if
		print("*****ploting "+var+" done!")
                delete(data)
                delete(data_regrid)
            end do
            print(ESM+" done!")
        end do
        print(rcp + "done!!!")
     end do
     print("ALL DONE!!!")
end
    
    
