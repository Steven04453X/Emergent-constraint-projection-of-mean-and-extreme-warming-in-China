; ===============================================================
; Author: Ziming Chen 
; This code is for:
;   - quantify how the constrained projections captures the spread correctly
;   1. indicate models, exp, var name, period 
;   2. for each model, read the raw projection which works as a pseudo warming at first,
; and then cal the standard deviation in this period 
;   3. read the constrained projection and then cal the inter-model spread 
;   4. cal how much the spread of constrained projection can cover the original spread
;   5. vis: 4 x 1
; ===============================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
load "$WACCEM/Work/Others/Special_Code/my_DrawFunction.ncl"
;
begin
;
; ===============================================================
;   1. indicate models, exp, var name, period 
s_ModelSelected           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models               = "_18Models"
;
if (isdefined("s_Exp") .eq. False) then
    s_Exp               = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
end if
;
i_yrBaseState           = (/1995, 2014/)
i_yrProj                = (/2080, 2099/)
i_Term                  = 2
;
s_VarName               = "TXx"
;
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Fig/06.PerfectModelTest/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FigOut                = s_DirOut + "07.Vis_ConstrainedSpread_PseudoSpread_" +\
                            s_VarName
wks                     = gsn_open_wks("pdf", s_FigOut)
;
; ===============================================================
;   2. for each model, read the raw projection which works as a pseudo warming at first,
; and then cal the standard deviation in this period 
do iExp = 0, dimsizes(s_Exp)-1
    do iModel = 0, dimsizes(s_ModelSelected)-1
        if (s_Exp(iExp) .eq. "ssp126" .and. \
            s_ModelSelected(iModel) .eq. "ACCESS-ESM1-5" ) then
            continue
        end if
        ;
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp(iExp), "day")
        s_ExpModel_Path:= l_list[0]
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp(iExp), "day")
        s_ExpModel_Path_:= l_list_[0]
        s_Real_         := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
            s_Real      := s_Real_
            s_ExpModel_Path := s_ExpModel_Path_
        end if
        ;
        ;; raw projection 
        s_DirRead        = "$WACCEM/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
        s_FileRead       = s_DirRead + s_VarName + "_AnMean_historical_" + \
                           s_ModelSelected(iModel) + "_CMIP6_AllReal.nc"
        ls              := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print("No historical " + s_ModelSelected(iModel))
            continue
        end if
        FileHist        := addfile(s_FileRead, "r")
        ;
        s_FileRead       = s_DirRead + s_VarName + "_AnMean_" + s_Exp(iExp) + \
                           "_" + s_ModelSelected(iModel) + "_CMIP6_AllReal.nc"
        ls              := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print("No " + s_Exp(iExp) + " " + s_ModelSelected(iModel))
            continue
        end if
        FileProj        := addfile(s_FileRead, "r")
        ;
        do iReal = 0, dimsizes(s_Real)-1
            s_VarName2   = s_VarName + "_" + s_ModelSelected(iModel) + \
                           "_historical_" + s_Real(iReal)
            if (my_VarNameExist(FileHist, s_VarName2) .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Hist      := FileHist->$s_VarName2$(\
                        {i_yrBaseState(0):i_yrBaseState(1)}, {-90:90}, {0:360})
            r_Hist      := dim_avg_n_Wrap(r_Hist, 0)
            ;
            s_VarName2   = s_VarName + "_" + s_ModelSelected(iModel) + "_" + \
                           s_Exp(iExp) + "_" + s_Real(iReal)
            if (my_VarNameExist(FileProj, s_VarName2) .eq. False) then
                print(s_VarName2 + "  Missing!")
                continue
            end if
            r_Proj      := FileProj->$s_VarName2$({i_yrProj(0):i_yrProj(1)}, \
                                                  {-90:90}, {0:360})
            ;
            r_Hist      := conform(r_Proj, r_Hist, (/1, 2/))
            copy_VarCoords(r_Proj, r_Hist)
            r_Proj       = r_Proj - r_Hist
            ;
            r_Proj      := my_SelChina([/r_Proj/])
            r_Proj      := my_AreaMean2([/r_Proj/], r_Proj&lat, r_Proj&lon, 0)
            ;
            if (isdefined("r_PseudoAll") .eq. False .or. iReal .eq. 0) then
                r_PseudoAll := new((/dimsizes(s_Real), 3/), float)
                r_PseudoAll!0= "Real"
                r_PseudoAll&Real = s_Real
                r_PseudoAll!1= "Mean_STD1_STD2"
                r_PseudoAll&Mean_STD1_STD2 = (/0, 1, 2/)
            end if
            ; r_PseudoAll(iReal, 0)= avg(r_Proj)
            ; r_PseudoAll(iReal, 1)= r_PseudoAll(iReal, 0) - stddev(r_Proj) * 1.64
            ; r_PseudoAll(iReal, 2)= r_PseudoAll(iReal, 0) + stddev(r_Proj) * 1.64
            i_Ind                := ind(ismissing(r_Proj) .eq. False)
            ; printVarSummary(r_Proj)
            ; printMinMax(i_Ind, False)
            if (all(ismissing(i_Ind)) ) then
                continue
            end if
            r_Proj               := r_Proj(i_Ind)
            r_PseudoAll(iReal, :) = my_QuaterNum_withMean3(r_Proj, 5, 95)
        end do
        r_PseudoAll      := dim_avg_n_Wrap(r_PseudoAll, 0)
        ;
        ; ===============================================================
        ;   3. read the constrained projection and then cal the inter-model spread 
        ;; constrained warming 
        s_DirRead := "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
        s_FileRead       := s_DirRead + "FurtherCorrectConstrainedWarming_PseModel_"\
                            + s_ModelSelected(iModel) + "_" + s_Exp(iExp) + \
                            s_4Models + ".nc"
        ls               := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print("No constrained in " + s_Exp(iExp) + " for " + \
                  s_ModelSelected(iModel))
            continue
        end if
        FileCon          := addfile(s_FileRead, "r")
        do iMFC = 0, dimsizes(s_ModelSelected)-1
            if (s_ModelSelected(iMFC) .eq. s_ModelSelected(iModel)) then
                continue
            end if
            s_VarName2   := "ConstrainedWarming_" + s_ModelSelected(iMFC) + "_"\
                            + s_Exp(iExp)
            if (my_VarNameExist(FileCon, s_VarName2) .eq. False) then
                print(s_VarName2 + " Missing!")
                continue
            end if
            r_ConWarm    := FileCon->$s_VarName2$(i_Term, {-90:90}, {0:360})
            ; printVarSummary(r_ConWarm)
            ;
            r_ConWarm    := my_SelChina([/r_ConWarm/])
            r_ConWarm    := my_AreaMean2([/r_ConWarm/], r_ConWarm&lat, \
                                           r_ConWarm&lon, 0)
            ;
            if (isdefined("r_ConWarmAll") .eq. False .or. iMFC .eq. 0) then
                r_ConWarmAll := new(dimsizes(s_ModelSelected), float)
                r_ConWarmAll!0= "Models"
                r_ConWarmAll&Models = s_ModelSelected
            end if
            r_ConWarmAll(iMFC)= r_ConWarm
            ;
        end do
        if (isdefined("r_ConWarmAll") .eq. False) then
            continue
        end if
        print(s_Exp(iExp) + ", " + s_ModelSelected(iModel))
        ; r_ConWarmMean    := avg(r_ConWarmAll)
        ; r_ConWarmSTD     := stddev(r_ConWarmAll) * 1.64
        ; r_ConWarmAll     := (/r_ConWarmMean, \
        ;                       r_ConWarmMean - r_ConWarmSTD, \
        ;                       r_ConWarmMean + r_ConWarmSTD/)
        i_Ind              := ind(ismissing(r_ConWarmAll) .eq. False)
        r_ConWarmAll       := r_ConWarmAll(i_Ind)
        r_ConWarmAll       := my_QuaterNum_withMean3(r_ConWarmAll, 5, 95)
        ; ===============================================================
        ;   4. cal how much the spread of constrained projection can cover the original spread
        if (r_ConWarmAll(1) .le. r_PseudoAll(1) .and. \
            r_ConWarmAll(2) .ge. r_PseudoAll(2)) then
            ;; all contain 
            r_Ratio      := (r_PseudoAll(2) - r_PseudoAll(1)) / \
                            (r_ConWarmAll(2) - r_ConWarmAll(1))
        else if (r_ConWarmAll(1) .gt. r_PseudoAll(1) .and. \
                 r_ConWarmAll(1) .lt. r_PseudoAll(2) .and. \
                 r_ConWarmAll(2) .gt. r_PseudoAll(2)) then
            ;; Constrained > Pseudo 
            r_Ratio      := r_PseudoAll(2) - r_ConWarmAll(1)
            r_Ratio       = r_Ratio / (r_ConWarmAll(2) - r_ConWarmAll(1))
        else if (r_ConWarmAll(2) .gt. r_PseudoAll(1) .and. \
                 r_ConWarmAll(2) .lt. r_PseudoAll(2) .and. \
                 r_ConWarmAll(1) .lt. r_PseudoAll(1)) then 
            ;; Constrained < Pseudo 
            r_Ratio      := r_ConWarmAll(2) - r_PseudoAll(1)
            r_Ratio       = r_Ratio / (r_ConWarmAll(2) - r_ConWarmAll(1))
        else if (r_ConWarmAll(1) .ge. r_PseudoAll(1) .and.\
                 r_ConWarmAll(2) .le. r_PseudoAll(2)) then
            ;; Pse spread contain constrained spread 
            r_Ratio       = 1
        else 
            r_Ratio       = 0
        end if 
        end if
        end if
        end if
        ;
        print("Ratio: " + r_Ratio + " of " + s_ModelSelected(iModel))
        print("Constrained: " + sprintf("%4.2f", r_ConWarmAll(1)) + " to " + \
              sprintf("%4.2f", r_ConWarmAll(2)) )
        print("Pseudo: " + sprintf("%4.2f", r_PseudoAll(1)) + " to " + \
              sprintf("%4.2f", r_PseudoAll(2)))
        print("")
        ;
        if (isdefined("r_RatioCon") .eq. False .or. iModel .eq. 0) then
            r_RatioCon   := new((/dimsizes(s_ModelSelected)/), float)
            r_RatioCon!0  = "Models"
            r_RatioCon&Models   = s_ModelSelected 
            ;
            r_ConWarmAllModels := new((/dimsizes(s_ModelSelected), 3/), float)
            copy_VarCoords(r_RatioCon, r_ConWarmAllModels(:, 0))
            r_ConWarmAllModels!1= "Mean_STD1_STD2"
            r_ConWarmAllModels&Mean_STD1_STD2 = (/0, 1, 2/)
            ;
            r_PseudoAllModels  := r_ConWarmAllModels
        end if
        r_RatioCon(iModel)= r_Ratio
        r_ConWarmAllModels(iModel, :) = r_ConWarmAll
        r_PseudoAllModels(iModel, :)  = r_PseudoAll
        ;
        delete([/r_ConWarmAll/])
        delete([/r_PseudoAll/])
        ;
    end do
    print(s_ModelSelected + ": " + sprintf("%4.2f", r_RatioCon))
    ;
    ;; Raw projection 
    r_RawProj                        := r_PseudoAllModels(:, 0)
    i_Ind                            := ind(ismissing(r_RawProj) .eq. False)
    r_RawProj                        := r_RawProj(i_Ind)
    r_RawProj                        := my_QuaterNum_withMean3(r_RawProj, 5, 95)
    r_RawProj                        := conform_dims((/2, 3/), r_RawProj, 1)
    ; ===============================================================
    ;   5. vis: 4 x 1
    sres                                 := True
    
    sres@trXMinF                          = 0.
    sres@trXMaxF                          = dimsizes(s_ModelSelected) + 2
    
    sres@gsnFrame                         = False
    sres@gsnDraw                          = False
    
    sres@gsnRightString                   = ""
    sres@gsnCenterString                  = ""
    sres@gsnLeftStringOrthogonalPosF      = 0.02
    sres@gsnRightStringOrthogonalPosF     = 0.01
    sres@gsnRightStringFontHeightF        = 0.015
    sres@gsnLeftStringFontHeightF         = 0.02
    sres@gsnYRefLine                      = 0.
    
    if (s_Exp(iExp) .eq. "ssp126") then
        s_ExpTitle                        = "SSP1-2.6"
        sres@trYMaxF                      = 4.5 ;
        sres@trYMinF                      = 0.1 ;
    else if (s_Exp(iExp) .eq. "ssp245") then
        s_ExpTitle                        = "SSP2-4.5"
        sres@trYMaxF                      = 7 ;
        sres@trYMinF                      = 1 ;
    else if (s_Exp(iExp) .eq. "ssp370") then
        s_ExpTitle                        = "SSP3-7.0"
        sres@trYMaxF                      = 10 ;
        sres@trYMinF                      = 1.5 ;
    else if (s_Exp(iExp) .eq. "ssp585") then
        s_ExpTitle                        = "SSP5-8.5"
        sres@trYMaxF                      = 12 ;
        sres@trYMinF                      = 2 ;
    end if
    end if
    end if
    end if
    s_OrderAll                            = (/"(a)", "(b)", "(c)", "(d)"/)
    sres@gsnLeftString                    = s_OrderAll(iExp) + \
                                            " Leave-one-out Test For TXx Warming"
    sres@gsnRightString                   = s_ExpTitle
    sres@tiYAxisString                    = "Warming Relative to 1995-2014"
    sres@tiYAxisFontHeightF               = 0.012
    sres@tmXTOn                           = False
    sres@tmYROn                           = False
    sres@tmXBMode                         = "Explicit"
    ;
    s_VarNameTitle                        = array_append_record("Raw", \
                                                            s_ModelSelected, 0)
    sres@tmXBLabels                       = s_VarNameTitle
    sres@tmXBLabelAngleF                  = 90.
    sres@tmXBLabelFontHeightF             = 0.011
    sres@tmYLLabelFontHeightF             = 0.012
    ;
    ; sres@xyMarkLineModes                  = "Markers"
    ; sres@xyMarkers                        = 4
    ; sres@xyMarkerThicknesses              = 1
    ; sres@xyMarkerSizeF                    = 0.004 

    txres                                 = True
    txres@txFontHeightF                   = 0.01
    
    sres@vpXF                             = 0.1 ; 0.3~0.95
    sres@vpYF                             = 0.8
    sres@vpWidthF                         = 0.7
    sres@vpHeightF                        = 0.25
    ;
    r_TempX                               = ispan(1, dimsizes(s_VarNameTitle), 1) * 1.
    r_TempY                               = r_TempX * 1.
    r_TempY                               = 0
    sres@tmXBValues                       = r_TempX
    ;
    if (iExp .eq. 0) then
        plot                              = new(dimsizes(s_Exp), graphic)
        plot1                             = new((/dimsizes(s_Exp), \
                                            dimsizes(s_ModelSelected)+1, 3/), \
                                            graphic)
        plot1_EC                          = plot1
        plot1_                            = plot1(:, 0, :)
        plotMarkers                       = new((/dimsizes(s_Exp), 3/), graphic)
    end if
    plot(iExp)                            = gsn_csm_xy(wks, r_TempX, r_TempY, \
                                                       sres)
    ;
    dxx                                   = 0.2
    resErrorBar                           = True
    resErrorBar@HorizontalLineLength      = 0.1
    resErrorBar@gsLineThicknessF          = 2
    ;
    r_PseudoAllModels                     = where(ismissing(r_PseudoAllModels),\
                                                  0, r_PseudoAllModels)
    r_ConWarmAllModels                    = where(ismissing(r_ConWarmAllModels),\
                                                  0, r_ConWarmAllModels)
    resErrorBar@gsLineColor               = (/0, 0, 0/) / 255.
    plot1(iExp, 1:, :)                    = add_yerrorbar_Function(wks, plot(iExp), \
                                    r_TempX(1:)+dxx, r_PseudoAllModels(:, 0), \
                                    r_PseudoAllModels(:, 1:), resErrorBar)
    resErrorBar@gsLineColor               = (/112, 160, 205/) / 255.
    plot1_EC(iExp, 1:, :)                 = add_yerrorbar_Function(wks, plot(iExp), \
                                    r_TempX(1:)-dxx, r_ConWarmAllModels(:, 0), \
                                    r_ConWarmAllModels(:, 1:), resErrorBar)
    ;
    ;; MME 
    resErrorBar@gsLineColor               = (/153, 0, 2/) / 255.
    plot1_Temp                           := add_yerrorbar_Function(wks, \
                                    plot(iExp), r_TempX(0:1), r_RawProj(0, 0), \
                                    r_RawProj(0:0, 1:), resErrorBar)
    plot1_(iExp, :)                       = plot1_Temp(0, :)
    ; printVarSummary(plot1_)
    ;
    ;; markers 
    res3                                 := sres
    res3@tmLabelAutoStride                = True
    res3@gsMarkerSizeF                    = 0.018
    res3@gsMarkerThicknessF               = 2.0
    res3@gsMarkerIndex                    = 1
    res3@gsMarkerColor                    = (/153, 0, 2/) / 255.
    plotMarkers(iExp, 0)                  = gsn_add_polymarker(wks, plot(iExp),\
                                            r_TempX(0), r_RawProj(0, 0), res3)
    res3@gsMarkerColor                    = (/0, 0, 0/) / 255.
    plotMarkers(iExp, 1)                  = gsn_add_polymarker(wks, plot(iExp),\
                                r_TempX(1:)+dxx, r_PseudoAllModels(:, 0), res3)
    res3@gsMarkerColor                    = (/112, 160, 205/) / 255.
    plotMarkers(iExp, 2)                  = gsn_add_polymarker(wks, plot(iExp),\
                                r_TempX(1:)-dxx, r_ConWarmAllModels(:, 0), res3)
    ; do iModel = 1, dimsizes(s_ModelSelected)-1  
    ;     res3@gsMarkerColor                = (/0, 0, 0/) / 255.
    ;     plotMarkers(iExp, iModel, 0)      = gsn_add_polymarker(wks, plot, \
    ;                 r_TempX(iModel)+dxx, r_PseudoAllModels(iModel-1, 0), sres)
    ;     res3@gsMarkerColor                = (/112, 160, 205/) / 255.
    ;     plotMarkers(iExp, iModel, 1)      = gsn_add_polymarker(wks, plot, \
    ;                 r_TempX(iModel)-dxx, r_ConWarmAllModels(iModel-1, 0), sres)
    ; end do 
    ; overlay(plot, plotMarkers(iExp, 0, 0))
    ; overlay(plot, plotMarkers(iExp, 1:, 0))
    ;
    ; Create legend
      lgres                     = True
      lgres@lgLineColors       := (/(/0, 0, 0/), (/112, 160, 205/), \
                                   (/153, 0, 2/)/) / 255.
      lgres@lgLineThicknessF    = (/5, 5, 5/) 
      lgres@lgLabelFontHeightF  = .15           ; set the legend label font thickness
      lgres@vpWidthF            = 0.2          ; width of legend (NDC)
      lgres@vpHeightF           = 0.08            ; height of legend (NDC)，调节Legend之间的垂直间距
      lgres@lgMonoDashIndex     = True           ; Same line dash patterns  ( 线型属性是lgDashIndexes )
      lgres@lgPerimOn           = False          ; draw the box perimeter in orange
      labels = "   " + (/ "Raw Projection", "Constrainted Projection", \
                          "Pseudo Warming" /)
      lbid = gsn_create_legend(wks, 3, labels(::-1), lgres)
    
      ; Add legend to bottom right corner of plot
      amres                     = True
      amres@amJust              = "TopLeft"  ; Corner of plot for positioning
      amres@amOrthogonalPosF    =  -0.48   ;  0.5 is the bottom edge of the plot.
      amres@amParallelPosF      =  -0.45   ;  0.5 is the right edge of the plot.
      annoid  = gsn_add_annotation(plot(iExp), lbid, amres)

    ;
    draw(plot(iExp))
    frame(wks)
    ;
end do
;
resP                            = True
resP@gsnMaximize                = False
gsn_panel(wks, plot, (/2, 2/), resP)
;
end


; (0) ACCESS-CM2: 1.00
; (1) ACCESS-ESM1-5: 0.40
; (2) AWI-CM-1-1-MR: 0.21
; (3) BCC-CSM2-MR: 0.67
; (4) CAMS-CSM1-0: 1.00
; (5) CNRM-CM6-1: 1.00
; (6) CNRM-ESM2-1: 0.32
; (7) CanESM5: 1.00
; (8) CanESM5-CanOE: 1.00
; (9) EC-Earth3: 0.00
; (10)    EC-Earth3-Veg: 0.97
; (11)    FGOALS-f3-L: 0.00
; (12)    FGOALS-g3: 0.21
; (13)    GFDL-CM4: 0.68
; (14)    GFDL-ESM4: 0.00
; (15)    GISS-E2-1-G: 1.00
; (16)    HadGEM3-GC31-LL: 0.11
; (17)    INM-CM4-8: 0.07
; (18)    IPSL-CM6A-LR: 1.00
; (19)    MIROC6: 0.74
; (20)    MPI-ESM1-2-LR: 0.40
; (21)    MRI-ESM2-0: 0.69
; (22)    UKESM1-0-LL: 0.97








