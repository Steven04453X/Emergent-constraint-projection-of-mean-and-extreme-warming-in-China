; =====================================================
; Author: Ziming Chen 
; This code is for:
;   - vis the relative contribution of uncertainty for each component 
;   1. read the var of uncertainty (mean and extreme)
;   2. cal the relative percentage 
;   3. vis, left for mean warming, right for extreme wamring 
; =====================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
begin
; =====================================================
;   1. read the var of uncertainty (mean and extreme)
s_DirOut = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Fig/07.UncerQualification/"
system("mkdir -p " + s_DirOut)
s_FigOut            = s_DirOut + "03.VisRelConForUncer"
wks                 = gsn_open_wks("X11", s_FigOut)
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/"
s_FileRead          = s_DirRead + "EachUncerContribution.nc"
FileMean            = addfile(s_FileRead, "r")
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/"
s_FileRead          = s_DirRead + "EachUncerContribution.nc"
FileExtreme         = addfile(s_FileRead, "r")
;
s_Exp               = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
r_TotalUncerMean    = new((/3, 4/), float)
r_TotalUncerMean!0  = "Terms"
r_TotalUncerMean&Terms = (/"Near", "Mid", "Long"/)
r_TotalUncerMean!1  = "SSP"
r_TotalUncerMean&SSP= (/"ssp126", "ssp245", "ssp370", "ssp585"/)
;
r_EC_UncerMean      = r_TotalUncerMean
r_Re_UncerMean      = r_TotalUncerMean
r_IV_UncerMean      = r_TotalUncerMean
r_IV_UncerProjMean  = r_TotalUncerMean
;
r_TotalUncerExtreme = r_TotalUncerMean
r_EC_UncerExtreme   = r_TotalUncerMean
r_Re_UncerExtreme   = r_TotalUncerMean
r_IV_UncerExtreme   = r_TotalUncerMean
r_IV_UncerProjExtreme = r_TotalUncerExtreme
;
do iExp = 0, dimsizes(r_TotalUncerMean(0, :))-1
    s_VarName2      = "TotalUncer_" + s_Exp(iExp)
    r_TotalUncerMean(:, iExp)    = FileMean->$s_VarName2$
    r_TotalUncerExtreme(:, iExp) = FileExtreme->$s_VarName2$
    ;
    s_VarName2      = "EC_Uncer_" + s_Exp(iExp)
    r_EC_UncerMean(:, iExp)      = FileMean->$s_VarName2$
    r_EC_UncerExtreme(:, iExp)   = FileExtreme->$s_VarName2$
    ;
    s_VarName2      = "Re_Uncer_" + s_Exp(iExp)
    r_Re_UncerMean(:, iExp)      = FileMean->$s_VarName2$
    r_Re_UncerExtreme(:, iExp)   = FileExtreme->$s_VarName2$
    if (s_Exp(iExp) .eq. "ssp370") then
        r_Re_UncerMean(:, iExp)  = where(r_Re_UncerMean(:, iExp) .lt. 0, \
                                         0, r_Re_UncerMean(:, iExp))
    end if
    ;
    s_VarName2      = "IV_Uncer_" + s_Exp(iExp)
    r_IV_UncerMean(:, iExp)      = FileMean->$s_VarName2$
    r_IV_UncerExtreme(:, iExp)   = FileExtreme->$s_VarName2$
    ;
    s_VarName2      = "IV_UncerInProj_" + s_Exp(iExp)
    r_IV_UncerProjMean(:, iExp)  = FileMean->$s_VarName2$
    r_IV_UncerProjExtreme(:,iExp)= FileExtreme->$s_VarName2$
    ;
    print(s_Exp(iExp) + "")
    print("Mean")
    print(sprintf("%4.2f", r_TotalUncerMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_EC_UncerMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerProjMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_Re_UncerMean(:, iExp)) )
    print("Extreme")
    print(sprintf("%4.2f", r_TotalUncerExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_EC_UncerExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerProjExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_Re_UncerExtreme(:, iExp)) )
    print("")
end do
;
; =====================================================
;   2. cal the relative percentage 
r_EC_UncerMean      = r_EC_UncerMean / r_TotalUncerMean * 100.
r_Re_UncerMean      = r_Re_UncerMean / r_TotalUncerMean * 100.
r_IV_UncerMean      = r_IV_UncerMean / r_TotalUncerMean * 100.
r_IV_UncerProjMean  = r_IV_UncerProjMean / r_TotalUncerMean * 100.
;
r_EC_UncerExtreme   = r_EC_UncerExtreme / r_TotalUncerExtreme * 100.
r_Re_UncerExtreme   = r_Re_UncerExtreme / r_TotalUncerExtreme * 100.
r_IV_UncerExtreme   = r_IV_UncerExtreme / r_TotalUncerExtreme * 100.
r_IV_UncerProjExtreme   = r_IV_UncerProjExtreme / r_TotalUncerExtreme * 100.
;
do iExp = 0, dimsizes(s_Exp)-1
    print(s_Exp(iExp) + ", from near to long")
    print("Mean Temp (%):")
    print("EC Coef Uncer: " + sprintf("%4.2f", r_EC_UncerMean(:, iExp)))
    print("PD IV Uncer: " + sprintf("%4.2f", r_IV_UncerMean(:, iExp)))
    print("Proj IV Uncer: " + sprintf("%4.2f", r_IV_UncerProjMean(:, iExp)))
    print("Re Uncer: " + sprintf("%4.2f", r_Re_UncerMean(:, iExp)))
    r_Temp          := r_IV_UncerMean(:, iExp) + r_IV_UncerProjMean(:, iExp)
    print("Total con of IV: " + sprintf("%4.2f", r_Temp))
    print("")
    print("Extreme Temp (%):")
    print("EC Coef Uncer: " + sprintf("%4.2f", r_EC_UncerExtreme(:, iExp)))
    print("PD IV Uncer: " + sprintf("%4.2f", r_IV_UncerExtreme(:, iExp)))
    print("Proj IV Uncer: " + sprintf("%4.2f", r_IV_UncerProjExtreme(:, iExp)))
    print("Re Uncer: " + sprintf("%4.2f", r_Re_UncerExtreme(:, iExp)))
    r_Temp          := r_IV_UncerExtreme(:, iExp) + r_IV_UncerProjExtreme(:, iExp)
    print("Total con of IV: " + sprintf("%4.2f", r_Temp))
    print("")
    print("")
end do
print("Ave across four emission")
r_EC_UncerMeanAve    = dim_avg_n_Wrap(r_EC_UncerMean, 1)
r_Re_UncerMeanAve    = dim_avg_n_Wrap(r_Re_UncerMean, 1)
r_IV_UncerMeanAve    = dim_avg_n_Wrap(r_IV_UncerMean, 1)
r_IV_UncerProjMeanAve= dim_avg_n_Wrap(r_IV_UncerProjMean, 1)

r_EC_UncerExtremeAve    = dim_avg_n_Wrap(r_EC_UncerExtreme, 1)
r_Re_UncerExtremeAve    = dim_avg_n_Wrap(r_Re_UncerExtreme, 1)
r_IV_UncerExtremeAve    = dim_avg_n_Wrap(r_IV_UncerExtreme, 1)
r_IV_UncerProjExtremeAve= dim_avg_n_Wrap(r_IV_UncerProjExtreme, 1)
;
print("Mean Temp (%):")
print("EC Coef Uncer: " + sprintf("%4.2f", r_EC_UncerMeanAve))
print("PD IV Uncer: " + sprintf("%4.2f", r_IV_UncerMeanAve))
print("Proj IV Uncer: " + sprintf("%4.2f", r_IV_UncerProjMeanAve))
print("Re Uncer: " + sprintf("%4.2f", r_Re_UncerMeanAve))
r_Temp              := r_IV_UncerMeanAve + r_IV_UncerProjMeanAve
print("Total con of IV: " + sprintf("%4.2f", r_Temp))
print("")
print("Extreme Temp (%):")
print("EC Coef Uncer: " + sprintf("%4.2f", r_EC_UncerExtremeAve))
print("PD IV Uncer: " + sprintf("%4.2f", r_IV_UncerExtremeAve))
print("Proj IV Uncer: " + sprintf("%4.2f", r_IV_UncerProjExtremeAve))
print("Re Uncer: " + sprintf("%4.2f", r_Re_UncerExtremeAve))
r_Temp              := r_IV_UncerExtremeAve + r_IV_UncerProjExtremeAve
print("Total con of IV: " + sprintf("%4.2f", r_Temp))
; exit()
; =====================================================
;   3. vis, left for mean warming, right for extreme wamring 
    sres                                  = True
    
    sres@trXMinF                          = -0.3
    sres@trXMaxF                          = 5.5
    sres@trYMaxF                          = 100. ;max(r_G11) + 0.1
    sres@trYMinF                          = 0. ;min(r_G11) - 0.1
    
    sres@gsnFrame                         = False
    sres@gsnDraw                          = False
    sres@gsnXYBarChart                    = True     ;; 设定柱状图
    sres@gsnXYBarChartBarWidth            = 0.25
    
    sres@gsnRightString                   = ""
    sres@gsnCenterString                  = ""
    sres@gsnLeftStringOrthogonalPosF      = 0.02
    sres@gsnRightStringFontHeightF        = 0.01
    sres@gsnLeftStringFontHeightF         = 0.013
    sres@gsnYRefLine                      = 0.
    
    sres@tiYAxisString                    = "Relative Contribution (%)"
    sres@tiYAxisFontHeightF               = 0.012
    sres@tmXTOn                           = False
    sres@tmYROn                           = False
    sres@tmXBMode                         = "Explicit"
    ; sres@tmXBValues                       = fspan(1., 6.5, dimsizes(s_VarName))
    s_VarNameTitle = (/"SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"/)
    sres@tmXBLabels                       = s_VarNameTitle
    sres@tmXBLabelAngleF                  = 90.
    sres@tmXBLabelFontHeightF             = 0.01
    sres@tmYLLabelFontHeightF             = 0.012

    txres                                 = True
    txres@txFontHeightF                   = 0.01
    txres@txJust                          = "CenterLeft"
    
    sres@vpXF                             = 0.1 ; 0.3~0.95
    sres@vpYF                             = 0.8
    sres@vpWidthF                         = 0.3
    sres@vpHeightF                        = 0.25
    ;
    plot1                                 = new(3, graphic)
    plot1_EC                              = plot1
    plot1_Re                              = plot1
    plot1_ProjIV                          = plot1
    xx_res2                               = (/0, 0.4, 0.8, 1.2/)
    dxx_res2                              = 2.
    s_Term                                = (/ "Near", "Middle", "Long" /)
    do iTerms = 0, 2
        sres@tmXBValues                   = xx_res2
        x                                 = sres@tmXBValues
        sres@gsnLeftString                = ""
        ;
        ;; Other Uncer 
        r_Temp                           := r_IV_UncerMean(iTerms, :)
        r_Temp                           := r_IV_UncerMean(iTerms, :) + \
                                            r_IV_UncerProjMean(iTerms, :) + \
                                            r_EC_UncerMean(iTerms, :) + \
                                            r_Re_UncerMean(iTerms, :) 
        ; print(r_Temp + "")
        ; print("")
        ; exit()
        sres@gsnXYBarChartColors          = "blue"
        sres@tmXBLabels                   = s_VarNameTitle
        sres@tmXBOn                       = True
        plot1_Re(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; EC Uncer
        r_Temp                            = r_IV_UncerMean(iTerms, :) + \
                                            r_IV_UncerProjMean(iTerms, :) + \
                                            r_EC_UncerMean(iTerms, :)
        ; print(r_Temp + "")
        ; print("")
        sres@gsnXYBarChartColors          = "green3"
        plot1_EC(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; Proj IV 
        r_Temp                            = r_IV_UncerMean(iTerms, :) + \
                                            r_IV_UncerProjMean(iTerms, :)
        sres@gsnXYBarChartColors          = "gold3"
        plot1_ProjIV(iTerms)              = gsn_csm_xy(wks, x, r_Temp, sres)
        ; print(r_Temp + "")
        ; print("")
        ;
        ;; PD IV 
        sres@gsnXYBarChartColors          = "orange"
        if (iTerms .eq. 0) then
            sres@gsnLeftString            = "(a) Annual Mean Warming"
        end if
        ;
        plot1(iTerms)                     = gsn_csm_xy(wks, x, \
                                            r_IV_UncerMean(iTerms, :), sres)
        ; print(r_IV_UncerMean(iTerms, :) + "")
        ; print("")
        ;
        r_xText                           = (/0.14, 0.24, 0.35/)
        r_yText                           = 0.47
        gsn_text_ndc(wks, s_Term(iTerms), r_xText(iTerms), r_yText, txres)
        ;
        overlay(plot1_Re(iTerms), plot1_Re(iTerms))
        overlay(plot1_Re(iTerms), plot1_EC(iTerms))
        overlay(plot1_Re(iTerms), plot1_ProjIV(iTerms))
        overlay(plot1_Re(iTerms), plot1(iTerms))
        draw(plot1_Re(iTerms))
        sres@tmYLOn                       = False
        sres@tiYAxisOn                    = False
        ;
        xx_res2                           = xx_res2 + dxx_res2 ;2.5 
    end do
    ;
    sres@vpXF                             = sres@vpXF + sres@vpWidthF + 0.02 ; 0.3~0.95
    plot1                                 = new(3, graphic)
    plot1_EC                              = plot1
    plot1_Re                              = plot1
    plot1_ProjIV                          = plot1
    xx_res2                               = (/0, 0.4, 0.8, 1.2/)
    dxx_res2                              = 2.
    do iTerms = 0, 2
        sres@tmXBValues                   = xx_res2
        x                                 = sres@tmXBValues
        sres@gsnLeftString                = ""
        ;
        ;; Other Uncer 
        r_Temp                           := r_IV_UncerExtreme(iTerms, :)
        r_Temp                           := r_IV_UncerExtreme(iTerms, :) + \
                                            r_IV_UncerProjExtreme(iTerms, :) + \
                                            r_EC_UncerExtreme(iTerms, :) + \
                                            r_Re_UncerExtreme(iTerms, :) 
        ; print(r_Temp + "")
        ; print("")
        ; exit()
        sres@gsnXYBarChartColors          = "blue"
        sres@tmXBLabels                   = s_VarNameTitle
        sres@tmXBOn                       = True
        plot1_Re(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; EC Uncer
        r_Temp                            = r_IV_UncerExtreme(iTerms, :) + \
                                            r_IV_UncerProjExtreme(iTerms, :) + \
                                            r_EC_UncerExtreme(iTerms, :)
        ; print(r_Temp + "")
        ; print("")
        sres@gsnXYBarChartColors          = "green3"
        plot1_EC(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        r_Temp                            = r_IV_UncerExtreme(iTerms, :) + \
                                            r_IV_UncerProjExtreme(iTerms, :) 
        ; print(r_Temp + "")
        ; print("")
        sres@gsnXYBarChartColors          = "gold3" ; "yellow3"
        plot1_ProjIV(iTerms)              = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; IV 
        sres@gsnXYBarChartColors          = "orange"
        if (iTerms .eq. 0) then
            sres@gsnLeftString            = "(b) Extreme Warming"
            sres@gsnLeftStringOrthogonalPosF = 0.05
        end if
        plot1(iTerms)                     = gsn_csm_xy(wks, x, \
                                            r_IV_UncerExtreme(iTerms, :), sres)
        ; print(r_IV_UncerExtreme(iTerms, :) + "")
        ; print("")
        ;
        r_xText                           = (/0.46, 0.56, 0.665/)
        r_yText                           = 0.47
        gsn_text_ndc(wks, s_Term(iTerms), r_xText(iTerms), r_yText, txres)
        ;
        overlay(plot1_Re(iTerms), plot1_Re(iTerms))
        overlay(plot1_Re(iTerms), plot1_EC(iTerms))
        overlay(plot1_Re(iTerms), plot1_ProjIV(iTerms))
        overlay(plot1_Re(iTerms), plot1(iTerms))
        draw(plot1_Re(iTerms))
        sres@tmYLOn                       = False
        sres@tiYAxisOn                    = False
        ;
        xx_res2                           = xx_res2 + dxx_res2 ;2.5 
    end do
    ;
    txres@txFontColor                     = "blue"
    s_TermTemp                            = "Var Unrelated to EC"
    gsn_text_ndc(wks, s_TermTemp, 0.73, 0.63, txres)
    ;
    txres@txFontColor                     = "green3"
    s_TermTemp                            = "Var Explained by EC"
    gsn_text_ndc(wks, s_TermTemp, 0.73, 0.605, txres)
    ;
    txres@txFontColor                     = "gold3" ; "yellow3"
    s_TermTemp                            = "Internal Variability in Projection"
    gsn_text_ndc(wks, s_TermTemp, 0.73, 0.585, txres)
    ;
    txres@txFontColor                     = "orange"
    s_TermTemp                            = "Internal Variability in Historical"
    gsn_text_ndc(wks, s_TermTemp, 0.73, 0.565, txres)
    ;
    ; drawNDCGrid(wks)
    frame(wks)
    
;
end








