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
s_FigOut            = s_DirOut + "04.1.VisRelConForUncer_InclReducedUncer_de"
wks                 = gsn_open_wks("X11", s_FigOut)
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/"
s_FileRead          = s_DirRead + "EachUncerContribution_InclReducedUncer.nc"
FileMean            = addfile(s_FileRead, "r")
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/"
s_FileRead          = s_DirRead + "EachUncerContribution_InclReducedUncer.nc"
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
r_EC_CoefUncerMean  = r_TotalUncerMean
r_Re_UncerMean      = r_TotalUncerMean
r_IV_UncerMean      = r_TotalUncerMean
r_IV_UncerProjMean  = r_TotalUncerMean
;
r_TotalUncerExtreme = r_TotalUncerMean
r_EC_UncerExtreme   = r_TotalUncerMean
r_EC_CoefUncerExtreme = r_TotalUncerMean
r_Re_UncerExtreme   = r_TotalUncerMean
r_IV_UncerExtreme   = r_TotalUncerMean
r_IV_UncerProjExtreme = r_TotalUncerExtreme
;
do iExp = 0, dimsizes(r_TotalUncerMean(0, :))-1
    s_VarName2      = "TotalUncer_" + s_Exp(iExp)
    r_TotalUncerMean(:, iExp)    = FileMean->$s_VarName2$
    r_TotalUncerExtreme(:, iExp) = FileExtreme->$s_VarName2$
    ;
    s_VarName2      = "EC_Redu_Uncer_" + s_Exp(iExp)
    r_EC_UncerMean(:, iExp)      = FileMean->$s_VarName2$
    r_EC_UncerExtreme(:, iExp)   = FileExtreme->$s_VarName2$
    ;
    s_VarName2      = "EC_Coef_Uncer_" + s_Exp(iExp)
    r_EC_CoefUncerMean(:, iExp)  = FileMean->$s_VarName2$
    r_EC_CoefUncerExtreme(:,iExp)= FileExtreme->$s_VarName2$
    ;
    s_VarName2      = "Re_Uncer_" + s_Exp(iExp)
    r_Re_UncerMean(:, iExp)      = FileMean->$s_VarName2$
    r_Re_UncerExtreme(:, iExp)   = FileExtreme->$s_VarName2$
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
r_EC_CoefUncerMean  = r_EC_CoefUncerMean / r_TotalUncerMean * 100.
;
r_EC_UncerExtreme   = r_EC_UncerExtreme / r_TotalUncerExtreme * 100.
r_Re_UncerExtreme   = r_Re_UncerExtreme / r_TotalUncerExtreme * 100.
r_IV_UncerExtreme   = r_IV_UncerExtreme / r_TotalUncerExtreme * 100.
r_IV_UncerProjExtreme   = r_IV_UncerProjExtreme / r_TotalUncerExtreme * 100.
r_EC_CoefUncerExtreme  = r_EC_CoefUncerExtreme / r_TotalUncerExtreme * 100.
;
do iExp = 0, dimsizes(s_Exp)-1
    print("Var percentage in " + s_Exp(iExp) + "")
    r_Temp         := r_EC_UncerMean(:, iExp) + r_IV_UncerMean(:, iExp) + \
                        r_IV_UncerProjMean(:, iExp) + r_Re_UncerMean(:, iExp)
    print("Mean")
    print(sprintf("%4.2f", r_EC_UncerMean(:, iExp)-r_EC_CoefUncerMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerProjMean(:, iExp)) + ", " + \
          sprintf("%4.2f", r_Re_UncerMean(:, iExp)) + " = " + r_Temp )
    print("Con Coef Var: " + sprintf("%4.2f", r_EC_CoefUncerMean(:, iExp)))
    print("Extreme")
    r_Temp         := r_EC_UncerExtreme(:, iExp) + r_IV_UncerExtreme(:, iExp) +\
                      r_IV_UncerProjExtreme(:, iExp) + r_Re_UncerExtreme(:, iExp)
    print(sprintf("%4.2f", r_EC_UncerExtreme(:, iExp) - r_EC_CoefUncerExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_IV_UncerProjExtreme(:, iExp)) + ", " + \
          sprintf("%4.2f", r_Re_UncerExtreme(:, iExp)) + " = " + r_Temp)
    print("Con Coef Var: " + sprintf("%4.2f", r_EC_CoefUncerExtreme(:, iExp)))
    print("")
end do
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
    sres@tmXBLabelFontHeightF             = 0.011
    sres@tmYLLabelFontHeightF             = 0.012

    txres                                 = True
    txres@txFontHeightF                   = 0.01
    
    sres@vpXF                             = 0.1 ; 0.3~0.95
    sres@vpYF                             = 0.8
    sres@vpWidthF                         = 0.3
    sres@vpHeightF                        = 0.25
    ;
    plot1                                 = new(3, graphic)
    plot1_EC                              = plot1
    plot1_Re                              = plot1
    plot1_ProjIV                          = plot1
    plot1_ConCoef                         = plot1
    xx_res2                               = (/0, 0.4, 0.8, 1.2/)
    dxx_res2                              = 2.
    s_Term                                = (/ "Near", "Middle", "Long" /)
    do iTerms = 0, 2
        sres@tmXBValues                   = xx_res2
        x                                 = sres@tmXBValues
        sres@gsnLeftString                = ""
        ;
        ;; Proj IV
        r_Temp                           := r_IV_UncerMean(iTerms, :)
        r_Temp                           := r_IV_UncerMean(iTerms, :) + \
                                            r_IV_UncerProjMean(iTerms, :) + \
                                            r_EC_UncerMean(iTerms, :) + \
                                            r_Re_UncerMean(iTerms, :) 
        sres@gsnXYBarChartColors          = "gold3" ; "blue"
        sres@tmXBLabels                   = s_VarNameTitle
        sres@tmXBOn                       = True
        plot1_Re(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; PD IV 
        r_Temp                            = r_IV_UncerMean(iTerms, :) + \
                                            r_EC_UncerMean(iTerms, :) + \
                                            r_Re_UncerMean(iTerms, :) 
        sres@gsnXYBarChartColors          = "orange" ; "green3"
        plot1_EC(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; Other Uncer 
        r_Temp                            = r_EC_UncerMean(iTerms, :) + \
                                            r_Re_UncerMean(iTerms, :) 
        sres@gsnXYBarChartColors          = "blue"
        plot1_ProjIV(iTerms)              = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; Con Coef Uncer 
        sres@gsnXYBarChartColors         = "gray"
        r_TempEC_Coef                   := r_EC_CoefUncerMean(iTerms, :)
        r_TempEC_Coef                    = where(r_EC_UncerMean(iTerms,:).lt.0,\
                                                 0, r_TempEC_Coef)
        plot1_ConCoef(iTerms)            = gsn_csm_xy(wks, x, r_TempEC_Coef, \
                                                      sres)
        ;; Reduced Uncer by EC 
        sres@gsnXYBarChartColors          = "green3"
        if (iTerms .eq. 0) then
            sres@gsnLeftString            = "(e) Annual Mean Warming"
        end if
        ;
        plot1(iTerms)                     = gsn_csm_xy(wks, x, \
                                            r_EC_UncerMean(iTerms, :), sres)
        ;
        r_xText                           = (/0.15, 0.25, 0.36/)
        r_yText                           = 0.46
        gsn_text_ndc(wks, s_Term(iTerms), r_xText(iTerms), r_yText, txres)
        ;
        overlay(plot1_Re(iTerms), plot1_Re(iTerms))
        overlay(plot1_Re(iTerms), plot1_EC(iTerms))
        overlay(plot1_Re(iTerms), plot1_ProjIV(iTerms))
        overlay(plot1_Re(iTerms), plot1(iTerms))
        overlay(plot1_Re(iTerms), plot1_ConCoef(iTerms))
        draw(plot1_Re(iTerms))
        sres@tmYLOn                       = False
        sres@tiYAxisOn                    = False
        ;
        xx_res2                           = xx_res2 + dxx_res2 ;2.5 
    end do
    ;
    sres@tmYLOn                           = True
    sres@vpXF                             = sres@vpXF + sres@vpWidthF + 0.08 ; 0.3~0.95
    plot1                                 = new(3, graphic)
    plot1_EC                              = plot1
    plot1_Re                              = plot1
    plot1_ProjIV                          = plot1
    plot1_ConCoef                         = plot1
    xx_res2                               = (/0, 0.4, 0.8, 1.2/)
    dxx_res2                              = 2.
    do iTerms = 0, 2
        sres@tmXBValues                   = xx_res2
        x                                 = sres@tmXBValues
        sres@gsnLeftString                = ""
        ;
        ;; Proj IV ; Other Uncer 
        r_Temp                           := r_IV_UncerExtreme(iTerms, :)
        r_Temp                           := r_IV_UncerExtreme(iTerms, :) + \
                                            r_IV_UncerProjExtreme(iTerms, :) + \
                                            r_EC_UncerExtreme(iTerms, :) + \
                                            r_Re_UncerExtreme(iTerms, :) 
        sres@gsnXYBarChartColors          = "gold3" ; "blue"
        sres@tmXBLabels                   = s_VarNameTitle
        sres@tmXBOn                       = True
        plot1_Re(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; PD IV; EC Uncer
        r_Temp                            = r_IV_UncerExtreme(iTerms, :) + \
                                            r_EC_UncerExtreme(iTerms, :) + \
                                            r_Re_UncerExtreme(iTerms, :) 
        sres@gsnXYBarChartColors          = "orange" ; "green3"
        plot1_EC(iTerms)                  = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; Other Uncer 
        r_Temp                            = r_EC_UncerExtreme(iTerms, :) + \
                                            r_Re_UncerExtreme(iTerms, :) 
        ; print(r_Temp + "")
        ; print("")
        sres@gsnXYBarChartColors          = "blue" ; "gold3" ; "yellow3"
        plot1_ProjIV(iTerms)              = gsn_csm_xy(wks, x, r_Temp, sres)
        ;
        ;; Con Coef Uncer 
        sres@gsnXYBarChartColors          = "gray"
        r_TempEC_Coef                   := r_EC_CoefUncerExtreme(iTerms, :)
        r_TempEC_Coef                    = where(r_EC_UncerExtreme(iTerms,:).lt.0,\
                                                 0, r_TempEC_Coef)
        plot1_ConCoef(iTerms)             = gsn_csm_xy(wks, x, r_TempEC_Coef, \
                                                       sres)
        ;; Reduced Uncer by EC 
        sres@gsnXYBarChartColors          = "green3" ; "orange"
        if (iTerms .eq. 0) then
            sres@gsnLeftString            = "(f) Extreme Warming"
            ;
        end if
        plot1(iTerms)                     = gsn_csm_xy(wks, x, \
                                            r_EC_UncerExtreme(iTerms, :), sres)
        ; print(r_IV_UncerExtreme(iTerms, :) + "")
        ; print("")
        ;
        r_xText                           = (/0.5, 0.6, 0.71/) + 0.03
        r_yText                           = 0.46
        gsn_text_ndc(wks, s_Term(iTerms), r_xText(iTerms), r_yText, txres)
        ;
        overlay(plot1_Re(iTerms), plot1_Re(iTerms))
        overlay(plot1_Re(iTerms), plot1_EC(iTerms))
        overlay(plot1_Re(iTerms), plot1_ProjIV(iTerms))
        overlay(plot1_Re(iTerms), plot1(iTerms))
        overlay(plot1_Re(iTerms), plot1_ConCoef(iTerms))
        draw(plot1_Re(iTerms))
        sres@tmYLOn                       = False
        sres@tiYAxisOn                    = False
        ;
        xx_res2                           = xx_res2 + dxx_res2 ;2.5 
    end do
    ;
    txres@txJust                          = "CenterLeft"
    txres@txFontColor                     = "blue"
    s_TermTemp                            = "Var(~F8~e~NF21~)" ; "Var Unrelated to EC"
    gsn_text_ndc(wks, s_TermTemp, 0.79, 0.7, txres)
    ;
    txres@txFontColor                     = "green3"
    s_TermTemp                            = "Var(EC)" ; "Reduced Var by EC" ; "Var Explained by Constraint Coef"
    gsn_text_ndc(wks, s_TermTemp, 0.79, 0.64, txres)
    ;
    txres@txFontColor                     = "gold3" ; "yellow3"
    s_TermTemp                            = "Var(IV~F10B~Proj~NF21~)" ; "Internal Variability in Projection"
    gsn_text_ndc(wks, s_TermTemp, 0.79, 0.79, txres)
    ;
    txres@txFontColor                     = "orange"
    s_TermTemp                            = "Var(IV~F10B~PD~NF21~)" ; "Internal Variability in Historical"
    gsn_text_ndc(wks, s_TermTemp, 0.79, 0.77, txres)
    ;
    txres@txFontColor                     = "gray"
    s_TermTemp                            = "Var(ConCoef)" ; "Var Explained by Constraint Coef"
    gsn_text_ndc(wks, s_TermTemp, 0.79, 0.55, txres)
    ;
    ; drawNDCGrid(wks)
    frame(wks)
print(s_FigOut + "")    
;
end








