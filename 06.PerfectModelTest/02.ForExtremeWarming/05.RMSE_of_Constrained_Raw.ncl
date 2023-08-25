; ==================================================================
; Author: Ziming Chen
; This code is for:
;   - cal the MME of RMSE btw constrained projection and pseudo obs, 
; and btw raw projection and pseudo obs 
;   - vis the results 
;   1. indicate the model, and then read the RMSE 
;   2. cal the MME and inter-model spread 
;   3. vis 
; ==================================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_DrawFunction.ncl"
;
undef("my_Legend")
function my_Legend(wks, plotUn, x[*]:numeric, y[*][*]:numeric)
local r_x, r_y
begin
r_x                        = x 
r_y                        = y 
;
  sspFill                  = (/(/ 252, 209, 197 /), (/ 146, 197, 222 /)/) /255.
  
  pgres1                   = True
  plres1                   = True
  plres1@gsLineDashPattern = 0
  plres1@gsLineThicknessF  = 4
  
  res_text                 = True
  res_text@txFontHeightF   = 0.02              ; change text size 
  res_text@txJust          = "CenterLeft"      ; text justification
  
  r_xx1                    = r_x
  r_yy1                    = r_y
  dum10                   := new(2, graphic)
  ; dum11                   := new(2, graphic)
  dum12                   := new(2, graphic)
  s_Text                   = (/ "Raw Projection", "Constrained Projection" /)
  do issp = 0, dimsizes(s_Text)-1
      pgres1@gsFillColor     := sspFill(issp, :) ;(/153, 0, 2/)/255.
      ; pgres1@gsFillOpacityF   = 0.4
      pgx                    := (/ r_xx1(0), r_xx1(0), r_xx1(1), \
                                        r_xx1(1), r_xx1(0)/)
      pgy                    := (/ r_yy1(issp,0), r_yy1(issp,1), r_yy1(issp,1), \
                                          r_yy1(issp,0), r_yy1(issp,0)  /)
      dum10(issp)              = gsn_add_polygon(wks, plotUn, pgx, pgy, pgres1)
      ; plres1@gsLineColor      := sspColor(issp, :)
      ; dum11(issp)              = gsn_add_polyline(wks, plotUn, \
      ;         (/r_xx1(0), r_xx1(1)/), (/avg(r_yy1(issp,:)), avg(r_yy1(issp,:))/), \
      ;         plres1)
      ;
      dum12(issp)              = gsn_add_text(wks, plotUn, s_Text(issp), \
                              r_xx1(1) + 0.1, avg(r_yy1(issp, :)), res_text)
  end do
  ;
  return([/dum10, dum12/])

;
end
;
begin
;
; ==================================================================
;   1. indicate the model, and then read the RMSE 
s_ModelSelected                 = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models                       = "_18Models"
s_4ModelsGSAT_Con               = "_23Models" ; Keep it as a largest ensemble
s_Exp                           = (/"ssp126", "ssp245", "ssp370", "ssp585"/)

s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Fig/06.PerfectModelTest/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FigOut                        = s_DirOut + "05.RMSE_of_Constrained_Raw" + s_4Models
wks                             = gsn_open_wks("X11", s_FigOut)
;
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
do iExp = 0, dimsizes(s_Exp)-1
    s_FileRead                  = s_DirRead + "CRPS_ForPMT_AllPseModel_" + \
                                  s_Exp(iExp) + s_4Models + ".nc"
    FileCRPS                    = addfile(s_FileRead, "r")
    do iModel = 0, dimsizes(s_ModelSelected)-1
        s_FileRead              = s_DirRead + \
                            "FurtherCorrectConstrainedWarming_PseModel_" + \
                            s_ModelSelected(iModel) + "_" + s_Exp(iExp) + \
                            s_4Models + ".nc"
        ls                     := systemfunc("ls " + s_FileRead)
        if (ismissing(ls(0))) then
            print(s_Exp(iExp) + ", " + s_ModelSelected(iModel) + " Missing!")
            continue
        end if
        FileRead               := addfile(s_FileRead, "r")
        s_VarName2             := "RMSE_China_ConstrainedWarming_" + s_Exp(iExp)
        r_RMSE                 := FileRead->$s_VarName2$
        s_VarName2             := "RMSE_China_RawProjectedWarming_"+ s_Exp(iExp)
        r_RMSE_Raw             := FileRead->$s_VarName2$
        ;
        s_VarName2              = "ConCRPS_ForPseudo_" + s_ModelSelected(iModel)\
                                + "_" + s_Exp(iExp)
        r_CRPS_Con             := FileCRPS->$s_VarName2$
        ;
        s_VarName2              = "RawCRPS_ForPseudo_" + s_ModelSelected(iModel)\
                                + "_" + s_Exp(iExp)
        r_CRPS_Raw             := FileCRPS->$s_VarName2$
        ;
        if (isdefined("r_RMSE_All") .eq. False) then
            i_dims             := array_append_record((/dimsizes(s_Exp), \
                                                dimsizes(s_ModelSelected)/), \
                                                dimsizes(r_RMSE), 0)
            r_RMSE_All          = new(i_dims, float)
            copy_VarCoords(r_RMSE, r_RMSE_All(0, 0, :))
            r_RMSE_All!0        = "Exp"
            r_RMSE_All&Exp      = s_Exp
            r_RMSE_All!1        = "Model"
            r_RMSE_All&Model    = s_ModelSelected
            ;
            r_RMSE_Raw_All      = r_RMSE_All
            ;
            r_CRPS_Con_All      = r_RMSE_All
            r_CRPS_Raw_All      = r_RMSE_Raw_All
        end if
        r_RMSE_All(iExp, iModel, :) = r_RMSE
        r_RMSE_Raw_All(iExp, iModel, :) = r_RMSE_Raw
        ;
        r_CRPS_Con_All(iExp, iModel, :) = r_CRPS_Con
        r_CRPS_Raw_All(iExp, iModel, :) = r_CRPS_Raw
    end do
end do
;
; ==================================================================
;   2. cal the MME and inter-model spread 
r_RMSE_AllMME                   = dim_avg_n_Wrap(r_RMSE_All, 1)
r_RMSE_AllSTD                   = dim_stddev_n_Wrap(r_RMSE_All, 1)
;
r_RMSE_Raw_AllMME               = dim_avg_n_Wrap(r_RMSE_Raw_All, 1)
r_RMSE_Raw_AllSTD               = dim_stddev_n_Wrap(r_RMSE_Raw_All, 1)
;
r_RMSE_Vis                      = r_RMSE_All(:, :2, :)
r_RMSE_Vis(:, 0, :)             = r_RMSE_AllMME
r_RMSE_Vis(:, 1, :)             = r_RMSE_AllMME - r_RMSE_AllSTD
r_RMSE_Vis(:, 2, :)             = r_RMSE_AllMME + r_RMSE_AllSTD
;
r_RMSE_Raw_Vis                  = r_RMSE_Vis
r_RMSE_Raw_Vis                  = r_RMSE_Raw_Vis@_FillValue
r_RMSE_Raw_Vis(:, 0, :)         = r_RMSE_Raw_AllMME
r_RMSE_Raw_Vis(:, 1, :)         = r_RMSE_Raw_AllMME - r_RMSE_Raw_AllSTD
r_RMSE_Raw_Vis(:, 2, :)         = r_RMSE_Raw_AllMME + r_RMSE_Raw_AllSTD
;
;; CRPS
r_CRPS_Con_AllMME               = dim_avg_n_Wrap(r_CRPS_Con_All, 1)
r_CRPS_Con_AllSTD               = dim_stddev_n_Wrap(r_CRPS_Con_All, 1)
;
r_CRPS_Raw_AllMME               = dim_avg_n_Wrap(r_CRPS_Raw_All, 1)
r_CRPS_Raw_AllSTD               = dim_stddev_n_Wrap(r_CRPS_Raw_All, 1)
;
r_CRPS_Vis                     := r_CRPS_Con_All(:, :2, :)
r_CRPS_Vis(:, 0, :)             = r_CRPS_Con_AllMME
r_CRPS_Vis(:, 1, :)             = r_CRPS_Con_AllMME - r_CRPS_Con_AllSTD
r_CRPS_Vis(:, 2, :)             = r_CRPS_Con_AllMME + r_CRPS_Con_AllSTD
;
r_CRPS_Raw_Vis                 := r_CRPS_Raw_All(:, :2, :)
r_CRPS_Raw_Vis(:, 0, :)         = r_CRPS_Raw_AllMME
r_CRPS_Raw_Vis(:, 1, :)         = r_CRPS_Raw_AllMME - r_CRPS_Raw_AllSTD
r_CRPS_Raw_Vis(:, 2, :)         = r_CRPS_Raw_AllMME + r_CRPS_Raw_AllSTD
;
;; Relative changes of RMSE 
r_RMSE_Relative_Vis             = r_RMSE_Vis
r_RMSE_Relative_Vis             = r_RMSE_Vis - r_RMSE_Raw_Vis
; r_RMSE_Relative_Vis             = r_RMSE_Relative_Vis / r_RMSE_Raw_Vis * 100.
r_Temp                         := conform(r_RMSE_Relative_Vis, \
                                          r_RMSE_Raw_Vis(:, 0, :), (/0, 2/))
r_RMSE_Relative_Vis             = r_RMSE_Relative_Vis / r_Temp * 100.
;
;; Relative changes of CRPS 
r_CRPS_Relative_Vis             = r_CRPS_Vis
r_CRPS_Relative_Vis             = r_CRPS_Vis - r_CRPS_Raw_Vis
; r_CRPS_Relative_Vis             = r_CRPS_Relative_Vis / r_CRPS_Raw_Vis * 100.
r_Temp                         := conform(r_CRPS_Relative_Vis, \
                                          r_CRPS_Raw_Vis(:, 0, :), (/0, 2/))
r_CRPS_Relative_Vis             = r_CRPS_Relative_Vis / r_Temp * 100.
;
do iExp = 0, dimsizes(s_Exp)-1
    print(s_Exp(iExp) + "")
    print("Relative Changes of RMSE")
    print(sprintf("%4.2f", r_RMSE_Relative_Vis(iExp, 0, :)) + " (" + \
          sprintf("%4.2f", r_RMSE_Relative_Vis(iExp, 1, :)) + "~" + \
          sprintf("%4.2f", r_RMSE_Relative_Vis(iExp, 2, :)) + ")")
    ;
    print("Relative Changes of CRPS")
    print(sprintf("%4.2f", r_CRPS_Relative_Vis(iExp, 0, :)) + " (" + \
          sprintf("%4.2f", r_CRPS_Relative_Vis(iExp, 1, :)) + "~" + \
          sprintf("%4.2f", r_CRPS_Relative_Vis(iExp, 2, :)) + ")")
    print("")
    ; print("RMSE")
    ; print("Constrained")
    ; print(sprintf("%4.2f", r_RMSE_Vis(iExp, 0, :))  + " (" + \
    ;       sprintf("%4.2f", r_RMSE_Vis(iExp, 1, :))  + "~" + \
    ;       sprintf("%4.2f", r_RMSE_Vis(iExp, 2, :))  + ")")
    ; print("Raw")
    ; print(sprintf("%4.2f", r_RMSE_Raw_Vis(iExp, 0, :))  + " (" + \
    ;       sprintf("%4.2f", r_RMSE_Raw_Vis(iExp, 1, :))  + "~"+\
    ;       sprintf("%4.2f", r_RMSE_Raw_Vis(iExp, 2, :))  + ")")
    ; print("")
    ; print("CRPS")
    ; print("Constrained")
    ; print(sprintf("%4.2f", r_CRPS_Vis(iExp, 0, :))  + " (" + \
    ;       sprintf("%4.2f", r_CRPS_Vis(iExp, 1, :))  + "~" + \
    ;       sprintf("%4.2f", r_CRPS_Vis(iExp, 2, :))  + ")")
    ; print("Raw")
    ; print(sprintf("%4.2f", r_CRPS_Raw_Vis(iExp, 0, :))  + " (" + \
    ;       sprintf("%4.2f", r_CRPS_Raw_Vis(iExp, 1, :))  + "~"+\
    ;       sprintf("%4.2f", r_CRPS_Raw_Vis(iExp, 2, :))  + ")")
    ; print("")
end do
;
; ===========================================================================
;   3. vis 
sres                                  = True
    
sres@trXMinF                          = 0.3
sres@trXMaxF                          = 5.2
r_trYF_L  = (/ (/0.15, 1.3/), (/0.05, 1.5/), (/0.35, 2.5/) /)
r_trYF_R  = (/ (/0.25, 5.95/),   (/0, 7.5/),    (/0.5, 10/) /)
; sres@trYMaxF                          = 1.6 ;max(r_G11) + 0.1
; sres@trYMinF                          = -0.9 ;min(r_G11) - 0.1
    
sres@gsnFrame                         = False
sres@gsnDraw                          = False
sres@gsnXYBarChart                    = True     ;; 设定柱状图
sres@gsnXYBarChartBarWidth            = 0.35
    
s_LeftString                          = (/(/"(a)", "(b)", "(c)"/) + " RMSE", \
                                          (/"(d)", "(e)", "(f)"/) + " CRPS"/)
; sres@gsnLeftString                    = "Global Land Summer Monsoon"
sres@gsnLeftStringFontHeightF         = 0.03
sres@gsnRightStringFontHeightF        = 0.025
s_RightString                         = (/"2021 to 2040", "2041 to 2060", \
                                          "2080 to 2099"/)
; sres@gsnRightString                   = i_yrRead(0) + "-" + i_yrRead(1)
sres@gsnCenterString                  = ""
; sres@gsnLeftStringOrthogonalPosF      = 0.02
; sres@gsnRightStringFontHeightF        = 0.01
; sres@gsnLeftStringFontHeightF         = 0.015
; sres@gsnYRefLine                      = 0.
    
s_tiYAxis                             = (/"RMSE (K)", "CRPS"/)
; sres@tiYAxisString                    = "RMSE (K)"
sres@tiYAxisFontHeightF               = 0.025
sres@tmXTOn                           = False
sres@tmYROn                           = False
sres@tmXBMode                         = "Explicit"
sres@tmXBValues                       = fspan(1., 4.5, 4)
s_VarNameTitle                        = (/"SSP1-2.6", "SSP2-4.5", "SSP3-7.0", \
                                          "SSP5-8.5"/)
sres@tmXBLabels                       = s_VarNameTitle
; sres@tmXBLabelAngleF                  = 50.
sres@tmXBLabelFontHeightF             = 0.025
sres@tmYLLabelFontHeightF             = 0.025

sres@vpWidthF                         = 0.9
;
; printVarSummary(r_RMSE_Raw_Vis) ; [Exp | 4] x [Model | 3] x [Terms | 3]
; printVarSummary(r_CRPS_Raw_Vis) ; [Exp | 4] x [Model | 3] x [Terms | 3]
x                                     = sres@tmXBValues
plot1                                 = new(6, graphic)
plot2                                 = plot1
plot_error                            = new((/6, 4, 3/), graphic)
plot_error2                           = plot_error
dum0_1                                = new((/3, 4/), graphic)
dum1_1                                = dum0_1
iplot                                 = 0
do iTerms = 0, 2
    ;; RMSE 
    sres@trYMaxF                      = r_trYF_L(iTerms, 1)
    sres@trYMinF                      = r_trYF_L(iTerms, 0)
    sres@gsnRightString               = s_RightString(iTerms)
    sres@gsnLeftString                = s_LeftString(0, iTerms)
    sres@tiYAxisString                = s_tiYAxis(0)
    sres@gsnXYBarChartColors          = (/ 252, 209, 197 /) /255.
    plot1(iplot)                      = gsn_csm_xy(wks, x-0.22, \
                                        r_RMSE_Raw_Vis(:, 0, iTerms), sres)
    plot_error(iplot, :, :)           = add_yerrorbar_Function(wks, \
                          plot1(iplot), x-0.22, r_RMSE_Raw_Vis(:, 0, iTerms), \
                          r_RMSE_Raw_Vis(:, 1:, iTerms), False)

    sres@gsnRightString               = ""
    sres@gsnLeftString                = ""
    sres@gsnXYBarChartColors          = (/ 146, 197, 222 /) /255.
    plot2(iplot)                      = gsn_csm_xy(wks, x+0.22, \
                                                r_RMSE_Vis(:, 0, iTerms), sres)
    plot_error2(iplot, :, :)          = add_yerrorbar_Function(wks, \
                          plot2(iplot), x+0.22, r_RMSE_Vis(:, 0, iTerms), \
                          r_RMSE_Vis(:, 1:, iTerms), False)
    overlay(plot1(iplot), plot2(iplot))
    ;; legend 
    r_Temp                           := (/x(0)-0.4, x(0)/)
    r_yy                              = (/ (/(/1.12, 1.2/),  (/1,    1.08/)/), \
                                           (/(/1.34, 1.44/), (/1.2,  1.3/)/), \
                                           (/(/2.24, 2.4/),  (/2,    2.16/)/) /)
    l_list                           := my_Legend(wks, plot1(iplot), r_Temp, \
                                                  r_yy(iTerms, :, :) )
    dum0_1(iTerms, :1)                = l_list[0]
    dum0_1(iTerms, 2:)                = l_list[1]
    ;; legend 
    ;
    iplot                             = iplot + 1
    ;
    ;; CRPS 
    sres@trYMaxF                      = r_trYF_R(iTerms, 1)
    sres@trYMinF                      = r_trYF_R(iTerms, 0)
    sres@gsnRightString               = s_RightString(iTerms)
    sres@gsnLeftString                = s_LeftString(1, iTerms)
    sres@tiYAxisString                = s_tiYAxis(1)
    sres@gsnXYBarChartColors          = (/ 252, 209, 197 /) /255.
    plot1(iplot)                      = gsn_csm_xy(wks, x-0.22, \
                                        r_CRPS_Raw_Vis(:, 0, iTerms), sres)
    plot_error(iplot, :, :)           = add_yerrorbar_Function(wks, \
                      plot1(iplot), x-0.22, r_CRPS_Raw_Vis(:, 0, iTerms), \
                      r_CRPS_Raw_Vis(:, 1:, iTerms), False)
    ;
    sres@gsnXYBarChartColors          = (/ 146, 197, 222 /) /255.
    sres@gsnRightString               = ""
    sres@gsnLeftString                = ""
    plot2(iplot)                      = gsn_csm_xy(wks, x+0.22, \
                                                r_CRPS_Vis(:, 0, iTerms), sres)
    plot_error2(iplot, :, :)          = add_yerrorbar_Function(wks, \
                              plot2(iplot), x+0.22, r_CRPS_Vis(:, 0, iTerms), \
                              r_CRPS_Vis(:, 1:, iTerms), False)
    overlay(plot1(iplot), plot2(iplot))
    ;; legend 
    r_Temp                           := (/x(0)-0.4, x(0)/)
    r_yy                              = (/ (/(/5.45, 5.85/), (/4.9, 5.3/)/), \
                                           (/(/6.65,7.2/), (/5.9, 6.45/)/), \
                                           (/(/8.9, 9.5/), (/7.9, 8.5/)/) /)
    l_list                           := my_Legend(wks, plot1(iplot), r_Temp, \
                                                  r_yy(iTerms, :, :) )
    dum1_1(iTerms, :1)                = l_list[0]
    dum1_1(iTerms, 2:)                = l_list[1]
    ;; legend 
    iplot                             = iplot + 1
end do
;
resP                                  = True
resP@gsnMaximize                      = False
gsn_panel(wks, plot1, (/3, 2/), resP)

print(s_FigOut + "")
;
end







