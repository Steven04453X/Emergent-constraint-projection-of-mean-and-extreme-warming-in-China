; =================================================================
; 本程序计算当前气候态的升温趋势，与未来指定时段内的中国地区极端升温幅度之间的简单相关
;   1. 指定模式，变量名，预估的情景，当前气候态的时段和不同预估时期的时段，当前求趋势的时段
;   2. 首先读取趋势数据，然后读取预估极端升温的数据，随后计算中国区域的平均
;   3. 求两者的简单相关，并绘图
;; 注意本程序附上了观测中GSAT的升温趋势，以及对应的未来升温幅度，作为约束的基础
; 本程序考虑观测资料之间的差异，以及内部变率的贡献
; 注意本程序是针对多样本的结果，即先对每一个模式的所有样本进行平均，然后再进行约束
; =================================================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
; ===================================================================
; 本程序用于绘制模式间的散点图，不同模式用不同的形状或颜色的标记绘制，并求两者的相关和线性拟合线
; 只需要输入两序列即可
undef("my_plot_Model_Scatter_Special_withConfidenceLines")
function my_plot_Model_Scatter_Special_withConfidenceLines(X[*]:numeric, \
                                                        Y[*]:numeric, wks, res)
;; 注意r_X和r_Y的坐标上，必须是对应的模式
;; 注意这个程序进行的部分修改，主要是给出观测中X的量值，如果没有指定Y，就通过拟合线计算，作为约束
local r_Coef, r_Reg
begin
r_X                     = X 
r_Y                     = Y

r_X!0                   = "Models"
r_Temp                  = avg(r_X)
r_Temp!0                = "Models"
r_Temp&Models           = "CMIP6 MME"
r_X                    := array_append_record(r_Temp, r_X, 0)

r_Y!0                   = "Models"
r_Temp                  = avg(r_Y)
r_Temp!0                = "Models"
r_Temp&Models           = "CMIP6 MME"
r_Y                    := array_append_record(r_Temp, r_Y, 0)

;; Obs Uncertainty
r_ObsX_STD              = res@ObsX_STD
r_ModelX_STD            = stddev(r_X(1:))
r_SNR                   = r_ModelX_STD^2 / r_ObsX_STD^2
print("SNR: " + r_SNR)
r_Correction            = 1. / r_SNR
r_Correction            = 1. / (1. + r_Correction)

s_ModelCMIP             = r_X&Models
; ---------------------------------------------------------
;   4. 计算两者的简单相关，然后绘图
r_Coef                  = escorc(r_X, r_Y)
; r_Coef                  = r_Coef * r_Correction
; print(r_Coef + "")
;
l_Model_X_Y             = [/ s_ModelCMIP, r_X, r_Y/]

;; 计算线形回归
r_Reg                  := regline_stats(r_X, r_Y)
r_Reg_Corr              = r_Reg * r_Correction

;;;; 计算相关系数的显著性水平
r_tval_Corr             = r_Reg@tval 
r_tval_Corr            := r_tval_Corr(1)
df                      = r_Reg@nptxy - 2
r_Sig                   = betainc(df/(df+r_tval_Corr^2), df/2.0, 0.5)
;;;; 计算相关系数的显著性水平
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; confidence intervals
r_b                     = r_Reg@b
r_yInter                = r_b(0)
r_Slope                 = r_Reg_Corr ; r_b(1)
r_yFit                  = r_yInter + r_Slope * r_X
r_x_Mean                = abs(avg(r_X))
r_yFit_x                = (/min(r_X)-r_x_Mean, max(r_X)+r_x_Mean/)
r_yFit_Draw             = r_yInter + r_Slope * r_yFit_x
;
;;;; 95%范围: yFit + or - tval(alpha/2) * std * sqrt(1/n + (xp-x_mean)^2/SSxx)
;;;; SSxx = sum(xp - x_mean)^2
;;;; std  = sqrt(SSE / (n-2)), SSE = SSyy - Beta1 * SSxy
;;;; Attention: this is one-tailed t value 
r_tval                  = 2.052     ; p = 0.05, n = 30
r_tval                  = 1.703     ; p = 0.1, n = 30 

r_SSxx                  = avg(r_X)
r_SSyy                  = avg(r_Y)
; print(r_SSxx + "  " + r_SSyy)
r_SSxy                  = sum((r_X - r_SSxx) * \
                              (r_Y - r_SSyy))
r_SSxx                 := sum((r_X - r_SSxx)^2)
r_SSyy                 := sum((r_Y - r_SSyy)^2)
; print(r_SSxx + "  " + r_SSyy + "  " + r_SSxy)
; print("")
r_SSE                   = sum( (r_yFit - r_Y)^2)
; print(r_SSE + "")
; print("")
r_Temp                 := r_SSE / (dimsizes(r_Y) - 2.)
r_std                   = sqrt(r_Temp)
; print(r_std + "")

n                       = dimsizes(r_Y)
r_X2                    = fspan(min(r_X)-(r_x_Mean*2), \
                                max(r_X)+(r_x_Mean*2), 40)
r_yFit2                 = r_yInter + r_Slope * r_X2
r_Yrange                = r_tval * r_std 
r_Yrange               := r_Yrange * sqrt(1.0/n + (r_X2 - r_x_Mean)^2 / r_SSxx)
r_Yrange               := (/ r_yFit2 - r_Yrange, r_yFit2 + r_Yrange /)
; printVarSummary(r_Yrange)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; prediction interval 
r_Yrange1              := r_Reg@YPI025
r_Yrange2              := r_Reg@YPI975
r_Yrange               := (/r_Yrange1, r_Yrange2/)
r_X2                   := r_X
;
if (any(ismissing(r_X))) then
    i_Ind              := ind(ismissing(r_X) .eq. False .and. \
                              ismissing(r_Y) .eq. False)
    r_X2               := r_X(i_Ind)
    r_Y2               := r_Y(i_Ind)
    ; print(r_X2 + ", " + r_Y2)
    r_Reg0             := regline_stats(r_X2, r_Y2)
    r_Yrange1          := r_Reg0@YPI025
    r_Yrange2          := r_Reg0@YPI975
    r_Yrange           := (/r_Yrange1, r_Yrange2/)
end if
i_pqsort                = dim_pqsort(r_X2, 2)
r_Yrange                = r_Yrange(:, i_pqsort)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; ===============================================
;   4.简单相关，然后进行绘图
l_Fit        = [/ r_yFit_x, r_X2, r_yFit_Draw, r_Yrange /]
;
if (res .eq. True) then
    res2     = res
else 
    res2     = True
end if

;;; just my_ObsX without my_ObsY, and then estimate the my_ObsY using fit line
; r_Reg, r_Reg@yintercept
if (ismissing(res2@my_ObsX) .eq. False .and. \
    ismissing(res2@my_ObsY) .eq. True) then 
    r_yintercept                = r_Reg@yintercept
    r_ObsX                      = tofloat(res2@my_ObsX)
    r_ObsY                      = r_ObsX * r_Reg_Corr + r_yintercept
    res2@my_ObsY                = r_ObsY
    print("Constrained Y: " + r_ObsY)
    print("Obs X: " + r_ObsX)
end if 
;
;;; Correct each model 
r_Y_Obs_EachModel               = r_Y(1:)
do iModel = 1, dimsizes(r_X)-1
    dx                          = r_X(iModel) - r_ObsX
    dy                          = r_Reg_Corr * dx ;+ r_yintercept
    r_Y_Obs_EachModel(iModel-1) = r_Y(iModel) - dy
end do
print("Corr Para: " + r_Correction)
print("Models: Constrained, Raw")
print(r_X&Models(1:) + ": " + sprintf("%4.2f", r_Y_Obs_EachModel) + ", " + \
      sprintf("%4.2f", r_Y(1:)) )
print("STD: " + stddev(r_Y_Obs_EachModel) + ", " + stddev(r_Y(1:)))
print("")
; exit()

res2@gsnDraw                    = False
res2@gsnFrame                   = False
res2@xyMarkLineMode             = "Markers"
res2@xyMarkerSizeF              = 0.01
;
;; x & y axis range 
if (ismissing(res2@trXMinF) .eq. True) then
    res2@trXMinF                = 0 ; 21 + 273.15
end if
if (ismissing(res2@trXMaxF) .eq. True) then 
    res2@trXMaxF                = 3 ; 25 + 273.15
end if 
if (ismissing(res2@trYMinF) .eq. True) then 
    res2@trYMinF                = min(r_Yrange) ; 4.5
end if 
if (ismissing(res2@trYMaxF) .eq. True) then 
    res2@trYMaxF                = max(r_Yrange) ; 8
end if 

if (ismissing(res2@tiYAxisString) .eq. True) then 
    res2@tiYAxisString          = "TXx Warming (~F35~J~F21~C)" ;"CEP SST ~F26~(~F35~J~F26~C)"
end if 
if (ismissing(res2@tiXAxisString) .eq. True) then
    res2@tiXAxisString          = "Present-day GSAT Trend (~F35~J~F21~C 10yr~SF21~-1~NF21~)" ; "LH Flux (W m~SF21~-2~NF21~)"
end if
res2@tiYAxisFontHeightF         = 0.014
res2@tiXAxisFontHeightF         = 0.014
res2@tiMainFontHeightF          = 0.018

if (ismissing(res2@gsnCenterString) .eq. True) then 
    res2@gsnCenterString        = ""
end if 
if (ismissing(res2@gsnRightString) .eq. True) then 
    res2@gsnRightString         = "r = " + sprintf("%4.2f", r_Coef)
end if 
print("Sig Level: " + r_Sig + "")

if (r_Sig .le. 0.01) then
    res2@gsnRightString         = res2@gsnRightString + " (p < 0.01)"
else if (r_Sig .le. 0.05 .and. r_Sig .gt. 0.01) then
    res2@gsnRightString         = res2@gsnRightString + " (p < 0.05)"
end if
end if

res2@gsnStringFontHeightF       = 0.017
res2@gsnLeftStringFontHeightF   = 0.018
;
res2@tmXTOn                     = False
res2@tmYROn                     = False
    
res2@vpHeightF                  = 0.4
res2@vpWidthF                   = 0.4

if (ismissing(res2@tmXBLabelFontHeightF) .eq. True) then 
    res2@tmXBLabelFontHeightF   = 0.014
end if 
if (ismissing(res2@tmYLLabelFontHeightF) .eq. True) then
    res2@tmYLLabelFontHeightF   = 0.014
end if 

;;
;
r_Markers                       = (/16, 2, 3, 5, 4, 6, 7, 9 /)
r_Color = (/ (/0, 0, 0/), (/153, 0, 2/), (/196, 121, 0/), (/112, 160, 205/), \
           (/0, 52, 102/), (/84, 39, 143/) /) / 255. ; , 

res2@xyMarkers                  = r_Markers(0)
res2@xyMarkerColors             = r_Color(0, :)

res3                            = res2
res3@tmLabelAutoStride          = True
res3@gsMarkerSizeF              = 0.009
res3@gsMarkerThicknessF         = 2.0
;
res2@vpYF                       = 0.9
res2@vpXF                       = 0.2
;
s_Model                        := l_Model_X_Y[0]
r_DrawX                         = l_Model_X_Y[1]
r_DrawY                         = l_Model_X_Y[2]
plot                            = gsn_csm_xy(wks, r_DrawX*0+99, r_DrawY*0+99, res2)

plotMarker                     := new(dimsizes(s_Model), graphic)

r_Color_2                       = new((/dimsizes(s_Model), 3/), float)
r_Markers_2                     = new(dimsizes(s_Model), integer)

res3@gsMarkerIndex              = r_Markers(0)
res3@gsMarkerColor             := (/153, 0, 2/) / 255. ; "red"
i_Ind                          := ind(s_Model .eq. "CMIP6 MME")
plotMarker(i_Ind)               = gsn_add_polymarker(wks, plot, \
                                r_DrawX(i_Ind), r_DrawY(i_Ind), res3)
r_Markers_2(i_Ind)              = res3@gsMarkerIndex
r_Color_2(i_Ind, :)             = res3@gsMarkerColor


iModel                          = 0
do iMarker = 1, dimsizes(r_Markers)-1
    res3@gsMarkerIndex          = r_Markers(iMarker)
    do iColor = 0, dimsizes(r_Color(:, 0))-1
        res3@gsMarkerColor     := r_Color(iColor, :)
        iModel                  = iModel + 1
        if (iModel .ge. dimsizes(s_Model)) then
            break
        end if
        ;
        print(s_Model(iModel) + "")
        if (s_Model(iModel) .eq. "CMIP5 MME" .or. \
            s_Model(iModel) .eq. "CMIP6 MME") then
            continue
        end if
        plotMarker(iModel)      = gsn_add_polymarker(wks, plot, \
                                r_DrawX(iModel), r_DrawY(iModel), res3)
        r_Markers_2(iModel)     = res3@gsMarkerIndex
        r_Color_2(iModel, :)    = res3@gsMarkerColor
    end do
    if (iModel .ge. dimsizes(s_Model)) then
        break
    end if
end do
;
;; 补充线形拟合线  
res4                            = res2 ;True
res4@xyMarkLineMode             = "Lines"
res4@xyLineColor                = "black"
res4@xyLineThicknessF           = 3
delete([/ r_yFit2, res4@xyMarkerSizeF, res4@xyMarkLineMode, \
       res4@xyMarkers, res4@xyMarkerColors/])

res4@tiYAxisString              = ""
res4@tiXAxisString              = ""
res4@gsnLeftString              = ""
res4@gsnCenterString            = ""
res4@gsnRightString             = ""

r_Contrast2                     = l_Fit[1]
r_yFit2                        := l_Fit[2]
r_Temp                         := l_Fit[0]

res4@xyDashPattern              = 0
plot_Reg                        = gsn_csm_xy(wks, r_Temp, r_yFit2, res4)
overlay(plot, plot_Reg)

res4@xyLineThicknesses          = 1
res4@xyLineColor                = "black"
res4@xyDashPattern              = 1
plot_Reg2                       = gsn_csm_xy(wks, r_Contrast2, r_Yrange(0, :), \
                                            res4)
overlay(plot, plot_Reg2)
plot_Reg3                       = gsn_csm_xy(wks, r_Contrast2, r_Yrange(1, :), \
                                             res4)
overlay(plot, plot_Reg3)
delete([/res2@xyLineThicknesses/])

;
;; 标记观测的水平
if (ismissing(res2@my_ObsX) .eq. False .and. ismissing(res2@my_ObsY) .eq. False) then 
    r_ObsX                      = res2@my_ObsX         ; W m-2
    r_ObsY                      = res2@my_ObsY          ; mm day-1
    resLine                     = True
    resLine@gsLineDashPattern   = 1
    resLine@gsLineColor         = "gray40"
    resLine@PostDraw            = True
    resLine@gsLineThicknessF    = 4
    ;;; Vertical 
    plx                         = (/ r_ObsX, r_ObsX /)
    ply                         = (/res2@trYMinF, r_ObsY/)
    dum                         = new(3, graphic)
    dum(0)                      = gsn_add_polyline(wks, plot, plx, ply, resLine)
    ;
    ;; Shading the Obs STD 
    plx:= (/ r_ObsX-r_ObsX_STD, r_ObsX+r_ObsX_STD, r_ObsX+r_ObsX_STD, r_ObsX-r_ObsX_STD, r_ObsX-r_ObsX_STD /)
    ply:= (/res2@trYMinF, res2@trYMinF, res2@trYMaxF, res2@trYMaxF, res2@trYMinF/)
    gsres                       = True
    gsres@gsFillColor           = "gray"
    gsres@gsFillOpacityF        = 0.3
    gsres@tfPolyDrawOrder       = "Predraw"
    dum(1)                      = gsn_add_polygon(wks, plot, plx, ply, gsres)

    ;;; Horizontal
    plx                        := (/ res2@trXMinF, r_ObsX /)
    ply                        := (/ r_ObsY, r_ObsY /)
    dum(2)                      = gsn_add_polyline(wks, plot, plx, ply, resLine)
end if 
;
draw(plot)
;
;----------------------------------------------------------------------
; Draw some markers and text.
;----------------------------------------------------------------------
  labels2  = s_Model ; array_append_record("MME", s_Model, 0)
  markers1 = r_Markers_2 
  Color2   = r_Color_2
  ypos2    = fspan(res2@vpYF, res2@vpYF - res2@vpHeightF -0.05, dimsizes(labels2)) 
  ypos3    = fspan(res2@vpYF, res2@vpYF - res2@vpHeightF*3 -0.05, \
                   dimsizes(labels2))

  dum2     = new(dimsizes(s_Model), graphic)
  dum2Text = dum2

  mkres               = True         ; Marker resources
  txres               = True         ; Text resources
  txres@txFontHeightF = 0.01
  txres@txJust        = "CenterLeft"
  
  do i = 0, dimsizes(s_Model)-1
     mkres@gsMarkerThicknessF = 2.
     mkres@gsMarkerSizeF      = 0.01
     mkres@gsMarkerIndex     := markers1(i)
     mkres@gsMarkerColor     := Color2(i, :)

     gsn_polymarker_ndc(wks,           0.65,  ypos2(i), mkres)
     gsn_text_ndc      (wks,labels2(i),0.68, ypos2(i),txres)
     ;
     dum2(i)          = gsn_add_polymarker(wks, plot, 0.65, ypos3(i), mkres)
     dum2Text(i)      = gsn_add_text(wks, plot, labels2(i), 0.68, ypos3(i), txres)
  end do
;
; if (ismissing(res2@my_ObsX) .eq. False .and. ismissing(res2@my_ObsY) .eq. False) then 
;     mkres@gsLineThicknessF = 5
;     mkres@gsLineColor      = resLine@gsLineColor
;     mkres@gsLineDashPattern= resLine@gsLineDashPattern
;     yy                     = abs(ypos2(2) - ypos2(1))
;     yy                    := ypos2(dimsizes(s_ModelCMIP5)-1) - i_temp * yy
;     gsn_polyline_ndc(wks, (/0.64, 0.67/), (/yy, yy/), mkres)
;     s_yy                   = "Obs Constraint"
;     gsn_text_ndc(wks, s_yy, 0.68, yy, txres)
; end if 

; drawNDCGrid(wks)
return([/r_ObsX, r_ObsY, plot, plotMarker, plot_Reg, plot_Reg2, plot_Reg3, \
       dum, dum2, dum2Text, r_Y_Obs_EachModel, r_Y(1:)/])
end 
; ==================================================================================

begin
;
; =================================================================
;   1. 指定模式，变量名，预估的情景，当前气候态的时段和不同预估时期的时段，当前求趋势的时段
; s_ModelSelected                 = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
;                             "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_ModelSelected                 = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models                       = "_18Models" ; "_30Models" ; "" ; "_4Models" ; 
s_Exp                           = "ssp370" ; "ssp585" ; "ssp245" ; "ssp126" ; 
s_ExpTitle                      = "SSP3-7.0" ; "SSP5-8.5" ; "SSP2-4.5" ; "SSP1-2.6" ; 

if (s_Exp .eq. "ssp126") then
    s_ModelSelected            := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                       "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                       "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if

s_VarName                       = "tas_Monthly"
i_yrHist                        = (/1995, 2014/)
i_yrProj                        = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)
i_Month                         = ispan(1, 12, 1)
s_Month                         = "AnMean"

i_yrTrend                       = (/1981, 2014/)

s_DirOut  = "/global/cfs/cdirs/m1867/zmchen/Work/OtherWorks/2020/ProjectionInChina/Fig/05.GSAT_scaled_RegionalTXx/02.ConstrainAllReal/"
system("mkdir -p " + s_DirOut)
s_FigOut  = s_DirOut + "05.TXx_inChina_vs_CurrGSAT_Trend_" + \
            s_Month + "_" + s_Exp
wks       = gsn_open_wks("X11", s_FigOut)

s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead                      = s_DirRead + \
                                "GSAT_AndRegionalTrend_1981-2014_23Models.nc"
FileGSAT_PD                     = addfile(s_FileRead, "r")
;
;; Obs Warming Trend 
s_DirOut2 = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/TEMP/POST/YearMon/"
s_FileRead                      = s_DirOut2 + \
                "SurfaceAirTemp_Post_YearMon_1900-2014.nc"
FileObs                         = addfile(s_FileRead, "r")

; s_VarNameObs                    = (/ "Cowtan_WayV2", "GISTEMP" /)
s_VarNameObs                    = (/ "Cowtan_WayV2", "BEST", \
                                     "GISTEMP", "NOAAGlobalTemp" /)
r_TrendObs                      = new(dimsizes(s_VarNameObs), float)
do iObs = 0, dimsizes(s_VarNameObs)-1
    s_VarName2                  = "tas_" + s_VarNameObs(iObs)
    r_HadCRUT4                 := FileObs->$s_VarName2$(\
                                {i_yrTrend(0):i_yrTrend(1)}, :, :, :)
    r_HadCRUT4                 := dim_avg_n_Wrap(r_HadCRUT4, 1)
    r_HadCRUT4!1                = "lat"
    r_HadCRUT4!2                = "lon"
    r_HadCRUT4                 := my_AreaMean2([/r_HadCRUT4/], \
                                            r_HadCRUT4&lat, r_HadCRUT4&lon, 0)
    ;
    r_Temp_Trend                = trend_manken(r_HadCRUT4, False, 0)
    r_TrendObs(iObs)            = r_Temp_Trend(1)
end do 
r_TrendObs                      = r_TrendObs * 10.
print(r_TrendObs + "")
r_TrendSTD                      = stddev(r_TrendObs)
print("Obs STD: " + r_TrendSTD)
; exit()
r_TrendObs                     := avg(r_TrendObs)
r_TrendObs@units                = "K per decade"
;
;; Blend Effect of surface temperature
r_TrendObs                      = r_TrendObs + 0.013
r_TrendSTD                      = r_TrendSTD + 0.005
;
; =================================================================
;   2. 首先读取趋势数据，然后读取预估升温的数据
r_GSAT_Trend                    = new(dimsizes(s_ModelSelected), float)
r_ProjWarming                   = new((/dimsizes(s_ModelSelected), 3/), float)
r_ProjWarming!0                 = "Models"
r_ProjWarming&Models            = s_ModelSelected
r_ProjWarming!1                 = "Terms"
r_ProjWarming&Terms             = (/ "Near", "Mid", "Long" /)
copy_VarCoords(r_ProjWarming(:, 0), r_GSAT_Trend)
do iModel = 0, dimsizes(s_ModelSelected)-1
    l_list                     := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "Amon")
    s_Real                     := l_list[1]
    l_list_                    := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelected(iModel), s_Exp, "Amon")
    s_Real_                    := l_list_[1]
    if (all(ismissing(s_Real)) .or. dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
        s_Real                 := s_Real_
    end if
    ;
    if (all(ismissing(s_Real)) ) then
        print(s_ModelSelected(iModel) + " Missing!")
        continue
    end if
    ;
    s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
    s_FileRead                  = s_DirRead + "TXx_AnMean_" + s_Exp + "_" + \
                                  s_ModelSelected(iModel) + "_CMIP6_AllReal.nc"
    FileGSAT_Proj              := addfile(s_FileRead, "r")
    s_FileRead                  = s_DirRead + "TXx_AnMean_historical_" + \
                                  s_ModelSelected(iModel) + "_CMIP6_AllReal.nc"
    FileGSAT_Hist              := addfile(s_FileRead, "r")
    ;
    r_GSAT_TrendTemp           := new(dimsizes(s_Real), float)
    r_ProjWarmingTemp          := new((/dimsizes(s_Real), 3/), float)
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2              = "GSAT_Trend_" + s_ModelSelected(iModel) + \
                                    "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileGSAT_PD, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_GSAT_TrendTemp(iReal) = FileGSAT_PD->$s_VarName2$(1)
        ;
        s_VarName2              = "TXx_" + s_ModelSelected(iModel) + "_" + \
                                  s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileGSAT_Proj, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_Temp                 := FileGSAT_Proj->$s_VarName2$
        r_Temp                 := my_SelChina([/r_Temp/])
        r_Temp                 := my_AreaMean2([/r_Temp/], r_Temp&lat, \
                                                 r_Temp&lon, 0)
        ;
        s_VarName2              = "TXx_" + s_ModelSelected(iModel) + \
                                  "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileGSAT_Hist, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_TempHist             := FileGSAT_Hist->$s_VarName2$(\
                                            {i_yrHist(0):i_yrHist(1)}, :, :)
        r_TempHist             := dim_avg_n_Wrap(r_TempHist, 0)
        r_TempHist             := my_SelChina([/r_TempHist/])
        r_TempHist             := my_AreaMean2([/r_TempHist/], r_TempHist&lat, \
                                               r_TempHist&lon, 0)
        ;
        r_Temp                  = r_Temp - r_TempHist
        r_ProjWarmingTemp(iReal, :)= my_3TermsMean(r_Temp, i_yrProj)
    end do
    r_GSAT_Trend(iModel)        = avg(r_GSAT_TrendTemp)
    r_ProjWarming(iModel, :)    = dim_avg_n_Wrap(r_ProjWarmingTemp, 0)
    ;
end do 
;
; print(s_ModelSelected + ": " + r_GSAT_Trend + ", " + r_ProjWarming(:, 2))
; exit()
; =================================================================
;   3. 求两者的简单相关，并绘图
r_Trend_All                   = r_GSAT_Trend
r_GSAT_All                    = r_ProjWarming
;
r_CorCoef                     = escorc_n(r_Trend_All, r_GSAT_All, 0, 0)
print(r_CorCoef + "")
print(s_ModelSelected + ": " + r_Trend_All + ", " + r_GSAT_All(:, 2))
; exit()
;
;; Vis
res                           = True
res@trXMinF                   = 0.14
res@trXMaxF                   = 0.4
if (s_Exp .eq. "ssp126") then
  res@trXMaxF                 = 0.42
end if
res@my_ObsX                   = r_TrendObs
res@ObsX_STD                  = r_TrendSTD

r_trY                         = (/ (/0.3, 1.5/), (/0.9, 2.8/), (/2., 6/) /)
if (s_Exp .eq. "ssp245") then
    r_trY                    := (/ (/0.3, 1.3/), (/0.6, 2.2/), (/1., 3.5/) /)
else if (s_Exp .eq. "ssp126") then
    r_trY                    := (/ (/0.3, 1.2/), (/0.3, 2./), (/.3, 2.2/) /)
else if (s_Exp .eq. "ssp370") then
    r_trY                     = (/ (/0.3, 1.4/), (/0.7, 2.5/), (/1.8, 5/) /)
end if
end if
end if
if (s_Exp .eq. "ssp126") then
    s_LeftString             := (/ "(a1) ", "(b1) ", "(c1) " /)
else if (s_Exp .eq. "ssp245") then
    s_LeftString             := (/ "(a2) ", "(b2) ", "(c2) " /)
else if (s_Exp .eq. "ssp370") then
    s_LeftString             := (/ "(a3) ", "(b3) ", "(c3) " /)
else if (s_Exp .eq. "ssp585") then
    s_LeftString             := (/ "(a4) ", "(b4) ", "(c4) " /)
end if
end if
end if
end if
r_ObsConstraint               = new((/3, 2/), float)
r_ObsConstraint!0             = "ProjectedTerms"
r_ObsConstraint&ProjectedTerms= (/2030, 2060, 2090/)
r_ObsConstraint!1             = "ObsX_ObsY"
r_ObsConstraint&ObsX_ObsY     = (/0, 1/)
r_ObsConstraint@ObsX          = "Multiple-obs average"
r_ObsConstraint_eachModel     = new((/3, dimsizes(s_ModelSelected)/), float)
r_RawGSAT_eachModel           = r_ObsConstraint_eachModel
;
plotAll                       = new(3, graphic)
plotMarkerAll                 = new((/3, dimsizes(s_ModelSelected)+1/), graphic)
plot_RegAll                   = new((/3, 3/), graphic)
dum                           = new((/3, 3/), graphic)
dum2                          = new((/3, dimsizes(s_ModelSelected)+1/), graphic)
dum2Text                      = dum2
do iSub = 0, 2
    ; res@trYMinF               = r_trY(iSub, 0)
    ; res@trYMaxF               = r_trY(iSub, 1)
    res@tiMainString          = s_ExpTitle
    if (iSub .gt. 0) then
        res@tiMainString      = ""
    end if
    res@gsnLeftString         = s_LeftString(iSub) + i_yrProj(iSub, 0) + "-" + \
                                i_yrProj(iSub, 1) ;+ ": " + s_ExpTitle
    l_list := my_plot_Model_Scatter_Special_withConfidenceLines(r_Trend_All, \
                                            r_GSAT_All(:,iSub), wks, res)
    r_ObsConstraint(iSub, 0)  = l_list[0]
    r_ObsConstraint(iSub, 1)  = l_list[1]
    frame(wks)
    ;
    ; return([/r_ObsX, r_ObsY, plot, plotMarker, plot_Reg, plot_Reg2, plot_Reg3, \
    ;    dum, dum2, dum2Text/])
    plotAll(iSub)             = l_list[2] ; plot
    plotMarkerAll(iSub, :)    = l_list[3] ; plotMarker
    plot_RegAll(iSub, 0)      = l_list[4] ; plot_Reg
    plot_RegAll(iSub, 1)      = l_list[5] ; plot_Reg2
    plot_RegAll(iSub, 2)      = l_list[6] ; plot_Reg3
    dum(iSub, :)              = l_list[7]
    dum2(iSub, :)             = l_list[8] ; dum2
    dum2Text(iSub, :)         = l_list[9] ; dum2Text
    ;
    ;; Correct each model 
    r_ObsConstraint_eachModel(iSub, :) = l_list[10]
    r_RawGSAT_eachModel(iSub, :) = l_list[11]
    ;
end do 
;
;; Print the constrained GSAT 
r_MME                         = dim_avg_n_Wrap(r_ObsConstraint_eachModel, 1)
do iSub = 0, 2
    i_Ind                    := ind(ismissing(r_ObsConstraint_eachModel(iSub,:))\
                                     .eq. False)
    r_Range                  := my_QuaterNum2(\
                            r_ObsConstraint_eachModel(iSub, i_Ind), 5, 95)
    print(sprintf("%4.2f", r_MME(iSub)) + " (" + sprintf("%4.2f", r_Range(0)) +\
          "~" + sprintf("%4.2f", r_Range(1)) + ")")
end do
;; panel 
resP                          = True
resP@gsnMaximize              = False
gsn_panel(wks, plotAll, (/3, 1/), resP)
;
end 



















