; ==================================================
; 本程序计算整个中国区域平均下，温度的变化序列，且在其右侧给出三个时期的不确定性和MME
;   1. 指定模式、试验、图形信息、文件信息
;   2. 首先处理温度数据，读取每个模式的结果后，计算MME和5-95th范围（注意历史模拟试验默认读取）
; 然后读取约束后的序列
;   3. 参考之前的脚本绘图，上方为序列图，下方为空间分布图
; 注意本程序参考AR6的Chap4中的检验方法，既考虑模式间的一致性，也考虑内部变率
; 注意本程序空间分布图的颜色，全部改成暖色调
; 注意：本程序计算不确定性范围时，通过样本间的百分位计算，而不是高斯分布假设下计算
; ==================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code/my_FunctionInPNNL.ncl"
;
undef("my_SelChina")
function my_SelChina(Var:list)
; 本程序挑选中国地区，并mask中国以外的区域
; 注意Var中的矩阵必须包含整个中国，至少是大陆地区，注意最后两位必须是lat x lon
;   1. 读取中国区域的mask场，并将其插值到与Var矩阵同样的维度上
;   2. 挑选包含整个中国的范围，然后基于Var中矩阵的维度，调整mask场的维度
;   3. mask中国以外的区域
local r_Var 
begin
    r_Var             = Var[0]
;
    s_DirRead         = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/CN05.1/"
    s_FileRead        = s_DirRead + "CN05.1_Pre_1961_2015_month_025x025.nc"
    FileMask          = addfile(s_FileRead, "r")
    r_Mask            = FileMask->pre(0, :, :)
    
    i_Mask            = where(ismissing(r_Mask) .eq. False, 1, r_Mask)
    copy_VarCoords(r_Mask, i_Mask)
    r_latM            = i_Mask&lat 
    r_lonM            = i_Mask&lon 
    
    i_dims            = dimsizes(r_Var)
    if (dimsizes(i_dims) .eq. 2) then 
        r_VarChina    = r_Var({min(r_latM):max(r_latM)}, \
                              {min(r_lonM):max(r_lonM)})
    else if (dimsizes(i_dims) .eq. 3) then
        r_VarChina    = r_Var(:, {min(r_latM):max(r_latM)}, \
                              {min(r_lonM):max(r_lonM)})
    else if (dimsizes(i_dims) .eq. 4) then
        r_VarChina    = r_Var(:, :, {min(r_latM):max(r_latM)}, \
                              {min(r_lonM):max(r_lonM)})
    end if
    end if
    end if
    ;
    r_latInter       := r_VarChina&lat
    r_lonInter       := r_VarChina&lon
    i_Mask2          := i_Mask({r_latInter}, {r_lonInter}) 
    ;
    if (dimsizes(i_dims) .eq. 2) then 
        i_Mask2_      = i_Mask2
    else if (dimsizes(i_dims) .eq. 3) then
        i_Mask2_      = conform(r_VarChina , i_Mask2, (/1, 2/))
    else if (dimsizes(i_dims) .eq. 4) then
        i_Mask2_      = conform(r_VarChina , i_Mask2, (/2, 3/))
    end if
    end if
    end if
    copy_VarCoords(r_VarChina, i_Mask2_)

    r_VarChina_       = where(i_Mask2_ .ne. 0, r_VarChina, \
                                      r_VarChina@_FillValue)
    copy_VarCoords(r_VarChina, r_VarChina_)
    
    ;; 台湾经纬度范围
    r_lonTaiwan       = (/119, 125/)
    r_latTaiwan       = (/20, 26/)
    r_Temp           := r_VarChina
    r_Temp            = my_mask([/r_Temp/], 0)
    if (dimsizes(i_dims) .eq. 2) then 
        r_VarChina_({r_latTaiwan(0):r_latTaiwan(1)}, \
                        {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(\
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
    else if (dimsizes(i_dims) .eq. 3) then
        r_VarChina_(:, {r_latTaiwan(0):r_latTaiwan(1)}, \
                       {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(:, \
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
    else if (dimsizes(i_dims) .eq. 4) then
        r_VarChina_(:, :, {r_latTaiwan(0):r_latTaiwan(1)}, \
                          {r_lonTaiwan(0):r_lonTaiwan(1)}) = r_Temp(:, :, \
                                            {r_latTaiwan(0):r_latTaiwan(1)}, \
                                            {r_lonTaiwan(0):r_lonTaiwan(1)})
    end if
    end if
    end if
    ;; 台湾经纬度范围
return(r_VarChina_)
;
end
;
undef("my_Add_SouthChinaSea")
function my_Add_SouthChinaSea(map, Draw:numeric, DrawSigLev:numeric, \
                              SigLev:numeric, wks)
; 本程序单独绘制指定图形右下角的南海区域
;; map表示已经绘制完成的plot变量
;; Draw表示被绘制的矩阵
;; DrawSigLev表示显著性水平的矩阵，其维度必须与Draw一致，如果没有，全部设为0即可
;; SigLev表示打点的显著性水平
;; wks表示画布
local r_Draw
begin
;
r_Draw                              = Draw
r_DrawSigLev                        = DrawSigLev
r_SigLev                            = SigLev 

; --------------- add South China Sea --------------
    nhres                           = True
    nhres@gsnFrame                  = False
    nhres@gsnDraw                   = False

    nhres@vpHeightF                 = 0.05
    nhres@vpWidthF                  = 0.05

    nhres@mpMinLatF                 = 2.
    nhres@mpMaxLatF                 = 23.
    nhres@mpMinLonF                 = 105.
    nhres@mpMaxLonF                 = 123.

    getvalues map
    "mpDataSetName"             : nhres@mpDataSetName
    "mpDataBaseVersion"         : nhres@mpDataBaseVersion
    "mpFillOn"                  : nhres@mpFillOn
    "mpFillDrawOrder"           : nhres@mpFillDrawOrder
    "mpAreaMaskingOn"           : nhres@mpAreaMaskingOn
    "mpMaskAreaSpecifiers"      : nhres@mpMaskAreaSpecifiers ;;;;;; very important
    "mpOutlineSpecifiers"       : nhres@mpOutlineSpecifiers
    "mpOutlineBoundarySets"     : nhres@mpOutlineBoundarySets
    "mpLandFillColor"           : nhres@mpLandFillColor
    "mpOceanFillColor"          : nhres@mpOceanFillColor
    "mpInlandWaterFillColor"    : nhres@mpInlandWaterFillColor
    "mpNationalLineColor"       : nhres@mpNationalLineColor
    "mpProvincialLineColor"     : nhres@mpProvincialLineColor
    end getvalues

    getvalues map@contour
    "cnFillOn"                  : nhres@cnFillOn
    "cnLevelSelectionMode"      : nhres@cnLevelSelectionMode
    "cnLevels"                  : nhres@cnLevels
    "cnFillColors"              : nhres@cnFillColors
    "cnFillDrawOrder"           : nhres@cnFillDrawOrder
    "cnLinesOn"                 : nhres@cnLinesOn
    "cnLineLabelsOn"            : nhres@cnLineLabelsOn
    end getvalues

    nhres@lbLabelBarOn              = False
    nhres@tmXBOn                    = False
    nhres@tmYLOn                    = False
    nhres@tmYROn                    = False
    nhres@tmXTOn                    = False
    nhres@gsnLeftString             = ""
    nhres@gsnRightString            = ""

    map_nanhai = gsn_csm_contour_map(wks, r_Draw, nhres)
    ;
    sres    = True
    sres@gsnDraw                    = False
    sres@gsnFrame                   = False
    sres@cnLinesOn                  = False
    sres@cnFillOn                   = False
    sres@cnLineLabelsOn             = False
    sres@cnInfoLabelOn              = False
    sres@gsnLeftString              = ""
    sres@gsnCenterString            = ""
    sres@gsnRightString             = ""
    sres@cnFillDrawOrder            = "PreDraw"
    opt1        = True
    opt1@gsnShadeFillType           = "pattern"
    ; opt1@gsnShadeHigh               = 3
    opt1@gsnShadeLow                = 17
    sres@gsnAddCyclic               = True
    opt1@gsnAddCyclic               = True

    map_stipple                     = gsn_csm_contour(wks, r_DrawSigLev,\
                                                     sres)
    plotSig                         = gsn_contour_shade(map_stipple, r_SigLev, \
                                                       1.6, opt1)
    overlay(map_nanhai, plotSig)
    ;
    ;;; 另一种overlay
    adres                           = True
    adres@amParallelPosF            = 0.495
    adres@amOrthogonalPosF          = 0.49
    adres@amJust                    = "BottomRight"

    plotnh = gsn_add_annotation(map, map_nanhai, adres)
    ; plotSig = gsn_add_annotation(map, plotSig, adres)
    return([/plotnh, plotSig/])
;
end
;
undef("add_y_errorbar_BarChart_DifColor_MME")
procedure add_y_errorbar_BarChart_DifColor_MME(wks, plot, x[*]:numeric, \
                    y[*]:numeric, yerror[*]:numeric, res, ii:numeric)
local dimz, xpoints, ypoints, gres
; 本程序不是挥之竖直线，而是空心框
begin
;
    gres    = res 
    if (gres .eq. False) then
        gres = True
    end if
    dimz    = dimsizes(y)

    gres@gsLineDashPattern = 1

    ; xpoints = (/x, x/)
    ; ypoints = yerror ; (/y-yerror, y+yerror/)
    ; ydot    = y 
    ; plot@$tostring("y_errorbar"+ii)$ = gsn_add_polyline(wks, plot, \
    ;                                     xpoints(:), ypoints(:), gres)
    ;
    r_BarWidth = res@myBarWidth
    ;; up 
    x_up    = (/x - r_BarWidth/2., x + r_BarWidth/2./)
    y_up    = max(yerror)
    y_up   := (/y_up, y_up/)
    plot@$tostring("y_up"+ii)$ = gsn_add_polyline(wks, plot, x_up, y_up, gres)
    ;
    ;; down 
    x_down  = (/x - r_BarWidth/2., x + r_BarWidth/2./)
    y_down  = min(yerror)
    y_down := (/y_down, y_down/)
    plot@$tostring("y_down"+ii)$ = gsn_add_polyline(wks, plot, x_down, y_down, \
                                                    gres)
    ;
    ;; left 
    x_left  = (/x - r_BarWidth/2., x - r_BarWidth/2./)
    y_left  = yerror
    plot@$tostring("y_left"+ii)$ = gsn_add_polyline(wks, plot, x_left, y_left, \
                                                    gres)
    ;
    ;; right
    x_right = (/x + r_BarWidth/2., x + r_BarWidth/2./)
    y_right = yerror
    plot@$tostring("y_right"+ii)$ = gsn_add_polyline(wks, plot, x_right, y_right, gres) 
    ;
    ;; Middle: MME 
    gres@gsLineDashPattern = 0
    gres@gsLineColor = gres@myMME_Color
    x_mid   = (/x - r_BarWidth/2., x + r_BarWidth/2./)
    y_mid   = (/y, y/)
    plot@$tostring("y_mid"+ii)$ = gsn_add_polyline(wks, plot, x_mid, y_mid, gres)
    ;
    ; ;; 补充MME
    ; if (ismissing(gres@gsMakerIndex) .eq. False) then 
    ;     gres@gsMarkerIndex           = 0
    ; end if 
    ; plot@$tostring("y_dot"+ii)$      = gsn_add_polymarker(wks, plot, \
    ;                                             xpoints(:), ydot(:), gres)
;
end
;
undef("my_InnerLegend")
function my_InnerLegend(wks, res, xx_res2, yy_res2)
local res2
begin
;
res2                = res 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Inner Legend 
res2@tmXBOn         = False
plot_res2          := gsn_csm_xy(wks, xx_res2, yy_res2, res2)
;
;; Constrained Projection 
res_text                 = True
res_text@txFontHeightF   = 0.01              ; change text size 
res_text@txJust          = "CenterLeft"      ; text justification
;
pgres               = True
pgres@gsFillColor  := "gray"
pgy                := (/6., 6., 6.5, 6.5, 6. /)
pgx                := (/0, 0.5, 0.5, 0, 0 /)
dum0               := gsn_add_polygon(wks, plot_res2, pgx, pgy, pgres)
plres               = True
plres@gsLineColor  := "black"
r_Temp             := (/min(pgx), max(pgx)/)
r_Temp             := avg(r_Temp)
dum1               := gsn_add_polyline(wks, plot_res2, (/r_Temp, r_Temp/), \
                             (/min(pgy), max(pgy)/), plres)
s_Text              = "Constrained Projection"
r_pgy               = (/min(pgy), max(pgy)/)
r_pgy              := avg(r_pgy)
dum1_text           = gsn_add_text(wks, plot_res2, s_Text, max(pgx)+0.1, \
                                   r_pgy, res_text)
;
;; Raw Projection 
rand1_tst          := toint(systemfunc(" date +%s"))
rand2_tst          := toint((54321l * rand1_tst) % 2147483398l) ;+ iModel
random_setallseed(rand1_tst, rand2_tst)
ii                  = random_uniform( 0, 100, 1 )

pgy                := (/7., 7., 7.5, 7.5, 7. /)
plres@gsLineDashPattern = 1
plres@gsLineColor = "gray"
    ;; up 
    x_up    = (/min(pgx), max(pgx)/)
    y_up    = (/max(pgy), max(pgy)/)
    plot_res2@$tostring("y_up"+ii)$ = gsn_add_polyline(wks, plot_res2, x_up, \
                                                       y_up, plres)
    ;
    ;; down 
    x_down  = (/min(pgx), max(pgx)/)
    y_down  = (/min(pgy), min(pgy)/)
    plot_res2@$tostring("y_down"+ii)$ = gsn_add_polyline(wks, plot_res2, \
                                                         x_down, y_down, plres)
    ;
    ;; left 
    x_left  = (/min(pgx), min(pgx)/)
    y_left  = (/min(pgy), max(pgy)/)
    plot_res2@$tostring("y_left"+ii)$ = gsn_add_polyline(wks, plot_res2, \
                                                    x_left, y_left, plres)
    ;
    ;; right
    x_right = (/max(pgx), max(pgx)/)
    y_right = (/min(pgy), max(pgy)/)
    plot_res2@$tostring("y_right"+ii)$ = gsn_add_polyline(wks, plot_res2, \
                                                    x_right, y_right, plres) 
    ;
    ;; Middle: MME 
    plres@gsLineDashPattern = 0
    plres@gsLineColor = "black"
    x_mid   = (/min(pgx), max(pgx)/)
    x_mid  := avg(x_mid)
    x_mid  := (/x_mid, x_mid/)
    y_mid   = (/min(pgy), max(pgy)/)
    plot_res2@$tostring("y_mid"+ii)$ = gsn_add_polyline(wks, plot_res2, x_mid, \
                                                        y_mid, plres)
s_Text              = "Raw Projection"
r_pgy              := (/min(pgy), max(pgy)/)
r_pgy              := avg(r_pgy)
dum2_text           = gsn_add_text(wks, plot_res2, s_Text, max(pgx)+0.1, \
                                   r_pgy, res_text)

draw(plot_res2)
return(True)
; return([/plot, dum1, dum1_text/])
;
end
;
begin
;
; ==================================================
;   1. 指定模式、试验、图形信息、文件信息
; s_ModelSelect = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", "BCC-CSM2-MR", \
;                 "CAMS-CSM1-0", "CanESM5", "CanESM5-CanOE", "CESM2-WACCM", \
;                 "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg", \
;                 "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
;                 "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", "INM-CM5-0", \
;                 "IPSL-CM6A-LR", "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
;                 "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NESM3", \
;                 "NorESM2-LM", "UKESM1-0-LL"/)
s_ModelSelect             = (/ "ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                         "BCC-CSM2-MR", "CAMS-CSM1-0", "CanESM5", \
                         "CanESM5-CanOE", "CNRM-CM6-1", "CNRM-ESM2-1", \
                         "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                         "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                         "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8", \
                         "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                         "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models                 = "_23Models"
;
s_Exp                     = (/ "ssp126", "ssp245", "ssp370", "ssp585" /)
s_ExpTitleAll             = (/"SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"/) ; 
s_Month                   = "AnMean" ; "JJA" ; "DJF" ; 
s_MonthTitle              = "Annual Mean" ; "JJA" ; "DJF" ; 
i_yrPro                   = (/2015, 2099/)
i_yrHist                  = (/1950, 2014/)
i_yrRun                   = 10
; r_RangePercentage         = 1.645       ; 5th - 95th 
; 1.285 ; 10th - 95th
; 1.645       ; 5th - 95th 

s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Fig/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/04.RunningConstraintAndScaling/"
system("mkdir -p " + s_DirOut)
s_FigOut                  = s_DirOut + \
        "01.Vis_ChinaAreaMean_SeriesWithConstrainedResults_RangeCalByPercentile" + \
        s_Month + s_4Models
wks                       = gsn_open_wks("X11", s_FigOut)
;
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/04.RunningConstraintAndScaling/"
s_FileRead                = s_DirRead + \
            "RunningConstrainedGSAT_ScaleCSAT_2015-2099_Runing" + i_yrRun + \
            "yr_30Models.nc"
FileConstrainedCSAT       = addfile(s_FileRead, "r")
; s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
; s_FileRead                = s_DirRead + "ConstrainedWarmingInChina" + s_4Models + ".nc"
; FileConstrainedCSAT_EachTerm = addfile(s_FileRead, "r")
; ==================================================
;   2. 首先处理温度数据，读取每个模式的结果后，计算MME和10-90th范围
;; Constrained CSAT 
; s_VarNameAll             := getfilevarnames(FileConstrainedCSAT_EachTerm)
;; 读取各个样本的约束结果，并读取每一个不确定性的贡献
do iExp = 0, dimsizes(s_Exp)-1
    s_VarName2            = "ConstrainedCSAT_" + s_Exp(iExp)
    r_TempSeries         := FileConstrainedCSAT->$s_VarName2$(\
                                                        {i_yrPro(0):i_yrPro(1)})
    ;
    ; s_VarName2            = "ConstrainedWarmingChina_" + s_Exp(iExp)
    ; r_Range              := FileConstrainedCSAT_EachTerm->$s_VarName2$
    ;
    ;
     ;; spread across all models and realizations
    s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
    s_FileRead           := s_DirRead + "ConstrainedWarmingInChina_" + \
                            s_Exp(iExp) + s_4Models + ".nc"
    FileConReal          := addfile(s_FileRead, "r")
    ;
    s_VarName2            = "ConstrainedWarmingChina_" + s_Exp(iExp)
    r_Temp               := FileConReal->$s_VarName2$(:, :)
    ;
    s_VarName2            = "RawWarmingChina_" + s_Exp(iExp)
    r_RangeRaw           := FileConReal->$s_VarName2$
    if (iExp .eq. 0) then
        r_ConstrainedCSAT = new((/dimsizes(s_Exp), dimsizes(r_TempSeries)/), \
                                  float)
        copy_VarCoords(r_TempSeries, r_ConstrainedCSAT(0, :))
        r_ConstrainedCSAT!0   = "Exp"
        r_ConstrainedCSAT&Exp = s_Exp
        ;
        ; i_dims           := dimsizes(r_Range)
        ; i_dims(1)         = 3
        i_dims           := array_append_record(dimsizes(s_Exp), \
                                                dimsizes(r_Temp), 0)
        r_ConstrainedCSAT_EachTerm   = new(i_dims, float)
        copy_VarCoords(r_Temp(:, :), r_ConstrainedCSAT_EachTerm(0, :, :))
        copy_VarCoords(r_ConstrainedCSAT(:, 0), \
                       r_ConstrainedCSAT_EachTerm(:, 0, 0))
        r_ConstrainedCSAT_EachTerm!2 = "MME_STD"
        ;
        r_tas_Range       = r_ConstrainedCSAT_EachTerm
    end if
    r_ConstrainedCSAT(iExp, :) = r_TempSeries + 0
    r_ConstrainedCSAT_EachTerm(iExp, :, :) = r_Temp + 0 ;r_Range(:, :) + 0
    r_tas_Range(iExp, :, :)    = r_RangeRaw(:, :) + 0
    ;
end do
;
;; Raw projection 
r_tas_Hist                = new((/dimsizes(s_ModelSelect), \
                                 (i_yrHist(1)-i_yrHist(0)+1)/), float)
r_tas_Hist!0              = "Models"
r_tas_Hist&Models         = s_ModelSelect
r_tas_Hist!1              = "year"
r_tas_Hist&year           = ispan(i_yrHist(0), i_yrHist(1), 1)

; r_pr_Hist                 = r_tas_Hist
;
r_tas_CMIP6               = new((/dimsizes(s_ModelSelect), dimsizes(s_Exp), \
                                 (i_yrPro(1)-i_yrPro(0)+1)/), float)
r_tas_CMIP6!0             = "Models"
r_tas_CMIP6&Models        = s_ModelSelect
r_tas_CMIP6!1             = "Exp"
r_tas_CMIP6&Exp           = s_Exp
r_tas_CMIP6!2             = "year"
r_tas_CMIP6&year          = ispan(i_yrPro(0), i_yrPro(1), 1)  

; r_pr_CMIP6                = r_tas_CMIP6 

s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/01.ChinaMean/"
s_FileRead                = s_DirRead + "tas_ProjectionInChina_" + s_Month + ".nc"
FileTas                   = addfile(s_FileRead, "r")
;
s_FileRead                = s_DirRead + "tas_HistoricalInChina_" + s_Month + ".nc"
FileTas_Hist              = addfile(s_FileRead, "r")
s_VarNameTasAll           = getfilevarnames(FileTas)
do iModel = 0, dimsizes(s_ModelSelect)-1
    ;
    do iExp = 0, dimsizes(s_Exp)-1
        ;
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        ; s_ExpModel_Path:= l_list[0]
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        ; s_ExpModel_Path_:= l_list_[0]
        s_Real_         := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
            s_Real      := s_Real_
            ; s_ExpModel_Path := s_ExpModel_Path_
        end if
        ;
        do iReal = 0, dimsizes(s_Real)-1
            ;; Projected Changes 
            s_VarName2   := "tas_ChinaArea_" + s_ModelSelect(iModel) + "_" + \
                            s_Exp(iExp) + "_" + s_Real(iReal)
            i_ind        := ind(s_VarName2 .eq. s_VarNameTasAll)
            if (ismissing(i_ind(0))) then
                print(s_VarName2 + " Missing!")
                continue
            end if
            r_Temp       := FileTas->$s_VarName2$({i_yrPro(0):i_yrPro(1)})
            ;
            ;; Historical Read 
            s_VarName2   := "tas_ChinaArea_" + s_ModelSelect(iModel) + \
                                    "_historical_" + s_Real(iReal)
            r_TempHist   := FileTas_Hist->$s_VarName2$({i_yrHist(0):i_yrHist(1)})
            ;
            if (isdefined("r_TempAll") .eq. False .or. iReal .eq. 0) then
                r_TempAll  := new((/dimsizes(s_Real), dimsizes(r_Temp)/), float)
                copy_VarCoords(r_Temp, r_TempAll(0, :))
                r_TempAll!0 = "Real"
                r_TempAll&Real = s_Real
                ;
                r_TempHistAll := new((/dimsizes(s_Real), dimsizes(r_TempHist)/), \
                                     float)
                copy_VarCoords(r_TempHist, r_TempHistAll(0, :))
                copy_VarCoords(r_TempAll(:, 0), r_TempHistAll(:, 0))
            end if 
            r_TempAll(iReal, :)= r_Temp
            r_TempHistAll(iReal, :) = r_TempHist
        ;
        end do 
        ;
        if (isdefined("r_TempAll") .eq. False) then
            print(s_ModelSelect(iModel) + " " + s_Exp(iExp) + " Missing!")
            continue
        end if
        ;
        r_tas_CMIP6(iModel, iExp, :) = dim_avg_n_Wrap(r_TempAll, 0)
        if (s_Exp(iExp) .eq. "ssp245") then
            ; printVarSummary(r_tas_Hist)
            ; printVarSummary(r_TempHistAll)
            ; exit()
            r_tas_Hist(iModel, :) = dim_avg_n_Wrap(r_TempHistAll, 0)
        end if
        delete([/r_TempAll, r_TempHistAll/])
        ;
    end do 
    print("")
end do
;
;; Running average 
r_tas_Hist_                := new((/dimsizes(s_ModelSelect), dimsizes(s_Exp),\
                                   i_yrPro(1)-i_yrHist(0)+1/), float)
copy_VarCoords(r_tas_CMIP6(:, :, 0), r_tas_Hist_(:, :, 0))
r_tas_Hist_!2               = "year"
r_tas_Hist_&year            = ispan(i_yrHist(0), i_yrPro(1), 1)
do iExp = 0, dimsizes(s_Exp)-1 
    r_tas_Hist_(:, iExp, {i_yrHist(0):i_yrHist(1)}) = r_tas_Hist + 0
end do 
r_tas_Hist_(:, :, {i_yrPro(0):i_yrPro(1)})    = r_tas_CMIP6(:, :, :)
r_tas_Hist                 := r_tas_Hist_ 

;;; 滑动平均前补全缺测值
r_tas_Hist_                 := dim_avg_n_Wrap(r_tas_Hist_, 0)
r_tas_Hist_                 := conform(r_tas_Hist, r_tas_Hist_, (/1, 2/))
copy_VarCoords(r_tas_Hist, r_tas_Hist_)
r_TempSign                  := r_tas_Hist
r_tas_Hist                   = where(ismissing(r_TempSign), r_tas_Hist_, r_tas_Hist)

r_tas_Hist                   = runave_n_Wrap(r_tas_Hist, i_yrRun, 1, 2)
r_tas_Hist                   = where(ismissing(r_TempSign), \
                                    r_tas_Hist@_FillValue, r_tas_Hist)
r_tas_CMIP6                 := r_tas_Hist(:, :, {i_yrPro(0):})
r_tas_Hist                  := r_tas_Hist(:, 3, {i_yrHist(0):i_yrHist(1)})

;;;;;
;
;; MME and 5th-95th 
r_tas_Hist_2                = r_tas_Hist(:2, :)
r_tas_Hist_2                = r_tas_Hist_2@_FillValue
r_tas_Hist_2(0, :)          = dim_avg_n_Wrap(r_tas_Hist, 0)
; r_tas_Hist_2(1, :)          = dim_stddev_n_Wrap(r_tas_Hist, 0)
; r_tas_Hist_2(2, :)          = r_tas_Hist_2(0, :) + r_tas_Hist_2(1, :) 
; r_tas_Hist_2(1, :)          = r_tas_Hist_2(0, :) - r_tas_Hist_2(1, :) 
; ; printVarSummary(r_tas_Hist_2) ; [Models | 3] x [year | 75]
do iyr = 0, dimsizes(r_tas_Hist_2(0, :))-1
    i_Ind                  := ind(ismissing(r_tas_Hist(:, iyr)) .eq. False)
    r_Temp                 := my_QuaterNum2(r_tas_Hist(i_Ind, iyr), 5, 95)
    r_tas_Hist_2(1:, iyr)   = r_Temp
end do
;
r_tas_CMIP6_2               = r_tas_CMIP6(:2, :, :)
r_tas_CMIP6_2               = r_tas_CMIP6@_FillValue
r_tas_CMIP6_2(0, :, :)      = dim_avg_n_Wrap(r_tas_CMIP6, 0)
; r_tas_CMIP6_2(1, :, :)      = dim_stddev_n_Wrap(r_tas_CMIP6, 0)
; r_tas_CMIP6_2(2, :, :)      = r_tas_CMIP6_2(0, :, :) + r_tas_CMIP6_2(1, :, :) 
; r_tas_CMIP6_2(1, :, :)      = r_tas_CMIP6_2(0, :, :) - r_tas_CMIP6_2(1, :, :) 
; printVarSummary(r_tas_CMIP6) ; [Models | 23] x [Exp | 4] x [year | 85]
do iExp = 0, dimsizes(s_Exp)-1
    do iyr = 0, dimsizes(r_tas_CMIP6_2(0, 0, :))-1
        i_Ind              := ind(ismissing(r_tas_CMIP6(:,iExp,iyr)) .eq. False)
        r_Temp             := my_QuaterNum2(r_tas_CMIP6(i_Ind,iExp,iyr), 5, 95)
        r_tas_CMIP6_2(1:, iExp, iyr)   = r_Temp
    end do 
end do
;
; =================================================
;   3. 参考之前的脚本绘图
; drawNDCGrid(wks)

sspColor                    = (/ (/153, 0, 2/)/255., (/196, 121, 0/)/255.,  \
                                (/112, 160, 205/)/255., (/0, 52, 102/)/255. /)
sspColor                    = sspColor(::-1, :)
sspFill                     = (/ (/252, 209, 197/), (/204, 174, 113/), \
                                   (/146, 197, 222/), (/67, 147, 195/) /) / 255.
sspFill                     = sspFill(::-1, :)
;
res2                        = True
res2@gsnDraw                = False
res2@gsnFrame               = False
res2@xyDashPattern          = 0
    ;
res2@tiYAxisFontHeightF     = 0.01
res2@tiXAxisFontHeightF     = 0.013

res2@gsnCenterString        = ""
res2@gsnRightString         = ""
res2@gsnLeftStringOrthogonalPosF= 0.010
;
res2@tmXTOn                 = False
res2@tmYROn                 = False
res2@tmYLOn                 = True
res2@tmYLLabelsOn           = True
;
res2@xyDashPattern          = 0
res2@vpHeightF              = 0.2
res2@vpWidthF               = 0.3

res2@vpXF                          = 0.7
res2@vpYF                          = 0.8
res2@gsnLeftStringFontHeightF      = 0.012
res2@gsnLeftStringOrthogonalPosF   = 0.01

res2@trXMaxF                = 5.5
res2@trXMinF                = -0.3
res2@trYMaxF                = 8 ; r_YMaxMinF(0,0) ;-0.1
res2@trYMinF                = 0 ; r_YMaxMinF(0,1) ;8
res2@tmXBOn                 = True
res2@tmXBMode               = "Explicit"

res2@tmXBLabelFontHeightF   = 0.011
res2@tmYLLabelFontHeightF   = 0.011

; res2@tiYAxisOn              = False
; res2@tiYAxisString          = "Relative to 1995-2014 ~F21~(~F35~J~F21~C)"
res2@gsnLeftString          = "(b) Mean in Three Projection Terms"
    ;
pgres                       = True
plres                       = True
plres@gsLineDashPattern     = 0
plres@gsLineThicknessF      = 2
    
s_Term                      = (/ "2021-2040     ", "2041-2060     ", \
                                 "2080-2099     " /)
s_Term                      = (/ "Near   ", "Middle     ", "Long   " /)
; printVarSummary(r_tas_Range) ; [Models | 3] x [Exp | 4] x [Projection_Term | 3]
; printVarSummary(r_ConstrainedCSAT_EachTerm) ; [Models | 3] x [Exp | 4] x [ProjectionTerm | 3]
; exit()
xx_res2                     = (/0, 0.4, 0.8, 1.2/)
dxx_res2                    = 2.
yy_res2                     = (/0, 0, 0, 0/)
ii                          = 0
do iTerm = 0, dimsizes(s_Term)-1
    res2@tmXBValues         = xx_res2
    res2@tmXBLabels         = (/ "", "", s_Term(iTerm), "" /)
    plot_res2               = gsn_csm_xy(wks, xx_res2, yy_res2, res2)
    res2@tmYLOn             = False
    s_Region                = ""
    res2@tiYAxisString      = ""
    res2@gsnLeftString      = ""
    ;
    r_yy_Unconstrained     := r_tas_Range(:, iTerm, :)
    ;
    r_yy                   := r_ConstrainedCSAT_EachTerm(:, iTerm, :) ; r_tas_Range(:,:,iTerm)

    r_xx                    = xx_res2 
    dum0                   := new(4, graphic)
    dum0_                  := new(4, graphic)
    dum1                   := new(4, graphic)
    dum1_                  := new(4, graphic)
    dum_CMIP5_0            := new((/4,5/), graphic)
    do issp = 0, dimsizes(s_Exp)-1
        pgres@gsFillColor  := sspFill(issp, :) ;(/153, 0, 2/)/255.
        ;
        ;; Left 
        resVerLines            := True
        resVerLines@gsLineColor:= sspFill(issp, :) ; 
        resVerLines@gsLineThicknessF = 2
        resVerLines@gsMarkerColor    = sspColor(issp, :)
        resVerLines@myMME_Color = sspColor(issp, :)
        resVerLines@myBarWidth  = 0.14
        add_y_errorbar_BarChart_DifColor_MME(wks, plot_res2, r_xx(issp)-0.1, \
                    r_yy_Unconstrained(issp, 0), r_yy_Unconstrained(issp, 1:), \
                    resVerLines, ii)
        ii                  = ii + 1
        ;; Right 
        ;;;;;;;;;;; 修改该部分，不同的不确定性用不同的颜色
        pgy                := (/ r_yy(issp, 1), r_yy(issp, 1), r_yy(issp, 2), \
                                      r_yy(issp, 2), r_yy(issp, 1)/)
        pgx                := (/ r_xx(issp)-0.07, r_xx(issp)+0.07, r_xx(issp)+0.07, \
                                        r_xx(issp)-0.07, r_xx(issp)-0.07  /) + 0.1
        dum0(issp)          = gsn_add_polygon(wks, plot_res2, pgx, pgy, pgres)
        plres@gsLineColor  := sspColor(issp, :)
        dum1(issp)          = gsn_add_polyline(wks, plot_res2, \
                             (/r_xx(issp)-0.07, r_xx(issp)+0.07/) + 0.1, \
                             (/r_yy(issp, 0), r_yy(issp, 0)/), plres)
    end do
        ;
    draw(plot_res2)
    xx_res2                 = xx_res2 + dxx_res2 ;2.5 
end do
L_Return                    = my_InnerLegend(wks, res2, xx_res2, yy_res2)
;;
delete([/res2@tmXBValues, res2@tmXBLabels/])
;
; frame(wks)
; exit()
;
;; 补充左侧的序列图
res2@tmYLLabelsOn           = True
res2@tmYLOn                 = True

;;; 原来的XF：0.7
res2@vpXF                   = 0.3
res2@vpWidthF               = 0.35

res2@trXMinF                = 2015 ; 1980
res2@trXMaxF                = 2095

if (s_MonthTitle .eq. "Annual Mean") then
    s_LeftString            = (/"(a) " + s_MonthTitle + " Temperature Anomaly",\
            "(b) " + s_MonthTitle + " Precipitation Anomaly"/)
    s_RightString           = (/"", ""/)
else 
    s_LeftString            = (/"Temperature Anomaly", "Precipitation Anomaly"/)
    s_RightString           = (/s_MonthTitle, s_MonthTitle/)
end if 

res2@tiYAxisOn              = True
s_tiYAxis                   = (/ "Relative to 1995-2014 ~F21~(~F35~J~F21~C)", \
                                 "Relative to 1995-2014 (%)"/)
l_DrawHist                  = [/r_tas_Hist_2, ""/]
l_DrawPro                   = [/r_tas_CMIP6_2, "" /]
do iFig = 0, 0
    res2@tiYAxisString      = s_tiYAxis(iFig)
    ; res2@vpYF               = res2@vpYF ; r_YF(iFig)
    res2@trYMinF            = 0;r_YMaxMinF(iFig, 0)
    res2@trYMaxF            = 8;r_YMaxMinF(iFig, 1)

    res2@gsnLeftString      = s_LeftString(iFig)
    res2@gsnRightString     = s_RightString(iFig)
    res2@tmYLOn             = True
    res2@tmXBOn             = True

    ; res2@xyLineColor       := -1
    ; res2@gsnXYFillOpacities = 0.2
    ; res2@gsnXYFillColors   := (/128, 128, 128/)/255.

    r_Draw0                := l_DrawHist[iFig]
    xHist                  := r_Draw0&year 
    r_Draw                 := r_Draw0
    ; plotCFill              := gsn_csm_xy(wks, xHist, \
    ;                 r_Draw(1:, {i_yrHist(0):i_yrHist(1)}), res2)
    plotCFill              := gsn_csm_xy(wks, xHist, \
                                    r_Draw(0, {i_yrHist(0):i_yrHist(1)}), res2)
;
    res2@gsnLeftString      = ""
    res2@gsnRightString     = ""
    res2@tmYLOn             = False
    res2@tmXBOn             = False
    
    r_DrawSSP0             := l_DrawPro[iFig]
    xSSP                   := r_DrawSSP0&year 

    ; ;; 只绘制SSP126和SSP585的spread
    ; res2@gsnXYFillColors   := sspFill(0, :)
    ; r_Draw                 := r_DrawSSP0(1:, {"ssp126"}, :)
    ; plot126                := gsn_csm_xy(wks, xSSP, \
    ;                 r_Draw(:, {i_yrPro(0):i_yrPro(1)}), res2)
    ; overlay(plotCFill, plot126)

    ; res2@gsnXYFillColors   := sspFill(3, :)
    ; r_Draw                 := r_DrawSSP0(1:, {"ssp585"}, :)
    ; plot585                := gsn_csm_xy(wks, xSSP, \
    ;                 r_Draw(:, {i_yrPro(0):i_yrPro(1)}), res2)
    ; overlay(plotCFill, plot585)
    ;
    ;; range of each scenario
    plotRange               = new(dimsizes(s_Exp), graphic)
    res2@xyLineThicknesses := (/1., 1./) * .5
    res2@xyDashPatterns    := (/1, 1/)
    do iExp = 0, dimsizes(s_Exp)-1, 3
        res2@xyLineColor   := sspColor(iExp, :)
        r_Draw             := r_DrawSSP0(1:, iExp, {i_yrPro(0):i_yrPro(1)})
        ;
        r_Draw(0, :)        = (r_DrawSSP0(0, iExp, {i_yrPro(0):i_yrPro(1)}) - \
                                r_Draw(0, :)) 
        r_Draw(0, :)        = r_DrawSSP0(0, iExp, {i_yrPro(0):i_yrPro(1)}) - \
                                r_Draw(0, :)

        r_Draw(1, :)        = (r_Draw(1, :) - \
            r_DrawSSP0(0, iExp, {i_yrPro(0):i_yrPro(1)})) 
        r_Draw(1, :)        = r_DrawSSP0(0, iExp, {i_yrPro(0):i_yrPro(1)}) + \
                                r_Draw(1, :)
        ;
        plotRange(iExp)     = gsn_csm_xy(wks, xSSP, r_Draw, res2)
        overlay(plotCFill, plotRange(iExp))
    end do
    delete([/res2@xyDashPatterns/])
;
    ;; 绘制每条序列
    res2@xyLineThicknesses := 0.5
    delete([/res2@gsnXYFillColors/])
    res2@xyLineColor       := (/0., 0., 0./)

    r_Draw                 := r_Draw0(0, :) ;r_tas_Hist_2(0, :)
    plotSerHist            := gsn_csm_xy(wks, xHist, r_Draw, res2)
    overlay(plotCFill, plotSerHist)
;
    delete([/res2@xyLineColor/])
    res2@xyLineThicknesses := conform_dims(dimsizes(s_Exp), 1., 0) ;(/5.0, 5.0, 5.0, 5.0/)
    res2@xyLineColors      := sspColor
    r_Draw                 := r_DrawSSP0(0, :, :) ;r_tas_CMIP6_2(0, :, :)
    plotSSP                := gsn_csm_xy(wks, xSSP, r_Draw, res2)
    overlay(plotCFill, plotSSP)
    ;
    ;; Constrained CSAT 
    i_Temp                 := r_ConstrainedCSAT&year
    res2@xyLineThicknesses  = conform_dims(dimsizes(s_Exp), 2.5, 0)
    plotSSP_Constrained     = gsn_csm_xy(wks, i_Temp, r_ConstrainedCSAT, res2)
    overlay(plotCFill, plotSSP_Constrained)
    ;
    delete([/res2@xyLineColors/])
;
    ;; Shading different periods
      gsres = True
      gsres@gsFillColor           = "gray80"
      gsres@gsFillOpacityF        = 0.2
      gsres@tfPolyDrawOrder       = "Predraw"
      txres                       = True
      txres@txFontHeightF         = 0.01
      txres@txFontThickness       = 5
      txres@txFontColor           = "black"
      txres@tfPolyDrawOrder       = "Postdraw"
      ; drawNDCGrid(wks)
      r_txY                       = (/ 0.3, -6 /)

      polygon1 = gsn_add_polygon(wks,plotCFill,(/2021, 2021, 2040, 2040, 2021/), \
                                 (/res2@trYMinF, res2@trYMaxF, res2@trYMaxF, res2@trYMinF, res2@trYMinF/), gsres)
      polytext1= gsn_add_text(wks, plotCFill, "Near", 2031, r_txY(iFig), txres)
      polygon2 = gsn_add_polygon(wks,plotCFill,(/2041, 2041, 2060, 2060, 2041/), \
                                 (/res2@trYMinF, res2@trYMaxF, res2@trYMaxF, res2@trYMinF, res2@trYMinF/), gsres)
      polytext2= gsn_add_text(wks, plotCFill, "Middle", 2051, r_txY(iFig), txres)
      polygon3 = gsn_add_polygon(wks,plotCFill,(/2080, 2080, 2099, 2099, 2080/), \
                                 (/res2@trYMinF, res2@trYMaxF, res2@trYMaxF, res2@trYMinF, res2@trYMinF/), gsres)
      polytext3= gsn_add_text(wks, plotCFill, "Long", 2087.5, r_txY(iFig), txres)
    ;
    ; Create legend
      lgres                    = True
      sspColor2                = new((/dimsizes(s_Exp)+1, 3/), float)
      sspColor2(1:,:)          = sspColor
      ; sspColor2(1,:)           = (/128, 128, 128/)/255.
      sspColor2(0,:)           = (/0, 0, 0/)/255.
      lgres@lgLineColors      := sspColor2(1:, :)
      lgres@lgLineThicknessF   = conform_dims(dimsizes(sspColor2(:,0)), 3, 0) ;(/3, 3, 3/) 
      lgres@lgLabelFontHeightF = .1           ; set the legend label font thickness
      lgres@vpWidthF           = 0.1          ; width of legend (NDC)
      lgres@vpHeightF          = 0.13            ; height of legend (NDC)，调节Legend之间的垂直间距
      lgres@lgMonoDashIndex    = True           ; Same line dash patterns  
      lgres@lgPerimOn          = False          ; draw the box perimeter in orange
      ; labels := " " + (/ "GPCP", "Historical" /)
      ; labels := array_append_record(labels, s_ExpTitle, 0)
      labels                   = " " + (/"Historical", "SSP1-2.6", "SSP2-4.5", \
                                    "SSP3-7.0", "SSP5-8.5" /)
      lbid    = gsn_create_legend(wks, dimsizes(labels)-1, labels(1:), lgres)
    
      ; Add legend to bottom right corner of plot
      amres                  = True
      amres@amJust           = "TopLeft"  ; Corner of plot for positioning
      amres@amOrthogonalPosF =  -0.46   ;  0.5 is the bottom edge of the plot.
      amres@amParallelPosF   =  -0.45   ;  0.5 is the right edge of the plot.
      annoid  = gsn_add_annotation(plotCFill, lbid, amres)
    ;
    draw(plotCFill)
end do 
;
frame(wks)
print(s_FigOut + "")
;
exit()
;; for Pattern 
;
res                 = True
res@gsnDraw         = False
res@gsnFrame        = False
; --------------------------------------------------------
; set for the map 
res@mpMinLatF       = 15.
res@mpMaxLatF       = 55.
res@mpMinLonF       = 70.
res@mpMaxLonF       = 140.

res@mpFillOn        = True
res@mpDataSetName   = "/WORK/chenzm/my_data/topo/NCL-Chinamap-master/database/Earth..4"
res@mpDataBaseVersion           = "MediumRes" ; or "Ncarg4_1"
res@mpAreaMaskingOn             = True
res@mpMaskAreaSpecifiers        = (/ "China" /)
res@mpOutlineSpecifiers         = (/ "China", "China:Provinces" /)

res@mpLandFillColor             = "white"
res@mpInlandWaterFillColor      = "white"
res@mpOceanFillColor            = "white"
res@mpOutlineBoundarySets       = "NoBoundaries"
res@mpFillBoundarySets          = "NoBoundaries"
res@mpNationalLineColor         = "black"
res@mpProvincialLineColor       = "gray" ;"black"
res@mpGeophysicalLineColor      = "black"
res@mpNationalLineThicknessF    = 2
res@mpProvincialLineThicknessF  = 1

; -------------------------------------------------------
; set for the polt 
res@cnFillOn                    = True
res@cnFillDrawOrder             = "PreDraw"
res@cnLinesOn                   = False
res@cnLineLabelsOn              = False

res@pmTickMarkDisplayMode       = "Always"

; cmap   := (/ (/103, 0, 31/), (/178, 24, 43/), (/214, 96, 77/), \
;              (/244, 165, 130/), (/253, 219, 219/), (/240, 240, 240/),\ ; (/247, 247, 247/), \
;              (/209, 229, 240/), (/146, 197, 222/), (/67, 147, 195/), \
;              (/33, 102, 172/), (/5, 48, 97/) /) / 255. ; 
; cmap                           := cmap(::-1, :)

r_cnLinesLevel_c               := (/1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6/)
i_dimColor                      = dimsizes(r_cnLinesLevel_c)

cmap                            = read_colormap_file("cmocean_matter")
i_dimColor                      = dimsizes(cmap(:, 0)) / i_dimColor
cmap                           := cmap(::i_dimColor, :)


res@cnFillColors               := cmap
res@cnLevelSelectionMode        = "ExplicitLevels"    
res@cnLevels                   := r_cnLinesLevel_c
res@tmXBLabelFontHeightF        = 0.009
res@tmYLLabelFontHeightF        = 0.009

; ---------------------------- Stipple ------------------------------
sres                            = True ; res
sres@cnFillDrawOrder            = "PreDraw"
sres@gsnDraw                    = False
sres@gsnFrame                   = False
sres@cnLinesOn                  = False
sres@cnFillOn                   = False
sres@cnLineLabelsOn             = False
sres@cnInfoLabelOn              = False
sres@gsnLeftString              = ""
sres@gsnRightString             = "" ;

opt1                            = True
opt1@cnFillDrawOrder            = "PreDraw"
opt1@gsnShadeFillType           = "pattern"

s_LeftString                   := (/ "(b)", "(c)", "(d)", "(e)" /)

res@vpWidthF                    = 0.3
res@vpHeightF                   = 0.2
r_vpXYF                         = (/ (/res2@vpXF, 0.62/), \
                                 (/res2@vpXF+0.35, 0.62/), \
                                 (/res2@vpXF, 0.38/), \
                                 (/res2@vpXF+0.35, 0.38/) /)
printVarSummary(r_VarAllExp)
printVarSummary(r_SigAllExp)
map_stipple                     = new(dimsizes(s_Exp), graphic)
map_stipple2                    = new(dimsizes(s_Exp), graphic)
map                             = map_stipple
plotnh                          = map 
plotnhSig                       = map 
do iFig = 0, dimsizes(s_Exp)-1
    res@vpXF                    = r_vpXYF(iFig, 0)
    res@vpYF                    = r_vpXYF(iFig, 1)
    res@gsnLeftString           = s_LeftString(iFig) + " " + s_ExpTitleAll(iFig)
    res@gsnRightString          = i_yrPattern(0) + "-" + i_yrPattern(1)
    ;
    res@lbLabelBarOn            = False
    if (iFig .eq. dimsizes(s_Exp)-1) then
        res@lbLabelBarOn        = True
        res@pmLabelBarOrthogonalPosF    = 0.2
        res@pmLabelBarParallelPosF      = -0.1
        res@pmLabelBarWidthF    = 0.6
        res@lbLabelFontHeightF  = 0.01
    end if
    map(iFig)     = gsn_csm_contour_map(wks, r_VarAllExp(iFig, :, :), res)

    res@gsnRightString          = ""
    sres@cnFillColor            = "black"
;
    ;; Model Agreements 
    opt1@gsnShadeLow            = 3
    r_Consistent                = 0.80
    r_Temp                     := where(r_SigAllExp(iFig, :, :).ge.r_Consistent, \
                                        1.1, 0)
    copy_VarCoords(r_SigAllExp(iFig, :, :), r_Temp)
    r_Temp(0, 0)                = 0
    map_stipple(iFig)           = gsn_csm_contour(wks, r_Temp, sres)
    map_stipple(iFig)           = gsn_contour_shade(map_stipple(iFig), 0, \
                                                    1, opt1)
    overlay(map(iFig), map_stipple(iFig))
    ;
    ;; Internal Variability 
    opt1@gsnShadeLow            = 12
    r_Consistent                = 0.66
    r_Temp                     := where(\
                        r_SigAllExp_InterVar(iFig, :, :) .ge. r_Consistent, \
                        1.1, 0)
    copy_VarCoords(r_SigAllExp_InterVar(iFig, :, :), r_Temp)
    map_stipple2(iFig)          = gsn_csm_contour(wks, r_Temp, sres)
    map_stipple2(iFig)          = gsn_contour_shade(map_stipple2(iFig), 0, \
                                                    1, opt1)
    overlay(map(iFig), map_stipple2(iFig))
    ;
    l_list                     := my_Add_SouthChinaSea(map(iFig), \
                                                    r_VarAllExp(iFig, :, :), \
                                                    r_SigAllExp(iFig, :, :), \
                                                    r_Consistent, wks)
    plotnh(iFig)                = l_list[0]
    plotnhSig(iFig)             = l_list[1]
    draw(map(iFig))
end do

; drawNDCGrid(wks)
;
frame(wks)
;
end











