; ==================================================
; 本程序基于GSAT与局地升温的回归系数，对中国地区的升温进行约束
;   1. 读取回归系数，然后读取约束后的GSAT，两者相乘和加上截距
;   2. 读取约束前的结果，然后约束前后作差
;   3. 进行可视化
; 注意本程序在旁边附上约束前后的差异，此外本程序是针对每个模式先多样本集合平均后，再约束的结果
; 注意本程序采用三个预估时期，GSAT与局地升温的回归系数，作为scaling的系数
; 注意本程序是用于完美模式检验，即需要制定虚假观测的模式，然后在计算scaling pattern
; 的时候，扣除该模式的回归分布
; 注意本程序首先订正每一个模式的tasmax，然后进行集合平均，由此避免截距的影响
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
    s_DirRead         = "$WACCEM/my_data/Obs/CN05.1/"
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
    r_lonTaiwan       = (/120, 125/)
    r_latTaiwan       = (/20, 25/)
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
begin
;
; ==================================================
;   1. 读取回归系数，然后读取约束后的GSAT，注意别用y = kx + b的方法，以免受到截距的影响
if (isdefined("s_PseudoModel") .eq. False) then
    s_PseudoModel         = "CanESM5" ; "MRI-ESM2-0"
end if
;
s_4Models                 = "_18Models" ; "_5Models" ; "_30Models" ; "" ; 
s_4ModelsGSAT_Con         = "_23Models" ; Keep it as a largest ensemble
; s_ModelSelected           = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
;                             "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_ModelSelected           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
;
s_Month                   = "AnMean"
if (isdefined("s_Exp") .eq. False) then
    s_Exp                 = "ssp245" ; "ssp126" ; "ssp370" ; "ssp585" ; 
end if
if (s_Exp .eq. "ssp126") then
    s_ModelSelected      := (/"ACCESS-CM2", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
;
i_yrProj                  = (/ (/2021, 2040/), (/2041, 2060/), \
                                   (/2080, 2099/) /)
s_yrProj                  = (/"near", "mid", "long"/)

; ;;;;; for Pattern 
; s_VarName                 = "tas" ; "pr" ; 
; s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/ProjectionOfAllReal/"
; s_FileRead                = s_DirRead + "ProjectionInDifTerms_AllReal" + \
;                             s_4Models + ".nc"
; FileRead_Raw              = addfile(s_FileRead, "r")
; ;;;;; for Pattern 
s_DirRead   = "$WACCEM/my_data/CMIP6/tasmax/POST/TXxRegGSAT/InTimeDim/"
s_FileRead              = s_DirRead + "GSAT_reg_TXxWarmingInTimeDim_AnMean_" + \
                            s_Exp + "_20yr" + s_4Models + ".nc"
FileRead                = addfile(s_FileRead, "r")
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    ; Reg Coef 
    s_VarName2          = "RegInTimDim_" + s_ModelSelected(iModel) + "_" + s_Exp
    if (my_VarNameExist(FileRead, s_VarName2) .eq. False) then
        print(s_VarName2 + "  Missing!") 
        continue
    end if
    r_RegTemp          := FileRead->$s_VarName2$(0, :, :, :)
    ; s_VarName2          = "yintercept_" + s_ModelSelected(iModel) + "_" + s_Exp
    ; r_yintercept_Temp  := FileRead->$s_VarName2$
    ;
    if (isdefined("r_Reg") .eq. False) then
        i_dims         := array_append_record(dimsizes(s_ModelSelected), \
                                              dimsizes(r_RegTemp), 0)
        r_Reg           = new(i_dims, float)
        copy_VarCoords(r_RegTemp, r_Reg(0, :, :, :))
        r_Reg!0         = "Models"
        r_Reg&Models    = s_ModelSelected
        ; r_yintercept    = r_Reg
    end if
    r_Reg(iModel,:,:,:) = r_RegTemp
end do
r_Reg                           := dim_avg_n(r_Reg, 0)
;
;; For GSAT and correction
s_DirRead = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tas/POST/GMSAT/MultiMembers/"
s_FileRead                       = s_DirRead + "GMSAT_Series_MultiMember_" + \
                                    s_Exp + ".nc"
FileGMSAT                        = addfile(s_FileRead, "r")
s_FileRead                       = s_DirRead + "GMSAT_Series_MultiMember_historical.nc"
FileGMSAT_Hist                   = addfile(s_FileRead, "r")
;
;; Constrained GSAT 
s_DirRead   = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
s_FileRead              = s_DirRead + "ConstrainedGSAT_inThreeTerms_PseModel_" \
                            + s_PseudoModel + s_4Models + ".nc"
FileConstrainedGSAT     = addfile(s_FileRead, "r")
s_VarName2              = "ConstrainedGSAT_MME_" + s_Exp
r_GSAT                  = FileConstrainedGSAT->$s_VarName2$(:)
print("Constrained GSAT: " + r_GSAT + "")
;
s_VarName                        = "TXx" ; "pr" ; 
s_DirRead = "$WACCEM/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    ;
    l_list                  := my_CMIP6Model_Exp_Path_ForAllReal(\
                                        s_ModelSelected(iModel), s_Exp, "day")
    s_Real                  := l_list[1]
    l_list_                 := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                        s_ModelSelected(iModel), s_Exp, "day")
    s_Real_                 := l_list_[1]
    if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
        s_Real              := s_Real_
    end if
        ;
    if (all(ismissing(s_Real)) ) then
        print(s_ModelSelected(iModel) + ", " + s_Exp + " Missing!")
        continue
    end if
    ;
    ;; Read TXx
    s_FileRead               = s_DirRead + s_VarName + "_" + s_Month + "_"+\
                                s_Exp + "_" + s_ModelSelected(iModel) + \
                                "_CMIP6_AllReal.nc"
    FileRaw                 := addfile(s_FileRead, "r")
    s_FileRead               = s_DirRead + s_VarName + "_" + s_Month + \
                                "_historical_" + s_ModelSelected(iModel) + \
                                "_CMIP6_AllReal.nc"
    FileRawHist             := addfile(s_FileRead, "r")
        ;
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2          := s_VarName + "_" + s_ModelSelected(iModel) + \
                                "_" + s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileRaw, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue 
        end if
        r_Temp              := FileRaw->$s_VarName2$(:, {-90:90}, {0:360})
            ;
        s_VarName2          := s_VarName + "_" + s_ModelSelected(iModel) + \
                                    "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileRawHist, s_VarName2)  .eq. False) then
            print(s_VarName2 + " Missing!")
            continue 
        end if
        r_TempHist          := FileRawHist->$s_VarName2$({1995:2014}, \
                                                             {-90:90}, {0:360})
        r_TempHist          := dim_avg_n_Wrap(r_TempHist, 0)
        r_TempHist          := conform(r_Temp, r_TempHist, (/1, 2/))
        copy_VarCoords(r_Temp, r_TempHist)
        r_Temp               = r_Temp - r_TempHist
            ;
        ;; GMSAT 
        s_VarName2           = "GMSAT_" + s_ModelSelected(iModel) + "_" + \
                                    s_Exp + "_" + s_Real(iReal)
        if (my_VarNameExist(FileGMSAT, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_TempGMSAT         := FileGMSAT->$s_VarName2$
        r_TempGMSAT         := dim_avg_n_Wrap(r_TempGMSAT, 1)
        s_VarName2           = "GMSAT_" + s_ModelSelected(iModel) + \
                                    "_historical_" + s_Real(iReal)
        r_TempGMSAT_Hist    := FileGMSAT_Hist->$s_VarName2$({1995:2014}, :)
        r_TempGMSAT_Hist    := dim_avg_n_Wrap(r_TempGMSAT_Hist, 0)
        r_TempGMSAT_Hist    := dim_avg_n_Wrap(r_TempGMSAT_Hist, 0)
        r_TempGMSAT          = r_TempGMSAT - r_TempGMSAT_Hist
            ;
        if (isdefined("r_TempAllReal") .eq. False .or. iReal .eq. 0) then
            i_dims          := array_append_record(dimsizes(s_Real), \
                                                    dimsizes(r_Temp(:2,:,:)), 0)
            r_TempAllReal   := new(i_dims, float)
            copy_VarCoords(r_Temp(0, :, :), r_TempAllReal(0, 0, :, :))
            r_TempAllReal!0  = "Real"
            r_TempAllReal&Real  = s_Real
            r_TempAllReal!1  = "Terms"
            r_TempAllReal&Terms = ispan(0, 2, 1) ; (/"near", "mid", "long"/)
                ;
            r_TempGMSAT_AllReal:= r_TempAllReal(:, :, 0, 0)
        end if
            ;
        do iTerm = 0, dimsizes(i_yrProj(:, 0))-1
            i_TermYr        := i_yrProj(iTerm, :)
            r_TempAllReal(iReal, iTerm, :, :) = dim_avg_n_Wrap(\
                                    r_Temp({i_TermYr(0):i_TermYr(1)}, :, :), 0)
                ;
            r_TempGMSAT_AllReal(iReal, iTerm) = dim_avg_n_Wrap(\
                                    r_TempGMSAT({i_TermYr(0):i_TermYr(1)}), 0)
        end do
            ;
    end do 
    ;
    if (isdefined("r_TempAllReal") ) then 
        if (isdefined("r_VarRaw") .eq. False) then
            i_dims          := array_append_record(dimsizes(s_ModelSelected), \
                                        dimsizes(r_TempAllReal(0, :2, :, :)), 0)
            r_VarRaw        := new(i_dims, float)
            copy_VarCoords(r_TempAllReal(0, :, :, :), r_VarRaw(0, :, :, :))
        ;
            r_VarRaw!0       = "Models"
            r_VarRaw&Models  = s_ModelSelected
        ;
            r_VarCorForAllModels    = r_VarRaw
        end if   
        ; r_VarRaw(iModel, :, :, :)   = dim_avg_n_Wrap(r_TempAllReal, 0)
            ;
        if (s_ModelSelected(iModel) .eq. s_PseudoModel) then
            print(s_ModelSelected(iModel) + " is pseudo-obs model!")
            r_PseudoWarming := dim_avg_n_Wrap(r_TempAllReal, 0)
        else
            r_TempGMSAT_AllReal        := dim_avg_n_Wrap(r_TempGMSAT_AllReal, 0) 
            r_VarRaw(iModel, :, :, :)   = dim_avg_n_Wrap(r_TempAllReal, 0)
            ;
            dx                         := r_TempGMSAT_AllReal - r_GSAT
            do iW = 0, 2
                r_VarCorForAllModels(iModel, iW, :, :) = \
                            r_VarRaw(iModel, iW, :, :) - r_Reg(iW, :, :) * dx(iW) 
            end do
        end if
        delete([/r_TempAllReal/])
    end if 
    ;
    ; r_yintercept(iModel, :, :, :) = r_yintercept_Temp
end do
print("")
;
r_Warming              := dim_avg_n_Wrap(r_VarCorForAllModels, 0)
;
;;; Area-Mean over China 
r_China                := my_SelChina([/r_Warming/])
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
print("3 Term Warming:")
print(r_China + " K")
;
;; For Output 
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/06.PerfectModelTest/02.ForExtremeWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut                        = s_DirOut + \
                "ConstrainedTXxWarmingInChina_RegMMEForScaling_PseModel_" + \
                s_PseudoModel + s_4Models + ".nc"
ls                               = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut                      = addfile(s_FileOut, "c")
    GAttr                        = True
    GAttr@history                = systemfunc("date")+";  "+systemfunc("pwd") +\
            "/02.CorTXxWarmingByScaling_withDif_RegMMEForScaling_ForPMT.ncl"
    GAttr@Modles                 = s_ModelSelected
    fileattdef(FileOut, GAttr)
else 
    FileOut                      = addfile(s_FileOut, "w")
end if
;
;; 单独计算中国地区的区域平均，以及spread
do iModel = 0, dimsizes(s_ModelSelected)-1
    if (s_ModelSelected(iModel) .eq. s_PseudoModel) then
        print(s_ModelSelected(iModel) + " is pseudo-obs model!")
        print("Skip")
        continue
    end if
    ;
    ; s_VarName2          = "ConstrainedGSAT_" + s_ModelSelected(iModel) + "_" + \
    ;                         s_Exp
    ; L_Return            = my_VarNameExist(FileConstrainedGSAT, s_VarName2)
    ; r_Reg_Warming      := FileConstrainedGSAT->$s_VarName2$
    ; print(s_ModelSelected(iModel) + ": " + r_Reg_Warming + "")
    ;
    ; r_Reg_Warming      := conform(r_Reg, r_Reg_Warming, (/0/))
    ; copy_VarCoords(r_Reg, r_Reg_Warming)
    ; r_Reg_Warming       = r_Reg * r_Reg_Warming + r_yintercept
    r_Reg_Warming      := r_VarCorForAllModels(iModel, :, :, :)
    ;
    s_VarName2          = "ConstrainedWarmingChina_" + \
                        s_ModelSelected(iModel) + s_Exp
    FileOut->$s_VarName2$ = r_Reg_Warming
    ;
    r_Reg_Warming      := my_SelChina([/r_Reg_Warming/])
    ;
    r_Reg_Warming      := my_AreaMean2([/r_Reg_Warming/], r_Reg_Warming&lat, \
                                       r_Reg_Warming&lon, 0)
    if (isdefined("r_Reg_Warming_All") .eq. False) then
        r_Reg_Warming_All        = new((/dimsizes(s_ModelSelected), 3/), float)
        copy_VarCoords(r_Reg_Warming, r_Reg_Warming_All(0, :))
        r_Reg_Warming_All!0      = "Models"
        r_Reg_Warming_All&Models = s_ModelSelected
    end if
    r_Reg_Warming_All(iModel, :) = r_Reg_Warming
end do
print(r_Reg_Warming_All&Models + ": " + r_Reg_Warming_All(:, 2))
r_Reg_Warming_Spread             = dim_stddev_n_Wrap(r_Reg_Warming_All, 0)
print("STD: " + r_Reg_Warming_Spread(2))
;
r_China_MME_Spread               = new((/3, 2/), float)
copy_VarMeta(r_China, r_China_MME_Spread(:, 0))
r_China_MME_Spread!1             = "MME_STD"
r_China_MME_Spread&MME_STD       = (/0, 1/)
do iPer = 0, 2
    r_China_MME_Spread(iPer, 0)  = r_China(iPer)
    r_China_MME_Spread(iPer, 1)  = r_Reg_Warming_Spread(iPer)
end do
; printVarSummary(r_China_MME_Spread) ; [ProjectionTerm | 3] x [MME_STD | 2]
;
;; Output the Constrained MME 
s_VarName2                       = "ConstrainedWarmingChina_" + s_Exp
FileOut->$s_VarName2$            = r_China_MME_Spread
s_VarName2                       = "ConstrainedWarmingChina_MME_" + s_Exp
FileOut->$s_VarName2$            = r_Warming
;
    r_VarRawMME        := dim_avg_n_Wrap(r_VarRaw, 0)
;
r_WarmingDif            = r_Warming
r_WarmingDif            = r_WarmingDif - r_PseudoWarming
r_VarRawMME0            = r_VarRawMME
r_VarRawMME             = r_VarRawMME - r_PseudoWarming
;
;; Area-mean of difference over China 
r_China                := my_SelChina([/r_WarmingDif/])
r_China                 = r_China^2
r_China                := my_AreaMean2([/r_China/], r_China&lat, r_China&lon, 0)
r_China                 = sqrt(r_China)
;
r_China_Raw            := my_SelChina([/r_VarRawMME/])
r_China_Raw             = r_China_Raw^2
r_China_Raw            := my_AreaMean2([/r_China_Raw/], r_China_Raw&lat, \
                                       r_China_Raw&lon, 0)
r_China_Raw             = sqrt(r_China_Raw)
;
;; PCC 
r_PCC                   = new(3, float)
r_PCC_Raw               = r_PCC
r_PseudoWarming        := my_SelChina([/r_PseudoWarming/])
r_Warming              := my_SelChina([/r_Warming/])
r_VarRawMME0           := my_SelChina([/r_VarRawMME0/])
do iTerm = 0, dimsizes(s_yrProj)-1
    l_list              = my_PCC_RMSE_for2D0(r_PseudoWarming(iTerm, :, :), \
                                             r_Warming(iTerm, :, :), 0)
    r_PCC(iTerm)        = l_list[0]
    l_list             := my_PCC_RMSE_for2D0(r_PseudoWarming(iTerm, :, :), \
                                             r_VarRawMME0(iTerm, :, :), 0)
    r_PCC_Raw(iTerm)    = l_list[0]
end do
print("RMSE of Constrained and Raw projection in 3 Term Warming:")
print(r_China + " K, " + r_China_Raw + " K")
;
print("PCC of Constrained and Raw projection in 3 Term Warming:")
print(r_PCC + ", " + r_PCC_Raw)

; s_FileTest = "./Test.nc"
; system("rm -rf " + s_FileTest)
; FileTest                = addfile(s_FileTest, "c")
; FileTest->Constrained   = r_Warming
; FileTest->Raw           = r_VarRawMME0
; FileTest->Pseudo        = r_PseudoWarming
;
print("")
print(s_FileOut + "")
;
end



































