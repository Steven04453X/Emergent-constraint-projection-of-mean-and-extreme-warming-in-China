; =================================================
; Author: Ziming Chen 
; This code is for:
;   - estimating the internal variability following Simpson (2021)
;   1. indicate models, LENS models with ensemble members, exp 
;   2. read the constrained projection of LENS for each realization
;   3. for each model, get the number of ensembel members (n_mem), and then 
; subsampling the n_mem members of each LENS and repeating 1000 times 
;   4. estimating the internal var by using the STD of 1000 sythentic members 
;   5. make inter-model ensemble of internal var and output 
; =================================================
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
; =================================================
;   1. indicate models, LENS models with ensemble members, exp 
s_ModelSelected         = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models               = "_18Models"
s_4ModelsGSAT_Con       = "_23Models" ; Keep it as a largest ensemble
s_ModelLENS             = (/"CanESM5", "MIROC6"/)
;
if (isdefined("s_Exp") .eq. False) then
    s_Exp               = "ssp585"
end if
if (s_Exp .eq. "ssp126") then
    s_ModelSelected    := (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                        "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                        "CNRM-ESM2-1", "CanESM5", "EC-Earth3", \
                        "EC-Earth3-Veg", "FGOALS-g3", "GFDL-CM4", \
                        "GFDL-ESM4", "IPSL-CM6A-LR", "INM-CM4-8", \
                        "MIROC6", "MRI-ESM2-0", "UKESM1-0-LL"/)
end if
;
i_yrProj                = (/(/2021, 2040/), (/2041, 2060/), (/2080, 2099/)/)
s_yrProj                = (/"near", "mid", "long"/)
s_VarName               = "TXx"
;
i_SubSamp               = 1000
;
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
s_FileOut               = s_DirOut + "EstInterVar_fromEachLENS_1000times.nc"
ls                      = systemfunc("ls " + s_FileOut)
if (ismissing(ls(0))) then
    FileOut             = addfile(s_FileOut, "c")
    GAttr               = True
    GAttr@history       = systemfunc("date") + "  " + systemfunc("pwd") + \
                        "/03.EstInterVar_fromEachLENS_1000times.ncl"
    GAttr@Des           = "For a CMIP model with nf ensemble members for the future," + \
                        " the internal variability is estimated by the subsmapling nf" + \
                        " members (with replacement) from the future periods of" + \
                        " each large ensemble (LEN) and repeating 1000 times. The LEN" + \
                        " ensemble mean is then subtracted from the mean of each subsample." + \
                        " Because the ensemble mean for each LE will be much closer to the true for" + \
                        " that model than the subsamples when nf is small, these 1000" + \
                        " anomalies can be considered to represent 1000 deviations from the ture value" + \
                        " that could arise due to sampling of internal var with only nf members."
    fileattdef(FileOut, GAttr)
else 
    FileOut             = addfile(s_FileOut, "w")
end if
; =================================================
;   2. read the constrained projection of LENS for each realization
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/02.ForExtremeWarming/04.ConEachRealForInterVarEst/"
s_FileRead              = s_DirRead + \
                    "FurtherCorrectConstrainedWarming_ForEachReal_" + s_Exp + \
                    s_4Models + ".nc"
FileLENS                = addfile(s_FileRead, "r")
do iLENS = 0, dimsizes(s_ModelLENS)-1
    ;; LENS 
    s_VarName2         := "ConstrainedWarming_" + s_ModelLENS(iLENS) + "_" + \
                            s_Exp + "_AllReal"
    r_LENS             := FileLENS->$s_VarName2$
    ;
    l_list             := my_CMIP6Model_Exp_Path_ForAllReal(s_ModelLENS(iLENS),\
                                                            s_Exp, "day")
    s_Real             := l_list[1]
    l_list_            := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelLENS(iLENS), s_Exp, "day")
    s_ExpModel_Path_   := l_list_[0]
    s_Real_            := l_list_[1]
    if (all(ismissing(s_Real)) .or. dimsizes(s_Real_).gt.dimsizes(s_Real)) then
        s_Real         := s_Real_
    end if
    ;
    s_DirRead2 = "/global/cfs/cdirs/m1867/zmchen/my_data/CMIP6/tasmax/POST/TXx/AllReal/"
    s_FileRead          = s_DirRead2 + "TXx_AnMean_" + s_Exp + "_" + \
                        s_ModelLENS(iLENS) + "_CMIP6_AllReal.nc"
    FileTXxProj        := addfile(s_FileRead, "r")
    s_FileRead          = s_DirRead2 + "TXx_AnMean_historical_" + \
                            s_ModelLENS(iLENS) + "_CMIP6_AllReal.nc"
    FileTXxHist        := addfile(s_FileRead, "r")
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2      = "TXx_" + s_ModelLENS(iLENS) + "_" + s_Exp + "_" + \
                        s_Real(iReal)
        if (my_VarNameExist(FileTXxProj, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_Temp         := FileTXxProj->$s_VarName2$
        r_Temp         := my_3TermsMean(r_Temp, i_yrProj)
        ;
        s_VarName2      = "TXx_" + s_ModelLENS(iLENS) + "_historical_" + \
                        s_Real(iReal)
        r_TempHist     := FileTXxHist->$s_VarName2$({1995:2014}, :, :)
        r_TempHist     := dim_avg_n_Wrap(r_TempHist, 0)
        r_TempHist     := conform(r_Temp, r_TempHist, (/1, 2/))
        copy_VarCoords(r_Temp, r_TempHist)
        ;
        r_Temp          = r_Temp - r_TempHist
        ;
        if (isdefined("r_ProjAll") .eq. False .or. iReal .eq. 0) then
            i_dims     := array_append_record(dimsizes(s_Real), \
                                              dimsizes(r_Temp), 0)
            r_ProjAll  := new(i_dims, float)
            copy_VarCoords(r_Temp, r_ProjAll(0, :, :, :))
            r_ProjAll!0 = "Real"
            r_ProjAll&Real = s_Real
        end if 
        r_ProjAll(iReal, :, :, :) = r_Temp + 0
    end do
    ; r_ProjAll          := my_SelChina([/r_ProjAll/])
    ; r_ProjAll          := my_AreaMean2([/r_ProjAll/], r_ProjAll&lat, \
    ;                                      r_ProjAll&lon, 0)
    ;
    ;; subsampling
    do iModel = 0, dimsizes(s_ModelSelected)-1  
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp, "day")
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp, "day")
        s_Real_        := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
            s_Real     := s_Real_
        end if
        ;
        i_dimReal      := dimsizes(s_Real)
        ;
        ; =========================================================
        ;   3. for each model, get the number of ensembel members (n_mem), and then 
        ; subsampling the n_mem members of each LENS and repeating 1000 times 
        rand1_tst      := toint(systemfunc(" date +%s"))
        rand2_tst      := toint((54321l * rand1_tst) % 2147483398l) ;+ iModel
        random_setallseed(rand1_tst, rand2_tst)
        ;
        i_N            := (/i_SubSamp, i_dimReal/)
        i_IndexRandom  := toint(random_uniform(0, dimsizes(r_LENS(:, 0, 0, 0)),\
                                               i_N))
        ;
        i_N            := array_append_record(i_SubSamp, \
                                              dimsizes(r_LENS(0, :, :, :)), 0)
        r_PDF          := new(i_N, float)
        copy_VarCoords(r_LENS(0, :, :, :), r_PDF(0, :, :, :))
        r_PDF_Proj     := r_PDF;(:, :, 0, 0)
        do ii = 0, i_SubSamp-1
            r_Temp     := r_LENS(i_IndexRandom(ii, :), :, :, :)
            r_TempProj := r_ProjAll(i_IndexRandom(ii, :), :, :, :)
            ; printVarSummary(r_Temp)
            ; printVarSummary(r_PDF)
            ; ;
            if (dimsizes(i_IndexRandom(ii, :)) .ge. 2) then
                r_PDF(ii, :, :, :) = dim_avg_n_Wrap(r_Temp, 0)
                r_PDF_Proj(ii, :, :, :)  = dim_avg_n_Wrap(r_TempProj, 0)
            else 
                r_PDF(ii, :, :, :) = r_Temp
                r_PDF_Proj(ii, :, :, :)  = r_TempProj
            end if
        end do
        ;
        ; ==========================================================
        ;   4. estimating the internal var by using the STD of 1000 sythentic members 
        ;; STD of PDF 
        if (iModel .eq. 0 .or. isdefined("r_STD") .eq. False) then
            i_dims     := array_append_record(dimsizes(s_ModelSelected), \
                                              dimsizes(r_PDF(0, :, :, :)), 0)
            r_STD      := new(i_dims, float)
            copy_VarCoords(r_PDF(0, :, :, :), r_STD(0, :, :, :))
            r_STD!0     = "Models"
            r_STD&Models= s_ModelSelected
            ;
            r_STD_Proj := r_STD;(:, :, 0, 0)
        end if
        r_STD(iModel, :, :, :) = dim_stddev_n_Wrap(r_PDF, 0)
        r_STD_Proj(iModel, :, :, :)  = dim_stddev_n_Wrap(r_PDF_Proj, 0)
    end do
    delete([/r_ProjAll/])
    ;
    ; =================================================
    ;   5. make inter-model ensemble of internal var and output 
    r_STD              := dim_avg_n_Wrap(r_STD, 0)
    r_STD_Proj         := dim_avg_n_Wrap(r_STD_Proj, 0)
    ;
    s_VarName2         := "InterVar_STD_ByLENS_" + s_ModelLENS(iLENS) + "_" + \
                        s_Exp
    FileOut->$s_VarName2$ = r_STD
    s_VarName2         := "InterVar_STD_Proj_ByLENS_" + s_ModelLENS(iLENS) + \
                        "_" + s_Exp
    FileOut->$s_VarName2$ = r_STD_Proj
    ;
    if (iLENS .eq. 0 .or. isdefined("r_STD_All") .eq. False) then
        i_dims         := array_append_record(dimsizes(s_ModelLENS), \
                                              dimsizes(r_STD), 0)
        r_STD_All      := new(i_dims, float)
        copy_VarCoords(r_STD, r_STD_All(0, :, :, :))
        r_STD_All!0     = "Models"
        r_STD_All&Models= s_ModelLENS
        ;
        r_STD_ProjAll   = r_STD_All;(:, :, 0, 0)
    end if
    r_STD_All(iLENS, :, :, :) = r_STD
    r_STD_ProjAll(iLENS, :, :, :)   = r_STD_Proj
end do
;
r_STD_China             = my_SelChina([/r_STD_All/])
r_STD_China             = r_STD_China^2
r_STD_China            := my_AreaMean2([/r_STD_China/], r_STD_China&lat, \
                                       r_STD_China&lon, 0)
r_STD_ProjAll          := my_SelChina([/r_STD_ProjAll/])
r_STD_ProjAll           = r_STD_ProjAll^2
r_STD_ProjAll          := my_AreaMean2([/r_STD_ProjAll/], r_STD_ProjAll&lat, \
                                       r_STD_ProjAll&lon, 0)
do iLENS = 0, dimsizes(s_ModelLENS)-1
    print(s_ModelLENS(iLENS) + "")
    print(r_STD_China(iLENS, :) + "")
    print("Proj IV:")
    print(r_STD_ProjAll(iLENS, :) + "")
    print("")
end do
;
end









