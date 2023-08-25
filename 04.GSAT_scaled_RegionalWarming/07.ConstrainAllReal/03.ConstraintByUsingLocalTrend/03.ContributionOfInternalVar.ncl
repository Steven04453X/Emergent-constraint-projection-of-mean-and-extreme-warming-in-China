; ===================================================================
; Author: Ziming Chen 
; This code is for:
;   - estimating the contribution of the internal var (IV)
;   1. indicate the models, exp, terms, var name
;   2. read the constrained results of each realization in each model 
;   3. estimating the internal variability by computing the inter-member standard deviation:
; CSAT and regional warming
;   4. make multi-model ensemble mean 
;   5. plot the results:
;   IV in near, IV in mid
;   IV in long, IV of CSAT in three terms
; ===================================================================
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
; 注意返回的矩阵就是挑选中国区域的矩阵了
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
        r_VarChina_({r_latTaiwan(0):r_latTiwan(1)}, \
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
; ===================================================================
;   1. indicate the models, exp, terms, var name
; s_ModelSelect           = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
;                             "BCC-CSM2-MR", "CAMS-CSM1-0", "CESM2-WACCM", \
;                             "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5", \
;                             "CanESM5-CanOE", "EC-Earth3", "EC-Earth3-Veg", \
;                             "FGOALS-f3-L", "FGOALS-g3", "GFDL-CM4", \
;                             "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", \
;                             "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", \
;                             "MCM-UA-1-0", "MIROC6", "MIROC-ES2L", \
;                             "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", \
;                             "NESM3", "NorESM2-LM", "UKESM1-0-LL"/)
s_ModelSelect           = (/"ACCESS-CM2", "CNRM-CM6-1", "EC-Earth3", \
                            "GFDL-CM4", "MPI-ESM1-2-LR"/)
s_4Models               = "_5Models"
;
s_Exp                   = (/"ssp126", "ssp245", "ssp370", "ssp585"/)
i_yrProj                = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
;
; ===================================================================
;   2. read the constrained results of each realization in each model
;; Constrained by regional warming 
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
s_FileRead              = s_DirRead + \
            "PD_LocalTrend_ConstrainLocalWarming_ForAllReal" + s_4Models + ".nc"
FileRegionalConstrained = addfile(s_FileRead, "r")
;
;; Constrained by GSAT 
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/"
s_FileRead              = s_DirRead + "Constrained_tas_InEachReal" + s_4Models \
                        + ".nc"
FileGSAT_Constrained    = addfile(s_FileRead, "r")
;
do iExp = 0, dimsizes(s_Exp)-1
    do iModel = 0, dimsizes(s_ModelSelect)-1
        l_list         := my_CMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        s_Real         := l_list[1]
        l_list_        := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                    s_ModelSelect(iModel), s_Exp(iExp), "Amon")
        s_Real_        := l_list_[1]
        if (all(ismissing(s_Real)) .or. \
            dimsizes(s_Real_) .gt. dimsizes(s_Real)) then
            s_Real     := s_Real_
        end if
        ;
        if (dimsizes(s_Real) .gt. 1) then
            do iReal = 0, dimsizes(s_Real)-1
                ;; Constrained By GSAT 
                s_VarName2 := "Constrained_tas_" + s_ModelSelect(iModel) + "_"+\
                        s_Exp(iExp) + "_" + s_Real(iReal)
                if (my_VarNameExist(FileGSAT_Constrained, s_VarName2) .eq. \
                    False) then
                    print(s_VarName2 + " Missing!")
                    continue
                end if
                r_GSAT_Constrained := FileGSAT_Constrained->$s_VarName2$
                r_CSAT_ConstrainedByGSAT := my_SelChina([/r_GSAT_Constrained/])
                ; s_FileTest = "./Test.nc"
                ; system("rm -rf " + s_FileTest)
                ; FileTest   = addfile(s_FileTest, "c")
                ; FileTest->Test = r_CSAT_ConstrainedByGSAT
                ; exit()
                r_CSAT_ConstrainedByGSAT := my_AreaMean2(\
                                            [/r_CSAT_ConstrainedByGSAT/],\
                                            r_CSAT_ConstrainedByGSAT&lat, \
                                            r_CSAT_ConstrainedByGSAT&lon, 0)
            ;
                ;; Constrained by regional warming trend 
                s_VarName2 := "ConstrainedRegionalWarming_" + s_ModelSelect(iModel)\
                        + "_" + s_Exp(iExp) + "_" + s_Real(iReal)
                if (my_VarNameExist(FileRegionalConstrained, s_VarName2) .eq. \
                    False) then
                    print(s_VarName2 + " Missing!")
                    continue
                end if
                r_RegionalConstrained := FileRegionalConstrained->$s_VarName2$
                r_CSAT_ConstrainedByRegionalWarming := \
                                        my_SelChina([/r_RegionalConstrained/])
                r_CSAT_ConstrainedByRegionalWarming := my_AreaMean2(\
                                    [/r_CSAT_ConstrainedByRegionalWarming/], \
                                    r_CSAT_ConstrainedByRegionalWarming&lat, \
                                    r_CSAT_ConstrainedByRegionalWarming&lon,0)
                ;
                if (isdefined("r_CSAT_ConstrainedByGSAT_AllReal") .eq. False \
                    .or. iReal .eq. 0) then
                    r_CSAT_ConstrainedByGSAT_AllReal := new(\
                                                (/dimsizes(s_Real),3/), float)
                    r_CSAT_ConstrainedByGSAT_AllReal!0= "Real"
                    r_CSAT_ConstrainedByGSAT_AllReal&Real = s_Real
                    ;
                    r_CSAT_ConstrainedByRegionalWarming_AllReal:= \
                                                r_CSAT_ConstrainedByGSAT_AllReal
                    ;
                    i_dims                           := array_append_record(\
                            dimsizes(s_Real), dimsizes(r_GSAT_Constrained), 0)
                    r_GSAT_Constrained_AllReal       := new(i_dims, float)
                    copy_VarCoords(r_GSAT_Constrained, \
                                   r_GSAT_Constrained_AllReal(0, :, :, :))
                    ;
                    i_dims                           := array_append_record(\
                            dimsizes(s_Real), dimsizes(r_RegionalConstrained), 0)
                    r_RegionalConstrained_AllReal    := new(i_dims, float)
                    copy_VarCoords(r_RegionalConstrained, \
                                   r_RegionalConstrained_AllReal(0, :, :, :))
                end if
                r_CSAT_ConstrainedByGSAT_AllReal(iReal, :)    = \
                                                r_CSAT_ConstrainedByGSAT + 0
                r_CSAT_ConstrainedByRegionalWarming_AllReal(iReal, :) = \
                                        r_CSAT_ConstrainedByRegionalWarming + 0
                r_GSAT_Constrained_AllReal(iReal, :, :, :)    = \
                                                        r_GSAT_Constrained + 0
                r_RegionalConstrained_AllReal(iReal, :, :, :) = \
                                                    r_RegionalConstrained + 0
            end do 
            ;
            ; ================================================================
            ;   3. estimating the internal variability by computing the inter-member standard deviation:
            ; CSAT and regional warming
            ; print(stddev(r_CSAT_ConstrainedByGSAT_AllReal(:, 2)) + ", " + \
            ;       stddev(r_CSAT_ConstrainedByRegionalWarming_AllReal(:, 2)))
            ; exit()
            r_CSAT_ConByGSAT_STD         := dim_stddev_n_Wrap(\
                                            r_CSAT_ConstrainedByGSAT_AllReal, 0)
            r_CSAT_ConByRegionalWarm_STD := dim_stddev_n_Wrap(\
                                r_CSAT_ConstrainedByRegionalWarming_AllReal, 0)
            r_PatternConByGSAT_STD       := dim_stddev_n_Wrap(\
                                            r_GSAT_Constrained_AllReal, 0)
            r_PatternConByRegionalWarm_STD:= dim_stddev_n_Wrap(\
                                            r_RegionalConstrained_AllReal, 0)
        end if
        ;
        if (isdefined("r_CSAT_ConByGSAT_STD") .eq. True) then
            ;; store up the STD for all models 
            if (isdefined("r_CSAT_ConByGSAT_STD_AllModels") .eq. False .or. \
                iModel .eq. 0) then
                i_dims                   := array_append_record(\
                                            dimsizes(s_ModelSelect), \
                                            dimsizes(r_CSAT_ConByGSAT_STD), 0)
                r_CSAT_ConByGSAT_STD_AllModels := new(i_dims, float)
                copy_VarCoords(r_CSAT_ConByGSAT_STD, \
                               r_CSAT_ConByGSAT_STD_AllModels(0, :))
                r_CSAT_ConByGSAT_STD_AllModels!0= "Models"
                r_CSAT_ConByGSAT_STD_AllModels&Models   = s_ModelSelect
                ;
                r_CSAT_ConByRegionalWarm_STD_AllModels := \
                                r_CSAT_ConByGSAT_STD_AllModels
                ;
                i_dims                    := array_append_record(\
                                            dimsizes(s_ModelSelect), \
                                            dimsizes(r_PatternConByGSAT_STD), 0)
                r_PatternConByGSAT_STD_AllModels       := new(i_dims, float)
                copy_VarCoords(r_PatternConByGSAT_STD, \
                               r_PatternConByGSAT_STD_AllModels(0, :, :, :))
                copy_VarCoords(r_CSAT_ConByGSAT_STD_AllModels(:, 0), \
                               r_PatternConByGSAT_STD_AllModels(:, 0, 0, 0))
                ;
                i_dims                    := array_append_record(\
                                    dimsizes(s_ModelSelect), \
                                    dimsizes(r_PatternConByRegionalWarm_STD), 0)
                r_PatternConByRegionalWarm_STD_AllModels:= new(i_dims, float)
                copy_VarCoords(r_PatternConByRegionalWarm_STD, \
                               r_PatternConByRegionalWarm_STD_AllModels(0,:,:,:))
                copy_VarCoords(r_CSAT_ConByGSAT_STD_AllModels(:, 0), \
                               r_PatternConByRegionalWarm_STD_AllModels(:,0,0,0))
            end if
            ;
            r_CSAT_ConByGSAT_STD_AllModels(iModel, :) = r_CSAT_ConByGSAT_STD
            r_CSAT_ConByRegionalWarm_STD_AllModels(iModel, :) = \
                                                    r_CSAT_ConByRegionalWarm_STD
            r_PatternConByGSAT_STD_AllModels(iModel, :, :, :) = \
                                                    r_PatternConByGSAT_STD
            r_PatternConByRegionalWarm_STD_AllModels(iModel, :, :, :) = \
                                                r_PatternConByRegionalWarm_STD
            ;
            delete([/r_CSAT_ConstrainedByGSAT_AllReal, \
                 r_CSAT_ConstrainedByRegionalWarming_AllReal, \
                 r_GSAT_Constrained_AllReal, r_RegionalConstrained_AllReal, \
                 r_CSAT_ConByGSAT_STD, r_CSAT_ConByRegionalWarm_STD, \
                 r_PatternConByGSAT_STD, r_PatternConByRegionalWarm_STD/])
        end if
        ;
    end do
    ;
    ; =========================================================
    ;   4. make multi-model ensemble mean 
    r_CSAT_ConByGSAT_STD_AllModels  := dim_avg_n_Wrap(\
                                            r_CSAT_ConByGSAT_STD_AllModels, 0)
    r_CSAT_ConByRegionalWarm_STD_AllModels   := dim_avg_n_Wrap(\
                                    r_CSAT_ConByRegionalWarm_STD_AllModels, 0)
    r_PatternConByGSAT_STD_AllModels:= dim_avg_n_Wrap(\
                                    r_PatternConByGSAT_STD_AllModels, 0)
    r_PatternConByRegionalWarm_STD_AllModels := dim_avg_n_Wrap(\
                                    r_PatternConByRegionalWarm_STD_AllModels, 0)
    ;
    r_AreaMean_ConByGSAT_STD        := my_SelChina([/r_PatternConByGSAT_STD_AllModels/])
    r_AreaMean_ConByGSAT_STD        := my_AreaMean2(\
                                    [/r_AreaMean_ConByGSAT_STD/], \
                                    r_AreaMean_ConByGSAT_STD&lat, \
                                    r_AreaMean_ConByGSAT_STD&lon, 0)
    ;
    r_AreaMean_ConByRegionalWarm_STD:= my_SelChina([/r_PatternConByRegionalWarm_STD_AllModels/])
    r_AreaMean_ConByRegionalWarm_STD:= my_AreaMean2(\
                                    [/r_AreaMean_ConByRegionalWarm_STD/], \
                                    r_AreaMean_ConByRegionalWarm_STD&lat, \
                                    r_AreaMean_ConByRegionalWarm_STD&lon, 0)
    ;
    print(s_Exp(iExp) + "")
    print(r_CSAT_ConByGSAT_STD_AllModels + ", " + \
          r_CSAT_ConByRegionalWarm_STD_AllModels)
    print("STD Area Mean: " + r_AreaMean_ConByGSAT_STD + ", " + \
          r_AreaMean_ConByRegionalWarm_STD)
    print("")
end do
;

;   5. plot the results:
;   IV in near, IV in mid
;   IV in long, IV of CSAT in three terms

;
end









