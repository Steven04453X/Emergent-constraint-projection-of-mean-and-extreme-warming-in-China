; ==============================================================
; Author: Ziming Chen 
; This code is for: 
;   - making the distribution of pattern scaling coef which are
; used to constrained waring by local residual trend 
;   - constrained the local warming by using the bootstraped 
; pattern scaling coef and cal the inter-member STD
;   1. indicate the models, exp, projection terms
;   2. read the local residual warming trend, the constrained local warming 
; by using the GSAT, and then read the inter-model regression btw 
; local residual trend and constrained lcoal warming 
;   3. for each bootstraped member of constrained local warming, 
; make 1000 bootstraped members for the regression btw local constrained 
; warming and present-day residual trend (attention: the bootstrape only is 
; done over the region with siginicant inter-model regression btw 
; local residual trend and constrained lcoal warming )
;   4. and then constrain the local warming for each bootstraped member 
; ==============================================================
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$WACCEM/Work/Others/Special_Code//my_FunctionInPNNL.ncl"
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
undef("my_FurtherConstrainedBootstrap")
function my_FurtherConstrainedBootstrap(FileYc1, rcBootAll[*]:numeric, \
                                ReTrendAllModel[*]:numeric, TrendObsMean:numeric,\
                                iTerm, ilat, ilon)
local r_rcBootAll, r_ReTrendAllModel, r_TrendObsMean
begin
r_rcBootAll             = rcBootAll
r_ReTrendAllModel       = ReTrendAllModel
r_TrendObsMean          = TrendObsMean
;
s_ModelSelected         = r_ReTrendAllModel&Models
;
r_ReTrendAllModel       = r_ReTrendAllModel - r_TrendObsMean
;
i_boot                  = FileYc1->boot
do iboot = 0, dimsizes(i_boot)-1    
    do iModel = 0, dimsizes(s_ModelSelected)-1
        s_VarName2      = "ConstrainedWarming_Yc1_" + s_ModelSelected(iModel) +\
                        "_In_pc_boot" + iboot
        r_Yc1          := FileYc1->$s_VarName2$(:, iTerm, ilat, ilon)
        ;
        do iboot2 = 0, dimsizes(r_Yc1(:))-1
            r_Yc2i     := r_Yc1(iboot2) - r_rcBootAll *r_ReTrendAllModel(iModel)
            if (iboot2 .eq. 0 .or. isdefined("r_Yc2i_AllBoot2") .eq. False) then
                i_dims := (/dimsizes(r_Yc1), dimsizes(r_Yc2i) /)
                r_Yc2i_AllBoot2 := new(i_dims, float)
                ; copy_VarCoords(r_Yc2i, r_Yc2i_AllBoot2(0, :))
                r_Yc2i_AllBoot2!0= "Boot1"
                r_Yc2i_AllBoot2&Boot1 = ispan(0, dimsizes(r_Yc1)-1, 1)
                r_Yc2i_AllBoot2!1= "Boot2"
                r_Yc2i_AllBoot2&Boot2 = ispan(0, dimsizes(r_Yc1)-1, 1)
            end if
            r_Yc2i_AllBoot2(iboot2, :)= r_Yc2i + 0
        end do
        ;
        if (isdefined("r_Yc2_AllModels") .eq. False .or. iModel .eq. 0) then
            i_dims     := array_append_record(dimsizes(s_ModelSelected), \
                                              dimsizes(r_Yc2i_AllBoot2), 0)
            r_Yc2_AllModels     := new(i_dims, float)
            copy_VarCoords(r_Yc2i_AllBoot2, r_Yc2_AllModels(0, :, :))
            r_Yc2_AllModels!0    = "Models"
            r_Yc2_AllModels&Models = s_ModelSelected
        end if
        r_Yc2_AllModels(iModel,:,:)= r_Yc2i_AllBoot2
        ;
        delete([/r_Yc2i_AllBoot2/])
    end do 
    ;
    ;; Multi-Model Ensemble
    r_Yc2_AllModels    := dim_avg_n_Wrap(r_Yc2_AllModels, 0)
    ;
    if (isdefined("r_Yc2_Boot0t02") .eq. False .or. iboot .eq. 0) then
        i_dims         := array_append_record(dimsizes(i_boot), \
                                              dimsizes(r_Yc2_AllModels), 0)
        r_Yc2_Boot0t02  = new(i_dims, float)
        copy_VarCoords(r_Yc2_AllModels, r_Yc2_Boot0t02(0, :, :))
        r_Yc2_Boot0t02!0= "Boot0"
        r_Yc2_Boot0t02&Boot0 = ispan(0, dimsizes(i_boot)-1, 1)
    end if
    r_Yc2_Boot0t02(iboot, :, :) = r_Yc2_AllModels
    delete([/r_Yc2_AllModels/])
end do
;
return(r_Yc2_Boot0t02)
;
end
;
undef("my_ConstrainedByGSAT_Bootstrap")
function my_ConstrainedByGSAT_Bootstrap(FileYc1, Model[*]:string, iTerm, ilat,\
                                        ilon)
local r_Yc1, s_VarName2, s_ModelSelected
begin
;
s_ModelSelected         = Model
i_boot                  = FileYc1->boot
do iboot = 0, dimsizes(i_boot)-1    
    do iModel = 0, dimsizes(s_ModelSelected)-1
        s_VarName2      = "ConstrainedWarming_Yc1_" + s_ModelSelected(iModel) +\
                        "_In_pc_boot" + iboot
        r_Yc1          := FileYc1->$s_VarName2$(:, iTerm, ilat, ilon)
        ;
        if (iModel .eq. 0 .or. isdefined("r_Yc1_AllModels") .eq. False) then
            i_dims     := array_append_record(dimsizes(s_ModelSelected), \
                                              dimsizes(r_Yc1), 0)
            r_Yc1_AllModels := new(i_dims, float)
            copy_VarCoords(r_Yc1, r_Yc1_AllModels(0, :))
            r_Yc1_AllModels!0= "Models"
            r_Yc1_AllModels&Models = s_ModelSelected
        end if
        r_Yc1_AllModels(iModel, :) = r_Yc1
    end do 
    r_Yc1_AllModels    := dim_avg_n_Wrap(r_Yc1_AllModels, 0)
    ;
    if (isdefined("r_Yc2_Boot0t02") .eq. False .or. iboot .eq. 0) then
        i_dims         := array_append_record(dimsizes(i_boot), \
                                              dimsizes(r_Yc1_AllModels), 0)
        r_Yc2_Boot0t02  = new(i_dims, float)
        copy_VarCoords(r_Yc1_AllModels, r_Yc2_Boot0t02(0, :))
        r_Yc2_Boot0t02!0= "Boot0"
        r_Yc2_Boot0t02&Boot0 = ispan(0, dimsizes(i_boot)-1, 1)
    end if
    r_Yc2_Boot0t02(iboot, :) = r_Yc1_AllModels
    delete([/r_Yc1_AllModels/])
end do
return(r_Yc2_Boot0t02)
;
end
;
begin
;
; ==============================================================
;   1. indicate the models, exp, projection terms
s_ModelSelected         = (/"ACCESS-CM2", "ACCESS-ESM1-5", "AWI-CM-1-1-MR", \
                            "BCC-CSM2-MR", "CAMS-CSM1-0", "CNRM-CM6-1", \
                            "CNRM-ESM2-1", "CanESM5", "CanESM5-CanOE",\
                            "EC-Earth3", "EC-Earth3-Veg", "FGOALS-f3-L", \
                            "FGOALS-g3", "GFDL-CM4", "GFDL-ESM4", \
                            "GISS-E2-1-G", "HadGEM3-GC31-LL", "INM-CM4-8",\
                            "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-LR", \
                            "MRI-ESM2-0", "UKESM1-0-LL"/)
s_4Models               = "_23Models"
if (isdefined("s_Exp") .eq. False) then
    s_Exp               = "ssp245"
end if
i_yrConstrained         = (/1981, 2014/)
i_yrProj                = (/ (/2021, 2040/), (/2041, 2060/), (/2080, 2099/) /)
s_yrProj                = (/"near", "mid", "long"/)
s_VarName               = "tas"
;
; ==============================================================
;   2. read the local residual warming trend, the constrained local warming 
; by using the GSAT, and then read the inter-model regression btw 
; local residual trend and constrained lcoal warming 
s_DirRead = "$WACCEM/my_data/CMIP6/tas/POST/GSAT_AndRegionalTrend_AllReal/"
s_FileRead              = s_DirRead + "GSAT_AndRegionalTrend_1981-2014" + \
                            s_4Models + ".nc"
FileResTrend            = addfile(s_FileRead, "r")
;
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/07.ConstrainAllReal/03.ConstraintByUsingLocalTrend/"
s_FileRead              = s_DirRead + "Regression_LocalTrendAndProjectedWarming_" \
                        + s_Exp + s_4Models + ".nc"
FileRegForResTrendAndCon= addfile(s_FileRead, "r")
;; Regression 
s_VarName2              = "RegInModelDim_" + s_Exp
r_RegResiTrend          = FileRegForResTrendAndCon->$s_VarName2$ ; Reg_SigLev x Terms x lat x lon
r_RegResiTrend         := my_SelChina([/r_RegResiTrend/])
;
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/04.GSAT_scaled_RegionalWarming/"
s_FileRead              = s_DirRead + \
                                "ConstrainedWarmingInChina_RegMMEForScaling" + \
                                s_4Models + ".nc"
FileConLocalWarming     = addfile(s_FileRead, "r")
;
do iModel = 0, dimsizes(s_ModelSelected)-1
    ;; constrained local warming 
    s_VarName2          = "ConstrainedWarmingChina_" + s_ModelSelected(iModel) + \
                            s_Exp
    if (my_VarNameExist(FileConLocalWarming, s_VarName2) .eq. False) then
        print(s_VarName2 + "")
        continue
    end if
    r_ConLocWarming    := FileConLocalWarming->$s_VarName2$ ; Terms x lat x lon
    ;
    ;; local residual warming trend 
    l_list                     := my_CMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp, "Amon")
    s_Real                     := l_list[1]
    l_list_                    := my_RawCMIP6Model_Exp_Path_ForAllReal(\
                                s_ModelSelected(iModel), s_Exp, "Amon")
    s_Real_                    := l_list_[1]
    if (all(ismissing(s_Real)) .or. \
        dimsizes(s_Real_) .gt. dimsizes(s_Real) ) then
        s_Real             := s_Real_
    end if
    ;
    do iReal = 0, dimsizes(s_Real)-1
        s_VarName2              = "ResidualTrend_" + s_ModelSelected(iModel) + \
                                "_historical_" + s_Real(iReal)
        if (my_VarNameExist(FileResTrend, s_VarName2) .eq. False) then
            print(s_VarName2 + " Missing!")
            continue
        end if
        r_Temp                 := FileResTrend->$s_VarName2$(1, :, :)
        ;
        if (isdefined("r_ReTrendAllReal") .eq. False .or. iReal .eq. 0) then
            i_dims             := array_append_record(dimsizes(s_Real), \
                                                      dimsizes(r_Temp), 0)
            r_ReTrendAllReal   := new(i_dims, float)
            copy_VarCoords(r_Temp, r_ReTrendAllReal(0, :, :))
            r_ReTrendAllReal!0  = "Real"
            r_ReTrendAllReal&Real = s_Real
        end if
        r_ReTrendAllReal(iReal, :, :) = r_Temp
    end do
    ;
    r_ReTrendAllReal           := dim_avg_n_Wrap(r_ReTrendAllReal, 0)
    ;
    if (isdefined("r_ReTrendAllModel") .eq. False .or. iModel .eq. 0) then
        i_dims                 := array_append_record(dimsizes(s_ModelSelected), \
                                                      dimsizes(r_ReTrendAllReal), 0)
        r_ReTrendAllModel      := new(i_dims, float)
        copy_VarCoords(r_ReTrendAllReal, r_ReTrendAllModel(0, :, :))
        r_ReTrendAllModel!0     = "Models"
        r_ReTrendAllModel&Models= s_ModelSelected
        ;
        i_dims                 := array_append_record(dimsizes(s_ModelSelected), \
                                                      dimsizes(r_ConLocWarming), 0)
        r_ConLocWarmingAllModel = new(i_dims, float)
        copy_VarCoords(r_ConLocWarming, r_ConLocWarmingAllModel(0, :, :, :))
        copy_VarCoords(r_ReTrendAllModel(:, 0, 0), \
                       r_ConLocWarmingAllModel(:, 0, 0, 0))
    end if
    r_ReTrendAllModel(iModel, :, :) = r_ReTrendAllReal
    r_ConLocWarmingAllModel(iModel, :, :, :) = r_ConLocWarming
end do
;
s_DirOut2 = "/global/cfs/cdirs/m1867/zmchen/my_data/Obs/TEMP/POST/YearMon/"
s_FileRead                      = s_DirOut2 + \
                "SurfaceAirTemp_Post_YearMon_1900-2014.nc"
FileObs                         = addfile(s_FileRead, "r")

s_VarNameObs                    = (/ "Cowtan_WayV2", "BEST", \
                                     "GISTEMP", "NOAAGlobalTemp" /)
;
opt                             = True
opt@DstGridLon                  = r_ReTrendAllModel&lon
opt@DstGridLat                  = r_ReTrendAllModel&lat
opt@DstReginal                  = True
opt@InterpMethod                = "conserve"
opt@ForceOverwrite              = True
;
do iObs = 0, dimsizes(s_VarNameObs)-1
    s_VarName2                  = "tas_" + s_VarNameObs(iObs)
    r_HadCRUT4                 := FileObs->$s_VarName2$(\
                                {i_yrConstrained(0):i_yrConstrained(1)}, :,:,:)
    r_HadCRUT4                  = where(abs(r_HadCRUT4) .lt. 100, r_HadCRUT4, \
                                        r_HadCRUT4@_FillValue)
    ;; Interpolation 
    opt@SrcGridLon             := r_HadCRUT4&lon
    opt@SrcGridLat             := r_HadCRUT4&lat
    system("rm -rf ./destination_grid_file.nc ./PET0.RegridWeightGen.Log ./source_grid_file.nc ./weights_file.nc")
    r_HadCRUT4                 := ESMF_regrid(r_HadCRUT4, opt)
    system("rm -rf ./destination_grid_file.nc ./PET0.RegridWeightGen.Log ./source_grid_file.nc ./weights_file.nc")
    ;
    r_HadCRUT4                 := dim_avg_n_Wrap(r_HadCRUT4, 1)
    r_HadCRUT4!1                = "lat"
    r_HadCRUT4!2                = "lon"
    r_HadCRUT4_Glo             := my_AreaMean2([/r_HadCRUT4/], \
                                            r_HadCRUT4&lat, r_HadCRUT4&lon, 0)
    r_HadCRUT4_Glo             := conform(r_HadCRUT4, r_HadCRUT4_Glo, (/0/))
    copy_VarCoords(r_HadCRUT4, r_HadCRUT4_Glo)
    r_HadCRUT4                  = r_HadCRUT4 - r_HadCRUT4_Glo
    ;
    r_Temp_Trend                = trend_manken(r_HadCRUT4, False, 0)
    copy_VarCoords(r_HadCRUT4(0, :, :), r_Temp_Trend(0, :, :))
    ;
    r_Temp_Trend                = where(abs(r_Temp_Trend) .lt. 100, r_Temp_Trend, \
                                        r_Temp_Trend@_FillValue)
    if (isdefined("r_TrendObs") .eq. False .or. iObs .eq. 0) then
        i_dims                 := array_append_record(dimsizes(s_VarNameObs), \
                                            dimsizes(r_Temp_Trend(0, :, :)), 0)
        r_TrendObs              = new(i_dims, float)
        copy_VarCoords(r_Temp_Trend(0, :, :), r_TrendObs(0, :, :))
        r_TrendObs!0            = "ObsData"
        r_TrendObs&ObsData      = ispan(1, dimsizes(s_VarNameObs), 1)
        r_TrendObs@ObsData      = s_VarNameObs
    end if
    r_TrendObs(iObs, :, :)      = r_Temp_Trend(1, :, :) * 10.
end do 
printMinMax(r_TrendObs, False)
r_TrendObsMean                  = dim_avg_n_Wrap(r_TrendObs, 0)
; printVarSummary(r_TrendObsMean)
; exit()
;
;; China region 
r_ReTrendAllModel               := my_SelChina([/r_ReTrendAllModel/])
r_ConLocWarmingAllModel         := my_SelChina([/r_ConLocWarmingAllModel/])
r_TrendObsMean                  := my_SelChina([/r_TrendObsMean/])
; ==============================================================
;   3. for each bootstraped member of constrained local warming, 
; make 1000 bootstraped members for the regression btw local constrained 
; warming and present-day residual trend (attention: the bootstrape only is 
; done over the region with siginicant inter-model regression btw 
; local residual trend and constrained lcoal warming )
s_DirRead = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/"
s_FileRead                       = s_DirRead + \
                "Dis_of_constrainedLocalWarming_Yc1_" + s_Exp + s_4Models +".nc"
FileYc1                          = addfile(s_FileRead, "r")
;
nBoot                            = 100
opt                              = True
opt@rseed3                       = "clock"
;; bootstraped between residual trend and constrained warming (Yc1)
; printVarSummary(r_RegResiTrend) ; [Reg_SigLev | 2] x [Terms | 3] x [lat | 16] x [lon | 28]
; printVarSummary(r_ReTrendAllModel) ; [Model] x [lat] x [lon]
; printVarSummary(r_ConLocWarmingAllModel) ; [Model] x [Terms] x [lat] x [lon]
i_dimsReg                       := dimsizes(r_RegResiTrend)
do iTerm = 0, i_dimsReg(1)-1
    do ilat = 0, i_dimsReg(2)-1
        do ilon = 0, i_dimsReg(3)-1
            if (ismissing(r_RegResiTrend(0, iTerm, ilat, ilon)) .or. \
                r_RegResiTrend(1, iTerm, ilat, ilon) .gt. 0.1) then
                continue
            end if
            ;
            BootStrap           := bootstrap_regcoef(\
                                r_ReTrendAllModel(:, ilat, ilon), \
                                r_ConLocWarmingAllModel(:, iTerm, ilat, ilon),\
                                nBoot, 0, opt)
            rcBoot              := BootStrap[0] ; Bootstrapped regression coefficients in ascending order
            ;
            if (isdefined("r_rcBootAll") .eq. False) then
                i_dims          := array_append_record(nBoot, \
                            dimsizes(r_ConLocWarmingAllModel(0, :, :, :)), 0)
                r_rcBootAll      = new(i_dims, float)
                copy_VarCoords(rcBoot, r_rcBootAll(:, 0, 0, 0))
                copy_VarCoords(r_ConLocWarmingAllModel(0, :, :, :), \
                               r_rcBootAll(0, :, :, :))
            end if
            r_rcBootAll(:, iTerm, ilat, ilon) = rcBoot
        end do
        ;
    end do
    ;
end do
; printVarSummary(r_rcBootAll) ; [boot | 100] x [Terms | 3] x [lat | 16] x [lon | 28]
; printMinMax(r_rcBootAll, False)
; exit()
;
; ==============================================================
;   4. and then constrain the local warming for each bootstraped member 
do iTerm = 0, i_dimsReg(1)-1
    do ilat = 0, i_dimsReg(2)-1
        do ilon = 0, i_dimsReg(3)-1
            i_Input                 := 1
            if (ismissing(r_RegResiTrend(0, iTerm, ilat, ilon)) .or. \
                r_RegResiTrend(1, iTerm, ilat, ilon) .gt. 0.1) then
                r_FurtherConstrained := my_ConstrainedByGSAT_Bootstrap(FileYc1,\
                                            s_ModelSelected, iTerm, ilat, ilon)
                r_FurtherConstrained_:= reshape(r_FurtherConstrained, (/10000/))
                i_Input              = 0
            ;
            else 
                r_FurtherConstrained:= my_FurtherConstrainedBootstrap(FileYc1, \
                                r_rcBootAll(:, iTerm, ilat, ilon), \
                                r_ReTrendAllModel(:, ilat, ilon), \
                                r_TrendObsMean(ilat, ilon), iTerm, ilat, ilon)
                printVarSummary(r_FurtherConstrained) ; [Boot0 | 100] x [Boot1 | 100] x [Boot2 | 100]
                printMinMax(r_FurtherConstrained, False)
            ;
                r_FurtherConstrained_ := reshape(r_FurtherConstrained, \
                                                 (/1000000/))
                ; printVarSummary(r_FurtherConstrained_)
                ; printMinMax(r_FurtherConstrained_, False)
            end if
            ;
            if (isdefined("r_FurtherConstrained_AllBoot") .eq. False) then
                i_dims          := array_append_record(1000000, \
                                        dimsizes(r_rcBootAll(0, :, :, :)), 0)
                r_FurtherConstrained_AllBoot          = new(i_dims, float)
                r_FurtherConstrained_AllBoot!0        = "AllBoots"
                r_FurtherConstrained_AllBoot&AllBoots = ispan(0, i_dims(0)-1, 1)
                copy_VarCoords(r_rcBootAll(0, :, :, :), \
                               r_FurtherConstrained_AllBoot(0, :, :, :))
            end if
            if (i_Input .eq. 1) then
                r_FurtherConstrained_AllBoot(:, iTerm, ilat, ilon) = \
                                                        r_FurtherConstrained_ +0
            else 
                r_FurtherConstrained_AllBoot(:10000-1, iTerm, ilat, ilon) = \
                                                        r_FurtherConstrained_ +0
            end if
            ;
        end do 
        ;
    end do 
    ;
end do 
printVarSummary(r_FurtherConstrained_AllBoot)
printMinMax(r_FurtherConstrained_AllBoot, False)
;
r_FurtherConstrained_STD         = dim_stddev_n_Wrap(\
                                                r_FurtherConstrained_AllBoot, 0)
printMinMax(r_FurtherConstrained_STD, False)
;
s_DirOut = "$WACCEM/Work/OtherWorks/2020/ProjectionInChina/Data/07.UncerQualification/01.ForMeanWarming/"
system("mkdir -p " + s_DirOut)
s_FileOut                        = s_DirOut + \
                    "Dis_of_ConstrainedLocalWarmingByLocalTrend_Yc2_" + s_Exp + \
                    s_4Models + ".nc"
system("rm -rf " + s_FileOut)
setfileoption("nc", "Format", "LargeFile")
FileOut                          = addfile(s_FileOut, "c")
GAttr                            = True
GAttr@history                    = systemfunc("date") + "  " + systemfunc("pwd")\
           + "/03.Dis_of_pc2_and_constrainedLocalWarmingByLocalTrend_Yc2.ncl"
fileattdef(FileOut, GAttr)
; s_VarName2                       = "PDF_FurConByLocalTrend_Yc2_" + s_Exp
; FileOut->$s_VarName2$            = r_FurtherConstrained_AllBoot
s_VarName2                       = "STD_FurConByLocalTrend_Yc2_" + s_Exp
FileOut->$s_VarName2$            = r_FurtherConstrained_STD
;
print(s_FileOut + "")
;
end

















