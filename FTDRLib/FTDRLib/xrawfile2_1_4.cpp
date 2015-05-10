// Created by Microsoft (R) C/C++ Compiler Version 14.00.50727.762 (137c9bb3).
//
// xrawfile2.tli
//
// Wrapper implementations for Win32 type library C:\\Xcalibur\\system\\programs\XRawfile2.dll
// compiler-generated file created 01/18/08 at 16:26:22 - DO NOT EDIT!

#pragma once

//
// interface IXRawfile wrapper method implementations
//

inline HRESULT IXRawfile::Open ( _bstr_t szFileName ) {
    HRESULT _hr = raw_Open(szFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::Close ( ) {
    HRESULT _hr = raw_Close();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFileName ( BSTR * pbstrFileName ) {
    HRESULT _hr = raw_GetFileName(pbstrFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetCreatorID ( BSTR * pbstrCreatorID ) {
    HRESULT _hr = raw_GetCreatorID(pbstrCreatorID);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetVersionNumber ( long * pnVersion ) {
    HRESULT _hr = raw_GetVersionNumber(pnVersion);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetCreationDate ( DATE * pCreationDate ) {
    HRESULT _hr = raw_GetCreationDate(pCreationDate);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::IsError ( long * pbIsError ) {
    HRESULT _hr = raw_IsError(pbIsError);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::IsNewFile ( long * pbIsNewFile ) {
    HRESULT _hr = raw_IsNewFile(pbIsNewFile);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetErrorCode ( long * pnErrorCode ) {
    HRESULT _hr = raw_GetErrorCode(pnErrorCode);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetErrorMessage ( BSTR * pbstrErrorMessage ) {
    HRESULT _hr = raw_GetErrorMessage(pbstrErrorMessage);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetWarningMessage ( BSTR * pbstrWarningMessage ) {
    HRESULT _hr = raw_GetWarningMessage(pbstrWarningMessage);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowNumber ( long * pnSeqRowNumber ) {
    HRESULT _hr = raw_GetSeqRowNumber(pnSeqRowNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowSampleType ( long * pnSampleType ) {
    HRESULT _hr = raw_GetSeqRowSampleType(pnSampleType);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowDataPath ( BSTR * pbstrDataPath ) {
    HRESULT _hr = raw_GetSeqRowDataPath(pbstrDataPath);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowRawFileName ( BSTR * pbstrRawFileName ) {
    HRESULT _hr = raw_GetSeqRowRawFileName(pbstrRawFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowSampleName ( BSTR * pbstrSampleName ) {
    HRESULT _hr = raw_GetSeqRowSampleName(pbstrSampleName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowSampleID ( BSTR * pbstrSampleID ) {
    HRESULT _hr = raw_GetSeqRowSampleID(pbstrSampleID);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowComment ( BSTR * pbstrComment ) {
    HRESULT _hr = raw_GetSeqRowComment(pbstrComment);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowLevelName ( BSTR * pbstrLevelName ) {
    HRESULT _hr = raw_GetSeqRowLevelName(pbstrLevelName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowUserText ( long nIndex, BSTR * pbstrUserText ) {
    HRESULT _hr = raw_GetSeqRowUserText(nIndex, pbstrUserText);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowInstrumentMethod ( BSTR * pbstrInstrumentMethod ) {
    HRESULT _hr = raw_GetSeqRowInstrumentMethod(pbstrInstrumentMethod);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowProcessingMethod ( BSTR * pbstrProcessingMethod ) {
    HRESULT _hr = raw_GetSeqRowProcessingMethod(pbstrProcessingMethod);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowCalibrationFile ( BSTR * pbstrCalibrationFile ) {
    HRESULT _hr = raw_GetSeqRowCalibrationFile(pbstrCalibrationFile);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowVial ( BSTR * pbstrVial ) {
    HRESULT _hr = raw_GetSeqRowVial(pbstrVial);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowInjectionVolume ( double * pdInjVol ) {
    HRESULT _hr = raw_GetSeqRowInjectionVolume(pdInjVol);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowSampleWeight ( double * pdSampleWt ) {
    HRESULT _hr = raw_GetSeqRowSampleWeight(pdSampleWt);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowSampleVolume ( double * pdSampleVolume ) {
    HRESULT _hr = raw_GetSeqRowSampleVolume(pdSampleVolume);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowISTDAmount ( double * pdISTDAmount ) {
    HRESULT _hr = raw_GetSeqRowISTDAmount(pdISTDAmount);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowDilutionFactor ( double * pdDilutionFactor ) {
    HRESULT _hr = raw_GetSeqRowDilutionFactor(pdDilutionFactor);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSeqRowUserLabel ( long nIndex, BSTR * pbstrUserLabel ) {
    HRESULT _hr = raw_GetSeqRowUserLabel(nIndex, pbstrUserLabel);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::InAcquisition ( long * pbInAcquisition ) {
    HRESULT _hr = raw_InAcquisition(pbInAcquisition);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumberOfControllers ( long * pnNumControllers ) {
    HRESULT _hr = raw_GetNumberOfControllers(pnNumControllers);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetControllerType ( long nIndex, long * pnControllerType ) {
    HRESULT _hr = raw_GetControllerType(nIndex, pnControllerType);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::SetCurrentController ( long nControllerType, long nControllerNumber ) {
    HRESULT _hr = raw_SetCurrentController(nControllerType, nControllerNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetCurrentController ( long * pnControllerType, long * pnControllerNumber ) {
    HRESULT _hr = raw_GetCurrentController(pnControllerType, pnControllerNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumSpectra ( long * pnNumberOfSpectra ) {
    HRESULT _hr = raw_GetNumSpectra(pnNumberOfSpectra);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumStatusLog ( long * pnNumberOfStatusLogEntries ) {
    HRESULT _hr = raw_GetNumStatusLog(pnNumberOfStatusLogEntries);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumErrorLog ( long * pnNumberOfErrorLogEntries ) {
    HRESULT _hr = raw_GetNumErrorLog(pnNumberOfErrorLogEntries);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumTuneData ( long * pnNumTuneData ) {
    HRESULT _hr = raw_GetNumTuneData(pnNumTuneData);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetMassResolution ( double * pdMassResolution ) {
    HRESULT _hr = raw_GetMassResolution(pdMassResolution);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetExpectedRunTime ( double * pdExpectedRunTime ) {
    HRESULT _hr = raw_GetExpectedRunTime(pdExpectedRunTime);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumTrailerExtra ( long * pnNumberOfTrailerExtraEntries ) {
    HRESULT _hr = raw_GetNumTrailerExtra(pnNumberOfTrailerExtraEntries);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetLowMass ( double * pdLowMass ) {
    HRESULT _hr = raw_GetLowMass(pdLowMass);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetHighMass ( double * pdHighMass ) {
    HRESULT _hr = raw_GetHighMass(pdHighMass);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStartTime ( double * pdStartTime ) {
    HRESULT _hr = raw_GetStartTime(pdStartTime);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetEndTime ( double * pdEndTime ) {
    HRESULT _hr = raw_GetEndTime(pdEndTime);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetMaxIntegratedIntensity ( double * pdMaxIntegIntensity ) {
    HRESULT _hr = raw_GetMaxIntegratedIntensity(pdMaxIntegIntensity);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetMaxIntensity ( long * pnMaxIntensity ) {
    HRESULT _hr = raw_GetMaxIntensity(pnMaxIntensity);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFirstSpectrumNumber ( long * pnFirstSpectrum ) {
    HRESULT _hr = raw_GetFirstSpectrumNumber(pnFirstSpectrum);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetLastSpectrumNumber ( long * pnLastSpectrum ) {
    HRESULT _hr = raw_GetLastSpectrumNumber(pnLastSpectrum);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstrumentID ( long * pnInstrumentID ) {
    HRESULT _hr = raw_GetInstrumentID(pnInstrumentID);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInletID ( long * pnInletID ) {
    HRESULT _hr = raw_GetInletID(pnInletID);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetErrorFlag ( long * pnErrorFlag ) {
    HRESULT _hr = raw_GetErrorFlag(pnErrorFlag);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSampleVolume ( double * pdSampleVolume ) {
    HRESULT _hr = raw_GetSampleVolume(pdSampleVolume);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSampleWeight ( double * pdSampleWeight ) {
    HRESULT _hr = raw_GetSampleWeight(pdSampleWeight);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetVialNumber ( long * pnVialNumber ) {
    HRESULT _hr = raw_GetVialNumber(pnVialNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInjectionVolume ( double * pdInjectionVolume ) {
    HRESULT _hr = raw_GetInjectionVolume(pdInjectionVolume);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFlags ( BSTR * pbstrFlags ) {
    HRESULT _hr = raw_GetFlags(pbstrFlags);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetAcquisitionFileName ( BSTR * pbstrFileName ) {
    HRESULT _hr = raw_GetAcquisitionFileName(pbstrFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstrumentDescription ( BSTR * pbstrInstrumentDescription ) {
    HRESULT _hr = raw_GetInstrumentDescription(pbstrInstrumentDescription);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetAcquisitionDate ( BSTR * pbstrAcquisitionDate ) {
    HRESULT _hr = raw_GetAcquisitionDate(pbstrAcquisitionDate);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetOperator ( BSTR * pbstrOperator ) {
    HRESULT _hr = raw_GetOperator(pbstrOperator);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetComment1 ( BSTR * pbstrComment1 ) {
    HRESULT _hr = raw_GetComment1(pbstrComment1);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetComment2 ( BSTR * pbstrComment2 ) {
    HRESULT _hr = raw_GetComment2(pbstrComment2);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSampleAmountUnits ( BSTR * pbstrSampleAmountUnits ) {
    HRESULT _hr = raw_GetSampleAmountUnits(pbstrSampleAmountUnits);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInjectionAmountUnits ( BSTR * pbstrInjectionAmountUnits ) {
    HRESULT _hr = raw_GetInjectionAmountUnits(pbstrInjectionAmountUnits);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetSampleVolumeUnits ( BSTR * pbstrSampleVolumeUnits ) {
    HRESULT _hr = raw_GetSampleVolumeUnits(pbstrSampleVolumeUnits);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFilters ( VARIANT * pvarFilterArray, long * pnArraySize ) {
    HRESULT _hr = raw_GetFilters(pvarFilterArray, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::ScanNumFromRT ( double dRT, long * pnScanNumber ) {
    HRESULT _hr = raw_ScanNumFromRT(dRT, pnScanNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::RTFromScanNum ( long nScanNumber, double * pdRT ) {
    HRESULT _hr = raw_RTFromScanNum(nScanNumber, pdRT);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFilterForScanNum ( long nScanNumber, BSTR * pbstrFilter ) {
    HRESULT _hr = raw_GetFilterForScanNum(nScanNumber, pbstrFilter);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFilterForScanRT ( double dRT, BSTR * pbstrFilter ) {
    HRESULT _hr = raw_GetFilterForScanRT(dRT, pbstrFilter);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetMassListFromScanNum ( long * pnScanNumber, _bstr_t bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double * pdCentroidPeakWidth, VARIANT * pvarMassList, VARIANT * pvarPeakFlags, long * pnArraySize ) {
    HRESULT _hr = raw_GetMassListFromScanNum(pnScanNumber, bstrFilter, nIntensityCutoffType, nIntensityCutoffValue, nMaxNumberOfPeaks, bCentroidResult, pdCentroidPeakWidth, pvarMassList, pvarPeakFlags, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetMassListFromRT ( double * pdRT, _bstr_t bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double * pdCentroidPeakWidth, VARIANT * pvarMassList, VARIANT * pvarPeakFlags, long * pnArraySize ) {
    HRESULT _hr = raw_GetMassListFromRT(pdRT, bstrFilter, nIntensityCutoffType, nIntensityCutoffValue, nMaxNumberOfPeaks, bCentroidResult, pdCentroidPeakWidth, pvarMassList, pvarPeakFlags, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNextMassListFromScanNum ( long * pnScanNumber, _bstr_t bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double * pdCentroidPeakWidth, VARIANT * pvarMassList, VARIANT * pvarPeakFlags, long * pnArraySize ) {
    HRESULT _hr = raw_GetNextMassListFromScanNum(pnScanNumber, bstrFilter, nIntensityCutoffType, nIntensityCutoffValue, nMaxNumberOfPeaks, bCentroidResult, pdCentroidPeakWidth, pvarMassList, pvarPeakFlags, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetPrevMassListFromScanNum ( long * pnScanNumber, _bstr_t bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double * pdCentroidPeakWidth, VARIANT * pvarMassList, VARIANT * pvarPeakFlags, long * pnArraySize ) {
    HRESULT _hr = raw_GetPrevMassListFromScanNum(pnScanNumber, bstrFilter, nIntensityCutoffType, nIntensityCutoffValue, nMaxNumberOfPeaks, bCentroidResult, pdCentroidPeakWidth, pvarMassList, pvarPeakFlags, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::IsProfileScanForScanNum ( long nScanNumber, long * pbIsProfileScan ) {
    HRESULT _hr = raw_IsProfileScanForScanNum(nScanNumber, pbIsProfileScan);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::IsCentroidScanForScanNum ( long nScanNumber, long * pbIsCentroidScan ) {
    HRESULT _hr = raw_IsCentroidScanForScanNum(nScanNumber, pbIsCentroidScan);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetScanHeaderInfoForScanNum ( long nScanNumber, long * pnNumPackets, double * pdStartTime, double * pdLowMass, double * pdHighMass, double * pdTIC, double * pdBasePeakMass, double * pdBasePeakIntensity, long * pnNumChannels, long * pbUniformTime, double * pdFrequency ) {
    HRESULT _hr = raw_GetScanHeaderInfoForScanNum(nScanNumber, pnNumPackets, pdStartTime, pdLowMass, pdHighMass, pdTIC, pdBasePeakMass, pdBasePeakIntensity, pnNumChannels, pbUniformTime, pdFrequency);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogForScanNum ( long nScanNumber, double * pdStatusLogRT, VARIANT * pvarLabels, VARIANT * pvarValues, long * pnArraySize ) {
    HRESULT _hr = raw_GetStatusLogForScanNum(nScanNumber, pdStatusLogRT, pvarLabels, pvarValues, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogForRT ( double * pdRT, VARIANT * pvarLabels, VARIANT * pvarValues, long * pnArraySize ) {
    HRESULT _hr = raw_GetStatusLogForRT(pdRT, pvarLabels, pvarValues, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogLabelsForScanNum ( long nScanNumber, double * pdStatusLogRT, VARIANT * pvarLabels, long * pnArraySize ) {
    HRESULT _hr = raw_GetStatusLogLabelsForScanNum(nScanNumber, pdStatusLogRT, pvarLabels, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogLabelsForRT ( double * pdRT, VARIANT * pvarLabels, long * pnArraySize ) {
    HRESULT _hr = raw_GetStatusLogLabelsForRT(pdRT, pvarLabels, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogValueForScanNum ( long nScanNumber, _bstr_t bstrLabel, double * pdStatusLogRT, VARIANT * pvarValue ) {
    HRESULT _hr = raw_GetStatusLogValueForScanNum(nScanNumber, bstrLabel, pdStatusLogRT, pvarValue);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogValueForRT ( double * pdRT, _bstr_t bstrLabel, VARIANT * pvarValue ) {
    HRESULT _hr = raw_GetStatusLogValueForRT(pdRT, bstrLabel, pvarValue);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTrailerExtraForScanNum ( long nScanNumber, VARIANT * pvarLabels, VARIANT * pvarValues, long * pnArraySize ) {
    HRESULT _hr = raw_GetTrailerExtraForScanNum(nScanNumber, pvarLabels, pvarValues, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTrailerExtraForRT ( double * pdRT, VARIANT * pvarLabels, VARIANT * pvarValues, long * pnArraySize ) {
    HRESULT _hr = raw_GetTrailerExtraForRT(pdRT, pvarLabels, pvarValues, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTrailerExtraLabelsForScanNum ( long nScanNumber, VARIANT * pvarLabels, long * pnArraySize ) {
    HRESULT _hr = raw_GetTrailerExtraLabelsForScanNum(nScanNumber, pvarLabels, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTrailerExtraLabelsForRT ( double * pdRT, VARIANT * pvarLabels, long * pnArraySize ) {
    HRESULT _hr = raw_GetTrailerExtraLabelsForRT(pdRT, pvarLabels, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTrailerExtraValueForScanNum ( long nScanNumber, _bstr_t bstrLabel, VARIANT * pvarValue ) {
    HRESULT _hr = raw_GetTrailerExtraValueForScanNum(nScanNumber, bstrLabel, pvarValue);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTrailerExtraValueForRT ( double * pdRT, _bstr_t bstrLabel, VARIANT * pvarValue ) {
    HRESULT _hr = raw_GetTrailerExtraValueForRT(pdRT, bstrLabel, pvarValue);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetErrorLogItem ( long nItemNumber, double * pdRT, BSTR * pbstrErrorMessage ) {
    HRESULT _hr = raw_GetErrorLogItem(nItemNumber, pdRT, pbstrErrorMessage);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTuneData ( long nSegmentNumber, VARIANT * pvarLabels, VARIANT * pvarValues, long * pnArraySize ) {
    HRESULT _hr = raw_GetTuneData(nSegmentNumber, pvarLabels, pvarValues, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumInstMethods ( long * pnNumInstMethods ) {
    HRESULT _hr = raw_GetNumInstMethods(pnNumInstMethods);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstMethod ( long nInstMethodItem, BSTR * pbstrInstMethod ) {
    HRESULT _hr = raw_GetInstMethod(nInstMethodItem, pbstrInstMethod);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetChroData ( long nChroType1, long nChroOperator, long nChroType2, _bstr_t bstrFilter, _bstr_t bstrMassRanges1, _bstr_t bstrMassRanges2, double dDelay, double * pdStartTime, double * pdEndTime, long nSmoothingType, long nSmoothingValue, VARIANT * pvarChroData, VARIANT * pvarPeakFlags, long * pnArraySize ) {
    HRESULT _hr = raw_GetChroData(nChroType1, nChroOperator, nChroType2, bstrFilter, bstrMassRanges1, bstrMassRanges2, dDelay, pdStartTime, pdEndTime, nSmoothingType, nSmoothingValue, pvarChroData, pvarPeakFlags, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::RefreshViewOfFile ( ) {
    HRESULT _hr = raw_RefreshViewOfFile();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTuneDataValue ( long nSegmentNumber, _bstr_t bstrLabel, VARIANT * pvarValue ) {
    HRESULT _hr = raw_GetTuneDataValue(nSegmentNumber, bstrLabel, pvarValue);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetTuneDataLabels ( long nSegmentNumber, VARIANT * pvarLabels, long * pnArraySize ) {
    HRESULT _hr = raw_GetTuneDataLabels(nSegmentNumber, pvarLabels, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstName ( BSTR * pbstrInstName ) {
    HRESULT _hr = raw_GetInstName(pbstrInstName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstModel ( BSTR * pbstrInstModel ) {
    HRESULT _hr = raw_GetInstModel(pbstrInstModel);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstSerialNumber ( BSTR * pbstrInstSerialNumber ) {
    HRESULT _hr = raw_GetInstSerialNumber(pbstrInstSerialNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstSoftwareVersion ( BSTR * pbstrInstSoftwareVersion ) {
    HRESULT _hr = raw_GetInstSoftwareVersion(pbstrInstSoftwareVersion);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstHardwareVersion ( BSTR * pbstrInstHardwareVersion ) {
    HRESULT _hr = raw_GetInstHardwareVersion(pbstrInstHardwareVersion);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstFlags ( BSTR * pbstrInstFlags ) {
    HRESULT _hr = raw_GetInstFlags(pbstrInstFlags);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstNumChannelLabels ( long * pnInstNumChannelLabels ) {
    HRESULT _hr = raw_GetInstNumChannelLabels(pnInstNumChannelLabels);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstChannelLabel ( long nChannelLabelNumber, BSTR * pbstrInstChannelLabel ) {
    HRESULT _hr = raw_GetInstChannelLabel(nChannelLabelNumber, pbstrInstChannelLabel);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetNumberOfControllersOfType ( long nControllerType, long * pnNumControllersOfType ) {
    HRESULT _hr = raw_GetNumberOfControllersOfType(nControllerType, pnNumControllersOfType);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetAverageMassList ( long * pnFirstAvgScanNumber, long * pnLastAvgScanNumber, long * pnFirstBkg1ScanNumber, long * pnLastBkg1ScanNumber, long * pnFirstBkg2ScanNumber, long * pnLastBkg2ScanNumber, _bstr_t bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double * pdCentroidPeakWidth, VARIANT * pvarMassList, VARIANT * pvarPeakFlags, long * pnArraySize ) {
    HRESULT _hr = raw_GetAverageMassList(pnFirstAvgScanNumber, pnLastAvgScanNumber, pnFirstBkg1ScanNumber, pnLastBkg1ScanNumber, pnFirstBkg2ScanNumber, pnLastBkg2ScanNumber, bstrFilter, nIntensityCutoffType, nIntensityCutoffValue, nMaxNumberOfPeaks, bCentroidResult, pdCentroidPeakWidth, pvarMassList, pvarPeakFlags, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::IsThereMSData ( long * pbMSData ) {
    HRESULT _hr = raw_IsThereMSData(pbMSData);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::HasExpMethod ( long * pbMethod ) {
    HRESULT _hr = raw_HasExpMethod(pbMethod);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetFilterMassPrecision ( long * pnFilterMassPrecision ) {
    HRESULT _hr = raw_GetFilterMassPrecision(pnFilterMassPrecision);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogForPos ( long nPos, VARIANT * pvarRT, VARIANT * pvarValue, long * pnArraySize ) {
    HRESULT _hr = raw_GetStatusLogForPos(nPos, pvarRT, pvarValue, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetStatusLogPlottableIndex ( VARIANT * pvarIndex, VARIANT * pvarValues, long * pnArraySize ) {
    HRESULT _hr = raw_GetStatusLogPlottableIndex(pvarIndex, pvarValues, pnArraySize);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetInstMethodNames ( long * pnNumInstMethods, VARIANT * pvarNames ) {
    HRESULT _hr = raw_GetInstMethodNames(pnNumInstMethods, pvarNames);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::SetMassTolerance ( long bUseUserDefined, double dMassTolerance, long nUnits ) {
    HRESULT _hr = raw_SetMassTolerance(bUseUserDefined, dMassTolerance, nUnits);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile::GetChros ( long nNumChros, double * pdStartTime, double * pdEndTime, VARIANT * pvarChroParamsArray, VARIANT * pvarChroDataSizeArray, VARIANT * pvarChroDataArray, VARIANT * pvarPeakFlagsArray ) {
    HRESULT _hr = raw_GetChros(nNumChros, pdStartTime, pdEndTime, pvarChroParamsArray, pvarChroDataSizeArray, pvarChroDataArray, pvarPeakFlagsArray);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

//
// interface IXRawfile2 wrapper method implementations
//

inline HRESULT IXRawfile2::GetLabelData ( VARIANT * pvarLabels, VARIANT * pvarFlags, long * pnScanNumber ) {
    HRESULT _hr = raw_GetLabelData(pvarLabels, pvarFlags, pnScanNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXRawfile2::GetNoiseData ( VARIANT * pvarNoisePacket, long * pnScanNumber ) {
    HRESULT _hr = raw_GetNoiseData(pvarNoisePacket, pnScanNumber);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

//
// interface IXVirMS wrapper method implementations
//

inline HRESULT IXVirMS::Create ( LPWSTR szFileName ) {
    HRESULT _hr = raw_Create(szFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::Close ( ) {
    HRESULT _hr = raw_Close();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::GetFileName ( BSTR * pbstrFileName ) {
    HRESULT _hr = raw_GetFileName(pbstrFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline long IXVirMS::GetIsError ( ) {
    long _result = 0;
    HRESULT _hr = get_IsError(&_result);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _result;
}

inline long IXVirMS::GetErrorCode ( ) {
    long _result = 0;
    HRESULT _hr = get_ErrorCode(&_result);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _result;
}

inline long IXVirMS::GetIsValid ( ) {
    long _result = 0;
    HRESULT _hr = get_IsValid(&_result);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _result;
}

inline HRESULT IXVirMS::WriteInstID ( LPWSTR szName, LPWSTR szModel, LPWSTR szSerialNumber, LPWSTR szSoftwareRev, LPWSTR ExpType ) {
    HRESULT _hr = raw_WriteInstID(szName, szModel, szSerialNumber, szSoftwareRev, ExpType);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteRunHeaderInfo ( double dExpectedRunTime, double dMassResolution, LPWSTR szComment1, LPWSTR szComment2 ) {
    HRESULT _hr = raw_WriteRunHeaderInfo(dExpectedRunTime, dMassResolution, szComment1, szComment2);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteInstData ( unsigned char * pcData, long nDataSize, enum MS_PacketTypes eType ) {
    HRESULT _hr = raw_WriteInstData(pcData, nDataSize, eType);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetTrailerHeaderNumFields ( long nFields ) {
    HRESULT _hr = raw_SetTrailerHeaderNumFields(nFields);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetTrailerHeaderField ( long nIdx, LPWSTR szLabel, enum MS_DataTypes eDataType, long nPrecision ) {
    HRESULT _hr = raw_SetTrailerHeaderField(nIdx, szLabel, eDataType, nPrecision);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteTrailerHeader ( ) {
    HRESULT _hr = raw_WriteTrailerHeader();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetStatusLogHeaderNumFields ( long nFields ) {
    HRESULT _hr = raw_SetStatusLogHeaderNumFields(nFields);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetStatusLogHeaderField ( long nIdx, LPWSTR szLabel, enum MS_DataTypes eDataType, long nPrecision ) {
    HRESULT _hr = raw_SetStatusLogHeaderField(nIdx, szLabel, eDataType, nPrecision);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteStatusLogHeader ( ) {
    HRESULT _hr = raw_WriteStatusLogHeader();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetTuneDataHeaderNumFields ( long nFields ) {
    HRESULT _hr = raw_SetTuneDataHeaderNumFields(nFields);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetTuneDataHeaderField ( long nIdx, LPWSTR szLabel, enum MS_DataTypes eDataType, long nPrecision ) {
    HRESULT _hr = raw_SetTuneDataHeaderField(nIdx, szLabel, eDataType, nPrecision);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteTuneDataHeader ( ) {
    HRESULT _hr = raw_WriteTuneDataHeader();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteTuneData ( unsigned char * pcData ) {
    HRESULT _hr = raw_WriteTuneData(pcData);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteStatusLog ( float fTime, unsigned char * pcData ) {
    HRESULT _hr = raw_WriteStatusLog(fTime, pcData);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteTrailer ( unsigned char * pcData ) {
    HRESULT _hr = raw_WriteTrailer(pcData);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::InitializeScanEvent ( struct MS_ScanEvent * pScanEvent ) {
    HRESULT _hr = raw_InitializeScanEvent(pScanEvent);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::InitMethodScanEvents ( ) {
    HRESULT _hr = raw_InitMethodScanEvents();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::SetMethodScanEvent ( long nSegment, long nScanEvent, struct MS_ScanEvent * pScanEvent ) {
    HRESULT _hr = raw_SetMethodScanEvent(nSegment, nScanEvent, pScanEvent);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteMethodScanEvents ( ) {
    HRESULT _hr = raw_WriteMethodScanEvents();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteScanIndex ( struct MS_ScanIndex * pScanIndex, struct MS_ScanEvent * pScanEvent ) {
    HRESULT _hr = raw_WriteScanIndex(pScanIndex, pScanEvent);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteInstData2 ( long nNumPkts, struct MS_DataPeak * pPackets ) {
    HRESULT _hr = raw_WriteInstData2(nNumPkts, pPackets);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::InitializeScanIndex ( long nScanIndexPosition, enum MS_PacketTypes eType ) {
    HRESULT _hr = raw_InitializeScanIndex(nScanIndexPosition, eType);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirMS::WriteScanIndex2 ( struct MS_ScanIndex * pScanIndex ) {
    HRESULT _hr = raw_WriteScanIndex2(pScanIndex);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

//
// interface IXVirUV wrapper method implementations
//

inline HRESULT IXVirUV::Create ( LPWSTR szFileName ) {
    HRESULT _hr = raw_Create(szFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirUV::Close ( ) {
    HRESULT _hr = raw_Close();
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirUV::GetFileName ( BSTR * pbstrFileName ) {
    HRESULT _hr = raw_GetFileName(pbstrFileName);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline long IXVirUV::GetIsError ( ) {
    long _result = 0;
    HRESULT _hr = get_IsError(&_result);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _result;
}

inline long IXVirUV::GetErrorCode ( ) {
    long _result = 0;
    HRESULT _hr = get_ErrorCode(&_result);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _result;
}

inline long IXVirUV::GetIsValid ( ) {
    long _result = 0;
    HRESULT _hr = get_IsValid(&_result);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _result;
}

inline HRESULT IXVirUV::WriteErrorLog ( float fTime, LPWSTR szError ) {
    HRESULT _hr = raw_WriteErrorLog(fTime, szError);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirUV::WriteInstID ( LPWSTR szName, LPWSTR szModel, LPWSTR szSerialNumber, LPWSTR szSoftwareRev, LPWSTR szLabel1, LPWSTR szLabel2, LPWSTR szLabel3, LPWSTR szLabel4 ) {
    HRESULT _hr = raw_WriteInstID(szName, szModel, szSerialNumber, szSoftwareRev, szLabel1, szLabel2, szLabel3, szLabel4);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirUV::WriteRunHeaderInfo ( double dExpectedRunTime ) {
    HRESULT _hr = raw_WriteRunHeaderInfo(dExpectedRunTime);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirUV::WriteInstData ( unsigned char * pcData, long nDataSize, enum MS_PacketTypes eType, long nDataLen ) {
    HRESULT _hr = raw_WriteInstData(pcData, nDataSize, eType, nDataLen);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}

inline HRESULT IXVirUV::WriteScanIndex ( struct MS_UVScanIndex * pScanIndex ) {
    HRESULT _hr = raw_WriteScanIndex(pScanIndex);
    if (FAILED(_hr)) _com_issue_errorex(_hr, this, __uuidof(this));
    return _hr;
}