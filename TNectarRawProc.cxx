#define CAL_XRAY90_A 0.0425282
#define CAL_XRAY90_B 0.0

#define CAL_XRAY145_A 0.037211
#define CAL_XRAY145_B 0.0


#define SCALER_RANGE 14400  //Previous value defined by Jan Glorius 7200 //in seconds

#include "TNectarRawProc.h"

#include "Riostream.h"

#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include "TGo4Analysis.h"
#include "TGo4Log.h"

#include "TGo4Picture.h"
#include "TGo4MbsEvent.h"
#include "TGo4UserException.h"
#include "TGo4WinCond.h"

#include "TNectarRawEvent.h"
#include "TNectarRawParam.h"
#include "TNectarDisplay.h"

#include "TGo4Fitter.h"
#include "TGo4FitDataHistogram.h"
#include "TGo4FitParameter.h"
#include "TGo4FitModelPolynom.h"

#include "TGo4UserException.h"

/** enable this definition to print out event sample data explicitely*/
//#define NECTAR_VERBOSE_PRINT 1
#ifdef NECTAR_VERBOSE_PRINT
#define printdeb( args... )  printf( args );
#else
#define printdeb( args...) ;
#endif


#define NECTAR_SKIP_EVENT( args...)\
  printf(args);\
  throw TGo4UserException(1,nullptr);

//GO4_SKIP_EVENT


static ULong_t skipped_events = 0;

/* helper macros for BuildEvent to check if payload pointer is still inside delivered region:*/
/* this one to be called at top data processing loop*/

#define  NECTAR_EVENT_CHECK_PDATA                                    \
if((pData - pSubevt->GetDataField()) > lWords ) \
{ \
  printf("############ unexpected end of event for subevent size :0x%x, skip event %ld\n", lWords, skipped_events++);\
  throw TGo4UserException(1,nullptr);\
  continue; \
}

//  GO4_SKIP_EVENT



//***********************************************************
TNectarRawProc::TNectarRawProc() :
    TGo4EventProcessor(), fPar(0), fNectarRawEvent(0), pSubevt(0), pData(0), lWords(0)
{
}

//***********************************************************
// this one is used in standard factory
TNectarRawProc::TNectarRawProc(const char* name) :
    TGo4EventProcessor(name), fPar(0), fNectarRawEvent(0), pSubevt(0), pData(0), lWords(0)
{
  TGo4Log::Info("TNectarRawProc: Create instance %s", name);
  fMdppDisplays.clear();
  fVmmrDisplays.clear();
  
  SetMakeWithAutosave(kTRUE);
//// init user analysis objects:

  fPar = dynamic_cast<TNectarRawParam*>(MakeParameter("NectarRawParam", "TNectarRawParam", "set_NectarRawParam.C"));
  if (fPar)
    fPar->SetConfigBoards();    // board configuration specified in parameter is used to build subevents and displays

  for (UInt_t i = 0; i < TNectarRawEvent::fgConfigMdppBoards.size(); ++i)
  {
    UInt_t uniqueid = TNectarRawEvent::fgConfigMdppBoards[i];
    Bool_t sixteenchannels = fPar->fMDPP_is16Channels[i];    // prepare displays for configured 16 or 32 channel MDPP
    fMdppDisplays.push_back(new TMdppDisplay(uniqueid, sixteenchannels));
  }

  for (UInt_t i = 0; i < TNectarRawEvent::fgConfigVmmrBoards.size(); ++i)
  {
    UInt_t uniqueid = TNectarRawEvent::fgConfigVmmrBoards[i];
    fVmmrDisplays.push_back(new TVmmrDisplay(uniqueid));
  }

  //add ADC/TDC/Scaler histos for Caen Modules by hand
  Int_t adc_range = 8192;
  Int_t tdc_range = 8192;
  uint32_t scaler_range = SCALER_RANGE;
  
  TString obname;
  TString obtitle;

  //////////////////
  // generic histos
  //////////////////

  for(int ch=0;ch<MADC_NUMCHANNELS;ch++)
      {
        obname.Form("Raw/MADC/ADC_Channel_%d", ch);
        obtitle.Form("mesytec MADC - Channel %d ", ch);
        hRaw_MADC[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), MADC_RANGE, 0.5, MADC_RANGE+0.5, "ADC value", "counts");
      }



  //ADC histo definition
  //for(int ch=0;ch<3;ch++)
  //  {
  //    obname.Form("Raw/ADC/ADC_Channel_%d", ch);
  //    obtitle.Form("Caen V785 ADC - Channel %d ", ch);
  //    hRawV785_ADC[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), adc_range, 0.5, adc_range+0.5, "ADC value", "counts");
  //  }
  ////TDC histo definition
  for(int ch=0;ch<32;ch++)
    {
      obname.Form("Raw/TDC/TDC_Channel_%d", ch);
      obtitle.Form("Caen V775 TDC - Channel %d ", ch);
      hRawV775_TDC[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), tdc_range, 0.5, tdc_range+0.5, "TDC value", "counts");
    }
  ////V830 scaler histo definition
  //for(int ch=0;ch<32;ch++)
  //  {
  //    obname.Form("Raw/Scaler/V830_Channel_%d", ch);
  //    obtitle.Form("Caen V830 Scaler - Channel %d ", ch);
  //    hRawV830_Scaler[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), tdc_range, 0.5, scaler_range+0.5, "scaler value", "counts");
  //  }
  // generic TRLO histos
  //for(int ch=0;ch<16;ch++)
  //  {
  //    obname.Form("TRLO/%d_bLMU", ch);
  //    obtitle.Form("TRLO before LMU - Channel %d ", ch);
  //    hRawTRLO_blmu[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), scaler_range, 0, scaler_range, "time", "rate");
  //
  //    obname.Form("TRLO/%d_bDT", ch);
  //    obtitle.Form("TRLO before Deadtime - Channel %d ", ch);
  //    hRawTRLO_bdt[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), scaler_range, 0, scaler_range, "scaler value", "rate");
  //
  //    obname.Form("TRLO/%d_aDT", ch);
  //    obtitle.Form("TRLO after Deadtime - Channel %d ", ch);
  //    hRawTRLO_adt[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), scaler_range, 0, scaler_range, "scaler value", "counts");
  //
  //    obname.Form("TRLO/%d_aRED", ch);
  //    obtitle.Form("TRLO after Reduction - Channel %d ", ch);
  //    hRawTRLO_ared[ch] = MakeTH1('I',obname.Data(), obtitle.Data(), scaler_range,0, scaler_range, "scaler value", "counts");
  //    
  //  }

  //////////////////
  // mapped histos
  /////////////////
  
  //TRLOii tpat & scaler histo definition
  hRawTRLO_tpat = MakeTH1('I',"TRLO/tpat", "trigger pattern number (2^tpat)", 256, 0.5, 256.5, "tpat", "counts");
  hRawTRLO_trigger = MakeTH1('I',"TRLO/trigger", "MBS trigger number", 16, 0, 16, "trigger", "counts");
  //Xray histos
  hRawV785_ADC[1] = MakeTH1('I',"XRAY/ADC_Xray_90", "ADC values Xray 90 degree",   adc_range, 0.5,   adc_range+0.5, "ADC values", "counts");
  hRawV785_ADC[2] = MakeTH1('I',"XRAY/ADC_Xray_145", "ADC values Xray 145 degree",   adc_range, 0.5,   adc_range+0.5, "ADC values", "counts");
  h_E_Xray[1] = MakeTH1('I',"XRAY/E_Xray_90", "Energy Xray 90 degree",   adc_range, 0.5*CAL_XRAY90_A + CAL_XRAY90_B,   (adc_range+0.5)*CAL_XRAY90_A + CAL_XRAY90_B, "Energy [keV]", "counts");
  h_E_Xray[2] = MakeTH1('I',"XRAY/E_Xray_145", "Energy Xray 145 degree", adc_range, 0.5*CAL_XRAY145_A + CAL_XRAY145_B, (adc_range+0.5)*CAL_XRAY145_A + CAL_XRAY145_B, "Energy [keV]", "counts");
  h_t_Xray[1] = MakeTH1('I',"XRAY/t_Xray_90", "time Xray 90 degree",   tdc_range, 0.5, tdc_range+0.5, "TDC value", "counts");
  h_t_Xray[2] = MakeTH1('I',"XRAY/t_Xray_145", "time Xray 145 degree",   tdc_range, 0.5, tdc_range+0.5, "TDC value", "counts");
  h_E_Xray_OFF[1] = MakeTH1('I',"XRAY/E_Xray_90_OFF", "Energy Xray 90 degree - jet OFF",   adc_range, 0.5*CAL_XRAY90_A + CAL_XRAY90_B,   (adc_range+0.5)*CAL_XRAY90_A + CAL_XRAY90_B, "Energy [keV]", "counts");
  h_E_Xray_OFF[2] = MakeTH1('I',"XRAY/E_Xray_145_OFF", "Energy Xray 145 degree - jet OFF", adc_range, 0.5*CAL_XRAY145_A + CAL_XRAY145_B, (adc_range+0.5)*CAL_XRAY145_A + CAL_XRAY145_B, "Energy [keV]", "counts");
  h_t_Xray_OFF[1] = MakeTH1('I',"XRAY/t_Xray_90_OFF", "time Xray 90 degree - jet OFF",   tdc_range, 0.5, tdc_range+0.5, "TDC value", "counts");
  h_t_Xray_OFF[2] = MakeTH1('I',"XRAY/t_Xray_145_OFF", "time Xray 145 degree - jet OFF",   tdc_range, 0.5, tdc_range+0.5, "TDC value", "counts");

  //beam scalers
  h_trafo = MakeTH1('I',"ESR/trafo", "ESR beam trafo", 7200, 0, 7200,"test", "test" );
  h_I_cooler = MakeTH1('I',"ESR/Icool", "ESR cooler current", scaler_range, 0, scaler_range, "time", "rate");
  h_U_cooler = MakeTH1('I',"ESR/Ucool", "ESR cooler voltage", scaler_range, 0, scaler_range, "time", "rate");
  h_jet_S1   = MakeTH1('I',"ESR/jetS1", "ESR target pressure S1", scaler_range, 0, scaler_range, "time", "rate");
  h_jet_S2   = MakeTH1('I',"ESR/jetS2", "ESR target pressure S2", scaler_range, 0, scaler_range, "time", "rate");
  h_pmt      = MakeTH1('I',"ESR/pmt", "ESR overlap PMT", scaler_range, 0, scaler_range, "time", "rate");

  //h_TRAFO->Fill(111.1,200);
  
  //TRLO scalers
  h_Tel_blmu = MakeTH1('I',"TRLO/Tel_bLMU", "Telescope before LMU", scaler_range, 0, scaler_range, "time", "rate");
  h_Tel_bdt =  MakeTH1('I',"TRLO/Tel_bDT",  "Telescope before DT" , scaler_range, 0, scaler_range, "time", "rate");
  h_Tel_adt =  MakeTH1('I',"TRLO/Tel_aDT",  "Telescope after DT"  , scaler_range, 0, scaler_range, "time", "rate");
  h_Tel_ared = MakeTH1('I',"TRLO/Tel_aRED", "Telescope after RED" , scaler_range, 0, scaler_range, "time", "rate");

  h_HRes_blmu = MakeTH1('I',"TRLO/HRes_bLMU", "HeavyRes. before LMU", scaler_range, 0, scaler_range, "time", "rate");
  h_HRes_bdt =  MakeTH1('I',"TRLO/HRes_bDT",  "HeavyRes. before DT" , scaler_range, 0, scaler_range, "time", "rate");
  h_HRes_adt =  MakeTH1('I',"TRLO/HRes_aDT",  "HeavyRes. after DT"  , scaler_range, 0, scaler_range, "time", "rate");
  h_HRes_ared = MakeTH1('I',"TRLO/HRes_aRED", "HeavyRes. after RED" , scaler_range, 0, scaler_range, "time", "rate");

  h_Xray1_blmu = MakeTH1('I',"TRLO/Xray1_bLMU", "Xray1 before LMU", scaler_range, 0, scaler_range, "time", "rate");
  h_Xray2_blmu = MakeTH1('I',"TRLO/Xray2_bLMU", "Xray2 before LMU", scaler_range, 0, scaler_range, "time", "rate");
  h_Xray3_blmu = MakeTH1('I',"TRLO/Xray3_bLMU", "Xray3 before LMU", scaler_range, 0, scaler_range, "time", "rate");
  h_Xray_bdt =  MakeTH1('I',"TRLO/Xray_bDT",  "Xray before DT" , scaler_range, 0, scaler_range, "time", "rate");
  h_Xray_adt =  MakeTH1('I',"TRLO/Xray_aDT",  "Xray after DT"  , scaler_range, 0, scaler_range, "time", "rate");
  h_Xray_ared = MakeTH1('I',"TRLO/Xray_aRED", "Xray after RED" , scaler_range, 0, scaler_range, "time", "rate");
  

  h_vulom_raw=MakeTH1('I',"Raw/VulomScalerRaw", "Raw content of VULOM" , VULOM_NUMSCALERS, 0, VULOM_NUMSCALERS, "scaler channel", "counts");

  InitDisplay(false);    // will init all subdisplays, creating histograms etc.

}

//***********************************************************
TNectarRawProc::~TNectarRawProc()
{
  std::cout << "**** TNectarRawProc dtor " << std::endl;

  for (UInt_t i = 0; i < fMdppDisplays.size(); ++i)
  {
    delete fMdppDisplays[i];
  }
  for (UInt_t i = 0; i < fVmmrDisplays.size(); ++i)
  {
    delete fVmmrDisplays[i];
  }
}

TMdppDisplay* TNectarRawProc::GetMdppDisplay(UInt_t uniqueid)
{
  for (UInt_t i = 0; i < fMdppDisplays.size(); ++i)
  {
    TMdppDisplay* theDisplay = fMdppDisplays[i];
    if (uniqueid == theDisplay->GetDevId())
      return theDisplay;
  }
  return 0;
}

TVmmrDisplay* TNectarRawProc::GetVmmrDisplay(UInt_t uniqueid)
{
  for (UInt_t i = 0; i < fVmmrDisplays.size(); ++i)
  {
    TVmmrDisplay* theDisplay = fVmmrDisplays[i];
    if (uniqueid == theDisplay->GetDevId())
      return theDisplay;
  }
  return 0;
}

void TNectarRawProc::InitDisplay(Bool_t replace)
{
  std::cout << "**** TNectarRawProc: Init Display with replace= " << replace << std::endl;

  for (UInt_t i = 0; i < fMdppDisplays.size(); ++i)
  {
    fMdppDisplays[i]->InitDisplay(0, replace);
  }
  for (UInt_t i = 0; i < fVmmrDisplays.size(); ++i)
  {
    fVmmrDisplays[i]->InitDisplay(0, replace);
  }

}

//-----------------------------------------------------------
// event function
Bool_t TNectarRawProc::BuildEvent(TGo4EventElement *target)
{
// called by framework from TNectarRawEvent to fill it
  fNectarRawEvent = (TNectarRawEvent*) target;
  fNectarRawEvent->SetValid(kFALSE);    // not store
  TGo4MbsEvent *source = (TGo4MbsEvent*) GetInputEvent();
  if (source == 0)
  {
    std::cout << "AnlProc: no input event !" << std::endl;
    return kFALSE;
  }
  UShort_t triggertype = source->GetTrigger();
  fNectarRawEvent->ftrigger = source->GetTrigger();
  
  if (triggertype > 11)
  {
    // frontend offset trigger can be one of these, we let it through to unpacking loop
    //cout << "**** TNectarRawProc: Skip trigger event" << endl;
    //GO4_SKIP_EVENT; - compiler warnings recently...
    //NECTAR_SKIP_EVENT("**** TNectarRawProc: Skip event of trigger type %d\n", triggertype); // works, but floods the terminal
    //throw TGo4UserException(1,""); // this will do it silently JAM24
    return kFALSE; // this would let the second step execute!
  }

// first we fill the TNectarRawEvent with data from MBS source
// we have up to two subevents, crate 1 and 2
// Note that one has to loop over all subevents and select them by
// crate number:   pSubevt->GetSubcrate(),
// procid:         pSubevt->GetProcid(),
// and/or control: pSubevt->GetControl()
// here we use only crate number

//  take general event number from mbs event header. Note that subsystem sequence may differ:
  fNectarRawEvent->fSequenceNumber = source->GetCount();

  source->ResetIterator();
  while ((pSubevt = source->NextSubEvent()) != 0)
  {    // loop over subevents

    // JAM here one can exclude data from other subsystem by these id numbers
    // see mbs setup.usf how this is defined!
    if (pSubevt->GetSubcrate() != 0)
      continue;
    if (pSubevt->GetControl() != 9)
      continue;
//  if (pSubevt->GetProcid() != 1)
//    continue;

    pData = pSubevt->GetDataField();
    lWords = pSubevt->GetIntLen();

    if ((UInt_t) *pData == 0xbad00bad)
    {
      NECTAR_SKIP_EVENT("**** TNectarRawProc: Found BAD mbs event (marked 0x%x), skip it.\n", (*pData));
    }

    // loop over single subevent data:
    int testiter = 0;
    while (pData - pSubevt->GetDataField() < lWords)
    {

      // first header word indicates board type:
      // 0xA - VMMR, 0xB - MDPP
      // also module id and firmware version is given
      Int_t header = *pData++;
      //Bool_t isVMMR = ((header >> 28) & 0xF) == 0xA;
      Bool_t isVMMR = ((header >> 24) & 0xFF) == 0xA0;    // mesytec VMMR
      Bool_t isMADC = ((header >> 24) & 0xFF) == 0xA1;    // mesytec MADC
      Bool_t isMDPP = ((header >> 24) & 0xFF) == 0xB0;    // mesytec MDPP
      Bool_t isMTDC = ((header >> 24) & 0xFF) == 0xB1;    // mesytec MTDC
      Bool_t isV775 = ((header >> 24) & 0xFF) == 0xC0;    //caen TDC v775
      Bool_t isV785 = ((header >> 24) & 0xFF) == 0xD0;    //caen ADC v785
      Bool_t isV830 = ((header >> 24) & 0xFF) == 0xE0;    //caen scaler v830
      Bool_t isTPAT = ((header >> 24) & 0xFF) == 0xcf;    //Vulom4 trloII tpat
      Bool_t isVSCA = ((header >> 24) & 0xFF) == 0xc7;    //Vulom4 trloII scalers
      Bool_t isPLAINVSCA = ((header >> 24) & 0xFF) == 0xc1;    //Vulom4 trloII scalers
      //  JAM todo: use following heaeder info somewhere?
      //Int_t module_nr = (header >> 16) & 0xF;
      //Int_t firmware = header & 0xFFFF;

      //printf("header: %x filter2: %x\n", header, ((header >> 24) & 0xFF));
      
      if (isVMMR)
      {
        if (!UnpackVmmr())
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of VMMR data failed! -  skip event  NEVER COME HERE!\n");
          // at the moment, any error in UnpackMdpp will throw exception, so you should never see this line
        }
      }
      else if (isMDPP)
      {
        if (!UnpackMdpp())
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of MDPP data failed! -  skip event  NEVER COME HERE!\n");
          // at the moment, any error in UnpackMdpp will throw exception, so you should never see this line
        }
      }
      else if (isV775)
      {
        //std::cout << "V775 header found!" << std::endl;
        if (!UnpackV775(triggertype))
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of V775 data failed! -  skip event\n");
        }
      }
      else if (isV785)
      {
        //std::cout << "V785 header found!" << std::endl;
        if (!UnpackV785(triggertype))
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of V785 data failed! -  skip event\n");
        }
      }
      else if (isV830)
      {
        //std::cout << "V830 header found!" << std::endl;
        if (!UnpackV830())
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of V830 data failed! -  skip event\n");
        }
      }
      else if (isTPAT)
      {
        //std::cout << "TPAT header found!" << std::endl;
        if (!UnpackTPAT(header))
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of TRLLOII TPAT data failed! -  skip event\n");
        }
      }
      else if (isVSCA)
      {
        //std::cout << "VScaler header found!" << std::endl;
        if (!UnpackVSCA(header))
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of TRLLOII Scaler data failed! -  skip event\n");
        }
      }
      else if (isPLAINVSCA)
      {
        //std::cout << "plain VULOM Scaler header found!" << std::endl;
        if (!UnpackPlainVSCA())
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of plain VULUM Scaler data failed! -  skip event\n");
        }
      }

      else if (isMADC)
      {
         //std::cout << "MADC header found!" << std::endl;
        if (!UnpackMadc(triggertype))
        {
          NECTAR_SKIP_EVENT("Data error: unpacking of MADC data failed! -  skip event\n");
        }
      }

      else if (isMTDC)
      {
        NECTAR_SKIP_EVENT("Data error: found unexpected MTDC data without implemented unpacker! - skip event!\n");
      }

      else
      {
        NECTAR_SKIP_EVENT("Data error: incompatible header - skip event! iter=%d word=0x%x \n", testiter, header);    // Removed and added printf command - 19 May 22
        return kFALSE;
      }

    }    // while pData - pSubevt->GetDataField() <lWords

  }    // while subevents

//

  UpdateDisplays();    // we fill the raw displays immediately, but may do additional histogramming later
  fNectarRawEvent->SetValid(kTRUE);    // to store

  if (fPar->fSlowMotion)
  {
    Int_t evnum = source->GetCount();
    GO4_STOP_ANALYSIS_MESSAGE(
        "Stopped for slow motion mode after MBS event count %d. Click green arrow for next event. please.", evnum);
    //printf("Stopped for slow motion mode after MBS event count %d. Click green arrow for next event please. \n");
  }

  return kTRUE;
}

Bool_t TNectarRawProc::UnpackMdpp()
{
// mdpp header should follow
  Int_t MDPP_head = *pData++;
  Bool_t isHeader = ((MDPP_head >> 28) & 0xF) == 0x4;
  if (!isHeader){
   NECTAR_SKIP_EVENT("Data error: invalid mdpp header 0x%x -  skip event!\n", MDPP_head);
  }
  Int_t sublen = MDPP_head & 0x3FF;    // number of following data words,including EOE
  UInt_t modid = (MDPP_head >> 16) & 0xFF;    //module id
  Int_t resolution = (MDPP_head >> 10) & 0x3F;    //3 bit TDC_resolution 0x604 , 3bit ADC_resolution 0x6046

// JAm 13-09-21: evaluate if MDPP has 16 or 32 channel format
  Bool_t mddpHas16Channels = kFALSE;
  Int_t ix = fPar->GetMDPPArrayindex(modid);
  if (ix >= 0)
    mddpHas16Channels = fPar->fMDPP_is16Channels[ix];

  TMdppBoard* theBoard = fNectarRawEvent->GetMdppBoard(modid);
  if (theBoard == 0)
  {
    GO4_STOP_ANALYSIS_MESSAGE(
        "Configuration error: MDPP module id %d does not exist as subevent. Please check TNectarRawParam setup", modid);
    return kFALSE;
  }
  theBoard->SetResolution(resolution);

  TMdppDisplay* boardDisplay = GetMdppDisplay(modid);
  if (boardDisplay == 0)
  {
    GO4_STOP_ANALYSIS_MESSAGE("Configuration error: MDPP module id %d does not exist as histogram display set!", modid);
    return kFALSE;
  }
  boardDisplay->ResetDisplay();

  boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_HEADER);    // always account header if found

// loop over following data words:
  Int_t* pdatastart = pData;    // remember begin of MDPP payload data section
  while ((pData - pdatastart) < sublen)
  {
    Int_t word = *pData++;
    Bool_t isTimestamp = ((word >> 28) & 0xF) == 0x2;
    Bool_t isData = ((word >> 28) & 0xF) == 0x1;
    Bool_t isEndmark = ((word >> 30) & 0x3) == 0x3;
    Bool_t isDummy = ((word >> 30) & 0x3) == 0x0;
    if (isTimestamp)
    {
      UShort_t ts = word & 0xFFFF;
      theBoard->fExtendedTimeStamp = ts;
      boardDisplay->hExtTimeStamp->Fill(ts);
      boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_TIMESTAMP);
    }
    else if (isData)
    {
      if (mddpHas16Channels)
      {
        // we expect MDPP 16 channel variant from configuration
        UShort_t maxchannel=16;
        UShort_t channel = (word >> 16) & 0x1F;
        UShort_t value = word & 0xFFFF;
        printdeb("UnpackMdpp 16 has data word 0x%x, channel=%d, value=0x%x\n", word, channel, value)
        Bool_t isTDC = channel > 15 ? kTRUE : kFALSE;
        if (isTDC)
        {
          channel -= 16;
          if (channel >= maxchannel)
          {
            printf("Warning: MDPP 16 TDC channel info %d out of bounds, maxchannel=%d\n", channel, maxchannel);
            continue;    // maybe better stop here ?
          }
          // put data to output event structure (tree file, next analysis steps...)
          theBoard->AddTdcMessage(new TMdppTdcData(value), channel);
          // now do histogramming:
          boardDisplay->hRawTDC[channel]->Fill(value);
          boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_TDC);
          boardDisplay->hTDC_ChannelScaler->Fill(channel);
        }
        else
        {
          // ADC:
          if (channel >= maxchannel)
          {
            printf("Warning: MDPP ADC channel info %d out of bounds, maxchannel=%d\n", channel, maxchannel);
            continue;    // maybe better stop here ?
          }

          theBoard->AddAdcMessage(new TMdppAdcData(value), channel);
          boardDisplay->hRawADC[channel]->Fill(value);
          boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_ADC);
          boardDisplay->hADC_ChannelScaler->Fill(channel);
        }

      }    //if(mddpHas16Channels)
      else
      {
        // we expect MDPP 32 channel variant from configuration
        UShort_t channel = (word >> 16) & 0x7F;
        UShort_t value = word & 0xFFFF;
        printdeb("UnpackMdpp 32 has data word 0x%x, channel=%d, value=0x%x\n", word, channel, value)
        Bool_t isTDC = channel > 31 ? kTRUE : kFALSE;
        if (isTDC)
        {
          // check here for extra trigger time difference information at highest channels:
          if (channel > 63)
          {
            UShort_t trigchan = channel - 64;
            if (trigchan >= MDPP_EXTDTCHANNELS)
            {
              printf("Warning: MDPP 32 external trigger input %d out of bounds, maxinputs=%d\n", trigchan, MDPP_EXTDTCHANNELS);
              continue;    // maybe better stop here ?
            }
            theBoard->AddExtDTMessage(new TMdppTdcData(value), trigchan);
            boardDisplay->hExtTrigTime[trigchan]->Fill(value);
          }

          else
          {
            // regular TDC
            channel -= 32;
            if (channel >= MDPP_CHANNELS)
            {
              printf("Warning: MDPP 32 TDC channel info %d out of bounds, maxchannel=%d\n", channel, MDPP_CHANNELS);
              continue;    // maybe better stop here ?
            }
            // put data to output event structure (tree file, next analysis steps...)
            theBoard->AddTdcMessage(new TMdppTdcData(value), channel);
            // now do histogramming:
            boardDisplay->hRawTDC[channel]->Fill(value);
            boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_TDC);
            boardDisplay->hTDC_ChannelScaler->Fill(channel);
          }
        }
        else
        {
          // ADC
          if (channel >= MDPP_CHANNELS)
          {
            printf("Warning: MDPP ADC channel info %d out of bounds, maxchannel=%d\n", channel, MDPP_CHANNELS);
            continue;    // maybe better stop here ?
          }
          theBoard->AddAdcMessage(new TMdppAdcData(value), channel);
          boardDisplay->hRawADC[channel]->Fill(value);
          boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_ADC);
          boardDisplay->hADC_ChannelScaler->Fill(channel);
        }

      }    //if(mddpHas16Channels)
    }
    else if (isEndmark)
    {
      Int_t count = (word & 0x3FFFFFFF);
      theBoard->fEventCounter = count;
      boardDisplay->hMsgTypes->Fill(TMdppMsg::MSG_EOE);

      // now check against global counter of all subsystems:
      if (fNectarRawEvent->fSequenceNumber < 0)
      {
        fNectarRawEvent->fSequenceNumber = count;    // init if not yet set by other board
        printdeb("Warning: MDPP board %d setting sequence number to %d \n", modid, count);

        // TODO: maybe we could put MBS header counter to fSequenceNumber, but this may differ by special trigger events
      }
      else
      {
        if (fNectarRawEvent->fSequenceNumber != count)    // check against global counter
        {
          printdeb("Warning: MDPP board %d eventcounter %d is not matching common counter %d \n", modid, count,
              fNectarRawEvent->fSequenceNumber);
          // JAM - optionally we may stop here in such case?
          //GO4_STOP_ANALYSIS_MESSAGE("Warning: MDPP board %d eventcounter %d is not matching common counter %d \n", modid, count, fNectarRawEvent->fSequenceNumber);
        }
      }

    }
    else if (isDummy)
    {
      // do nothing, just proceed.
    }

    else
    {
      GO4_STOP_ANALYSIS_MESSAGE("Data error: unknown format in MDPP word 0x%x from module id %d-  stopped!", word,
          modid);
      return kFALSE;    // redundant, above will do exception
    }
    // following check will avoid that sublen of frontend is wrong..
    NECTAR_EVENT_CHECK_PDATA
  }    //  while ((pData - pdatastart) < sublen)

  return kTRUE;
}

Bool_t TNectarRawProc::UnpackVmmr()
{
  Int_t VMMR_head = *pData++;
  Bool_t isHeader = ((VMMR_head >> 28) & 0xF) == 0x4;
  if (!isHeader){
    NECTAR_SKIP_EVENT("Data error: invalid vmmr header 0x%x -  skip event!\n", VMMR_head);
  }
  Int_t sublen = VMMR_head & 0xFFF;    // number of following data words,including EOE
  Int_t modid = (VMMR_head >> 16) & 0xFF;    //module id

  TVmmrBoard* theBoard = fNectarRawEvent->GetVmmrBoard(modid);
  if (theBoard == 0)
  {
    GO4_STOP_ANALYSIS_MESSAGE(
        "Configuration error: VMMR module id %d does not exist as subevent. Please check TNectarRawParam setup", modid);
    return kFALSE;
  }

  TVmmrDisplay* boardDisplay = GetVmmrDisplay(modid);
  if (boardDisplay == 0)
  {
    GO4_STOP_ANALYSIS_MESSAGE("Configuration error: VMMR module id %d does not exist as histogram display set!", modid);
    return kFALSE;
  }
  boardDisplay->ResetDisplay();

  boardDisplay->hMsgTypes->Fill(TVmmrMsg::MSG_HEADER);    // always account header if found
// loop over following data words:
  Int_t* pdatastart = pData;    // remember begin of MDPP payload data section
  while ((pData - pdatastart) < sublen)
  {
    Int_t word = *pData++;
    Bool_t isTimestamp = ((word >> 28) & 0xF) == 0x2;
    Bool_t isAdc = ((word >> 28) & 0xF) == 0x1;
    //Bool_t isDeltaT = ((word >> 28) & 0xF) == 0x3;  // old definition  JAM
    Bool_t isEndmark = ((word >> 30) & 0x3) == 0x3;
    Bool_t isDummy = ((word >> 30) & 0x3) == 0x0;
    Bool_t isDeltaT=kFALSE;
    Bool_t isExtendedTS=kFALSE;


    if(isTimestamp)
    {
      isExtendedTS= ((word >> 23) & 0x1) == 0x1;
      isDeltaT = !isExtendedTS;
    }

    if (isExtendedTS)
    {
      UShort_t ts = word & 0xFFFF;
      theBoard->fExtendedTimeStamp = ts;
      boardDisplay->hExtTimeStamp->Fill(ts);
      boardDisplay->hMsgTypes->Fill(TVmmrMsg::MSG_TIMESTAMP);
    }
    else if (isAdc)
    {

      UChar_t slave = (word >> 24) & 0xF;
      UShort_t channel = (word >> 12) & 0xFFF;
      UShort_t value = word & 0xFFF;

      // first put extracted data to output event structure:
      TVmmrSlave* theSubBoard = theBoard->GetSlave(slave, kTRUE); // will create subcomponent if not exisiting
      if (theSubBoard == 0)
      {
        GO4_STOP_ANALYSIS_MESSAGE("NEVER  COME HERE: Could not access Vmmr subevent for FE slave %d", slave);
        //  JAM -maybe paranoid, but you never know...
      }
      theSubBoard->AddAdcMessage(new TVmmrAdcData(value, channel));

      // then fill histograms of corresponding container:
      boardDisplay->hMsgTypes->Fill(TVmmrMsg::MSG_ADC);

      TVmmrSlaveDisplay* slaveDisplay = boardDisplay->GetSlaveDisplay(slave);
      if (slaveDisplay == 0)
      {
        GO4_STOP_ANALYSIS_MESSAGE("NEVER  COME HERE: Could not access Vmmr display container for  FE slave %d", slave);

      }
      if (channel < VMMR_CHANNELS)
      {
        slaveDisplay->hRawADC[channel]->Fill(value);
      }
      else
      {
        // possible range of FE subaddress may be bigger than number of channels to plot:
        printf("VMMR %d slave %d warning: supress FE subadress %d for histogramming! please check setup.\n", modid,
            slave, channel);
      }

    }
    else if (isDeltaT)
    {
      UChar_t slave = (word >> 24) & 0xF;
      UShort_t value = word & 0xFFFF;

      // first put extracted data to output event structure:
      TVmmrSlave* theSubBoard = theBoard->GetSlave(slave, kTRUE);
      if (theSubBoard == 0)
      {
        GO4_STOP_ANALYSIS_MESSAGE("NEVER  COME HERE: Could not access Vmmr subevent for FE slave %d", slave);
        //  JAM -maybe paranoid, but you never know...
      }
      theSubBoard->fDeltaTGate = value;
      //GO4_SKIP_EVENT_MESSAGE("DT: ", value);

      // then fill histograms of corresponding container:
      boardDisplay->hMsgTypes->Fill(TVmmrMsg::MSG_DELTA_T);

      TVmmrSlaveDisplay* slaveDisplay = boardDisplay->GetSlaveDisplay(slave);
      if (slaveDisplay == 0)
      {
        GO4_STOP_ANALYSIS_MESSAGE("NEVER  COME HERE: Could not access Vmmr display container for  FE slave %d", slave);

      }
      slaveDisplay->hDeltaTime->Fill(value);

    }
    else if (isEndmark)
    {
      Int_t count = (word & 0x3FFFFFFF);
      theBoard->fEventCounter = count;
      boardDisplay->hMsgTypes->Fill(TVmmrMsg::MSG_EOE);

      // now check against global counter of all subsystems:
      if (fNectarRawEvent->fSequenceNumber < 0)
      {
        fNectarRawEvent->fSequenceNumber = count;    // init if not yet set by other board
	printdeb("Warning: VMMR board %d setting sequence number to %d \n", modid, count);
        // TODO: maybe we could put MBS header counter to fSequenceNumber, but this may differ by special trigger events
      }
      else
      {
        if (fNectarRawEvent->fSequenceNumber != count)    // check against global counter
        {
          printdeb("Warning: VMMR board %d eventcounter %d is not matching common counter %d \n", modid, count,
              fNectarRawEvent->fSequenceNumber);
          // JAM - optionally we may stop here in such case?
          //GO4_STOP_ANALYSIS_MESSAGE("Warning: VMMR board %d eventcounter %d is not matching common counter %d \n", modid, count, fNectarRawEvent->fSequenceNumber);
        }
      }

    }
    else if (isDummy)
    {
      // do nothing, just proceed.
    }
    else
    {
      GO4_STOP_ANALYSIS_MESSAGE("Data error: unknown format in VMMR word 0x%x from module id %d-  stopped!", word,
          modid);
      return kFALSE;    // redundant, above will do exception
    }

    // following check will avoid that sublen of frontend is wrong..
    NECTAR_EVENT_CHECK_PDATA
  }    // while ((pData - pdatastart) < sublen)
  return kTRUE;
}

Bool_t TNectarRawProc::UnpackV775(UShort_t triggertype)
{
  //read header
  Int_t V775_head = *pData++;
  //Bool_t isHeader = ((V775_head >> 28) & 0xF) == 0x4;



  uint32_t count =(uint32_t)(V775_head & 0x3f00) >> 8;
#ifdef NECTAR_VERBOSE_PRINT
  uint32_t crate =(uint32_t)(V775_head & 0xff0000) >> 16;
  uint32_t geom  =(uint32_t)(V775_head & 0xf8000000) >> 27;
#endif
  
//  std::cout << "UnpackV775: count: " << count << " crate: " << crate << " geom: " << geom << std::endl;
  printdeb("UnpackV775: count:%d crate %d geom:%d \n", count, crate, geom)

  uint32_t ch_val, channel;
#ifdef NECTAR_VERBOSE_PRINT
  uint32_t ch_valid; //ch_geom, ch_overflow, ch_underflow;
#endif
  uint32_t word;
  uint32_t counter = 0;
  // loop over following data words:
  while (counter < count)
    {
      counter++;
      //      std::cout << "Data: " << pData << std::endl;
      word = *pData++;
      ch_val       = (uint32_t)(word & 0x00000fff);
#ifdef NECTAR_VERBOSE_PRINT
      //ch_overflow  = (uint32_t)(word & 0x1000)       >> 12;
      //ch_underflow = (uint32_t)(word & 0x2000)       >> 13;

      ch_valid     = (uint32_t)(word & 0x4000)       >> 14;
#endif
      channel      = (uint32_t)(word & 0x1f0000)     >> 16;
      //ch_geom      = (uint32_t)(word & 0xf8000000)   >> 27;
//      std::cout << "ch: " << channel << " val: " << ch_val << " ok: " << ch_valid << std::endl;
      printdeb("775 ch:%d, value:0x%x, valid:%d\n", channel, ch_val, ch_valid)
      hRawV775_TDC[channel]->Fill(ch_val);
      switch (channel)
      {
        case 0:
          if (triggertype == 1)
            h_t_Xray[1]->Fill(ch_val);
          else
            h_t_Xray_OFF[1]->Fill(ch_val);
          break;
        case 1:
          if (triggertype == 1)
            h_t_Xray[2]->Fill(ch_val);
          else
            h_t_Xray_OFF[2]->Fill(ch_val);
          break;

        default:
          printdeb("Warning: found wrong channel number %d in V775!\n", channel);
          break;

      }
    }
  
  pData++; //skip EOB
  return kTRUE;
}

Bool_t TNectarRawProc::UnpackV785(UShort_t triggertype)
{
  uint32_t counter = 0;
  //read header
  Int_t V785_head = *pData++;

  uint32_t count = (uint32_t) (V785_head & 0x3f00) >> 8;
  //uint32_t crate = (uint32_t) (V785_head & 0xff0000) >> 16;
  //uint32_t geom = (uint32_t) (V785_head & 0xf8000000) >> 27;

  //std::cout << "count: " << count << " crate: " << crate << " geom: " << geom << std::endl;

  uint32_t ch_val, channel;
  //uint32_t ch_geom, ch_overflow, ch_underflow, ch_valid;
  uint32_t word;
  // loop over following data words:  
  while (counter < count)
  {
    counter++;
    //      std::cout << "Data: " << pData << std::endl;
    word = *pData++;
    ch_val = (uint32_t) (word & 0x00000fff);
    //ch_overflow = (uint32_t) (word & 0x1000) >> 12;
    //ch_underflow = (uint32_t) (word & 0x2000) >> 13;
    //ch_valid = (uint32_t) (word & 0x4000) >> 14;
    channel = (uint32_t) (word & 0x1f0000) >> 16;
    //ch_geom = (uint32_t) (word & 0xf8000000) >> 27;
    //std::cout << "ch: " << channel << " val: " << ch_val << " ok: " << ch_valid << std::endl;
    //hRawV785_ADC[channel]->Fill(ch_val);

    if (!fPar->fUseSetup2024)
    {
      switch (channel)
      {
        case 0:
          hRawV785_ADC[1]->Fill(ch_val);
          if (triggertype == 1)
            h_E_Xray[1]->Fill(ch_val * CAL_XRAY90_A + CAL_XRAY90_B);
          else
            h_E_Xray_OFF[1]->Fill(ch_val * CAL_XRAY90_A + CAL_XRAY90_B);
          break;
        case 1:
          hRawV785_ADC[2]->Fill(ch_val);
          if (triggertype == 1)
            h_E_Xray[2]->Fill(ch_val * CAL_XRAY145_A + CAL_XRAY145_B);
          else
            h_E_Xray_OFF[2]->Fill(ch_val * CAL_XRAY145_A + CAL_XRAY145_B);
          break;
        default:
          //printf("Warning: found wrong channel number %d in V785!\n",channel);
          break;
      }
    }
    else
    {
      //printf("Warning: found V785 data although new setup is enabled! Please check fUseSetup2024 in parameter\n");
      GO4_STOP_ANALYSIS_MESSAGE("Stopped: found V785 data although new setup is enabled! Please check fUseSetup2024 in parameter\n");
    }

  }
  
  pData++;    //skip EOB
  
  return kTRUE;
}

Bool_t TNectarRawProc::UnpackV830()
{
  uint32_t counter = 0;
  //read header
  Int_t V830_head = *pData++;

  uint32_t count =(uint32_t)(V830_head & 0xfc0000) >> 18;
  //uint32_t geom  =(uint32_t)(V830_head & 0xf8000000) >> 27;

  //std::cout << "count: " << count << " geom: " << geom << std::endl;

  uint32_t ch_val, channel, word; //, ch_diff;
  
  // loop over following data words:
  while (counter < count)
    {
      counter++;
      //      std::cout << "Data: " << pData << std::endl;
      word = *pData++;

      //26 bit mode
      //ch_val       = (uint32_t)(word & 0x3ffffff);
      //channel      = (uint32_t)(word & 0xf8000000) >> 27;
            
      //32 bit mode
      ch_val       = (uint32_t)(word & 0xffffffff);
      channel      = counter;
      if(!fPar->fUseSetup2024)
      {
        V830_diff[channel] = ch_val - V830_old[channel];

        //std::cout << "ch: " << channel << " diff: " << V830_diff[channel] << std::endl;
        V830_old[channel] = ch_val;
      }
      else
      {
          //printf("Warning: found V830 data although new setup is enabled! Please check fUseSetup2024 in parameter\n");
          GO4_STOP_ANALYSIS_MESSAGE("Stopped: found V830 data although new setup is enabled! Please check fUseSetup2024 in parameter\n");
      }
      NECTAR_EVENT_CHECK_PDATA
    }

  return kTRUE;
}


Bool_t TNectarRawProc::UnpackTPAT(Int_t header)
{
  uint32_t counter = 0;
  uint32_t count =(uint32_t)(header   & 0xfff);
  //uint32_t ev_num  =(uint32_t)(header & 0xfff000) >> 12;
  //uint32_t id  =(uint32_t)(header &     0xff000000) >> 24;

  //std::cout << "count: " << count << " ev: " << ev_num << " id " << id << std::endl;
  
  uint32_t ts_lo, ts_hi, tpat, trig, word;

  while (counter < count/3)
    {
      counter++;
      //std::cout << "header: " <<  << std::endl;
      ts_lo = *pData++;
      //std::cout << "ts_lo: " << ts_lo << std::endl;
      ts_hi = *pData++;
      //std::cout << "ts_hi: " << ts_hi << std::endl;
      if (TRLO_TS_old == 0) TRLO_TS_start = ((uint64_t) ts_hi << 32 ) | ts_lo;
      
      TRLO_TS = ((uint64_t) ts_hi << 32 ) | ts_lo;

      if ( (TRLO_TS - TRLO_TS_start - (uint64_t)SCAoffset*100000000) >= (uint64_t)SCALER_RANGE*100000000 )
	{
	  SCAoffset+=SCALER_RANGE;
	  ResetScalers();
	  //std::cout << "SCA offset: " << SCAoffset << std::endl;
	}
      TRLO_runtime = (((double_t)TRLO_TS - (double_t)TRLO_TS_start)/100000000.);
      SCAtime = TRLO_runtime - (double)SCAoffset;
      TRLO_TS_old = TRLO_TS;
      //FillBeamScalers();
      
      word = *pData++;
      tpat = (word & 0xffff);
      trig = (word & 0xf000000) >> 24;
      //std::cout << "----- tpat: " << tpat << " trig: " << trig << " -----" << std::endl;
      hRawTRLO_tpat->Fill(tpat);
      hRawTRLO_trigger->Fill(trig);
      NECTAR_EVENT_CHECK_PDATA
    }
  return kTRUE;
}

Bool_t TNectarRawProc::UnpackVSCA(Int_t header)
{
  uint32_t counter = 0;
  
  uint32_t b_lmu =(uint32_t)(header   & 0x3f);
  uint32_t b_lmu_mux  = (uint32_t)(header & 0x7c0) >> 6;
  uint32_t b_lmu_aux  = (uint32_t)(header & 0xf800) >> 11;
  uint32_t a_lmu  = (uint32_t)(header & 0xff0000) >> 16;
  //uint32_t id  =(uint32_t)(header &     0xff000000) >> 24;

  //std::cout << "b_lmu: " << b_lmu << " mux: " << b_lmu_mux << " aux: " << b_lmu_aux << " a_lmu: " << a_lmu << " id: " << id << std::endl;
  
  uint32_t before_lmu[b_lmu];
  uint32_t d_blmu;	
  while (counter < b_lmu)
    {
      counter++;
      before_lmu[counter] = *pData++;

      d_blmu = before_lmu[counter]-blmu_old[counter];
      if ( blmu_old[counter] > 0 && d_blmu > 0)
        {
          //printf( "counter: %d\n", counter);
          switch (counter)
            {
            case 1:
              h_Tel_blmu->Fill(SCAtime, d_blmu);
              break;
            case 2:
              h_HRes_blmu->Fill(SCAtime, d_blmu);
              break;
            case 3:
              h_Xray1_blmu->Fill(SCAtime, d_blmu);
              break;
            case 4:
              h_Xray2_blmu->Fill(SCAtime, d_blmu);
              break;
            case 5:
              h_Xray3_blmu->Fill(SCAtime, d_blmu);
              break;
            default:
              break;
            }
        }
      //printf("blmu_Data ch %d: %x\n",counter, before_lmu[counter]);
      blmu_old[counter] = before_lmu[counter];
      NECTAR_EVENT_CHECK_PDATA
    }
  counter = 0;
  // JAM 5-24: just skip the non used trloii scalers, suppress warnings
  //uint32_t before_lmu_mux[b_lmu_mux];
  while (counter < b_lmu_mux)
    {
      counter++;
      //before_lmu_mux[counter] = *pData++;
      pData++;
      //printf("mux_Data ch %d: %x\n", counter, before_lmu_mux[counter] );
      NECTAR_EVENT_CHECK_PDATA
    }
  counter = 0;
  // JAM 5-24: just skip the non used trloii scalers, suppress warnings
  //uint32_t before_lmu_aux[b_lmu_aux];
  while (counter < b_lmu_aux)
    {
      counter++;
      //before_lmu_aux[counter] = *pData++;
      pData++;
      //printf("aux_Data ch %d: %x\n",counter,before_lmu_aux[counter] );
      NECTAR_EVENT_CHECK_PDATA
    }

  counter = 0;

  uint32_t before_dt[a_lmu];
  uint32_t after_dt[a_lmu];
  uint32_t after_red[a_lmu];

  uint32_t d_bdt = 0;
  uint32_t d_adt = 0;
  uint32_t d_ared = 0;
  
  while (counter < a_lmu)
    {
      counter++;
      before_dt[counter] = *pData++;
      after_dt[counter] = *pData++;
      after_red[counter] = *pData++;
      d_bdt = before_dt[counter]-bdt_old[counter];
      d_adt = after_dt[counter]-adt_old[counter];
      d_ared = after_red[counter]-ared_old[counter];
      if ( bdt_old[counter] != 0 && d_bdt > 0)
        {
          //printf( "counter: %d\n", counter);
          switch (counter)
            {
            case 1:
            case 4:
              //printf( "filling Tel");
              h_Tel_bdt->Fill(SCAtime, d_bdt);
              h_Tel_adt->Fill(SCAtime, d_adt);
              h_Tel_ared->Fill(SCAtime, d_ared);
              break;
            case 2:
            case 5:
              //printf( "filling HRes");
              h_HRes_bdt->Fill(SCAtime, d_bdt);
              h_HRes_adt->Fill(SCAtime, d_adt);
              h_HRes_ared->Fill(SCAtime, d_ared);
              break;
            case 3:
            case 6:
              //printf( "filling Xray");
              h_Xray_bdt->Fill(SCAtime, d_bdt);
              h_Xray_adt->Fill(SCAtime, d_adt);
              h_Xray_ared->Fill(SCAtime, d_ared);
              break;
            }

        }
      //printf("ch %d BDT: %x\t ADT: %x\t ARED: %x\n",counter, before_dt[counter], after_dt[counter], after_red[counter] );
	bdt_old[counter] = before_dt[counter];
	adt_old[counter] = after_dt[counter];
	ared_old[counter] = after_red[counter];
	NECTAR_EVENT_CHECK_PDATA
    }
	
  return kTRUE;
}

Bool_t TNectarRawProc::UnpackPlainVSCA()
{
  for (Int_t channel = 0; channel < VULOM_NUMSCALERS; ++channel)
  {
    uint32_t ch_val = (uint32_t) ((*pData++) & 0xFFFFFFFF);
    h_vulom_raw->SetBinContent(channel+1,ch_val);

    if (fPar->fUseSetup2024)
    {
      V830_diff[channel] = ch_val - V830_old[channel];
      //std::cout << "UnpackPlainVSCA: ch: " << channel << " diff: " << V830_diff[channel] <<", value: "<< ch_val << std::endl;
      printdeb("UnpackPlainVSCA: ch: %d diff:%d , value:%d\n",channel, V830_diff[channel], ch_val);
      V830_old[channel] = ch_val;
    }
    else
    {
      //printf("Warning: found VULOM data although old setup is enabled! Please check fUseSetup2024 in parameter\n");
      GO4_STOP_ANALYSIS_MESSAGE("Stopped: found VULOM data although old setup is enabled! Please check fUseSetup2024 in parameter\n");
    }
    NECTAR_EVENT_CHECK_PDATA
  }      // for

  return kTRUE;
}

Bool_t TNectarRawProc::UnpackMadc(UShort_t triggertype)
{

  uint32_t counter = 0;
  uint32_t ch_val, channel, ch_overflow;    //data
  uint32_t sig = (uint32_t) (*pData & 0xC0000000) >> 30;
  uint32_t sub = (uint32_t) (*pData & 0x3F000000) >> 24;
#ifdef NECTAR_VERBOSE_PRINT
  uint32_t geom       = (uint32_t)(*pData & 0x00FF0000) >> 16;
  //uint32_t format   = (uint32_t)(*pData & 0x00008000) >> 15;
  uint32_t adc_resol  = (uint32_t)(*pData & 0x00007C00) >> 10;
#endif
  uint32_t word_num = (uint32_t) (*pData & 0x000003FF);    // this gives the # of following words

  //check module header
  if ((sig == 1) && (sub == 0))
  {
    printdeb("header MADC: 0x%x\n", *pData);
    printdeb("adc=%d; geom=%d; word_num=%d\n", adc_resol, geom, word_num);
  }
  else
  {
    NECTAR_SKIP_EVENT("expected madc32 header not found: sig=0x%x sub=0x%x ! Skipping event!\n", sig, sub);
    //return kFALSE;
  }

  pData++;

  while (counter < word_num)
  {
    sig = (uint32_t) (*pData & 0xC0000000) >> 30;
    sub = (uint32_t) (*pData & 0x3FE00000) >> 21;

    printdeb("MADC32: sig=%d; sub=%d\n", sig, sub);

    //dummy - data - ts
    if (sig == 0)
    {

      if (sub == 0)    //dummy -> skip
      {
        pData++;
      }
      else if (sub == 32)    //data -> read
      {
        ch_val = (uint32_t) (*pData & 0x00001FFF);    //channel value
        ch_overflow = (uint32_t) (*pData & 0x00004000) >> 14;
        channel = (uint32_t) (*pData & 0x001F0000) >> 16;

        printdeb("MADC32: channel=%d; ch_val=%d; overflow=%d\n", channel, ch_val, ch_overflow);

        if (ch_overflow)
            ch_val = MADC_RANGE;// -1;    //fill overflow bin if channel in overlow
          hRaw_MADC[channel]->Fill(ch_val);
        //////////////////
        if (fPar->fUseSetup2024)
           {
             switch (channel)
             {
               case 0:
                 hRawV785_ADC[1]->Fill(ch_val);
                 if (triggertype == 1)
                   h_E_Xray[1]->Fill(ch_val * CAL_XRAY90_A + CAL_XRAY90_B);
                 else
                   h_E_Xray_OFF[1]->Fill(ch_val * CAL_XRAY90_A + CAL_XRAY90_B);
                 break;
               case 1:
                 hRawV785_ADC[2]->Fill(ch_val);
                 if (triggertype == 1)
                   h_E_Xray[2]->Fill(ch_val * CAL_XRAY145_A + CAL_XRAY145_B);
                 else
                   h_E_Xray_OFF[2]->Fill(ch_val * CAL_XRAY145_A + CAL_XRAY145_B);
                 break;
               default:
                 //printf("Warning: found wrong channel number %d in MADC!\n",channel);
                 break;
             }
           }
           else
           {
             //printf("Warning: found MADC data although old setup is enabled! Please check fUseSetup2024 in parameter\n");
             GO4_STOP_ANALYSIS_MESSAGE("Stopped: found MADC data although old setup is enabled! Please check fUseSetup2024 in parameter.\n");
           }



        /////////////////

        pData++;
      }
      else if (sub == 36)    // ext_ts -> read
      {
        //_ext_ts       = (uint32_t)(*(*p_data) & 0x0000FFFF); //timestamp
        //ext_ts = _ext_ts;
        pData++;
      }
      else
      {
        NECTAR_SKIP_EVENT("MADC32: sub signature of word not identified -> skipping event");
        //return kFALSE;
      }
    }
    else if (sig == 3)    //event counter
    {
      //evt_counter = (uint32_t) (*pData++ & 0x3FFFFFFF);
      pData++;

    }
    else
    {
      NECTAR_SKIP_EVENT("MADC32: signature of word not identified -> skipping event");
      //return kFALSE;
    }

    counter++;
    NECTAR_EVENT_CHECK_PDATA
  }// while







   return kTRUE;
 }



void TNectarRawProc::ScalerReset(TH1* histo)
{
  if ( histo->GetEntries() > 0) histo->Reset();
}

void TNectarRawProc::ResetScalers()
{
  h_Tel_blmu->Reset();
  ScalerReset(h_Tel_bdt);
  ScalerReset(h_Tel_adt);
  h_Tel_ared->Reset();
  h_HRes_blmu->Reset();
  h_HRes_bdt->Reset();
  h_HRes_adt->Reset();
  h_HRes_ared->Reset();
  h_Xray1_blmu->Reset();
  h_Xray2_blmu->Reset();
  h_Xray3_blmu->Reset();
  ScalerReset(h_Xray_bdt);
  ScalerReset(h_Xray_adt);
  ScalerReset(h_Xray_ared);
  ScalerReset(h_trafo);
  h_I_cooler->Reset();
  h_U_cooler->Reset();
  h_jet_S1->Reset();
  h_jet_S2->Reset();
  h_pmt->Reset();
  
  h_vulom_raw->Reset();

}

void TNectarRawProc::FillBeamScalers()
{

  for (int i=0; i<32; i++)
    {
      if ( V830_old[i] > 0 && V830_diff[i] > 0 )
	{
	  switch (i)
	    {
	    case 0: h_trafo->Fill(SCAtime, V830_diff[i]);
	    case 1: h_I_cooler->Fill(SCAtime, V830_diff[i]);
	    case 2: h_U_cooler->Fill(SCAtime, V830_diff[i]);
	    case 3: h_jet_S1->Fill(SCAtime, V830_diff[i]);
	    case 4: h_jet_S2->Fill(SCAtime, V830_diff[i]);
	    case 5: h_pmt->Fill(SCAtime, V830_diff[i]);
	    }
	}
    }
}

Bool_t TNectarRawProc::UpdateDisplays()
{

// place for some advanced analysis from output event data here, if not in the second go4 analysis step...
  
  
  // access now collected data of mdpp boards to evaluate delta t histograms:
  for (UInt_t i = 0; i < TNectarRawEvent::fgConfigMdppBoards.size(); ++i)
  {
    UInt_t uniqueid = TNectarRawEvent::fgConfigMdppBoards[i];

    TMdppBoard* theMDPP = fNectarRawEvent->GetMdppBoard(uniqueid);
    if (theMDPP == 0)
    {
      GO4_STOP_ANALYSIS_MESSAGE(
          "Configuration error from UpdateDisplays: MDPP module id %d does not exist as subevent. Please check TNectarRawParam setup",
          uniqueid);
      return kFALSE;
    }

//////////////////// JAM check if there is any data:
#ifdef  NECTAR_VERBOSE_PRINT
    for(UChar_t c=0; c<MDPP_CHANNELS; ++c)
    {
      UInt_t maxmessages=theMDPP->NumTdcMessages(c);
      printdeb("Channel %d has %d maxmessages\n", c, maxmessages);
      for(UInt_t i=0; i<maxmessages; ++i)
      {
        TMdppTdcData* msg= theMDPP->GetTdcMessage(c, i);
        if(msg)
        {
          printdeb(" - data in channel %d (i=%d)  is %d\n", c, i, msg->fData);
        }
      }
    }

#endif

    TMdppDisplay* boardDisplay = GetMdppDisplay(uniqueid);
    if (boardDisplay == 0)
    {
      GO4_STOP_ANALYSIS_MESSAGE(
          "Configuration error from UpdateDisplays: MDPP module id %d does not exist as histogram display set!",
          uniqueid);
      return kFALSE;
    }
    Int_t refchannel = fPar->fMDPP_ReferenceChannel[i];
    if (refchannel < 0 || refchannel >= MDPP_CHANNELS)
      refchannel = 0;

    printdeb("UpdateDisplays - refchannel %d for board %d, index %d\n", refchannel, uniqueid, i);

    // first get all reference time messages of this event:
    UInt_t maxrefs = theMDPP->NumTdcMessages(refchannel);
    printdeb("UpdateDisplays - maxrefs=%d\n", maxrefs);

    if (maxrefs == 0)
      continue;
    Int_t reftime[maxrefs];
    for (UInt_t r = 0; r < maxrefs; ++r)
    {
      TMdppTdcData* refdata = theMDPP->GetTdcMessage(refchannel, r);
      if (refdata)
        reftime[r] = refdata->fData;

      printdeb(" - reftime[%d]=%d\n", r, reftime[r]);
    }
    for (UChar_t c = 0; c < MDPP_CHANNELS; ++c)
    {
      UInt_t maxmessages = theMDPP->NumTdcMessages(c);
      printdeb("Channel %d has %d maxmessages\n", c, maxmessages);
      for (UInt_t i = 0; i < maxmessages; ++i)
      {
        TMdppTdcData* msg = theMDPP->GetTdcMessage(c, i);
        if (msg)
        {
          UInt_t r = i;    // always use reference time message of same buffer index
          if (r && r >= maxrefs)
            r = maxrefs - 1;    // (...if existing, otherwise use most recent one)
          Int_t deltaT = msg->fData - reftime[r];
          printdeb(" - deltaT channel %d (i=%d=  is %d\n", c, i, deltaT);

          boardDisplay->hDeltaTDC[c]->Fill(deltaT);
        }
      }
    }
  }
    
    //////////////////////////////////////////////
  /* for (UInt_t i = 0; i < TNectarRawEvent::fgConfigVmmrBoards.size(); ++i)
  {
    UInt_t uniqueid = TNectarRawEvent::fgConfigVmmrBoards[i];

    TVmmrBoard* theVmmr = fNectarRawEvent->GetVmmrBoard(uniqueid);
    if (theVmmr == 0)
    {
      GO4_STOP_ANALYSIS_MESSAGE(
          "Configuration error from UpdateDisplays: VMMR module id %d does not exist as subevent. Please check TNectarRawParam setup",
          uniqueid);
      return kFALSE;
    }
    //uniqueid=0;
    TVmmrDisplay* boardDisplay = GetVmmrDisplay(uniqueid);
    if (boardDisplay == 0)
    {
      GO4_STOP_ANALYSIS_MESSAGE(
          "Configuration error from UpdateDisplays: VMMR module id %d does not exist as histogram display set!",
          uniqueid);
      return kFALSE;
    }

    for (int slid=0; slid<VMMR_CHAINS; ++slid)
    {

    TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE); // will return 0 if not existing, do not create new element!
     if (theslave==0)
    {
      continue; //  this means there is just empty slot in event
    }
    
    TVmmrSlaveDisplay* sldisplay= boardDisplay->GetSlaveDisplay(slid);
    if (sldisplay == 0)
    {
      GO4_STOP_ANALYSIS_MESSAGE(
          "Configuration error from UpdateDisplays: VMMR slave id %d does not exist as histogram display set!", slid);
      return kFALSE;
    }
    
    // now get channel data from event structure for channels 1 and 2:
    
     UInt_t maxvmmrmessages = theslave->NumAdcMessages();
      printdeb("slave %d has %d maxmessages\n", slid,  maxvmmrmessages);
       //UShort_t ch1[10];
       //UShort_t ch2[10];
       //Modifications - cobus and michele 7 feb 22
       UShort_t chS[64][10];
       //UShort_t chH[16][10];
       
       int ix1A[64]={0};
       int ix1=0;
       
       
     for(int a=0;a<64;a++){ 
     ix1=0; 
      for (UInt_t i = 0; i < maxvmmrmessages; ++i)
      {
          TVmmrAdcData* dat = theslave->GetAdcMessage(i);
           if(dat->fFrontEndSubaddress==a)
           {
                chS[a][ix1++]=dat->fData;
                if(ix1>+10) break;
            }
      }
      ix1A[a]=ix1;
     } 
     //Filling New histogram :: hRawADCVvsH (1 ch Vertical vs all chs horizontal) 
     for(int v=16;v<32;v++){
      for(int a=0;a<16;a++){
      int maxi=ix1A[v];
      if(ix1A[a]<ix1A[v])maxi=ix1A[a];      
      for(int i=0; i<maxi;i++)
      {
            sldisplay->hRawADCVvsH[v-16]->Fill(chS[a][i],chS[v][i]);
      	}
      }	
   }  
     //Filling New histogram :: hRawADCHvsV (1 ch Horizontal vs all chs vertical) 
    for(int v=0;v<16;v++){
      for(int a=16;a<32;a++){
      int maxi=ix1A[v];
      if(ix1A[a]<ix1A[v])maxi=ix1A[a];      
      for(int i=0; i<maxi;i++)
      {
            sldisplay->hRawADCHvsV[v]->Fill(chS[a][i],chS[v][i]);
      	}
      }	
   }  
    
   int RstripP[32] = {15,11,16,14,12,9,13,8,10,5,6,3,7,2,4,1,15,16,10,14,13,12,9,11,4,8,7,6,3,5,1,2};    
     //Filling new histogram :: hSTRIPmap (BB8 MMR to strips conversion)
    for(int h=0;h<16;h++){
      for(int v=16;v<32;v++){
     int maxi = ix1A[h];
      if(ix1A[v]<ix1A[h])maxi=ix1A[v];
     for(int i=0; i<maxi;i++){  
    if((chS[h][i]>230)&&(chS[v][i]>230)){
    sldisplay->hSTRIPmap->Fill(RstripP[h],RstripP[v]);
   }
   }
    }
    }
  }		
  } */	
  return kTRUE;	
}	

