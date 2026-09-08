#include "TNectarProfileProc.h"

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
#include "TNectarRawProc.h"
#include "TNectarProfileEvent.h"


#include "TGo4Fitter.h"
#include "TGo4FitDataHistogram.h"
#include "TGo4FitParameter.h"
#include "TGo4FitModelPolynom.h"

#include "TGo4UserException.h"


//*********************************************************
TNectarProfileProc::TNectarProfileProc() :
  TGo4EventProcessor(), fNectarProfileEvent(nullptr)
{
  TGo4Log::Info("TNectarProfileProc: Create");

  // BuildNEvent();
}
//
//
//*********************************************************
TNectarProfileProc::TNectarProfileProc(const char* name) :
  TGo4EventProcessor(name), fNectarProfileEvent(nullptr)
{

  TGo4Log::Info("TNectarProfileProc: Create1");
  ////init user analysis object;

  //// init user analysis objects:
  //fPar1 = dynamic_cast<TNectarRawParam*>(MakeParameter("NectarRawParam", "TNectarRawParam", "set_NectarRawParam.C"));
  //if (fPar1)
  //xfPar1->SetConfigBoards();    // board configuration specified in parameter is used to build subevents and displays

  // Creation of histograms (or take them from autosave)
     TString obname;
     TString obtitle;
     TString foldername;
  
   //******************** 
   //Init BB8 Histograms
   //*******************

     for(int i=0;i<16;i++){
       obname.Form("Profile/BB8/ADC_Vertical_Strip_%d", i+1);
       obtitle.Form("ADC_Vertical_Strip_%d", i+1);
       histoVStrip[i]= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");
     } 
     for(int i=0;i<16;i++){
       obname.Form("Profile/BB8/ADC_Horizontal_Strip_%d", i+1);
       obtitle.Form("ADC_Horizontal_Strip_%d", i+1);  
       histoHStrip[i]= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");
     } 

   //******************** 
   //Init MSXO4-1 & 2  Histograms
   //*******************
       obname.Form("Profile/MSX04/ADC_MSX04_1");
       obtitle.Form("ADC_MSX04_1");
       h_ADC_MSX04_1= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");
      
       obname.Form("Profile/MSX04/ADC_MSX04_2");
       obtitle.Form("ADC_MSX04_2");  
       h_ADC_MSX04_2= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");  
     
   //******************** 
   //Init MSXO4-7 & 8  Histograms
   //*******************

     obname.Form("Profile/MSX04/ADC_MSX04_7");
     obtitle.Form("ADC_MSX04_7");
     h_ADC_MSX04_7= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");
     
     obname.Form("Profile/MSX04/ADC_MSX04_8");
     obtitle.Form("ADC_MSX04_8");  
     h_ADC_MSX04_8= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");  
     
     //******************** 
     //Init MSXO4-9 & 10   Histograms
     //*******************

     obname.Form("Profile/MSX04/ADC_MSX04_9");
     obtitle.Form("ADC_MSX04_9");
     h_ADC_MSX04_9= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");
     
     obname.Form("Profile/MSX04/ADC_MSX04_10");
     obtitle.Form("ADC_MSX04_10");  
     h_ADC_MSX04_10= MakeTH1('I', obname.Data(), obtitle.Data(), VMMR_ADC_RANGE, 0, VMMR_ADC_RANGE, "ADC value", "counts");  
     

   //******************** 
   //Init BB29  Histograms
   //*******************
     
}

//********************************************************

Bool_t TNectarProfileProc::BuildEvent(TGo4EventElement *target)
{
  //TNectarBB8* DSSD = (TNectarBB8*) target;
  TNectarTelescope* Telescope = (TNectarTelescope*) target;
  TNectarRawEvent* rawEvent = (TNectarRawEvent*) GetInputEvent(); 
  if(!rawEvent->IsValid()) return kTRUE; // ignore invalid 1st step JAM 2026
  //Bool_t fOutput = kFALSE; //not store the output
  Double_t ADC1=0,ADC2=0,ADC3=0,ADC4=0,ADC5=0,ADC6=0;
  for (UInt_t i = 0; i < 1; ++i)
  {
    UInt_t uniqueid = rawEvent->fgConfigVmmrBoards[i];
    TVmmrBoard* theVmmr = rawEvent->GetVmmrBoard(uniqueid);
    
    for (int slid=1; slid<2; ++slid)
    {
      TVmmrSlave* theslave=theVmmr->GetSlave(slid, kFALSE);
      UInt_t maxvmmrmessages = theslave->NumAdcMessages();
      TGo4Log::Info("maxvmmrmessages - %d",maxvmmrmessages);
          
      //For cycle around the strips values to fill new variables 
      //maxvmmrmessages;
      //Int_t ix1=0;
      //UInt_t maxvmmrmessages = theslave->NumAdcMessages();
      int Vmmr_strip_C[32] = {15,11,16,14,12,9,13,8,10,5,6,3,7,2,4,1,15,16,10,14,13,12,9,11,4,8,7,6,3,5,1,2};    //Convension Array VMMR channel -> Strips number 
      for(UInt_t a=0;a<64;a++){    
        for (UInt_t i = 0; i < maxvmmrmessages; ++i)
        {
          //TGo4Log::Info("1");
          TVmmrAdcData* dat = theslave->GetAdcMessage(i);
          //TGo4Log::Info("2");
          if(dat->fFrontEndSubaddress==a){
            if(a<=15){  // Filling BB8 ProfileEvent with Convertion (VMMR channel -> Strips number) 
              Int_t id = Vmmr_strip_C[a]; 
              Telescope->HStrip[id] = dat->fData;
            }
          }
          if(dat->fFrontEndSubaddress==a){
            if ((a>15)&&(a<32))
            {
              Int_t id = Vmmr_strip_C[a];
              Telescope->VStrip[id] = dat->fData;
            }
          }
          // TGo4Log::Info("3");
          if ((dat->fFrontEndSubaddress==a)&&(a==42)) ADC1 = dat->fData;//MSX04->ADC_MSX04_1 = dat->fData;
          if ((dat->fFrontEndSubaddress==a)&&(a==43)) ADC2 = dat->fData;//MSX04->ADC_MSX04_2 = dat->fData; 
          if ((dat->fFrontEndSubaddress==a)&&(a==44)) ADC3 = dat->fData;//MSX04->ADC_MSX04_7 = dat->fData;
          if ((dat->fFrontEndSubaddress==a)&&(a==45)) ADC4 = dat->fData;//MSX04->ADC_MSX04_8 = dat->fData; 
          if ((dat->fFrontEndSubaddress==a)&&(a==46)) ADC5 = dat->fData;//MSX04->ADC_MSX04_9 = dat->fData;
          if ((dat->fFrontEndSubaddress==a)&&(a==47)) ADC6 = dat->fData;//MSX04->ADC_MSX04_10 = dat->fData; 
          //if(ix1>+10) break;
          
        }
      }   
      // TNectarMSX04* MSX04 = (TNectarMSX04*) target;
      
      Telescope->ADC_MSX04_1 = ADC1;
      Telescope->ADC_MSX04_2 = ADC2;
      Telescope->ADC_MSX04_7 = ADC3;
      Telescope->ADC_MSX04_8 = ADC4;
      Telescope->ADC_MSX04_9 = ADC5;
      Telescope->ADC_MSX04_10 = ADC6;
     
      //fNectarProfileEvent->SetValid(kTRUE);
      //MSX04->SetValid(kTRUE);
      //DSSD->SetValid(kTRUE);
      // Telescope->SetValid(kTRUE);
      
      Int_t Adc=0;
      //Fill the vertical Strips Histograms::  BB8
      for(Int_t i=0; i<16; i++){
        Adc=Telescope->VStrip[i+1];
        histoVStrip[i]->Fill(Adc);
      }
      //Fill the horizontal Strips Histograms:: BB8
      for(Int_t i=0; i<16; i++){
        Adc=Telescope->HStrip[i+1]; 
        histoHStrip[i]->Fill(Adc);
      }

  
      
     
      //Fill the ADC_histogram:  MSX04 Detectors
      // MSX04_1 & 2
        Adc =Telescope->ADC_MSX04_1;
      h_ADC_MSX04_1->Fill(Adc);
      Adc = Telescope->ADC_MSX04_2;
      h_ADC_MSX04_2->Fill(Adc);
      // MSX04_7 & 8
      Adc = Telescope->ADC_MSX04_7;
      h_ADC_MSX04_7->Fill(Adc);
      Adc =  Telescope->ADC_MSX04_8;
      h_ADC_MSX04_8->Fill(Adc);
      // MSX04_9 & 10
      Adc =  Telescope->ADC_MSX04_9;
      h_ADC_MSX04_9->Fill(Adc);
      Adc =  Telescope->ADC_MSX04_10;
      h_ADC_MSX04_10->Fill(Adc);     

      Telescope->SetValid(kTRUE);
    }
  }	
  return kTRUE;	
}	



