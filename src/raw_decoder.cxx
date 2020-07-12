#include <arpa/inet.h>
#include <assert.h>
#include <utility>

#include "raw_decoder.h"

#include <algorithm> 
#include <fstream>
#include <iostream>
#include <utility>

#include "TPaveText.h"
#include "TLatex.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>

using namespace std;

double sigmacut=5.0;
//==========================================================================
RawDecoder::RawDecoder(unsigned int * buffer, int n)
{ 
  mAPVRawSingleEvent.clear();
  mAPVRawHisto.clear();

  fBuf = n;
  buf = new unsigned int[fBuf];
  for(int i=0;i<fBuf;i++)
  {
    buf[i] = buffer[i];
  }
  Decode();
}

//==========================================================================
RawDecoder::RawDecoder(const vector<uint32_t> &buffer, int start, int end)
{
  mAPVRawSingleEvent.clear();
  mAPVRawHisto.clear();

  fBuf = end - start;
  buf = new unsigned int[fBuf];
  int bufp = 0;
  for(int i=start;i<end;i++)
    {

      buf[bufp] = buffer[i];
      bufp++;
    }
  if(bufp!=fBuf){cout<<"vector passed and vector doesnt match"<<endl;}
  if(fBuf <= 0)
    {
      cout<<"empty vector passed in..."<<endl;
      return;
    }
  
  //std::cout<<__FUNCTION__<<": Data size get : "<<  bufp<<std::endl;
  Decode();
  
};

//==========================================================================
RawDecoder::~RawDecoder()
{
  //free buf
  delete[] buf;

  //clear maps
  map<int, vector<int> >::iterator it;
  for(it=mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
      it->second.clear();
  }
  mAPVRawSingleEvent.clear();

  //clear APV raw histos
  map<int, TH1F*>::iterator it_apv;
  for(it_apv=mAPVRawHisto.begin(); it_apv!=mAPVRawHisto.end(); ++it_apv)
  {
      TH1F *h = it_apv->second;
      h->Delete();
  }
  mAPVRawHisto.clear();
  //
}

//==========================================================================
void RawDecoder::Decode()
{
  unsigned int word32bit;

  int mpdid;
  int adc_ch;
  int hybridID;

  //map<int, vector<int> > mpd_event;
  vector<int> apv_event;
  vector<int> apv_margin;
  vector<int> mpd_margin;

  mpd_margin.clear();
  apv_margin.clear();
  
  //find MPD margin, find apv margin

#define BLOCK_HEADER    0x0	// {3'h0, MODULE_ID, EVENT_PER_BLOCK, BLOCK_CNT[7:0]}
#define BLOCK_TRAILER   0x1	// {3'h1, 1'b0, BlockWordCounter}
#define EVENT_HEADER    0x2	// {3'h2, 1'b0, EventCounterFifo_Data}
#define TRIGGER_TIME   0x3	// {3'h3, 1'b0, TimeCounterFifo_Data[39:20]}
//#define TRIGGER_TIME2   0x3	// {3'h3, 1'b1, TimeCounterFifo_Data[19:0]}
#define APV_CH_DATA     0x4	// {3'h4, ChannelData[20:0]}
#define EVENT_TRAILER   0x5	// {3'h5, 1'b0, LoopDataCounter[11:0], TRIGGER_TIME_FIFO}
#define DATA_NOT_VALID  0x6	// {3'h6, 21'b0}
#define FILLER_WORD     0x7	// {3'h7, 21'b0}

  //find apv margin
    for(int i=0;i<fBuf;i++)//skipping ssp data
  {
    uint32_t data = buf[i];
    uint32_t header;
    uint32_t apv_header;
    header = (data & 0x00e00000)>>21; 
    switch(header)
      {
      case BLOCK_HEADER:
	mpdid=(data&0x001F0000) >> 16;
	//	if(mpdid == 9)	cout<<"MPDID: "<<mpdid<<endl;
	break;
      case EVENT_HEADER:
	break;
      case TRIGGER_TIME:
	break;
      case APV_CH_DATA:
	switch((data& 0x00180000)>>19)
	  {
	  case 0: //apv header  -------  {1'b0, 1'b0, 1'b0, MEAN[11], DATA_IN[12:0], CH_ID[3:0]};
	    adc_ch=(data&0xf);
	    hybridID=(mpdid<<12)|(adc_ch<<8);
	    //   cout<<"adc_ch: "<<adc_ch<<endl;
	    //cout<<"Update apv header"<<std::endl;
	    break;
	  case 1: //data  -------  {1'b0, 1'b1, THRESHOLD_ADDRESS[6:0], data_minus_baseline[11:0]};
	    mAPVRawSingleEvent[hybridID].push_back(data & 0x00000fff);
	    //if((mpdid==9)&&(adc_ch==13))
	    //cout<<" sample count 1: "<<"  MPD "<< mpdid<<"  apvID:"<<adc_ch <<"  hybridID:"<<hybridID<<"  adc:"<<(data & 0x00000fff)<<"  Data:"<<std::hex<<data<<std::dec<<std::endl;
	    break;
	  case 2: //apv trailer  -------  {1'b1, 1'b0, 2'b0, MODULE_ID[4:0], DATA_IN[11:0]};
	                                                                   //DATA_IN[11:0] = {ApvSampleCounterMinusOne[3:0], frame_counter[7:0]};
	    mAPVRawSingleEvent[hybridID].push_back((data&0xf00)>>8);

	   // 	    cout<<" sample count 2: "<<"  MPD "<< mpdid<<"  apvID:"<<adc_ch <<"  "<<((data&0xf00)>>8)<<endl<<endl<<endl;
	    break;
	  case 3: //Trailer  -------  {1'b1, 1'b1, MEAN[10:0], word_count[7:0]};
	    break;
	  default:
	    break;
	  }
	break;
      case EVENT_TRAILER:
	break;
      case BLOCK_TRAILER:
	break;
      case DATA_NOT_VALID:
	break;
      case FILLER_WORD:
	break;
      default:
	break;
      }

  }
}

//==========================================================================
map<int, vector<int> > RawDecoder::GetDecoded()
{
  return mAPVRawSingleEvent;
}


//=============================================================
// ARGS:
//      C; the canvas that used for draw the plot
map<int, TH1F* > RawDecoder::DrawCorrectedRawHisto(map<int,vector<int> > mMapping,TCanvas *c) {
    static int idx_ec=0;
    static int mpd_off=9999; // first mpd id
    int mpd_id=0;
    int adc_ch=0;
    int hybridID=0;
    int nbAPVs = 0 ;

    map<int, vector<int> >::iterator it;

    if (mpd_off == 9999) { // compute for first event only
        for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it) {
            hybridID=it->first;
            mpd_id = GetMPD_ID(hybridID);
            if (mpd_off>mpd_id) mpd_off=mpd_id;
        }
    }

    std::map<int, std::map<int, TH1F *>>planeRawHist;

    std::cout<<"Map Size"<<mMapping.size()<<std::endl;
    //check how many GEM modules and CREATE
    std::map<int, std::map<int, int>> maxStrips;
    for (auto iter = mMapping.begin(); iter !=mMapping.end(); iter++){
        // load the pdestal to read
        hybridID = iter->first;
        int detID = mMapping[hybridID][0];
        int planeID = mMapping[hybridID][1];
        int RstripNb=128;
        int RstripPos=RstripNb+128*mMapping[hybridID][2];

        if((!(maxStrips.find(detID)!=maxStrips.end() || maxStrips[detID].find(planeID)!=maxStrips[detID].end()))|| maxStrips[detID][planeID] < RstripPos){
            maxStrips[detID][planeID]=RstripPos;
        }
    }

    // print out the detector information
    for (auto iter = maxStrips.begin(); iter!=maxStrips.end(); iter++){
        for (auto itter = iter->second.begin(); itter!=iter->second.end() ; ++itter) {
            std::cout <<"Detector::"<<iter->first<<"    PlaneID::"<<itter->first<<"    value::"<<itter->second<<std::endl;
            planeRawHist[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_rawHist",iter->first,itter->first),Form("Detector%d_dimension%d_rawHist",iter->first,itter->first),itter->second,0,itter->second);
        }
    }

    // write the data in to the orded raw histo
    for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
    {
        hybridID = it->first;
        mpd_id = GetMPD_ID(hybridID);
        adc_ch = GetADC_ch(hybridID);

        vector<int> adc_temp = it->second;
        int N = adc_temp.size();//774
        int TSsize=N/129;

        // check the macth of the frame size
        if(adc_temp.size()!=129*3){
            if((adc_temp.size()>129*3)&&(adc_temp[386]==2)&&(adc_temp[128]==0)&&(adc_temp[257]==1)){
                adc_temp.resize(129*3);
            }else{
                std::cout<<"  MisMatch::  expect 387   get : "<<adc_temp.size()<<"   skip!!("<<__FUNCTION__<<"@"<<__LINE__<<")"<<std::endl;
                return  mAPVRawHisto;
            }
        }
        for(int j=0; j<128;j++){
            int adcSum_temp=0;
            for(int i=0;i<TSsize;i++){
                adcSum_temp = adcSum_temp+adc_temp[j+129*i];
            }
            adcSum_temp = adcSum_temp/TSsize;

            //get the position
            int RstripPos=j;
            int RstripNb = ChNb[j];

            // check the existance of HIOD
            if(mMapping.find(hybridID)!=mMapping.end()){
                RstripNb=RstripNb+(127-2*RstripNb)*mMapping[hybridID][3];                   //re-matching for inverted strips Nb
                RstripPos=RstripNb+128*mMapping[hybridID][2];                               // calculate position
                int detID = mMapping[hybridID][0];
                int planeID = mMapping[hybridID][1];
                planeRawHist[detID][planeID]->Fill(RstripPos, adcSum_temp);
            } else{
                std::cout<<"[WORNING]::"<<__FUNCTION__ <<" CAN NOT FIND "<<hybridID<<" in the Mapping file!!!!!!!(MPD->"<<mpd_id<<",  apv->"<<adc_ch<<")"<<std::endl;
            }

        }
    }


    c->Clear();
    c->Divide(2,planeRawHist.size());
    double counter_temp=1;
    for (auto iter = planeRawHist.begin();iter!=planeRawHist.end(); iter++){
        c->cd(counter_temp);
        if(iter->second.find(0)!=iter->second.end()) (iter->second)[0]->Draw("hist");
        c->cd(counter_temp+1);
        if(iter->second.find(1)!=iter->second.end()) (iter->second)[1]->Draw("hist");
        counter_temp+=2;
    }
    c->Update();

    int typc = getchar();
    // int typc=-1;
    //  c->SaveAs(Form("Result/rawhisto_e%04d.png",idx_ec));
    idx_ec++;

    if ((typc==83) || (typc==115)) { // if s or S has pressed, it save the histos in a png file
        c->SaveAs("Result/rawhisto.png");
        getchar(); // return key
    }

    if ((typc==80) || (typc==112)) { // if p or P has pressed, it save the histos in a pdf file
        c->SaveAs("Result/rawhisto.pdf");
        getchar(); // return key
    }

    for(auto iter = planeRawHist.begin(); iter!=planeRawHist.end();iter++){
        for(auto itter=iter->second.begin(); itter!=iter->second.end(); itter++){
            itter->second->Delete();
        }
    }
    return mAPVRawHisto;
}


//===========================================================================
map<int, TH1F* > RawDecoder::DrawRawHisto(TCanvas *c)
{
  static int idx_ec=0;
  static int mpd_off=9999; // first mpd id
  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;
  int nbAPVs = 0 ;

  map<int, vector<int> >::iterator it;

  if (mpd_off == 9999) { // compute for first event only
    for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it) {
      hybridID=it->first;
      mpd_id = GetMPD_ID(hybridID);     
      if (mpd_off>mpd_id) mpd_off=mpd_id;
    }
  }

  int mpd_count = -1, last_mpd_id = -1, draw_index;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
    {
      hybridID=it->first;
      mpd_id = GetMPD_ID(hybridID);     
      adc_ch = GetADC_ch(hybridID);
      vector<int> adc_temp = it->second;
      int N = adc_temp.size();//cout<<"adc_tempsize:"<<N<<endl;
 
     // std::cout<<"MPD:: "<< mpd_id<<std::endl;
      TH1F* h = new TH1F(Form("mpd_%d_ch_%d",mpd_id, adc_ch), Form("mpd_%d_ch_%d_raw_data",mpd_id, adc_ch), 780, 0, 779);
 
      //cout<<"EventNb = "<<idx_ec<<"  mpdid: "<<mpd_id<<"  adcCh: "<<adc_ch<<"  histo: "<<h->GetName()<<" nbAPVs: "<<nbAPVs<<endl;
 
      for(int i=0;i<N;i++) h->Fill(i+1, (Float_t) adc_temp[i]);
       
      mAPVRawHisto[hybridID] = h;

      if(mpd_id != last_mpd_id){
	mpd_count++;
      }

      draw_index = mpd_count * 16 + adc_ch + 1;

      //c->cd((mpd_id-mpd_off)*15+adc_ch+1)->SetLogy();
      c->cd(draw_index);//->SetLogy();
      mAPVRawHisto[hybridID]->SetMaximum(3000);
      mAPVRawHisto[hybridID]->SetMinimum(100);
      mAPVRawHisto[hybridID]->Draw("HISTO");
      nbAPVs++ ;
      last_mpd_id = mpd_id;
      
    }
  map<int, TH1F*>::iterator itRaw;

  /*
    int i;
  for(itRaw = mAPVRawHisto.begin();itRaw!=mAPVRawHisto.end();itRaw++)
    {
      c->cd(i);
      itRaw->second->SetTitleSize(0.2);
      itRaw->second->SetMaximum(3500);
      itRaw->second->SetMinimum(300);
      itRaw->second->Draw();
      i++;
    }
  */
  c->Update();
 
  int typc = getchar();
  // int typc=-1;
  //  c->SaveAs(Form("Result/rawhisto_e%04d.png",idx_ec));
  idx_ec++;

  if ((typc==83) || (typc==115)) { // if s or S has pressed, it save the histos in a png file
    c->SaveAs("Result/rawhisto.root");
    getchar(); // return key
  }

  if ((typc==80) || (typc==112)) { // if p or P has pressed, it save the histos in a pdf file
    c->SaveAs("Result/rawhisto.pdf");
    getchar(); // return key
  }

  return mAPVRawHisto;
}

map<int, vector<int> > RawDecoder::GetStripTsAdcMap()
{
  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;

  map<int, vector<int> > nullmap;
  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    
    hybridID=it->first;mpd_id= GetMPD_ID(hybridID);
    adc_ch = GetADC_ch(hybridID);
    vector<int> adc_temp = it->second;

    //std::cout<<"[Test]:: MPD->"<<GetMPD_ID(hybridID)<<"  ADCid->"<<GetADC_ch(hybridID);//<<std::endl; //774 for 6 time samples 387 for 3 time samples

    int TSsize=adc_temp.size()/129;

    if(adc_temp.size()!=129*3){
       
       if((adc_temp.size()>129*3)&&(adc_temp[386]==2)&&(adc_temp[128]==0)&&(adc_temp[257]==1)){
         adc_temp.resize(129*3);
       }else{
         std::cout<<"  MisMatch::  expect 387   get : "<<adc_temp.size()<<"   skip!!"<<std::endl;
    	
	 //for (int i =0; i < adc_temp.size(); i++){
	//	 std::cout<< i << "("<< "MPD "<<mpd_id <<", adcID "<<adc_ch<<", " <<adc_temp[i]<<")";
	 //}
	 //std::cout<<std::endl;
   	 //getchar(); 

	 return nullmap;
       }	
   }
// check the dataframe miss match 
   
    if(!((adc_temp[386]==2)&&(adc_temp[128]==0)&&(adc_temp[257]==1))){
    std::cout<<"  MisMatch::  the data strcuture is wrong!!! "<<__FUNCTION__<<" @"<<__LINE__<<std::endl;
    }

    for(int i=0; i<TSsize;i++)
      {
	      vector<int> singleTSadc_temp;
	      singleTSadc_temp.insert(singleTSadc_temp.end(),&adc_temp[129*i],&adc_temp[129*(i+1)]);
	      //cout<<"singleTSADCSIZE: "<<singleTSadc_temp.size()<<"     ";
	      vector<int> singleTSadc_temp_sorting; 
	      //	for(int j=0; j<128;j++)
	      //	  {
	      //	if(mpd_id==8&&adc_ch==0){cout<<dec<<"i:"<<i<<"singleTSadc_temp:"<<singleTSadc_temp[j]<<endl;}
	      //	  }

	      singleTSadc_temp_sorting.insert(singleTSadc_temp_sorting.end(),singleTSadc_temp.begin(),singleTSadc_temp.end());
	      sort(singleTSadc_temp_sorting.begin(),singleTSadc_temp_sorting.end()-1);
	      int iCommonMode=0;

	      for ( int k=28; k <100; k++) //omitting largest 4 channels for common mode, necessary to calculate this from the middle
	        {
	          iCommonMode+=singleTSadc_temp_sorting[k];
	          //if((i==0|i==1)&&adc_ch==0)cout <<i<<"  "<<iCommonMode << " ";//if(k==singleTSadc_temp.size()-1){cout<<endl;}
	        }

	      iCommonMode = iCommonMode/72; //cout<<"i "<<i<<"commonmode: "<<iCommonMode<<endl;
	      for ( int k=0; k <singleTSadc_temp.size()-1; k++) 
	        {
	            singleTSadc_temp[k]-=iCommonMode;
	            //cout<<"commonmode: "<<singleTSadc_temp[k]<<endl;
	        }

	      int temphybridID;
	      for(int j=0; j<128;j++)
	        {
	          //cout<<"hybridid:"<<hex<<hybridID;
	          temphybridID=hybridID|j;     //if(adc_ch==0){cout<<hex<<"temphybridid:"<<temphybridID<<endl;}
	          mPedestalTsAdc[temphybridID].push_back(singleTSadc_temp[j]);//singleTSadc_temp[j];
	          //cout <<singleTSadc_temp[j]<<endl;
	          //   if(mpd_id==8&&adc_ch==0){cout<<dec<<"i:"<<i<<"singleTSadc_temp:"<<singleTSadc_temp[j]<<endl;}
	          //if((i==0|i==1)&&adc_ch==0)cout <<i<<" XXX  "<<singleTSadc_temp[j] << " ";
	        }
          
      }
  }// 
  return mPedestalTsAdc;
}

// used for get the hit mode result
map<int, vector<int> > RawDecoder::ZeroSup(map<int,vector<int> > mMapping, map<int,vector<int> > mPedestalMean, map<int,vector<int> > mPedestalRMS)
{

  map<int, vector<int> > mmHit;
  map<int, vector<int> > nullmap;
  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;




    //===================================================================================================================
    std::map<int, std::map<int, TH1F *>>PlaneRawHistRaw;  // strips　
    std::map<int, std::map<int, TH1F *>>PlaneCutThresHist;  // strips　
    std::map<int, std::map<int, TH1F *>>PlaneRawHist;  // strips that pass the 5-sigma cut
    std::map<int, std::map<int, TH1F *>>PlaneRawHist_NonZeroCut;  // strips that pass the 5-sigma cut
//    std::cout<<"Map Size"<<mMapping.size()<<std::endl;
    //check how many GEM modules and CREATE
    std::map<int, std::map<int, int>> maxStrips;
    for (auto iter = mMapping.begin(); iter !=mMapping.end(); iter++){
        // load the pdestal to read
        int hybridID = iter->first;
        int detID = mMapping[hybridID][0];
        int planeID = mMapping[hybridID][1];
        int RstripNb=128;
        int RstripPos=RstripNb+128*mMapping[hybridID][2];

        if((!(maxStrips.find(detID)!=maxStrips.end() || maxStrips[detID].find(planeID)!=maxStrips[detID].end()))|| maxStrips[detID][planeID] < RstripPos){
            maxStrips[detID][planeID]=RstripPos;
        }
    }

    // print out the detector information
    for (auto iter = maxStrips.begin(); iter!=maxStrips.end(); iter++){
        for (auto itter = iter->second.begin(); itter!=iter->second.end() ; ++itter) {
            PlaneRawHistRaw[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_rawHist",iter->first,itter->first),Form("Detector%d_dimension%d_rawHist",iter->first,itter->first),itter->second,0,itter->second);
            PlaneRawHist[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_HitHist",iter->first,itter->first),Form("Detector%d_dimension%d_HitHist",iter->first,itter->first),itter->second,0,itter->second);
            PlaneRawHist_NonZeroCut[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_before0",iter->first,itter->first),Form("Detector%d_dimension%d_before0",iter->first,itter->first),itter->second,0,itter->second);
            PlaneCutThresHist[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_thrHist",iter->first,itter->first),Form("Detector%d_dimension%d_thrHist",iter->first,itter->first),itter->second,0,itter->second);
            PlaneCutThresHist[iter->first][itter->first]->SetLineColor(6);
            PlaneRawHistRaw[iter->first][itter->first]->SetLineColor(7);
            PlaneRawHist[iter->first][itter->first]->SetLineColor(2);
        }
    }
//===================================================================================================================


  double gemPlaneMaxRMS=0.0;
  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    hybridID = it->first;
    mpd_id = GetMPD_ID(hybridID);     
    adc_ch = GetADC_ch(hybridID);
    
    vector<int> adc_temp = it->second;
    int N = adc_temp.size();//774
    int TSsize=N/129;

    //if ((N!=388)&&(N!=387)&&(N!=386))
    //{
    //	std::cout<<"[ERROR]:: expected apv frame size is :"<<N<<"  ("<<TSsize*129+1<<")"<<std::endl;
    //	continue;
    //}
    
    // check the macth of the frame size
    if(adc_temp.size()!=129*3){
   	if((adc_temp.size()>129*3)&&(adc_temp[386]==2)&&(adc_temp[128]==0)&&(adc_temp[257]==1)){
		adc_temp.resize(129*3);
	}else{
		std::cout<<"  MisMatch::  expect 387   get : "<<adc_temp.size()<<"   skip!!("<<__FUNCTION__<<"@"<<__LINE__<<")"<<std::endl;
		return nullmap;
	}
    }

    //cout<<TSsize<<endl;//774
    int CommonMode[TSsize];


    for(int i=0;i<N;i++)
	{
	  if((i%129)!=128) adc_temp[i]-=mPedestalMean[hybridID][i%129];
	}


    for(int i=0; i<TSsize;i++)
      { 
	CommonMode[i]=0;
	vector<int> singleTSadc_temp_sorting;
	singleTSadc_temp_sorting.insert(singleTSadc_temp_sorting.end(),&adc_temp[129*i],&adc_temp[129*(i+1)]);
	  	 
	//cout<<(singleTSadc_temp_sorting.size()-5)<<endl;
	sort(singleTSadc_temp_sorting.begin(),singleTSadc_temp_sorting.end()-1);

	for ( int k=28; k <100; k++) //omitting largest 4 channels for common mode, necessary to calculate this from the middle
	  {  
	    CommonMode[i]+=singleTSadc_temp_sorting[k];
	    //CommonMode[i]+=500;

	  }

	CommonMode[i] = CommonMode[i]/72; 
	//	cout<<"Commonmode: "<<CommonMode[i];
      }
 
    for(int j=0; j<128;j++)
      {
	int adcSum_temp=0;
	for(int i=0;i<TSsize;i++)
	  {//cout<<"ADC:"<<adc_temp[j+129*i]<<"    ";
	    adcSum_temp = adcSum_temp+adc_temp[j+129*i]-CommonMode[i];
	    //cout<<"ADC:"<<adcSum_temp<<"    ";
	  }
	// cout<<endl<<endl;
	adcSum_temp = adcSum_temp/TSsize; 
	 
	//cout<<j<<"Mean  "<<mPedestalMean[mpd_id][adc_ch][j]<<"  "<<mPedestalRMS[mpd_id][adc_ch][j]<<endl;

	//========================================================================================================
	{
	    if((mPedestalRMS.find(hybridID)!=mPedestalRMS.end()))
              {
                  int RstripPos=j;
                  int RstripNb = ChNb[j];

                  RstripNb=RstripNb+(127-2*RstripNb)*mMapping[hybridID][3];                   //re-matching for inverted strips Nb
                  RstripPos=RstripNb+128*mMapping[hybridID][2];                               // calculate position

                  int detID = mMapping[hybridID][0];
                  int planeID = mMapping[hybridID][1];
                  int HitHybridID = (detID<<13)|(planeID<<12)|RstripPos;

                  PlaneRawHist_NonZeroCut[detID][planeID]->Fill(RstripPos,adcSum_temp);
                  if((adcSum_temp>sigmacut*mPedestalRMS[hybridID][j])){
                      PlaneRawHist[detID][planeID]->Fill(RstripPos,adcSum_temp);
                  }
                  // Get all the hit, and draw it into the files
              }

    }
//========================================================================================================



	if(adcSum_temp>sigmacut*mPedestalRMS[hybridID][j])
	  { 
	  
	    int RstripPos=j;	  
	    int RstripNb = ChNb[j];
	    //int RstripNb=32*(j%4)+8*(int)(j/4)-31*(int)(j/16);                        //channel re-matching for apv25 chip
	    ////stripNb=(8*(int)(stripNb/4)+3-stripNb)*((int)(stripNb/4)%2)+stripNb*(1-((int)(stripNb/4)%2));
	    //RstripNb=RstripNb+1+RstripNb%4-5*(((int)(RstripNb/4))%2);                 //channel re-matching for INFN type APV front-end card
	    //cout<<RstripNb<<", ";
	    RstripNb=RstripNb+(127-2*RstripNb)*mMapping[hybridID][3];                   //re-matching for inverted strips Nb
	    RstripPos=RstripNb+128*mMapping[hybridID][2];                               // calculate position

	    int detID = mMapping[hybridID][0];
	    int planeID = mMapping[hybridID][1];
	    int HitHybridID = (detID<<13)|(planeID<<12)|RstripPos;
	    for(int i=0; i<TSsize;i++)
	      { 
		mmHit[HitHybridID].push_back(adc_temp[j+129*i]-CommonMode[i]);//if stop here, reduce 200ms/10k event
		//	cout<<(adc_temp[j+129*i]-CommonMode[i]-mPedestalMean[hybridID][j])<<"  ";
	      }     
	    //	    cout<<endl;
	    //  j++;
	  } 

      }  
  }//
  
  // get the maximum RMS, AND apply cut
  for(auto iter = PlaneRawHist.begin(); iter!=PlaneRawHist.end(); iter++){
      for (auto itter=iter->second.begin();itter!=iter->second.end();itter++){
          if(gemPlaneMaxRMS< itter->second->GetRMS()) gemPlaneMaxRMS=itter->second->GetRMS();
      }
  }

  if(gemPlaneMaxRMS > 280) mmHit.clear();

  return mmHit;
}



map<int, vector<int> > RawDecoder::DrawHits(map<int,vector<int> > mMapping, map<int,vector<int> > mPedestalMean, map<int,vector<int> > mPedestalRMS, TCanvas *c)
{
  std::map<int,std::map<int,std::map<int,int>>> plotData;

    std::map<int, std::map<int, TH1F *>>PlaneRawHistRaw;  // strips　
    std::map<int, std::map<int, TH1F *>>PlaneCutThresHist;  // strips　
    std::map<int, std::map<int, TH1F *>>PlaneRawHist;  // strips that pass the 5-sigma cut
    std::map<int, std::map<int, TH1F *>>PlaneRawHist_NonZeroCut;  // strips that pass the 5-sigma cut


//    std::cout<<"Map Size"<<mMapping.size()<<std::endl;
    //check how many GEM modules and CREATE
    std::map<int, std::map<int, int>> maxStrips;
    for (auto iter = mMapping.begin(); iter !=mMapping.end(); iter++){
        // load the pdestal to read
        int hybridID = iter->first;
        int detID = mMapping[hybridID][0];
        int planeID = mMapping[hybridID][1];
        int RstripNb=128;
        int RstripPos=RstripNb+128*mMapping[hybridID][2];

        if((!(maxStrips.find(detID)!=maxStrips.end() || maxStrips[detID].find(planeID)!=maxStrips[detID].end()))|| maxStrips[detID][planeID] < RstripPos){
            maxStrips[detID][planeID]=RstripPos;
        }
    }

    // print out the detector information
    for (auto iter = maxStrips.begin(); iter!=maxStrips.end(); iter++){
        for (auto itter = iter->second.begin(); itter!=iter->second.end() ; ++itter) {
//            std::cout <<"Detector::"<<iter->first<<"    PlaneID::"<<itter->first<<"    value::"<<itter->second<<std::endl;
            PlaneRawHistRaw[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_rawHist",iter->first,itter->first),Form("Detector%d_dimension%d_rawHist",iter->first,itter->first),itter->second,0,itter->second);
            PlaneRawHist[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_HitHist",iter->first,itter->first),Form("Detector%d_dimension%d_HitHist",iter->first,itter->first),itter->second,0,itter->second);
            PlaneRawHist_NonZeroCut[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_before0",iter->first,itter->first),Form("Detector%d_dimension%d_before0",iter->first,itter->first),itter->second,0,itter->second);
//            PlaneRawHist_NonZeroCut[iter->first][itter->first]->SetLineColor()
            PlaneCutThresHist[iter->first][itter->first]=new TH1F(Form("Detector%d_dimension%d_thrHist",iter->first,itter->first),Form("Detector%d_dimension%d_thrHist",iter->first,itter->first),itter->second,0,itter->second);
            PlaneCutThresHist[iter->first][itter->first]->SetLineColor(6);
            PlaneRawHistRaw[iter->first][itter->first]->SetLineColor(7);
            PlaneRawHist[iter->first][itter->first]->SetLineColor(2);
        }
    }


  map<int, vector<int> > mmHit;
  int mpd_id=0;
  int adc_ch=0;
  int hybridID=0;

  map<int, vector<int> >::iterator it;
  for(it = mAPVRawSingleEvent.begin(); it!=mAPVRawSingleEvent.end(); ++it)
  {
    hybridID = it->first;
    mpd_id = GetMPD_ID(hybridID);     
    adc_ch = GetADC_ch(hybridID);
    
    vector<int> adc_temp = it->second;
    int N = adc_temp.size();//774
    int TSsize=N/129;

    // check the size of the frame
    if(adc_temp.size()!=129*3){
          if((adc_temp.size()>129*3)&&(adc_temp[386]==2)&&(adc_temp[128]==0)&&(adc_temp[257]==1)){
              adc_temp.resize(129*3);
          }else{
              std::cout<<"  MisMatch::  expect 387   get : "<<adc_temp.size()<<"   skip!!("<<__FUNCTION__<<"@"<<__LINE__<<")"<<std::endl;
              continue;
          }
      }
    //cout<<TSsize<<endl;//774
    int CommonMode[TSsize];

    // load the pedestals
    for(int i=0;i<N;i++)
	  {
          if(mPedestalMean.find(hybridID)!=mPedestalMean.end())
          {
            if((i%129)!=128) adc_temp[i]-=mPedestalMean[hybridID][i%129];
            //std::cout<<"Loading MPD="<<GetMPD_ID(hybridID)<<" ADC="<<GetADC_ch(hybridID)<<" in the pedestal file"<<std::endl;
          }else{
            std::cout<<"ERROR : can not find MPD="<<GetMPD_ID(hybridID)<<" ADC="<<GetADC_ch(hybridID)<<" in the pedestal file"<<std::endl;
          }
	  }

    for(int i=0; i<TSsize;i++)
      { 
	      CommonMode[i]=0;
	      vector<int> singleTSadc_temp_sorting;
	      singleTSadc_temp_sorting.insert(singleTSadc_temp_sorting.end(),&adc_temp[129*i],&adc_temp[129*(i+1)]);
	  	 
	      //cout<<(singleTSadc_temp_sorting.size()-5)<<endl;
	      sort(singleTSadc_temp_sorting.begin(),singleTSadc_temp_sorting.end()-1);

	      for ( int k=28; k <100; k++) //omitting largest 4 channels for common mode, necessary to calculate this from the middle
	      {  
	        CommonMode[i]+=singleTSadc_temp_sorting[k];
	        //CommonMode[i]+=500;
	      }
	      CommonMode[i] = CommonMode[i]/72;
	      //	cout<<"Commonmode: "<<CommonMode[i];
      }
 
    for(int j=0; j<128;j++)
      {

        int rawADC_aver=0;
          for(int i=0;i<TSsize;i++){
              rawADC_aver = rawADC_aver+adc_temp[j+129*i];
          }
          rawADC_aver = rawADC_aver/TSsize;

	      int adcSum_temp=0;
	      for(int i=0;i<TSsize;i++)
	        {//cout<<"ADC:"<<adc_temp[j+129*i]<<"    ";
	          adcSum_temp = adcSum_temp+adc_temp[j+129*i]-CommonMode[i];
	          //cout<<"ADC:"<<adcSum_temp<<"    ";
	        }
	      // cout<<endl<<endl;
	      adcSum_temp = adcSum_temp/TSsize; 
  	    //cout<<j<<"Mean  "<<mPedestalMean[mpd_id][adc_ch][j]<<"  "<<mPedestalRMS[mpd_id][adc_ch][j]<<endl;

	      if((mPedestalRMS.find(hybridID)!=mPedestalRMS.end()))
	        {
//            std::cout<<"Working on  MPD="<<GetMPD_ID(hybridID)<<" ADC="<<GetADC_ch(hybridID)<<" , get the graph  "<<__LINE__<<std::endl;
            
	          int RstripPos=j;	  
	          int RstripNb = ChNb[j];
	          
	          RstripNb=RstripNb+(127-2*RstripNb)*mMapping[hybridID][3];                   //re-matching for inverted strips Nb
	          RstripPos=RstripNb+128*mMapping[hybridID][2];                               // calculate position

	          int detID = mMapping[hybridID][0];
	          int planeID = mMapping[hybridID][1];
	          int HitHybridID = (detID<<13)|(planeID<<12)|RstripPos;

	          PlaneRawHistRaw[detID][planeID]->Fill(RstripPos,rawADC_aver);
	          PlaneRawHist_NonZeroCut[detID][planeID]->Fill(RstripPos,adcSum_temp);
	          PlaneCutThresHist[detID][planeID]->Fill(RstripPos,sigmacut*mPedestalRMS[hybridID][j]);
              if((adcSum_temp>sigmacut*mPedestalRMS[hybridID][j])){
                  PlaneRawHist[detID][planeID]->Fill(RstripPos,adcSum_temp);
	          }
	          plotData[detID][planeID][RstripNb]=adcSum_temp;
            // Get all the hit, and draw it into the files
	      } 
      }  
  }// 

  // draw the canvas
    c->Clear();
    c->Divide(2,PlaneRawHist.size());
    double counter_temp=1;
    for (auto iter = PlaneRawHist.begin();iter!=PlaneRawHist.end(); iter++){
        c->cd(counter_temp);
//        PlaneRawHistRaw[iter->first][0]->Draw("hist");
        PlaneRawHist_NonZeroCut[iter->first][0]->Draw("hist same");
//        PlaneCutThresHist[iter->first][0]->Draw("hist same");
        if(iter->second.find(0)!=iter->second.end()) {
            (iter->second)[0]->Draw("hist same");
            TPaveText *text0=new TPaveText(0.1,0.6,0.4,0.9,"NDC");
            text0->AddText(Form("RMS:%f",(iter->second)[0]->GetRMS()));
            text0->SetTextColor(kRed);
            text0->Draw("same");
        }



        c->cd(counter_temp+1);
//        PlaneRawHistRaw[iter->first][1]->Draw("hist");
        PlaneRawHist_NonZeroCut[iter->first][1]->Draw("hist same");
//        PlaneCutThresHist[iter->first][1]->Draw("hist same");
        if(iter->second.find(1)!=iter->second.end()) {
            (iter->second)[1]->Draw("hist same");
            TPaveText *text1=new TPaveText(0.1,0.6,0.4,0.9,"NDC");
            text1->AddText(Form("RMS:%f",(iter->second)[1]->GetRMS()));
            text1->SetTextColor(kRed);
            text1->Draw("same");
        }
        counter_temp+=2;
    }
    c->Update();

    int typc = getchar();
    // int typc=-1;
    //  c->SaveAs(Form("Result/rawhisto_e%04d.png",idx_ec));

    if ((typc==83) || (typc==115)) { // if s or S has pressed, it save the histos in a png file
        c->SaveAs("Result/histoplot.root");
        getchar(); // return key
    }

    if ((typc==80) || (typc==112)) { // if p or P has pressed, it save the histos in a pdf file
        c->SaveAs("Result/rawhisto.pdf");
        getchar(); // return key
    }

    for(auto iter = PlaneRawHist.begin(); iter!=PlaneRawHist.end();iter++){
        for(auto itter=iter->second.begin(); itter!=iter->second.end(); itter++){
            itter->second->Delete();
        }
    }

  return mmHit;
}
