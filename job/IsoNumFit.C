#include  <iostream>
#include  <dirent.h>//DIR
#include  <TCanvas.h>
#include  <TMath.h>
#include  <TLegend.h>
#include  <TH1F.h>
#include  <TFile.h>
#include  <fstream>
#include  <sstream>
#include  <TStyle.h>
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include    "RooPlot.h"
#include "RooFitResult.h"
#include  <string>
#include  <TChain.h>
#include  <TROOT.h>
#include  <vector>
using namespace RooFit;

int fitHisto(int siteNum,string dataVer)
{

//===>initialize variable
	bool anaLiHe=true;//He8/Li9 bg
    int ADNumOfSite[3]={0};
    int daqHistNum=5;
    if( dataVer.find("13")!=string::npos || dataVer.find("14")!=string::npos || dataVer.find("15")!=string::npos || dataVer.find("16")!=string::npos)
    {
        ADNumOfSite[0]=2;
        ADNumOfSite[1]=2;
        ADNumOfSite[2]=4;
        daqHistNum=5;
    } else
    {
        if( dataVer.find("11")!=string::npos || dataVer.find("12")!=string::npos )
        {
            ADNumOfSite[0]=2;
            ADNumOfSite[1]=1;
            ADNumOfSite[2]=3;
            daqHistNum=4;
        }
    }
	string site;
	if( siteNum==1 )
	{
		site="EH1";
	} else if(siteNum==2)
	{
		site="EH2";
	}else if (siteNum==3)
	{
		site="EH3";
	}else
	{
		site="EH1";
	}

//===>initialize histograms
	//livetime
	string runnum;
	TH1F* h[daqHistNum];//the fifth is for daqtime
	Double_t totalTime[5]={0.};

	//Li9 He8
	double n98total=0.;
	//double in98total=0.;
    double ifitn98total=0.;
    double n98Rate=0.;
    double in98Rate=0.;
    double tlivetime=0.;
	TH1F* showermuonNum[6]; 
	TH1F* time2lastshowermuon[6];
    TH1F* LiResult[5];
    TH1F* hh[6];

//===>get histograms from .root file for analyse
	string filename;
    filename=site;
    filename+="totalHisto_";
    filename+=dataVer;
    filename+=".root";
    TFile* f=new TFile(filename.c_str());

		//livetime
		for( int i=0 ; i<daqHistNum ; i++ )
		{

			stringstream hname;
			if( i==daqHistNum-1 )
			{
				hname<<"LiveTime/DaqTime";
			} else
			{
				hname<<"LiveTime/AD";
				hname<<i+1;
				hname<<"LiveTime";	
			}
			h[i] = (TH1F*)f->Get(hname.str().c_str());
			if( !(h[i]) )
			{
				cout<<"Can not get Hist : "<<hname.str()<<" from "<<site<<"/"<<runnum<<" ."<<endl;
				//return true;
                continue;
			}
		}

		//Li9 He8
        if(anaLiHe)
		{
            for( int i=0 ; i<6 ; i++ )
            {
                TString hnameLi;
                hnameLi="lidj/showermuonNum";
                hnameLi+=i+1;
                showermuonNum[i]=(TH1F*)f->Get(hnameLi);
                hnameLi="lidj/time2lastshowermuon";
                hnameLi+=i+1;
                //hnameLi+="4Li";
                time2lastshowermuon[i]=(TH1F*)f->Get(hnameLi);
            }
		}
	
//===>analyse 
	TAxis *xaxis = h[0]->GetXaxis();
	int binNum = xaxis->GetNbins();
	//livetime
	for( int i=0 ; i<daqHistNum ; i++ )
	{
		totalTime[i]=0.;
		for( int j=1 ; j<=binNum ; j++ )
		{
			totalTime[i]+=h[i]->GetBinContent(j);
		}
		
	}

	//Li He
    std::cout<<"begin to analyse Li "<<endl;
	if(anaLiHe)
	{
        TString xTitle[6]={"20~60MeV","60~500MeV","0.5~1.5GeV","1.5~2.5GeV",">2.5GeV",">2.5GeV"};
        double showerTh[6] = {0.02, 0.06, 0.5, 1.5, 2.5, 5.0};
        LiResult[0]= new TH1F("B12Yield", "All", 5, showerTh);
        LiResult[1]= new TH1F("B12YieldAD1", "AD1", 5, showerTh);
        LiResult[2]= new TH1F("B12YieldAD2", "AD2", 5, showerTh);
        LiResult[3]= new TH1F("B12YieldAD3", "AD3", 5, showerTh);
        LiResult[4]= new TH1F("B12YieldAD4", "AD4", 5, showerTh);
		double NumMuon[6]={0.};
		double RateMuon[6]={0.};
		double NumIbd[6]={0.};
        
		RooRealVar x("x","x",0.001,0.5, "s");
		RooRealVar tauLi9("tauLi9", "tauLi9", -0.0291, "s");
		RooRealVar rateMu("rateMu","rate of showermuon",-0.1,-10., 0.,"Hz");
		RooRealVar N98("N98","total number of Li9 and He8",500.,0.,1.e7);
		RooFormulaVar lambdaLi9("lambdaLi9","lambdaLi9","1/@0 + @1",RooArgList(tauLi9, rateMu));
		RooRealVar NIbd("NIbd","number of background",3e2, 1e1, 2.e7);
		RooExponential expLi9("expLi9", "Li9 distribution", x, lambdaLi9);
		RooExponential expIbd("expIbd","Ibd distribution",x,rateMu);
		
		RooFitResult* fitres;
		RooAddPdf* sum;
		RooDataHist* data;
        //RooPlot* mesframe[6];
		double nIbd[6]={0.};
		double inIbd[6]={0.};
        double rateMuValue[6]={0.};
        double rateMuErr[6]={0.};
		double n98[6]={0.};
		double in98[6]={0.};
		double binwidth=0.;
        TH1F* lh[5][6];
        TString lname[5]={"","AD1","AD2","AD3","AD4"};
        int lcolor[5]={4,3,2,6,5};
        TCanvas* c = new TCanvas("c","c",2000,800) ; 
        c->Divide(3,2);
        gStyle->SetEndErrorSize(5.0);
        gStyle->SetMarkerSize(0.1);
        //gStyle->SetHistLineWidth(1);
                int hnum=time2lastshowermuon[0]->FindBin(100);
        for( int k=0 ; k<6 ; k++ )
        {
            TString hhname="slice";
            hhname+=k+1;
            hh[k]=new TH1F(hhname,"hh",hnum-2,0.001,100);
                for( int i=1 ; i<hnum-2 ; i++ )
                {
                    hh[k]->SetBinContent(i,time2lastshowermuon[k]->GetBinContent(i+1));
                }
                hh[k]->SetOption("E1");
        }
            if( 0 )
            {
                // it's strange ,if you can't draw following six figs in the big 3*2 pad,you can fist use these code,run once,don't exit the ROOT,then you can draw six figs!!!
	            hh[0]->Rebin(5);
				data = new RooDataHist("data", "data", x, hh[0]);
				sum  = new RooAddPdf("sum","sum",RooArgList(expLi9, expIbd),RooArgList(N98, NIbd));
	            RooPlot* mesframe = x.frame() ;
	            data->plotOn(mesframe) ;
	            c->cd(1);
	            mesframe->Draw();
                return true;
                
            }
        for( int ihist=0 ; ihist<1 ; ihist++ )
        {
            
			for( int j=0 ; j<5 ; j++ )
			{
	            std::cout<<"now is  "<< j+1<<endl;
	            NumIbd[j]=hh[j]->Integral(1,hh[j]->FindBin(0.5));
	            hh[j]->Rebin(5);
				data = new RooDataHist("data", "data", x, hh[j]);
				sum  = new RooAddPdf("sum","sum",RooArgList(expLi9, expIbd),RooArgList(N98, NIbd));
				fitres = sum->fitTo((*data),Save(),PrintLevel(-1),Extended(kTRUE));
	            //fitres->Print();
				n98[j]=N98.getVal(0);
				in98[j]=N98.getError();
	            nIbd[j]=NIbd.getVal(0);
	            inIbd[j]=NIbd.getError();
	            rateMuValue[j]=rateMu.getVal(0);
	            rateMuErr[j]=rateMu.getError();
	
	            RooPlot* mesframe = x.frame() ;
	            data->plotOn(mesframe) ;
	            sum->plotOn(mesframe);
	            sum->plotOn(mesframe,Components(expIbd),LineStyle(kDashed),LineColor(kGreen)) ;
	            sum->plotOn(mesframe,Components(expLi9),LineStyle(kDashed),LineColor(kRed)) ;
	            xTitle[j]+="        time since last muon (s)";
	            mesframe->GetXaxis()->SetTitle(xTitle[j]);
	            mesframe->GetYaxis()->SetTitle("Entries");
	            c->cd(j+1);
	            mesframe->Draw();
                gPad->SetLogy();
			}
	        
			n98total=n98[0]+n98[1]+n98[2]+(n98[3]+n98[4]+n98[5]);
			ifitn98total=sqrt(in98[0]*in98[0]+in98[1]*in98[1]+in98[2]*in98[2]+in98[3]*in98[3]+in98[4]*in98[4]+in98[5]*in98[5]);
	        if( ihist==0 )
	        {
	            for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
	            {
	                tlivetime+=totalTime[i];
	            }
	        }else
	        {
	            tlivetime=totalTime[ihist-1];
	        }
	        n98Rate=n98total/(tlivetime/(24*3600));
	        in98Rate=ifitn98total/(tlivetime/(24*3600));
	        for( int j=0 ; j<6 ; j++ )
	        {
			    for( int jbin=0 ; jbin<binNum ; jbin++ )
	            {
	                NumMuon[j]+=showermuonNum[j]->GetBinContent(jbin);
	            }
	            LiResult[ihist]->SetBinContent(j+1,n98[j]/NumMuon[j]);
	            LiResult[ihist]->SetBinError(j+1,in98[j]/NumMuon[j]);
				RateMuon[j]=NumMuon[j]/totalTime[daqHistNum-1]/ADNumOfSite[siteNum-1];
	            std::cout<<"   n98["<<j<<"] "<<n98[j]<<"   in98["<<j<<"] "<< in98[j]<<"   fitRateMu["<<j<<"] "<<rateMuValue[j]<<"  realRateMu "<<RateMuon[j]<<"   NumTol "<<NumIbd[j]<<endl;
	        }
	        std::cout<<"n98Num   : "<<n98total<<" +- "<<ifitn98total<<" Rate(/day/AD)  : "<<n98Rate<<" +- "<< in98Rate<<endl;
        }
        
        c->cd(6);
        LiResult[0]->SetMarkerStyle(20);
		LiResult[0]->SetMarkerColor(lcolor[0]);
		LiResult[0]->SetLineColor(lcolor[0]);
		LiResult[0]->SetMarkerSize(1.0);
		LiResult[0]->SetMinimum(0);
		LiResult[0]->SetMaximum(LiResult[0]->GetMaximum()*1.5);
		LiResult[0]->GetXaxis()->SetTitle("Muon visible energy (GeV)");
		LiResult[0]->GetYaxis()->SetTitle("^{12}B/^{12}N yield per muon");
		LiResult[0]->SetTitle("");
        LiResult[0]->SetStats(kFALSE);
		LiResult[0]->Draw("EP");
            /*
        TLegend *legend=new TLegend(.6,.65,.79,.89);
        TString legendLabel=dataVer+" "+site+" "+"All";
        legend->AddEntry(LiResult[0],legendLabel,"lp");

        for( int ihist=1 ; ihist<ADNumOfSite[siteNum-1]+1; ihist++ )
        {
            LiResult[ihist]->SetMarkerStyle(20);
		    LiResult[ihist]->SetMarkerColor(lcolor[ihist]);
		    LiResult[ihist]->SetLineColor(lcolor[ihist]);
		    LiResult[ihist]->SetMarkerSize(1.0);
            LiResult[ihist]->SetStats(kFALSE);
		    LiResult[ihist]->Draw("EPsame");
            legend->AddEntry(LiResult[ihist],lname[ihist],"lp");
        }
        legend->SetFillColor(0);
        legend->Draw();
        */
	}	
	
    //f->Close();
//===>print result
	std::cout<<""<<endl;
	std::cout<<site <<"'s infomation : "<<endl;
	std::cout<<""<<endl;
	for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
	{
		std::cout<<"Total AD"<<i+1<<"LiveTime                         : "<<totalTime[i]/(24*3600)<<endl;
	}
	std::cout<<"Total DaqTime : "<<totalTime[daqHistNum-1]/(24*3600)<<" day" <<endl;
	std::cout<<""<<endl;

	if(anaLiHe)
	{
	    std::cout<<" "<<endl;
        /*
	    std::cout<<"Li9/He8 "<<endl;
	    std::cout<<"n98total  : "<<n98total <<" +- "<<in98total <<" Rate:"<<n98total/totalTime[0]*24*3600 <<" +- "<<in98total/totalTime[0]*24*3600 <<endl;
	    std::cout<<"n98total0(without rpc)  : "<<n98total0 <<" +- "<<in98total0 <<" Rate:"<<n98total0/totalTime0[0]*24*3600 <<" +- "<<in98total0/totalTime0[0]*24*3600 <<endl;
        */
	}

//===>write into .root file
	string rootname=site;
	rootname+="FitResult_"+dataVer+".root";
	TFile* file = new TFile(rootname.c_str(),"RECREATE");
    file->cd();
    for( int i=0 ; i<daqHistNum ; i++ )
    {
        h[i]->Write();
    }
    
    if( anaLiHe )
    {
		//for( int i=0 ; i<ADNumOfSite[siteNum-1] ; i++ )
		for( int i=0 ; i<6 ; i++ )
		{
			hh[i]->Write();
		}
		LiResult[0]->Write();
    }
	file->Close();
	f->Close();
	return 0;
}

/*
int anaTree(int siteNum,vector<string> dataVer)
{
	gROOT->ProcessLine(".L LiTree.C+");
	TChain chain("Tree/LiTree");
	string site;
    if( siteNum==1 )
	{   
	    site="EH1";
	} else if(siteNum==2)
	{   
	    site="EH2";
	}else if (siteNum==3)
	{   
	    site="EH3";
	}else
	{   
	    site="EH1";
	} 

	string chainFile;
	
    chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/";
    chainFile+=dataVer;
    chainFile+="/";
    chainFile+=site;
    chainFile+="/*LiAna.root";
    chain.Add(chainFile.c_str());
    if( dataVer=="P12E" )
    {
        chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/";
        chainFile+="P12A";
        chainFile+="/";
        chainFile+=site;
        chainFile+="/*LiAna.root";
        chain.Add(chainFile.c_str());
    }
    string  siteAndDataVer=site+dataVer;
    std::cout<<"site And DataVer  : "<<siteAndDataVer<<endl;
	chain.Process("LiTree",siteAndDataVer.c_str());
	return 0;
}
*/
vector<string> checkdata(string dataVer)
{
    string dataIsGood="1";
    string verNum[6]={"11","12","13","14","15","16"};
    vector<string> dataVerReal;
    for( int i=0 ; i<6 ; i++ )
    {
        int pos=0;
        string TmpVer;
        TmpVer=dataVer.substr(pos);
        while( pos!=-1 )
        {
            pos=TmpVer.find(verNum[i]);
            if( pos!=-1 )
            {
                string subdataVer;
                subdataVer="P"+verNum[i];
                subdataVer+=toupper(*(TmpVer.substr(pos+2,1).c_str()));
                for( vector<string>::iterator it=dataVerReal.begin() ; it!=dataVerReal.end() ; ++it  )
                {
                    if( *it==subdataVer )
                    {
                        dataVerReal.erase(it);
                        it--;
                    }
                }
                dataVerReal.push_back(subdataVer);
                TmpVer=TmpVer.substr(pos+3);
            }
        }
    }
    sort(dataVerReal.begin(),dataVerReal.end());

    std::cout<<" "<<endl;
    std::cout<<"Find following data version will be analysed : "<<endl;
    std::cout<<" "<<endl;
    string runlistSiteNum[3]={"EH1","EH2","EH3"};
    vector<string> runlistFileContent[3];
    for( int i=0 ; i<3 ; i++ )
    {
        runlistFileContent[i].clear();
    }
    string verSuf;
    for( vector<string>::iterator it=dataVerReal.begin() ; it!=dataVerReal.end() ; ++it )
    {
        std::cout<<*it<<endl;
        verSuf+="_"+*it;
        //check data path
        string datadir="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/";
        DIR *dir=NULL;
        datadir+=*it;
        dir = opendir(datadir.c_str());
        if( dir==NULL )
        {
            std::cout<<" !!! data dir doesn't exist : "<<datadir<<endl;
            dataIsGood="0";
        } else
        {
            std::cout<<" data dir exist : "<<datadir<<endl;
            closedir(dir);
        }

        
        std::cout<<" "<<endl;
    }
        //check merged file 
        for( int i=0 ; i<3 ; i++ )
        {
            string runlistName=runlistSiteNum[i]+"totalHisto";
            runlistName+=verSuf;
            //runlistName+=*it;
            runlistName+=".root";
            TFile *ff = new TFile(runlistName.c_str());
            if( ff->IsZombie() )
            {
                std::cout<<" !!! merged file doesn't exist  : "<<runlistName<<endl;
                /*
                if( *it!="P12A" )
                {
                    dataIsGood="0";
                }
                */
            } else
            {
                std::cout<<" merged file exist  : "<<runlistName<<endl;
            }
            ff->Close();
        }
    verSuf=verSuf.substr(1);
    std::cout<<"data version : "<<verSuf<<endl;

    dataVerReal.push_back(verSuf);
    dataVerReal.push_back(dataIsGood);
    
    return dataVerReal; 
}
int IsoNumFit(string dataVer,int siteNum=0)
{
    vector<string> dataVerVec=checkdata(dataVer);
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    std::cout<<" "<<endl;
    dataVer=dataVerVec[dataVerVec.size()-2];
    if( dataVerVec[dataVerVec.size()-1]=="0" )
    {
        std::cout<<"  The dataVerion["<<dataVer<<"] or data is wrong ,please check! "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        return 0;
    }

    if( siteNum==0 )
    {
        for( int i=1 ; i<=3 ; i++ )
        {
            
            std::cout<<"====> begin to analyse EH"<<i<<"'s DaqTime, He8/Li9 "<<endl;
            std::cout<<"dataVersion  : "<<dataVer<<endl;
            fitHisto(i,dataVer);
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;

            std::cout<<"====> begin to analyse EH"<<i<<"'s IbdNum ,fast neutron"<<endl;
            //anaTree(i,dataVerVec);
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
            std::cout<<" "<<endl;
        }
        
    }else if(siteNum>=1 && siteNum<=3)
    {
    
        std::cout<<"====> begin to analyse EH"<<siteNum<<"'s DaqTime,He8/Li9 "<<endl;
        std::cout<<"dataVersion  : "<<dataVer<<endl;
        fitHisto(siteNum,dataVer);
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;

        std::cout<<"====> begin to analyse EH"<<siteNum<<"'s IbdNum ,fast neutron"<<endl;
        //anaTree(siteNum,dataVerVec);
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
        std::cout<<" "<<endl;
    }
    std::cout<<"ALL DONE!!! "<<endl;
    return 0;
}
