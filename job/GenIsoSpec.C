#include  <iostream>
#include  <TChain.h>
#include  <TROOT.h>
#include  <string>
int anaTree(int siteNum,string dataVer)
{
	gROOT->ProcessLine(".L SingleTree.C+");
	TChain chain("Tree/SingleTree");
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
    chainFile+="bak7/";
    chainFile+=site;
    chainFile+="/*LiAna.root";
    chain.Add(chainFile.c_str());
    if( dataVer=="P12E" )
    {
        chainFile="/afs/ihep.ac.cn/users/l/lidj/largedata/LiAna/";
        chainFile+="P12A";
        chainFile+="bak7/";
        chainFile+=site;
        chainFile+="/*LiAna.root";
        chain.Add(chainFile.c_str());
    }
    string  siteAndDataVer=site+dataVer;
    std::cout<<"site And DataVer  : "<<siteAndDataVer<<endl;
	chain.Process("SingleTree",siteAndDataVer.c_str());
	return 0;
}
int GenIsoSpec(string dataVer,int sitenum=0)
{
    //string dataVer="P13A";
    if( !(dataVer=="P13A"||dataVer=="P12E") )
    {
        std::cout<<"  The dataVerion is wrong ,please check! "<<endl;
    }
    if( sitenum==0 )
    {
        anaTree(1,dataVer);
        anaTree(2,dataVer);
        anaTree(3,dataVer);
    }
    if( sitenum==1 )
    {
        anaTree(1,dataVer);
    }
    if( sitenum==2 )
    {
        anaTree(2,dataVer);
    }
    if( sitenum==3 )
    {
        anaTree(3,dataVer);
    }
    return 0;
}
