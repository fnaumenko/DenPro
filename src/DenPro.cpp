/*
 * DensPro
 * Copyright (C) 2014 Fedor Naumenko
*/

#include "DenPro.h"
#include <fstream>	// for dout

using namespace std;

const string Product::Title = "DenPro";
const string Product::Version = "1.0";
const string Product::Descr = "Density Profile";

const string OutFile = string(Product::Title) +  "_out.txt";
const string HelpOutFile = "duplicate standart output to " + OutFile + " file";
const string Alignment = "alignment";

enum eOptGroup	{ oINPUT, oTREAT, oOUTPUT, oOTHER };	// oOTHER should be the last 
const BYTE	Options::_GroupCount = oOTHER + 1;

const char* Options::_OptGroups [] = {
	"Input", "Treatment", "Ambig output", "Other"
};

Options::Option Options::_Options [] = {
	{ 'g', "gen",	1, true, tNAME, oINPUT, vUNDEF, 0, 0, NULL,
	"genome size file, or genome library, or single nucleotide sequence" },
	{ HPH, "gap-len",0,true, tINT,	oINPUT, 1000, 10, 100000, NULL,
	"minimal length of undefined nucleotides region in genome\nwhich is declared as a gap.\nIgnored for genome size file" },
	{ 'd', "dupl",	0, true, tENUM,	oINPUT, TRUE,	0, 2, (char*)Options::Booleans,
	"accept duplicate reads" },
	{ HPH, "diff-sz",	0, true, tENUM,	oINPUT, FALSE,	0, 2, (char*)Options::Booleans,
	"allow to ignore reads with different size" },
	{ 'c',Chrom::Abbr,	0,true,tCHAR, oTREAT, vUNDEF, 0, 0, NULL,
	"treat stated chromosome only (all)" },
	{ HPH, "min-scr",	0, true, tINT, oTREAT, vUNDEF, 0, 1000, NULL, "score threshold for treated reads (lack)" },
	{ HPH, "cons",	0, true, tINT, oTREAT, vUNDEF, 2, 500, NULL, "step of number of consolidated reads" },
	{ 'f', "fbed",	0, true, tNAME, oTREAT, vUNDEF,	0, 0, NULL,
	"'template' bed file which features define given regions" },
	{ 'e',"exp-len",0, true, tINT,	oTREAT, 0, 0, 1e3, NULL,
	"length of expanding features in 'template' bed file" },
	{ 's',"space",	0, true, tINT,	oTREAT, 50, 1, 1e4, NULL,
	"resolution: span in bp in which reads will be counted\nto define a density" },
	{ HPH, "alarm",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL,
	"output features ambiguities, if they exist" },
	{ HPH, "stat",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL,
	"output features ambiguities statistics, if they exist" },
	//{ 'o', "out",	0, false,tENUM, oOUTPUT,FALSE,	0, 2, NULL, HelpOutFile.c_str() },
	{ 't', "time",	0, false,tENUM, oOTHER,	FALSE,	0, 2, NULL, "output run time" },
	{ 'v', Version,	0, false,tVERS,	oOTHER, vUNDEF, 0, 0, NULL, "print program's version and quit" },
	{ 'h', "help",	0, false,tHELP,	oOTHER, vUNDEF, 0, 0, NULL, "print usage information and quit" }
};

const BYTE	Options::_OptCount = oHELP + 1;
const BYTE	Options::_UsageCount = 1;		// count of 'Usage' variants in help
const Options::Usage Options::_Usages[] = {	// content of 'Usage' variants in help
	{	vUNDEF,	" sequence"	}
};

ofstream outfile;			// file ostream duplicated cout; inizialised by file in code
dostream dout(cout, outfile);	// stream's duplicator

/*****************************************/
int main(int argc, char* argv[])
{
	if (argc < 2)	
	{ Options::PrintUsage(false);	return 0; }		// output tip

	short fileInd = Options::Tokenize(argc, argv, Alignment.c_str());
	if( fileInd < 0 )	return 1;		// wrong otpion
	int ret = 0;						// main() return code
	
	BedF *bedF = NULL;
	bool getTime = Options::GetBVal(oTIME);
	outfile.open(OutFile.c_str());
	Timer timer(getTime);
	timer.Start();
	try {
		// check file names first of all
		const char *genName  = FS::CheckedFileDirName(oGFILE);
		const char* rBedName = FS::CheckedFileName(argv[fileInd]);
		const char* fBedName = Options::GetSVal(oFBED);
		bool alarm	= Options::GetBVal(oALARM);		// print ambiguity messages
		bool stats	= Options::GetBVal(oSTATS);		// print ambiguity statistics
		chrid cID	= Chrom::ID(Options::GetSVal(oCHROM));
		Timer timer = Timer(Options::GetBVal(oTIME));

		if( fBedName ) {
			dout << Template << MSGSEP_BLANK;
			bedF = new BedF(FS::CheckedFileName(fBedName), cID, true, true, alarm, stats);
			bedF->Expand(Options::GetIVal(oEXPLEN), false);
		}
		BedR bedR((Alignment + MSGSEP_BLANK).c_str(), rBedName, cID, true, getTime, true, alarm, stats,
			Options::GetBVal(oDUPL), Options::GetBVal(oDIFFSZ), Options::GetIVal(oMINSCR));

		if( bedF && !SetCommonChroms(bedR, *bedF, false) )
			Err("no common chromosomes", "").Throw();
		dout << EOL;
		GenomeRegions grgn(genName, &cID, Options::GetIVal(oGAPLEN));
		DenPro(cID, bedR, grgn, bedF);
	}
	catch(Err &e)		{ ret = 1;	dout << e.what() << EOL; }
	catch(exception &e)	{ ret = 1;	dout << e.what() << EOL; }
	catch(...)			{ ret = 1;	dout << "Unregistered error\n"; }
	if( bedF )	delete bedF;

	timer.Stop(true);
	if( outfile.is_open() )		outfile.close();	// in case of holding execution by user
//#ifdef OS_Windows
//	system("pause");
//#endif
	return ret;
}

/************************ class DenPro ************************/
chrlen	DenPro::WinLen = 0;
//const char* sRatio = "ratio";
//const string sDensUnit = sReads + "/kbs";

DenPro::DenPro(chrid cID, BedR &bedR, GenomeRegions &gRegn, BedF *bedF)
{
	WinLen = Options::GetIVal(oSPACE);
	GenomeReadDistrib grDist(cID, bedR, gRegn, bedF);
	grDist.Scan(WinLen);
	grDist.PrintDensity();

	// fill dens distribution
	ULONG rCnt = 0;
	if( bedF )  rCnt += Init(grDist, PairReadDistrib::IN_P);
	rCnt += Init(grDist, PairReadDistrib::OUT_P);
	dout << "\nDistributed " << BedR::ReadTitle(true) << MSGSEP_BLANK
			<< rCnt << sPercent(rCnt, ULLONG(bedR.ReadsCount(cID)), 3) << EOL;

	Consolidate(Options::GetIVal(oCONS));
	Print(outfile);
}

// Initializes precise Read's distribution
//	@grDistr: genome reads distribution
//	@loc: In-|Out- features location
//	return: count of distributed Reads
ULONG DenPro::Init(GenomeReadDistrib& grDistr, PairReadDistrib::eLocating loc)
{
	ULONG rCnt, sumrCnt = 0;
	const ReadDistrib& rDistr = grDistr.Distrib(loc);

	_densDistr[loc].reserve(rDistr.Length());
	for(chrlen ri=0; ri<rDistr.Length(); ri++) {
		AddVal(loc, rCnt=rDistr[ri]);
		sumrCnt += rCnt;
	}
	sort(_densDistr[loc].begin(), _densDistr[loc].end());
	return sumrCnt;
}

// Adds value to the distribution: increase counter if value exists or add new counter.
//	@densDistr: adding distribution
//	@cntReads: added value
void DenPro::AddVal(BYTE loc, chrlen rCnt)
{
	vector<DensDistr>::iterator iter = _densDistr[loc].begin(); 
	for(; iter != _densDistr[loc].end(); iter++)
		if( iter->CountReads == rCnt ) {
			(iter->CountWins)++;
			return;
		}
	_densDistr[loc].push_back( DensDistr(rCnt) );		// if val not found
}

// Split sets of pairs <ReadCnt>-<WinCnt> into subset with defined length of ReadCnt 
//	@consLen: step of number of consolidated reads or vUNDEF if no consolidation
void DenPro::Consolidate(int consLen)
{
	if( consLen == vUNDEF )	return;

	ULONG sumR;		// total count of Reads
	//float sumR;
	chrlen size;

	//dout << "consolidate\n";
	for(BYTE i=0; i<2; i++) {
		size = _densDistr[i].size();
		sumR = 0;
		for(chrlen k=1; k<size; k++)
			sumR += _densDistr[i][k].CountReads * _densDistr[i][k].CountWins;
		if( consLen ) {
			ULONG sumW = 0;
			chrlen range;
			_consDensDistr[i].reserve(size/consLen);
			for(chrlen j=1, x=1; j<size; x++) {
				range = x*consLen;
				for(sumW=sumR=0; true; j++) {
					if( j>=size )			break;
					sumR += _densDistr[i][j].CountReads;
					if( sumR >= range )		break;
					sumW += _densDistr[i][j].CountWins;
				}
				_consDensDistr[i].push_back(DensDistr(sumW, x*consLen));
			}
		}
	}
	//Print(dout);
}

// Outputs to stream Read's distribution
//	@stream: outputted stream. Currently not used.
//	@printTwoDistribs: true if template is setting
void DenPro::Print(ofstream& stream)	//, bool printTwoDistribs)
{
	chrlen	cnt1, j;
	
	// print precise distribution
	if( !_densDistr[PairReadDistrib::IN_P].size() ) {	// no template?
		//stream << "Precise distribution:\n";
		//stream << "ReadCnt\tWinCnt\n";
		dout << "Precise distribution:\n";
		dout << "ReadCnt\tWinCnt\n";
		cnt1 = _densDistr[PairReadDistrib::OUT_P].size();
		for(j=0; j<cnt1; j++)
			//stream	<< _densDistr[PairReadDistrib::OUT_P][j].CountReads << TAB
			dout	<< _densDistr[PairReadDistrib::OUT_P][j].CountReads << TAB
					<< _densDistr[PairReadDistrib::OUT_P][j].CountWins << EOL;
	}

	// print consolidate distribution
	cnt1 = _consDensDistr[PairReadDistrib::OUT_P].size();
	if( !cnt1 )	return;			// no consolidation
	chrlen cnt0 = _consDensDistr[PairReadDistrib::IN_P].size();
	//stream << "Consolidate distribution:\n";
	//stream << "ReadCnt\t";
	//if( cnt0 )	stream << "inWinCnt\toutWinCnt\n";
	//else		stream <<  "WinCnt\n";
	dout << "Consolidate distribution:\nReadCnt\t";
	if( cnt0 )	dout << "inWinCnt\toutWinCnt\n";
	else		dout <<  "WinCnt\n";
	//stream << cnt0 ? "inWinCnt\toutWinCnt\n" : "WinCnt\n";
	for(j=0; j<cnt0; j++) {
		//stream	<< _consDensDistr[PairReadDistrib::IN_P][j].CountReads << TAB 
		dout	<< _consDensDistr[PairReadDistrib::IN_P][j].CountReads << TAB 
				<< _consDensDistr[PairReadDistrib::IN_P][j].CountWins << TAB;
		if( j<cnt1 )
			//stream << _consDensDistr[PairReadDistrib::OUT_P][j].CountWins;
			dout << _consDensDistr[PairReadDistrib::OUT_P][j].CountWins;
		//stream << EOL;
		dout << EOL;
	}
	for(; j<cnt1; j++)
		//stream	<< _consDensDistr[PairReadDistrib::OUT_P][j].CountReads << TAB 
		dout	<< _consDensDistr[PairReadDistrib::OUT_P][j].CountReads << TAB 
				//<< TAB 
				<< _consDensDistr[PairReadDistrib::OUT_P][j].CountWins << EOL;
}

//void DenPro::Write()
//{
//	ofstream outFile;
//
//	outFile.open(OutFile.c_str());
//	Print(outFile);
//	outFile.close();
//}

//void DenPro::Write(const char * bedFfName, chrlen fCnt, const char * bedRfName, chrlen rCnt)
//{
//	ofstream outFile;
//	const char* fNm = strrchr(bedRfName, SLASH)+1;	// short file name
//	string outfname = string(fNm).substr(0, strchr(fNm,'.')-fNm) + "_dens.txt";;
//
//	cout << "write to " << outfname << EOL;
//	outFile.open(outfname.c_str());
//	outFile << "align: " << bedRfName << MSGSEP_TAB << rCnt << " reads\n";
//	if( bedFfName )		outFile << "templ: " << bedFfName << MSGSEP_TAB << fCnt << " features\n";
//	outFile << "winLen: " << WinLen << "\twinShift: " << WinShift << EOL;
//	ULLONG rTotalCnt = TotalReadCount();
//	outFile << "total recorded reads: " << rTotalCnt << TAB << 100 * rTotalCnt / rCnt << "%\n";
//	PrintDenRatio(outFile);
//	Print(outFile);
//	outFile.close();
//}

/************************ end of class DenPro ************************/

//
///************************ class ReadDistrib ************************/
//
//// Prepares to scanning: reserves capacity and zero counters.
////	@regns: Regions for wich distribution is defined
////	@owner: owner (holder) of this distribution
//void ReadDistrib::Init(Regions* const regns, PairReadDistrib* const owner)
//{
//	_rtotalCnt = _ri = _wi = 0;
//	_regns = regns;
//	_owner = owner;
//	_wCurrLen = _owner->WinLength();
//	Reserve(regns->Length()/_wCurrLen); // + 1);
//}
//
//// Returns count of Reads in window
////	@start: window's start position
////	return: number of Reads in window
//chrlen ReadDistrib::ScanWindow(chrlen start)
//{
//	chrlen	rCentre, resCnt = 0, rCnt = _owner->ChromReadsCount();
//
//	_outWin = false;
//	for(; _ri<rCnt; _ri++) {		// loop through window for Reads
//		rCentre = _owner->ReadCentrePos(_ri);
//		if( rCentre >= start ) {				// pass left-of-window Reads
//			if( rCentre >= start + _wCurrLen ) {
//				_outWin = true;
//				return resCnt;					// exit on first right-of-window Read
//			}
//			resCnt++;
//			_rtotalCnt++;
//		}
//	}
//	return resCnt;
//}
//
//// Returns count of Reads in Region
//void ReadDistrib::ScanRegion(const Region& rgn)
//{
//	chrlen start = rgn.Start;
//	chrlen end = rgn.End;
//	chrlen len = static_cast<chrlen>(Length());
//
//	chrlen res;
//
//	_outWin = false;
//	if(start + _wCurrLen <= end) {		// current window belong to region entirely
//		do {
//			if( _wi+1 == len )	return;
//			res = ScanWindow(start);
//			_data[_wi] += res;
//			_wi++;
//			//if( _wi == static_cast<chrlen>(Length()) )
//			//	Err("window's index is out of range", "ReadDistrib::ScanRegion").Throw();
//			if( _outWin ) {
//				_data[_wi]++;
//				_rtotalCnt++;
//				_ri++;
//			}
//			start += _wCurrLen;
//			_wCurrLen = _owner->WinLength();	// set user's win length
//		}
//		while(start + _wCurrLen <= end);
//		
//		if(start + _wCurrLen == end)	// window's- and region's ends are equal
//			return;
//	}
//	// current window exceeds region
//	_wCurrLen -= end - start;			// decrease current window by the rest of region
//	_data[_wi] = ScanWindow(start);		// scan part of current window
//}
//
//// Counts Reads in defined Regions through the chromosome
//void ReadDistrib::Scan()
//{
//	for(Regions::Iter it=_regns->Begin(); it!=_regns->End(); it++)
//		ScanRegion(*it);
//}
//
//#ifdef DEBUG
//void ReadDistrib::Print(ostream& stream) const
//{
//	bool zeroflag = false;
//	stream << "Count of reads in cutting windows:\n";
//	stream << "window\treads\n";
//	long wCnt = Length();
//	for(long i=0; i<wCnt; i++) {
//		if( _data[i] == 0 )
//			stream << i << TAB << _data[i] << EOL;
//		//if( _data[i] != 0 ) {
//		//	stream << i << TAB << _data[i] << EOL;
//		//	zeroflag = false;
//		//}
//		//else
//		//	if( !zeroflag )	{
//		//		stream << "...\n";
//		//		zeroflag = true;
//		//	}
//	}
//}
//#endif
//
///************************ end of class ReadDistrib ************************/
//
///************************ class PairReadDistrib ************************/
//
//// Creates new instance by chrom, BedR, and regions
////	@cID: chromosome's ID
////	@bedR: BedR for what distributions should be builded
////	@templRegns: external additional template Regions from BedF or from whole chromosome.
////	For In-peak distribution has been used directly,
////	for Out-peak has been inverted and ovarlaped with chrom defined Regions.
//PairReadDistrib::PairReadDistrib(chrid cID, const BedR &bedR, const Regions &defRegns, const BedF *bedF) :
//	_halfRLen( bedR.ReadLen()/2 ),
//	_beginIt( bedR.ReadsBegin(cID) ),
//	_rCnt( bedR.ReadsCount(cID) )
//{
//	if( bedF ) {
//		// fill _templRegns[IN_P] by Regions-features directly
//		bedF->FillRegions(cID, _templRegns[IN_P]);
//		// fill _templRegns[OUT_P] by Regions-feature, inverting and and overlaping with defRegns 
//		Regions tmpRegns;
//		tmpRegns.FillInvert(_templRegns[IN_P], defRegns.LastEnd());
//		_templRegns[OUT_P].FillOverlap(defRegns, tmpRegns);
//	}
//	else// fill _templRegns[IN_P] by defRegns
//		_templRegns[OUT_P].Copy(defRegns);
//}
//
///************************ end of class PairReadDistrib ************************/
//
///************************ class GenomeReadDistrib ************************/
//
//#define DENS_BASE	1000	// base on wich density is defined
//const char* sRatio = "ratio";
//const string sDensUnit = "rd/kbs";
//
//// Initializes In-|Out-Read's distribution container:
//// for each chromosome initializes default Regions and PairReadDistrib.
//// Only chromosomes marked as 'Treated' would be treated.
////	@cID: chromosome's ID
////	@bedR: Reads wich are distribeted
////	@gName: chrom.sizes file or genome directory, determining define Regions
////	@bedF: features determining peaks or NULL
//GenomeReadDistrib::GenomeReadDistrib (
//	chrid cID, const BedR & bedR, 
//	GenomeRegions& gRegn,
//	const BedF* bedF) 
//	: _twoDistrs(bedF ? true : false)
//{
//	for(BedR::cIter it = bedR.cBegin(); it != bedR.cEnd(); it++)
//		if(TREATED(it))
//			AddClass(CID(it), PairReadDistrib(CID(it), bedR, gRegn[CID(it)], bedF));
//}
//
//// Initializes and fills In-|Out-Read's distribution for each chromosome in genome
////	@winLen: length of layout (cutting) window
//GenomeReadDistrib& GenomeReadDistrib::Scan(chrlen winLen)
//{
//	_wLen = winLen;
//	// loop through chromosomes
//	for(Iter it = Begin(); it != End(); it++) {
//		if( _twoDistrs ) {
//			it->second.Init(PairReadDistrib::IN_P, winLen);
//			it->second.Scan(PairReadDistrib::IN_P);
//		}
//		it->second.Init(PairReadDistrib::OUT_P, winLen);
//		it->second.Scan(PairReadDistrib::OUT_P);
//	}
//	return *this;
//}
//
//// Returns density, averaged for whole genome
////	@loc: In-|Out- location
//float	GenomeReadDistrib::Density (PairReadDistrib::eLocating loc) const
//{
//	if( loc == PairReadDistrib::IN_P && !_twoDistrs )
//		return -1;	// only out_peak distribution is defined
//		float res = 0;
//		for(cIter it = cBegin(); it != cEnd(); it++)
//			res += it->second.Density(loc);
//		return DENS_BASE * res;
//}
//
//// Returns total number of Reads in n-|Out- distribution
////	@loc: In-|Out- location
//ULONG	GenomeReadDistrib::ReadsCount(PairReadDistrib::eLocating loc) const
//{
//	if( loc == PairReadDistrib::IN_P && _twoDistrs )
//		return 0;		// only out_peak distribution is defined
//	ULONG res = 0;
//	for(cIter it = cBegin(); it != cEnd(); it++)
//		res += it->second.ReadsCount(loc);
//	return res;
//}
//
//// Returns Read's distribution after scanning.
////	@loc: In-|Out- location
//ReadDistrib & GenomeReadDistrib::Distrib(PairReadDistrib::eLocating loc)
//{
//	chrlen len = 0;
//	Iter it;
//	for(it = Begin(); it != End(); it++)
//		len += it->second.Distrib(loc).Length();
//	_totalReadDistr.Reserve(len);
//	len = 0;
//	for(it = Begin(); it != End(); it++) {
//		_totalReadDistr.ConcatDistrib(it->second.Distrib(loc), len);
//		len += it->second.Distrib(loc).Length();
//	}
//	return _totalReadDistr;
//}
//
//void GenomeReadDistrib::PrintDensity()
//{
//	float	inDens = 0, outDens;
//	ULONG	rInCnt = 0, rOutCnt = 0,	// total count of reads in locations
//			wInCnt = 0, wOutCnt = 0;	// total count of windows in locations
//
//	// header
//	dout << "Mean density, " << sDensUnit << EOL;
//	if(_twoDistrs)	dout << "\tIn-peaks\tOut-peaks\t" << sRatio << EOL;
//
//	for(cIter it=cBegin(); it!=cEnd(); it++) {
//		dout << Chrom::AbbrName(CID(it)) << TAB;
//		if(_twoDistrs) {
//			inDens = DENS_BASE * it->second.Density(PairReadDistrib::IN_P);
//			dout << inDens << TAB << TAB;
//			rInCnt += it->second.ReadsCount(PairReadDistrib::IN_P);
//			wInCnt += it->second.WinsCount(PairReadDistrib::IN_P);
//		}
//		outDens = DENS_BASE * it->second.Density(PairReadDistrib::OUT_P);
//		dout << outDens;
//
//		//dout << EOL;
//		//it->second._rDistribs[PairReadDistrib::OUT_P].Print();
//
//		rOutCnt += it->second.ReadsCount(PairReadDistrib::OUT_P);
//		wOutCnt += it->second.WinsCount(PairReadDistrib::OUT_P);
//		if(_twoDistrs)	dout << TAB << TAB << inDens/outDens;
//		dout << EOL;
//	}
//	if( ChromsCount() > 1 ) {
//		dout << Total << ":\t";
//		if(_twoDistrs) {
//			inDens = DENS_BASE * float(rInCnt) / wInCnt / _wLen;
//			dout << inDens << TAB << TAB;
//		}
//		outDens = DENS_BASE * float(rOutCnt) / wOutCnt / _wLen;
//		dout << outDens;
//		if(_twoDistrs)	dout << TAB << TAB << inDens/outDens;
//		dout << EOL;
//	}
//}
//
///************************ end of class GenomeReadDistrib ************************/
//
