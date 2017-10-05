#include "Distrib.h"

/************************ class ReadDistrib ************************/

// Prepares to scanning: reserves capacity and zero counters.
//	@regns: Regions for wich distribution is defined
//	@owner: owner (holder) of this distribution
void ReadDistrib::Init(Regions* const regns, PairReadDistrib* const owner)
{
	_rtotalCnt = _ri = _wi = 0;
	_regns = regns;
	_owner = owner;
	_wCurrLen = _owner->WinLength();
	Reserve(regns->Length()/_wCurrLen); // + 1);
}

// Returns count of Reads in window
//	@start: window's start position
//	@outWin: returned true if current Read passed left boundary of window, otherwise false
//	return: number of Reads in window
chrlen ReadDistrib::ScanWindow(chrlen start, bool* outWin)
{
	chrlen	rCentre, resCnt = 0, rCnt = _owner->ChromReadsCount();

	*outWin = false;
	for(; _ri<rCnt; _ri++) {		// loop through window for Reads
		rCentre = _owner->ReadCentrePos(_ri);
		if( rCentre >= start ) {				// pass left-of-window Reads
			if( rCentre >= start + _wCurrLen ) {
				*outWin = true;
				return resCnt;					// exit on first right-of-window Read
			}
			resCnt++;
			_rtotalCnt++;
		}
	}
	return resCnt;
}

// Returns count of Reads in Region
void ReadDistrib::ScanRegion(const Region& rgn)
{
	chrlen start = rgn.Start;
	chrlen end = rgn.End;
	chrlen len = static_cast<chrlen>(Length());
	bool outWin;
	bool* poutWin = &outWin;
	//chrlen res;

	while(start + _wCurrLen <= end) {	// current window belong to region entirely
		if( _wi+1 == len )	return;
		//res = ScanWindow(start);
		_data[_wi] += ScanWindow(start, poutWin);
		_wi++;
		//if( _wi == static_cast<chrlen>(Length()) )
		//	Err("window's index is out of range", "ReadDistrib::ScanRegion").Throw();
		if( outWin ) {
			_data[_wi]++;
			_rtotalCnt++;
			_ri++;
		}
		//else if( res == 0 )		dout << "0:\t" << _wi << TAB << start << EOL;
		start += _wCurrLen;
		_wCurrLen = _owner->WinLength();	// set user's win length
	}
	if(start + _wCurrLen > end) {	// current window exceeds region
		_wCurrLen -= end - start;			// decrease current window by the rest of region
		_data[_wi] += ScanWindow(start, poutWin);	// scan part of current window
	}
}

// Counts Reads in defined Regions through the chromosome
void ReadDistrib::Scan()
{
	for(Regions::Iter it=_regns->Begin(); it!=_regns->End(); it++)
		ScanRegion(*it);
}

#ifdef DEBUG
void ReadDistrib::Print(ostream& stream) const
{
	bool zeroflag = false;
	stream << "Count of reads in cutting windows:\n";
	stream << "window\treads\n";
	long wCnt = Length();
	for(long i=0; i<wCnt; i++) {
		if( _data[i] == 0 )
			stream << i << TAB << _data[i] << EOL;
		//if( _data[i] != 0 ) {
		//	stream << i << TAB << _data[i] << EOL;
		//	zeroflag = false;
		//}
		//else
		//	if( !zeroflag )	{
		//		stream << "...\n";
		//		zeroflag = true;
		//	}
	}
}
#endif

/************************ end of class ReadDistrib ************************/

/************************ class PairReadDistrib ************************/

// Creates new instance by chrom, BedR, and regions
//	@cID: chromosome's ID
//	@bedR: BedR for what distributions should be builded
//	@templRegns: external additional template Regions from BedF or from whole chromosome.
//	For In-peak distribution has been used directly,
//	for Out-peak has been inverted and ovarlaped with chrom defined Regions.
PairReadDistrib::PairReadDistrib(chrid cID, const BedR &bedR, const Regions &defRegns, const BedF *bedF) :
	_halfRLen( bedR.ReadLen()/2 ),
	_beginIt( bedR.ReadsBegin(cID) ),
	_rCnt( bedR.ReadsCount(cID) )
{
	if( bedF ) {
		// fill _templRegns[IN_P] by Regions-features directly
		bedF->FillRegions(cID, _templRegns[IN_P]);
		// fill _templRegns[OUT_P] by Regions-feature, inverting and and overlaping with defRegns 
		Regions tmpRegns;
		tmpRegns.FillInvert(_templRegns[IN_P], defRegns.LastEnd());
		_templRegns[OUT_P].FillOverlap(defRegns, tmpRegns);
	}
	else// fill _templRegns[IN_P] by defRegns
		_templRegns[OUT_P].Copy(defRegns);
}

/************************ end of class PairReadDistrib ************************/

/************************ class GenomeReadDistrib ************************/

#define DENS_BASE	1000	// base on wich density is defined
const char* sRatio = "ratio";
const string sDensUnit = "rd/kbs";

// Initializes In-|Out-Read's distribution container:
// for each chromosome initializes default Regions and PairReadDistrib.
// Only chromosomes marked as 'Treated' would be treated.
//	@bedR: Reads wich are distribeted
//	@gName: chrom.sizes file or genome directory, determining define Regions
//	@bedF: features determining peaks or NULL
GenomeReadDistrib::GenomeReadDistrib (
	const BedR & bedR, 
	GenomeRegions& gRegn,
	const BedF* bedF) 
	: _twoDistrs(bedF ? true : false)
{
	for(BedR::cIter it = bedR.cBegin(); it != bedR.cEnd(); it++)
		if(TREATED(it))
			AddClass(CID(it), PairReadDistrib(CID(it), bedR, gRegn[CID(it)], bedF));
}

// Initializes and fills In-|Out-Read's distribution for each chromosome in genome
//	@winLen: length of layout (cutting) window
GenomeReadDistrib& GenomeReadDistrib::Scan(chrlen winLen)
{
	_wLen = winLen;
	// loop through chromosomes
	for(Iter it = Begin(); it != End(); it++) {
		if( _twoDistrs ) {
			it->second.Init(PairReadDistrib::IN_P, winLen);
			it->second.Scan(PairReadDistrib::IN_P);
		}
		it->second.Init(PairReadDistrib::OUT_P, winLen);
		it->second.Scan(PairReadDistrib::OUT_P);
	}
	return *this;
}

// Returns density, averaged for whole genome
//	@loc: In-|Out- location
float	GenomeReadDistrib::Density (PairReadDistrib::eLocating loc) const
{
	if( loc == PairReadDistrib::IN_P && !_twoDistrs )
		return -1;	// only out_peak distribution is defined
		float res = 0;
		for(cIter it = cBegin(); it != cEnd(); it++)
			res += it->second.Density(loc);
		return DENS_BASE * res;
}

// Returns total number of Reads in n-|Out- distribution
//	@loc: In-|Out- location
ULONG	GenomeReadDistrib::ReadsCount(PairReadDistrib::eLocating loc) const
{
	if( loc == PairReadDistrib::IN_P && _twoDistrs )
		return 0;		// only out_peak distribution is defined
	ULONG res = 0;
	for(cIter it = cBegin(); it != cEnd(); it++)
		res += it->second.ReadsCount(loc);
	return res;
}

// Returns Read's distribution after scanning.
//	@loc: In-|Out- location
ReadDistrib & GenomeReadDistrib::Distrib(PairReadDistrib::eLocating loc)
{
	chrlen len = 0;
	Iter it;
	for(it = Begin(); it != End(); it++)
		len += it->second.Distrib(loc).Length();
	_totalReadDistr.Reserve(len);
	len = 0;
	for(it = Begin(); it != End(); it++) {
		_totalReadDistr.ConcatDistrib(it->second.Distrib(loc), len);
		len += it->second.Distrib(loc).Length();
	}
	return _totalReadDistr;
}

void GenomeReadDistrib::PrintDensity()
{
	float	inDens = 0, outDens;
	ULONG	rInCnt = 0, rOutCnt = 0,	// total count of reads in locations
			wInCnt = 0, wOutCnt = 0;	// total count of windows in locations

	// header
	dout << "Mean density, " << sDensUnit << EOL;
	if(_twoDistrs)	dout << "\tIn-peaks\tOut-peaks\t" << sRatio << EOL;

	for(cIter it=cBegin(); it!=cEnd(); it++) {
		dout << Chrom::AbbrName(CID(it)) << TAB;
		if(_twoDistrs) {
			inDens = DENS_BASE * it->second.Density(PairReadDistrib::IN_P);
			dout << inDens << TAB << TAB;
			rInCnt += it->second.ReadsCount(PairReadDistrib::IN_P);
			wInCnt += it->second.WinsCount(PairReadDistrib::IN_P);
		}
		outDens = DENS_BASE * it->second.Density(PairReadDistrib::OUT_P);
		dout << outDens;

		//dout << EOL;
		//it->second._rDistribs[PairReadDistrib::OUT_P].Print();

		rOutCnt += it->second.ReadsCount(PairReadDistrib::OUT_P);
		wOutCnt += it->second.WinsCount(PairReadDistrib::OUT_P);
		if(_twoDistrs)	dout << TAB << TAB << inDens/outDens;
		dout << EOL;
	}
	if( ChromsCount() > 1 ) {
		dout << Total << ":\t";
		if(_twoDistrs) {
			inDens = DENS_BASE * float(rInCnt) / wInCnt / _wLen;
			dout << inDens << TAB << TAB;
		}
		outDens = DENS_BASE * float(rOutCnt) / wOutCnt / _wLen;
		dout << outDens;
		if(_twoDistrs)	dout << TAB << TAB << inDens/outDens;
		dout << EOL;
	}
}

/************************ end of class GenomeReadDistrib ************************/

