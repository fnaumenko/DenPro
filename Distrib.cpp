#include "Distrib.h"

/************************ class ReadDistrib ************************/

// Returns count of Reads in window
//	@wStart: current window's start position
//	@wLen: current window's length
//	@outWin: returned true if current Read passed left boundary of window, otherwise false
//	return: number of Reads in window
chrlen ReadDistrib::ScanWindow(chrlen wStart, chrlen wLen, bool* outWin)
{
	chrlen	rCentre, resCnt = 0, rCnt = _owner->ChromReadsCount();

	*outWin = false;
	for(; _ri<rCnt; _ri++) {		// loop through window for Reads
		rCentre = _owner->ReadCentrePos(_ri);
		if( rCentre >= wStart ) {				// pass left-of-window Reads
			if( rCentre >= wStart + wLen ) {
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
//	@rgn: current region
//	@wLen: length of window
void ReadDistrib::ScanRegion(const Region& rgn, chrlen wLen)
{
	chrlen start = rgn.Start;
	chrlen end = rgn.End;
	chrlen len = static_cast<chrlen>(Length());
	bool outWin;
	bool* poutWin = &outWin;

	while(start + wLen <= end) {	// current window belong to region entirely
		if( _wi+1 == len )	return;
		_data[_wi] += ScanWindow(start, wLen, poutWin);
		_wi++;
		//if( _wi == static_cast<chrlen>(Length()) )
		//	Err("window's index is out of range", "ReadDistrib::ScanRegion").Throw();
		if( outWin ) {
			_data[_wi]++;
			_rtotalCnt++;
			_ri++;
		}
		start += wLen;
	}
	if(start + wLen > end)		// current window exceeds region
		_data[_wi] += ScanWindow(start, wLen - end + start, poutWin);	// scan part of current window
}

// Counts Reads in given Regions through the chromosome
void ReadDistrib::Scan(const Regions& rgns)
{
	chrlen wLen = _owner->WinLength();
	_rtotalCnt = _ri = _wi = 0;

	Reserve(rgns.Length()/wLen); // + 1);
	for(Regions::Iter it=rgns.Begin(); it!=rgns.End(); it++)
		ScanRegion(*it, wLen);
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

// Initializes regions
//	@cID: chromosome's ID
//	@cRegns: defined reference chrom Regions
//	For In-peak distribution has been used directly,
//	for Out-peak has been inverted and ovarlaped with chrom defined Regions.
//	@bedF: optional additional regions
void PairReadDistrib::InitRegions(chrid cID, const Regions &cRegns, const BedF *bedF)
{
	if( bedF ) {
		// fill _regns[IN_P] by Regions-features directly
		bedF->FillRegions(cID, _regns[IN_P]);
		// fill _regns[OUT_P] by Regions-feature, inverting and and overlaping with defRegns 
		Regions tmpRegns;
		tmpRegns.FillInvert(_regns[IN_P], cRegns.LastEnd());
		_regns[OUT_P].FillOverlap(cRegns, tmpRegns);
		_rDistribs[IN_P].SetOwner(this);	// duplicates copy constructor; needed for _NO_UNODMAP build
	}
	else// fill _regns[IN_P] by cRegns
		_regns[OUT_P] = cRegns;
	_rDistribs[OUT_P].SetOwner(this);		// duplicates copy constructor; needed for _NO_UNODMAP build

}

float PairReadDistrib::GetResults(eLocating loc, ULONG* rCnt, ULONG* wCnt) const
{
	rCnt[loc] += ReadsCount(loc);
	wCnt[loc] += WinsCount(loc);
	return DENS_BASE * Density(loc);;
}

/************************ end of class PairReadDistrib ************************/

/************************ class GenomeReadDistrib ************************/

const char* sRatio = "ratio";
const string sDensUnit = "rd/kbs";

// Initializes In-|Out-Read's distribution container:
// for each chromosome initializes default Regions and PairReadDistrib.
// Only chromosomes marked as 'Treated' would be treated.
//	@bedR: Reads wich are distribeted
//	@gRegn: defined reference genome Regions
//	@bedF: features determining peaks or NULL
GenomeReadDistrib::GenomeReadDistrib (
	const BedR& bedR, 
	GenomeRegions& gRegn,
	const BedF* bedF) 
	: _twoDistrs(bedF ? true : false)
{
	Reserve(bedR.ChromsCount());
	for(BedR::cIter it = bedR.cBegin(); it != bedR.cEnd(); it++)
		if(TREATED(it))
			AddClass(CID(it), PairReadDistrib(CID(it), bedR)).
				InitRegions(CID(it), gRegn[CID(it)], bedF);		// initialize Regions after adding to collection,
																// to avoid regions copying
}

// Initializes and fills In-|Out-Read's distribution for each chromosome in genome
//	@winLen: length of layout (cutting) window
GenomeReadDistrib& GenomeReadDistrib::Scan(chrlen winLen)
{
	_wLen = winLen;
	// loop through chromosomes
	for(Iter it = Begin(); it != End(); it++) {
		if( _twoDistrs )
			it->second.Scan(PairReadDistrib::IN_P, winLen);
		it->second.Scan(PairReadDistrib::OUT_P, winLen);
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
	chrid	cID;
	float	inDens, outDens;
	ULONG	rCnt[2], wCnt[2];	// total count of reads in locations, total count of windows in locations
	memset(rCnt, 0, 2*sizeof(ULONG));
	memset(wCnt, 0, 2*sizeof(ULONG));

	// header
	dout << "Mean density, " << sDensUnit << EOL;
	if(_twoDistrs)	dout << "\tIn-peaks\tOut-peaks\t" << sRatio << EOL;

#ifdef _NO_UNODMAP
	for(cIter it=cBegin(); it!=cEnd(); it++) {
		cID = CID(it);
		const PairReadDistrib& prd = it->second;
#else
	// sort printed chroms
	vector<chrid> cids;
	cids.reserve(ChromsCount());
	for(cIter it=cBegin(); it!=cEnd(); it++)	cids.push_back(CID(it));
	sort(cids.begin(), cids.end());

	for(vector<chrid>::iterator it=cids.begin(); it!=cids.end(); it++) {
		const PairReadDistrib& prd = At(cID = *it);
#endif	// _NO_UNODMAP
		dout << Chrom::AbbrName(cID) << TAB;
		if(_twoDistrs)
			dout << (inDens = prd.GetResults(PairReadDistrib::IN_P, rCnt, wCnt)) << TAB << TAB;
		dout << (outDens = prd.GetResults(PairReadDistrib::OUT_P, rCnt, wCnt));
		if(_twoDistrs)	dout << TAB << TAB << inDens/outDens;
		dout << EOL;

	}
	if( ChromsCount() > 1 ) {
		dout << Total << SepClTab;
		if(_twoDistrs)	dout << (inDens = GetTotalDensity(0, rCnt, wCnt)) << TAB << TAB;
		dout << (outDens = GetTotalDensity(1, rCnt, wCnt));
		if(_twoDistrs)	dout << TAB << TAB << inDens/outDens;
		dout << EOL;
	}
}

/************************ end of class GenomeReadDistrib ************************/

