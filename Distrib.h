#pragma once
#include "Data.h"

#define DENS_BASE	1000	// base on wich density is defined

class PairReadDistrib;		// forward declaration used in ReadDistrib

class ReadDistrib : public Array<chrlen>
/*
 * 'ReadDistrib' represents simple distribution of Reads,
 * that's a number of Reads in each cutting window. Windows are indexed from 0 to maxWinCnt
 */
{
	// members to keep values during discontinuous scanning
	chrlen	_ri,			// current read's index in BedR
			_wi;			// current window's index
	// summerised members
	ULONG	_rtotalCnt;		// total number of reads in distribution
	// members - constants in this instance 
	PairReadDistrib *_owner;// Container of this instance
	
	// Returns count of Reads in window
	//	@wStart: current window's start position
	//	@wLen: current window's length
	//	@outWin: returned true if current Read passed left boundary of window, otherwise false
	//	Return: count of Reads in window
	chrlen ScanWindow(chrlen wStart, chrlen wLen, bool* outWin);

	// Returns count of Reads in Region
	//	@rgn: current region
	//	@wLen: length of window
	//	Return: count of Reads in Region
	void ScanRegion(const Region&, chrlen wLen);

public:
	//inline ReadDistrib() { _rtotalCnt = _wi = _ri = 0; _owner = NULL; }

	// Set owner (holder) of this distribution
	inline void SetOwner(PairReadDistrib* owner) { _owner = owner; }

	// Counts Reads in given Regions through the chromosome
	void Scan(const Regions& rgns);
	
	// Gets number of filled windows
	//inline chrlen WinsCount() const	{ return _wi+1; }

	// Gets average count of Reads for one widow in this distribution
	float WinDensity() const	{ return (float)_rtotalCnt / Length(); }	// WinsCount(); }
	
	// Gets total number of Reads in distribution
	inline ULONG ReadTotalCount() const { return _rtotalCnt; }
	
	// Adds Read distribution beginning from given position
	//	@rDistr: added distribution
	//	@start: position from which distribution will be added
	void ConcatDistrib(const ReadDistrib &rDistr, chrlen start)	{
		Concat(rDistr, start);
		_wi += rDistr._wi;
		//_rtotalCnt += rDistr._rtotalCnt;
	}

#ifdef DEBUG
	void Print(ostream& stream=cout) const;
#endif
};

class PairReadDistrib
/*
 * Class 'PairReadDistrib' represents two distribution of Reads: In- and Out of peaks,
 * for single chromosome
 */
{
private:
	chrlen	_wLen;	// user-defined length of window; duplicated from GenomeReadDistrib for the efficiency
	// members - constants in this instance 
	readlen	_halfRLen;			// half-length of Read
	chrlen	_rCnt;				// number of Reads in BedR for current chromosome
	BedR::cItemsIter _firstRit;	// iterator to first Read in BedR for current chromosome
	Regions	_regns[2];			// Regions in which it is required to calculate the density
	// For In-peaks distribution it has been used directly,
	// for Out-peak it has been inverted and ovarlaped with defined Regions.
	ReadDistrib	_rDistribs[2];	// Read density distributions

public:
	// sign of In-|Out- features location
	enum eLocating {
		IN_P = 0,	// in regions (peaks)
		OUT_P = 1	// out of regions (peaks)
	};

	// default constructor needed for base Chroms collection
	inline PairReadDistrib() : _rCnt(0) {}

#ifdef _NO_UNODMAP
	// copy constructor used in Chroms::AddClass() and by automatic reallocation
	PairReadDistrib(const PairReadDistrib& prd) {
		_halfRLen = prd._halfRLen;
		_rCnt = prd._rCnt;
		_firstRit = prd._firstRit;
		// these initializations are needed when moving basic vector by automatic reallocation
		_regns[IN_P] = prd._regns[IN_P];
		_regns[OUT_P] = prd._regns[OUT_P];
		_rDistribs[IN_P].SetOwner(this);
		_rDistribs[OUT_P].SetOwner(this);
	}
#endif

	// Creates new instance by chrom, BedR, and regions
	//	@cID: chromosome's ID
	//	@bedR: BedR for what distribution should be builded
	PairReadDistrib(chrid cID, const BedR &bedR) :
		_halfRLen( bedR.ReadLen()/2 ),
		_firstRit( bedR.ReadsBegin(cID) ),
		_rCnt( bedR.ReadsCount(cID) )
		{}

	// Initializes regions
	//	@cID: chromosome's ID
	//	@cRegns: defined reference chrom Regions
	//	For In-peak distribution has been used directly,
	//	for Out-peak has been inverted and ovarlaped with chrom defined Regions.
	//	@bedF: optional additional regions
	void InitRegions(chrid cID, const Regions &cRegns, const BedF *bedF);

	// Creates distributions for givem @location
	void Scan(eLocating loc, chrlen winLen)
	{
		_wLen = winLen;
		_rDistribs[loc].Scan(_regns[loc]);
	}

	// Gets number of filled windows
	inline chrlen WinsCount(eLocating loc) const { return _rDistribs[loc].Length(); }

	// Gets density, averaged for whole genome
	inline float Density(eLocating loc) const { return _rDistribs[loc].WinDensity() / _wLen; }

	// Gets the centre of given Read for current chromosome
	// Used once in ReadDistrib.ScanWindow(..) only.
	// Public because of denied 'friend ReadDistrib' for g++ 4.1.2
	inline chrlen ReadCentrePos(chrlen rRelInd) const {
		return (_firstRit + rRelInd)->Pos + _halfRLen;
	}

	// Gets number of reads in BedR for current chromosome
	// Used once in ReadDistrib::ScanWindow(..) only.
	// Public because of denied 'friend ReadDistrib' for g++ 4.1.2
	inline chrlen ChromReadsCount() const	{ return _rCnt; }

	// Gets length of cutting windows
	inline chrlen WinLength() const	{ 
		return _wLen; 
	}

	// Gets total number of Reads in distribution
	inline ULONG ReadsCount(eLocating loc) const	{ 
		return _rDistribs[loc].ReadTotalCount();
	}

	// Gets Read's distribution after scanning.
	//	@loc: In-|Out- location
	inline  ReadDistrib & Distrib(eLocating loc) {
		return _rDistribs[loc];
	}

	float GetResults(eLocating loc, ULONG* rCnt, ULONG* wCnt) const;

	//friend ReadDistrib;	// g++ 4.1.2 doesn't allow to declare 'friend' for class that used predifined name of this class
};

class GenomeReadDistrib : public Chroms<PairReadDistrib>
/*
 * Class 'GenomeReadDistrib' generates and keeps a collection of Read's distributions 
 * for each chromosome.
 */
{
private:
	ReadDistrib _totalReadDistr;	// ChromReadDistr for total gemone
	bool	_twoDistrs;				// true if both distributiona are defined
	chrlen	_wLen;					// user-defined length of window

	// Gets total density
	inline float GetTotalDensity(BYTE loc, ULONG* const rCnt, ULONG* const wCnt) const {
		return (double(rCnt[loc]) / wCnt[loc] / _wLen) * DENS_BASE;
	}

public:
	// Initializes In-|Out-Read's distribution container:
	// for each chromosome initializes default Regions and PairReadDistrib.
	// Only chromosomes marked as 'Treated' would be treated.
	//	@bedR: Reads wich are distribeted
	//	@gRegn: define genome Regions
	//	@bedF: features determining peaks or NULL
	GenomeReadDistrib(const BedR& bedR, GenomeRegions& gRegn, const BedF* bedF);
	
	// Initializes and fills In-|Out-Read's distribution for each chromosome in genome
	//	@winLen: length of layout (cutting) window
	GenomeReadDistrib& Scan	(chrlen winLen);
	
	// Returns density, averaged for whole genome
	//	@loc: In-|Out- location
	float	Density (PairReadDistrib::eLocating loc) const;

	// Returns total number of Reads in n-|Out- distribution
	//	@loc: In-|Out- location
	ULONG ReadsCount	(PairReadDistrib::eLocating loc) const;

	// Returns Read's distribution after scanning.
	//	@loc: In-|Out- location
	ReadDistrib& Distrib(PairReadDistrib::eLocating loc);

	inline const PairReadDistrib & operator[] (chrid cID) const { return At(cID); }

	void PrintDensity();
};
