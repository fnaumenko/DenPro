#pragma once
#include "Data.h"

class PairReadDistrib;		// forward declaration used in ReadDistrib

class ReadDistrib : public Array<chrlen>
/*
 * 'ReadDistrib' represents simple distribution of Reads,
 * that's a number of Reads in each cutting window. Windows are indexed from 0 to maxWinCnt
 */
{
	// members are needed to keep their values during discontinuous scanning
	chrlen	_ri,			// current read's index in BedR
			_wCurrLen,		// length of uncsanning rest of current window
			_wi;			// current window's index
	// summerised members
	ULONG	_rtotalCnt;		// total number of reads in distribution
	// members - constants in this instance 
	Regions *_regns;		// Regions on which distribution is defined
	PairReadDistrib *_owner;// Container of this instance
	
	// Returns count of Reads in window
	//	@start: window's start position
	//	@outWin: returned true if current Read passed left boundary of window, otherwise false
	//	Return: count of Reads in window
	chrlen ScanWindow(chrlen start, bool* outWin);

	// Returns count of Reads in Region
	//	Return: count of Reads in Region
	void ScanRegion(const Region&);

public:
	inline ReadDistrib() { _rtotalCnt = _wi = _ri = _wCurrLen = 0; }

	//	@regns: Regions for wich distribution is defined
	// Prepares to scanning: reserves capacity and zero counters.
	//	@owner: owner (holder) of this distribution
	void Init(Regions* const regns, PairReadDistrib* const owner);

	// Counts Reads in defined Regions through the chromosome
	void Scan();
	
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
	chrlen	_wLen;		// user-defined length of window
	// members - constants in this instance 
	readlen	_halfRLen;	// half-length of Read
	chrlen	_rCnt;		// number of Reads in BedR for current chromosome
	BedR::cItemsIter _beginIt;	// iterator to first Read in BedR for current chromosome
	Regions	_templRegns[2];		// external additional template Regions at BedF or at chromosome.
	// For In-peaks distribution it has been used directly,
	// for Out-peak it has been inverted and ovarlaped with defined Regions.

public:
	ReadDistrib	_rDistribs[2];

	// sign of In-|Out- features location
	enum eLocating {
		IN_P = 0,	// in peaks
		OUT_P = 1	// out of peaks
	};

	// default constructor needed for base Chroms collection
	inline PairReadDistrib() {};

	// default constructor needed for base Chroms collection
	PairReadDistrib(const PairReadDistrib& d) {
		_wLen = d._wLen;
		_halfRLen = d._halfRLen;
		_rCnt = 0;
	};

	// Creates new instance by chrom, BedR, and regions
	//	@cID: chromosome's ID
	//	@bedR: BedR for what distribution should be builded
	//	@templRegns: external additional template Regions from BedF or whole chromosome.
	//	For In-peak distribution has been used directly,
	//	for Out-peak has been inverted and ovarlaped with chrom defined Regions.
	PairReadDistrib(chrid cID, const BedR &bedR, const Regions &defRegns, const BedF *bedF);

	// Prepares to scanning
	//	@loc: location
	//	@winLen: length of cutting window
	inline void Init(eLocating loc, chrlen winLen)
	{
		_wLen = winLen;
		_rDistribs[loc].Init(&_templRegns[loc], this);
	}

	// Creates distributions for givem @location
	inline void Scan(eLocating loc) { _rDistribs[loc].Scan(); }

	// Gets number of filled windows
	inline chrlen WinsCount(eLocating loc) const { return _rDistribs[loc].Length(); }

	// Gets density, averaged for whole genome
	inline float Density(eLocating loc) const { return _rDistribs[loc].WinDensity() / _wLen; }

	// Gets the centre of given Read for current chromosome
	// Used once in ReadDistrib.ScanWindow(..) only.
	// Public because of denied 'friend ReadDistrib' for g++ 4.1.2
	inline chrlen ReadCentrePos(chrlen rRelInd) const {
		return (_beginIt + rRelInd)->Pos + _halfRLen;
	}

	// Gets number of reads in BedR for current chromosome
	// Used once in ReadDistrib::ScanWindow(..) only.
	// Public because of denied 'friend ReadDistrib' for g++ 4.1.2
	inline chrlen ChromReadsCount() const	{ return _rCnt; }

	// Gets length of cutting windows
	inline chrlen WinLength() const	{ return _wLen; }

	// Gets total number of Reads in distribution
	inline ULONG ReadsCount(PairReadDistrib::eLocating loc) const	{ 
		return _rDistribs[loc].ReadTotalCount();
	}

	// Gets Read's distribution after scanning.
	//	@loc: In-|Out- location
	inline  ReadDistrib & Distrib(PairReadDistrib::eLocating loc) {
		return _rDistribs[loc];
	}

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
	ReadDistrib & Distrib(PairReadDistrib::eLocating loc);

	inline const PairReadDistrib & operator[] (chrid cID) const { return At(cID); }

	void PrintDensity();
};
