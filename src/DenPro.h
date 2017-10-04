#pragma once
#include "Distrib.h"

enum optValue {
	oGFILE,
	oGAPLEN,
	oDUPL,
	oDIFFSZ,
	oCHROM,
	oMINSCR,
	oCONS,
	oFBED,
	oEXTLEN,
	oSPACE,
	oINFO,
	oALARM,
	oOUTFILE,
	oTIME,
	oVERSION,
	oHELP
};

class DenPro
/*
 * 'DenPro' represents density profile
 */
{
private:
	struct DensDistr {
		chrlen	CountWins;
		chrlen	CountReads;

		inline DensDistr(chrlen cntReads) : CountWins(1),  CountReads(cntReads) {};
		inline DensDistr(chrlen cntWins, chrlen cntReads) : CountWins(cntWins), CountReads(cntReads) {};
		// for ascending sorting by default
		inline bool operator < (const DensDistr& denDistr) const {	return (CountReads < denDistr.CountReads); }
	};

	vector<DensDistr>	_densDistr[2];		// precise density distribution
	vector<DensDistr>	_consDensDistr[2];	// consolidated density distribution

	// Adds value to the distribution: increase counter if value exists or add new counter.
	//	@densDistr: adding distribution
	//	@cntReads: added value
	void AddVal(BYTE loc, chrlen cntReads);

public:
	static chrlen	WinLen;		// user-defined length of sliding window

	DenPro(BedR &bed, GenomeRegions &gRegn, BedF *fbed);
	
	// Initializes precise Read's distribution
	//	@grDistr: genome reads distribution
	//	@loc: In-|Out- features location
	//	return: count of distributed Reads
	ULONG Init(GenomeReadDistrib& grDists, PairReadDistrib::eLocating loc);

	// Split sets of pairs <ReadCnt>-<WinCnt> into subset with defined length of ReadCnt 
	//	@consLen: step of number of consolidated reads or vUNDEF if no consolidation
	void Consolidate(int consLen);

	// Outputs to stream Read's distribution
	//	@stream: outputted stream. Currently not used.
	//	@printTwoDistribs: true if template is setting
	void Print(ofstream& stream);	//, bool printTwoDistribs);
	
	// Writes output to output file
	//void Write();
	//void Write(const char * bedFfName, chrlen fCnt, const char * bedRfName, chrlen rCnt);
};


/*
class Test
{
private:
	Array<int> arr;
public:
	inline Test() {
		arr.Reserve(3);
		arr[0] = 10; arr[1] = 20; arr[2] = 30;
	}
	inline void Print() { 
		for(int i=0; i<arr.Length(); i++) 
			cout << i << ": " << arr[i] << EOL;
	}
	inline void Reinit() {
		//arr.~Array();
		arr = Array<int>(5);
		for(int i=0; i<arr.Length(); i++)
			arr[i] = i+1;
	}
};
*/
