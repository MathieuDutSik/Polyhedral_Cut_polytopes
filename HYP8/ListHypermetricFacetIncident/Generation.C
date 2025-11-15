#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <math.h>

using namespace std;


// We have some function for iterating over all elements of 
// the symmetric group. This is useful for combinatorial
// enumeration
//
// very good timings
//  1s for iterating over Sym(9)
//  5s for iterating over Sym(10)
// 52s for iterating over Sym(11)


void PrintVector(vector<int> & TheVect)
{
  int i;
  cout << "(";
  for (i=1; i<=TheVect.size(); i++)
    {
      cout << TheVect[i-1];
      if (i < TheVect.size())
	{
	  cout << ",";
	}
    }
  cout << ")\n";
}


void FirstElementSymmetricGroup(int& n, vector<int> & TheElt, vector<vector<int> > & ListListSet)
{
  int i, j;
  vector<int> TheList;
  for (i=1; i<=n; i++)
    {
      TheElt.push_back(i);
    }
  for (i=1; i<=n; i++)
    {
      for (j=i; j<=n; j++)
	{
	  TheList.push_back(j);
	}
      ListListSet.push_back(TheList);
      TheList.clear();
    }
}






void NextElementSymmetricGroup(int& n, vector<int> & TheElt, vector<vector<int> > &ListListSet, int & result)
{
  int i, iCol, jCol, MatchedPos, MatchedCol, pos, siz;
  set<int> ListMatched;
  vector<int> Residue, TheSet, TheList;
  MatchedCol=-1;
  for (i=n; i>1; i--)
    {
      if (TheElt[i-1] > TheElt[i-2])
	{
	  MatchedCol=i-1;
	  siz=n+1-MatchedCol;
	  break;
	}
    }
  if (MatchedCol == -1)
    {
      result=0;
      return;
    }
  TheSet=ListListSet[MatchedCol-1];
  MatchedPos=-4;
  for (pos=1; pos<=siz; pos++)
    {
      if (TheSet[pos-1] == TheElt[MatchedCol-1])
	{
	  MatchedPos=pos;
	}
      if (MatchedPos+1 != pos)
	{
	  Residue.push_back(TheSet[pos-1]);
	}
    }
  TheElt[MatchedCol-1]=TheSet[MatchedPos];
  for (iCol=MatchedCol+1; iCol<=n; iCol++)
    {
      TheElt[iCol-1]=Residue[iCol-MatchedCol-1];
    }
  for (iCol=MatchedCol+1; iCol<=n; iCol++)
    {
      for (jCol=iCol; jCol<=n; jCol++)
	{
	  TheList.push_back(Residue[jCol-MatchedCol-1]);
	}
      ListListSet[iCol-1]=TheList;
      TheList.clear();
    }
  result=1;
}



void ReadMetaData(int& n,
		  int& nbOrbHyp,
		  int& nbCuts,
		  FILE *META)
{
  fscanf(META, "%ld", &n);
  fprintf(stderr, "Dimension = %ld\n", n);
  fscanf(META, "%ld", &nbOrbHyp);
  fprintf(stderr, "|Orbit hypermetric facet|= %ld\n", nbOrbHyp);
  fscanf(META, "%ld", &nbCuts);
  fprintf(stderr, "|Incident cuts| = %ld\n", nbCuts);
}


void ReadListINTS(vector<vector<int> >& ListHypermet, int& n, int& nbOrbHyp, FILE *f)
{
  int iOrb, iCol, eVal;
  vector<int> eList;
  for (iOrb=1; iOrb<=nbOrbHyp; iOrb++)
    {
      for (iCol=1; iCol<=n; iCol++)
	{
	  fscanf(f, "%ld", &eVal);
	  eList.push_back(eVal);
	}
      ListHypermet.push_back(eList);
      eList.clear();
    }
}

void PermutedVector(vector<int> & TheElt, vector<int> & TheVect, 
		    vector<int> & TheImage)
{
  int iCol;
  TheImage.clear();
  for (iCol=1; iCol<=TheVect.size(); iCol++)
    {
      TheImage.push_back(TheVect[TheElt[iCol-1]-1]);
    }
}



int IsIncident(int & nbCuts, int & n, 
	       vector<vector<int> > & ListCut,
	       vector<int> & TheHyp)
{
  int iCut, TheSum, iCol;
  vector<int> eCut;
  for (iCut=1; iCut<=nbCuts; iCut++)
    {
      eCut=ListCut[iCut-1];
      TheSum=0;
      for (iCol=1; iCol<=n; iCol++)
	{
	  if (eCut[iCol-1] == 1)
	    {
	      TheSum=TheSum+TheHyp[iCol-1];
	    }
	}
      if (TheSum != 0 && TheSum != 1)
	{
	  return 0;
	}
    }
  return 1;
}




void ComputeIncidence(int& n, int& nbOrbHyp, int& nbCuts,
		      vector<set<vector<int> > > & ListIncident, 
		      vector<vector<int> > &ListCut, 
		      vector<vector<int> > &ListHypermet)
{
  int result, iOrb, test;
  vector<int> TheImage, TheHyp;
  vector<int> TheElt;
  vector<vector<int> > ListListSet;
  set<vector<int> > TheEmpty;
  TheEmpty.clear();
  for (iOrb=1; iOrb<=nbOrbHyp; iOrb++)
    {
      ListIncident.push_back(TheEmpty);
    }
  FirstElementSymmetricGroup(n, TheElt, ListListSet);
  while(1)
    {
      for (iOrb=1; iOrb<=nbOrbHyp; iOrb++)
	{
	  TheHyp=ListHypermet[iOrb-1];
	  PermutedVector(TheElt, TheHyp, TheImage);
	  test=IsIncident(nbCuts, n, ListCut, TheImage);
	  if (test == 1)
	    {
	      ListIncident[iOrb-1].insert(TheImage);
	    }
	}
      NextElementSymmetricGroup(n, TheElt, ListListSet, result);
      if (result == 0)
	{
	  break;
	}
    }
}

void PrintListIncident(int& n, int& nbOrbHyp,
		       vector<set<vector<int> > > & ListIncident, 
		       FILE *OUTPUT)
{
  int IsFirstPrint, iOrb, iCol, eVal;
  set<vector<int> >::iterator iter;
  vector<int> eHyp;

  IsFirstPrint=1;
  fprintf(OUTPUT, "return [");
  for (iOrb=1; iOrb<=nbOrbHyp; iOrb++)
    {
      iter=ListIncident[iOrb-1].begin();
      while(iter != ListIncident[iOrb-1].end())
	{
	  eHyp=*iter;
	  if (IsFirstPrint == 0)
	    {
	      fprintf(OUTPUT, ",\n");
	    }
	  IsFirstPrint=0;
	  fprintf(OUTPUT, "[");
	  for (iCol=1; iCol<=n; iCol++)
	    {
	      if (iCol > 1)
		{
		  fprintf(OUTPUT, ",");
		}
	      eVal=eHyp[iCol-1];
	      fprintf(OUTPUT, "%ld", eVal);
	    }
	  fprintf(OUTPUT, "]");
	  iter++;
	}
    }
  fprintf(OUTPUT, "];\n");
}






int main(int argc, char *argv[])
{
  FILE *METADATA=NULL;
  FILE *HYPORBS=NULL;
  FILE *LISTCUT=NULL;
  FILE *OUTPUT=NULL;
  int n, nbOrbHyp, nbCuts;
  vector<vector<int> > ListHypermet;
  vector<vector<int> > ListCut;
  vector<set<vector<int> > > ListIncident;
  if (argc !=5)
    {
      fprintf(stderr, "Number of argument is = %ld\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "Generation [METADATA] [HYPORBS] [LISTCUT] [OUTPUT]\n");
      fprintf(stderr, "METADATA: n |OrbHyps| |LISTCUTS|\n");
      fprintf(stderr, "HYPORBS: the list of hypermetric vectors\n");
      fprintf(stderr, "LISTCUT: The List of selecting cuts\n");
      fprintf(stderr, "OUTPUT: list of incident hypermetrics\n");
      return -1;
    }
  METADATA=fopen(argv[1], "r");
  if (METADATA == NULL)
    {
      fprintf(stderr,"The file %s not found\n",argv[2]);
      return -1;
    }
  ReadMetaData(n, nbOrbHyp, nbCuts, METADATA);
  //
  HYPORBS=fopen(argv[2], "r");
  if (HYPORBS == NULL)
    {
      fprintf(stderr,"The file %s not found\n",argv[3]);
      goto _END;
    }
  ReadListINTS(ListHypermet, n, nbOrbHyp, HYPORBS);
  //
  LISTCUT=fopen(argv[3], "r");
  if (LISTCUT == NULL)
    {
      fprintf(stderr,"The file %s not found\n",argv[4]);
      goto _END;
    }
  ReadListINTS(ListCut, n, nbCuts, LISTCUT);
  //
  ComputeIncidence(n, nbOrbHyp, nbCuts,
		   ListIncident, ListCut, ListHypermet);
  //
  OUTPUT=fopen(argv[4], "w");
  PrintListIncident(n, nbOrbHyp, ListIncident, OUTPUT);
  fclose(OUTPUT);

 _END: ;
  return 0;
}
