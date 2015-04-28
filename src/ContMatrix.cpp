#include "ContMatrix.h"

namespace CM
{

//////////////////////
// CopyTH2DtoTMatrix()
void CopyTH2DtoTMatrix( TH2D *th2d_matrix, TMatrix *&tmatrix, int indices )
{

	int nSize = GetNdigits( indices );

	int rowcolno[nSize];
	FillIndices( rowcolno, indices );

	tmatrix = new TMatrix(nSize, nSize);
	
	for(int i=0; i<nSize; i++)
	for(int j=0; j<nSize; j++)
	{
		(*tmatrix)(i,j) = th2d_matrix->GetBinContent( rowcolno[i], rowcolno[j] );
	}
}

///////////////////
// CopyTH2DtoTH2D()
void CopyTH2DtoTH2D( TH2D *matrix_old, TH2D *&matrix_new, const char matrixname[], int indices )
{

	int nSize = GetNdigits( indices );

	int rowcolno[nSize];
	FillIndices( rowcolno, indices );

	matrix_new = new TH2D( matrixname, ";RECO;MC", nSize, 0.0, ((double)nSize), nSize, 0.0, ((double)nSize));
	
	for(int i=0; i < nSize; i++)
	for(int j=0; j < nSize; j++)
	{
		matrix_new -> SetBinContent( (i+1) , (j+1) , matrix_old->GetBinContent( rowcolno[i], rowcolno[j] ) );
	}
}

//////////////////////////
// normalizeColoumn_TH2D()
void normalizeColoumn_TH2D (TH2D *matrix, int colno )
{
	double entries = 0;

	int nSize = matrix->GetNbinsX();

	// To Be Finished
	for( int j = 1; j < nSize+1; j++ )
	{
		entries = entries + matrix->GetBinContent(colno,j);
	}

	for( int j = 1; j < nSize+1; j++ )
	{
		double r = matrix->GetBinContent(colno,j)/entries;
		matrix->SetBinContent(colno, j, r);
	}	
	
}

//////////////////////////
// normalizeColoumn_TH2D()
void normalizeMatrix_TH2D (TH2D *matrix )
{

	int nSize = matrix->GetNbinsX();

	for( int i = 1; i < nSize+1; i++ )
	{
		normalizeColoumn_TH2D(matrix, i);
	}	
	
}


/////////////////////
// displayMatrix_TH2D
void displayMatrix_TH2D( TH2D *matrix )
{
	int nSize = matrix->GetNbinsX();	

	std::cout << std::endl;
	for(int i = 1; i < (nSize+1); i++)
	{
		for(int j = 1; j < (nSize+1); j++)
		{ std::cout << Form("%2.2f", matrix->GetBinContent(i,j) ) << " "; }
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

////////////////////////
// displayMatrix_TMatrix
void displayMatrix_TMatrix( TMatrix *matrix )
{
	int nSize = matrix->GetNrows();	

	std::cout << std::endl;
	for(int i = 0; i < nSize; i++)
	{
		for(int j = 0; j < nSize; j++)
		{ std::cout << Form("%2.2f", (*matrix)(i,j) ) << " "; }
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

////////////////
// FillIndices()
void FillIndices( int *indices, int number )
{
	int nSize = GetNdigits(number);
	
	for(int i=0; i < nSize; i++)
	{ indices[i] = GetNthDigit(number, (nSize-1-i)); }
}

////////////////
// GetNthDigit()
int GetNthDigit(int number, int n)
{
	return ( ((int)(number/pow(10, n))) % 10);
}

///////////////
// GetNdigits()
int GetNdigits( int number )
{
	int digits = 0;
	 while (number) {
	            number /= 10;
	  	        digits++;
	  	     }
	return digits;	  
};

///////////////////
// plotContMatrix()
void plotContMatrix(TH2D *matrix, int indices, const char figbasename[])
{

	int nSize = matrix->GetNbinsX();
	std::cout << "plotContMatrix_nSize: " << nSize << std::endl;

	//
	
	TCanvas canvas_matrix ("contmatrix", ";RECO;GEN", canvas_res_x, canvas_res_y);

	gStyle->SetOptStat(0);
	gStyle->SetPaintTextFormat("2.2f");

	matrix->SetMarkerSize(2);
	matrix->Draw("TEXT");

 	matrix->GetXaxis()->SetLabelSize(0.05);
 	matrix->GetYaxis()->SetLabelSize(0.05);

	int indx[nSize];

	FillIndices( indx, indices );

	for (int i=1; i < (nSize+1); i++)
	{

		matrix->GetXaxis()->SetBinLabel( i, GetTypeLatex(indx[i-1]).c_str() );
		matrix->GetYaxis()->SetBinLabel( i, GetTypeLatex(indx[i-1]).c_str() );
	}

	std::string outPNG = Form("%s.png", figbasename);
	std::string outPDF = Form("%s.pdf", figbasename);

	canvas_matrix.SaveAs(outPDF.c_str());
	canvas_matrix.SaveAs(outPNG.c_str());

}

std::string GetTypeLatex (int ID)
{
	std::string out;

	switch (ID)
   { 
		case 1: out = "#pi"; return out;
		case 2: out = "K"; return out;
		case 3: out = "p"; return out;
		case 4: out = "non-p"; return out;
   }
}

}
