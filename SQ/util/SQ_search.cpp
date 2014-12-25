#include <algorithm>
#include <mex.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>

using namespace std;

/* Asymmetric Distance computation (ADC), for stacked quantizers. */
 
/* Input Arguments */
#define	CODES       prhs[0]  //
#define CODEBOOKS   prhs[1]  //
#define QUERIES     prhs[2]  //
#define DBNORMS     prhs[3]  //
#define NNEIGHBOURS prhs[4]  //

/* Output Arguments */
#define	DISTS	plhs[0]  // The (sorted) products from each query to the database points.
#define	IDX 	plhs[1]  // The indices in the database of the distances above.

/* Define C types to actually match matlab types. */
typedef unsigned long  UINT64;
typedef unsigned int   UINT32;
typedef unsigned short UINT16;
typedef unsigned char  UINT8;
typedef double         REAL;

// MATLAB ENTRY POINT
void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray*prhs[] )

{
    // ===== Check for proper number of arguments. =====
    if (nrhs != 5)
    mexErrMsgTxt("5 input arguments required.");

    if (nlhs != 2)
    mexErrMsgTxt("2 output arguments required.");

    // ===== Check input types. =====
    if (!mxIsUint8( CODES ))
	mexErrMsgTxt("CODES must be uint8");

    // Codebooks have to be checked individually.

    if (!mxIsSingle( QUERIES ))
	mexErrMsgTxt("QUERIES must be single");

    if (!mxIsSingle( DBNORMS ))
    mexErrMsgTxt("DBNORMS must be single");

    // ===== Check input dimensions =====
    long nlevels = (long) mxGetM( CODES );
    long ncodes  = (long) mxGetN( CODES );

    // DBNORMS must have ncodes entries.
    long mdbnorms = (long) mxGetM( DBNORMS );
    long ndbnorms = (long) mxGetN( DBNORMS );

    if (mdbnorms == 1 && ndbnorms == 1 && ncodes != 1)
    mexErrMsgTxt("Number of DBNORMS seems too small.");

    long normspassed = mdbnorms > ndbnorms ? mdbnorms : ndbnorms;

    if ( normspassed != ncodes )
    mexErrMsgTxt("Number of DBNORMS not equal to number of codes.");

    long d        = mxGetM( QUERIES );
    long nqueries = mxGetN( QUERIES );

    long mcodebooks = (long) mxGetM( CODEBOOKS );
    long ncodebooks = (long) mxGetN( CODEBOOKS );

    long passedcodebooks = mcodebooks > ncodebooks ? mcodebooks : ncodebooks;

    // Assuming that the codebooks have all the same number of elements.
    long codebooksz = (long) mxGetN( mxGetCell( CODEBOOKS, 1 ));

    // ===== Get pointers to the input. =====
    unsigned char  *codesp     = (unsigned char*)  mxGetPr( CODES );
    float *queriesp   = (float*) mxGetPr( QUERIES );
    float *dbnormsp   = (float*) mxGetPr( DBNORMS );
    int   nneighbours = (int) mxGetScalar( NNEIGHBOURS );

    // ===== Create output arrays. =====
    DISTS = mxCreateNumericMatrix(nneighbours, nqueries, mxSINGLE_CLASS, mxREAL);
    IDX   = mxCreateNumericMatrix(nneighbours, nqueries, mxINT32_CLASS, mxREAL);

    // ===== Get pointers to the output. =====
    float *distsp = (float*) mxGetPr( DISTS );
    int   *idxp   = (int*)   mxGetPr( IDX );

    // ===== Convert the codebook cell array to a big matrix. =====
    float *pcodebook = new float[ codebooksz * passedcodebooks * d ]();

    for (int i=0; i<passedcodebooks; i++) {

        // Check codebooks type.
        if (!mxIsSingle( mxGetCell( CODEBOOKS, i ) ))
        mexErrMsgTxt("CODEBOOKS must be single");

        float *ithcodebook = (float *) mxGetPr( mxGetCell( CODEBOOKS, i ));

        for (int j=0; j < codebooksz; j++) { // Loop over entries.
            for (int k=0; k < d; k++) {      // Loop over dimensions.
                pcodebook[ i*codebooksz*d + j*d + k ] = ithcodebook[ j*d + k ];
            }
        }
    }

    // ===== Compute the distance from each query to the database =====
    int totalcentries = codebooksz*passedcodebooks;

    #pragma omp parallel for
    for( int i=0; i<nqueries; i++ ) {                           // Loop over queries.
        float* query  = queriesp + i*d;      // Create a pointer for this query.
        float* tentry = new float[ totalcentries ]();

        // ===== Create table entries =====
        for( int j=0; j<codebooksz * passedcodebooks; j++ ) {   // Loop over codebook entries.
            float* centry = pcodebook + j*d; // Create a pointer for this codebook entry.

            for ( int k=0; k<d; k++) {                          // Loop over dimensions.
                tentry[ j ] += query[k]*centry[k]; // Compute the dot product.
            }
        }

        pair<float,int> * pairs = new pair<float,int>[ ncodes ](); // We will store the distances here.
		unsigned char *code;

        // ===== Distance computation =====
        for (long j=0; j<ncodes; j++) {      // Loop through database vectors.
            pairs[j].second = j+1; // Save the index.

			code = codesp + j*nlevels;   // Pointer to this code.

            for (int k=0; k<nlevels; k++) {  // Loop through vector codes.
                // int code = codesp[ j*nlevels + k ];

                // IN THE NEXT LINE, WE WOULD NORMALLY SUBSTRACT ONE, BUT
                // WE'VE TAKEN CARE OF THAT IN MATLAB.
                pairs[j].first -= tentry[ 256*k + code[k] ]; // Table lookup and sum.
            }
            pairs[j].first *= 2;
            pairs[j].first += dbnormsp[ j ]; // Add the database norm.
        }

        partial_sort(pairs, pairs + nneighbours, pairs + ncodes); // Sort the distances.

        for (long j=0; j<nneighbours; j++) {
            distsp[i*nneighbours + j] = pairs[j].first;
            idxp  [i*nneighbours + j] = pairs[j].second;
        }
        delete [] pairs;
        delete [] tentry;
    }

    delete [] pcodebook;

    return;
}
