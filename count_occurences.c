#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = " <history_vector> <seqfile>";
static char banner[] = "px, the first parallel tests of H4";

int
main(int argc, char **argv)
{

  /* This is some Easel convenience code for managing arguments to functions */
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 2, argc, argv, banner, usage);

  int status, i;
  int num_sequences = 0;
  ESL_DSQDATA *dd = NULL;

  /* History_vector and seqfile are the two inputs.  Seqfile is the sequence database (in dsq format).  History_vector specifies which previous amino acids you're conditioning on, using the binary format you've defined. */
  uint32_t        history_vector = atoi(esl_opt_GetArg(go, 1));
  char           *seqfile = esl_opt_GetArg(go, 2);
  
  /* Code to open the sequence database */
  ESL_ALPHABET *abc = NULL;
  ESL_DSQDATA_CHUNK *chu = NULL;
  /* Open sequence database */
  status = esl_dsqdata_Open(&abc, seqfile, 1, &dd);
  if      (status == eslENOTFOUND) esl_fail("Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   esl_fail("Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        esl_fail("Unexpected error in opening dsqdata \n %s", dd->errbuf );

/* Ok, here's where the real program starts. */
/* Step one: allocate an array of unsigned 64-bit integers that will hold the counts
   of amino patterns that you see */

   int num_past_residues = 0; // number of prior residues we look at
   int count_array_size = 20;  // size of pattern_counts array that holds counts
   int furthest_back = 0;     // how many residues back is the furthest we have to look
   int look_back_indeces[32]; // which of the prior residues do we look at
   int look_back_idx = 0;     // counter for loop below

/* Count the number of ones in history_vector, which is the number of resdues you're
  conditioning on, and compute the size of the count array as 20 ^ number of residues
  tracked */ 

   for(i = 0; i < 32; i++){
     if(history_vector & 1){ // low bit is set
      num_past_residues++; // increase number of residues we need to look at
      count_array_size = count_array_size * 20; // increase size of array
      furthest_back = i + 1; // update furthest lookback index
      look_back_indeces[look_back_idx] = i+1;
      look_back_idx++; 
    }
    history_vector = history_vector >> 1;
   }
  
  printf("furthest_back = %d \n", furthest_back);
  printf("look_back_idx = %d \n", look_back_idx);
  printf("look_back_indeces: \n");
  for(i = 0; i < look_back_idx; i++) {
      printf("%d",look_back_indeces[i]);
  }
  printf("\n");
  uint64_t *pattern_counts = NULL;  // This is the array that will hold the counts of 
  //patterns seen

  // Allocate memory for the array. Will signal an error if this fails
  ESL_ALLOC(pattern_counts, count_array_size*8);

  /* Ok, this is the main loop. The dsq file reader returns chunks of sequences in
  data structures that look like: 
typedef struct esl_dsqdata_chunk_s {
  int64_t   i0;           // Chunk contains sequences i0..i0+N-1 from the database, 0-offset
  int       N;            // Chunk contains N sequences

  ESL_DSQ **dsq;          // Pointers to each of the N sequences
  char    **name;         // Names, \0 terminated.  Ptr into <metadata> buffer.
  char    **acc;          // Optional accessions, \0 terminated;   "\0" if none. 
  char    **desc;         // Optional descriptions, \0 terminated; "\0" if none
  int32_t  *taxid;        // NCBI taxonomy identifiers. (>=1 is a taxid; -1 means none)
  int64_t  *L;            // Sequence lengths, in residues. The unpacker figures these out.

  So, there's an outer loop that iterates over the chunks in the file, and an  inner loop that iterates over the sequences within each chuk.  Right now, it just counts the
  number of sequences and residues in the file, but this is where your code will go.

  */

  int m; // index into pattern_counts

  uint64_t num_residues = 0;
  uint8_t *the_sequence; // Sequences are arrays of 8-bit unsigned integers,
  //one byte/residue

  // Iterate over chunks in file
  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK) {
     
      // Iterate over sequences in chunk
	for (i = 0; i < chu->N; i++) {
		the_sequence = chu->dsq[i]; // Get the sequence 
       
         	// Iterate through residues in sequence
	 	for(int j = furthest_back; j <= chu->L[i]; j++) {
			num_residues++;
          		// j increases by one until it gets to the end of the length of the chunk/sequence
			int m_last = the_sequence[j];
			m = 0;
                        if(m > 19){continue;} // if residue is not one of 20 amino acids, skip
			
			int bad_letter = 0; // keep track of if we find a non-amino acid letter

			// iterate over lookback indeces and update m
          		for(int k = look_back_idx - 1; k >= 0; k--) {
	     		        if(the_sequence[j-look_back_indeces[k]] > 19){
				    bad_letter = 1; // found non-amino acid letter
				    break;} // if residue is not one of 20 amino acids skip it
				// m is the idx into pattern counts for the specific lookback of this residue
				m = m*20 + the_sequence[j-look_back_indeces[k]];
			} 

			// add final residue to m
			m = m*20 + m_last;

			// increment pattern count for this pattern as long as no bad letters are in midst
        		if(bad_letter == 0){pattern_counts[m]++;}
		}
		num_sequences++;
	
	} 
    	esl_dsqdata_Recycle(dd, chu);
    }

    int modulo = 0; // keep track of conditionals
    uint64_t num_conditionals;
    double *pattern_probabs = NULL;
    ESL_ALLOC(pattern_probabs, count_array_size*8);

    //for printing purposes
    int meep0 = 29;
    int meep1 = 27;
    int meep2 = 28;

    // cycle through each pattern
    for(int i = 0; i < count_array_size; i++) {

	num_conditionals = 0;	
	modulo = i - (i % 20);
	
	for (int j = 0; j < 20; j++) {
	    num_conditionals += pattern_counts[modulo + j];
	}
	
	pattern_probabs[i] = (float) pattern_counts[i] / (float) num_conditionals;

	// for printing purposes
	if(i == 0){meep0 = num_conditionals;}
	if(i == 1){meep1 = num_conditionals;}
	if(i == 2){meep2 = num_conditionals;}
    }

    printf("%f %llu %d \n", pattern_probabs[0], pattern_counts[0], meep0);
    printf("%lf %llu %d \n", pattern_probabs[1], pattern_counts[1], meep1);
    printf("%f %llu %d \n", pattern_probabs[2], pattern_counts[2], meep2);

//    printf("%llu %llu %llu %llu \n", pattern_counts[0], pattern_counts[1], pattern_counts[2], pattern_counts[3]); 
//    printf("%llu %llu %llu %llu \n", pattern_counts[4], pattern_counts[5], pattern_counts[6], pattern_counts[7]);
//    printf("%llu %llu %llu %llu \n", pattern_counts[8], pattern_counts[9], pattern_counts[10], pattern_counts[11]);
//    printf("%llu %llu %llu %llu \n", pattern_counts[12], pattern_counts[13], pattern_counts[14], pattern_counts[15]);
//    printf("%llu %llu %llu %llu \n", pattern_counts[16], pattern_counts[17], pattern_counts[18], pattern_counts[19]);
//    printf("%llu %llu %llu %llu ", pattern_counts[20], pattern_counts[21], pattern_counts[22], pattern_counts[23]);
//    printf("%llu %llu %llu %llu ", pattern_counts[24], pattern_counts[25], pattern_counts[26], pattern_counts[27]);
//    printf("%llu \n", pattern_counts[28]);
    printf("count_array_size = %d \n", count_array_size);
    printf("Conditioning on %d residues\n", num_past_residues);
    printf("Saw %d sequences, %llu residues\n", num_sequences,num_residues);
  esl_dsqdata_Close(dd);
  esl_getopts_Destroy(go);
  exit(0);

  ERROR:
  esl_fail("Unable to allocate memory\n", NULL);
}

