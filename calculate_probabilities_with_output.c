#include <string.h>
#include <stdio.h>
#include <stdlib.h>
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

  int status, i, j, k;
  int num_sequences = 0;
  ESL_DSQDATA *dd = NULL;

  /* History_vector and seqfile are the two inputs.  Seqfile is the sequence database (in dsq format).
  History_vector specifies which previous amino acids you're conditioning on, using the binary format you've defined. */
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

  int num_past_residues = 0;  // number of prior residues we look at
  int count_array_size = 20;  // size of pattern_counts array that holds counts
  int furthest_back = 0;      // how many residues back is the furthest we have to look
  int look_back_indeces[32];  // which of the prior residues do we look at
  int look_back_idx = 0;      // counter for loop below

  /* Count the number of ones in history_vector, which is the number of resdues you're
  conditioning on, and compute the size of the count array as 20 ^ number of residues
  tracked */

  // get output file name for conditional probabilities output file
  char output_file_name[256];
  char training_set_name[256];
  char folder[] = "cond_probab_tables/";
  char *arg = strrchr(seqfile, '/');
  if (arg != NULL) {
      arg++;
      sprintf(output_file_name, "%scond_probabs_%d_%s.txt", folder, history_vector, arg);
      sprintf(training_set_name, "%s", arg);
  }
  else {
      sprintf(output_file_name, "%scond_probabs_%d_%s.txt", folder, history_vector, seqfile);
      sprintf(training_set_name, "%s", seqfile);
  }
  printf("output file name: %s\n", output_file_name);

  // get lookback index into formats needed
  for(i = 0; i < 32; i++){
    if(history_vector & 1){                     // low bit is set
      num_past_residues++;                      // increase number of residues we need to look at
      count_array_size = count_array_size * 20; // increase size of array
      furthest_back = i + 1;                    // update furthest lookback index
      look_back_indeces[look_back_idx] = i+1;   // add lookback to index
      look_back_idx++;                          // update number of lookbacks
    }
    history_vector = history_vector >> 1;
   }

  // open file to write the conditional probabilities
  FILE *fptr;
  fptr = fopen(output_file_name,"w");
  // fprintf(fptr,"training dataset: %s\n", training_set_name);
  // fprintf(fptr,"look back idx: ");

  // // ----------------------------------------------
  // // open file to write pattern counts for checking
  // FILE *pcountsfile;
  // pcountsfile = fopen("../databases/count_occurs_rando.txt","w");
  // fprintf(pcountsfile,"training dataset: swissprot\n");
  // fprintf(pcountsfile,"look back idx: ");
  // for(i = 0; i < look_back_idx; i++) {
  //     fprintf(pcountsfile,"%d",look_back_indeces[i]);
  // }
  // fprintf(pcountsfile,"\n");
  // // ----------------------------------------------

  // add lookback index to file to keep track
  printf("furthest_back = %d \n", furthest_back);
  printf("look_back_idx = %d \n", look_back_idx);
  printf("look_back_indeces: \n");
  for(i = 0; i < look_back_idx; i++) {
      printf("%d",look_back_indeces[i]);
      // fprintf(fptr,"%d",look_back_indeces[i]);
  }
  if(look_back_idx == 0){printf("None");}
  // fprintf(fptr,"\n");
  printf("\n");

  // this is the array that will hold the counts of patterns seen
  uint64_t *pattern_counts = NULL;

  // this array holds marginal counts
  uint64_t *pattern_counts_marginal = NULL;

  // Allocate memory for the array. Will signal an error if this fails
  ESL_ALLOC(pattern_counts, count_array_size*8);
  ESL_ALLOC(pattern_counts_marginal, 20*8);

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

  // int m; // index into pattern_counts
  // int bad_letter; // bool to track if idx includes bad letters

  uint64_t num_residues = 0;
  uint8_t *the_sequence; // sequences are arrays of 8-bit unsigned integers, one byte/residue

  int num_bad_letters = 0; // keep track of how many non-amino-acid residues there are in dataset

  // Iterate over chunks in file
  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK) {

      // Iterate over sequences in chunk
    	for (i = 0; i < chu->N; i++) {
      		the_sequence = chu->dsq[i]; // Get the sequence

          // keep track of residues at the beginning
      		for(j = 1; j < furthest_back+1 && j <= chu->L[i]; j++){
              int m_last = the_sequence[j];
              if(m_last > 20 && m_last < 27){ num_bad_letters++; } // if it's a bad letter increase the bad letter count
              if(m_last != 20 && m_last < 27){ num_residues++; } // if its a residue, increase number of residues seen
              if(m_last < 20){ pattern_counts_marginal[m_last]++; } // if it's an amino acid, add to marginal counts
      		}

          // Iterate through residues in sequence
      	 	for(j = furthest_back + 1; j <= chu->L[i]; j++) {
        			int m_last = the_sequence[j]; // get residue of interest
              if(m_last > 20 && m_last < 27){ num_bad_letters++; } // count number of non-amino acids
        			if(m_last != 20 && m_last < 27){ num_residues++; } // if residue increase count
              if(m_last < 20){ pattern_counts_marginal[m_last]++; } // if amino acid add to marginals

        			if(m_last > 19){continue;} // if residue is not one of 20 amino acids, skip

              int m = 0; // index into pattern_counts
              int bad_letter = 0; // bool to track if idx includes bad letters

        			// iterate over lookback indeces and update m
              for(k = look_back_idx - 1; k >= 0; k--) {

                  // if residue is not one of 20 amino acids skip it
      	     		  if(the_sequence[j-look_back_indeces[k]] > 19){
          				    bad_letter = 1; // found non-amino acid letter
          				    break;}

          				// m is the idx into pattern counts for the specific lookback of this residue
          				m = m*20 + the_sequence[j-look_back_indeces[k]];
        			}

        			// add final residue to m
        			m = m*20 + m_last;

        			// increment pattern count for this pattern as long as no bad letters are in midst
              if(bad_letter == 0){
                pattern_counts[m]++;
              }
      		}
    		  num_sequences++; // count sequences
      }
    	esl_dsqdata_Recycle(dd, chu);
  }

  /* Start conditional probability calculation */

  // keep track of conditionals for denominator calculation
  int modulo = 0;
  uint64_t num_conditionals;
  double *pattern_probabs = NULL;
  ESL_ALLOC(pattern_probabs, count_array_size*8);

  //for printing purposes
  int denom_len = 10;
  int denominators[10];

  // cycle through each pattern
  for(i = 0; i < count_array_size; i++) {

    	num_conditionals = 0;
    	modulo = i - (i % 20);

      // calculate denominator
    	for (j = 0; j < 20; j++) {
    	    num_conditionals += pattern_counts[modulo + j];
    	}

      // include Laplace smoothing and caluclate conditional probability
    	pattern_probabs[i] = ((float) pattern_counts[i] + 1.0) / ((float) num_conditionals + 20.0);
    	fprintf(fptr,"%f\n",pattern_probabs[i]);

    	// for printing purposes
    	if(i < denom_len){denominators[i] = num_conditionals;}
  }

  // add marginal probabilities to end of file and laplace smooth for consistency
  for(i = 0; i < 20; i++){
      fprintf(fptr,"%f\n",((float) pattern_counts_marginal[i] + 1.0) / ((float) num_residues - num_bad_letters + 20.0));
  }

  // // --------------------------------------------------
  // // output pattern counts to file for checking
  // for(i = 0; i < count_array_size; i++) {
  //     fprintf(pcountsfile,"%llu\n", pattern_counts[i]);
  // }
  // fclose(pcountsfile);
  // // --------------------------------------------------

  // print first few probabilities for spot checking
  printf("\nhead probabilities\n");
  for(k = 0; k < denom_len; k++){
      printf("%d: %f %llu %d \n", k, pattern_probabs[k], pattern_counts[k], denominators[k]);
  }

  // print rest of stats
  printf("\ncount_array_size = %d \n", count_array_size);
  printf("Conditioning on %d residues\n", num_past_residues);
  printf("Saw %d sequences, %llu amino acid residues\n", num_sequences, num_residues);
  printf("Number of non-amino-acid residues in dataset: %d\n", num_bad_letters);

  // close printing file
  fclose(fptr);

  esl_dsqdata_Close(dd);
  esl_getopts_Destroy(go);
  exit(0);

  ERROR:
  esl_fail("Unable to allocate memory\n", NULL);
}
