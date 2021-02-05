#include <string.h>
#include "easel.h"
#include "esl_dsqdata.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                               docgroup*/
  { "-h",        eslARG_NONE,  FALSE,  NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",  0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "<test_set_to_score> <conditional_probabilities_file> <history_vector>";
static char banner[] = "px, the first parallel tests of H4";

int
main(int argc, char **argv)
{

  /* easel convenience code for managing arguments to functions */
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 3, argc, argv, banner, usage);
  int status, i;
  int num_sequences = 0;
  ESL_DSQDATA *dd = NULL;

  /* history_vector is the binary lookback index, seqfile is the test dataset file, and trainingfile is the cond probabs */
  uint32_t        history_vector = atoi(esl_opt_GetArg(go, 3));
  char           *seqfile = esl_opt_GetArg(go, 1);
  char           *trainingfile = esl_opt_GetArg(go, 2);

  /* open the test sequence file */
  ESL_ALPHABET *abc = NULL;
  ESL_DSQDATA_CHUNK *chu = NULL;
  status = esl_dsqdata_Open(&abc, seqfile, 1, &dd);
  if      (status == eslENOTFOUND) esl_fail("Failed to open dsqdata files:\n  %s",    dd->errbuf);
  else if (status == eslEFORMAT)   esl_fail("Format problem in dsqdata files:\n  %s", dd->errbuf);
  else if (status != eslOK)        esl_fail("Unexpected error in opening dsqdata \n %s", dd->errbuf);

  /* get details on lookback index */
  int num_past_residues = 0; // number of prior residues we look at
  int count_array_size = 20;  // size of pattern_counts array that holds counts
  int furthest_back = 0;     // how many residues back is the furthest we have to look
  int look_back_indeces[32]; // which of the prior residues do we look at
  int look_back_idx = 0;     // counter for loop below, same thing as num_past_residues
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

  // /* print variables to terminal for checking */
  // printf("count_array_size  = %d\n", count_array_size);
  // printf("num_past_residues = %d\n", num_past_residues);
  // printf("furthest_back     = %d\n", furthest_back);
  // printf("look_back_indeces:  ");
  // for(i = 0; i < look_back_idx; i++) {
  //     printf("%d",look_back_indeces[i]);
  // }
  // printf("\n\n");
  // printf("seqfile: %s\n",seqfile);
  // printf("trainingfile: %s\n\n",trainingfile);

  /* open the conditional probabilities file */
  FILE *cond_probab_file = fopen(trainingfile, "r");
  if(cond_probab_file == NULL){
      fprintf(stderr,"Error opening file.\n");
      exit(2);
  }

  /* start array to hold conditional probabilities */
  double *pattern_probabs = NULL;
  ESL_ALLOC(pattern_probabs, count_array_size*8);
  double *marginal_probabs = NULL;
  ESL_ALLOC(marginal_probabs, count_array_size*8);

  /* read each conditional probability in the cond_probabs file */
  for(int i = 0; i < count_array_size; i++){
      fscanf(cond_probab_file,"%lf\n", &(pattern_probabs[i]));
  }

  // // printing check
  // for(int i = 0; i < count_array_size; i++){
  //     printf("%d: %lf\n",i,pattern_probabs[i]);
  // }

  /* read each marginal probability in the cond_probabs file */
  for(int i = 0; i < 20; i++){
      fscanf(cond_probab_file,"%lf\n", &(marginal_probabs[i]));
  }

  // // printing check
  // for(int i = 0; i < 20; i++){
  //     printf("%d: %lf\n",i,marginal_probabs[i]);
  // }

  // close file now that everything has been read from it
  fclose(cond_probab_file);

  /* The dsq file reader returns chunks of sequences in data structures that look like:

  typedef struct esl_dsqdata_chunk_s {
  int64_t   i0;           // Chunk contains sequences i0..i0+N-1 from the database, 0-offset
  int       N;            // Chunk contains N sequences
  ESL_DSQ **dsq;          // Pointers to each of the N sequences
  char    **name;         // Names, \0 terminated.  Ptr into <metadata> buffer.
  char    **acc;          // Optional accessions, \0 terminated;   "\0" if none.
  char    **desc;         // Optional descriptions, \0 terminated; "\0" if none
  int32_t  *taxid;        // NCBI taxonomy identifiers. (>=1 is a taxid; -1 means none)
  int64_t  *L;            // Sequence lengths, in residues. The unpacker figures these out.

  So, there's an outer loop that iterates over the chunks in the file, and an  inner loop that
  iterates over the sequences within each chuk.  Right now, it just counts the number of sequences
  and residues in the file, but this is where your code will go.
  */

  /* start variables for while loop below */
  int m; // index into pattern_counts
  int bad_letter; // bool to track if idx includes bad letters
  uint64_t num_residues = 0; // total number of residues in test set
  uint8_t *the_sequence; // sequences are arrays of 8-bit unsigned integers, one byte/residue
  int num_bad_letters = 0; // keep track of how many non-amino-acid residues there are in dataset
  int num_aa_in_seq = 0; // number of probabilities used in a sequence (not exactly equal to aa but close)
  int num_aa_in_testset = 0; // for normalization of total score
  double test_set_score = 0;

  /* start counters to keep track of scores */
  double *seq_scores = NULL;
  int scores_array_size = 600000;
  ESL_ALLOC(seq_scores, scores_array_size*8); // can do alloc realloc for this later

  /*cycle through each chunk in the file, each seq in chunk, and get score */
  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK) {

      // iterate over sequences in chunk
    	for (int i = 0; i < chu->N; i++) {
      		the_sequence = chu->dsq[i]; // Get the sequence
          num_aa_in_seq = 0;

          if(num_sequences >= scores_array_size){
              scores_array_size = scores_array_size * 2;
              ESL_REALLOC(seq_scores, scores_array_size*8);
            }

          // keep track of residues at the beginning
      		for(int j = 1; j < furthest_back+1 && j <= chu->L[i]; j++){
            // printf("hi %d\n",the_sequence[j]);
              if(the_sequence[j] > 20 && the_sequence[j] < 27){num_bad_letters++;}
              if(the_sequence[j] != 20 && the_sequence[j] < 27){num_residues++;}
              if(the_sequence[j] < 20){
                  num_aa_in_seq++;
                  m = the_sequence[j];
                  seq_scores[num_sequences] += log(marginal_probabs[m])/log(2);
                  test_set_score += log(marginal_probabs[m])/log(2);
                  // printf("inside %d: %lf and log: %lf and m=%d\n",num_sequences,marginal_probabs[m],log(marginal_probabs[m])/log(2),m);
              }
      		}

          // iterate through residues in sequence
      	 	for(int j = furthest_back+1; j <= chu->L[i]; j++) {
        			num_residues++;
        			int m_last = the_sequence[j];
        			m = 0;
        			if(m_last > 20 && m_last < 27){num_bad_letters++;}
              if(m_last > 19){continue;} // if residue is not one of 20 amino acids, skip

        			int bad_letter = 0; // keep track of if we find a non-amino acid letter

        			// iterate over lookback indeces and update m
              for(int k = look_back_idx - 1; k >= 0; k--) {
                  if(the_sequence[j-look_back_indeces[k]] > 19){
                      bad_letter = 1; // found non-amino acid letter
                      break; // if residue is not one of 20 amino acids skip it
                  }
                  // m is the idx into pattern counts for the specific lookback of this residue
                  m = m*20 + the_sequence[j-look_back_indeces[k]];
              }

        			// add final residue to m
        			m = m*20 + m_last;

        			// use m to get probability from table
              if(bad_letter == 0){
                  num_aa_in_seq++;
                  seq_scores[num_sequences] += log(pattern_probabs[m])/log(2);
                  test_set_score += log(pattern_probabs[m])/log(2);
                  // printf("inside %d: %lf and log: %lf and m=%d\n",num_sequences,pattern_probabs[m],log(pattern_probabs[m])/log(2),m);
              }
      		}
          seq_scores[num_sequences] = seq_scores[num_sequences] / num_aa_in_seq;
          // printf("sequence length: %d\n", num_aa_in_seq);
      		num_sequences++;
          num_aa_in_testset += num_aa_in_seq;

      }
      esl_dsqdata_Recycle(dd, chu);
  }

  // printf("total test set score: %lf\n", test_set_score / num_aa_in_testset);
  printf("%lf\n", test_set_score / num_aa_in_testset);

  // for(int i = 0; i < 2; i++){
  //     printf("%f\n", seq_scores[i]);
  // }

  // printf("Saw %d sequences, %llu amino acid residues\n", num_sequences, num_residues);
  // printf("Number of non-amino-acid residues in dataset: %d\n", num_bad_letters);
  // printf("Number of amino-acid residues in dataset: %d\n", num_aa_in_testset);

  esl_dsqdata_Close(dd);
  esl_getopts_Destroy(go);
  exit(0);

  ERROR:
  esl_fail("Unable to allocate memory\n", NULL);
}
