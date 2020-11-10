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

   int num_past_residues = 0;
   int count_array_size = 1;
   int furthest_back = 0;

/* Count the number of ones in history_vector, which is the number of resdues you're
  conditioning on, and compute the size of the count array as 20 ^ number of residues
  tracked */ 

   for(i = 0; i < 32; i++){
     if(history_vector & 1){ // low bit is set
      num_past_residues++;
      count_array_size = count_array_size * 20;
      furthest_back = i + 1;
    }
    history_vector = history_vector >> 1;
   }

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

  uint64_t num_residues = 0;
  uint8_t *the_sequence; // Sequences are arrays of 8-bit unsigned integers,
  //one byte/residue
  while (( status = esl_dsqdata_Read(dd, &chu)) == eslOK)  
    { // Iterate over chunks in file
      for (i = 0; i < chu->N; i++) { // Iterate over sequences in chunk
    	  the_sequence = chu->dsq[i]; // Get the sequence 
        for(int j = 0; j < chu->L[i]; j++){ // Iterate through residues in sequence
          num_residues++;    
      }
    num_sequences++;
	  } 
      esl_dsqdata_Recycle(dd, chu);
    }
    printf("Conditioning on %d residues\n", num_past_residues);
    printf("Saw %d sequences, %lu residues\n", num_sequences,num_residues);
  esl_dsqdata_Close(dd);
  esl_getopts_Destroy(go);
  exit(0);

  ERROR:
  esl_fail("Unable to allocate memory\n", NULL);
}





